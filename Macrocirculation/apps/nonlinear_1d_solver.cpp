////////////////////////////////////////////////////////////////////////////////
//  Copyright (c) 2021 Andreas Wagner.
//
//  Distributed under the Boost Software License, Version 1.0. (See accompanying
//  file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
////////////////////////////////////////////////////////////////////////////////

#include <chrono>
#include <cxxopts.hpp>
#include <memory>

#include "macrocirculation/communication/mpi.hpp"
#include "macrocirculation/csv_vessel_tip_writer.hpp"
#include "macrocirculation/dof_map.hpp"
#include "macrocirculation/embedded_graph_reader.hpp"
#include "macrocirculation/explicit_nonlinear_flow_solver.hpp"
#include "macrocirculation/graph_csv_writer.hpp"
#include "macrocirculation/graph_partitioner.hpp"
#include "macrocirculation/graph_pvd_writer.hpp"
#include "macrocirculation/graph_storage.hpp"
#include "macrocirculation/interpolate_to_vertices.hpp"
#include "macrocirculation/quantities_of_interest.hpp"
#include "macrocirculation/vessel_formulas.hpp"
#include "macrocirculation/rcr_estimator.hpp"
#include <nlohmann/json.hpp>

namespace mc = macrocirculation;

constexpr std::size_t degree = 2;

void output_flows(const std::string &filepath, const mc::GraphStorage &graph, const mc::FlowData &flows, double t_end, double t_start_averaging) {
  using json = nlohmann::json;

  json j;

  auto vertices_list = json::array();

  for (auto t : flows.flows) {
    auto vertex = graph.get_vertex(t.first);
    auto edge = graph.get_edge(vertex->get_edge_neighbors()[0]);
    const auto r = edge->get_physical_data().radius;
    auto flow = t.second;
    json vessel_obj = {
      {"vertex_name", vertex->get_name()},
      {"average_flow", flow / (t_end - t_start_averaging)},
      {"radius", r}
    };
    vertices_list.push_back(vessel_obj);
  }

  j["vertex_flow_data"] = vertices_list;
  j["unit"] = {
    {"average_flow", "cm^3/s"},
    {"radius", "cm"}
  };
  j["comments"] = "Formula for flows: Q_cap = (r_cap / r_art)^gamma * Q_art";

  std::ofstream f(filepath, std::ios::out);
  f << j.dump(1);
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  {
    cxxopts::Options options(argv[0], "Nonlinear 1D solver");
    options.add_options()                                                                                                                                           //
      ("mesh-file", "path to the input file", cxxopts::value<std::string>()->default_value("./data/1d-meshes/33-vessels.json"))                               //
      ("boundary-file", "path to the file for the boundary conditions", cxxopts::value<std::string>()->default_value(""))                                           //
      ("output-directory", "directory for the output", cxxopts::value<std::string>()->default_value("./output/"))                                                   //
      ("inlet-name", "the name of the inlet", cxxopts::value<std::string>()->default_value("cw_in"))                                                                //
      ("heart-amplitude", "the amplitude of a heartbeat", cxxopts::value<double>()->default_value("485.0"))                                                         //
      ("tau", "time step size", cxxopts::value<double>()->default_value(std::to_string(2.5e-4 / 16.)))                                                              //
      ("tau-out", "time step size for the output", cxxopts::value<double>()->default_value("1e-2"))                                                                 //
      ("t-start-averaging", "Time when to start averaging flows", cxxopts::value<double>()->default_value("0"))                                                                             //
      ("t-end", "Endtime for simulation", cxxopts::value<double>()->default_value("0.01"))                                                                             //
      ("h,help", "print usage");
    options.allow_unrecognised_options(); // for petsc
    auto args = options.parse(argc, argv);
    if (args.count("help")) {
      std::cout << options.help() << std::endl;
      exit(0);
    }
    if (!args.unmatched().empty()) {
      std::cout << "The following arguments were unmatched: " << std::endl;
      for (auto &it : args.unmatched())
        std::cout << " " << it;
      std::cout << "\nAre they part of petsc or a different auxillary library?" << std::endl;
    }

    // create_for_node the ascending aorta
    auto graph = std::make_shared<mc::GraphStorage>();

    mc::EmbeddedGraphReader graph_reader;
    graph_reader.append(args["mesh-file"].as<std::string>(), *graph);

    // read in other data
    auto boundary_file_path = args["boundary-file"].as<std::string>();
    if (!boundary_file_path.empty()) {
      std::cout << "Using separate file at " << boundary_file_path << " for boundary conditions." << std::endl;
      graph_reader.set_boundary_data(boundary_file_path, *graph);
    }

    graph->find_vertex_by_name(args["inlet-name"].as<std::string>())->set_to_inflow_with_fixed_flow(mc::heart_beat_inflow(args["heart-amplitude"].as<double>()));

    // set_0d_tree_boundary_conditions(graph, "bg_");
    graph->finalize_bcs();

    // mc::naive_mesh_partitioner(*graph, MPI_COMM_WORLD);
    mc::flow_mesh_partitioner(MPI_COMM_WORLD, *graph, degree);

    auto dof_map_flow = std::make_shared<mc::DofMap>(graph->num_vertices(), graph->num_edges());
    dof_map_flow->create(MPI_COMM_WORLD, *graph, 2, degree, false);

    const double t_start_averaging = args["t-start-averaging"].as<double>();

    const double t_end = args["t-end"].as<double>();
    const std::size_t max_iter = 160000000;

    const auto tau = args["tau"].as<double>();
    const auto tau_out = args["tau-out"].as<double>();

    const auto output_interval = static_cast<std::size_t>(tau_out / tau);
    std::cout << "tau = " << tau << ", tau_out = " << tau_out << ", output_interval = " << output_interval << std::endl;

    // configure solver
    auto flow_solver = std::make_shared<mc::ExplicitNonlinearFlowSolver>(MPI_COMM_WORLD, graph, dof_map_flow, degree);
    flow_solver->use_ssp_method();

    std::vector<mc::Point> points;
    std::vector<double> Q_vertex_values;
    std::vector<double> A_vertex_values;
    std::vector<double> p_total_vertex_values;
    std::vector<double> p_static_vertex_values;
    std::vector<double> c_vertex_values;
    std::vector<double> vessel_ids;

    // vessels ids do not change, thus we can precalculate them
    mc::fill_with_vessel_id(MPI_COMM_WORLD, *graph, points, vessel_ids);

    mc::GraphCSVWriter csv_writer(MPI_COMM_WORLD, args["output-directory"].as<std::string>(), "abstract_vessels", graph);
    csv_writer.add_setup_data(dof_map_flow, flow_solver->A_component, "a");
    csv_writer.add_setup_data(dof_map_flow, flow_solver->Q_component, "q");
    csv_writer.setup();

    mc::GraphPVDWriter pvd_writer(MPI_COMM_WORLD, args["output-directory"].as<std::string>(), "abstract_vessels");
    mc::CSVVesselTipWriter vessel_tip_writer(MPI_COMM_WORLD, "output", "abstract_vessels_tips", graph, dof_map_flow);

    // output for 0D-Model:
    std::vector<double> list_t;
    std::map<size_t, std::vector<mc::Values0DModel>> list_0d;
    for (auto v_id : graph->get_active_vertex_ids(mc::mpi::rank(MPI_COMM_WORLD))) {
      if (graph->get_vertex(v_id)->is_windkessel_outflow())
        list_0d[v_id] = std::vector<mc::Values0DModel>{};
    }

    for (auto &v_id : graph->get_active_vertex_ids(mc::mpi::rank(MPI_COMM_WORLD))) {
      auto &vertex = *graph->get_vertex(v_id);
      if (vertex.is_vessel_tree_outflow()) {
        const auto &vertex_dof_map = dof_map_flow->get_local_dof_map(vertex);
        const auto &vertex_dofs = vertex_dof_map.dof_indices();
        auto &u = flow_solver->get_solution();

        for (size_t k = 0; k < vertex_dofs.size(); k += 1)
          u[vertex_dofs[k]] = 0. * 1.33333;
      }
    }

    double t = 0;

    const auto write_output = [&](){
      csv_writer.add_data("a", flow_solver->get_solution());
      csv_writer.add_data("q", flow_solver->get_solution());
      csv_writer.write(t);

      mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph, *dof_map_flow, 0, flow_solver->get_solution(), points, Q_vertex_values);
      mc::interpolate_to_vertices(MPI_COMM_WORLD, *graph, *dof_map_flow, 1, flow_solver->get_solution(), points, A_vertex_values);
      mc::calculate_total_pressure(MPI_COMM_WORLD, *graph, *dof_map_flow, flow_solver->get_solution(), points, p_total_vertex_values);
      mc::calculate_static_pressure(MPI_COMM_WORLD, *graph, *dof_map_flow, flow_solver->get_solution(), points, p_static_vertex_values);

      pvd_writer.set_points(points);
      pvd_writer.add_vertex_data("Q", Q_vertex_values);
      pvd_writer.add_vertex_data("A", A_vertex_values);
      pvd_writer.add_vertex_data("p_static", p_static_vertex_values);
      pvd_writer.add_vertex_data("p_total", p_total_vertex_values);
      pvd_writer.add_vertex_data("c", c_vertex_values);
      pvd_writer.add_vertex_data("vessel_id", vessel_ids);
      pvd_writer.write(t);

      vessel_tip_writer.write(t, flow_solver->get_solution());
    };

    write_output();

    double flow_solution_time = 0;
    size_t num_iteration = 0;

    const auto begin_t = std::chrono::steady_clock::now();

    mc::FlowIntegrator flow_integrator(graph);

    for (std::size_t it = 0; it < max_iter; it += 1) {
      auto start = std::chrono::high_resolution_clock::now();
      flow_solver->solve(tau, t);
      auto end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> elapsed = (end - start);
      flow_solution_time += elapsed.count();
      num_iteration += 1;

      t += tau;

      if (it % output_interval == 0) {
        std::cout << "iter = " << it << ", t = " << t << std::endl;

        write_output();
      }

      // break
      if (t > t_end + 1e-12)
        break;

      // add total flows
      if (t >= t_start_averaging)
        flow_integrator.update_flow(*flow_solver, tau);
    }

    const auto end_t = std::chrono::steady_clock::now();
    const auto elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end_t - begin_t).count();
    if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
    {
      std::cout << "time = " << elapsed_ms * 1e-6 << " s" << std::endl;
      std::cout << "total time flow solver = " << flow_solution_time << ", average = " << flow_solution_time / num_iteration << ", iteration = " << num_iteration << std::endl;
    }

    auto flows = flow_integrator.get_windkessel_outflow_data();

    output_flows("flows.json", *graph, flows, t_end, t_start_averaging);

    if (mc::mpi::rank(MPI_COMM_WORLD) == 0)
    {
      for (auto t : flows.flows) {
        auto vertex = graph->get_vertex(t.first);
        auto flow = t.second;

        std::cout << "vertex name = " << vertex->get_name()
                  << ", id = " << vertex->get_id()
                  << ", flow = " << flow
                  << ", average flow = " << flow / (t_end - t_start_averaging)
                  << std::endl;
      }
    }
  }

  MPI_Finalize();
}
