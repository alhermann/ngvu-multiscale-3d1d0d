cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# ****************************************************************************
# Project information
# ****************************************************************************
set(ProjectName "Macrocirculation")
set(ProjectLib "Lib${ProjectName}")
project(${ProjectName} CXX)

# Create executible in main directory of build directory
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/macrocirculation/)

# ****************************************************************************
# Create executable
# ****************************************************************************

# nonlinear flow solver line:
add_executable(${ProjectName}NonlinearFlowLine nonlinear_flow_line.cpp)
target_link_libraries(${ProjectName}NonlinearFlowLine ${ProjectLib})

# nonlinear flow solver bifurcation:
add_executable(${ProjectName}NonlinearFlowBifurcation nonlinear_flow_bifurcation.cpp)
target_link_libraries(${ProjectName}NonlinearFlowBifurcation ${ProjectLib})

# test the convergence:
add_executable(${ProjectName}ConvergenceStudy convergence_study.cpp)
target_link_libraries(${ProjectName}ConvergenceStudy ${ProjectLib})

# test the 1d network:
add_executable(${ProjectName}Nonlinear1DSolver nonlinear_1d_solver.cpp)
target_link_libraries(${ProjectName}Nonlinear1DSolver ${ProjectLib})
target_link_libraries(${ProjectName}Nonlinear1DSolver cxxopts)
target_link_libraries(${ProjectName}Nonlinear1DSolver nlohmann_json::nlohmann_json)

add_executable(${ProjectName}Calibration calibration.cpp)
target_link_libraries(${ProjectName}Calibration ${ProjectLib})
target_link_libraries(${ProjectName}Calibration cxxopts)

foreach (TargetName ${ProjectName}NonlinearFlowLine
        ${ProjectName}NonlinearFlowBifurcation
        ${ProjectName}ConvergenceStudy
        ${ProjectName}Nonlinear1DSolver
        )

    foreach (FolderName output
            output_darcy3d
            output_darcy_hetero
            output_part_perfusion_test
            output_voxels_to_mesh
            output_darcy_flow_reg_outlets
            output_full_1d0d3d_pkj
            output_3d_model
            output_multiscale_time_3d1d0d)
        add_custom_command(
                TARGET ${TargetName} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E make_directory
                $<TARGET_FILE_DIR:${TargetName}>/${FolderName})
    endforeach ()

    add_custom_command(
            TARGET ${TargetName} PRE_BUILD
            COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/../tools/visualization $<TARGET_FILE_DIR:${TargetName}>/visualization
            COMMENT "Symlinking folder: ${CMAKE_CURRENT_SOURCE_DIR}/../tools/visualization")

    add_custom_command(
            TARGET ${TargetName} PRE_BUILD
            COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/../data $<TARGET_FILE_DIR:${TargetName}>/data
            COMMENT "Copying symlinking folder: data")

    file(GLOB visualization_scripts RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}/../tools/visualization ${CMAKE_CURRENT_SOURCE_DIR}/../tools/visualization/*)
    foreach (f ${visualization_scripts})
        MESSAGE(${f})
        add_custom_command(
                TARGET ${TargetName} PRE_BUILD
                COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/../tools/visualization/${f} $<TARGET_FILE_DIR:${TargetName}>/${f}
                COMMENT "Symlinking : ${CMAKE_CURRENT_SOURCE_DIR}/../tools/visualization/${f} to $<TARGET_FILE_DIR:${TargetName}>/${f}")
    endforeach ()
endforeach ()
