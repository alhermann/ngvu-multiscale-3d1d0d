#ifndef __CREATECLOUDS_H__
#define __CREATECLOUDS_H__

////////////////////////////////////////////////////////////////////////////////
//
// kd_tree_interface.h
//
// Interface to libkdtree++ library
//
// Copyright (c) Farshid Mossaiby, 2016, 2017
//
////////////////////////////////////////////////////////////////////////////////

// Includes, system

#include <algorithm>
#include <cmath>
#include "kdtree++/kdtree.hpp"
#include <stdexcept>

typedef void (*build_cloud_callback)(size_t i, size_t &min_node_node, double &r,
                                     double &dr);

class Node;

Node *ref_node;

class Node {
 public:
  typedef double value_type;

  double xyz[3];
  int index;

  double operator[](int n) const { return xyz[n]; }

  double distance(const Node &node) const {
    double x = xyz[0] - node.xyz[0];
    double y = xyz[1] - node.xyz[1];
    double z = xyz[2] - node.xyz[2];

    return std::max(std::max(std::abs(x), std::abs(y)), std::abs(z));
  }

  double distance_euclidean(const Node &node) const {
    double x = xyz[0] - node.xyz[0];
    double y = xyz[1] - node.xyz[1];
    double z = xyz[2] - node.xyz[2];

    return std::sqrt(x * x + y * y + z * z);
  }
};

bool NodeComp(const Node &i, const Node &j) {
  return i.distance_euclidean(*ref_node) < j.distance_euclidean(*ref_node);
}

template <typename double_matrix_type, typename int_matrix_type,
          typename int_vector_type>
void build_clouds(int min_node_no, int max_node_no, double r, double dr,
                  double_matrix_type &nodes, double_matrix_type &nodes2,
                  int_matrix_type &clouds, int_vector_type &cloud_lengths,
                  int &max_cloud_length) {
  KDTree::KDTree<3, Node> tree;

  ref_node = new Node;

  int size = nodes.N();
  int size2 = nodes2.N();

  clouds.resize(size, max_node_no);
  cloud_lengths.resize(size);

  for (int i = 0; i < size2; i++) {
    Node node;

    node.xyz[0] = nodes2[i][0];
    node.xyz[1] = nodes2[i][1];
    node.xyz[2] = nodes2[i][2];

    node.index = i;

    tree.insert(node);
  }

  tree.optimize();

  std::vector<Node> found_nodes;

  max_cloud_length = 0;

  for (int i = 0; i < size; i++) {
    ref_node->xyz[0] = nodes[i][0];
    ref_node->xyz[1] = nodes[i][1];
    ref_node->xyz[2] = nodes[i][2];

    double current_r = r;

    found_nodes.reserve(max_node_no);
    int found_node_no = 0;

    do {
      found_nodes.clear();

      tree.find_within_range(
          *ref_node, current_r,
          std::back_insert_iterator<std::vector<Node> >(found_nodes));

      found_node_no = 0;

      for (int j = 0; j < found_nodes.size(); j++) {
        if (found_nodes[j].distance_euclidean(*ref_node) <= current_r) {
          found_node_no++;
        }
      }

      if (found_node_no > max_node_no) {
        throw std::runtime_error(
            "found_node_no > max_node_no; increase max_node_no.");
      }

      if (found_node_no < min_node_no) {
        current_r += dr;
      }

    } while (found_node_no < min_node_no);

    std::sort(found_nodes.begin(), found_nodes.end(), NodeComp);

    int current_node = 0;

    for (int j = 0; j < found_nodes.size(); j++) {
      if (found_nodes[j].distance_euclidean(*ref_node) <= current_r) {
        clouds[i][current_node++] = found_nodes[j].index;
      }
    }

    for (int j = found_node_no; j < max_node_no; j++) {
      clouds[i][j] = -1;
    }

    // std::cout << found_node_no << std::endl;

    cloud_lengths[i] = found_node_no;
    max_cloud_length = std::max(max_cloud_length, found_node_no);
  }
}


#endif  // __CREATECLOUDS_H__
