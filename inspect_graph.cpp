#include "Snap.h"
#include "walker/base.h"
#include "walker/consts.h"
#include "walker/io.h"
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <unistd.h>

namespace fs = std::filesystem;

struct WalkConfig {
  int sigma_max = 10000;
  int nr_of_walkers = 1;
  int walker_min_id = 0;
  int interval = 500;
  int limit_length = -1;
  int export_prob_dist = false;
  double diffusion_constant = 0.5;
  std::string graph;
  std::string graph_file;
};

void print_usage() {
  std::cout << "usage: ./inspect_graph.x  "
               "graph"
            << std::endl
            << "Displays all relevant info on a graph" << std::endl;
}

template <class Graph> int process_network_or_graph(WalkConfig config) {
  std::string walk_dirname = walker::WALK_DIR + config.graph + '/';
  TPt<Graph> G;
  try {
    TFIn FIn(config.graph_file.c_str());
    G = Graph::Load(FIn);
  } catch (...) {
    std::cerr << "Error reading file " << config.graph_file << std::endl;
    return 1;
  }

  std::cout << "  Nodes: " << G->GetNodes() << std::endl;
  std::cout << "  Edges: " << G->GetEdges() << std::endl;
  return 0;
}

int main(int argc, char *argv[]) {
  WalkConfig config;
  // parse arguments
  int option;
  while ((option = getopt(argc, argv, "")) != -1) { // get option from the getopt() method
    switch (option) {
    case '?':
      std::cout << "Unknown option: " << optopt << std::endl;
      print_usage();
      return 1;
    }
  }

  int nr_of_files = 0;
  for (; optind < argc; optind++) {
    config.graph = argv[optind];
    nr_of_files++;
  }

  if (nr_of_files != 1) {
    std::cout << "Please specify exactly one graph" << std::endl;
    print_usage();
    return 1;
  }

  config.graph_file = walker::GRAPH_DIR + config.graph + ".dat";
  if (fs::exists(config.graph_file)) {
    std::cout << "Info for graph " << config.graph_file << std::endl;
    return process_network_or_graph<TUNGraph>(config);
  }

  config.graph_file = walker::NETWORK_DIR + config.graph + ".dat";
  if (fs::exists(config.graph_file)) {
    std::cout << "Info for network " << config.graph_file << std::endl;
    return process_network_or_graph<TIntNEDNet>(config);
  }

  std::cout << "- Neither graph or network found - please check that the graph " << config.graph << " exists!"
            << std::endl;
  return 1;
}
