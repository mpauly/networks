#include "Snap.h"
#include "walker/base.h"

#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>

namespace fs = std::filesystem;

// gives different shades of blue as a function of value in [0,1]
TStr get_color(double value) {
  char hex[9];
  sprintf(hex, "#%02x%02x%02x%02x", 0, 0, 255, (int)(255 * value));
  return TStr(hex);
}

void print_usage() {
  std::cout
      << "usage: ./visualize_walk.x [-s start_node] [-l length] [-d diffusion_constant] graph_name" << std::endl
      << "  -s start_node: The node to start at - can only be given for exactly one walker" << std::endl
      << "  -l length: Length of the random walk, defaults to 100" << std::endl
      << "  -d diffusion_constant delta: Diffusion constant to regularize osciallations in the spectral dimension, "
         "default: 0.5, (1-delta) is the probability to stay at the current node"
      << std::endl
      << "  graph_filename: The graph to run the random walker on in format roadnet_pa" << std::endl;
}

int main(int argc, char *argv[]) {

  int frames = 20;
  int sigma_max = 20;
  int start_node = 0;
  double diffusion_constant = 0.5;
  std::string graph_name;

  int option;
  while ((option = getopt(argc, argv, "f:s:l:d:?")) != -1) { // get option from the getopt() method
    switch (option) {
    case 'f':
      frames = atoi(optarg);
      break;
    case 's':
      start_node = atoi(optarg);
      break;
    case 'l':
      sigma_max = atoi(optarg);
      break;
    case 'd':
      diffusion_constant = atof(optarg);
      break;
    case '?':
      std::cout << "Unknown option: " << optopt << std::endl;
      print_usage();
      return 1;
    }
  }

  int nr_of_files = 0;
  std::string graph_filename;
  for (; optind < argc; optind++) {
    graph_name = argv[optind];
    nr_of_files++;
  }

  if (nr_of_files != 1) {
    std::cout << "Please specify exactly one file" << std::endl;
    print_usage();
    return 1;
  }

  typedef PUNGraph PGraph; // undirected graph
  PUNGraph G;
  try {
    std::string graph_filename = "graphs/" + graph_name + ".dat";
    TFIn FIn(graph_filename.c_str());
    G = TUNGraph::Load(FIn);
  } catch (...) {
    std::cerr << "Error reading file " << std::endl;
    return 1;
  }

  walker::RandomWalk walk = walker::setupRandomWalk(start_node, diffusion_constant);

  std::string dirname = "plots/frames/" + graph_name;
  if (!fs::is_directory(dirname) || !fs::exists(dirname)) {
    fs::create_directory(dirname);
  }

  std::cout << " == Walking on graph " << graph_name << " == " << std::endl;
  std::cout << "  Graph has " << G->GetNodes() << " nodes and " << G->GetEdges() << " edges" << std::endl;

  for (int frame = 0; frame < frames; frame++) {
    // make a snapshot
    std::cout << " == Frame " << frame << " == " << std::endl;
    TIntV visisted_nodes;
    TIntStrH NIdColorH;

    walk.lvl_probabilities.GetKeyV(visisted_nodes);

    std::cout << "  " << visisted_nodes.Len() << " Nodes on this level" << std::endl;

    // here we are just rescaling the values to the interval [0,1] in order to have more vivid colors
    double min_val = 2.0;
    double max_val = 0.0;

    for (auto it = walk.lvl_probabilities.BegI(); it < walk.lvl_probabilities.EndI(); it++) {
      if (it.GetDat() > max_val)
        max_val = it.GetDat();
      if (it.GetDat() < min_val)
        min_val = it.GetDat();
    }

    double scaling = 1.0 / (max_val - min_val);

    for (auto it = walk.lvl_probabilities.BegI(); it < walk.lvl_probabilities.EndI(); it++) {
      NIdColorH.AddDat(it.GetKey(), get_color(scaling * (it.GetDat() - min_val)));
    }

    PGraph SubG = TSnap::GetSubGraph(G, visisted_nodes);
    std::string framename = dirname + "/frame_" + std::to_string(frame) + ".png";
    TSnap::DrawGViz<PGraph>(SubG, gvlNeato, framename.c_str(), "", true, NIdColorH);
    // walk on
    walker::progressRandomWalk(G, walk, sigma_max / frames);
  }

  return 0;
}