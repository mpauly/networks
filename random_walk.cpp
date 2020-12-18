#include "Snap.h"
#include "walker.h"
#include <fstream>
#include <functional>
#include <iostream>
#include <unistd.h>

const double DIFFUSION_CONSTANT = 0.5;

int main(int argc, char *argv[]) {
  int start_node = 0;
  int sigma_max = 100;
  std::string dimension_file;

  // parse arguments
  int option;
  while ((option = getopt(argc, argv, ":s:l:d:")) != -1) { // get option from the getopt() method
    switch (option) {
    case 's':
      start_node = atoi(optarg);
      break;
    case 'l':
      sigma_max = atoi(optarg);
      break;
    case 'd':
      dimension_file = optarg;
      break;
    case '?':
      printf("unknown option: %c\n", optopt);
      break;
    }
  }

  int nr_of_files = 0;
  std::string graph_filename;
  for (; optind < argc; optind++) { // when some extra arguments are passed
    graph_filename = argv[optind];
    nr_of_files++;
  }

  if (nr_of_files != 1) {
    std::cout << "Please specify exactly one file" << std::endl;
    return 1;
  }

  if (dimension_file.empty()) {
    dimension_file = graph_filename;
    dimension_file = dimension_file.replace(dimension_file.find("graphs/"), sizeof("graphs/") - 1, "data/dim_");
    std::cout << "- no output file given - will write to " << dimension_file << std::endl;
  }

  // =====  load the graph  ========
  typedef PUNGraph PGraph; // undirected graph
  std::cout << "- Loading graph file " << graph_filename << std::endl;

  PUNGraph G;
  try {
    TFIn FIn(graph_filename.c_str());
    G = TUNGraph::Load(FIn);
  } catch (...) {
    std::cerr << "Error reading file " << graph_filename << std::endl;
    return 1;
  }

  // ==== do the walk ====

  std::cout << "- Starting to walk at node " << start_node << " for " << sigma_max << " steps" << std::endl;
  std::function<void(int)> progress_monitor = [](int sigma) {
    if (sigma % 50 == 0)
      std::cout << "\t - sigma = " << sigma << std::endl;
  };

  auto walk_dimensions = spectralDimensionAtNode(G, start_node, sigma_max, progress_monitor, DIFFUSION_CONSTANT);

  // ===  write to file  ==
  std::cout << "- Writing results to file " << dimension_file << std::endl;

  std::ofstream dimfile;
  dimfile.open(dimension_file.c_str());
  dimfile << "# dimensions for graph " << graph_filename << std::endl;
  dimfile << "# walk length: " << sigma_max << std::endl;
  dimfile << "# start node: " << start_node << std::endl;
  for (int sigma = 0; sigma < walk_dimensions.size(); sigma++) {
    dimfile << sigma << "\t";
    dimfile << std::fixed << std::setprecision(12);
    dimfile << walk_dimensions[sigma] << "\n";
  }
  dimfile.close();

  return 0;
}
