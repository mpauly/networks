#include "Snap.h"
#include "walker.h"
#include <fstream>
#include <functional>
#include <iostream>
#include <unistd.h>

const double DIFFUSION_CONSTANT = 0.5;

void print_usage() {
  std::cout << "usage: ./random_walk.x [-w walker_nr] [-s start_node] [-l length] [-d dimensionfile] graph_filename"
            << std::endl
            << "\twalker_nr: Number of walkers that will be run, defaults to 1" << std::endl
            << "\tstart_node: The node to start at - can only be given for exactly one walker" << std::endl
            << "\tlength: Length of the random walk, defaults to 100" << std::endl
            << "\tdimensionfile: Output file to write the result to, defaults to an automatically generated filename "
               "in data/dim_"
            << std::endl
            << "\tgraph_filename: The graph to run the random walker on" << std::endl;
}

int main(int argc, char *argv[]) {
  int start_node = 0;
  int sigma_max = 100;
  int nr_of_walkers = 1;
  std::string dimension_file;

  // parse arguments
  int option;
  while ((option = getopt(argc, argv, ":s:l:w:d:")) != -1) { // get option from the getopt() method
    switch (option) {
    case 's':
      start_node = atoi(optarg);
      break;
    case 'l':
      sigma_max = atoi(optarg);
      break;
    case 'w':
      nr_of_walkers = atoi(optarg);
      break;
    case 'd':
      dimension_file = optarg;
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
    graph_filename = argv[optind];
    nr_of_files++;
  }

  if (nr_of_files != 1) {
    std::cout << "Please specify exactly one file" << std::endl;
    print_usage();
    return 1;
  }

  if (dimension_file.empty()) {
    dimension_file = graph_filename;
    dimension_file = dimension_file.replace(dimension_file.find("graphs/"), sizeof("graphs/") - 1, "data/dim_");
    std::cout << "- no output file given - will write to " << dimension_file << std::endl;
  }

  if (nr_of_walkers > 1 && start_node > 0) {
    std::cerr << "Can't give options for number of walkers and start nodes at the same time" << '\n';
    print_usage();
    return 1;
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
  std::cout << "- Writing results to file " << dimension_file << std::endl;

  std::ofstream dimfile;
  dimfile.open(dimension_file.c_str());
  dimfile << "# dimensions for graph " << graph_filename << std::endl;
  dimfile << "# nr of walks: " << nr_of_walkers << std::endl;
  dimfile << "# walk length: " << sigma_max << std::endl;
  dimfile << "# format: start_node sigma d_spec" << std::endl;

#pragma omp parallel for
  for (int walk = 0; walk < nr_of_walkers; walk++) {
    int walker_start_node;
    if (nr_of_walkers > 1) {
      walker_start_node = G->GetRndNId();
    } else {
      walker_start_node = start_node;
    }

#pragma omp critical
    {
      std::cout << "- Starting walk " << walk << " at node " << start_node << " for " << sigma_max << " steps"
                << std::endl;
    }
    std::function<void(int)> progress_monitor = [walk](int sigma) {
      if (sigma % 50 == 0)
        std::cout << "\t - w" << walk << " sigma = " << sigma << std::endl;
    };

    auto walk_dimensions =
        spectralDimensionAtNode(G, walker_start_node, sigma_max, progress_monitor, DIFFUSION_CONSTANT);

// ===  write to file  ==
#pragma omp critical
    {
      for (int sigma = 0; sigma < walk_dimensions.size(); sigma++) {
        dimfile << walker_start_node << "\t" << sigma << "\t";
        dimfile << std::fixed << std::setprecision(12);
        dimfile << walk_dimensions[sigma] << "\n";
      }
      std::cout << "- Writing walk " << walk << " to file" << std::endl;
    }
  }
  dimfile.close();
  std::cout << "- Finished walking..." << std::endl;
  return 0;
}
