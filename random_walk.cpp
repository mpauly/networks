#include "Snap.h"
#include "walker.h"
#include <fstream>
#include <functional>
#include <iostream>
#include <unistd.h>

void print_usage() {
  std::cout
      << "usage: ./random_walk.x [-a] [-w walker_nr] [-s start_node] [-l length] [-d diffusion_constant] [-o outfile] "
         "graph_filename"
      << std::endl
      << "  -a appending: Append to dimension file instead of overwriting it" << std::endl
      << "  -w walker_nr: Number of walkers that will be run, defaults to 1" << std::endl
      << "  -s start_node: The node to start at - can only be given for exactly one walker" << std::endl
      << "  -l length: Length of the random walk, defaults to 100" << std::endl
      << "  -d diffusion_constant delta: Diffusion constant to regularize osciallations in the spectral dimension, "
         "default: 0.5, (1-delta) is the probability to stay at the current node"
      << std::endl
      << "  -o outfile: Output file to write the result to, defaults to an automatically generated filename "
         "in data/dim_"
      << std::endl
      << "  graph_filename: The graph to run the random walker on" << std::endl;
}

int main(int argc, char *argv[]) {
  bool file_appending = false;
  int start_node = 0;
  int sigma_max = 100;
  int nr_of_walkers = 1;
  double diffusion_constant = 0.5;
  std::string dimension_file;

  // parse arguments
  int option;
  while ((option = getopt(argc, argv, "as:l:w:d:o:")) != -1) { // get option from the getopt() method
    switch (option) {
    case 'a':
      file_appending = true;
      break;
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
      diffusion_constant = atof(optarg);
      break;
    case 'o':
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

  auto filemode = std::ofstream::out;
  if (file_appending) {
    filemode = std::ofstream::out | std::ofstream::app;
    std::cout << "- Appending to existing file" << std::endl;
  }

  std::ofstream dimfile;
  dimfile.open(dimension_file.c_str(), filemode);
  if (!file_appending) {
    dimfile << "# dimensions for graph " << graph_filename << std::endl;
    dimfile << "# nr of walks: " << nr_of_walkers << std::endl;
    dimfile << "# walk length: " << sigma_max << std::endl;
    dimfile << "# diffusion const: " << diffusion_constant << std::endl;
    dimfile << "# format: start_node sigma d_spec" << std::endl;
  }

  // Determine starting positions - this is outside of the parallel section to avoid complications regarding the
  // generation of random numbers in different threads
  std::vector<int> walker_start_nodes(nr_of_walkers);
  for (int walk = 0; walk < nr_of_walkers; walk++) {
    if (nr_of_walkers > 1) {
      walker_start_nodes[walk] = G->GetRndNId();
    } else {
      walker_start_nodes[walk] = start_node;
    }
  }

#pragma omp parallel for
  for (int walk = 0; walk < nr_of_walkers; walk++) {
    int walker_start_node;
    walker_start_node = walker_start_nodes[walk];

#pragma omp critical
    {
      std::cout << "- Starting walk " << walk << " at node " << walker_start_node << " for " << sigma_max
                << " steps with diffusion constant " << diffusion_constant << std::endl;
    }
    std::function<void(int)> progress_monitor = [walk](int sigma) {
      if (sigma % 50 == 0)
        std::cout << "\t - w" << walk << " sigma = " << sigma << std::endl;
    };

    auto walk_dimensions =
        spectralDimensionAtNode(G, walker_start_node, sigma_max, progress_monitor, diffusion_constant);

// ===  write to file  ==
// We drop the first data point because one cannot compute a derivative at that data point
#pragma omp critical
    {
      for (int sigma = 1; sigma < walk_dimensions.size(); sigma++) {
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
