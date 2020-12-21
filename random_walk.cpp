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

  if (nr_of_walkers > 1 && start_node > 0) {
    std::cerr << "Can't give options for number of walkers and start nodes at the same time" << '\n';
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

  for (int walk = 0; walk < nr_of_walkers; walk++) {
    int walker_start_node;
    if (nr_of_walkers > 1) {
      walker_start_node = G->GetRndNId();
    } else {
      walker_start_node = start_node;
    }

    std::cout << "- Starting walk " << walk << " at node " << start_node << " for " << sigma_max << " steps"
              << std::endl;
    std::function<void(int)> progress_monitor = [walk](int sigma) {
      if (sigma % 50 == 0)
        std::cout << "\t - w" << walk << " sigma = " << sigma << std::endl;
    };

    auto walk_dimensions =
        spectralDimensionAtNode(G, walker_start_node, sigma_max, progress_monitor, DIFFUSION_CONSTANT);

    // ===  write to file  ==
    for (int sigma = 0; sigma < walk_dimensions.size(); sigma++) {
      dimfile << walker_start_node << "\t" << sigma << "\t";
      dimfile << std::fixed << std::setprecision(12);
      dimfile << walk_dimensions[sigma] << "\n";
    }
  }
  dimfile.close();

  return 0;
}
