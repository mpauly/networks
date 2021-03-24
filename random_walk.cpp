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
#include <string>
#include <unistd.h>

namespace fs = std::filesystem;

void print_usage() {
  std::cout
      << "usage: ./random_walk.x [-a] [-e] [-c] [-w walker_nr] [-W walker_min_id] [-s start_node] [-l length] [-d "
         "diffusion_constant] [-o "
         "outfile] graphname"
      << std::endl
      << "  -a appending: Append to dimension file instead of overwriting it" << std::endl
      << "  -e export: Export final state of the walker such that the walk can be continued in the future" << std::endl
      << "  -c continue: Continue -w existing walks for -l steps - ignores -s and -d - this also implies -e"
      << std::endl
      << "  -w walker_nr: Number of walkers that will be run, defaults to 1" << std::endl
      << "  -W walker_min_id: the numbering of the walkers starts at this id - this allows to add random walks to an "
         "existing ensemble"
      << std::endl
      << "  -s start_node: The node to start at - can only be given for exactly one walker" << std::endl
      << "  -l length: Length of the random walk, defaults to 100" << std::endl
      << "  -d diffusion_constant delta: Diffusion constant to regularize osciallations in the spectral dimension, "
         "default: 0.5, (1-delta) is the probability to stay at the current node"
      << std::endl
      << "  -o outfile: Output file to write the result to, defaults to an automatically generated filename "
      << std::endl
      << "  graphname: The graph to run the random walker on" << std::endl;
}

struct WalkConfig {
  bool file_appending = false;
  bool export_final_state = false;
  bool continue_walk = false;
  bool random_start_nodes = true;
  std::vector<int> start_node;
  int sigma_max = 100;
  int nr_of_walkers = 1;
  int walker_min_id = 0;
  double diffusion_constant = 0.5;
  std::string dimension_file;
  std::string graph;
  std::string graph_file;
  std::string walk_dir;
};

template <class Graph> int process_network_or_graph(WalkConfig config) {
  TPt<Graph> G;
  try {
    TFIn FIn(config.graph_file.c_str());
    G = Graph::Load(FIn);
  } catch (...) {
    std::cerr << "Error reading file " << config.graph_file << std::endl;
    return 1;
  }

  std::cout << "  Input graph has " << G->GetNodes() << " nodes and " << G->GetEdges() << " edges";
  if (!TSnap::IsWeaklyConn(G)) {
    std::cout << std::endl;
    std::cout << "============================================" << std::endl;
    std::cout << "== WARNING: Graph is not weakly connected ==" << std::endl;
    std::cout << "============================================" << std::endl;
    TIntPrV counts;
    TSnap::GetWccSzCnt(G, counts);
    for (int i = 0; i < counts.Len(); i++) {
      std::cout << " - " << counts[i].GetVal2() << " component(s) with " << counts[i].GetVal1() << " nodes"
                << std::endl;
    }
    return 1;
  } else {
    std::cout << " and is connected" << std::endl;
  }

  std::cout << "- Writing results to file " << config.dimension_file << std::endl;

  auto filemode = std::ofstream::out;
  if (config.file_appending) {
    filemode = std::ofstream::out | std::ofstream::app;
    std::cout << "- Appending to existing file" << std::endl;
  }

  std::ofstream dimfile;
  dimfile.open(config.dimension_file.c_str(), filemode);
  if (!config.file_appending) {
    dimfile << "# dimensions for graph " << config.graph << std::endl;
    dimfile << "# nr of walks: " << config.nr_of_walkers << std::endl;
    dimfile << "# walk length: " << config.sigma_max << std::endl;
    dimfile << "# diffusion const: " << config.diffusion_constant << std::endl;
    dimfile << "# format: start_node sigma d_spec" << std::endl;
  }

  // Determine starting positions - this is outside of the parallel section to avoid complications regarding the
  // generation of random numbers in different threads
  // in case walker_min_id > 0 we are generating more than the needed starting positions - the snap RNG is
  // deterministic, i.e. we are first reproducing the start positions of existing walks and then generating new starting
  // positions
  std::vector<int> walker_start_nodes;
  walker_start_nodes.reserve(config.walker_min_id + config.nr_of_walkers);
  if (config.random_start_nodes && !config.continue_walk) {
    while (walker_start_nodes.size() < config.walker_min_id + config.nr_of_walkers) {
      int candidate = G->GetRndNId();
      // is that start node already in our pool of start nodes?
      if (std::find(walker_start_nodes.begin(), walker_start_nodes.end(), candidate) == walker_start_nodes.end())
        walker_start_nodes.push_back(candidate);
    }
  } else {
    for (int i = 0; i < config.nr_of_walkers; i++)
      walker_start_nodes[i] = config.start_node[i];
  }

#pragma omp parallel for
  for (int walk = config.walker_min_id; walk < config.walker_min_id + config.nr_of_walkers; walk++) {
    walker::RandomWalk random_walk;

    std::function<void(const walker::RandomWalk &)> progress_monitor = [walk](const walker::RandomWalk &walk_in) {
      if (walk_in.sigma % 50 == 0)
        std::cout << "\t - w" << walk << " sigma = " << walk_in.sigma << std::endl;
    };

#pragma omp critical
    {
      if (config.continue_walk) {
        std::cout << "- Continuing Walk ";
        random_walk = walker::importRandomWalkFromBinaryFile(config.walk_dir + std::to_string(walk) + ".bin.dat");
        // random_walk = walker::importRandomWalkFromFile(walk_dirname + std::to_string(walk) + ".dat");
      } else {
        std::cout << "- Starting Walk ";
        random_walk = walker::setupRandomWalk(walker_start_nodes[walk], config.diffusion_constant);
      }

      std::cout << walk << " started at node " << random_walk.start_node << " and will walk for " << config.sigma_max
                << " steps with diffusion constant " << config.diffusion_constant << std::endl;
    }

    progressRandomWalk(G, random_walk, config.sigma_max, progress_monitor);

// ===  write to file  ==
// We drop the first and the last point because one cannot compute a derivative at these data point
#pragma omp critical
    {
      if (config.export_final_state) {
        std::cout << "- Exporting final state for walk " << walk << std::endl;
        if (!fs::is_directory(config.walk_dir) || !fs::exists(config.walk_dir)) {
          fs::create_directory(config.walk_dir);
        }
        std::string comment = "Random walk of graph " + config.graph;

        walker::exportRandomWalkToBinaryFile(random_walk, config.walk_dir + "/" + std::to_string(walk) + ".bin.dat");
      }

      for (int sigma = 1; sigma < random_walk.dimension.size() - 1; sigma++) {
        dimfile << random_walk.start_node << "\t" << sigma << "\t";
        dimfile << std::fixed << std::setprecision(12);
        dimfile << random_walk.dimension[sigma] << "\n";
      }
      std::cout << "- Writing walk " << walk << " to file" << std::endl;
    }
  }
  dimfile.close();
  std::cout << "- Finished walking..." << std::endl;
  return 0;
}

int main(int argc, char *argv[]) {
  WalkConfig config;
  // parse arguments
  int option;
  while ((option = getopt(argc, argv, "aecs:l:W:w:d:o:")) != -1) { // get option from the getopt() method
    switch (option) {
    case 'a':
      config.file_appending = true;
      break;
    case 'e':
      config.export_final_state = true;
      break;
    case 'c':
      config.continue_walk = true;
      config.export_final_state = true;
      break;
    case 's': {
      std::stringstream input(optarg);
      std::string int_str;
      while (std::getline(input, int_str, ',')) {
        config.start_node.push_back(std::stoi(int_str));
      }
    }
      config.random_start_nodes = false;
      break;
    case 'l':
      config.sigma_max = atoi(optarg);
      break;
    case 'W':
      config.walker_min_id = atoi(optarg);
      break;
    case 'w':
      config.nr_of_walkers = atoi(optarg);
      break;
    case 'd':
      config.diffusion_constant = atof(optarg);
      break;
    case 'o':
      config.dimension_file = optarg;
      break;
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

  if (config.dimension_file.empty()) {
    config.dimension_file = walker::DIMENSION_DIR + config.graph + ".dat";
    std::cout << "- no output file given - will write to " << config.dimension_file << std::endl;
  }

  if (!config.random_start_nodes && config.start_node.size() != config.nr_of_walkers) {
    std::cerr << "The number of start nodes does not match the number of walks you requested" << '\n';
    print_usage();
    return 1;
  }

  config.walk_dir = walker::WALK_DIR + config.graph;
  config.graph_file = walker::GRAPH_DIR + config.graph + ".dat";

  if (fs::exists(config.graph_file)) {
    std::cout << "- Loading graph file " << config.graph_file << std::endl;
    return process_network_or_graph<TUNGraph>(config);
  }

  config.graph_file = walker::NETWORK_DIR + config.graph + ".dat";
  if (fs::exists(config.graph_file)) {
    std::cout << "- Loading network file " << config.graph_file << std::endl;
    return process_network_or_graph<TIntNEDNet>(config);
  }

  std::cout << "- Neither graph or network found - please check that the graph " << config.graph << " exists!"
            << std::endl;
  return 1;
}
