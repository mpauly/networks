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
  double diffusion_constant = 0.5;
  std::string graph;
};

void print_usage() {
  std::cout
      << "usage: ./queue_walk.x [-w walker_nr] [-W walker_min_id] [-l length] [-d diffusion_constant] [-i interval] "
         "graph"
      << std::endl
      << "This walker is similar to random_walk.x but optimized for long running queue jobs" << std::endl
      << "  -w walker_nr: Number of walkers that will be run, defaults to 1" << std::endl
      << "  -W walker_min_id: the numbering of the walkers starts at this id - this allows to add random walks to an "
         "existing ensemble"
      << std::endl
      << "  -l maximum_length: Length of the random walk, defaults to 10000" << std::endl
      << "  -L limit: target length for all walks - when this length is reached a walker is stopped" << std::endl
      << "  -d diffusion_constant delta: Diffusion constant to regularize osciallations in the spectral dimension, "
         "default: 0.5, (1-delta) is the probability to stay at the current node"
      << std::endl
      << "  -i interval: at multiples of interval the walker will write a snapshot to disk" << std::endl
      << "  graph: The graph to run the random walker on" << std::endl;
}

template <class Graph> int process_network_or_graph(WalkConfig config) {
  std::string walk_dirname = walker::WALK_DIR + config.graph + '/';
  std::string graph_filename = walker::GRAPH_DIR + config.graph + ".dat";
  TPt<Graph> G;
  try {
    TFIn FIn(graph_filename.c_str());
    G = Graph::Load(FIn);
  } catch (...) {
    std::cerr << "Error reading file " << graph_filename << std::endl;
    return 1;
  }

  std::cout << "  Input graph has " << G->GetNodes() << " nodes and " << G->GetEdges() << " edges" << std::endl;

  // Determine starting positions - this is outside of the parallel section to avoid complications regarding the
  // generation of random numbers in different threads
  // in case walker_min_id > 0 we are generating more than the needed starting positions - the snap RNG is
  // deterministic, i.e. we are first reproducing the start positions of existing walks and then generating new starting
  // positions
  std::vector<int> walker_start_nodes;
  walker_start_nodes.reserve(config.walker_min_id + config.nr_of_walkers);
  while (walker_start_nodes.size() < config.walker_min_id + config.nr_of_walkers) {
    int candidate = G->GetRndNId();
    // is that start node already in our pool of start nodes?
    if (std::find(walker_start_nodes.begin(), walker_start_nodes.end(), candidate) == walker_start_nodes.end())
      walker_start_nodes.push_back(candidate);
  }

#pragma omp parallel for
  for (int walk = config.walker_min_id; walk < config.walker_min_id + config.nr_of_walkers; walk++) {
    walker::RandomWalk random_walk;
    std::string walk_filename = walk_dirname + std::to_string(walk) + ".bin.dat";

    std::function<void(walker::RandomWalk)> progress_monitor = [walk, config,
                                                                walk_filename](walker::RandomWalk walk_in) {
      if (walk_in.sigma % 50 == 0)
        std::cout << "\t - w" << walk << " sigma = " << walk_in.sigma << std::endl;
      if (walk_in.sigma % config.interval == 0) {
#pragma omp critical
        {
          std::cout << "- Exporting walk " << walk << std::endl;
          walker::exportRandomWalkToBinaryFile(walk_in, walk_filename);
        }
      }
    };

    int steps = config.sigma_max;

#pragma omp critical
    {
      if (fs::exists(walk_filename)) {
        std::cout << "- Continuing Walk ";
        random_walk = walker::importRandomWalkFromBinaryFile(walk_filename);
      } else {
        std::cout << "- Starting Walk ";
        random_walk = walker::setupRandomWalk(walker_start_nodes[walk], config.diffusion_constant);
      }

      if (config.limit_length > 0 && random_walk.sigma + steps > config.limit_length) {
        steps = config.limit_length - random_walk.sigma;
      }

      std::cout << walk << " at node " << random_walk.start_node << " and will walk for " << steps
                << " steps with diffusion constant " << config.diffusion_constant << std::endl;
    }

    progressRandomWalk(G, random_walk, steps, progress_monitor);

    if (random_walk.sigma % config.interval != 0) {
#pragma omp critical
      {

        std::cout << "- Exporting final state for walk " << walk << std::endl;
        if (!fs::is_directory(walk_dirname) || !fs::exists(walk_dirname)) {
          fs::create_directory(walk_dirname);
        }

        walker::exportRandomWalkToBinaryFile(random_walk, walk_filename);
      }
    }
    std::cout << "- Walk " << walk << " finished" << std::endl;
  }
  std::cout << "- Finished walking..." << std::endl;
  return 0;
}

int main(int argc, char *argv[]) {
  WalkConfig config;
  // parse arguments
  int option;
  while ((option = getopt(argc, argv, "l:L:W:w:d:i:")) != -1) { // get option from the getopt() method
    switch (option) {
    case 'l':
      config.sigma_max = atoi(optarg);
      break;
    case 'L':
      config.limit_length = atoi(optarg);
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
    case 'i':
      config.interval = atoi(optarg);
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

  std::string graph_filename = walker::GRAPH_DIR + config.graph + ".dat";
  if (fs::exists(graph_filename)) {
    std::cout << "- Loading graph file " << graph_filename << std::endl;
    return process_network_or_graph<TUNGraph>(config);
  }

  graph_filename = walker::NETWORK_DIR + config.graph + ".dat";
  if (fs::exists(graph_filename)) {
    std::cout << "- Loading network file " << graph_filename << std::endl;
    return process_network_or_graph<TIntNEDNet>(config);
  }

  std::cout << "- Neither graph or network found - please check that the graph " << config.graph << " exists!"
            << std::endl;
  return 1;
}
