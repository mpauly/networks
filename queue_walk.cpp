#include "Snap.h"
#include "walker/base.h"
#include "walker/io.h"
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <unistd.h>

namespace fs = std::filesystem;

void print_usage() {
  std::cout
      << "usage: ./queue_walk.x [-w walker_nr] [-W walker_min_id] [-l length] [-d diffusion_constant] [-i interval] "
         "graph_filename"
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
      << "  graph_filename: The graph to run the random walker on" << std::endl;
}

int main(int argc, char *argv[]) {
  int sigma_max = 10000;
  int nr_of_walkers = 1;
  int walker_min_id = 0;
  int interval = 500;
  int limit_length = -1;
  double diffusion_constant = 0.5;

  // parse arguments
  int option;
  while ((option = getopt(argc, argv, "l:L:W:w:d:i:")) != -1) { // get option from the getopt() method
    switch (option) {
    case 'l':
      sigma_max = atoi(optarg);
      break;
    case 'L':
      limit_length = atoi(optarg);
      break;
    case 'W':
      walker_min_id = atoi(optarg);
      break;
    case 'w':
      nr_of_walkers = atoi(optarg);
      break;
    case 'd':
      diffusion_constant = atof(optarg);
      break;
    case 'i':
      interval = atoi(optarg);
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
    std::cout << "Please specify exactly one graph" << std::endl;
    print_usage();
    return 1;
  }

  std::string walk_dirname = graph_filename;
  walk_dirname = walk_dirname.replace(walk_dirname.find("graphs/"), sizeof("graphs/") - 1, "walks/");
  walk_dirname = walk_dirname.replace(walk_dirname.find(".dat"), sizeof(".dat") - 1, "/");

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

  std::cout << "  Input graph has " << G->GetNodes() << " nodes and " << G->GetEdges() << " edges" << std::endl;

  // Determine starting positions - this is outside of the parallel section to avoid complications regarding the
  // generation of random numbers in different threads
  // in case walker_min_id > 0 we are generating more than the needed starting positions - the snap RNG is
  // deterministic, i.e. we are first reproducing the start positions of existing walks and then generating new starting
  // positions
  std::vector<int> walker_start_nodes;
  walker_start_nodes.reserve(walker_min_id + nr_of_walkers);
  while (walker_start_nodes.size() < walker_min_id + nr_of_walkers) {
    int candidate = G->GetRndNId();
    // is that start node already in our pool of start nodes?
    if (std::find(walker_start_nodes.begin(), walker_start_nodes.end(), candidate) == walker_start_nodes.end())
      walker_start_nodes.push_back(candidate);
  }

#pragma omp parallel for
  for (int walk = walker_min_id; walk < walker_min_id + nr_of_walkers; walk++) {
    walker::RandomWalk random_walk;
    std::string walk_filename = walk_dirname + std::to_string(walk) + ".bin.dat";

    std::function<void(walker::RandomWalk)> progress_monitor = [walk, interval,
                                                                walk_filename](walker::RandomWalk walk_in) {
      if (walk_in.sigma % 50 == 0)
        std::cout << "\t - w" << walk << " sigma = " << walk_in.sigma << std::endl;
      if (walk_in.sigma % interval == 0) {
#pragma omp critical
        {
          std::cout << "- Exporting walk " << walk << std::endl;
          walker::exportRandomWalkToBinaryFile(walk_in, walk_filename);
        }
      }
    };

    int steps = sigma_max;

#pragma omp critical
    {
      if (fs::exists(walk_filename)) {
        std::cout << "- Continuing Walk ";
        random_walk = walker::importRandomWalkFromBinaryFile(walk_filename);
      } else {
        std::cout << "- Starting Walk ";
        random_walk = walker::setupRandomWalk(walker_start_nodes[walk], diffusion_constant);
      }

      if (limit_length > 0 && random_walk.sigma + steps > limit_length) {
        steps = limit_length - random_walk.sigma;
      }

      std::cout << walk << " at node " << random_walk.start_node << " and will walk for " << steps
                << " steps with diffusion constant " << diffusion_constant << std::endl;
    }

    progressRandomWalk(G, random_walk, steps, progress_monitor);

    if (random_walk.sigma % interval != 0) {
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
