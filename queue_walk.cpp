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
  std::cout
      << "usage: ./queue_walk.x [-w walker_nr] [-W walker_min_id] [-l length] [-d diffusion_constant] [-i interval] "
         "[-p] "
         "graph"
      << std::endl
      << "This walker is similar to random_walk.x but optimized for long running queue jobs" << std::endl
      << "  -w walker_nr: Number of walkers that will be run, defaults to 1" << std::endl
      << "  -W walker_min_id: the numbering of the walkers starts at this id - this allows to add random walks to an "
         "existing ensemble"
      << std::endl
      << "  -l maximum_length: Length of the random walk, defaults to 10000" << std::endl
      << "  -p probability distribution: Export the probability distribution p(r) at every interval" << std::endl
      << "  -L limit: target length for all walks - when this length is reached a walker is stopped" << std::endl
      << "  -d diffusion_constant delta: Diffusion constant to regularize osciallations in the spectral dimension, "
         "default: 0.5, (1-delta) is the probability to stay at the current node"
      << std::endl
      << "  -i interval: at multiples of interval the walker will write a snapshot to disk" << std::endl
      << "  graph: The graph to run the random walker on" << std::endl;
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

  std::ofstream radius_file;
  std::ofstream radius_nodes_file;
  if (config.export_prob_dist) {
    radius_file.open(walker::RADII_DIR + config.graph + ".dat", std::ios::out);
    radius_file << "# format: start_node sigma prob(0) prob(1) ..." << std::endl;
    radius_file.close();
    std::ofstream radius_nodes_file(walker::RADII_DIR + config.graph + "_nodes.dat", std::ios::out);
    radius_nodes_file << "# format: start_node nr_of_nodes_at(0) nr_of_nodes_at(1) ..." << std::endl;
    radius_nodes_file.close();
  }

#pragma omp parallel for
  for (int walk = config.walker_min_id; walk < config.walker_min_id + config.nr_of_walkers; walk++) {
    walker::RandomWalk random_walk;
    std::string walk_filename = walk_dirname + std::to_string(walk) + ".bin.dat";

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

    int max_radius = 0;
    TIntH node_radii(G->GetNodes());

    if (config.export_prob_dist) {
      // compute the shortest distance from startnode to every node
      TSnap::GetShortPath(G, random_walk.start_node, node_radii);
      for (auto it = node_radii.BegI(); it < node_radii.EndI(); it++) {

        if (it.GetDat() > max_radius)
          max_radius = it.GetDat();
      }
      std::vector<int> counts(max_radius + 1, 0);
      for (auto it = node_radii.BegI(); it < node_radii.EndI(); it++) {
        counts[it.GetDat()]++;
      }
      std::ofstream radius_file(walker::RADII_DIR + config.graph + "_nodes.dat", std::ios::out | std::ios::app);
      radius_file << random_walk.start_node;
      for (int count : counts) {
        radius_file << "\t" << count;
      }
      radius_file << std::endl;
      radius_file.close();
    }

    std::function<void(walker::RandomWalk)> progress_monitor = [walk, config, walk_filename, max_radius,
                                                                node_radii](walker::RandomWalk walk_in) {
      if (walk_in.sigma % 50 == 0)
        std::cout << "\t - w" << walk << " sigma = " << walk_in.sigma << std::endl;
      if (walk_in.sigma % config.interval == 0) {
        if (config.export_prob_dist) {
          std::vector<double> prob_of_radius(max_radius + 1, 0.0);
          for (auto it = walk_in.lvl_probabilities.BegI(); it < walk_in.lvl_probabilities.EndI(); it++) {
            const int radius = node_radii.GetDat(it.GetKey());
            prob_of_radius[radius] += it.GetDat();
          }
#pragma omp critical
          {
            std::ofstream radius_file(walker::RADII_DIR + config.graph + ".dat", std::ios::out | std::ios::app);
            radius_file << walk_in.start_node << "\t" << walk_in.sigma;
            for (double prob : prob_of_radius) {
              radius_file << "\t" << prob;
            }
            radius_file << std::endl;
            radius_file.close();
          }
        }
#pragma omp critical
        {
          std::cout << "- Exporting walk " << walk << std::endl;
          walker::exportRandomWalkToBinaryFile(walk_in, walk_filename);
        }
      }
    };

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
  while ((option = getopt(argc, argv, "l:L:W:w:d:i:p")) != -1) { // get option from the getopt() method
    switch (option) {
    case 'l':
      config.sigma_max = atoi(optarg);
      break;
    case 'L':
      config.limit_length = atoi(optarg);
      break;
    case 'p':
      config.export_prob_dist = true;
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
