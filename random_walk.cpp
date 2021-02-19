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
      << "usage: ./random_walk.x [-a] [-e] [-c] [-w walker_nr] [-W walker_min_id] [-s start_node] [-l length] [-d "
         "diffusion_constant] [-o "
         "outfile] graph_filename"
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
         "in data/dim_"
      << std::endl
      << "  graph_filename: The graph to run the random walker on" << std::endl;
}

int main(int argc, char *argv[]) {
  bool file_appending = false;
  bool export_final_state = false;
  bool continue_walk = false;
  int start_node = 0;
  int sigma_max = 100;
  int nr_of_walkers = 1;
  int walker_min_id = 0;
  double diffusion_constant = 0.5;
  std::string dimension_file;

  // parse arguments
  int option;
  while ((option = getopt(argc, argv, "aecs:l:W:w:d:o:")) != -1) { // get option from the getopt() method
    switch (option) {
    case 'a':
      file_appending = true;
      break;
    case 'e':
      export_final_state = true;
      break;
    case 'c':
      continue_walk = true;
      export_final_state = true;
      break;
    case 's':
      start_node = atoi(optarg);
      break;
    case 'l':
      sigma_max = atoi(optarg);
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

  std::cout << "  Input graph has " << G->GetNodes() << " nodes and " << G->GetEdges() << " edges";
  if (!TSnap::IsWeaklyConn<PUNGraph>(G)) {
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
  // in case walker_min_id > 0 we are generating more than the needed starting positions - the snap RNG is
  // deterministic, i.e. we are first reproducing the start positions of existing walks and then generating new starting
  // positions
  std::vector<int> walker_start_nodes;
  walker_start_nodes.reserve(walker_min_id + nr_of_walkers);
  if (nr_of_walkers > 1 && !continue_walk) {
    while (walker_start_nodes.size() < walker_min_id + nr_of_walkers) {
      int candidate = G->GetRndNId();
      // is that start node already in our pool of start nodes?
      if (std::find(walker_start_nodes.begin(), walker_start_nodes.end(), candidate) == walker_start_nodes.end())
        walker_start_nodes.push_back(candidate);
    }
  } else {
    walker_start_nodes[0] = start_node;
  }

#pragma omp parallel for
  for (int walk = walker_min_id; walk < walker_min_id + nr_of_walkers; walk++) {
    walker::RandomWalk random_walk;

    std::function<void(const walker::RandomWalk &)> progress_monitor = [walk](const walker::RandomWalk &walk_in) {
      if (walk_in.sigma % 50 == 0)
        std::cout << "\t - w" << walk << " sigma = " << walk_in.sigma << std::endl;
    };

#pragma omp critical
    {
      if (continue_walk) {
        std::cout << "- Continuing Walk ";
        random_walk = walker::importRandomWalkFromBinaryFile(walk_dirname + std::to_string(walk) + ".bin.dat");
        // random_walk = walker::importRandomWalkFromFile(walk_dirname + std::to_string(walk) + ".dat");
      } else {
        std::cout << "- Starting Walk ";
        random_walk = walker::setupRandomWalk(walker_start_nodes[walk], diffusion_constant);
      }

      std::cout << walk << " started at node " << random_walk.start_node << " and will walk for " << sigma_max
                << " steps with diffusion constant " << diffusion_constant << std::endl;
    }

    progressRandomWalk(G, random_walk, sigma_max, progress_monitor);

// ===  write to file  ==
// We drop the first and the last point because one cannot compute a derivative at these data point
#pragma omp critical
    {
      if (export_final_state) {
        std::cout << "- Exporting final state for walk " << walk << std::endl;
        if (!fs::is_directory(walk_dirname) || !fs::exists(walk_dirname)) {
          fs::create_directory(walk_dirname);
        }
        std::string comment = "Random walk of graph " + graph_filename;

        walker::exportRandomWalkToBinaryFile(random_walk, walk_dirname + std::to_string(walk) + ".bin.dat");
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
