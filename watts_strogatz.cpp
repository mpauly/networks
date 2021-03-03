#include "Snap.h"
#include "walker/base.h"
#include "walker/consts.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

int main(int argc, char *argv[]) {

  const int nr_of_nodes = 1000;
  const std::vector<int> degrees = {2, 4, 6};
  const int nr_of_graphs = 5;
  const int nr_of_startnodes = 25;
  const int walk_length = 200;
  const std::vector<double> diffusion_constants = {0.25, 0.5, 0.75};
  const std::vector<double> rewiring_probs = {0.0, 0.05, 0.1, 0.15, 0.2, 0.4, 0.6, 0.8};

  std::ofstream dimfile;
  dimfile.open(walker::DIMENSION_DIR + "watts_strogatz.dat", std::ofstream::out);
  dimfile << "rewiring_prob\tdiffusion_const\tdegree\tsigma\tdim" << std::endl;

  const std::function<void(walker::RandomWalk)> progress_monitor = [](const walker::RandomWalk &walk) {
    // if (walk.sigma % 50 == 0)
    //  std::cout << "." << std::flush;
  };

  for (int degree : degrees) {
    for (double diffusion_constant : diffusion_constants) {
      for (double rewiring : rewiring_probs) {
        std::cout << "== Rewiring prob. " << rewiring << " == " << std::endl;
        std::vector<double> dimension(walk_length, 0);
        // we are averaging over nr_of_graphs different graphs
        for (int graph = 0; graph < nr_of_graphs; graph++) {
          PUNGraph G;
          G = TSnap::GenSmallWorld(nr_of_nodes, degree, rewiring);
          std::vector<int> start_nodes = {};
          start_nodes.reserve(nr_of_startnodes);
          while (start_nodes.size() < nr_of_startnodes) {
            const int candidate = G->GetRndNId();
            if (std::find(start_nodes.begin(), start_nodes.end(), candidate) == start_nodes.end())
              start_nodes.emplace_back(candidate);
          }

          // ... and nr_of_startnodes different starting positions
          for (int startnode : start_nodes) {
            std::vector<double> this_walk_dimensions =
                walker::spectralDimensionAtNode(G, startnode, walk_length, progress_monitor, diffusion_constant);
            for (int i = 0; i < walk_length; i++) {
              dimension[i] += this_walk_dimensions[i] / nr_of_startnodes / nr_of_graphs;
            }
          }
        }
        for (int i = 1; i < walk_length; i++) {
          dimfile << rewiring << "\t" << diffusion_constant << "\t" << degree << "\t" << i << "\t" << dimension[i]
                  << "\n";
        }
      }
    }
  }
}