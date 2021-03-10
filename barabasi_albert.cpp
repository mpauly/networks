#include "Snap.h"
#include "walker/base.h"
#include "walker/consts.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

void process_barabasi_albert(std::ofstream &dimfile, const std::vector<int> &degrees,
                             const std::vector<double> &diffusion_constants, const int &walk_length,
                             const int &nr_of_graphs, const int &nr_of_startnodes, const int &nr_of_nodes) {
  dimfile << "diffusion_const\tdegree\tsigma\tdim" << std::endl;

  const std::function<void(walker::RandomWalk)> progress_monitor = [](const walker::RandomWalk &walk) {
    // if (walk.sigma % 50 == 0)
    //  std::cout << "." << std::flush;
  };

  for (int degree : degrees) {
    for (double diffusion_constant : diffusion_constants) {
      std::cout << "== Degree " << degree << " diffusion_const " << diffusion_constant << " == " << std::endl;
      std::vector<double> dimension(walk_length, 0);
// we are averaging over nr_of_graphs different graphs
#pragma omp parallel for
      for (int graph = 0; graph < nr_of_graphs; graph++) {
        PUNGraph G;
        G = TSnap::GenPrefAttach(nr_of_nodes, degree);
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
#pragma omp critical
      {
        for (int i = 1; i < walk_length; i++) {
          dimfile << diffusion_constant << "\t" << degree << "\t" << i << "\t" << dimension[i] << "\n";
        }
      }
    }
  }
}

int main(int argc, char *argv[]) {

  {
    const int nr_of_nodes = 1000;
    const std::vector<int> degrees = {2, 3, 4, 5, 6};
    const int nr_of_graphs = 5;
    const int nr_of_startnodes = 25;
    const int walk_length = 200;
    const std::vector<double> diffusion_constants = {0.25, 0.5, 0.75};

    std::ofstream dimfile;
    dimfile.open(walker::DIMENSION_DIR + "barabasi_albert.dat", std::ofstream::out);

    process_barabasi_albert(dimfile, degrees, diffusion_constants, walk_length, nr_of_graphs, nr_of_startnodes,
                            nr_of_nodes);

    dimfile.close();
  }
  {
    const int nr_of_nodes = 100000;
    const std::vector<int> degrees = {2, 3, 4, 5, 6};
    const int nr_of_graphs = 5;
    const int nr_of_startnodes = 25;
    const int walk_length = 200;
    const std::vector<double> diffusion_constants = {0.25, 0.5, 0.75};

    std::ofstream dimfile;
    dimfile.open(walker::DIMENSION_DIR + "barabasi_albert_large.dat", std::ofstream::out);

    process_barabasi_albert(dimfile, degrees, diffusion_constants, walk_length, nr_of_graphs, nr_of_startnodes,
                            nr_of_nodes);

    dimfile.close();
  }
}