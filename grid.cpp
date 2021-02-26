#include "Snap.h"
#include "walker/base.h"
#include "walker/consts.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#define EDGE_LENGTH 100
#define STARTNODE 5050
#define WALKLENGTH 10000
#define INDEXAT2D(x, y) ((x)*EDGE_LENGTH + y)

using namespace std;

int main(int argc, char *argv[]) {

  typedef PUNGraph PGraph; // undirected graph

  PGraph G = PGraph::TObj::New();

  std::ofstream dimfile;
  dimfile.open(walker::DIMENSION_DIR + "grid.dat", std::ofstream::out);

  dimfile << "# dimensions for grid with length " << EDGE_LENGTH << std::endl;
  dimfile << "# start_node: " << STARTNODE << std::endl;
  dimfile << "# walk length: " << WALKLENGTH << std::endl;
  dimfile << "# format: nr_of_random_conn diffusion_const sigma d_spec" << std::endl;

  const std::function<void(walker::RandomWalk)> progress_monitor = [](walker::RandomWalk walk) {
    if (walk.sigma % 250 == 0)
      std::cout << "." << std::flush;
  };

  printf("== Making 2D grid with length %d ==\n", EDGE_LENGTH);
  for (int n = 0; n < EDGE_LENGTH * EDGE_LENGTH; n++) {
    G->AddNode(); // if no parameter is given, node ids are 0,1,...
  }
  // Build a square
  for (int n = 0; n < EDGE_LENGTH - 1; n++) {
    for (int m = 0; m < EDGE_LENGTH - 1; m++) {
      G->AddEdge(INDEXAT2D(n, m), INDEXAT2D(n + 1, m));
      G->AddEdge(INDEXAT2D(n, m), INDEXAT2D(n, m + 1));
    }
  }

  for (int j = 0; j < 7; j++) {
    double diffusion_constant = 0.05 + 0.15 * j;
    std::cout << "delta = " << diffusion_constant << " ";
    auto walk_dimensions =
        walker::spectralDimensionAtNode(G, STARTNODE, WALKLENGTH, progress_monitor, diffusion_constant);

    // We drop the first data point because one cannot compute a derivative at that data point
    for (int sigma = 1; sigma < walk_dimensions.size(); sigma++) {
      dimfile << diffusion_constant << "\t" << sigma << "\t";
      dimfile << std::fixed << std::setprecision(12);
      dimfile << walk_dimensions[sigma] << "\n";
    }
    std::cout << std::endl;
  }

  dimfile.close();

  return 0;
}
