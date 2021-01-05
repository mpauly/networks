#include "Snap.h"
#include "walker.h"
#include <fstream>
#include <string>

#define STARTNODE 0
#define WALKLENGTH 1000

using namespace std;

int main(int argc, char *argv[]) {

  typedef PUNGraph PGraph; // undirected graph

  PGraph G = PGraph::TObj::New();
  const int ring_length = 100;

  std::ofstream dimfile;
  dimfile.open("data/dim_ring.dat", std::ofstream::out);

  dimfile << "# dimensions for ring with length " << ring_length << std::endl;
  dimfile << "# start_node: " << STARTNODE << std::endl;
  dimfile << "# walk length: " << WALKLENGTH << std::endl;
  dimfile << "# format: nr_of_random_conn diffusion_const sigma d_spec" << std::endl;

  std::function<void(int)> progress_monitor = [](int sigma) {
    if (sigma % 50 == 0)
      std::cout << "." << std::flush;
  };

  for (int nr_of_random_connections = 0; nr_of_random_connections <= 1000; nr_of_random_connections += 200) {
    std::cout << "== Making 1D ring with length " << ring_length << " and " << nr_of_random_connections
              << " random connections ==" << std::endl;
    for (int n = 0; n < ring_length; n++) {
      G->AddNode(); // if no parameter is given, node ids are 0,1,...
    }

    for (int n = 0; n < ring_length; n++) {
      G->AddEdge(n, (n + 1) % ring_length);
    }

    // We are adding random edges to our graph
    int edges_added = 0;
    while (edges_added < nr_of_random_connections) {
      const int NId1 = G->GetRndNId();
      const int NId2 = G->GetRndNId();
      // this if condition checks if an edge already exists
      if (G->AddEdge(NId1, NId2) != -2) {
        edges_added++;
      }
    }

    for (int j = 0; j < 5; j++) {
      double diffusion_constant = 0.05 + 0.15 * j;
      std::cout << "delta = " << diffusion_constant << " ";
      auto walk_dimensions = spectralDimensionAtNode(G, STARTNODE, WALKLENGTH, progress_monitor, diffusion_constant);

      // We drop the first data point because one cannot compute a derivative at that data point
      for (int sigma = 1; sigma < walk_dimensions.size(); sigma++) {
        dimfile << nr_of_random_connections << "\t" << diffusion_constant << "\t" << sigma << "\t";
        dimfile << std::fixed << std::setprecision(12);
        dimfile << walk_dimensions[sigma] << "\n";
      }
      std::cout << std::endl;
    }
  }
  dimfile.close();

  return 0;
}
