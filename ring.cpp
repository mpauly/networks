#include "Snap.h"
#include "walker/base.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#define STARTNODE_REGULAR 0
#define WALKLENGTH_REGULAR 15000
#define WALKLENGTH_RANDOM 500
#define WALKLENGTH_RANDOM2 100000

using namespace std;

struct Parameter_pair {
  int nr_of_random_connections;
  double diffusion_constant;
};

int main(int argc, char *argv[]) {

  const std::function<void(int)> progress_monitor = [](int sigma) {
    // if (sigma % 50 == 0)
    //  std::cout << "." << std::flush;
  };
  // =========== without random connections ===========
  {
    typedef PUNGraph PGraph; // undirected graph

    PGraph G = PGraph::TObj::New();
    const int ring_length = 100;
    double diffusion_consts[7] = {0.05, 0.2, 0.35, 0.5, 0.65, 0.8, 0.95};

    std::ofstream dimfile;
    dimfile.open("data/dim_ring.dat", std::ofstream::out);
    dimfile << "# dimensions for ring with length " << ring_length << std::endl;
    dimfile << "# start_node: " << STARTNODE_REGULAR << std::endl;
    dimfile << "# walk length: " << WALKLENGTH_REGULAR << std::endl;
    dimfile << "# format: diffusion_const sigma d_spec" << std::endl;

    for (int n = 0; n < ring_length; n++) {
      G->AddNode(); // if no parameter is given, node ids are 0,1,...
    }
    for (int n = 0; n < ring_length; n++) {
      G->AddEdge(n, (n + 1) % ring_length);
    }

    for (const double &diffusion_constant : diffusion_consts) {
      std::cout << "== Making 1D regular ring with length " << ring_length << " and delta = " << diffusion_constant
                << " ==" << std::endl;

      auto walk_dimensions = walker::spectralDimensionAtNode(G, STARTNODE_REGULAR, WALKLENGTH_REGULAR, progress_monitor,
                                                             diffusion_constant);

      // We drop the first data point because one cannot compute a derivative at that data point
      for (int sigma = 1; sigma < walk_dimensions.size(); sigma++) {
        dimfile << diffusion_constant << "\t" << sigma << "\t";
        dimfile << std::fixed << std::setprecision(12);
        dimfile << walk_dimensions[sigma] << "\n";
      }
      std::cout << std::endl;
    }
    dimfile.close();
  }
  // =========== with random connections ==============
  {
    typedef PUNGraph PGraph; // undirected graph

    const int ring_length = 100;

    Parameter_pair parameter_pairs[12] = {{200, 0.05}, {200, 0.35}, {200, 0.5}, {200, 0.65}, {200, 0.8},  {200, 0.95},
                                          {200, 0.2},  {400, 0.2},  {600, 0.2}, {800, 0.2},  {1000, 0.2}, {2000, 0.2}};

    std::ofstream dimfile;
    dimfile.open("data/dim_ring_random.dat", std::ofstream::out);
    dimfile << "# average dimensions for ring with length " << ring_length << " and random connections" << std::endl;
    dimfile << "# averaging over all start_nodes " << std::endl;
    dimfile << "# walk length: " << WALKLENGTH_RANDOM << std::endl;
    dimfile << "# format: nr_of_random_conn diffusion_const sigma d_spec" << std::endl;

    for (const Parameter_pair &params : parameter_pairs) {
      std::cout << "== Making 1D ring with length " << ring_length << " and " << params.nr_of_random_connections
                << " random connections ==" << std::endl;

      PGraph G = PGraph::TObj::New();
      for (int n = 0; n < ring_length; n++) {
        G->AddNode(); // if no parameter is given, node ids are 0,1,...
      }

      for (int n = 0; n < ring_length; n++) {
        G->AddEdge(n, (n + 1) % ring_length);
      }

      // We are adding random edges to our graph
      int edges_added = 0;
      while (edges_added < params.nr_of_random_connections) {
        const int NId1 = G->GetRndNId();
        const int NId2 = G->GetRndNId();
        // this if condition checks if an edge already exists
        if (G->AddEdge(NId1, NId2) != -2) {
          edges_added++;
        }
      }

      // Average over all start positions
      std::vector<double> walk_dimensions(WALKLENGTH_RANDOM);
      for (int start = 0; start < ring_length; start++) {
        if (start % 10 == 0)
          std::cout << " start node = " << start << std::endl;
        std::vector<double> this_walk_dimensions =
            walker::spectralDimensionAtNode(G, start, WALKLENGTH_RANDOM, progress_monitor, params.diffusion_constant);
        for (int i = 0; i < WALKLENGTH_RANDOM; i++)
          walk_dimensions[i] += this_walk_dimensions[i];
      }
      for (int i = 0; i < WALKLENGTH_RANDOM; i++)
        walk_dimensions[i] = walk_dimensions[i] / ring_length;

      // We drop the first data point because one cannot compute a derivative at that data point
      for (int sigma = 1; sigma < walk_dimensions.size(); sigma++) {
        dimfile << params.nr_of_random_connections << "\t" << params.diffusion_constant << "\t" << sigma << "\t";
        dimfile << std::fixed << std::setprecision(12);
        dimfile << walk_dimensions[sigma] << "\n";
      }
      std::cout << std::endl;
    }
    dimfile.close();
  }
  // =========== with random connections and maximum random length link ==============
  {
    typedef PUNGraph PGraph; // undirected graph

    const int ring_length = 1000;
    const int max_length = 10;
    const double diffusion_constant = 0.5;

    std::vector<int> random_connections = {0, 50, 100, 200};

    std::ofstream dimfile;
    dimfile.open("data/dim_ring_random_maxlength.dat", std::ofstream::out);
    dimfile << "# average dimensions for ring with length " << ring_length
            << " and random connections of maximum length " << max_length << std::endl;
    dimfile << "# averaging over all start_nodes " << std::endl;
    dimfile << "# walk length: " << WALKLENGTH_RANDOM2 << std::endl;
    dimfile << "# diffusion const: " << diffusion_constant << std::endl;
    dimfile << "# format: nr_of_random_conn sigma d_spec" << std::endl;

#pragma omp parallel for
    for (int ran_con_count = 0; ran_con_count < random_connections.size(); ran_con_count++) {
      int nr_of_random_connections = random_connections.at(ran_con_count);
      std::cout << "== Making 1D ring with length " << ring_length << " and " << nr_of_random_connections
                << " random connections ==" << std::endl;

      PGraph G = PGraph::TObj::New();
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
        // check if this link is too long
        if (abs(NId1 - NId2) > max_length)
          continue;
        // this if condition checks if an edge already exists
        if (G->AddEdge(NId1, NId2) != -2) {
          edges_added++;
        }
      }

      // Average over all start positions
      std::vector<double> walk_dimensions(WALKLENGTH_RANDOM2);
      for (int start = 0; start < ring_length; start++) {
        if (start % 100 == 0)
          std::cout << " start node = " << start << std::endl;
        std::vector<double> this_walk_dimensions =
            walker::spectralDimensionAtNode(G, start, WALKLENGTH_RANDOM2, progress_monitor, diffusion_constant);
        for (int i = 0; i < WALKLENGTH_RANDOM2; i++)
          walk_dimensions[i] += this_walk_dimensions[i];
      }
      for (int i = 0; i < WALKLENGTH_RANDOM2; i++)
        walk_dimensions[i] = walk_dimensions[i] / ring_length;

// We drop the first data point because one cannot compute a derivative at that data point
#pragma omp critical
      {
        for (int sigma = 1; sigma < walk_dimensions.size(); sigma++) {
          dimfile << nr_of_random_connections << "\t" << sigma << "\t";
          dimfile << std::fixed << std::setprecision(12);
          dimfile << walk_dimensions[sigma] << "\n";
        }
        std::cout << std::endl;
      }
    }
    dimfile.close();
  }
  return 0;
}
