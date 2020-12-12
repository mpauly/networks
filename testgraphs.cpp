#include "Snap.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#define EDGE_LENGTH 100
#define INDEXAT2D(x, y) ((x)*EDGE_LENGTH + y)
#define INDEXAT3D(x, y, z) ((x)*EDGE_LENGTH * EDGE_LENGTH + (y)*EDGE_LENGTH + z)
#define DEBUG false

// For a given startNode this method returns a hashmap of vectors, such that the keys of the hashmap
// corresond to the Ids of the nodes that one wants to study and the vectors correspond to the probability
// probability_trail for one node. The implementation probably is not super efficent, complexity grows exponentially.
// The implementation roughly corresponds to a breadth-first search tree.
// Inspired by DoBfs
template <class PGraph>
std::vector<double> spectralDimensionAtNode(const PGraph &Graph, const int &start_node, const int &max_depth,
                                            std::ofstream &dimfile) {
  // setup data structures
  TSnapQueue<int> queue;
  std::vector<THash<TInt, double>> levelProbabilities(max_depth + 1);
  // std::vector<long> totalCounts(max_depth + 1);
  std::vector<double> return_probability(max_depth + 1);
  std::vector<double> dimension(max_depth);

  double probability_derivative, mean_probability;
  double node_probability;

  IAssert(Graph->IsNode(start_node));

  // setup the queue
  queue.Clr(false);
  queue.Push(start_node);
  int v, MaxDist = 0;
  // Set the level zero quantities by hand before we start walking
  levelProbabilities[0].AddDat(start_node, 1);

  for (int sigma = 1; sigma <= max_depth; sigma++) {
    if (DEBUG)
      printf("\n Sigma %d\n", sigma);
    // deal with all the nodes that are in the queue at the moment
    // this relies on the loop initialization begin run exactly once, as the queue grows while the loop is executed
    for (int node = queue.Len(); node > 0; node--) {
      // take a node out of the queue
      const int nodeId = queue.Top();
      if (DEBUG)
        printf("   Node %d   ", nodeId);
      queue.Pop();
      const typename PGraph::TObj::TNodeI NodeI = Graph->GetNI(nodeId);
      // loop over child nodes
      const int nodeDegree = NodeI.GetOutDeg();
      for (v = 0; v < nodeDegree; v++) {
        const int childId = NodeI.GetOutNId(v);
        // did we already touch this node? If yes just update existing entry
        // if no then create a new entry and add the node to the queue
        if (levelProbabilities[sigma].IsKey(childId)) {
          node_probability = levelProbabilities[sigma].GetDat(childId) +
                             levelProbabilities[sigma - 1].GetDat(nodeId) / ((double)nodeDegree);
        } else {
          node_probability = levelProbabilities[sigma - 1].GetDat(nodeId) / ((double)nodeDegree);
          queue.Push(childId);
        }
        if (DEBUG)
          printf(" %d : %f ", childId, node_probability);
        // AddDat also overwrites an existing entry
        levelProbabilities[sigma].AddDat(childId, node_probability);
      }
    }
  }
  // Compute the return probability
  for (int sigma = 0; sigma <= max_depth; sigma++) {
    if (levelProbabilities[sigma].IsKey(start_node)) {
      return_probability[sigma] = levelProbabilities[sigma].GetDat(start_node);
    } else {
      return_probability[sigma] = 0.0;
    }
  }

  // Extract the dimensionality
  int nodeId, startLevel;
  for (int sigma = 0; sigma < max_depth; sigma++) {
    // we are using a symmetric definition of the first derivative here and additional using the mean of
    // neighbouring sites for the probability
    mean_probability = 0.5 * (return_probability[sigma - 1] + return_probability[sigma + 1]);
    probability_derivative = 0.5 * (return_probability[sigma + 1] - return_probability[sigma - 1]);
    // compute the spectral dimension
    dimension[sigma] = -2.0 * ((double)sigma) / mean_probability * probability_derivative;
    if (DEBUG) {
      printf("Dimension: \n");
      printf("%d %f \n", sigma, dimension[sigma]);
    }
    // write to file
    dimfile << sigma << "\t";
    dimfile << std::fixed << std::setprecision(12);
    dimfile << dimension[sigma] << "\n";
  }

  return dimension;
}

// ------------------ the main function -----------------------
int main(int argc, char *argv[]) {
  typedef PUNGraph PGraph; // undirected graph

  // =============== 1D ring ===============================
  PGraph G = PGraph::TObj::New();
  int ring_length = 26;
  for (int n = 0; n < ring_length; n++) {
    G->AddNode(); // if no parameter is given, node ids are 0,1,...
  }

  for (int n = 0; n < ring_length; n++) {
    G->AddEdge(n, (n + 1) % ring_length);
  }

  std::ofstream dimfile;
  dimfile.open("data/dimension_1d_short.dat");
  dimfile << "# dimensions for a ring of length " << ring_length << "\n";
  auto treeCounts = spectralDimensionAtNode(G, 7, 150, dimfile);
  dimfile.close();

  // =============== 1D ring ===============================
  G = PGraph::TObj::New();
  ring_length = 100;
  for (int n = 0; n < ring_length; n++) {
    G->AddNode(); // if no parameter is given, node ids are 0,1,...
  }

  for (int n = 0; n < ring_length; n++) {
    G->AddEdge(n, (n + 1) % ring_length);
  }

  dimfile.open("data/dimension_1d.dat");
  dimfile << "# dimensions for a ring of length " << ring_length << "\n";
  treeCounts = spectralDimensionAtNode(G, 7, 150, dimfile);
  dimfile.close();

  // =============== 1D ring with random connections ===============================
  PGraph G1a = PGraph::TObj::New();
  ring_length = 100;
  for (int n = 0; n < ring_length; n++) {
    G1a->AddNode(); // if no parameter is given, node ids are 0,1,...
  }

  for (int n = 0; n < ring_length; n++) {
    G1a->AddEdge(n, (n + 1) % ring_length);
  }

  const int nr_of_random_connections = 50;

  // Adding random edges
  int edges_added = 0;
  while (edges_added < nr_of_random_connections) {
    const int NId1 = G->GetRndNId();
    const int NId2 = G->GetRndNId();
    // this if condition checks if an edge already exists
    if (G->AddEdge(NId1, NId2) != -2) {
      edges_added++;
    }
  }

  TSnap::DrawGViz(G, gvlDot, "graph.png");

  std::ofstream dimfile1a;
  dimfile1a.open("data/dimension_1d_connected.dat");
  dimfile1a << "# dimensions for a ring with " << nr_of_random_connections
            << " additional random connections of length " << ring_length << "\n";
  auto treeCounts1a = spectralDimensionAtNode(G, 7, 150, dimfile1a);
  dimfile1a.close();

  // ================== 2D lattice ===================================
  PGraph G2 = PGraph::TObj::New();
  for (int n = 0; n < EDGE_LENGTH * EDGE_LENGTH; n++) {
    G2->AddNode(); // if no parameter is given, node ids are 0,1,...
  }
  // Build a square
  for (int n = 0; n < EDGE_LENGTH - 1; n++) {
    for (int m = 0; m < EDGE_LENGTH - 1; m++) {
      G2->AddEdge(INDEXAT2D(n, m), INDEXAT2D(n + 1, m));
      G2->AddEdge(INDEXAT2D(n, m), INDEXAT2D(n, m + 1));
    }
  }
  std::ofstream dimfile2;
  dimfile2.open("data/dimension_2d.dat");
  dimfile2 << "# dimensions for a square \n";
  auto treeCounts2 = spectralDimensionAtNode(G2, INDEXAT2D(50, 50), 100, dimfile2);
  dimfile2.close();

  // ================== 3D lattice ===================================
  PGraph G3 = PGraph::TObj::New();
  for (int n = 0; n < EDGE_LENGTH * EDGE_LENGTH * EDGE_LENGTH; n++) {
    G3->AddNode(); // if no parameter is given, node ids are 0,1,...
  }

  // Build a cube
  for (int n1 = 0; n1 < EDGE_LENGTH - 1; n1++) {
    for (int n2 = 0; n2 < EDGE_LENGTH - 1; n2++) {
      for (int n3 = 0; n3 < EDGE_LENGTH - 1; n3++) {
        G3->AddEdge(INDEXAT3D(n1, n2, n3), INDEXAT3D(n1 + 1, n2, n3));
        G3->AddEdge(INDEXAT3D(n1, n2, n3), INDEXAT3D(n1, n2 + 1, n3));
        G3->AddEdge(INDEXAT3D(n1, n2, n3), INDEXAT3D(n1, n2, n3 + 1));
      }
    }
  }

  std::ofstream dimfile3;
  dimfile3.open("data/dimension_3d.dat");
  dimfile3 << "# dimensions for a cube \n";
  auto treeCounts3 = spectralDimensionAtNode(G3, INDEXAT3D(50, 50, 50), 100, dimfile3);
  dimfile3.close();

  // the following are approximate neighbourhood functions - might be
  // interesting at some point printf("== TestAnf ==\n");
  // TSnap::TestAnf<PUNGraph>();
  return 0;
}
