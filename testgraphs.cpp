#include "stdafx.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#define EDGE_LENGTH 100
#define INDEXAT2D(x, y) ((x)*EDGE_LENGTH + y)
#define INDEXAT3D(x, y, z) ((x)*EDGE_LENGTH * EDGE_LENGTH + (y)*EDGE_LENGTH + z)

// This method just computes the probability weight per node - here we assume a uniform distribution
double get_node_weight(int outDegree) { return 1.0 / outDegree; }

// For a given startNode this method returns a hashmap of vectors, such that the keys of the hashmap
// corresond to the Ids of the nodes that one wants to study and the vectors correspond to the probability
// probability_trail for one node. The implementation probably is not super efficent, complexity grows exponentially.
// The implementation roughly corresponds to a breadth-first search tree.
// Inspired by DoBfs
template <class PGraph>
THash<TInt, std::vector<double>> traverseBfsTree(const PGraph &Graph, const int &start_node, const int &max_depth,
                                                 std::ofstream &dimfile) {
  // setup data structures
  TSnapQueue<int> queue;
  std::vector<THash<TInt, TInt>> levelCounts(max_depth + 1);
  std::vector<int> totalCounts(max_depth + 1);
  THash<TInt, std::vector<double>> levelProbabilities;
  THash<TInt, std::vector<double>> nodeDimensionality;
  THash<TInt, int> nodeFirstVisited;
  std::vector<double> thisLvlProbs;
  IAssert(Graph->IsNode(start_node));
  int visits;

  // setup the queue
  queue.Clr(false);
  queue.Push(start_node);
  int v, MaxDist = 0;
  levelCounts[0].AddDat(start_node, 1);
  std::vector<double> zero_probs(max_depth + 1);
  zero_probs[0] = 1.0;
  levelProbabilities.AddDat(start_node, zero_probs);
  nodeFirstVisited.AddDat(start_node, 0);
  double node_probability_weight, probability_derivative, mean_probability;

  // first we keep around a vector of length max_depth for each of the sites - in a second step we then trim this vector
  for (int lvl = 1; lvl <= max_depth; lvl++) {
    // deal with all the nodes that are in the queue at the moment
    // this relies on the loop initialization begin run exactly once, as the queue grows while the loop is executed
    for (int node = queue.Len(); node > 0; node--) {
      // take a node out of the queue
      const int nodeId = queue.Top();
      queue.Pop();
      const typename PGraph::TObj::TNodeI NodeI = Graph->GetNI(nodeId);
      // loop over child nodes
      const int noteDegree = NodeI.GetOutDeg();
      for (v = 0; v < noteDegree; v++) {
        const int childId = NodeI.GetOutNId(v);
        // compute the probability weight for this node
        node_probability_weight = get_node_weight(noteDegree);
        // did we already touch this node? If yes just update existing entry
        // if no then create a new entry and add the node to the queue

        if (levelCounts[lvl].IsKey(childId)) {
          visits = levelCounts[lvl].GetDat(childId) + 1;
        } else {
          visits = 1;
          queue.Push(childId);
        }

        if (levelProbabilities.IsKey(childId)) {
          thisLvlProbs = levelProbabilities.GetDat(childId);
          thisLvlProbs[lvl] = thisLvlProbs[lvl] + levelProbabilities.GetDat(nodeId)[lvl - 1] * node_probability_weight;
        } else {
          nodeFirstVisited.AddDat(childId, lvl);
          thisLvlProbs = std::vector<double>(max_depth + 1);
          thisLvlProbs[lvl] = levelProbabilities.GetDat(nodeId)[lvl - 1] * node_probability_weight;
        }
        // AddDat also overwrites an existing entry
        levelCounts[lvl].AddDat(childId, visits);
        levelProbabilities.AddDat(childId, thisLvlProbs);
      }
    }
  }

  // Extract the dimensionality

  int nodeId, startLevel, sigma;
  for (int lvl = 0; lvl < max_depth; lvl++) {
    sigma = lvl;

    // we are using a symmetric definition of the first derivative here and additional using the mean of
    // neighbouring sites for the probability
    mean_probability = 0.5 * (probability_trail[lvl - 1] + probability_trail[lvl + 1]);
    probability_derivative = 0.5 * (probability_trail[lvl + 1] - probability_trail[lvl - 1]);
    // compute the spectral dimension
    dimensionality[sigma] = -2.0 * sigma / mean_probability * probability_derivative;
    // printf("%f \t", dimensionality[sigma]);
    // write to file
    dimfile << std::fixed << std::setprecision(12);
    dimfile << dimensionality[sigma] << "\t";
  }
  dimfile << "\n";

  return levelProbabilities;
}

// ------------------ the main function -----------------------
int main(int argc, char *argv[]) {
  typedef PUNGraph PGraph; // undirected graph

  // Create the graph
  PGraph G = PGraph::TObj::New();
  int ring_length = 1000;
  for (int n = 0; n < ring_length; n++) {
    G->AddNode(); // if no parameter is given, node ids are 0,1,...
  }

  // Build a ring
  for (int n = 0; n < ring_length; n++) {
    G->AddEdge(n, (n + 1) % ring_length);
  }

  // Adding random edges
  // for (int e = 0; e < 10; e++) {
  //  const int NId1 = G->GetRndNId();
  //  const int NId2 = G->GetRndNId();
  //  if (G->AddEdge(NId1, NId2) != -2) {
  //    printf("  Edge %d -- %d added\n", NId1, NId2);
  //  } else {
  //    printf("  Edge %d -- %d already exists\n", NId1, NId2);
  //  }
  //}

  TSnap::DrawGViz(G, gvlDot, "graph.png");

  std::ofstream dimfile;
  dimfile.open("dimension.dat");
  dimfile << "# dimensions for a chain of length \n";
  // arguments here are start_node and depth
  auto treeCounts = traverseBfsTree(G, 7, 100, dimfile);
  dimfile.close();

  // ================== 2D lattice ===================================
  // Create the graph
  PGraph G2 = PGraph::TObj::New();
  for (int n = 0; n < EDGE_LENGTH * EDGE_LENGTH; n++) {
    G2->AddNode(); // if no parameter is given, node ids are 0,1,...
  }

  // Build a ring
  for (int n = 0; n < EDGE_LENGTH - 1; n++) {
    for (int m = 0; m < EDGE_LENGTH - 1; m++) {
      if (n == 50 && m == 50) {
        printf("  Edge %d -- %d added\n", INDEXAT2D(n, m), INDEXAT2D(n + 1, m));
        printf("  Edge %d -- %d added\n", INDEXAT2D(n, m), INDEXAT2D(n, m + 1));
      }
      G2->AddEdge(INDEXAT2D(n, m), INDEXAT2D(n + 1, m));
      G2->AddEdge(INDEXAT2D(n, m), INDEXAT2D(n, m + 1));
    }
  }

  const typename PGraph::TObj::TNodeI NodeI = G2->GetNI(INDEXAT2D(50, 50));
  // loop over child nodes
  const int noteDegree = NodeI.GetOutDeg();
  printf("  Start node %d has degree %d\n", INDEXAT2D(50, 50), noteDegree);
  for (int v = 0; v < noteDegree; v++) {
    printf("\t %d", NodeI.GetOutNId(v));
  }
  printf("Should be %d \t %d \t %d \t %d", INDEXAT2D(49, 50), INDEXAT2D(51, 50), INDEXAT2D(50, 51), INDEXAT2D(50, 49));
  printf("\n");
  // Adding random edges
  // for (int e = 0; e < 10; e++) {
  //  const int NId1 = G->GetRndNId();
  //  const int NId2 = G->GetRndNId();
  //  if (G->AddEdge(NId1, NId2) != -2) {
  //    printf("  Edge %d -- %d added\n", NId1, NId2);
  //  } else {
  //    printf("  Edge %d -- %d already exists\n", NId1, NId2);
  //  }
  //}

  std::ofstream dimfile2;
  dimfile2.open("dimension_2D.dat");
  dimfile2 << "# dimensions for a square \n";
  // arguments here are start_node and depth
  printf("Starting Walker \n");
  auto treeCounts2 = traverseBfsTree(G2, INDEXAT2D(50, 50), 50, dimfile2);
  dimfile2.close();

  // ================== 3D lattice ===================================
  // Create the graph
  PGraph G3 = PGraph::TObj::New();
  for (int n = 0; n < EDGE_LENGTH * EDGE_LENGTH * EDGE_LENGTH; n++) {
    G3->AddNode(); // if no parameter is given, node ids are 0,1,...
  }

  // Build a ring
  for (int n1 = 0; n1 < EDGE_LENGTH - 1; n1++) {
    for (int n2 = 0; n2 < EDGE_LENGTH - 1; n2++) {
      for (int n3 = 0; n3 < EDGE_LENGTH - 1; n3++) {
        G3->AddEdge(INDEXAT3D(n1, n2, n3), INDEXAT3D(n1 + 1, n2, n3));
        G3->AddEdge(INDEXAT3D(n1, n2, n3), INDEXAT3D(n1, n2 + 1, n3));
        G3->AddEdge(INDEXAT3D(n1, n2, n3), INDEXAT3D(n1, n2, n3 + 1));
      }
    }
  }

  // Adding random edges
  // for (int e = 0; e < 10; e++) {
  //  const int NId1 = G->GetRndNId();
  //  const int NId2 = G->GetRndNId();
  //  if (G->AddEdge(NId1, NId2) != -2) {
  //    printf("  Edge %d -- %d added\n", NId1, NId2);
  //  } else {
  //    printf("  Edge %d -- %d already exists\n", NId1, NId2);
  //  }
  //}

  std::ofstream dimfile3;
  dimfile3.open("dimension_3D.dat");
  dimfile3 << "# dimensions for a cube \n";
  // arguments here are start_node and depth
  printf("Starting Walker \n");
  auto treeCounts3 = traverseBfsTree(G3, INDEXAT3D(50, 50, 50), 50, dimfile3);
  dimfile3.close();

  // the following are approximate neighbourhood functions - might be
  // interesting at some point printf("== TestAnf ==\n");
  // TSnap::TestAnf<PUNGraph>();
  return 0;
}
