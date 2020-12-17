#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>

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
  std::vector<THash<TInt, double>> levelProbabilities(max_depth + 2);
  // std::vector<long> totalCounts(max_depth + 1);
  std::vector<double> return_probability(max_depth + 2);
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

  for (int sigma = 1; sigma < max_depth + 2; sigma++) {
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
  for (int sigma = 0; sigma < max_depth + 2; sigma++) {
    if (levelProbabilities[sigma].IsKey(start_node)) {
      return_probability[sigma] = levelProbabilities[sigma].GetDat(start_node);
    } else {
      return_probability[sigma] = 0.0;
    }
  }

  // Extract the dimensionality
  int nodeId, startLevel;
  for (int sigma = 2; sigma < max_depth; sigma++) {
    // we to something like the five-point stencil for the averaging
    mean_probability = 0.5 * (return_probability[sigma - 1] + return_probability[sigma + 1]);
    // mean_probability =
    //    0.25 * (return_probability[sigma - 1] + 2 * return_probability[sigma] + return_probability[sigma + 1]);
    // symmetric five-point stencil - cf. wikpedia "Finite difference coefficents"
    probability_derivative = 0.5 * (-return_probability[sigma - 1] + return_probability[sigma + 1]);
    // probability_derivative = (return_probability[sigma - 2] - 8 * return_probability[sigma - 1] +
    //                          8 * return_probability[sigma + 1] - return_probability[sigma + 2]) /
    //                         12.0;
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
