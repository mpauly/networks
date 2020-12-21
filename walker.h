#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <vector>

#define WALKER_DEBUG false

// For a given startNode this method returns a hashmap of vectors, such that the keys of the hashmap
// corresond to the Ids of the nodes that one wants to study and the vectors correspond to the probability
// probability_trail for one node. The implementation probably is not super efficent, complexity grows exponentially.
// The implementation roughly corresponds to a breadth-first search tree.
// Inspired by DoBfs
template <class PGraph>
std::vector<double> spectralDimensionAtNode(const PGraph &Graph, const int &start_node, const int &max_depth,
                                            std::function<void(int)> progress_monitor,
                                            const double diffusion_constant) {
  // setup data structures
  TSnapQueue<int> queue;
  // std::vector<THash<TInt, double>> levelProbabilities(max_depth + 2);
  THash<TInt, double> last_lvl_probabilities;
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
  last_lvl_probabilities.AddDat(start_node, 1.0);
  return_probability[0] = 1.0;

  for (int sigma = 1; sigma < max_depth + 2; sigma++) {
    THash<TInt, double> lvl_probabilities;

    if (WALKER_DEBUG)
      printf("\n == Sigma %d ==\n", sigma);
    progress_monitor(sigma);
    // deal with all the nodes that are in the queue at the moment
    // this relies on the loop initialization begin run exactly once, as the queue grows while the loop is executed
    for (int node = queue.Len(); node > 0; node--) {
      // take a node out of the queue
      const int nodeId = queue.Top();
      if (WALKER_DEBUG)
        printf("   Node %d   ", nodeId);
      queue.Pop();
      const typename PGraph::TObj::TNodeI NodeI = Graph->GetNI(nodeId);

      const int nodeDegree = NodeI.GetOutDeg();
      // the diffusion constant term at the same node
      if (lvl_probabilities.IsKey(nodeId)) {
        node_probability = lvl_probabilities.GetDat(nodeId);
      } else {
        node_probability = 0;
        queue.Push(nodeId);
      }
      node_probability += (1.0 - diffusion_constant) * last_lvl_probabilities.GetDat(nodeId);
      lvl_probabilities.AddDat(nodeId, node_probability);
      // loop over child nodes
      for (v = 0; v < nodeDegree; v++) {
        const int childId = NodeI.GetOutNId(v);
        // did we already touch this node? If yes just update existing entry
        // if no then create a new entry and add the node to the queue
        if (lvl_probabilities.IsKey(childId)) {
          node_probability = lvl_probabilities.GetDat(childId);
        } else {
          node_probability = 0.0;
          queue.Push(childId);
        }
        node_probability += diffusion_constant * last_lvl_probabilities.GetDat(nodeId) / ((double)nodeDegree);
        if (WALKER_DEBUG)
          printf(" %d : %f ", childId, node_probability);
        // AddDat also overwrites an existing entry
        lvl_probabilities.AddDat(childId, node_probability);
      }

      if (lvl_probabilities.IsKey(start_node)) {
        return_probability[sigma] = lvl_probabilities.GetDat(start_node);
      } else {
        return_probability[sigma] = 0.0;
      }
    }
    last_lvl_probabilities = lvl_probabilities;
  }

  // Extract the dimensionality
  if (WALKER_DEBUG) {
    printf(" == Dimension == \n");
  }
  int nodeId, startLevel;
  for (int sigma = 1; sigma < max_depth; sigma++) {
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
    if (WALKER_DEBUG)
      printf("%d %f \n", sigma, dimension[sigma]);
  }

  return dimension;
}
