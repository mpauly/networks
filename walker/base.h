#ifndef WALKER_BASE_H
#define WALKER_BASE_H

#include "Snap.h"
#include <functional>
#include <vector>

#define WALKER_DEBUG false

namespace walker {

struct RandomWalk {
  int start_node;
  int sigma;
  double diffusion_constant;
  std::vector<double> return_probability;
  std::vector<double> dimension;
  THash<TInt, double> lvl_probabilities;
};

RandomWalk setupRandomWalk(const int &start_node, const double diffusion_constant) {
  RandomWalk walk;
  walk.return_probability.resize(1);
  walk.dimension.resize(1);
  walk.start_node = start_node;
  walk.sigma = 0;
  walk.diffusion_constant = diffusion_constant;
  walk.return_probability[0] = 1.0;
  walk.lvl_probabilities.AddDat(start_node, 1.0);
  return walk;
}

// legacy wrapper
template <class PGraph>
std::vector<double> spectralDimensionAtNode(const PGraph &Graph, const int &start_node, const int &max_depth,
                                            std::function<void(const RandomWalk)> progress_monitor,
                                            const double diffusion_constant) {
  // setup the walk
  RandomWalk walk = setupRandomWalk(start_node, diffusion_constant);
  progressRandomWalk(Graph, walk, max_depth, progress_monitor);
  return walk.dimension;
}

// The implementation roughly corresponds to a breadth-first search tree.
// Inspired by SNAP's DoBfs
template <class PGraph>
void progressRandomWalk(const PGraph &Graph, RandomWalk &walk, int nr_of_steps,
                        std::function<void(const RandomWalk &)> progress_monitor) {

  walk.return_probability.resize(walk.sigma + nr_of_steps + 1);
  walk.dimension.resize(walk.sigma + nr_of_steps + 1);

  TSnapQueue<int> queue;
  THash<TInt, double> last_lvl_probabilities = walk.lvl_probabilities;
  int childnode;
  double probability_derivative, node_probability;

  // setup the queue
  queue.Clr(false);
  for (auto it = last_lvl_probabilities.BegI(); it < last_lvl_probabilities.EndI(); it++) {
    IAssert(Graph->IsNode(it.GetKey()));
    queue.Push(it.GetKey());
  }

  int target_size = walk.sigma + nr_of_steps + 1;
  for (walk.sigma = walk.sigma + 1; walk.sigma < target_size; walk.sigma++) {
    walk.lvl_probabilities.Clr();

    if (WALKER_DEBUG)
      printf("\n == Sigma %d ==\n", walk.sigma);

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
      if (walk.lvl_probabilities.IsKey(nodeId)) {
        node_probability = walk.lvl_probabilities.GetDat(nodeId);
      } else {
        node_probability = 0;
        queue.Push(nodeId);
      }
      node_probability += (1.0 - walk.diffusion_constant) * last_lvl_probabilities.GetDat(nodeId);
      walk.lvl_probabilities.AddDat(nodeId, node_probability);
      // loop over child nodes
      for (childnode = 0; childnode < nodeDegree; childnode++) {
        const int childId = NodeI.GetOutNId(childnode);
        // did we already touch this node? If yes just update existing entry
        // if no then create a new entry and add the node to the queue
        if (walk.lvl_probabilities.IsKey(childId)) {
          node_probability = walk.lvl_probabilities.GetDat(childId);
        } else {
          node_probability = 0.0;
          queue.Push(childId);
        }
        node_probability += walk.diffusion_constant * last_lvl_probabilities.GetDat(nodeId) / ((double)nodeDegree);
        if (WALKER_DEBUG)
          printf(" %d : %f ", childId, node_probability);
        // AddDat also overwrites an existing entry
        walk.lvl_probabilities.AddDat(childId, node_probability);
      }

      if (walk.lvl_probabilities.IsKey(walk.start_node)) {
        walk.return_probability[walk.sigma] = walk.lvl_probabilities.GetDat(walk.start_node);
      } else {
        walk.return_probability[walk.sigma] = 0.0;
      }
    }
    last_lvl_probabilities = walk.lvl_probabilities;

    // finally extract the dimension
    if (walk.sigma > 1) {
      // This is the probability derivative and dimension at position sigma -1
      probability_derivative = 0.5 * (-walk.return_probability[walk.sigma - 2] + walk.return_probability[walk.sigma]);
      // this is the spectral dimension at sigma - 1
      walk.dimension[walk.sigma - 1] =
          -2.0 * ((double)walk.sigma - 1) / walk.return_probability[walk.sigma - 1] * probability_derivative;
      if (WALKER_DEBUG)
        printf("\n d=%f", walk.dimension[walk.sigma - 1]);
    }

    progress_monitor(walk);
  }
  walk.sigma -= 1; // take minus one, because the last step was never taken
}

// Convenience function without progress monitor
template <class PGraph> void progressRandomWalk(const PGraph &Graph, RandomWalk &walk, int nr_of_steps) {
  std::function<void(const RandomWalk)> progress_monitor = [](const RandomWalk &walk) {};
  progressRandomWalk(Graph, walk, nr_of_steps, progress_monitor);
}
} // namespace walker

#endif