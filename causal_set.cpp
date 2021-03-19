#include "Snap.h"
#include "walker/base.h"
#include "walker/consts.h"
#include "walker/graphs/save_to_file.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>

double minkowski_line_element(const std::vector<double> &x, const std::vector<double> &y) {
  return (x[0] - y[0]) * (x[0] - y[0]) - (x[1] - y[1]) * (x[1] - y[1]);
}

double euclidean_distance(const std::vector<double> &x, const std::vector<double> &y) {
  return (x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1]);
}

bool accept_all(const double &t, const double &x) { return true; }

bool connected(const std::vector<double> &x, const std::vector<double> &y) {
  if (minkowski_line_element(x, y) < 0)
    return false;
  return true;
}
bool x_before_y(const std::vector<double> &x, const std::vector<double> &y) { return x[0] < y[0]; }

int main(int argc, char *argv[]) {

  std::ofstream nodefile;
  nodefile.open("nodes.csv", std::ofstream::out);
  std::ofstream edgefile;
  edgefile.open("edges.csv", std::ofstream::out);

  // const int nV = 100;
  std::default_random_engine generator;
  // std::poisson_distribution<int> poisson_distribution(nV);
  std::uniform_real_distribution<double> uni_dist(0.0, 1.0);

  const int numberPoints = 1e6; // poisson_distribution(generator);
  const double euclidean_cutoff1 = 30. / sqrt(numberPoints);
  const double euclidean_cutoff2 = 100. / sqrt(numberPoints);

  std::cout << "Making causal set with " << numberPoints << " points" << std::endl;
  std::vector<std::vector<double>> positions(numberPoints, std::vector<double>(2, 0));
  PUNGraph Graph = TUNGraph::New();
  PUNGraph GraphCutoff1 = TUNGraph::New();
  PUNGraph GraphCutoff2 = TUNGraph::New();
  // sprinkle the points - simple rejection sampling for now
  int point = 0;
  while (point < numberPoints) {
    const double t = uni_dist(generator);
    const double x = uni_dist(generator);
    if (accept_all(t, x)) {
      IAssert(Graph->AddNode(point) == point);
      IAssert(GraphCutoff1->AddNode(point) == point);
      IAssert(GraphCutoff2->AddNode(point) == point);
      positions[point][0] = t;
      positions[point][1] = x;
      nodefile << point << "\t" << t << "\t" << x << "\n";
      point++;
    }
  }
  // one could probably speed up the next step by sorting the points in a clever way
  const int ten_percent = numberPoints * 10 / 100;
  double distance_squared;
  int edge_count = 0;
  bool edge;
  std::cout << "- Transitive reduction" << std::endl;
  // the following loops are O(N^3) and could probably be done slightly more efficently
  for (int i = 0; i < numberPoints; i++) {
    if (i % ten_percent == 0) {
      std::cout << i * 100 / numberPoints << "% done" << std::endl;
    }
#pragma omp parallel for
    for (int j = 0; j < i; j++) {
      // these two points are space-like assuming a mostly negative metric
      if (!connected(positions[i], positions[j])) {
        continue;
      }
      // sort the two points
      const bool i_before_j = x_before_y(positions[i], positions[j]);
      const int in = (i_before_j) ? i : j;
      const int out = (i_before_j) ? j : i;
      // first assume they are directly connected
      edge = true;
      for (int k = 0; k < numberPoints; k++) {
        // the relation actually is transitive
        // here the order of the calls is optimised for performance
        // this is why it looks slightly weird
        if (x_before_y(positions[in], positions[k]) && x_before_y(positions[k], positions[out]) &&
            connected(positions[in], positions[k]) && connected(positions[k], positions[out])) {
          edge = false;
          break;
        }
      }
      if (edge) {
#pragma omp critical
        {
          Graph->AddEdge(i, j);
          const double euclidean_dist = euclidean_distance(positions[i], positions[j]);
          if (euclidean_dist < euclidean_cutoff1) {
            GraphCutoff1->AddEdge(i, j);
          }
          if (euclidean_dist < euclidean_cutoff2) {
            GraphCutoff2->AddEdge(i, j);
          }
          edgefile << positions[in][0] << "\t" << positions[in][1] << "\t" << positions[out][0] << "\t"
                   << positions[out][1] << "\t" << euclidean_dist << "\n";
          edge_count++;
        }
      }
    }
  }
  std::cout << "- Writing graph with " << edge_count << " edges" << std::endl;
  Graph->Defrag();
  GraphCutoff1->Defrag();
  GraphCutoff2->Defrag();
  save_graph_to_file(Graph, walker::GRAPH_DIR + "causal_set.dat");
  save_graph_to_file(Graph, walker::GRAPH_DIR + "causal_set_l30.dat");
  save_graph_to_file(Graph, walker::GRAPH_DIR + "causal_set_l100.dat");
}