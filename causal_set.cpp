#include "Snap.h"
#include "walker/base.h"
#include "walker/consts.h"
#include "walker/graphs/save_to_file.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <string>
#define _USE_MATH_DEFINES
#include <cmath>

inline double minkowski_line_element(const std::vector<double> &x, const std::vector<double> &y) {
  return (x[0] - y[0]) * (x[0] - y[0]) - (x[1] - y[1]) * (x[1] - y[1]);
}

inline double euclidean_distance(const std::vector<double> &x, const std::vector<double> &y) {
  return (x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1]);
}

inline bool accept_all(const double &t, const double &x) { return true; }

inline bool flat_connected(const std::vector<double> &x, const std::vector<double> &y) {
  if (minkowski_line_element(x, y) < 0)
    return false;
  return true;
}

inline bool x_before_y(const std::vector<double> &x, const std::vector<double> &y) { return x[0] < y[0]; }

std::vector<std::vector<double>> sprinkle_minkowski(const int numberPoints) {
  // const int nV = 1e6;
  std::default_random_engine generator;
  // std::poisson_distribution<int> poisson_distribution(nV);
  std::uniform_real_distribution<double> uni_dist(0.0, 1.0);

  // const int numberPoints = 1e6; // poisson_distribution(generator);

  std::cout << "Making causal set with " << numberPoints << " points" << std::endl;
  std::vector<std::vector<double>> positions(numberPoints, std::vector<double>(2, 0));
  // sprinkle the points - simple rejection sampling for now
  int point = 0;
  while (point < numberPoints) {
    const double t = uni_dist(generator);
    const double x = uni_dist(generator);
    if (accept_all(t, x)) {
      positions[point][0] = t;
      positions[point][1] = x;
      point++;
    }
  }
  // sort based on time
  std::sort(positions.begin(), positions.end(), x_before_y);
  return positions;
}

std::vector<std::vector<double>> sprinkle_desitter(const int numberPoints) {
  std::default_random_engine generator;
  std::uniform_real_distribution<double> uni_dist(0.0, 1.0);
  std::vector<std::vector<double>> positions(numberPoints, std::vector<double>(2, 0));

  int point = 0;
  // this allows to scale eta
  const double scaling = M_PI / 2.0;
  while (point < numberPoints) {
    const double eta = std::atan(scaling * uni_dist(generator));
    const double theta = 2.0 * M_PI * uni_dist(generator);
    if (accept_all(eta, theta)) {
      positions[point][0] = eta;
      positions[point][1] = theta;
      point++;
    }
  }
  // sort based on time
  std::sort(positions.begin(), positions.end(), x_before_y);
  return positions;
}

std::vector<std::vector<double>> sprinkle_anti_desitter(const int numberPoints) {
  std::default_random_engine generator;
  std::uniform_real_distribution<double> uni_dist(0.0, 1.0);
  std::vector<std::vector<double>> positions(numberPoints, std::vector<double>(2, 0));

  int point = 0;
  // this allows to scale eta
  const double scaling = M_PI / 2.0;
  while (point < numberPoints) {
    const double rho = std::atan(scaling * uni_dist(generator));
    const double eta = 2.0 * M_PI * uni_dist(generator);
    if (accept_all(eta, rho)) {
      positions[point][0] = eta;
      positions[point][1] = rho;
      point++;
    }
  }
  // sort based on time
  std::sort(positions.begin(), positions.end(), x_before_y);
  return positions;
}

int construct_flat_edges(std::vector<std::vector<double>> &positions, std::function<void(int, int)> add_edge) {
  const int numberPoints = positions.size();
  const int ten_percent = numberPoints * 10 / 100;
  double distance_squared;
  bool edge;
  int edge_count = 0;
  std::cout << "- Transitive reduction" << std::endl;

  for (int i = 0; i < numberPoints; i++) {
    if (i % ten_percent == 0) {
      std::cout << "  " << i * 100 / numberPoints << "% done" << std::endl;
    }
#pragma omp parallel for
    for (int j = 0; j < i; j++) {
      // these two points are space-like
      if (!flat_connected(positions[i], positions[j])) {
        continue;
      }
      // sort the two points
      const bool i_before_j = x_before_y(positions[i], positions[j]);
      const int in = (i_before_j) ? i : j;
      const int out = (i_before_j) ? j : i;
      // first assume they are directly connected
      edge = true;
      // check all points in between the two points - this assumes a sorted input
      for (int k = in + 1; k < out - 1; k++) {
        // the relation actually is transitive
        // here the order of the calls is optimised for performance
        // this is why it looks slightly weird
        if (flat_connected(positions[in], positions[k]) && flat_connected(positions[k], positions[out])) {
          edge = false;
          break;
        }
      }
      if (edge) {
#pragma omp critical
        {
          add_edge(i, j);
          edge_count++;
        }
      }
    }
  }
  return edge_count;
}

// simple Minkowski space
void generate_minkowski() {
  PUNGraph Graph = TUNGraph::New();
  PUNGraph GraphCutoff1 = TUNGraph::New();
  PUNGraph GraphCutoff2 = TUNGraph::New();
  std::ofstream nodefile;
  nodefile.open("nodes.csv", std::ofstream::out);
  std::ofstream edgefile;
  edgefile.open("edges.csv", std::ofstream::out);

  const int numberPoints = 1e6;
  const double euclidean_cutoff1 = 30. / sqrt(numberPoints);
  const double euclidean_cutoff2 = 100. / sqrt(numberPoints);

  std::vector<std::vector<double>> positions = sprinkle_minkowski(numberPoints);

  for (int point = 0; point < numberPoints; point++) {
    IAssert(Graph->AddNode(point) == point);
    IAssert(GraphCutoff1->AddNode(point) == point);
    IAssert(GraphCutoff2->AddNode(point) == point);
    nodefile << point << "\t" << positions[point][0] << "\t" << positions[point][1] << "\n";
  }
  std::function<void(int, int)> add_edge = [&Graph, &GraphCutoff1, &GraphCutoff2, &positions, euclidean_cutoff1,
                                            euclidean_cutoff2, &edgefile](int i, int j) {
    Graph->AddEdge(i, j);
    const double euclidean_dist = euclidean_distance(positions[i], positions[j]);
    if (euclidean_dist < euclidean_cutoff1) {
      GraphCutoff1->AddEdge(i, j);
    }
    if (euclidean_dist < euclidean_cutoff2) {
      GraphCutoff2->AddEdge(i, j);
    }
    edgefile << positions[i][0] << "\t" << positions[i][1] << "\t" << positions[j][0] << "\t" << positions[j][1] << "\t"
             << euclidean_dist << "\n";
  };
  int nr_of_edges = construct_flat_edges(positions, add_edge);
  Graph->Defrag();
  GraphCutoff1->Defrag();
  GraphCutoff2->Defrag();
  std::cout << "- Writing graph with " << nr_of_edges << " edges" << std::endl;
  save_graph_to_file(Graph, walker::GRAPH_DIR + "causal_set.dat");
  save_graph_to_file(GraphCutoff1, walker::GRAPH_DIR + "causal_set_l30.dat");
  save_graph_to_file(GraphCutoff2, walker::GRAPH_DIR + "causal_set_l100.dat");
}

// (anti-)desitter space
void generate_hyperbolic(bool anti) {
  PUNGraph Graph = TUNGraph::New();
  std::ofstream nodefile;
  nodefile.open("nodes.csv", std::ofstream::out);
  std::ofstream edgefile;
  edgefile.open("edges.csv", std::ofstream::out);

  const int numberPoints = 1e3;
  std::vector<std::vector<double>> positions;

  if (anti) {
    positions = sprinkle_anti_desitter(numberPoints);
  } else {
    positions = sprinkle_desitter(numberPoints);
  }
  for (int point = 0; point < numberPoints; point++) {
    IAssert(Graph->AddNode(point) == point);
    nodefile << point << "\t" << positions[point][0] << "\t" << positions[point][1] << "\n";
  }
  std::function<void(int, int)> add_edge = [&Graph, &positions, &edgefile](int i, int j) {
    Graph->AddEdge(i, j);
    // to be implemented
    const double euclidean_dist = 0;
    edgefile << positions[i][0] << "\t" << positions[i][1] << "\t" << positions[j][0] << "\t" << positions[j][1] << "\t"
             << euclidean_dist << "\n";
  };
  int nr_of_edges = construct_flat_edges(positions, add_edge);
  Graph->Defrag();
  std::cout << "- Writing graph with " << nr_of_edges << " edges" << std::endl;
  if (anti) {
    save_graph_to_file(Graph, walker::GRAPH_DIR + "causal_set_antidesitter.dat");
  } else {
    save_graph_to_file(Graph, walker::GRAPH_DIR + "causal_set_desitter.dat");
  }
}

void generate_desitter() { generate_hyperbolic(false); }
void generate_anti_desitter() { generate_hyperbolic(true); }

void average_shortest_path() {
  std::ofstream outfile;
  outfile.open("data/average_path.tsv", std::ofstream::out);

  const std::vector<int> numberPointSet = {500, 1000, 5000, 10000, 50000, 100000};
  const int nr_of_sets = 10;

  for (int numberPoints : numberPointSet) {
    std::cout << "- Nr. of points: " << numberPoints << std::endl;

#pragma omp parallel for
    for (int set = 0; set < nr_of_sets; set++) {
      PUNGraph Graph = TUNGraph::New();
      std::vector<std::vector<double>> positions = sprinkle_minkowski(numberPoints);

      for (int point = 0; point < numberPoints; point++) {
        IAssert(Graph->AddNode(point) == point);
      }
      std::function<void(int, int)> add_edge = [&Graph, &positions](int i, int j) { Graph->AddEdge(i, j); };
      int nr_of_edges = construct_flat_edges(positions, add_edge);
      Graph->Defrag();

      std::cout << "Computing distances" << std::endl;
      double average_distance = 0;
      int nr_of_paths = 0;
      const int nr_of_paths_expected = (numberPoints * numberPoints - numberPoints) / 2;
      const int ten_percent = nr_of_paths_expected / 10;
      for (int i = 1; i < numberPoints; i++) {
        // we do this per point in order to avoid problems with averages of small/large numbers
        std::vector<double> distances(i, 0.0);
        for (int j = 0; j < i; j++) {
          nr_of_paths++;
          if (i % ten_percent == 0) {
            std::cout << "  " << i * 100 / numberPoints << "% done" << std::endl;
          }
          distances[j] = TSnap::GetShortPath(Graph, i, j);
        }
        average_distance += std::accumulate(distances.begin(), distances.end(), 0.0) / i;
      }
#pragma omp critical
      {
        std::cout << "  Finished set with " << nr_of_paths << " distances computed" << std::endl;
        outfile << numberPoints << "\t" << set << "\t" << average_distance << std::endl;
      }
    }
  }
}

int main(int argc, char *argv[]) {
  std::map<std::string, std::function<void()>> known_functions = {
      {"minkowski", generate_minkowski},
      {"desitter", generate_desitter},
      {"antidesitter", generate_anti_desitter},
      {"shortestPath", average_shortest_path},
  };

  int count = 0;
  for (; optind < argc; optind++) {
    count++;
    std::string key = argv[optind];
    if (key == "all") {
      std::cout << "Making all graphs" << std::endl;
      for (std::pair<std::string, std::function<void()>> element : known_functions) {
        std::cout << "- making " << element.first << std::endl;
        element.second();
      }
    } else if (key == "help" || key == "list") {
      std::cout << "Available targets are: " << std::endl;
      for (std::pair<std::string, std::function<void()>> element : known_functions) {
        std::cout << "- " << element.first << std::endl;
      }
      return 0;
    } else {
      auto element = known_functions.find(key);
      if (element == known_functions.end()) {
        std::cout << " - Target " << key << " not found - SKIPPING!" << std::endl;
        std::cout << " - use ./causal_set.x help to list all available targets" << std::endl;
      } else {
        std::cout << "Making " << key << std::endl;
        element->second();
      }
    }
  }
  if (count == 0) {
    std::cout << " No target given" << std::endl;
    std::cout << " - use ./causal_set.x help to list all available sets" << std::endl;
    return 1;
  } else {
    std::cout << std::endl << "Done writing sets" << std::endl;
    return 0;
  }
}