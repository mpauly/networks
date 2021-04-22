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

// The parametrization and the formulars correspond to
// 2010 Krioukov et.al. - Hyperbolic geometry of complex networks

// returns the cosh(d(x,y)) of the distance between two points x and y
inline double hyperbolic_distance(const std::vector<double> &x, const std::vector<double> &y) {
  return cosh(x[0]) * cosh(y[0]) - sinh(x[0]) * sinh(y[0]) * cos(M_PI - abs(M_PI - abs(x[1] - y[1])));
}

inline double rho_of_r(const double &r, const double &alpha, const double &radius) {
  return alpha * sinh(alpha * r) / (cosh(alpha * radius) - 1);
}

std::vector<std::vector<double>> sprinkle_hyperbolic(const int numberPoints, const double radius) {
  std::default_random_engine generator;
  std::uniform_real_distribution<double> uni_dist(0.0, 1.0);

  std::cout << "Making hyperbolic geometry with " << numberPoints << " points" << std::endl;
  std::vector<std::vector<double>> positions(numberPoints, std::vector<double>(2, 0));
  // sprinkle the points - simple rejection sampling for now
  int point = 0;
  const double alpha = 1.0;
  const double rho_max = rho_of_r(radius, alpha, radius);
  while (point < numberPoints) {
    const double r = radius * uni_dist(generator);
    const double theta = 2. * M_PI * uni_dist(generator);
    const double rand = rho_max * uni_dist(generator);
    const double rho = rho_of_r(r, alpha, radius);
    if (rand < rho) {
      positions[point][0] = r;
      positions[point][1] = theta;
      point++;
    }
  }
  return positions;
}

int construct_edges(std::vector<std::vector<double>> &positions, std::function<void(int, int)> add_edge,
                    double cutoff_distance) {
  const int numberPoints = positions.size();
  const int ten_percent = numberPoints * 10 / 100;
  double distance_squared;
  bool edge;
  int edge_count = 0;
  static double cutoff_distance_acosh = cosh(cutoff_distance);
  std::cout << "- Constructing edges" << std::endl;
  std::cout << "- cutoff is cosh(" << cutoff_distance << ")=" << cutoff_distance_acosh << std::endl;

  for (int i = 0; i < numberPoints; i++) {
    if (i % ten_percent == 0) {
      std::cout << i * 100 / numberPoints << "% done" << std::endl;
    }
#pragma omp parallel for
    for (int j = 0; j < i; j++) {
      if (hyperbolic_distance(positions[i], positions[j]) < cutoff_distance_acosh) {
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

// hyperbolic disk
void generate_hyperbolic_disk() {
  PUNGraph Graph = TUNGraph::New();
  std::ofstream nodefile;
  nodefile.open("nodes_hyperbolic.csv", std::ofstream::out);
  std::ofstream edgefile;
  edgefile.open("edges_hyperbolic.csv", std::ofstream::out);

  const int numberPoints = 1e6;
  const double average_degree = 4.0;
  const double radius = 2.0 * log(8 * numberPoints / M_PI / average_degree);
  const double edge_cutoff = radius;

  std::cout << "- Sprinkling on a disk of radius " << radius << std::endl;
  std::vector<std::vector<double>> positions = sprinkle_hyperbolic(numberPoints, radius);

  for (int point = 0; point < numberPoints; point++) {
    IAssert(Graph->AddNode(point) == point);
    nodefile << point << "\t" << positions[point][0] << "\t" << positions[point][1] << "\n";
  }
  std::function<void(int, int)> add_edge = [&Graph, &positions, &edgefile](int i, int j) {
    Graph->AddEdge(i, j);
    edgefile << positions[i][0] << "\t" << positions[i][1] << "\t" << positions[j][0] << "\t" << positions[j][1]
             << "\n";
  };
  int nr_of_edges = construct_edges(positions, add_edge, edge_cutoff);
  Graph->Defrag();
  std::cout << "- Writing graph with " << nr_of_edges << " edges" << std::endl;
  save_graph_to_file(Graph, walker::GRAPH_DIR + "hyperbolic_disk.dat");
}

// network cosmology
void generate_growing_disk() {
  PUNGraph Graph = TUNGraph::New();
  std::ofstream nodefile;
  nodefile.open("nodes_hyperbolic.csv", std::ofstream::out);
  std::ofstream edgefile;
  edgefile.open("edges_hyperbolic.csv", std::ofstream::out);

  const int numberPoints = 1e6;
  const double r_cutoff = 3.0;
  const double r_cutoff_cosh = cosh(r_cutoff);

  std::default_random_engine generator;
  std::uniform_real_distribution<double> uni_dist(0.0, 1.0);

  std::vector<std::vector<double>> positions(numberPoints, {0.0, 0.0});
  const int ten_percent = numberPoints / 10;

  for (int point = 0; point < numberPoints; point++) {
    if (point % ten_percent == 0) {
      std::cout << point * 100 / numberPoints << "% done" << std::endl;
    }
    positions[point][0] = log(point);
    positions[point][1] = 2.0 * M_PI * uni_dist(generator);
    IAssert(Graph->AddNode(point) == point);
    nodefile << point << "\t" << positions[point][0] << "\t" << positions[point][1] << "\n";

    for (int j = 0; j < point; j++) {
      if (hyperbolic_distance(positions[point], positions[j]) < r_cutoff_cosh) {
        Graph->AddEdge(point, j);
        edgefile << positions[point][0] << "\t" << positions[point][1] << "\t" << positions[j][0] << "\t"
                 << positions[j][1] << "\n";
      }
    }
  }
  Graph->Defrag();
  Graph = TSnap::GetMxWcc(Graph);
  std::cout << "- Writing graph with " << Graph->GetNodes() << " nodes and " << Graph->GetEdges() << " edges"
            << std::endl;
  save_graph_to_file(Graph, walker::GRAPH_DIR + "growing_disk.dat");
}

int main(int argc, char *argv[]) {
  std::map<std::string, std::function<void()>> known_functions = {
      {"hyperbolic_disk", generate_hyperbolic_disk},
      {"growing_disk", generate_growing_disk},
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
        std::cout << " - use ./hyperbolic.x help to list all available targets" << std::endl;
      } else {
        std::cout << "Making " << key << std::endl;
        element->second();
      }
    }
  }
  if (count == 0) {
    std::cout << " No target given" << std::endl;
    std::cout << " - use ./hyperbolic.x help to list all available sets" << std::endl;
    return 1;
  } else {
    std::cout << std::endl << "Done writing graphs" << std::endl;
    return 0;
  }
}