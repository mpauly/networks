#include "walker/base.h"
#include "walker/consts.h"
#include "walker/io.h"
#include <filesystem>
#include <iostream>
#include <string>
namespace fs = std::filesystem;

int main(int argc, char *argv[]) {
  std::string graph_name;
  for (; optind < argc; optind++) {
    graph_name = argv[optind];
  }

  std::string dimension_file = walker::DIMENSION_DIR + graph_name + ".dat";
  std::cout << "- will write to " << dimension_file << std::endl;

  std::ofstream dimfile;
  dimfile.open(dimension_file.c_str(), std::ofstream::out);
  dimfile << "# merged dimensions for graph " << graph_name << std::endl;
  dimfile << "# format: start_node sigma d_spec" << std::endl;

  std::string path = walker::WALK_DIR + graph_name;
  walker::RandomWalk walk;
  for (const auto &entry : fs::directory_iterator(path)) {
    std::cout << "Reading file " << entry.path() << " with length ";
    walk = walker::importRandomWalkFromBinaryFile(entry.path(), true);
    // we are ignoring the first and the last data point
    for (int sig = 1; sig < walk.sigma; sig++) {
      dimfile << walk.start_node << "\t" << sig << "\t";
      dimfile << std::fixed << std::setprecision(8);
      dimfile << walk.dimension[sig] << "\n";
    }
    std::cout << walk.sigma << std::endl;
  }
}