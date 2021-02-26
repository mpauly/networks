#include "walker/base.h"
#include "walker/consts.h"
#include "walker/io.h"
#include <filesystem>
#include <iostream>
#include <string>
namespace fs = std::filesystem;

void print_usage() {
  std::cout << "usage: ./inspect_walk.x [-f walk_id] graph_name" << std::endl
            << "Displays information on all walks for a given graph" << std::endl
            << "  -f walk_id: Filter for the walk with the given walk_id" << std::endl
            << "  graph_filename: The graph to study walks for in the format roadnet_pa" << std::endl;
}

int main(int argc, char *argv[]) {
  bool filter_walks = false;
  int walk_id = 0;
  // parse arguments
  int option;
  while ((option = getopt(argc, argv, "f:")) != -1) { // get option from the getopt() method
    switch (option) {
    case 'f':
      filter_walks = true;
      walk_id = atoi(optarg);
      break;
    case '?':
      std::cout << "Unknown option: " << optopt << std::endl;
      print_usage();
      return 1;
    }
  }

  std::string graph_name = "";
  for (; optind < argc; optind++) {
    graph_name = argv[optind];
  }

  if (graph_name.empty()) {
    print_usage();
    return 1;
  }

  const std::vector<std::string> headers = {"Walk ID", "Start Node", "Sigma", "Delta"};
  const std::string separator = "      ";

  size_t headerWidths[headers.size()];
  for (int i = 0; i < headers.size(); i++)
    headerWidths[i] = headers[i].size();

  std::cout << "== Walks for Graph " << graph_name << " == \n" << std::endl;
  if (filter_walks) {
    std::cout << " Filtering output to walk " << walk_id << " only" << std::endl << std::endl;
  }

  // Table header
  int table_width = 0;
  for (int i = 0; i < headers.size(); i++) {
    std::cout << headers[i] << separator;
    table_width += headers[i].size() + separator.size();
  }
  table_width -= separator.size();
  std::cout << std::endl;
  for (int i = 0; i < table_width; i++)
    std::cout << "-";
  std::cout << std::endl;

  // Table content
  std::string path = walker::WALK_DIR + graph_name;
  walker::RandomWalk walk;
  if (filter_walks) {
    walk = walker::importRandomWalkFromBinaryFile(path + "/" + std::to_string(walk_id) + ".bin.dat", true);
    std::cout << std::left << std::setw(headerWidths[0]) << walk_id << separator << std::setw(headerWidths[1])
              << walk.start_node << separator << std::setw(headerWidths[2]) << walk.sigma << separator
              << std::setw(headerWidths[3]) << walk.diffusion_constant << std::endl;
    std::cout << std::endl;
  } else {
    int nr_of_walks = 0;
    for (const auto &entry : fs::directory_iterator(path)) {
      walk = walker::importRandomWalkFromBinaryFile(entry.path(), true);
      std::string filename = fs::path(entry).filename().string();
      // the 8 just cuts away the .bin.dat
      std::cout << std::left << std::setw(headerWidths[0]) << filename.substr(0, filename.size() - 8) << separator
                << std::setw(headerWidths[1]) << walk.start_node << separator << std::setw(headerWidths[2])
                << walk.sigma << separator << std::setw(headerWidths[3]) << walk.diffusion_constant << std::endl;
      nr_of_walks++;
    }
    std::cout << std::endl << nr_of_walks << " walks in total" << std::endl;
  }
}
