#include "walker/graphs/datagraphs.h"
#include "walker/graphs/testgraphs.h"
#include <map>

int main(int argc, char *argv[]) {
  std::map<std::string, std::function<void()>> known_functions = {
      // data graphs
      {"as_skitter", make_as_skitter},
      {"internet_caida_large", make_internet_caida_large},
      {"roadnet_pa", make_roadnet_pa},
      {"metabolism", make_metabolism},
      {"drosophila", make_drosophila},
      {"drosophila_large", make_drosophila_large},
      {"drosophila_large_network", make_drosophila_large_network},
      {"rat_voxel_brain", make_rat_voxel_brain},
      {"europe_osm", make_europe_osm},
      {"brain", make_brain},
      {"human_brain_regions_300", make_human_brain_regions_300},
      // test graphs
      {"1d_ring_26", make_1d_ring_26},
      {"1d_ring_100", make_1d_ring_100},
      {"1d_ring_100_random", make_1d_ring_100_random},
      {"2d_lattice_10", make_2d_lattice_10},
      {"2d_lattice_100", make_2d_lattice_100},
      {"3d_lattice_100", make_3d_lattice_100},
      {"3d_lattice_7", make_3d_lattice_7},
      {"watts_strogatz_2d", make_ws_2d},
      {"watts_strogatz_3d", make_ws_3d},
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
        std::cout << " - use ./make_graph.x help to list all available targets" << std::endl;
      } else {
        std::cout << "Making " << key << std::endl;
        element->second();
      }
    }
  }
  if (count == 0) {
    std::cout << " No target given" << std::endl;
    std::cout << " - use ./make_graph.x help to list all available targets" << std::endl;
    return 1;
  } else {
    std::cout << std::endl << "Done writing graphs" << std::endl;
    return 0;
  }
}
