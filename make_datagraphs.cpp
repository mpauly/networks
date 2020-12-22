#include "Snap.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;

void save_graph_to_file(PUNGraph G, std::string filename) {
  TFOut FOut(filename.c_str());
  G->Save(FOut);
  FOut.Flush();
}

int main(int argc, char *argv[]) {
  const string graph_dir = "graphs/";

  typedef PUNGraph PGraph; // undirected graph

  // =============== PA streets ===============================
  {
    std::cout << "== Writing Pensylvania street network graph ==" << std::endl;
    PNGraph G = TSnap::LoadEdgeList<PNGraph>("graphs/data/roadNet-PA.txt", 0, 1);
    // Here we are converting to an undirected graph for now
    save_graph_to_file(TSnap::ConvertGraph<PUNGraph>(G), graph_dir + "roadnet_pa.dat");
  }
  // =============== Internet topology ===============================
  {
    std::cout << "== Writing AS Skitter Internet topology graph ==" << std::endl;
    PUNGraph G = TSnap::LoadEdgeList<PUNGraph>("graphs/data/as-skitter.txt", 0, 1);
    // Here we are converting to an undirected graph for now
    save_graph_to_file(G, graph_dir + "as-skitter.dat");
  }
  // =============== Drosophila brain ===============================
  {
    std::cout << "== Writing Drosophila Brain graph ==" << std::endl;
    PUNGraph G = TSnap::LoadEdgeList<PUNGraph>("graphs/data/fly-drosophila_edges.txt", 0, 1);
    // Here we are converting to an undirected graph for now
    save_graph_to_file(G, graph_dir + "fly-drosophila.dat");
  }
  // =============== Europe OSM map ===============================
  {
    std::cout << "== Writing Europe OSM graph ==" << std::endl;
    fstream osmfile;
    int rows, columns, entries;
    osmfile.open("graphs/data/europe_osm/europe_osm.mtx", ios::in);
    if (osmfile.is_open()) {
      std::string line;
      int nr_of_nodes = 0, nr_of_edges = 0;
      while (std::getline(osmfile, line)) {
        if (line[0] == '%')
          continue;
        std::istringstream iss(line);
        iss >> rows, columns, entries;
        break;
      }

      PGraph G = PUNGraph::TObj::New();
      for (int n = 0; n < rows; n++) {
        G->AddNode();
        nr_of_nodes++;
      }
      int in_id, out_id;
      while (osmfile >> in_id >> out_id) {
        G->AddEdge(in_id - 1, out_id - 1);
        nr_of_edges++;
      }
      std::cout << "  # Nodes: " << nr_of_nodes << "  # Edges: " << nr_of_edges << std::endl;
      save_graph_to_file(G, graph_dir + "europe_osm.dat");
      osmfile.close();
    }
  }

  printf("\n Done writing graphs to ");
  printf(graph_dir.c_str());
  printf("\n");
  return 0;
}
