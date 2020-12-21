#include "Snap.h"
#include <fstream>
#include <iostream>
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
    PNGraph G = TSnap::LoadEdgeList<PNGraph>("graphs/roadNet-PA.txt", 0, 1);
    // Here we are converting to an undirected graph for now
    save_graph_to_file(TSnap::ConvertGraph<PUNGraph>(G), graph_dir + "roadnet_pa.dat");
  }
  // =============== Internet topology ===============================
  {
    std::cout << "== Writing AS Skitter Internet topology graph ==" << std::endl;
    PUNGraph G = TSnap::LoadEdgeList<PUNGraph>("graphs/as-skitter.txt", 0, 1);
    // Here we are converting to an undirected graph for now
    save_graph_to_file(G, graph_dir + "as-skitter.dat");
  }
  // =============== Drosophila brain ===============================
  {
    std::cout << "== Writing Drosophila Brain graph ==" << std::endl;
    PUNGraph G = TSnap::LoadEdgeList<PUNGraph>("graphs/fly-drosophila_edges.txt", 0, 1);
    // Here we are converting to an undirected graph for now
    save_graph_to_file(G, graph_dir + "fly-drosophila.dat");
  }

  printf("\n Done writing graphs to ");
  printf(graph_dir.c_str());
  printf("\n");
  return 0;
}
