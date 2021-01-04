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
    int nr_of_nodes = 0, nr_of_edges = 0;
    PGraph G = PUNGraph::TObj::New();
    if (osmfile.is_open()) {
      std::string line;

      while (std::getline(osmfile, line)) {
        if (line[0] == '%')
          continue;
        std::istringstream iss(line);
        iss >> rows, columns, entries;
        break;
      }

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

    std::cout << "== Reducing Europe OSM graph ==" << std::endl;

    int nodes_reduced = 0;
    int node0, node1;
    bool nodes_removed;

    // We are reducing the graph because it still contains many long roads composed of many nodes
    // 1) Looping over all edges
    // 2) If both nodes leading into an edge have degree two we can get rid of one of the nodes
    // 3) Rewiring then requires to add one edge between node j and the other child node of i
    for (int i = 0; i < nr_of_nodes; i++) {
      // did we already delete the node in a previous step?
      if (!G->IsNode(i))
        continue;
      const typename PGraph::TObj::TNodeI NodeI = G->GetNI(i);
      const int nodeIDegree = NodeI.GetOutDeg();
      if (nodeIDegree != 2)
        continue;

      do {
        const int j1 = NodeI.GetOutNId(0);
        const typename PGraph::TObj::TNodeI NodeJ1 = G->GetNI(j1);
        const int nodeJ1Degree = NodeJ1.GetOutDeg();
        const int j2 = NodeI.GetOutNId(1);
        const typename PGraph::TObj::TNodeI NodeJ2 = G->GetNI(j2);
        const int nodeJ2Degree = NodeJ2.GetOutDeg();

        // check if one of the two child nodes has degree two, if so remove it
        int removalNode;
        if (nodeJ1Degree == 2) {
          nodes_removed = true;
          removalNode = j1;
          const int nodeJ1Degree = NodeJ1.GetOutDeg();
        } else if (nodeJ2Degree == 2) {
          nodes_removed = true;
          removalNode = j2;
        } else {
          nodes_removed = false;
        }

        // if one of the child nodes are degree two we can get rid of that node
        if (nodes_removed) {
          const typename PGraph::TObj::TNodeI RemovedNode = G->GetNI(removalNode);
          const int childId0 = RemovedNode.GetOutNId(0);
          const int childId1 = RemovedNode.GetOutNId(1);
          // ... except if we already form a triangle
          if (childId0 == j1 || childId1 == j1 || childId0 == j2 || childId1 == j2) {
            nodes_removed = false;
            continue;
          }
          G->DelNode(removalNode);
          nodes_reduced++;
          nodes_removed = true;

          // rewire the edges
          if (childId0 != i && !G->IsEdge(childId0, i)) {
            G->AddEdge(childId0, i);
          }
          if (childId1 != i && !G->IsEdge(childId1, i)) {
            G->AddEdge(childId1, i);
          }
        }
      } while (nodes_removed);
    }
    std::cout << "  Reduced by " << nodes_reduced << " nodes - saving" << std::endl;
    save_graph_to_file(G, graph_dir + "europe_osm_reduced.dat");
  }

  // =============== Human Brain ===============================
  {
    std::cout << "== Writing Human brain graph ==" << std::endl;
    PUNGraph G = TSnap::LoadEdgeList<PUNGraph>("graphs/data/brain_edgeslist.txt", 0, 1);
    // Here we are converting to an undirected graph for now
    save_graph_to_file(G, graph_dir + "brain.dat");
  }

  std::cout << std::endl << "Done writing graphs to " << graph_dir << std::endl;
  return 0;
}
