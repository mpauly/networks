#include "../consts.h"
#include "Snap.h"
#include "save_to_file.h"
#include <fstream>
#include <string>
#define _USE_MATH_DEFINES
#include <cmath>

#define INDEXAT2D(x, y, edge) ((x) * (edge) + y)
#define INDEXAT3D(x, y, z, edge) ((x) * (edge) * (edge) + (y) * (edge) + z)

using namespace std;

// =============== 1D ring short ===============================
void make_1d_ring_26() {
  PUNGraph G = PUNGraph::TObj::New();
  const int ring_length = 26;
  printf("== Making 1D ring with length %d ==\n", ring_length);
  for (int n = 0; n < ring_length; n++) {
    G->AddNode(); // if no parameter is given, node ids are 0,1,...
  }

  for (int n = 0; n < ring_length; n++) {
    G->AddEdge(n, (n + 1) % ring_length);
  }

  save_graph_to_file(G, walker::GRAPH_DIR + "1d_ring_" + to_string(ring_length) + ".dat");
}
// =============== 1D ring long ===============================
void make_1d_ring_100() {
  PUNGraph G = PUNGraph::TObj::New();
  const int ring_length = 100;
  printf("== Making 1D ring with length %d ==\n", ring_length);

  for (int n = 0; n < ring_length; n++) {
    G->AddNode(); // if no parameter is given, node ids are 0,1,...
  }

  for (int n = 0; n < ring_length; n++) {
    G->AddEdge(n, (n + 1) % ring_length);
  }
  save_graph_to_file(G, walker::GRAPH_DIR + "1d_ring_" + to_string(ring_length) + ".dat");
}
// =============== 1D ring with random connections ===============================
void make_1d_ring_100_random() {
  PUNGraph G = PUNGraph::TObj::New();
  const int ring_length = 100;

  const int nr_of_random_connections = 50;
  printf("== Making 1D ring with length %d and %d random connections ==\n", ring_length, nr_of_random_connections);
  for (int n = 0; n < ring_length; n++) {
    G->AddNode(); // if no parameter is given, node ids are 0,1,...
  }

  for (int n = 0; n < ring_length; n++) {
    G->AddEdge(n, (n + 1) % ring_length);
  }

  // We are adding random edges to our graph
  int edges_added = 0;
  while (edges_added < nr_of_random_connections) {
    const int NId1 = G->GetRndNId();
    const int NId2 = G->GetRndNId();
    // this if condition checks if an edge already exists
    if (G->AddEdge(NId1, NId2) != -2) {
      edges_added++;
    }
  }

  save_graph_to_file(G, walker::GRAPH_DIR + "1d_ring_" + to_string(ring_length) + "_with_random_" +
                            to_string(nr_of_random_connections) + ".dat");
}
// ================== 2D small lattice ===================================
void make_2d_lattice_10() {
  printf("== Making 2D grid with length %d ==\n", 10);
  PUNGraph G = PUNGraph::TObj::New();
  const int small_edge = 10;
  for (int n = 0; n < small_edge * small_edge; n++) {
    G->AddNode(); // if no parameter is given, node ids are 0,1,...
  }
  // Build a square
  for (int n = 0; n < small_edge - 1; n++) {
    for (int m = 0; m < small_edge - 1; m++) {
      G->AddEdge(INDEXAT2D(n, m, small_edge), INDEXAT2D(n + 1, m, small_edge));
      G->AddEdge(INDEXAT2D(n, m, small_edge), INDEXAT2D(n, m + 1, small_edge));
    }
  }
  save_graph_to_file(G, walker::GRAPH_DIR + "2d_lattice_" + to_string(small_edge) + ".dat");
}
// ================== 2D lattice ===================================
void make_2d_lattice_100() {
  const int edge_length = 100;
  printf("== Making 2D grid with length %d ==\n", edge_length);
  PUNGraph G = PUNGraph::TObj::New();
  for (int n = 0; n < edge_length * edge_length; n++) {
    G->AddNode(); // if no parameter is given, node ids are 0,1,...
  }
  // Build a square
  for (int n = 0; n < edge_length; n++) {
    for (int m = 0; m < edge_length; m++) {
      if (n + 1 < edge_length)
        G->AddEdge(INDEXAT2D(n, m, edge_length), INDEXAT2D(n + 1, m, edge_length));
      if (m + 1 < edge_length)
        G->AddEdge(INDEXAT2D(n, m, edge_length), INDEXAT2D(n, m + 1, edge_length));
    }
  }
  // and repair the surface
  const int last_ind = edge_length - 1;
  for (int n = 0; n < edge_length - 1; n++) {
    G->AddEdge(INDEXAT2D(last_ind, n, edge_length), INDEXAT2D(last_ind, n + 1, edge_length));
    G->AddEdge(INDEXAT2D(n, last_ind, edge_length), INDEXAT2D(n + 1, last_ind, edge_length));
  }
  save_graph_to_file(G, walker::GRAPH_DIR + "2d_lattice_" + to_string(edge_length) + ".dat");
}
// ================== 3D lattice ===================================
void make_3d_lattice(const int edge_length) {
  printf("== Making 3D grid with length %d ==\n", edge_length);
  PUNGraph G = PUNGraph::TObj::New();
  for (int n = 0; n < edge_length * edge_length * edge_length; n++) {
    G->AddNode(); // if no parameter is given, node ids are 0,1,...
  }

  // Build a cube
  for (int n1 = 0; n1 < edge_length; n1++) {
    for (int n2 = 0; n2 < edge_length; n2++) {
      for (int n3 = 0; n3 < edge_length; n3++) {
        if (n1 + 1 < edge_length)
          G->AddEdge(INDEXAT3D(n1, n2, n3, edge_length), INDEXAT3D(n1 + 1, n2, n3, edge_length));
        if (n2 + 1 < edge_length)
          G->AddEdge(INDEXAT3D(n1, n2, n3, edge_length), INDEXAT3D(n1, n2 + 1, n3, edge_length));
        if (n3 + 1 < edge_length)
          G->AddEdge(INDEXAT3D(n1, n2, n3, edge_length), INDEXAT3D(n1, n2, n3 + 1, edge_length));
      }
    }
  }
  save_graph_to_file(G, walker::GRAPH_DIR + "3d_lattice_" + to_string(edge_length) + ".dat");
}
void make_3d_lattice_100() { make_3d_lattice(100); }
void make_3d_lattice_7() { make_3d_lattice(7); }
// =================== multi-dimensional generalization of WS ===========
void make_ws_2d() {
  const int edge_length = 100;
  const int total_nodes = edge_length * edge_length;
  const double radius_cutoff = 2.1;
  const int shift = ((int)radius_cutoff) + 1;
  const double RewireProb = 0.001;

  const double PI = 3.141592653589793;
  TRnd Rnd = TInt::Rnd;
  const int estimate_nr_edges = (int)(PI * radius_cutoff * radius_cutoff / 4.0) + shift;
  std::cout << "== Making 2D Watts strogatz with edge length" << edge_length << std::endl;
  std::cout << "  Edges per node estimated to be " << estimate_nr_edges << std::endl;

  THashSet<TIntPr> EdgeSet(total_nodes * estimate_nr_edges);

  for (int node_x = 0; node_x < edge_length; node_x++) {
    for (int node_y = 0; node_y < edge_length; node_y++) {
      const int src_x = node_x;
      const int src_y = node_y;
      const int src = INDEXAT2D(node_x, node_y, edge_length);
      for (int shift_x = 0; shift_x <= shift; shift_x++) {
        for (int shift_y = 0; shift_y <= shift; shift_y++) {
          if (shift_x == 0 && shift_y == 0)
            continue;
          // Euclidean distance for now
          if (shift_x * shift_x + shift_y * shift_y < radius_cutoff * radius_cutoff) {
            const int dst_x = (src_x + shift_x) % edge_length;
            const int dst_y = (src_y + shift_y) % edge_length;
            int dst = INDEXAT2D(dst_x, dst_y, edge_length);
            if (Rnd.GetUniDev() < RewireProb) { // random edge
              do {
                dst = Rnd.GetUniDevInt(total_nodes);
              } while (dst == src || EdgeSet.IsKey(TIntPr(src, dst)));
            }
            EdgeSet.AddKey(TIntPr(src, dst));
          }
        }
      }
    }
  }

  PUNGraph GraphPt = TUNGraph::New();
  TUNGraph &Graph = *GraphPt;
  Graph.Reserve(total_nodes, EdgeSet.Len());
  int node;
  for (node = 0; node < total_nodes; node++) {
    IAssert(Graph.AddNode(node) == node);
  }
  for (int edge = 0; edge < EdgeSet.Len(); edge++) {
    Graph.AddEdge(EdgeSet[edge].Val1, EdgeSet[edge].Val2);
  }
  Graph.Defrag();

  save_graph_to_file(GraphPt, walker::GRAPH_DIR + "watts_strogatz_2d_" + to_string(edge_length) + ".dat");
}

void make_ws_3d() {
  const int edge_length = 100;
  const int total_nodes = edge_length * edge_length * edge_length;
  const double radius_cutoff = 2.1;
  const int shift = ((int)radius_cutoff) + 1;
  std::vector<double> rewiringProbs = {0.001, 0.005, 0.01, 0.05};

  for (double RewireProb : rewiringProbs) {
    TRnd Rnd = TInt::Rnd;
    // the following is just a rough upper bound on the number of edges
    const int estimate_nr_edges = (int)(M_PI * radius_cutoff * radius_cutoff * radius_cutoff / 6.0) + shift;
    std::cout << "== Making 3D Watts strogatz with edge length " << edge_length << std::endl;
    std::cout << "  Edges per node estimated to be " << estimate_nr_edges << std::endl;

    THashSet<TIntPr> EdgeSet(total_nodes * estimate_nr_edges);

    for (int node_x = 0; node_x < edge_length; node_x++) {
      for (int node_y = 0; node_y < edge_length; node_y++) {
        for (int node_z = 0; node_z < edge_length; node_z++) {
          const int src_x = node_x;
          const int src_y = node_y;
          const int src_z = node_z;
          const int src = INDEXAT3D(node_x, node_y, node_z, edge_length);
          for (int shift_x = 0; shift_x <= shift; shift_x++) {
            for (int shift_y = 0; shift_y <= shift; shift_y++) {
              for (int shift_z = 0; shift_z <= shift; shift_z++) {
                if (shift_x == 0 && shift_y == 0 && shift_z == 0)
                  continue;
                // Euclidean distance for now
                if (shift_x * shift_x + shift_y * shift_y + shift_z * shift_z < radius_cutoff * radius_cutoff) {
                  const int dst_x = (src_x + shift_x) % edge_length;
                  const int dst_y = (src_y + shift_y) % edge_length;
                  const int dst_z = (src_z + shift_z) % edge_length;
                  int dst = INDEXAT3D(dst_x, dst_y, dst_z, edge_length);
                  if (Rnd.GetUniDev() < RewireProb) { // random edge
                    do {
                      dst = Rnd.GetUniDevInt(total_nodes);
                    } while (dst == src || EdgeSet.IsKey(TIntPr(src, dst)));
                  }
                  EdgeSet.AddKey(TIntPr(src, dst));
                }
              }
            }
          }
        }
      }
    }

    PUNGraph GraphPt = TUNGraph::New();
    TUNGraph &Graph = *GraphPt;
    Graph.Reserve(total_nodes, EdgeSet.Len());
    int node;
    for (node = 0; node < total_nodes; node++) {
      IAssert(Graph.AddNode(node) == node);
    }
    for (int edge = 0; edge < EdgeSet.Len(); edge++) {
      Graph.AddEdge(EdgeSet[edge].Val1, EdgeSet[edge].Val2);
    }
    Graph.Defrag();

    save_graph_to_file(GraphPt, walker::GRAPH_DIR + "watts_strogatz_3d_" + to_string(edge_length) + "_" +
                                    to_string(RewireProb) + ".dat");
  }
}