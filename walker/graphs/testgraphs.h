#include "../consts.h"
#include "Snap.h"
#include "save_to_file.h"
#include <fstream>
#include <string>

#define INDEXAT2D(x, y, edge) ((x)*edge + y)
#define INDEXAT3D(x, y, z) ((x)*edge_length * edge_length + (y)*edge_length + z)

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
  for (int n = 0; n < edge_length - 1; n++) {
    for (int m = 0; m < edge_length - 1; m++) {
      G->AddEdge(INDEXAT2D(n, m, edge_length), INDEXAT2D(n + 1, m, edge_length));
      G->AddEdge(INDEXAT2D(n, m, edge_length), INDEXAT2D(n, m + 1, edge_length));
    }
  }
  save_graph_to_file(G, walker::GRAPH_DIR + "2d_lattice_" + to_string(edge_length) + ".dat");
}
// ================== 3D lattice ===================================
void make_3d_lattice_100() {
  const int edge_length = 100;
  printf("== Making 3D grid with length %d ==\n", edge_length);
  PUNGraph G = PUNGraph::TObj::New();
  for (int n = 0; n < edge_length * edge_length * edge_length; n++) {
    G->AddNode(); // if no parameter is given, node ids are 0,1,...
  }

  // Build a cube
  for (int n1 = 0; n1 < edge_length - 1; n1++) {
    for (int n2 = 0; n2 < edge_length - 1; n2++) {
      for (int n3 = 0; n3 < edge_length - 1; n3++) {
        G->AddEdge(INDEXAT3D(n1, n2, n3), INDEXAT3D(n1 + 1, n2, n3));
        G->AddEdge(INDEXAT3D(n1, n2, n3), INDEXAT3D(n1, n2 + 1, n3));
        G->AddEdge(INDEXAT3D(n1, n2, n3), INDEXAT3D(n1, n2, n3 + 1));
      }
    }
  }
  save_graph_to_file(G, walker::GRAPH_DIR + "3d_lattice_" + to_string(edge_length) + ".dat");
}