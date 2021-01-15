#include "Snap.h"
#include <fstream>
#include <string>

#define EDGE_LENGTH 100
#define INDEXAT2D(x, y, edge) ((x)*edge + y)
#define INDEXAT3D(x, y, z) ((x)*EDGE_LENGTH * EDGE_LENGTH + (y)*EDGE_LENGTH + z)

using namespace std;

void save_graph_to_file(PUNGraph G, std::string filename) {
  TFOut FOut(filename.c_str());
  G->Save(FOut);
  FOut.Flush();
}

// ------------------ the main function -----------------------
int main(int argc, char *argv[]) {
  const string graph_dir = "graphs/";

  typedef PUNGraph PGraph; // undirected graph

  // =============== 1D ring short ===============================
  {
    PGraph G = PGraph::TObj::New();
    const int ring_length = 26;
    printf("== Making 1D ring with length %d ==\n", ring_length);
    for (int n = 0; n < ring_length; n++) {
      G->AddNode(); // if no parameter is given, node ids are 0,1,...
    }

    for (int n = 0; n < ring_length; n++) {
      G->AddEdge(n, (n + 1) % ring_length);
    }

    save_graph_to_file(G, graph_dir + "1d_ring_" + to_string(ring_length) + ".dat");
  }
  // =============== 1D ring long ===============================
  {
    PGraph G = PGraph::TObj::New();
    const int ring_length = 100;
    printf("== Making 1D ring with length %d ==\n", ring_length);

    for (int n = 0; n < ring_length; n++) {
      G->AddNode(); // if no parameter is given, node ids are 0,1,...
    }

    for (int n = 0; n < ring_length; n++) {
      G->AddEdge(n, (n + 1) % ring_length);
    }
    save_graph_to_file(G, graph_dir + "1d_ring_" + to_string(ring_length) + ".dat");
  }
  // =============== 1D ring with random connections ===============================
  {
    PGraph G = PGraph::TObj::New();
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

    save_graph_to_file(G, graph_dir + "1d_ring_" + to_string(ring_length) + "_with_random_" +
                              to_string(nr_of_random_connections) + ".dat");
  }
  // ================== 2D small lattice ===================================
  {
    printf("== Making 2D grid with length %d ==\n", 10);
    PGraph G = PGraph::TObj::New();
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
    save_graph_to_file(G, graph_dir + "2d_lattice_" + to_string(small_edge) + ".dat");
  }
  // ================== 2D lattice ===================================
  {
    printf("== Making 2D grid with length %d ==\n", EDGE_LENGTH);
    PGraph G = PGraph::TObj::New();
    for (int n = 0; n < EDGE_LENGTH * EDGE_LENGTH; n++) {
      G->AddNode(); // if no parameter is given, node ids are 0,1,...
    }
    // Build a square
    for (int n = 0; n < EDGE_LENGTH - 1; n++) {
      for (int m = 0; m < EDGE_LENGTH - 1; m++) {
        G->AddEdge(INDEXAT2D(n, m, EDGE_LENGTH), INDEXAT2D(n + 1, m, EDGE_LENGTH));
        G->AddEdge(INDEXAT2D(n, m, EDGE_LENGTH), INDEXAT2D(n, m + 1, EDGE_LENGTH));
      }
    }
    save_graph_to_file(G, graph_dir + "2d_lattice_" + to_string(EDGE_LENGTH) + ".dat");
  }
  // ================== 3D lattice ===================================
  {
    printf("== Making 3D grid with length %d ==\n", EDGE_LENGTH);
    PGraph G = PGraph::TObj::New();
    for (int n = 0; n < EDGE_LENGTH * EDGE_LENGTH * EDGE_LENGTH; n++) {
      G->AddNode(); // if no parameter is given, node ids are 0,1,...
    }

    // Build a cube
    for (int n1 = 0; n1 < EDGE_LENGTH - 1; n1++) {
      for (int n2 = 0; n2 < EDGE_LENGTH - 1; n2++) {
        for (int n3 = 0; n3 < EDGE_LENGTH - 1; n3++) {
          G->AddEdge(INDEXAT3D(n1, n2, n3), INDEXAT3D(n1 + 1, n2, n3));
          G->AddEdge(INDEXAT3D(n1, n2, n3), INDEXAT3D(n1, n2 + 1, n3));
          G->AddEdge(INDEXAT3D(n1, n2, n3), INDEXAT3D(n1, n2, n3 + 1));
        }
      }
    }
    save_graph_to_file(G, graph_dir + "3d_lattice_" + to_string(EDGE_LENGTH) + ".dat");
  }
  printf("\n Done writing graphs to ");
  printf(graph_dir.c_str());
  printf("\n");
  return 0;
}
