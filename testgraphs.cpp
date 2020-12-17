#include "Snap.h"
#include "walker.h"
#include <fstream>

#define EDGE_LENGTH 100
#define INDEXAT2D(x, y) ((x)*EDGE_LENGTH + y)
#define INDEXAT3D(x, y, z) ((x)*EDGE_LENGTH * EDGE_LENGTH + (y)*EDGE_LENGTH + z)

// ------------------ the main function -----------------------
int main(int argc, char *argv[]) {
  typedef PUNGraph PGraph; // undirected graph

  // =============== 1D ring ===============================
  PGraph G = PGraph::TObj::New();
  int ring_length = 26;
  for (int n = 0; n < ring_length; n++) {
    G->AddNode(); // if no parameter is given, node ids are 0,1,...
  }

  for (int n = 0; n < ring_length; n++) {
    G->AddEdge(n, (n + 1) % ring_length);
  }

  std::ofstream dimfile;
  dimfile.open("data/dimension_1d_short.dat");
  dimfile << "# dimensions for a ring of length " << ring_length << "\n";
  auto treeCounts = spectralDimensionAtNode(G, 7, 150, dimfile);
  dimfile.close();

  // =============== 1D ring ===============================
  G = PGraph::TObj::New();
  ring_length = 100;
  for (int n = 0; n < ring_length; n++) {
    G->AddNode(); // if no parameter is given, node ids are 0,1,...
  }

  for (int n = 0; n < ring_length; n++) {
    G->AddEdge(n, (n + 1) % ring_length);
  }

  dimfile.open("data/dimension_1d.dat");
  dimfile << "# dimensions for a ring of length " << ring_length << "\n";
  treeCounts = spectralDimensionAtNode(G, 7, 150, dimfile);
  dimfile.close();

  // =============== 1D ring with random connections ===============================
  PGraph G1a = PGraph::TObj::New();
  ring_length = 100;
  for (int n = 0; n < ring_length; n++) {
    G1a->AddNode(); // if no parameter is given, node ids are 0,1,...
  }

  for (int n = 0; n < ring_length; n++) {
    G1a->AddEdge(n, (n + 1) % ring_length);
  }

  const int nr_of_random_connections = 50;

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

  // TSnap::DrawGViz(G, gvlDot, "graph.png");

  std::ofstream dimfile1a;
  dimfile1a.open("data/dimension_1d_connected.dat");
  dimfile1a << "# dimensions for a ring with " << nr_of_random_connections
            << " additional random connections of length " << ring_length << "\n";
  auto treeCounts1a = spectralDimensionAtNode(G, 7, 150, dimfile1a);
  dimfile1a.close();

  // ================== 2D lattice ===================================
  PGraph G2 = PGraph::TObj::New();
  for (int n = 0; n < EDGE_LENGTH * EDGE_LENGTH; n++) {
    G2->AddNode(); // if no parameter is given, node ids are 0,1,...
  }
  // Build a square
  for (int n = 0; n < EDGE_LENGTH - 1; n++) {
    for (int m = 0; m < EDGE_LENGTH - 1; m++) {
      G2->AddEdge(INDEXAT2D(n, m), INDEXAT2D(n + 1, m));
      G2->AddEdge(INDEXAT2D(n, m), INDEXAT2D(n, m + 1));
    }
  }
  std::ofstream dimfile2;
  dimfile2.open("data/dimension_2d.dat");
  dimfile2 << "# dimensions for a square \n";
  auto treeCounts2 = spectralDimensionAtNode(G2, INDEXAT2D(50, 50), 100, dimfile2);
  dimfile2.close();

  // ================== 3D lattice ===================================
  PGraph G3 = PGraph::TObj::New();
  for (int n = 0; n < EDGE_LENGTH * EDGE_LENGTH * EDGE_LENGTH; n++) {
    G3->AddNode(); // if no parameter is given, node ids are 0,1,...
  }

  // Build a cube
  for (int n1 = 0; n1 < EDGE_LENGTH - 1; n1++) {
    for (int n2 = 0; n2 < EDGE_LENGTH - 1; n2++) {
      for (int n3 = 0; n3 < EDGE_LENGTH - 1; n3++) {
        G3->AddEdge(INDEXAT3D(n1, n2, n3), INDEXAT3D(n1 + 1, n2, n3));
        G3->AddEdge(INDEXAT3D(n1, n2, n3), INDEXAT3D(n1, n2 + 1, n3));
        G3->AddEdge(INDEXAT3D(n1, n2, n3), INDEXAT3D(n1, n2, n3 + 1));
      }
    }
  }

  std::ofstream dimfile3;
  dimfile3.open("data/dimension_3d.dat");
  dimfile3 << "# dimensions for a cube \n";
  auto treeCounts3 = spectralDimensionAtNode(G3, INDEXAT3D(50, 50, 50), 100, dimfile3);
  dimfile3.close();

  // the following are approximate neighbourhood functions - might be
  // interesting at some point printf("== TestAnf ==\n");
  // TSnap::TestAnf<PUNGraph>();
  return 0;
}
