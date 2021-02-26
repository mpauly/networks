#include <string>
#ifndef WALKER_GRAPHS_SAVE_H
#define WALKER_GRAPHS_SAVE_H

template <class Graph> void save_graph_to_file(Graph G, std::string filename) {
  TFOut FOut(filename.c_str());
  G->Save(FOut);
  FOut.Flush();
}
#endif