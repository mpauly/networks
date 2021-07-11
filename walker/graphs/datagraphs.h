#include "../consts.h"
#include "Snap.h"
#include "save_to_file.h"
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

// =============== Internet topology ===============================
void make_as_skitter() {
  std::cout << "== Writing AS Skitter Internet topology graph ==" << std::endl;
  PUNGraph G = TSnap::LoadEdgeList<PUNGraph>((walker::DATA_DIR + "as-skitter.txt").c_str(), 0, 1);
  save_graph_to_file(G, walker::GRAPH_DIR + "as-skitter.dat");
}
// =============== Internet topology large ===============================
void make_internet_caida_large() {
  std::cout << "== Writing Internet topology large graph ==" << std::endl;
  std::ifstream file(walker::DATA_DIR + "midar-iff.links.bz2", std::ios_base::in | std::ios_base::binary);
  boost::iostreams::filtering_stream<boost::iostreams::input> in;
  in.push(boost::iostreams::bzip2_decompressor());
  in.push(file);
  std::string line;
  std::string word;
  std::vector<int> nodes;

  PUNGraph G = PUNGraph::TObj::New();

  int count = 0;

  while (std::getline(in, line)) {
    // skip comments
    if (line[0] == '#')
      continue;
    std::istringstream iss(line);
    // this reads "link"
    iss >> word;
    // this reads "L42:"
    iss >> word;
    nodes.clear();
    while (iss >> word) {
      // get rid of the ip
      word = word.substr(0, word.find(':'));
      // get rid of the leading N
      word = word.substr(1);
      nodes.push_back(std::stoi(word));
    }
    for (int i = 0; i < nodes.size(); i++) {
      if (!G->IsNode(nodes[i])) {
        G->AddNode(nodes[i]);
      }
    }
    for (int i = 0; i < nodes.size(); i++) {
      for (int j = i + 1; j < nodes.size(); j++) {
        G->AddEdge(nodes[i], nodes[j]);
      }
    }
    count++;
    if (count % 10000 == 0)
      std::cout << "Read " << count << " links" << std::endl;
  }
  std::cout << "Total: " << count << " links" << std::endl;
  G = TSnap::GetMxWcc(G);
  save_graph_to_file(G, walker::GRAPH_DIR + "internet_caida_large.dat");
}

// =============== PA streets ===============================
void make_roadnet_pa() {
  std::cout << "== Writing Pensylvania street network graph ==" << std::endl;
  PUNGraph G = TSnap::LoadEdgeList<PUNGraph>((walker::DATA_DIR + "roadNet-PA.txt").c_str(), 0, 1);
  // Here we are converting to an undirected graph for now
  save_graph_to_file(TSnap::ConvertGraph<PUNGraph>(G), walker::GRAPH_DIR + "roadnet_pa.dat");
}

// =============== Metabolism network ===============================
void make_metabolism() {
  std::cout << "== Writing Metabolism network graph ==" << std::endl;
  PUNGraph G = TSnap::LoadEdgeList<PUNGraph>((walker::DATA_DIR + "metabolism.csv").c_str(), 0, 1, ',');
  G = TSnap::GetMxWcc(G);
  save_graph_to_file(G, walker::GRAPH_DIR + "metabolism.dat");
}
// =============== Drosophila brain ===============================
void make_drosophila() {
  std::cout << "== Writing Drosophila Brain graph - small version ==" << std::endl;
  PUNGraph G = TSnap::LoadEdgeList<PUNGraph>((walker::DATA_DIR + "fly-drosophila_edges.txt").c_str(), 0, 1);
  save_graph_to_file(G, walker::GRAPH_DIR + "fly-drosophila.dat");
}
// =============== Drosophila brain large ==========================
void make_drosophila_large() {
  std::cout << "== Writing Drosophila Brain graph - large version ==" << std::endl;
  PUNGraph G = TSnap::LoadEdgeList<PUNGraph>(
      (walker::DATA_DIR + "exported-traced-adjacencies-v1.2/traced-total-connections-renumbered.csv").c_str(), 0, 1,
      ',');
  save_graph_to_file(G, walker::GRAPH_DIR + "fly-drosophila-large.dat");
}
// =============== Drosophila brain large network ==========================
void make_drosophila_large_network() {
  std::cout << "== Writing Drosophila Brain network with weights - large version ==" << std::endl;
  auto G = TNodeEDatNet<TInt, TInt>::New();

  std::fstream osmfile;
  int rows, columns, entries;
  osmfile.open(walker::DATA_DIR + "exported-traced-adjacencies-v1.2/traced-total-connections-renumbered.csv",
               std::ios::in);

  if (osmfile.is_open()) {
    std::string line;
    std::string token;

    int in_id, out_id, weight;

    int id_counter;

    int node_count = 0;
    int edge_count = 0;
    while (std::getline(osmfile, line)) {
      if (line[0] == 'b')
        continue;
      std::istringstream iss(line);
      std::getline(iss, token, ',');
      in_id = std::stoi(token);
      std::getline(iss, token, ',');
      out_id = std::stoi(token);
      std::getline(iss, token, ',');
      weight = std::stoi(token);
      if (!G->IsNode(in_id)) {
        G->AddNode(in_id);
        node_count++;
      }
      if (!G->IsNode(out_id)) {
        G->AddNode(out_id);
        node_count++;
      }

      G->AddEdge(in_id, out_id, weight);
      G->AddEdge(out_id, in_id, weight);
      edge_count++;
    }

    std::cout << " Built network with " << node_count << " nodes and " << edge_count
              << " edges - precomputing nodecount next" << std::endl;
    // loop over nodes and compute summed weights for this node
    int childnode, nodecount;
    for (TNodeEDatNet<TInt, TInt>::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
      nodecount = 0;
      const int nodeDegree = NI.GetOutDeg();
      for (childnode = 0; childnode < nodeDegree; childnode++) {
        const int childId = NI.GetOutNId(childnode);
        nodecount += NI.GetOutEDat(childnode);
      }
      G->SetNDat(NI.GetId(), nodecount);
    }
    osmfile.close();

    save_graph_to_file(G, walker::NETWORK_DIR + "fly-drosophila-weighted.dat");
  }
}
// =============== Rat Voxel Brain ===============================
void make_rat_voxel_brain() {
  std::cout << "== Writing Rat Voxel Brain Network ==" << std::endl;
  auto G = TNodeEDatNet<TInt, TInt>::New();
  // we are scaling the floats to ints by multiplying by the following factor
  const int upscaling_factor = 1000000;

  std::fstream ratfile;
  int rows, columns, entries;
  ratfile.open(walker::DATA_DIR + "rat_voxel_brain.csv", std::ios::in);

  if (ratfile.is_open()) {
    std::string line;
    std::string token;

    int in_id, out_id, weight;

    int id_counter;

    int node_count = 0;
    int edge_count = 0;
    while (std::getline(ratfile, line)) {
      // skip lines that are too short
      if (line.length() < 5)
        continue;
      std::istringstream iss(line);
      try {
        std::getline(iss, token, ',');
        in_id = std::stoi(token);
        std::getline(iss, token, ',');
        out_id = std::stoi(token);
        std::getline(iss, token, ',');
        weight = (int)(upscaling_factor * std::stof(token));
      } catch (const std::exception &ex) {
        std::cerr << ex.what() << '\n';
        std::cerr << "line: " << line << std::endl;
      }
      if (!G->IsNode(in_id)) {
        G->AddNode(in_id);
        node_count++;
      }
      if (!G->IsNode(out_id)) {
        G->AddNode(out_id);
        node_count++;
      }

      G->AddEdge(in_id, out_id, weight);
      G->AddEdge(out_id, in_id, weight);
      edge_count++;
    }

    std::cout << " Built network with " << node_count << " nodes and " << edge_count
              << " edges - precomputing nodecount next" << std::endl;
    // loop over nodes and compute summed weights for this node
    int childnode, nodecount;
    for (TNodeEDatNet<TInt, TInt>::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
      nodecount = 0;
      const int nodeDegree = NI.GetOutDeg();
      for (childnode = 0; childnode < nodeDegree; childnode++) {
        const int childId = NI.GetOutNId(childnode);
        nodecount += NI.GetOutEDat(childnode);
      }
      G->SetNDat(NI.GetId(), nodecount);
    }
    ratfile.close();

    save_graph_to_file(G, walker::NETWORK_DIR + "rat_voxel_brain.dat");
  }
}
// =============== Human Connectome Project (HCP) netmat - regions of the brain ==========================
void make_human_brain_regions_300() {
  std::cout << "== Writing human_brain_regions_300 ==" << std::endl;
  auto G = TNodeEDatNet<TInt, TInt>::New();
  const int matrix_size = 300;
  std::array<std::array<int, matrix_size>, matrix_size> matrix;

  std::fstream matrixfile;
  std::string line, token;
  const int upscaling_factor = 100000;
  int weight;
  matrixfile.open(walker::DATA_DIR + "hcp1200_netmat_d300.tsv", std::ios::in);
  for (int i = 0; i < matrix_size; i++)
    G->AddNode(i);
  for (int i = 0; i < matrix_size; i++) {
    std::getline(matrixfile, line);
    std::istringstream iss(line);
    for (int j = 0; j < i; j++) {
      std::getline(iss, token, '\t');
      weight = abs((int)(upscaling_factor * std::stof(token)));
      G->AddEdge(i, j, weight);
      G->AddEdge(j, i, weight);
    }
  }

  int childnode, nodecount;
  for (TNodeEDatNet<TInt, TInt>::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
    nodecount = 0;
    const int nodeDegree = NI.GetOutDeg();
    for (childnode = 0; childnode < nodeDegree; childnode++) {
      const int childId = NI.GetOutNId(childnode);
      nodecount += NI.GetOutEDat(childnode);
    }
    G->SetNDat(NI.GetId(), nodecount);
  }

  save_graph_to_file(G, walker::NETWORK_DIR + "human_brain_regions_300.dat");
}
// =============== Europe OSM map ===============================
void make_europe_osm() {
  std::cout << "== Writing Europe OSM graph ==" << std::endl;
  std::fstream osmfile;
  int rows, columns, entries;
  osmfile.open(walker::DATA_DIR + "europe_osm/europe_osm.mtx", std::ios::in);
  int nr_of_nodes = 0, nr_of_edges = 0;
  PUNGraph G = PUNGraph::TObj::New();
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
    save_graph_to_file(G, walker::GRAPH_DIR + "europe_osm.dat");
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
    const typename PUNGraph::TObj::TNodeI NodeI = G->GetNI(i);
    const int nodeIDegree = NodeI.GetOutDeg();
    if (nodeIDegree != 2)
      continue;

    do {
      const int j1 = NodeI.GetOutNId(0);
      const typename PUNGraph::TObj::TNodeI NodeJ1 = G->GetNI(j1);
      const int nodeJ1Degree = NodeJ1.GetOutDeg();
      const int j2 = NodeI.GetOutNId(1);
      const typename PUNGraph::TObj::TNodeI NodeJ2 = G->GetNI(j2);
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
        const typename PUNGraph::TObj::TNodeI RemovedNode = G->GetNI(removalNode);
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
  save_graph_to_file(G, walker::GRAPH_DIR + "europe_osm_reduced.dat");
}
// =============== Human Brain ===============================
void make_brain() {
  std::cout << "== Writing Human brain graph ==" << std::endl;
  PUNGraph G = TSnap::LoadEdgeList<PUNGraph>((walker::DATA_DIR + "brain_edgeslist.txt").c_str(), 0, 1);
  // Here we are converting to an undirected graph for now
  save_graph_to_file(G, walker::GRAPH_DIR + "brain.dat");
}