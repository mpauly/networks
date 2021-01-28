#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <vector>

#define WALKER_DEBUG false

struct RandomWalk {
  int start_node;
  int sigma;
  double diffusion_constant;
  std::vector<double> return_probability;
  std::vector<double> dimension;
  THash<TInt, double> lvl_probabilities;
};

RandomWalk setupRandomWalk(const int &start_node, const double diffusion_constant) {
  RandomWalk walk;
  walk.return_probability.resize(1);
  walk.dimension.resize(1);
  walk.start_node = start_node;
  walk.sigma = 0;
  walk.diffusion_constant = diffusion_constant;
  walk.return_probability[0] = 1.0;
  walk.lvl_probabilities.AddDat(start_node, 1.0);
  return walk;
}

// legacy wrapper
template <class PGraph>
std::vector<double> spectralDimensionAtNode(const PGraph &Graph, const int &start_node, const int &max_depth,
                                            std::function<void(int)> progress_monitor,
                                            const double diffusion_constant) {
  // setup the walk
  RandomWalk walk = setupRandomWalk(start_node, diffusion_constant);
  progressRandomWalk(Graph, walk, max_depth, progress_monitor);
  return walk.dimension;
}

// The implementation roughly corresponds to a breadth-first search tree.
// Inspired by SNAP's DoBfs
template <class PGraph>
void progressRandomWalk(const PGraph &Graph, RandomWalk &walk, int nr_of_steps,
                        std::function<void(int)> progress_monitor) {

  walk.return_probability.resize(walk.sigma + nr_of_steps + 1);
  walk.dimension.resize(walk.sigma + nr_of_steps + 1);

  TSnapQueue<int> queue;
  THash<TInt, double> last_lvl_probabilities = walk.lvl_probabilities;
  int childnode;
  double probability_derivative, node_probability;

  // setup the queue
  queue.Clr(false);
  for (auto it = last_lvl_probabilities.BegI(); it < last_lvl_probabilities.EndI(); it++) {
    IAssert(Graph->IsNode(it.GetKey()));
    queue.Push(it.GetKey());
  }

  for (int sigma = walk.sigma + 1; sigma < walk.sigma + nr_of_steps + 1; sigma++) {
    walk.lvl_probabilities.Clr();

    if (WALKER_DEBUG)
      printf("\n == Sigma %d ==\n", sigma);
    progress_monitor(sigma);
    // deal with all the nodes that are in the queue at the moment
    // this relies on the loop initialization begin run exactly once, as the queue grows while the loop is executed
    for (int node = queue.Len(); node > 0; node--) {
      // take a node out of the queue
      const int nodeId = queue.Top();
      if (WALKER_DEBUG)
        printf("   Node %d   ", nodeId);
      queue.Pop();
      const typename PGraph::TObj::TNodeI NodeI = Graph->GetNI(nodeId);

      const int nodeDegree = NodeI.GetOutDeg();
      // the diffusion constant term at the same node
      if (walk.lvl_probabilities.IsKey(nodeId)) {
        node_probability = walk.lvl_probabilities.GetDat(nodeId);
      } else {
        node_probability = 0;
        queue.Push(nodeId);
      }
      node_probability += (1.0 - walk.diffusion_constant) * last_lvl_probabilities.GetDat(nodeId);
      walk.lvl_probabilities.AddDat(nodeId, node_probability);
      // loop over child nodes
      for (childnode = 0; childnode < nodeDegree; childnode++) {
        const int childId = NodeI.GetOutNId(childnode);
        // did we already touch this node? If yes just update existing entry
        // if no then create a new entry and add the node to the queue
        if (walk.lvl_probabilities.IsKey(childId)) {
          node_probability = walk.lvl_probabilities.GetDat(childId);
        } else {
          node_probability = 0.0;
          queue.Push(childId);
        }
        node_probability += walk.diffusion_constant * last_lvl_probabilities.GetDat(nodeId) / ((double)nodeDegree);
        if (WALKER_DEBUG)
          printf(" %d : %f ", childId, node_probability);
        // AddDat also overwrites an existing entry
        walk.lvl_probabilities.AddDat(childId, node_probability);
      }

      if (walk.lvl_probabilities.IsKey(walk.start_node)) {
        walk.return_probability[sigma] = walk.lvl_probabilities.GetDat(walk.start_node);
      } else {
        walk.return_probability[sigma] = 0.0;
      }
    }
    last_lvl_probabilities = walk.lvl_probabilities;

    // finally extract the dimension
    if (sigma > 1) {
      // This is the probability derivative and dimension at position sigma -1
      probability_derivative = 0.5 * (-walk.return_probability[sigma - 2] + walk.return_probability[sigma]);
      // this is the spectral dimension at sigma - 1
      walk.dimension[sigma - 1] =
          -2.0 * ((double)sigma - 1) / walk.return_probability[sigma - 1] * probability_derivative;
      if (WALKER_DEBUG)
        printf("\n d=%f", walk.dimension[sigma - 1]);
    }
  }
  walk.sigma += nr_of_steps;
}

// Convenience function without progress monitor
template <class PGraph> void progressRandomWalk(const PGraph &Graph, RandomWalk &walk, int nr_of_steps) {
  std::function<void(int)> progress_monitor = [](int sigma) {};
  progressRandomWalk(Graph, walk, nr_of_steps, progress_monitor);
}

void exportRandomWalkToFile(RandomWalk walk, std::string filename, std::string comment) {
  std::ofstream outfile(filename, std::ios::out);
  outfile << comment << std::endl;
  outfile << "start_node\t" << walk.start_node << "\nsigma\t" << walk.sigma << "\ndiffusion_constant\t"
          << walk.diffusion_constant << std::endl;
  outfile << "\n\ndimension" << std::endl;
  for (int i = 0; i < walk.sigma + 1; i++) {
    outfile << walk.dimension[i] << "\t";
  }
  outfile << std::fixed << std::setprecision(20);
  outfile << std::scientific;
  outfile << "\n\nreturn_probability" << std::endl;
  for (int i = 0; i < walk.sigma + 1; i++) {
    outfile << walk.return_probability[i] << "\t";
  }
  outfile << "\n\ndistribution\n";
  outfile << std::setprecision(50);
  for (auto it = walk.lvl_probabilities.BegI(); it < walk.lvl_probabilities.EndI(); it++) {
    outfile << it.GetKey() << "\t" << it.GetDat() << "\n";
  }
  outfile.close();
}

void exportRandomWalkToBinaryFile(RandomWalk walk, std::string filename) {
  std::ofstream outfile(filename, std::ios::out | std::ios::binary);
  outfile.write(reinterpret_cast<const char *>(&walk.start_node), sizeof(int));
  outfile.write(reinterpret_cast<const char *>(&walk.sigma), sizeof(int));
  outfile.write(reinterpret_cast<const char *>(&walk.diffusion_constant), sizeof(double));

  for (int i = 0; i < walk.sigma + 1; i++) {
    outfile.write(reinterpret_cast<const char *>(&walk.dimension[i]), sizeof(double));
  }
  for (int i = 0; i < walk.sigma + 1; i++) {
    outfile.write(reinterpret_cast<const char *>(&walk.return_probability[i]), sizeof(double));
  }
  int length = walk.lvl_probabilities.Len();
  outfile.write(reinterpret_cast<const char *>(&length), sizeof(int));
  for (auto it = walk.lvl_probabilities.BegI(); it < walk.lvl_probabilities.EndI(); it++) {
    outfile.write(reinterpret_cast<const char *>(&it.GetKey()), sizeof(int));
    outfile.write(reinterpret_cast<const char *>(&it.GetDat()), sizeof(double));
  }
  outfile.close();
}

RandomWalk importRandomWalkFromFile(std::string filename) {
  std::ifstream infile(filename, std::ios::in);
  RandomWalk walk;
  std::string line;
  std::string key;
  std::string value;

  std::getline(infile, line);

  for (int i = 0; i < 3; i++) {
    std::getline(infile, line);
    std::stringstream ss(line);
    std::getline(ss, key, '\t');
    std::getline(ss, value, '\t');
    if (key == "start_node")
      walk.start_node = std::stoi(value);
    if (key == "diffusion_constant")
      walk.diffusion_constant = std::stod(value);
    if (key == "sigma")
      walk.sigma = std::stoi(value);
  }
  walk.return_probability.resize(walk.sigma + 1);
  walk.dimension.resize(walk.sigma + 1);
  {
    while (line != "dimension")
      std::getline(infile, line);

    std::getline(infile, line);
    std::stringstream ss(line);
    int i = 0;
    while (std::getline(ss, value, '\t')) {
      walk.dimension[i] = std::stod(value);
      i++;
    }
    if (i != walk.sigma + 1) {
      throw std::invalid_argument("Inconsistent input file");
    }
  }
  {
    while (line != "return_probability")
      std::getline(infile, line);

    std::getline(infile, line);
    std::stringstream ss(line);
    int i = 0;
    while (std::getline(ss, value, '\t')) {
      walk.return_probability[i] = std::stod(value);
      i++;
    }
    if (i != walk.sigma + 1) {
      throw std::invalid_argument("Inconsistent input file");
    }
  }
  while (line != "distribution")
    std::getline(infile, line);

  while (std::getline(infile, line)) {
    std::stringstream ss(line);
    std::getline(ss, key, '\t');
    std::getline(ss, value, '\t');
    walk.lvl_probabilities.AddDat(std::stoi(key), std::stod(value));
  }

  infile.close();
  return walk;
}

RandomWalk importRandomWalkFromBinaryFile(std::string filename) {
  std::ifstream infile(filename, std::ios::in | std::ios::binary);
  RandomWalk walk;

  infile.read(reinterpret_cast<char *>(&(walk.start_node)), sizeof(int));
  infile.read(reinterpret_cast<char *>(&(walk.sigma)), sizeof(int));
  infile.read(reinterpret_cast<char *>(&(walk.diffusion_constant)), sizeof(double));

  walk.return_probability.resize(walk.sigma + 1);
  walk.dimension.resize(walk.sigma + 1);

  double val_d;

  for (int i = 0; i < walk.sigma + 1; i++) {
    infile.read(reinterpret_cast<char *>(&(walk.dimension[i])), sizeof(double));
  }
  for (int i = 0; i < walk.sigma + 1; i++) {
    infile.read(reinterpret_cast<char *>(&(walk.return_probability[i])), sizeof(double));
  }
  int length;
  infile.read(reinterpret_cast<char *>(&length), sizeof(int));
  int key;
  double value;
  for (int i = 0; i < length; i++) {
    infile.read(reinterpret_cast<char *>(&key), sizeof(int));
    infile.read(reinterpret_cast<char *>(&value), sizeof(double));
    walk.lvl_probabilities.AddDat(key, value);
  }

  return walk;
}