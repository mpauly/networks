#ifndef WALKER_IO_H
#define WALKER_IO_H

#include "base.h"
#include <fstream>
#include <iomanip>
#include <iostream>

namespace walker {

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
  infile.close();
  return walk;
}

} // namespace walker
#endif