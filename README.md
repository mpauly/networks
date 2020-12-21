## Random Walk on Networks
This repository contains code to implement random walks and the computation of the spectral dimension on various networks utilizing [SNAP](http://snap.stanford.edu/).

### Installation
Clone this repository with `git clone --recurse-submodules git@github.com:mpauly/networks.git`. This will not only clone the content of this repository, but also a copy of snap.
To install first run `make snap`. This should build all parts of snap that you need.

After that you should be able to run `make`. This should build all executables. All graph files are stored in `graphs/`, `data/` contains intermediate files such as e.g. the spectral dimension extracted from a set of random walks. `plots/` then contains the resulting plots. To generate a few example graphs just execute `./make_testgraphs.x`. To run a random walk on one of the graphs run `./random_walk.x graphs/1d_ring_100.dat`, and explore the options for the random walker with `./random_walk.x -?`.

### Downloading graph data
To download some real-world graph data just run `make graph_data`. This downloads a few graph datasets. If you then run `./make_datagraphs.x`, this generates the corresponding graph files that you can then run the walker on.
