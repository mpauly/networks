## Random Walk on Networks
This repository contains code to implement random walks and the computation of the spectral dimension on various networks utilizing [SNAP](http://snap.stanford.edu/). See [the arxiv](https://arxiv.org/abs/2107.07325) for the corresponding paper.

### Installation
Clone this repository with `git clone --recurse-submodules git@github.com:mpauly/networks.git`. This will not only clone the content of this repository, but also a copy of snap.
To install first run `make snap`. This should build all parts of snap that you need.

After that you should be able to run `make`. This should build all executables. All graph files are stored in `data/graphs/` and `data/networks/` (where networks have a weight for their edges), `data/dimension` contains the spectral dimension extracted from a set of random walks and `data/walks/` contains the random walks as binary files. The directory `plots/` contains the resulting plots. To generate graphs execute `./make_graphs.x`. To run a random walk on one of the graphs run `./random_walk.x 1d_ring_100`, and explore the options for the random walker with `./random_walk.x -?`. There also is a `./queue_walk.x` that is more suitable for use in queuing systems.


#### Work flow
The standard work flow are the three following steps:
  1.) Generating/Downloading a graph via `./make_graph.x`. This will generate graph files in the `data/graphs/` subdirectory.
  2.) Executing a random walk via `./random_walk.x` or `./queue_walk.x` on one of the graph files. See `./random_walk.x -?` for more details on the usage. 
  If you use `./random_walk.x`, the result will be a file containing the spectral dimension as a function of the random walker steps in the subfolder `data/dimension`.
  If you use `./queue_walk.x`, the resulting walks are exported to `data/walks`. To analyze the spectral dimension, the walks can be merged using `./merge_walks_to_data.x`, and then also give a file containing the spectral dimension in the subfolder `data/dimension`.
  3.) Plotting the dimension using `python plots/plot_dimension.py data/dimension/SomeDimFile.dat` or any other plotting utility of your choice. Also see the various python files in `plots/` for more detailed analysis.

#### Artificial graph data
Artifical graphs can be generated (and in some cases analyzed) using the commands `grid.x`, `hyperbolic.x`, `ring.x` and `watts_strogatz.x`. 

#### Data handling
Existing graphs can be checked with `./inspect_graph.x`.
If you perform a random walk with `./queue_walk.x` on multiple cores then you typically end up with a set of binary files in `data/walks/`.
To inspect these walks use `./inspect_walks.x`. To merge the various walks into a single binary file use `./merge_walks_to_data.x`.

#### Downloading real-world graph data
To download some real-world graph data just run `make graph_data`. This downloads a few graph datasets. If you then run the corresponding subcommand for `./make_datagraphs.x`, this generates the corresponding graph files that you can then run the walker on.
For an overview of the corresponding datasets see [the list of datasets](Datasets.md).