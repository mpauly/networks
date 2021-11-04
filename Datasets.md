## PA Roadnet
This network is a network of roads in Pennsylvania from the [SNAP project](http://snap.stanford.edu/data/).

## OSM
The OSM dataset is from the '10 DIMACS challenge, see [this webpage](https://www.cc.gatech.edu/dimacs10/archive/streets.shtml). They are graph representations of street data from the OpenStreetMap project. This implies that roads are modelled as a series of individual nodes connected by edges (to get their 2d shape right). We remove these intermediate edges in a post-processing step. This is the difference between the `osm` and `osm_reduced`.

## Drosophila Brain
The large drosophila brain dataset is the one described in [this preprint](https://www.biorxiv.org/content/10.1101/2020.01.21.911859v1). The corresponding download can be found [here](https://www.janelia.org/project-team/flyem/hemibrain) after registration. It needs to be downloaded and unpacked to `data/raw/`. We are using `v1.2` of the dataset for now.

## Brain - unused
The brain dataset is from [neurodata.io](https://neurodata.io/mri/) and the dataset BNU1, with the download here [Download](https://mrneurodata.s3.amazonaws.com/data/BNU1/ndmg_0-0-48/graphs/DS72784/sub-0025864_ses-1_dwi_DS72784.gpickle).
BNU1 refers to Beijing Normal University, scans are of healthy young adults. We convert the gpickle file by taking it apart (there are various issues related to python2 and python3), and writing a new file. Note that the edges here come with a weight — for now we simply ignore this weight when converting to a graph.

## Rat brain
The rat network is obtained by taking the correlation matrix obtained [here](https://direct.mit.edu/netn/article/3/1/217/2194/High-resolution-data-driven-model-of-the-mouse), that measures the correlations between various voxels (three-dimensional volumes) in a mouse brain. We model each voxel by a graph node. We then introduce a weighted edge between two nodes if they are linked by a correlation that is larger than the cutoff $10^-3$. See also [here](https://github.com/AllenInstitute/mouse_connectivity_models).

## Brain small (300 regions) - unused
The network of human brain areas is derived from a network matrix by the Human Connectome Project. To obtain such a matrix N=300 distinct regions are distinguished in the brain. For each participant a time series of measurements is converted into a matrix that measures connectivity between these areas. The resulting matrices are averaged over 1003 participants. We use this group-averaged full correlation network matrix and build a weighted, undirected graph from it, by assigning a node to each region, and weighting the edge between two regions by the absolute value of the corresponding correlation. The corresponding file is `3T_HCP1200_MSMAll_d300_ts2/Mnet1.pconn.nii`, which can be converted to a `.tsv` with command-line tools from the HCP. 
The resulting graph is fairly small (approximately 300 nodes) and hence is not ideal for our analysis.

## Metabolism
The metabolism dataset is a reaction network from [KEGG](https://www.genome.jp/kegg/pathway.html). We are just considering reactions, and only take into account the largest connected component of the resulting graph. The graph that we are considering is [this one](https://www.kegg.jp/kegg-bin/show_pathway?rn01100). See the corresponding script in `utils/convert_metabolism_edges.py`.

## Internet
Files can be downloaded from CAIDA after [registration](https://www.caida.org/data/request_user_info_forms/ark.xml).
Then just extract the AS ids like so: 
`cat cycle-aslinks.l7.t1.c008218.20200228.txt | grep "^D" | awk '{print $2,$3}' - > internet_caida_as.dat`.