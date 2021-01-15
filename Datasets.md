## PA Roadnet
This network is a network of roads in Pennsylvania from the [SNAP project](http://snap.stanford.edu/data/).

## OSM
The OSM dataset is from the '10 DIMACS challenge, see [this webpage](https://www.cc.gatech.edu/dimacs10/archive/streets.shtml). They are graph representations of street data from the OpenStreetMap project. This implies that roads are modelled as a series of individual nodes connected by edges (to get their 2d shape right). We remove these intermediate edges in a post-processing step.

## Drosophila Brain
The large drosophila brain dataset is the one described in [this preprint](https://www.biorxiv.org/content/10.1101/2020.01.21.911859v1). The corresponding download can be found [here](https://www.janelia.org/project-team/flyem/hemibrain) after registration. It needs to be downloaded and unpacked to `graphs/data`. We are using `v1.2` of the dataset for now.

## Brain
The brain dataset is from [neurodata.io](https://neurodata.io/mri/) and the dataset BNU1, with the download here [Download](https://mrneurodata.s3.amazonaws.com/data/BNU1/ndmg_0-0-48/graphs/DS72784/sub-0025864_ses-1_dwi_DS72784.gpickle).
BNU1 refers to Beijing Normal University, scans are of healthy young adults. We convert the gpickle file by taking it apart (there are various issues related to python2 and python3), and writing a new file. Note that the edges here come with a weight - for now we simply ignore this weight when converting to a graph.