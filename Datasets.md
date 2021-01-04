## PA Roadnet
This network is a network of roads in Pensylvania from the [SNAP project](http://snap.stanford.edu/data/).

## OSM
The OSM dataset is from the '10 DIMACS challenge, see [this webpage]([)https://www.cc.gatech.edu/dimacs10/archive/streets.shtml). They are graph representations of street data from the OpenStreetmap project. This implies that roads are modelled as a series of individual nodes connected by edges (to get their 2d shape right). We remove these intermediate edges in a post-processing step.

## Brain
The brain dataset is from [neurodata.io](https://neurodata.io/mri/) and the dataset BNU1, with the download here [Download](https://mrneurodata.s3.amazonaws.com/data/BNU1/ndmg_0-0-48/graphs/DS72784/sub-0025864_ses-1_dwi_DS72784.gpickle).
BNU1 refers to Beijing Normal University, scans are of healthy young adults. We convert the gpickle file by taking it apart (there are various issues related to python2 and python3), and writing a new file. 