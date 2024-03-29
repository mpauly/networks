SHELL := /bin/bash
CC = g++
CXXFLAGS = -std=c++17 -O3 -g -fopenmp
LIBS = -lrt -lsnap -lboost_iostreams
SNAP_DIR = snap
LDLIBS = -lstdc++fs
LIB_DIR = -L ${SNAP_DIR}/snap-core/
INCLUDE = -I ${SNAP_DIR}/glib-core/ -I ${SNAP_DIR}/snap-core/
objects := $(wildcard *.cpp)

all: ${objects:.cpp=.x}

%.x: %.cpp
	${CC} ${INCLUDE} ${LIB_DIR} ${LIBS} ${CXXFLAGS} -o $(@) $< ${SNAP_DIR}/snap-core/Snap.o ${LDLIBS}

snap:
	cd snap && $(MAKE)
	cd snap/snap-core && $(MAKE) lib

graph_data:
	 wget -O - http://snap.stanford.edu/data/roadNet-PA.txt.gz | gunzip -c > data/raw/roadNet-PA.txt
	 wget -O - http://snap.stanford.edu/data/as-skitter.txt.gz | gunzip -c > data/raw/as-skitter.txt
	 wget -O tmp_network.zip http://nrvis.com/download/data/bn/bn-fly-drosophila_medulla_1.zip && unzip -j -d data/tmp tmp_network.zip && mv data/tmp/bn-fly-drosophila_medulla_1.edges data/raw/fly-drosophila_edges.txt && rm -r data/tmp && rm tmp_network.zip
	 wget -O - https://suitesparse-collection-website.herokuapp.com/MM/DIMACS10/europe_osm.tar.gz | tar -xz -C data/raw/
	 wget -O data/raw/sub-0025864_ses-1_dwi_DS72784.gpickle https://mrneurodata.s3.amazonaws.com/data/BNU1/ndmg_0-0-48/graphs/DS72784/sub-0025864_ses-1_dwi_DS72784.gpickle
	 # the following does not work - one needs to authenticate in order to download that dataset
	 # wget -O - https://storage.cloud.google.com/hemibrain/v1.2/exported-traced-adjacencies-v1.2.tar.gz | tar -xz -C data/raw/
	 wget -O data/raw/metabolism.kgml http://rest.kegg.jp/get/rn01100/kgml

diffusion_test:
	for i in {1..9}; do ./random_walk.x -l 200 -d 0.$$i -s 5050 -o data/2dtest/d0$$i.dat graphs/2d_lattice_100.dat; done
	python plots/plot_2dtest.py

oneDimPlateau:
	./random_walk.x -s 10 -l 1000 -d 0.25 graphs/1d_ring_26.dat

clean:
	rm -f *.x
	rm -f graphs/*.dat
	rm -f plots/*.png
