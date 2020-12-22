SHELL := /bin/bash
CC = g++
CXXFLAGS = -std=c++17 -O3 -g -fopenmp
LIBS = -lrt -lsnap
SNAP_DIR = snap
LIB_DIR = -L ${SNAP_DIR}/snap-core/
INCLUDE = -I ${SNAP_DIR}/glib-core/ -I ${SNAP_DIR}/snap-core/
objects := $(wildcard *.cpp)

all: ${objects:.cpp=.x}

%.x: %.cpp
	${CC} ${INCLUDE} ${LIB_DIR} ${LIBS} ${CXXFLAGS} -o $(@) $< ${SNAP_DIR}/snap-core/Snap.o

snap:
	cd snap && $(MAKE)
	cd snap/snap-core && $(MAKE) lib

graph_data:
	 wget -O - http://snap.stanford.edu/data/roadNet-PA.txt.gz | gunzip -c > graphs/data/roadNet-PA.txt
	 wget -O - http://snap.stanford.edu/data/as-skitter.txt.gz | gunzip -c > graphs/data/as-skitter.txt
	 wget -O tmp_network.zip http://nrvis.com/download/data/bn/bn-fly-drosophila_medulla_1.zip && unzip -j -d graphs/tmp tmp_network.zip && mv graphs/tmp/bn-fly-drosophila_medulla_1.edges graphs/data/fly-drosophila_edges.txt && rm -r graphs/tmp && rm tmp_network.zip
	 wget -O - https://suitesparse-collection-website.herokuapp.com/MM/DIMACS10/europe_osm.tar.gz | tar -xz -C graphs/data/

diffusion_test:
	for i in {1..9}; do ./random_walk.x -l 200 -d 0.$$i -s 5050 -o data/2dtest/d0$$i.dat graphs/2d_lattice_100.dat; done
	python plots/plot_2dtest.py

clean:
	rm -f *.x
	rm -f graphs/*.dat
	rm -f plots/*.png
