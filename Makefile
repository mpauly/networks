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

graph_data:
	 wget -O - http://snap.stanford.edu/data/roadNet-PA.txt.gz | gunzip -c > graphs/roadNet-PA.txt
	 wget -O - http://snap.stanford.edu/data/as-skitter.txt.gz | gunzip -c > graphs/as-skitter.txt
	 wget -O tmp_network.zip http://nrvis.com/download/data/bn/bn-fly-drosophila_medulla_1.zip && unzip -j -d graphs/tmp tmp_network.zip && mv graphs/tmp/bn-fly-drosophila_medulla_1.edges graphs/fly-drosophila_edges.txt && rm -r graphs/tmp && rm tmp_network.zip

clean:
	rm -f *.x
	rm -f graphs/*.dat
	rm -f data/*.dat
	rm -f plots/*.dat
