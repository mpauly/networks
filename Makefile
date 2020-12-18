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


clean:
	rm -f *.x
	rm -f graphs/*.dat
	rm -f data/*.dat
	rm -f plots/*.dat
