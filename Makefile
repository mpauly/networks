CC = g++
CXXFLAGS = -std=c++11 -O3 -g -fopenmp
LIBS = -lrt -lsnap
SNAP_DIR = ../Snap
LIB_DIR = -L ${SNAP_DIR}/snap-core/
INCLUDE = -I ${SNAP_DIR}/glib-core/ -I ${SNAP_DIR}/snap-core/

make_testgraphs:
	${CC} ${INCLUDE} ${LIB_DIR} ${LIBS} ${CXXFLAGS} -o make_testgraphs.x make_testgraphs.cpp ${SNAP_DIR}/snap-core/Snap.o

clean:
	rm -f *.x
	rm -f graphs/*.dat
