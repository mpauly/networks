CC = g++
CXXFLAGS = -O3 -g -fopenmp
LIBS = -lrt -lsnap
SNAP_DIR = ../Snap
LIB_DIR = -L ${SNAP_DIR}/snap-core/
INCLUDE = -I ${SNAP_DIR}/glib-core/ -I ${SNAP_DIR}/snap-core/

testgraphs:
	${CC} ${INCLUDE} ${LIB_DIR} ${LIBS} ${CXXFLAGS} -o testgraphs.x testgraphs.cpp ${SNAP_DIR}/snap-core/Snap.o

clean:
	rm -f *.x
