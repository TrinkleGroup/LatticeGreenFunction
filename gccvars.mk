CC=gcc
CXX=g++
DEBUGFLAG=-g
CFLAGS=-O3 -ftree-vectorize -march=native -flto -floop-parallelize-all -fgraphite-identity -floop-block -floop-strip-mine -floop-interchange
#CFLAGS=-O0 -g -DDEBUG
CXXFLAGS=$(CFLAGS)
CXXPROGFLAGS=-fuse-linker-plugin -fwhole-program
#CXXPROGFLAGS=
INCLUDES=-I/sw/include -I../matrixlib
LIBS=-L/sw/lib -lgsl -lgslcblas
AR=gcc-ar
RANLIB=gcc-ranlib
