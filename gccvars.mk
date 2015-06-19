CC=gcc
CXX=g++
DEBUGFLAG=-g
CFLAGS=-O2 -march=core2
CXXFLAGS=$(CFLAGS)
INCLUDES=-I/sw/include -I../matrixlib
LIBS=-L/sw/lib -lgsl -lgslcblas
AR=ar
