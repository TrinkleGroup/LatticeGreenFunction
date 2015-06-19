CC=icc
CXX=icpc
CFLAGS=-O2 -xT -fp-model source
CXXFLAGS=$(CFLAGS)
INCLUDES=-I/usr/local/include -I../matrixlib
LIBS=-L/usr/local/lib -lgsl -lgslcblas -lm
AR=xiar
