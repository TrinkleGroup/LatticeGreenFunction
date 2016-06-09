include ../compilervars.mk

all: build2dlgf genKgrid build3dlgf build2dlgf2d 2d-kptmesh

.cc.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ -c $<

LIBMATLIB=-L../matrixlib -lmatlib

build2dlgf.o: ../matrixlib/libmatlib.a build2dlgf.cc
build2dlgf2d.o: ../matrixlib/libmatlib.a build2dlgf2d.cc
build3dlgf.o: ../matrixlib/libmatlib.a build3dlgf.cc
genKgrid.o: ../matrixlib/libmatlib.a genKgrid.cc
2d-kptmesh.o: 2dkgrid/cell.H 2dkgrid/dcomp.H 2dkgrid/dislocation.H 2dkgrid/elastic.H 2dkgrid/gcd.H 2dkgrid/io.H 2dkgrid/kpts.H 2dkgrid/matrix.H 2d-kptmesh.cc

B2DLGFOBJS=build2dlgf.o
B2DLGF2DOBJS=build2dlgf2d.o
B3DLGFOBJS=build3dlgf.o
GKGRIDOBJS=genKgrid.o
KPT2DOBJS=2d-kptmesh.o

build2dlgf: $(B2DLGFOBJS)
	$(CXX) $(CXXFLAGS) $(CXXPROGFLAGS) -o build2dlgf $(B2DLGFOBJS) $(LIBMATLIB) $(LIBS)

build2dlgf2d: $(B2DLGF2DOBJS)
	$(CXX) $(CXXFLAGS) $(CXXPROGFLAGS) -o build2dlgf2d $(B2DLGF2DOBJS) $(LIBMATLIB) $(LIBS)

build3dlgf: $(B3DLGFOBJS)
	$(CXX) $(CXXFLAGS) $(CXXPROGFLAGS) -o build3dlgf $(B3DLGFOBJS) $(LIBMATLIB) $(LIBS)

genKgrid: $(GKGRIDOBJS)
	$(CXX) $(CXXFLAGS) $(CXXPROGFLAGS) -o genKgrid $(GKGRIDOBJS) $(LIBMATLIB) $(LIBS)

2d-kptmesh: $(KPT2DOBJS)
	$(CXX) $(CXXFLAGS) $(CXXPROGFLAGS) -o 2d-kptmesh $(KPT2DOBJS)


.PHONY: clean

clean:
	-rm *.o build2dlgf genKgrid build3dlgf build2dlgf2d 2d-kptmesh
