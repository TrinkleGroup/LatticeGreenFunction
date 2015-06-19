include ../compilervars.mk

all: libmatlib.a

.cc.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ -c $<

Matrix.o: Matrix.hh
DynMat.o: Matrix.o DynMat.hh
ZMatrix.o: ZMatrix.hh Matrix.hh Matrix.o
allmats.o: Matrix.o DynMat.o allmats.cc
gsl_complex_op.o: gsl_complex_op.hh
SphericalIntegrator.o: SphericalIntegrator.hh CutOffFcn.hh Matrix.o DynMat.o ZMatrix.o
PolarIntegrator.o: PolarIntegrator.hh CutOffFcn.hh Matrix.o DynMat.o ZMatrix.o
SemiCont.o: Matrix.o ZMatrix.o SemiCont.hh gsl_complex_op.o
UnitCell.o: Matrix.o UnitCell.hh

OBJS=Matrix.o DynMat.o allmats.o ZMatrix.o SphericalIntegrator.o SemiCont.o gsl_complex_op.o PolarIntegrator.o UnitCell.o

.PHONY: clean

libmatlib.a: $(OBJS)
	-rm libmatlib.a
	$(AR) vq libmatlib.a $(OBJS)
	ranlib libmatlib.a

clean:
	-rm *.o libmatlib.a
