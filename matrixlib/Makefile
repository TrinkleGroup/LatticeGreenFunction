include ../compilervars.mk

all: libmatlib.a

.cc.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ -c $<

SystemDimension.o: SystemDimension.hh
Matrix.o: Matrix.hh
DynMat.o: Matrix.o DynMat.hh SystemDimension.hh
ZMatrix.o: ZMatrix.hh Matrix.hh Matrix.o
allmats.o: Matrix.o DynMat.o allmats.cc SystemDimension.hh
gsl_complex_op.o: gsl_complex_op.hh
SphericalIntegrator.o: SphericalIntegrator.hh CutOffFcn.hh Matrix.o DynMat.o ZMatrix.o
PolarIntegrator.o: PolarIntegrator.hh CutOffFcn.hh Matrix.o DynMat.o ZMatrix.o SystemDimension.hh
SemiCont.o: Matrix.o ZMatrix.o SemiCont.hh gsl_complex_op.o SystemDimension.hh
UnitCell.o: Matrix.o UnitCell.hh SystemDimension.hh

OBJS=Matrix.o DynMat.o allmats.o ZMatrix.o SphericalIntegrator.o SemiCont.o gsl_complex_op.o PolarIntegrator.o UnitCell.o SystemDimension.o

.PHONY: clean

libmatlib.a: $(OBJS)
	-rm libmatlib.a
	$(AR) vq libmatlib.a $(OBJS)
	$(RANLIB) libmatlib.a

clean:
	-rm *.o libmatlib.a
