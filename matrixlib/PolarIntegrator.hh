#ifndef _POLAR_INTEGRATOR_H_
#define _POLAR_INTEGRATOR_H_
#include <gsl/gsl_integration.h>
#include "Matrix.hh"
#include "ZMatrix.hh"
#include "DynMat.hh"
#include "CutOffFcn.hh"

typedef struct {
	double Rperp;    // |\vec{R}_{\perp}|
	double kmax; // k_{max}
	int    n;    // n of J_n
} intparam2d;


class PolarIntegrator {
private:
	static const size_t INTWORKSZ = 1024l; /* Numerical Integration Workspace Size */
	static const double TOL_A = 0.0; /* absolute integration tolerance */
	static const double TOL_R = 1e-9; /* relative integration tolerance */
	static const int INT_TYPE = GSL_INTEG_GAUSS61; /* Integration Type */
	DynMat *dynmat;
	Matrix (*gexpan[3])(DynMat &, double[3]);
	double tnorm[3];
	double tmag;
	double mnorm[3];
	double mmag;
	double nnorm[3];
	double nmag;

	static double integrandm2(double x, void *params);
	static double integrandm2cut(double x, void *params);
	static double integrandm1(double x, void *params);
	static double integrandm1cut(double x, void *params);
	static double integrandm0(double x, void *params);
	static double integrandm0cut(double x, void *params);
public:
	PolarIntegrator(DynMat &dynmat, double t[3], double m[3], Matrix (*one_on_k2)(DynMat &, double[3]), Matrix (*i_on_k)(DynMat &, double[3]), Matrix (*k0)(DynMat &, double[3]));
	void calcGnExp(int b, unsigned int nmax, Matrix *gn);
	Matrix calcG_b_R(Matrix *Gn, int b, int nmax, double kmax, double R[3], double V);
	double radialBesselIntegral(double kmax, double R[3], int n, int b);
};

#endif
