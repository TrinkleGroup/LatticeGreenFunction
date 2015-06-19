#ifndef _SPHERICAL_INTEGRATOR_H_
#define _SPHERICAL_INTEGRATOR_H_

#include <gsl/gsl_complex.h>
#include <gsl/gsl_integration.h>
#include "DynMat.hh"

typedef struct {
	double R;    // |\vec{R}|
	double kmax; // k_{max}
	int    l;    // l of j_l
} intparam;

class SphericalIntegrator {
private:
	DynMat *dynmat;
	Matrix (*gexpan[3])(DynMat &, double[3]);

	static const size_t INTWORKSZ = 1024l; /* Numerical Integration Workspace Size */
	static const double TOL_A = 0.0; /* absolute integration tolerance */
	static const double TOL_R = 1e-10; /* relative integration tolerance */
	static const double TOL_ROOT = 3e-13; /* Legendre Root tolerance */
	static const int INT_TYPE = GSL_INTEG_GAUSS61; /* Integration Type */
	/* Gauss-Legendre integration points from Numerical Recipes */
	static void gauleg(double x[], double w[], int n);
	static double calcYlmR(int l, int m, double Rvec[3]);
	static double integrandm2(double x, void *params);
	static double integrandm2cut(double x, void *params);
	static double integrandm1(double x, void *params);
	static double integrandm1cut(double x, void *params);
	static double integrandm0(double x, void *params);
	static double integrandm0cut(double x, void *params);
	static Matrix calcG_b_R_even(Matrix *Glm, int b, unsigned int lmax, double kmax, double R[3], double V, unsigned int natoms);
	static Matrix calcG_b_R_odd(Matrix *Glm, int b, unsigned int lmax, double kmax, double R[3], double V, unsigned int natoms);
public:
	SphericalIntegrator(DynMat &dynmat, Matrix (*one_on_k2)(DynMat &, double[3]), Matrix (*i_on_k)(DynMat &, double[3]), Matrix (*k0)(DynMat &, double[3]));
	Matrix one_on_k2_angle(double theta, double phi);
	Matrix i_on_k_angle(double theta, double phi);
	Matrix k0_angle(double theta, double phi);
	void calcGlmExp(int b, unsigned int lmax, Matrix *glm);
	Matrix calcGlm(int b, unsigned int l, int m, unsigned int lmax);
	static Matrix calcG_b_R(Matrix *Glm, int b, unsigned int lmax, double kmax, double R[3], double V, unsigned int natoms);
	static double radialBesselIntegral(double kmax, double R, int l, int b);
};

#endif
