/* File: PolarIntegrator.cc
 * 2D lattice Green function analytic integration routine
 */
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "gsl_complex_op.hh"
#include "Matrix.hh"
#include "DynMat.hh"
#include "PolarIntegrator.hh"
#include "CutOffFcn.hh"
#include "VecMath.hh"
#include "SystemDimension.hh"

const size_t PolarIntegrator::INTWORKSZ; /* Numerical Integration Workspace Size */
const double PolarIntegrator::TOL_A; /* absolute integration tolerance */
const double PolarIntegrator::TOL_R; /* relative integration tolerance */
const int PolarIntegrator::INT_TYPE; /* Integration Type */

/* dynmat: Dynamical Matrix
 * t: Dislocation line direction vector
 * m: Dislocation cut direction vector
 * (one_on_k2): function pointer to the 1/k^2 dynamical matrix expansion
 * (i_on_k): function pointer to the i/k dynamical matrix expansion
 * (k0): function pointer to the k^0 dynamical matrix expansion
 */
PolarIntegrator::PolarIntegrator(DynMat &dynmat, double t[3], double *m, Matrix (*one_on_k2)(DynMat &, double *), Matrix (*i_on_k)(DynMat &, double *), Matrix (*k0)(DynMat &, double *)) {
	if(CARTDIM != 3 && CARTDIM != 2) {
		std::cerr << "This only works in 2 and 3 dimensions." << std::endl;
		std::exit(1);
	}
	this->gexpan[0] = one_on_k2;
	this->gexpan[1] = i_on_k;
	this->gexpan[2] = k0;
	this->dynmat = &dynmat;
	mnorm = new double[CARTDIM];
	nnorm = new double[CARTDIM];
	if(CARTDIM == 3) {
		this->tmag = vecmag(t);
		for(unsigned int idx = 0u; idx < 3u; ++idx) {
			this->tnorm[idx] = t[idx]/this->tmag;
		}
	} else {
		this->tmag = 1.0;
		this->tnorm[0] = 0.0;
		this->tnorm[1] = 0.0;
		this->tnorm[2] = 1.0;
	}

	this->mmag = vecmag(m);
	for(unsigned int idx = 0u; idx < CARTDIM; ++idx) {
		this->mnorm[idx] = m[idx]/this->mmag;
	}

	if(CARTDIM == 3) {
		crossprod(t,m,this->nnorm);
		this->nmag = vecmag(this->nnorm);
		this->nnorm[0] /= this->nmag;
		this->nnorm[1] /= this->nmag;
		this->nnorm[2] /= this->nmag;
	} else if(CARTDIM == 2) {
		this->nnorm[0] = -mnorm[1];
		this->nnorm[1] = mnorm[0];
	}

}

PolarIntegrator::~PolarIntegrator() {
	if(mnorm != NULL) {
		delete [] mnorm;
	}
	if(nnorm != NULL) {
		delete [] nnorm;
	}
}

/* Calculates the angular Fourier coefficients for the k^b term
 * b: the Laurent expansion term (-2,-1,0) for the angular fourier expansion
 *    of the Dynamical matrix
 * nmax: maximum term of the Fourier expansion
 * Do a real polar expansion since G^{(b)}(k_phi) are real
 * Thus the output gn Matrix array is packed as:
 * gn[0] = G_0
 * gn[n=1..N/2-1] = Re(G_n)
 * gn[N/2] = G_{N/2}
 * gn[n=N/2+1..N-1] = Im(G_{N-n})
 */
void PolarIntegrator::calcGnExp(int b, unsigned int nmax, Matrix *gn) {
	/* nmax must be a power of 2 for the fft */

	Matrix (*G)(DynMat &, double *);

	switch (b) {
	case -2:
		G = gexpan[0];
		break;
	case -1:
		G = gexpan[1];
		break;
	case 0:
		G = gexpan[2];
		break;
	default:
		exit(-1);
	}

	double Gpoints[CARTDIM*dynmat->getNions()][CARTDIM*dynmat->getNions()][nmax];

	/* phi_m = 2*pi*n/Nmax 
	 * expansion is in the plane \vec{k} \cdot \vec{t} = 0
	 * thus we start from reference direction m and
	 * rotate about it
	 */
	double dphi = 2.0*M_PI/nmax;
	double k[CARTDIM];
	for(int m=0; m < nmax; ++m) {
		double cosmdphi = std::cos(m*dphi);
		double sinmdphi = std::sin(m*dphi);
		for(unsigned int idx = 0u; idx < CARTDIM; ++idx) {
			k[idx] = cosmdphi*mnorm[idx] + sinmdphi*nnorm[idx];
		}
		Matrix gval = (*G)(*dynmat, k);
		/*
		if(b==-1) {
			std::cerr << "Gk_b=-1(k={"<< k[0] << "," << k[1] << "," << k[2] << "})::\n"<<gval<<"\n";
		}
		*/
		for(int i = 0; i<CARTDIM*dynmat->getNions(); ++i) {
			for(int j = 0; j<CARTDIM*dynmat->getNions(); ++j) {
				Gpoints[i][j][m] = gval.val(i,j);
			}
		}
	}
	for(int i = 0; i<CARTDIM*dynmat->getNions(); ++i) {
		for(int j = 0; j<CARTDIM*dynmat->getNions(); ++j) {
			if(gsl_fft_real_radix2_transform(Gpoints[i][j], 1, nmax) != 0) {
				std::cerr << "FFT ERROR!" << std::endl;
				exit(-1);
			}
		}
	}
	/*
	if(b==-1) {
		for(int m=0; m < nmax; ++m) {
			std::cerr << "Gfft_b=-1_m="<<m<<":\n"<<Gpoints[0][CARTDIM][m]<<"\n";
		}
	}
	*/
	/* for b==-2 && b==0, only even n are non-zero */
	if(b==-2 || b==0) {
		for(int q=0; q < nmax/2; ++q) {
			gn[q] = Matrix(CARTDIM*dynmat->getNions(),CARTDIM*dynmat->getNions());
			for(int i = 0; i<CARTDIM*dynmat->getNions(); ++i) {
				for(int j = 0; j<CARTDIM*dynmat->getNions(); ++j) {
					gn[q].val(i,j) = Gpoints[i][j][2*q]/(nmax);
				}
			}
		}
	/* for b==-1, only odd n are non-zero */
	} else {
		for(int q=0; q < nmax/2; ++q) {
			gn[q] = Matrix(CARTDIM*dynmat->getNions(),CARTDIM*dynmat->getNions());
			for(int i = 0; i<CARTDIM*dynmat->getNions(); ++i) {
				for(int j = 0; j<CARTDIM*dynmat->getNions(); ++j) {
					gn[q].val(i,j) = Gpoints[i][j][2*q+1]/(nmax);
				}
			}
		}
	}
}

/* Calculate the real space integrated Green function expansion term
 * for vector R from the polar Fourier expansion pieces Gn for k^b term
 * V is the unit cell volume
 * kmax is the cutoff k value
 * nmax is the number of Fourier expansion pieces
 */
Matrix PolarIntegrator::calcG_b_R(Matrix *Gn, int b, int nmax, double kmax, double *R, double V) {
	Matrix sum(CARTDIM*dynmat->getNions(), CARTDIM*dynmat->getNions());

	double RdotM = vecdot(R,mnorm);
	double RdotN = vecdot(R,nnorm);

	double phiR = atan2(RdotN,RdotM);

	double V_t2pi = V/(tmag*2.0*M_PI);

	if(abs(b)%2 == 0) {
		double coeff0 = V_t2pi*radialBesselIntegral(kmax, R, 0, b);
#ifdef DEBUG
		std::cerr << "coeff0(R={" << R[0] << "," << R[1] << "," << R[2]
		          << "},kmax=" << kmax << ",b=" << b << ",V=" << V << "):\n"
		          << coeff0 << std::endl;
#endif
		sum = coeff0 * Gn[0];
		double coeffN2 = V_t2pi*((nmax/4)%2 == 0?1.0:-1.0)*radialBesselIntegral(kmax, R, nmax/2, b);
		sum += (coeffN2 * cos(phiR*nmax/2)) * Gn[nmax/4];
		for(int q = 1; q<nmax/4; ++q) {
			int n = 2*q;
			double coeff = 2.0*V_t2pi*(q%2 == 0?1.0:-1.0)*radialBesselIntegral(kmax, R, n, b);
			sum += (coeff*cos(phiR*n))*Gn[q];
			sum -= (coeff*sin(phiR*n))*Gn[nmax/2-q];
		}
	} else {
		for(int q = 0; q<nmax/4; ++q) {
			int n = 2*q+1;
			double coeff = 2.0*V_t2pi*(q%2 == 0?1.0:-1.0)*radialBesselIntegral(kmax, R, n, b);
			sum += (coeff*cos(phiR*n))*Gn[q];
			sum -= (coeff*sin(phiR*n))*Gn[nmax/2-q-1];
		}
	}
	return sum;
}

/* This does not converge for n == 0, must treat it as a special case */
double PolarIntegrator::integrandm2(double x, void *params) {
	intparam2d p = *(intparam2d *)params;
	return gsl_sf_bessel_Jn(p.n,x*p.kmax*p.Rperp)/x;
}
double PolarIntegrator::integrandm2cut(double x, void *params) {
	intparam2d p = *(intparam2d *)params;
	return fcutspline(x)*gsl_sf_bessel_Jn(p.n,x*p.kmax*p.Rperp)/x;
}

double PolarIntegrator::integrandm1(double x, void *params) {
	intparam2d p = *(intparam2d *)params;
	return p.kmax*gsl_sf_bessel_Jn(p.n,x*p.kmax*p.Rperp);
}
double PolarIntegrator::integrandm1cut(double x, void *params) {
	intparam2d p = *(intparam2d *)params;
	return p.kmax*fcutspline(x)*gsl_sf_bessel_Jn(p.n,x*p.kmax*p.Rperp);
}

double PolarIntegrator::integrandm0(double x, void *params) {
	intparam2d p = *(intparam2d *)params;
	return p.kmax*p.kmax*x*gsl_sf_bessel_Jn(p.n,x*p.kmax*p.Rperp);
}
double PolarIntegrator::integrandm0cut(double x, void *params) {
	intparam2d p = *(intparam2d *)params;
	return p.kmax*p.kmax*x*fcutspline(x)*gsl_sf_bessel_Jn(p.n,x*p.kmax*p.Rperp);
}


/* Bessel function integrals to handle the cut-off function properly */
double PolarIntegrator::radialBesselIntegral(double kmax, double *R, int nsigned, int b) {
	assert(b == -2 || b == -1 | b == 0);
	double Rperp = std::sqrt(fabs(vecmag2(R) - vecdot(R,tnorm)*vecdot(R,tnorm)));
	int n = abs(nsigned);
	/* special value: Rperp==0 */
	if(Rperp < 1.0e-14) {
		double val = 0.0;
		if(n!=0) {
			return 0.0;
		}
		switch(b) {
		case -2:
			val = std::log(kmax) + M_EULER - M_LN2
			    + (6.0*ALPHA_CUT*ALPHA_CUT*(3.0-ALPHA_CUT)*std::log(ALPHA_CUT)
			    - (1.0-ALPHA_CUT)*(5.0*ALPHA_CUT*ALPHA_CUT - 22.0*ALPHA_CUT + 5.0))
			    / (6.0*(1.0-ALPHA_CUT)*(1.0-ALPHA_CUT)*(1.0-ALPHA_CUT));
			break;
		case -1:
			val = kmax*(ALPHA_CUT+1.0)/2.0;
			break;
		case 0:
			val = kmax*kmax*(0.5 - (1.0-ALPHA_CUT)*(3.0*ALPHA_CUT+7.0)/20.0);
			break;
		default:
			break;
		}
		return val;
	}

	/* Adaptive Gauss-Kronrod routine for the general case: */
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(INTWORKSZ);
	
	int errcode;
	double res;
	double abserr;
	gsl_function fun;
	intparam2d p;
	p.Rperp = Rperp;
	p.kmax = kmax;
	p.n = n;
	fun.params = &p;
	switch(b) {
		case -2: /* Special case for n==0.  Must rewrite integral for convergence*/
			fun.function = &integrandm2;
			break;
		case -1:
			fun.function = &integrandm1;
			break;
		case 0:
			fun.function = &integrandm0;
			break;
		default:
			exit(-1);
	}

	if(n>0 || b!=-2) {
		errcode = gsl_integration_qag(&fun, 0.0, ALPHA_CUT, TOL_A, TOL_R, INTWORKSZ, INT_TYPE, w, &res, &abserr);
#ifdef DEBUG
		std::cerr << "Integration Error Estimate 1: " << abserr << std::endl;
#endif
	} else { /*b == -2 && n == 0 special case*/
		res = std::log(ALPHA_CUT*kmax/2.0) + M_EULER;
		double sumval = 0.0;
		double term = 0.0;
		int q = 0;
		do {
			double fact = gsl_sf_fact(q+1u);
			term = (q%2 == 1 ? 0.5 : -0.5)*gsl_sf_pow_int(ALPHA_CUT*kmax*Rperp/2.0,2*q+2)/(q+1.0)/fact/fact;
			sumval += term;
			++q;
		} while(fabs(term/sumval) > TOL_R);
#ifdef DEBUG
		std::cerr << "Sumval = " << sumval << "\nlog+euler = " << res << "\ntotal = " << res+sumval << std::endl;
#endif
		res += sumval;
	}

	double res2;
	switch(b) {
		case -2:
			fun.function = &integrandm2cut;
			break;
		case -1:
			fun.function = &integrandm1cut;
			break;
		case 0:
			fun.function = &integrandm0cut;
			break;
		default:
			exit(-1);
	}
	errcode = gsl_integration_qag(&fun, ALPHA_CUT, 1.0, TOL_A, TOL_R, INTWORKSZ, INT_TYPE, w, &res2, &abserr);

#ifdef DEBUG
	std::cerr << "Integration Error Estimate 2: " << abserr << std::endl;
#endif
	gsl_integration_workspace_free(w);

	return (res + res2);
}
