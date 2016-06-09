#include <cmath>
#include <cassert>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
using std::size_t;
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_fft_real.h>
#include "SphericalIntegrator.hh"
#include "CutOffFcn.hh"
#include "VecMath.hh"

//#define DEBUG 1

#include <iostream>

const size_t SphericalIntegrator::INTWORKSZ; /* Numerical Integration Workspace Size */
const double SphericalIntegrator::TOL_A; /* absolute integration tolerance */
const double SphericalIntegrator::TOL_R; /* relative integration tolerance */
const double SphericalIntegrator::TOL_ROOT; /* Legendre Root tolerance */
const int SphericalIntegrator::INT_TYPE; /* Integration Type */


void SphericalIntegrator::gauleg(double x[], double w[], int n)
{
	gsl_integration_glfixed_table *gaulegtab = gsl_integration_glfixed_table_alloc(n);
	for(int i=0; i<n; ++i) {
		x[i] = gaulegtab->x[i];
		w[i] = gaulegtab->w[i];
	}
	gsl_integration_glfixed_table_free(gaulegtab);
}


/* radix-2 FFT for m*phi integrations, 2*lmax + 2 must be a power of 2
 * -> lmax = 2^n - 1 for some integer n
 */
void SphericalIntegrator::calcGlmExp(int b, unsigned int lmax, Matrix *glm)
{
	int N = lmax+1;  //Theta integration: up to order 2*N - 1
	                     //Poly of order 2L -> L + 1 points
	int M = 2*(lmax+1);     //phi integration needs 2L+1 points
	                        //even number of M for nice cancelation
	unsigned int natoms = dynmat->getNions();
	double w[N];
	double x[N];
	Matrix Gtmp(3u*natoms,3u*natoms);
	/* Doing an FFT in m phi for the phi integration,
	 * so we need an array of doubles rather than an array of Matrix classes.
	 * We'll convert to Matrix[] after integrating.
	 */
	double Gpoints[3u*natoms][3u*natoms][N][M];
	/*
	Matrix *glm_ret = new Matrix[lmax + ((lmax*(lmax+1u))>>1u) + 1u];
	*/
	Matrix (SphericalIntegrator::*G)(double, double);

	switch (b) {
	case -2:
		G = &SphericalIntegrator::one_on_k2_angle;
		break;
	case -1:
		G = &SphericalIntegrator::i_on_k_angle;
		break;
	case 0:
		G = &SphericalIntegrator::k0_angle;
		break;
	default:
		exit(-1);
	}
	gauleg(x, w, N);

	for(int i = 0; i<N; ++i) {
		double theta = acos(x[i]);
		for(int j = 0; j < M; ++j) {
			double phi = 2.0*M_PI*((double)j)/(double)M;
			Gtmp = (this->*G)(theta,phi);
			for(int p=0; p<3u*natoms; ++p) {
				for(int q=0; q<3u*natoms; ++q) {
					Gpoints[p][q][i][j] = Gtmp.val(p,q);
				}
			}
		}
	}

	for(int p=0; p<3u*natoms; ++p) {
		for(int q=0; q<3u*natoms; ++q) {
			for(int i = 0; i<N; ++i) {
				if(gsl_fft_real_radix2_transform(Gpoints[p][q][i], 1, M) != 0) {
					std::cerr << "FFT ERROR!" << std::endl;
					exit(-1);
				}
			}
		}
	}

	/* ignore the m=lmax+1 term, we really only need 2l+1 points
	 * for integration, but we'd like an even number of points
	 * for symmetry, plus radix-2 fft is faster and has better
	 * memory localization.  possible to convert to mixed radix
	 * FFT if necessary.
	 * M/2 = lmax+1 
	 * Gpoints[p][q][i][j = m = 0..lmax+1] = int dphi g_{pq}(theta_i, phi)cos(m phi)
	 * Gpoints[p][q][i][j = M - m = lmax+2 .. 2*lmax+1,m=1..lmax] = int dphi g_{pq}(theta_i, phi)sin(m phi)
	 */

	/* allocate zero matrixes for return values */
	for(int i = 0; i<=(lmax + ((lmax*(lmax+1u))>>1u)); ++i) {
		glm[i] = Matrix(3u*natoms,3u*natoms);
	}



	double *plm = new double[lmax+1];
	/* m = 0: special case for real spherical harmonics */
	for(int i = 0; i<N; ++i) {
		if(int err = gsl_sf_legendre_sphPlm_array(lmax, 0, x[i], plm) != 0) {
			std::cerr << "GSL error calculating associated legendre polynomials: err = " << err << " lmax = " << lmax << " m = " << 0 << std::endl;
			exit(-1);
		}
		/* only need even ls for b even, and odd ls for b odd */
		/* must be careful with integer division and boundaries in
		 * this case */
		for(int n = 0; n <= (lmax-(abs(b)%2))/2; ++n) {
			int l = 2*n + (abs(b)%2);
			double weight = 2.0*M_PI/M * w[i] * plm[l];
			for(int p=0;p<3*natoms;++p) {
				for(int q=0; q<3*natoms; ++q) {
					glm[(l*(l+1u))>>1u].val(p,q) += weight*Gpoints[p][q][i][0];
				}
			}
		}
	}

	/* m < 0 and m > 0 */
	for(int i = 0; i<N; ++i) {
		for(int m = 1; m<=lmax; ++m) {
			if(int err = gsl_sf_legendre_sphPlm_array(lmax, m, x[i], plm) != 0) {
				std::cerr << "GSL error calculating associated legendre polynomials: err = " << err << " lmax = " << lmax << " m = " << m << std::endl;
				exit(-1);
			}
			/* only need even ls for b even, and odd ls for b odd */
			/* must be careful with integer division and boundaries in
			 * this case */
			for(int n=(m+1-(abs(b)%2))/2; n<=(lmax-(abs(b)%2))/2; ++n) {
				int l = 2*n + (abs(b)%2);
				double weight = M_SQRT2*2.0*M_PI/M * w[i] * plm[l-m];
				for(int p=0;p<3*natoms;++p) {
					for(int q=0; q<3*natoms; ++q) {
						glm[m+((l*(l+1u))>>1u)].val(p,q) += weight*Gpoints[p][q][i][m];
						glm[-m+((l*(l+1u))>>1u)].val(p,q) -= weight*Gpoints[p][q][i][M-m];
					}
				}
			}
		}

	}
	free(plm);
}

Matrix SphericalIntegrator::calcGlm(int b, unsigned int l, int m, unsigned int lmax) {
	int N = lmax+1;  //Theta integration: up to order 2*N - 1
	                     //Poly of order 2L -> L + 1 points
	int M = 2*(lmax+1);     //phi integration needs 2L+1 points
	                        //even number of M for nice cancelation
	unsigned int natoms = dynmat->getNions();
	double w[N];
	double x[N];
	Matrix glm(3u*natoms, 3u*natoms);
	Matrix (SphericalIntegrator::*G)(double, double);
	switch (b) {
	case -2:
		G = &SphericalIntegrator::one_on_k2_angle;
		break;
	case -1:
		G = &SphericalIntegrator::i_on_k_angle;
		break;
	case 0:
		G = &SphericalIntegrator::k0_angle;
		break;
	default:
		exit(-1);
	}

	gauleg(x, w, N);
#ifdef DEBUG
	std::cerr << "x points : weights\n";
	for(int j=0; j<N; ++j) {
		std::cerr << x[j] << " : " << w[j] << "\n";
	}
	std::cerr.flush();
#endif

/* Use the real basis spherical harmonics:
 * Ylm = Yl0  m = 0
 *     = sqrt(2)Nlm Plm cos m\phi  m>0
 *     = sqrt(2)Nl-m Pl-m sin -m\phi  m<0
 */
	if(m == 0) {
		/* group + and - k together for numerical cancellation */
		for(int j=0; j<N/2; ++j) {
			double weight = 2.0*M_PI/M * w[j] * gsl_sf_legendre_sphPlm(l,0,x[j]);
			double weight2 = 2.0*M_PI/M * w[N-1-j] * gsl_sf_legendre_sphPlm(l,0,x[N-1-j]);
			double theta_j = acos(x[j]);
			double theta_j2 = acos(x[N-1-j]);
			for(int i=0; i<M/2; ++i) {
				double phi_i = 2.0*M_PI*((double)i)/(double)M;
				double phi_i2 = 2.0*M_PI*((double)(M/2+i))/(double)M;
				glm += (weight*(this->*G)(theta_j, phi_i) + weight2*(this->*G)(theta_j2,phi_i2))
				     + (weight2*(this->*G)(theta_j2, phi_i) + weight*(this->*G)(theta_j,phi_i2));
			}
		}
		if(N%2==1) {
			double weight = 2.0*M_PI/M * w[N/2] * gsl_sf_legendre_sphPlm(l,0,x[N/2]);
			double theta_j = acos(x[N/2]);
			for(int i=0; i<M/2; ++i) {
				double phi_i = 2.0*M_PI*((double)i)/(double)M;
				double phi_i2 = 2.0*M_PI*((double)(M/2+i))/(double)M;
				glm += weight*((this->*G)(theta_j, phi_i) + (this->*G)(theta_j,phi_i2));
			}
		}
	} else if(m < 0) {
		for(int j=0; j<N/2; ++j) {
			double weight = 2.0*M_PI/M * w[j] * gsl_sf_legendre_sphPlm(l,-m,x[j]) * M_SQRT2;
			double theta_j = acos(x[j]);
			double weight2 = 2.0*M_PI/M * w[N-1-j] * gsl_sf_legendre_sphPlm(l,-m,x[N-1-j]) * M_SQRT2;
			double theta_j2 = acos(x[N-1-j]);
			for(int i=0; i<M/2; ++i) {
				double phi_i = 2.0*M_PI*((double)i)/(double)M;
				double phi_i2 = 2.0*M_PI*((double)(M/2+i))/(double)M;
				glm += (sin(-m*phi_i)*weight)*(this->*G)(theta_j, phi_i) + (sin(-m*phi_i2)*weight2)*(this->*G)(theta_j2, phi_i2)
			             + (sin(-m*phi_i)*weight2)*(this->*G)(theta_j2, phi_i) + (sin(-m*phi_i2)*weight)*(this->*G)(theta_j, phi_i2);
			}
		}
		if(N%2==1) {
			double weight = 2.0*M_PI/M * w[N/2] * gsl_sf_legendre_sphPlm(l,-m,x[N/2]) * M_SQRT2;
			double theta_j = acos(x[N/2]);
			for(int i=0; i<M/2; ++i) {
				double phi_i = 2.0*M_PI*((double)i)/(double)M;
				double phi_i2 = 2.0*M_PI*((double)(M/2+i))/(double)M;
				glm += (sin(-m*phi_i)*weight)*(this->*G)(theta_j, phi_i) + (sin(-m*phi_i2)*weight)*(this->*G)(theta_j, phi_i2);
			}
		}
	} else { /* m > 0 */
		for(int j=0; j<N/2; ++j) {
			double weight = 2.0*M_PI/M * w[j] * gsl_sf_legendre_sphPlm(l,m,x[j]) * M_SQRT2;
			double theta_j = acos(x[j]);
			double weight2 = 2.0*M_PI/M * w[N-1-j] * gsl_sf_legendre_sphPlm(l,m,x[N-1-j]) * M_SQRT2;
			double theta_j2 = acos(x[N-1-j]);
			for(int i=0; i<M/2; ++i) {
				double phi_i = 2.0*M_PI*((double)i)/(double)M;
				double phi_i2 = 2.0*M_PI*((double)(M/2+i))/(double)M;
				glm += (cos(m*phi_i)*weight)*(this->*G)(theta_j, phi_i) + (cos(m*phi_i2)*weight2)*(this->*G)(theta_j2, phi_i2)
				     + (cos(m*phi_i)*weight2)*(this->*G)(theta_j2, phi_i) + (cos(m*phi_i2)*weight)*(this->*G)(theta_j, phi_i2);
			}
		}
		if(N%2==1) {
			double weight = 2.0*M_PI/M * w[N/2] * gsl_sf_legendre_sphPlm(l,m,x[N/2]) * M_SQRT2;
			double theta_j = acos(x[N/2]);
			for(int i=0; i<M/2; ++i) {
				double phi_i = 2.0*M_PI*((double)i)/(double)M;
				double phi_i2 = 2.0*M_PI*((double)(M/2+i))/(double)M;
				glm += (cos(m*phi_i)*weight)*(this->*G)(theta_j, phi_i) + (cos(m*phi_i2)*weight)*(this->*G)(theta_j, phi_i2);
			}
		}
	}

	return glm;
}

/* real spherical harmonics */
double SphericalIntegrator::calcYlmR(int l, int m, double Rvec[3]) {
	double R = vecmag(Rvec);
	double cos_theta = ((R == 0.0) ? 1.0 : (Rvec[2]/R));
	double phi = atan2(Rvec[1],Rvec[0]);

	double plm = gsl_sf_legendre_sphPlm(l,abs(m),cos_theta);
	/* real spherical harmonics */
	if(m == 0) {
		return plm;
	} else if(m < 0) {
		return M_SQRT2 * plm * sin(-m*phi);
	} else {
		return M_SQRT2 * plm * cos(m*phi);
	}
}

/* Glm is a packed array containing the even l up to lmax Spherical
 * Harmonic componenents of G^{(b)} for even b, where G^{(b)} is the bth
 * angular coefficient of the Laurent Expansion of the Lattice Green Function.
 * Glm contains the odd l up to lmax for odd b.
 * i.e. b = -2 -> Elastic Green Function
 *      b = -1 -> i/k Term
 *      b =  0 -> Discontinuity Correction
 * For each l there are 2l + 1 components
 *
 * For even b:
 * Indexed by l = 2n since only even l are non-zero (Time reversal symmetry).
 * Packing is as follows:
 * Glm[m + (l*(l+1))>>1]
 *
 * For odd b: 
 * Indexed by l = 2n+1 since only odd l are non-zero (Time reversal symmetry).
 * Packing is as follows:
 * Glm[m + (l*(l+1))>>1]
 */

Matrix SphericalIntegrator::calcG_b_R(Matrix *Glm, int b, unsigned int lmax,
                                      double kmax, double R[3], double V, unsigned int natoms)
{
	switch (abs(b)%2) {
	case 0:
		return calcG_b_R_even(Glm, b, lmax, kmax, R, V, natoms);
	case 1:
		return calcG_b_R_odd(Glm, b, lmax, kmax, R, V, natoms);
	default: /* Can't get here... */
		std::cerr << "Can't get here?" << std::endl;
		return Matrix(3u*natoms,3u*natoms);
	}
}

Matrix SphericalIntegrator::calcG_b_R_even(Matrix *Glm, int b, unsigned int lmax, double kmax, double R[3], double V, unsigned int natoms) {
	Matrix G_R(3u*natoms,3u*natoms);
	double Rmag = vecmag(R);

	/* special case: R == 0.0
	 * only need l=0, all other integrals are zero
	 */
	if(Rmag == 0.0) {
		unsigned int l = 0;
		double cweight = V/(2.0*M_PI*M_PI);
		double radialval = radialBesselIntegral(kmax, Rmag, l, b);
		double Ylm = M_2_SQRTPI/4.0; //Y00
		double term = Ylm*radialval*cweight;
		G_R = Glm[0]*term;
		return G_R;
	}

	double lweight[(lmax>>1)+1];
	/* l = 2n */
	for(unsigned int n=0; n<=(lmax>>1); ++n) {
		unsigned int l = 2u*n;
		double cweight = ((n%2) == 0 ? 1.0 : -1.0) * V/(2.0*M_PI*M_PI);
		double radialval = radialBesselIntegral(kmax, Rmag, l, b);
		lweight[n] = cweight*radialval;
#ifdef DEBUG
		std::cerr << "radial val (kmax=" << kmax << ",Rmag=" << Rmag
		          << ",l=" << l << ",b=" << b << "): " << radialval
		          << std::endl;
#endif
	}

	double sphPlm[lmax+1];
	double cos_theta = ((Rmag == 0.0) ? 1.0 : (R[2]/Rmag));
	double phi = atan2(R[1],R[0]);
	// m == 0:
	gsl_sf_legendre_sphPlm_array(lmax, 0, cos_theta, sphPlm);
	for(int n=0; n<=(lmax>>1); ++n) {
		int l = 2*n;
		double term = sphPlm[l]*lweight[n];
		G_R += Glm[((l*(l+1))>>1)]*term;
	}

	for(int m = 1; m <= (int)lmax; ++m) {
		gsl_sf_legendre_sphPlm_array(lmax, m, cos_theta, sphPlm);
		for(int n=(m + (m%2))/2; n<=(lmax>>1); ++n) {
			int l = 2*n;
			double Ylm_plus = M_SQRT2 * sphPlm[l-m] * cos(m*phi)*lweight[n];
			double Ylm_minus = M_SQRT2 * sphPlm[l-m] * sin(m*phi)*lweight[n];
#ifdef DEBUG
			std::cerr << "Ylm_plus (l=" << l << ",m=" << m << ",R={"
			          << R[0] << "," << R[1] << "," << R[2]
			          << "}): " << Ylm_plus << std::endl;
#endif
			G_R += Glm[m + ((l*(l+1))>>1)]*Ylm_plus;
			G_R += Glm[((l*(l+1))>>1) - m]*Ylm_minus;
		}
	}
	return G_R;
}

Matrix SphericalIntegrator::calcG_b_R_odd(Matrix *Glm, int b, unsigned int lmax, double kmax, double R[3], double V, unsigned int natoms) {
	Matrix G_R(3u*natoms,3u*natoms);
	double Rmag = vecmag(R);

	/* special case: R == 0.0
	 * only need l=0, all other integrals are zero
	 * Since only odd l are non-zero due to symmetry,
	 * the odd G_b_R is zero for R == 0.0
	 */
	if(Rmag == 0.0) {
		return G_R;
	}
	/* l = 2n+1 */
	double lweight[(((int)lmax-1)>>1)+1];
	for(int n=0; n<=(((int)lmax-1)>>1); ++n) {
		int l = 2*n + 1;
		double cweight = ((n%2) == 0 ? 1.0 : -1.0) * V/(2.0*M_PI*M_PI);
		double radialval = radialBesselIntegral(kmax, Rmag, l, b);
		lweight[n] = cweight*radialval;
#ifdef DEBUG
		std::cerr << "radial val (kmax=" << kmax << ",Rmag=" << Rmag
		          << ",l=" << l << ",b=" << b << "): " << radialval
		          << std::endl;
#endif
	}

	double sphPlm[lmax+1];
	double cos_theta = ((Rmag == 0.0) ? 1.0 : (R[2]/Rmag));
	double phi = atan2(R[1],R[0]);
	// m == 0:
	gsl_sf_legendre_sphPlm_array(lmax, 0, cos_theta, sphPlm);
	for(int n=0; n<=(((int)lmax-1)>>1); ++n) {
		int l = 2*n + 1;
		double term = sphPlm[l]*lweight[n];
		G_R += Glm[(n+1)*(2*n+1)]*term;
	}

	for(int m = 1; m <= (int)lmax; ++m) {
		gsl_sf_legendre_sphPlm_array(lmax, m, cos_theta, sphPlm);
		for(int n=(m - (m%2))/2; n<=(((int)lmax-1)>>1); ++n) {
			int l = 2*n + 1;
			double Ylm_plus = M_SQRT2 * sphPlm[l-m] * cos(m*phi)*lweight[n];
			double Ylm_minus = M_SQRT2 * sphPlm[l-m] * sin(m*phi)*lweight[n];
#ifdef DEBUG
			std::cerr << "Ylm_plus (l=" << l << ",m=" << m << ",R={"
			          << R[0] << "," << R[1] << "," << R[2]
			          << "}): " << Ylm_plus << std::endl;
			std::cerr << "Ylm_minus (l=" << l << ",m=" << m << ",R={"
			          << R[0] << "," << R[1] << "," << R[2]
			          << "}): " << Ylm_minus << std::endl;
			std::cerr << "Glm*Ylm_plus\n" << Glm[m + (n+1)*(2*n+1)]*Ylm_plus;
			std::cerr << "Glm*Ylm_minus\n" << Glm[-m + (n+1)*(2*n+1)]*Ylm_minus;
			std::cerr << "Glm_minus\n" << Glm[-m + (n+1)*(2*n+1)];
			std::cerr << "Glm_minus_idx\n" << (-m + (n+1)*(2*n+1)) << "\n";
			std::cerr << "Glm_minus\n" << Glm[(n+1)*(2*n+1) - m];
#endif
			G_R += Glm[m + (n+1)*(2*n+1)]*Ylm_plus;
			G_R += Glm[(n+1)*(2*n+1) - m]*Ylm_minus;
		}
	}
	return G_R;
}

double SphericalIntegrator::integrandm2(double x, void *params) {
	intparam p = *(intparam *)params;
	return p.kmax*gsl_sf_bessel_jl(p.l,x*p.kmax*p.R);
}
double SphericalIntegrator::integrandm2cut(double x, void *params) {
	intparam p = *(intparam *)params;
	return p.kmax*fcutspline(x)*gsl_sf_bessel_jl(p.l,x*p.kmax*p.R);
}

double SphericalIntegrator::integrandm1(double x, void *params) {
	intparam p = *(intparam *)params;
	return p.kmax*p.kmax*x*gsl_sf_bessel_jl(p.l,x*p.kmax*p.R);
}
double SphericalIntegrator::integrandm1cut(double x, void *params) {
	intparam p = *(intparam *)params;
	return p.kmax*p.kmax*x*fcutspline(x)*gsl_sf_bessel_jl(p.l,x*p.kmax*p.R);
}

double SphericalIntegrator::integrandm0(double x, void *params) {
	intparam p = *(intparam *)params;
	return p.kmax*p.kmax*p.kmax*x*x*gsl_sf_bessel_jl(p.l,x*p.kmax*p.R);
}
double SphericalIntegrator::integrandm0cut(double x, void *params) {
	intparam p = *(intparam *)params;
	return p.kmax*p.kmax*p.kmax*x*x*fcutspline(x)*gsl_sf_bessel_jl(p.l,x*p.kmax*p.R);
}

double SphericalIntegrator::radialBesselIntegral(double kmax, double R, int l, int b) {
	assert(b==-2 || b==-1 || b==0);
	/* special value: R==0 */
	if(R == 0.0) {
		if(l!=0) return 0.0;

		switch(b) {
			case -2:
				return kmax*(0.5*(1.0 + ALPHA_CUT));
			case -1:
				return kmax*kmax*(0.15*(1.0 + ALPHA_CUT*ALPHA_CUT) + 0.2*ALPHA_CUT);
			case 0:
				return kmax*kmax*kmax*(1.0/3.0 - (1.0 - ALPHA_CUT)*(2.0*ALPHA_CUT*ALPHA_CUT + 5.0*ALPHA_CUT + 8.0)/30.0);
			default:
				return 0.0;
		}
	}

	/* Adaptive Gauss-Kronrod routine for the general case: */
	gsl_integration_workspace *w = gsl_integration_workspace_alloc(INTWORKSZ);

	int errcode;
	double res;
	double abserr;
	gsl_function fun;
	intparam p;
	p.R = R;
	p.kmax = kmax;
	p.l = l;
	fun.params = &p;
	switch(b) {
		case -2:
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
	errcode = gsl_integration_qag(&fun, 0.0, ALPHA_CUT, TOL_A, TOL_R, INTWORKSZ, INT_TYPE, w, &res, &abserr);

#ifdef DEBUG
	std::cerr << "Integration Error Estimate 1: " << abserr << std::endl;
#endif
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

SphericalIntegrator::SphericalIntegrator(DynMat &dynmat, Matrix (*one_on_k2)(DynMat &, double[3]), Matrix (*i_on_k)(DynMat &, double[3]), Matrix (*k0)(DynMat &, double[3])) {
	this->gexpan[0] = one_on_k2;
	this->gexpan[1] = i_on_k;
	this->gexpan[2] = k0;
	this->dynmat = &dynmat;
}

Matrix SphericalIntegrator::one_on_k2_angle(double theta, double phi) {
	double kvec[3];
	kvec[0] = sin(theta)*cos(phi);
	kvec[1] = sin(theta)*sin(phi);
	kvec[2] = cos(theta);
	return (*gexpan[0])(*dynmat, kvec);
}

Matrix SphericalIntegrator::i_on_k_angle(double theta, double phi) {
	double kvec[3];
	kvec[0] = sin(theta)*cos(phi);
	kvec[1] = sin(theta)*sin(phi);
	kvec[2] = cos(theta);
	return (*gexpan[1])(*dynmat, kvec);
}

Matrix SphericalIntegrator::k0_angle(double theta, double phi) {
	double kvec[3];
	kvec[0] = sin(theta)*cos(phi);
	kvec[1] = sin(theta)*sin(phi);
	kvec[2] = cos(theta);
	return (*gexpan[2])(*dynmat, kvec);
}
