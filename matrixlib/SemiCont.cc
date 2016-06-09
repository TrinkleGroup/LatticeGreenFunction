/* File: SemiCont.cc
 * Semi-continuum correction for the lattice Green function
 * Numerical calculated for the non-divergent pieces
 */
#include <cmath>
#include <fstream>
#include <iostream>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "SemiCont.hh"
#include "SphericalIntegrator.hh"
#include "CutOffFcn.hh"
#include "UnitCell.hh"
#include "VecMath.hh"
#include "SystemDimension.hh"

#define MAX(x,y) ((x>y)?x:y)

const double SemiCont::SMALLK_TOL;
const double SemiCont::LOG_SMALLK_TOL;
const int SemiCont::SMALLK_NMAX;

SemiCont::SemiCont() { }


SemiCont::SemiCont(DynMat &dmat, InvExpPtrs fc, double *kpoints, int nk, UnitCell &uc)
	: dynmat(&dmat), ucell(uc), nkpt(nk), fcptrs(fc)
{
	this->kpoints = new double[CARTDIM*nk];
	for(int i=0; i<nk*CARTDIM; ++i) {
		this->kpoints[i] = kpoints[i];
	}
	kmax = uc.getKmax();
	ksmall = 0.1*kmax;
	//std::cerr << "kmax\n" << kmax << std::endl;
}

SemiCont::SemiCont(DynMat &dmat, InvExpPtrs fc, std::ifstream &kpt_file, UnitCell &uc)
	: dynmat(&dmat), ucell(uc), fcptrs(fc)
{
	load_kpoints(kpt_file);
	kmax = ucell.getKmax();
	ksmall = 0.1*kmax;
	//std::cerr << "kmax\n" << kmax << std::endl;
}

SemiCont::~SemiCont() {
	delete [] kpoints;
}

/* K points format:
 * num_kpoints(N)
 * k1_1 k1_2 .. k1_CARTDIM
 * ...
 * kN_1 kN_2 .. kN_CARTDIM
 */
void SemiCont::load_kpoints(std::ifstream &ifs) {
	ifs >> nkpt;
	kpoints = new double[CARTDIM*nkpt];
	for(unsigned int i=0u; i<nkpt; ++i) {
		for(unsigned int j=0u; j<CARTDIM; ++i) {
			ifs >> kpoints[CARTDIM*i + j];
		}
	}
}

/* calculate the complex semicontinuum matrix
 * handle both small k values and direct subtraction for large k
 */
void SemiCont::calcSemiCont(ZMatrix *gsc) {
	//gsc = new ZMatrix[nkpt];

	for(int i=0; i < nkpt; ++i) {
		double *curr_k = kpoints + CARTDIM*i;
		double kmag = vecmag(curr_k);

		/* Gamma point -> SemiCont is zero by definition */
		if(kmag == 0.0) {
			gsc[i] = ZMatrix(CARTDIM*dynmat->getNions(), CARTDIM*dynmat->getNions());
		} else if(kmag < ksmall) {
			gsc[i] = smallKexp(curr_k);
		} else { 
			gsc[i] = directSub(curr_k);
		}
	}

}

/* calculate semicontinuum piece at kpoint kpt
 * via direct subtraction of the analytic divergent acoustic pieces
 */
ZMatrix SemiCont::directSub(double *kpt) {
	double cutoff = 1.0;
	double kmag = vecmag(kpt);
	double kdir[CARTDIM];
	for(unsigned int i = 0u; i<CARTDIM ; ++i) {
		kdir[i] = kpt[i]/kmag;
	}
	ZMatrix ft(dynmat->getNions()*CARTDIM, dynmat->getNions()*CARTDIM);
	Matrix zero(dynmat->getNions()*CARTDIM, dynmat->getNions()*CARTDIM);
	dynmat->getFourierTransformRot(kpt, ft);
	ZMatrix ret = ft.inverse();
	if(kmag >= kmax) { //cutoff = 0.0, outside sphere
		return ret;
	} else if(kmag > kmax*ALPHA_CUT) {
		cutoff = fcutspline(kmag/kmax);
	}
	ret -= ZMatrix((*(fcptrs.one_on_k2))(*dynmat, kdir) * (cutoff/kmag/kmag));
	if(dynmat->getNions() != 1) {
		ret -= ZMatrix(zero, (*(fcptrs.i_on_k))(*dynmat, kdir)*(cutoff/kmag));
	}
	ret -= ZMatrix((*(fcptrs.k0))(*dynmat, kdir)*cutoff);

	return ret;
}

/* calculate semicontinuum piece at kpoint kpt
 * via a small k expansion the analytic divergent acoustic pieces
 * XiInvK is the leading order pieces in each quadrant of the Green function
 * matrix
 */
/* currently not used
ZMatrix SemiCont::smallKtaylor(double *kpt) {
	double kmag = vecmag(kpt);
	double kdir[CARTDIM];
	for(unsigned int i = 0u; i < CARTDIM; ++i) {
		kdir[i] = kpt[i]/kmag;
	}

	//std::cerr << "Small K Exp called: kmag = " << kmag << std::endl;

	int msz = dynmat->getNions()*CARTDIM;
	ZMatrix smft(msz,msz);
	ZMatrix XiInvK(msz,msz);
	dynmat->getSmallFourierTransform(kpt, smft);
	XiInvK.setReal((*(fcptrs.XiInvAA))(*dynmat,kdir)*(1.0/kmag/kmag) + (*(fcptrs.XiInvOO))(*dynmat,kdir));
	XiInvK.setImag((*(fcptrs.XiInvAO))(*dynmat,kdir)*(1.0/kmag));

	gsl_complex i_on_kmag;
	GSL_SET_COMPLEX(&i_on_kmag, 0.0, 1.0/kmag);
	ZMatrix B = XiInvK * smft;

	double condest = B.cond1est();
	double condinv = 1.0/condest;

	// Little bit of an issue here. the max eigenvalue doesn't
	// really determine stability as in the single atom case
	// because of those pesky optical modes.  We really want
	// just the eigenvalues of the acoustic piece for comparison

	//std::cerr << "condest = " << condest << std::endl;
	//std::cerr << "1/condest = " << condinv << std::endl;


	int nseries; 
	if(condinv > 0.0 && condinv < 1.0) {
		nseries = lround(-LOG_SMALLK_TOL/std::log(condest));
	} else if(condinv > 0.0) {
		nseries = SMALLK_NMAX + 1;
	} else {
		nseries = 1;
	}

	ZMatrix eye = ZMatrix::eye(msz);

	ZMatrix pk(msz,msz);

	if(nseries > SMALLK_NMAX) {
		//std::cerr << "Direct Inverse. nseries = " << nseries << std::endl;
		// Direct inverse
		ZMatrix A = eye - B;
		A = A.inverse() - eye;
		pk = A*XiInvK;
		pk -= (*(fcptrs.i_on_k_sub))(*dynmat,kdir)*i_on_kmag;
		pk -= (*(fcptrs.k0_sub))(*dynmat,kdir)*GSL_COMPLEX_ONE;
	} else {
		//std::cerr << "Power Inverse. nseries = " << nseries << std::endl;
		// Power Series Inversion
		ZMatrix tmpA = eye+B;
		for(int n=0; n<nseries+CARTDIM; ++n) {
			tmpA = eye + B*tmpA;
		}
		pk = ((B*XiInvK - (*(fcptrs.i_on_k_sub))(*dynmat,kdir)*i_on_kmag) - (*(fcptrs.k0_sub))(*dynmat,kdir)*GSL_COMPLEX_ONE);
		pk += B*B*tmpA*XiInvK;
	}

	return pk;
}
*/

/* calculate semicontinuum piece at kpoint kpt
 * via a small k expansion the analytic divergent acoustic pieces
 * for very small K
 * using the small K fourier transform.
 * XiInvK is the leading order pieces in each quadrant of the Green function
 * matrix
 */
ZMatrix SemiCont::smallKexp(double *kpt) {
	double kmag = vecmag(kpt);
	double kdir[CARTDIM];
	for(unsigned int i = 0u; i < CARTDIM; ++i) {
		kdir[i] = kpt[i]/kmag;
	}


	int msz = dynmat->getNions()*CARTDIM;
	ZMatrix smft(msz,msz);
	ZMatrix XiInvK(msz,msz);
	dynmat->getSmallFourierTransform(kpt, smft);
	//std::cerr << "Small K Exp called: kmag = " << kmag << "\n" << "smft = \n" << smft << std::endl;
	gsl_complex i_on_kmag;

	if(dynmat->getNions() == 1) {
		XiInvK.setReal((*(fcptrs.XiInvAA))(*dynmat,kdir)*(1.0/kmag/kmag));
	} else {
		XiInvK.setReal((*(fcptrs.XiInvAA))(*dynmat,kdir)*(1.0/kmag/kmag) + (*(fcptrs.XiInvOO))(*dynmat,kdir));
		XiInvK.setImag((*(fcptrs.XiInvAO))(*dynmat,kdir)*(1.0/kmag));

		GSL_SET_COMPLEX(&i_on_kmag, 0.0, 1.0/kmag);
	}

	ZMatrix B = XiInvK * smft;

	//std::cerr << "XiInvK =\n" << XiInvK << std::endl;

	/*
	double condest = B.cond1est();
	double condinv = 1.0/condest;
	*/

	/* Little bit of an issue here. the max eigenvalue doesn't
	 * really determine stability as in the single atom case
	 * because of those pesky optical modes.  We really want
	 * just the eigenvalues of the acoustic piece for comparison
	 */

	/*
	std::cerr << "condest = " << condest << std::endl;
	std::cerr << "1/condest = " << condinv << std::endl;
	*/


	/*
	int nseries; 
	if(condinv > 0.0 && condinv < 1.0) {
		nseries = lround(-LOG_SMALLK_TOL/std::log(condest));
	} else if(condinv > 0.0) {
		nseries = SMALLK_NMAX + 1;
	} else {
		nseries = 1;
	}
	*/

	ZMatrix eye = ZMatrix::eye(msz);

	ZMatrix pk(msz,msz);

	//std::cerr << "Direct Inverse. nseries = " << nseries << std::endl;
	// Direct inverse
	ZMatrix A = eye - B;
	A = A.inverse() - eye;
	//std::cerr << "A =\n" << A << std::endl;
	//std::cerr << "B =\n" << B << std::endl;
	pk = A*XiInvK;
	//std::cerr << "pk1 =\n" << pk << std::endl;
	if(dynmat->getNions() > 1) {
		pk -= (*(fcptrs.i_on_k_sub))(*dynmat,kdir)*i_on_kmag;
	}
	pk -= (*(fcptrs.k0_sub))(*dynmat,kdir)*GSL_COMPLEX_ONE;

	//std::cerr << "pk0 =\n" << pk << std::endl;
	return pk;
}
