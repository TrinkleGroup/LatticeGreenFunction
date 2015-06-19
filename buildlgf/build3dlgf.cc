#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "gsl_complex_op.hh"
#include "UnitCell.hh"
#include "SemiCont.hh"
#include "Matrix.hh"
#include "ZMatrix.hh"
#include "DynMat.hh"
#include "SphericalIntegrator.hh"
#include "VecMath.hh"
#include "allmats.hh"

#define L_MAX_E 31u
#define L_MAX_O 31u

using namespace std;

void load_Rpoints(istream &ifs, unsigned int &n, double * &Rpoints, UnitCell &uc) {
	double scale;
	ifs >> n;
	ifs >> scale;
	string s;
	getline(ifs, s);
	int idx = 0;
	while(s.length() > idx && (s[idx] == ' ' || s[idx] == '\t'))
		++idx;
	char c = 'L';
	if(s.length() > idx) {
		c = s[idx];
	}
	Rpoints = new double[3u*n];
	if(c == 'L') { //Reciprocal lattice coords
		Matrix avec = uc.getAvec();
		double lsc = uc.getScale();
		for(unsigned int i=0u; i<n; ++i) {
			double rpttmp[3];
			ifs >> rpttmp[0]
			    >> rpttmp[1]
			    >> rpttmp[2];
			Rpoints[3u*i] = (rpttmp[0]*avec.val(0,0)
			              +  rpttmp[1]*avec.val(1,0)
			              +  rpttmp[2]*avec.val(2,0))*lsc*scale;
			Rpoints[3u*i+1u] = (rpttmp[0]*avec.val(0,1)
			                 +  rpttmp[1]*avec.val(1,1)
			                 +  rpttmp[2]*avec.val(2,1))*lsc*scale;
			Rpoints[3u*i+2u] = (rpttmp[0]*avec.val(0,2)
			                 +  rpttmp[1]*avec.val(1,2)
			                 +  rpttmp[2]*avec.val(2,2))*lsc*scale;
		}
	} else {
		for(unsigned int i=0u; i<n; ++i) {
			ifs >> Rpoints[3u*i]
			    >> Rpoints[3u*i + 1u]
			    >> Rpoints[3u*i + 2u];
			Rpoints[3u*i] *= scale;
			Rpoints[3u*i+1u] *= scale;
			Rpoints[3u*i+2u] *= scale;
		}
	}
}

/* K points format: (Cartesian coordinates)
 * num_kpoints(N) scale, coord
 * k1_1 k1_2 k1_3 weight
 * ...
 * kN_1 kN_2 kN_3 weight
 */
void load_kpoints(istream &ifs, unsigned int &n, double * &kpoints, double * &weights, UnitCell &uc) {
	double scale;
	ifs >> n;
	ifs >> scale;
	string s;
	getline(ifs, s);
	int idx = 0;
	while(s.length() > idx && (s[idx] == ' ' || s[idx] == '\t'))
		++idx;
	char c = 'L';
	if(s.length() > idx) {
		c = s[idx];
	}
	kpoints = new double[3u*n];
	weights = new double[n];
	if(c == 'L') { //Reciprocal lattice coords
		Matrix bvec = uc.getBvec();
		double rscale = uc.getRscale();
		//cerr << "bvec =\n" << bvec;
		//cerr << "rscale =\n" << rscale << "\n";
		//cerr << "scale =\n" << scale << endl;
		for(unsigned int i=0u; i<n; ++i) {
			double kpttmp[3];
			ifs >> kpttmp[0]
			    >> kpttmp[1]
			    >> kpttmp[2]
			    >> weights[i]; //weight
			kpoints[3u*i] = (kpttmp[0]*bvec.val(0,0)
			              +  kpttmp[1]*bvec.val(1,0)
			              +  kpttmp[2]*bvec.val(2,0))*rscale*scale;
			kpoints[3u*i+1u] = (kpttmp[0]*bvec.val(0,1)
			                 +  kpttmp[1]*bvec.val(1,1)
			                 +  kpttmp[2]*bvec.val(2,1))*rscale*scale;
			kpoints[3u*i+2u] = (kpttmp[0]*bvec.val(0,2)
			                 +  kpttmp[1]*bvec.val(1,2)
			                 +  kpttmp[2]*bvec.val(2,2))*rscale*scale;
		}
	} else {
		for(unsigned int i=0u; i<n; ++i) {
			ifs >> kpoints[3u*i]
			    >> kpoints[3u*i + 1u]
			    >> kpoints[3u*i + 2u]
			    >> weights[i]; //weight
			kpoints[3u*i] *= scale;
			kpoints[3u*i+1u] *= scale;
			kpoints[3u*i+2u] *= scale;
		}
	}
}

static Matrix glm_m2[L_MAX_E + ((L_MAX_E*(L_MAX_E+1u))>>1) + 1u];
static Matrix glm_m1[L_MAX_O + ((L_MAX_O*(L_MAX_O+1u))>>1) + 1u];
static Matrix glm_0[L_MAX_E + ((L_MAX_E*(L_MAX_E+1u))>>1) + 1u];

void calcGlms(DynMat &dynmat) {
	SphericalIntegrator si(dynmat, &one_on_k_2_total_rot, &i_on_k_total_rot, &k0_total_rot);

	double ktest[3] = {0.0, 0.0, 1.0};

	/* check one_on_k2...*/
	//cerr << "ktest-2\n" << one_on_k_2_total_rot(dynmat, ktest) << endl;
	//cerr << "ktest0\n" << k0_total_rot(dynmat, ktest) << endl;

	si.calcGlmExp(-2, L_MAX_E, glm_m2);
	si.calcGlmExp(-1, L_MAX_O, glm_m1);
	si.calcGlmExp(0, L_MAX_E, glm_0);

	/*
	cerr << "G^{-1}(k) l=1,m=-1:\n" << glm_m1[0];
	cerr << "G^{-1}(k) l=1,m=0:\n" << glm_m1[1];
	cerr << "G^{-1}(k) l=1,m=1:\n" << glm_m1[2];
	cerr << "G^{-1}(k) l=3,m=-3:\n" << glm_m2[3];
	cerr << "G^{-1}(k) l=3,m=-2:\n" << glm_m2[4];
	cerr << "G^{-1}(k) l=3,m=-1:\n" << glm_m2[5];
	cerr << "G^{-1}(k) l=3,m=0:\n" << glm_m2[6];
	cerr << "G^{-1}(k) l=3,m=1:\n" << glm_m2[7];
	cerr << "G^{-1}(k) l=3,m=2:\n" << glm_m2[8];
	cerr << "G^{-1}(k) l=3,m=3:\n" << glm_m2[9];
	cerr << "G^{-1}(k) l=5,m=-5:\n" << glm_m2[10];
	cerr << "G^{-1}(k) l=5,m=-4:\n" << glm_m2[11];
	cerr << "G^{-1}(k) l=5,m=-3:\n" << glm_m2[12];
	cerr << "G^{-1}(k) l=5,m=-2:\n" << glm_m2[13];
	cerr << "G^{-1}(k) l=5,m=-1:\n" << glm_m2[14];
	cerr << "G^{-1}(k) l=5,m=0:\n" << glm_m2[15];
	*/

}

ZMatrix *Gsc;// = new ZMatrix[nkpt];

void doSemiCont(DynMat &dynmat, UnitCell &uc, double *kpoints, unsigned int nkpt) {
	InvExpPtrs fcptrs;

	fcptrs.one_on_k2 = &one_on_k_2_total_rot;
	fcptrs.i_on_k = &i_on_k_total_rot;
	fcptrs.k0 = &k0_total_rot;

	fcptrs.XiInvAA = &XiInvAA;
	fcptrs.XiInvAO = &XiInvAO;
	fcptrs.XiInvOO = &XiInvOO;

	fcptrs.i_on_k_sub = &i_on_k_sub;
	fcptrs.k0_sub = &k0_sub;

	SemiCont sc(dynmat, fcptrs, kpoints, nkpt, uc);
	sc.calcSemiCont(Gsc);
	cerr << "SemiCont done..." << endl;

	/*
	for(int i = 0; i < nkpt; ++i) {
		cerr << "k = (" << kpoints[3u*i+0u] << ","
		     << kpoints[3u*i+1u] << ","
		     << kpoints[3u*i+2u] << ")\n"
		     << Gsc[i].getReal()
		     << Gsc[i].getImag();
	}
	cerr.flush();
	*/
}

Matrix getGbr(DynMat &dynmat, double R[3], UnitCell &uc, double *kpoints, double *weights, unsigned int nkpt) {
	Matrix gbrm2 = SphericalIntegrator::calcG_b_R(glm_m2, -2, L_MAX_E, uc.getKmax(), R, uc.getVolume(), dynmat.getNions());
	Matrix gbrm1 = SphericalIntegrator::calcG_b_R(glm_m1, -1, L_MAX_O, uc.getKmax(), R, uc.getVolume(), dynmat.getNions());
	Matrix gbr0 = SphericalIntegrator::calcG_b_R(glm_0, 0, L_MAX_E, uc.getKmax(), R, uc.getVolume(), dynmat.getNions());
	cerr << "R=(" << R[0] << "," << R[1] << "," << R[2] << ")\n";
	//cerr << "rbi=" << SphericalIntegrator::radialBesselIntegral(uc.getKmax(), vecmag(R), 1, -1) << "\n";
	cerr << "gbrm2=\n" << gbrm2;
	cerr << "gbrm1=\n" << gbrm1;
	cerr << "gbr0=\n" << gbr0;


	/* calculate the contribution from the in plane semi-continuum piece
	 * Assume time reversal symmetry k -> -k. G_{(sc)} -> G^{*}_{(sc)}.
	 */
	Matrix gR_sc(3u*dynmat.getNions(), 3u*dynmat.getNions());
	for(unsigned int i=0u; i<nkpt; ++i) {
		double *curr_k = kpoints + 3u*i;
		double kdotR = vecdot(curr_k,R);
		gR_sc += (cos(kdotR)*(Gsc[i].getReal()) + sin(kdotR)*(Gsc[i].getImag()))*weights[i];
	}

	cerr << "Semicontinuum:\n" << gR_sc;


	return gbrm2 + gbrm1 + gbr0 + gR_sc;
}

int main(int argc, char *argv[]) {
	if(argc != 5) {
		cerr << "Usage: " << argv[0] << " dynamical_matrix kpoints cellvec rvecfile"
		          << endl;
		return(-1);
	}
	double *kpoints, *weights;
	unsigned int num_kpoints;
	DynMat dynmat;
	ifstream dmstream(argv[1]);
	ifstream kpstream(argv[2]);
	ifstream cellstream(argv[3]);
	ifstream rvecstream(argv[4]);

	Matrix avec(3u,3u);

	cout.precision(16);
	cout << scientific;
	cout.setf(cout.showpos);
	cerr.precision(16);
	cerr << scientific;
	cerr.setf(cerr.showpos);

	dynmat.load(dmstream);

	UnitCell uc(cellstream);

	load_kpoints(kpstream, num_kpoints, kpoints, weights, uc);

	double *Rpoints;
	unsigned int num_rpoints;
	load_Rpoints(rvecstream, num_rpoints, Rpoints, uc);


	rvecstream.close();
	cellstream.close();
	dmstream.close();
	kpstream.close();

	calcGlms(dynmat);
	Gsc = new ZMatrix[num_kpoints];

	doSemiCont(dynmat, uc, kpoints, num_kpoints);

	for(double *R = Rpoints; R < Rpoints + 3u*num_rpoints; R += 3) {
		Matrix gbr = getGbr(dynmat, R, uc, kpoints, weights, num_kpoints);

		cout << R[0] << " " << R[1] << " " << R[2] << "\n" << gbr;
	}
	cout.flush();

	delete [] Gsc;
	delete [] Rpoints;
	delete [] kpoints;
	delete [] weights;

	return 0;
}
