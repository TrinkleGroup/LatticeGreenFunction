/* Calculate the 2D lattice Green function for dislocations */
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
#include "PolarIntegrator.hh"
#include "VecMath.hh"
#include "allmats.hh"
#include "SystemDimension.hh"

#define L_MAX_N 16u
#define TOL_PLANE 1e-12

using namespace std;

/* load in the R vectors to calculate the Green function */
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
	Rpoints = new double[CARTDIM*n];
	if(c == 'L') { //Reciprocal lattice coords
		Matrix avec = uc.getAvec();
		double lsc = uc.getScale();
		for(unsigned int i=0u; i<n; ++i) {
			double rpttmp[CARTDIM];
			for(int jdx = 0; jdx < CARTDIM; ++jdx) {
				ifs >> rpttmp[jdx];
			}
			for(int kdx = 0; kdx < CARTDIM; ++kdx) {
				for(int jdx = 0; jdx < CARTDIM; ++jdx) {
					Rpoints[CARTDIM*i + kdx] += rpttmp[jdx]*avec.val(jdx,kdx);
				}
				Rpoints[CARTDIM*i + kdx] *= lsc*scale;
			}
		}
	} else {
		for(unsigned int i=0u; i<n; ++i) {
			for(int jdx = 0; jdx < CARTDIM; ++jdx) {
				ifs >> Rpoints[CARTDIM*i + jdx];
				Rpoints[CARTDIM*i + jdx] *= scale;
			}
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
	kpoints = new double[CARTDIM*n];
	weights = new double[n];
	if(c == 'L') { //Reciprocal lattice coords
		Matrix bvec = uc.getBvec();
		double rscale = uc.getRscale();
		for(unsigned int i=0u; i<n; ++i) {
			double kpttmp[CARTDIM];
			for(int idx = 0; idx < CARTDIM; ++idx) {
				ifs >> kpttmp[idx];
			}
			ifs >> weights[i]; //weight
			for(int jdx = 0; jdx < CARTDIM; ++jdx) {
				kpoints[CARTDIM*i + jdx] = 0.0;
				for(int idx = 0; idx < CARTDIM; ++idx) {
					kpoints[CARTDIM*i + jdx] += kpttmp[idx]*bvec.val(idx,jdx);
				}
				kpoints[CARTDIM*i + jdx] *= rscale*scale;
			}
		}
	} else {
		for(unsigned int i=0u; i<n; ++i) {
			for(int idx = 0; idx < CARTDIM; ++idx) {
				ifs >> kpoints[CARTDIM*i + idx];
				kpoints[CARTDIM*i + idx] *= scale;
			}
			ifs >> weights[i]; //weight
		}
	}
}

/* Polar Fourier expansion matrices for the G^-2, G^-1 and G^0 terms */
static Matrix gn_m2[L_MAX_N/2];
static Matrix gn_m1[L_MAX_N/2];
static Matrix gn_0[L_MAX_N/2];

/* Calculate the polar expansions of the divergent Green function pieces
 * dynmat: Dynamical Matrix
 * t: dislocation line direction
 * m: dislocation cut direction
 */
void calcGns(DynMat &dynmat, double t[3], double *m) {
/*	cout << dynmat.getNions() << " " << n << "\n";
	Matrix rot;
	Matrix rotT;
	getAcousticRots(dynmat, rot, rotT);
*/
	PolarIntegrator *pi = NULL;
	if(dynmat.getNions() > 1) {
		pi = new PolarIntegrator(dynmat, t, m, &one_on_k_2_total_rot, &i_on_k_total_rot, &k0_total_rot);
	} else {
		pi = new PolarIntegrator(dynmat, t, m, &one_on_k_2_1a, NULL, &k0_1a);
	}

	pi->calcGnExp(-2, L_MAX_N, gn_m2);
#ifdef DEBUG
	cerr << "k^-2 expansion nmax=" << L_MAX_N << endl;
	for(unsigned int n=0u; n<L_MAX_N/2; ++n) {
		cerr << "n=" << 2*n << "\n";
		cerr << gn_m2[n];
	}
	cerr.flush();
#endif
	if(dynmat.getNions() > 1) {
		pi->calcGnExp(-1, L_MAX_N, gn_m1);
#ifdef DEBUG
		cerr << "k^-1 expansion nmax=" << L_MAX_N << endl;
		for(unsigned int n=0u; n<L_MAX_N/2; ++n) {
			cerr << "n=" << 2*n+1 << "\n";
			cerr << gn_m1[n];
		}
		cerr.flush();
#endif
	}
	pi->calcGnExp(0, L_MAX_N, gn_0);
	delete pi;
#ifdef DEBUG
	cerr << "k^0 expansion nmax=" << L_MAX_N << endl;
	for(unsigned int n=0u; n<L_MAX_N/2; ++n) {
		cerr << "n=" << 2*n << "\n";
		cerr << gn_0[n];
	}
	cerr.flush();
#endif
}

/* Semi-continuum matrix */
static ZMatrix *Gsc;

/* calculate the semi-continuum correction
 * dynmat: Dynamical matrix
 * uc: crystal unit cell
 * inkpoints: kpoints in the dislocation slab plane perp to t
 * ninkpt: number of kpoints
 */
void doSemiCont(DynMat &dynmat, UnitCell &uc, double *inkpoints, unsigned int ninkpt) {
	InvExpPtrs fcptrs;

	if(dynmat.getNions() > 1) {
		fcptrs.one_on_k2 = &one_on_k_2_total_rot;
		fcptrs.i_on_k = &i_on_k_total_rot;
		fcptrs.k0 = &k0_total_rot;
		fcptrs.XiInvAA = &XiInvAA;
		fcptrs.XiInvAO = &XiInvAO;
		fcptrs.XiInvOO = &XiInvOO;
		fcptrs.i_on_k_sub = &i_on_k_sub;
		fcptrs.k0_sub = &k0_sub;
	} else {
		fcptrs.one_on_k2 = &one_on_k_2_1a;
		fcptrs.k0 = &k0_1a;
		fcptrs.XiInvAA = &XiInv_1a;
		fcptrs.k0_sub = &k0_sub_1a;
	}

	SemiCont sc(dynmat, fcptrs, inkpoints, ninkpt, uc);
	sc.calcSemiCont(Gsc);
	/*
	cerr << "SemiCont done..." << endl;
	for(int i = 0; i < ninkpt; ++i) {
		cerr << Gsc[i];
	}
	cerr.flush();
	*/
}

/* Calculate the real space lattice Green function pieces
 * dynmat: dynamical matrix
 * t: dislocation line direction
 * m: dislocation cut direction
 * R: real space vector of KGF
 * uc: crystal unit cell
 * inkpoints: kpoints perp to t
 * ninkpt: number of kpoints perp to t
 * outkpoints: kpoints with component along t
 * noutkpt: number of kpoints with component along t
 */
Matrix getGbr(DynMat &dynmat, double t[3], double *m, double *R, UnitCell &uc, double *inkpoints, unsigned int ninkpt, double *outkpoints, unsigned int noutkpt) {
	PolarIntegrator *pi = NULL;
	if(dynmat.getNions() > 1) {
		pi = new PolarIntegrator(dynmat, t, m, &one_on_k_2_total_rot, &i_on_k_total_rot, &k0_total_rot);
	} else {
		pi = new PolarIntegrator(dynmat, t, m, &one_on_k_2_1a, NULL, &k0_1a);
	}

	Matrix gbrm2(pi->calcG_b_R(gn_m2, -2, L_MAX_N, uc.getKmax(), R, uc.getVolume()));
	Matrix gbrm1(CARTDIM*dynmat.getNions(),CARTDIM*dynmat.getNions());
	if(dynmat.getNions() > 1) {
		gbrm1 = pi->calcG_b_R(gn_m1, -1, L_MAX_N, uc.getKmax(), R, uc.getVolume());
	}
	Matrix gbr0(pi->calcG_b_R(gn_0, 0, L_MAX_N, uc.getKmax(), R, uc.getVolume()));
	delete pi;
#ifdef DEBUG
	cerr << "R = (" << R[0] << "," << R[1] << ")\n";
	cerr << "g(-2)(R) =\n" << gbrm2;
	if(dynmat.getNions() > 1) {
		cerr << "g(-1)(R) =\n" << gbrm1;
	}
	cerr << "g( 0)(R) =\n" << gbr0;
#endif



	/* calculate the contribution from the in plane semi-continuum piece
	 * Assume time reversal symmetry k -> -k. G_{(sc)} -> G^{*}_{(sc)}.
	 */
	Matrix gR_sc(CARTDIM*dynmat.getNions(), CARTDIM*dynmat.getNions());
	for(unsigned int i=0u; i<ninkpt; ++i) {
		double *curr_k = inkpoints + CARTDIM*i;
		double kdotR = vecdot(curr_k,R);
		gR_sc += cos(kdotR)*(Gsc[i].getReal()) + sin(kdotR)*(Gsc[i].getImag());
	}
	gR_sc *= 1.0/(ninkpt+noutkpt);

#ifdef DEBUG
	cerr << "Semicontinuum:\n" << gR_sc;

	cerr.flush();
#endif

	return gbrm2 + gbrm1 + gbr0 + gR_sc;
}

/* load in the dislocation description file from dislstream */
void load_disl(istream &dislstream, const UnitCell &uc, double t[3], double *m) {
	string s;
	do {
		getline(dislstream,s);
	} while(s[0] == '#' && !dislstream.eof());

	istringstream iss(s);
	int t_unit[3];
	iss >> t_unit[0] >> t_unit[1] >> t_unit[2];
	iss.clear();
	//Burger's vector
	do {
		getline(dislstream,s);
	} while(s[0] == '#' && !dislstream.eof());
	//cut vector
	do {
		getline(dislstream,s);
	} while(s[0] == '#' && !dislstream.eof());
	iss.str(s);
	int m_unit[CARTDIM];
	for(int idx = 0; idx < CARTDIM; ++idx) {
		iss >> m_unit[idx];
	}
	Matrix avec = uc.getAvec();
	double scale = uc.getScale();
	for(int idx = 0; idx < CARTDIM; ++idx) {
		m[idx] = 0.0;
		for(int jdx = 0; jdx < CARTDIM; ++jdx) {
			m[idx] += m_unit[jdx]*avec.val(jdx,idx);
		}
		m[idx] *= scale;
	}

	if(CARTDIM == 2) {
		t[0] = 0.0;
		t[1] = 0.0;
		t[2] = 1.0;
	}
}

int main(int argc, char *argv[]) {
	if(argc != 6) {
		cerr << "Usage: " << argv[0] << " dynamical_matrix kpoints cellvec dislocfile rvecfile"
		     << endl;
		return(-1);
	}
	double *kpoints, *weights;
	unsigned int num_kpoints;
	SystemDimension::initialize(2u);
	DynMat dynmat;
	ifstream dmstream(argv[1]);
	ifstream kpstream(argv[2]);
	ifstream cellstream(argv[3]);
	ifstream dislstream(argv[4]);
	ifstream rvecstream(argv[5]);

	Matrix avec(CARTDIM,CARTDIM);

	cout.precision(16);
	cout << scientific;
	cout.setf(cout.showpos);
	cerr.precision(16);
	cerr << scientific;
	cerr.setf(cerr.showpos);

	/* load our input data */

	UnitCell uc(cellstream);

	dynmat.load(dmstream, uc);


	load_kpoints(kpstream, num_kpoints, kpoints, weights, uc);

	double *Rpoints;
	unsigned int num_rpoints;
	load_Rpoints(rvecstream, num_rpoints, Rpoints, uc);


	double t[3];
	double m[CARTDIM];


	load_disl(dislstream, uc, t, m);

	rvecstream.close();
	dislstream.close();
	cellstream.close();
	dmstream.close();
	kpstream.close();

	//unsigned int n_inplane_kpt, n_outplane_kpt;
	//double *inplane_kpt, *outplane_kpt, *inplane_weights, *outplane_weights;

	// split kpoints into in threading plane, and out of threading plane
	//inplane_kpoints(kpoints, weights, num_kpoints, inplane_kpt, inplane_weights, n_inplane_kpt, t);
	//outplane_kpoints(kpoints, weights, num_kpoints, outplane_kpt, outplane_weights, n_outplane_kpt, t);

	//calculate polar expansion of the lattice Green function divergent
	//terms
	calcGns(dynmat, t, m);
	Gsc = new ZMatrix[num_kpoints];

	//calculate semi-continuum piece
	doSemiCont(dynmat, uc, kpoints, num_kpoints);

	Matrix rot(dynmat.getRot());
	Matrix rotT(dynmat.getRotT());

	//rotate back from acoustic/optical to cartesian basis
	for(double *R = Rpoints; R < Rpoints + CARTDIM*num_rpoints; R += CARTDIM) {
		Matrix gbr = rotT*getGbr(dynmat, t, m, R, uc, kpoints, num_kpoints, NULL, 0)*rot;

		cout << R[0];
		for(int idx = 1; idx < CARTDIM; ++idx) {
			cout << " " << R[idx];
		}
		cout << "\n" << gbr;
	}
	cout.flush();

	delete [] Gsc;
	delete [] Rpoints;
	delete [] kpoints;
	delete [] weights;
	/*
	delete [] inplane_kpt;
	delete [] outplane_kpt;
	delete [] inplane_weights;
	delete [] outplane_weights;
	*/

	return 0;
}
