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

#define L_MAX_N 32u
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

/* calculate k points the are in the plane of the dislocation */
void inplane_kpoints(double *kpoints, double *weights, unsigned int nkpt, double * &plane_kpt, double * &plane_weights, unsigned int &n_plane_kpt, double *t) {
	n_plane_kpt = 0;
	for(unsigned int i=0u; i<nkpt; ++i) {
		if(abs(vecdot(&kpoints[3u*i],t)) < TOL_PLANE) {
			n_plane_kpt++;
		}
	}
	plane_kpt = new double[3u*n_plane_kpt];
	plane_weights = new double[n_plane_kpt];
	unsigned int j = 0u;
	for(unsigned int i=0u; i<nkpt; ++i) {
		if(abs(vecdot(&kpoints[3u*i],t)) < TOL_PLANE) {
			plane_kpt[3u*j] = kpoints[3u*i];
			plane_kpt[3u*j+1u] = kpoints[3u*i+1u];
			plane_kpt[3u*j+2u] = kpoints[3u*i+2u];
			plane_weights[j] = weights[i];
			j++;
		}
	}
}

/* calculate k points along the line direction */
void outplane_kpoints(double *kpoints, double *weights, unsigned int nkpt, double * &outplane_kpt, double * &outplane_weights, unsigned int &n_outplane_kpt, double *t) {
	n_outplane_kpt = 0;
	for(unsigned int i=0u; i<nkpt; ++i) {
		if(abs(vecdot(&kpoints[3u*i],t)) >= TOL_PLANE) {
			n_outplane_kpt++;
		}
	}
	outplane_kpt = new double[3u*n_outplane_kpt];
	outplane_weights = new double[n_outplane_kpt];
	unsigned int j = 0u;
	for(unsigned int i=0u; i<nkpt; ++i) {
		if(abs(vecdot(&kpoints[3u*i],t)) >= TOL_PLANE) {
			outplane_kpt[3u*j] = kpoints[3u*i];
			outplane_kpt[3u*j+1u] = kpoints[3u*i+1u];
			outplane_kpt[3u*j+2u] = kpoints[3u*i+2u];
			outplane_weights[j] = weights[i];
			j++;
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
void calcGns(DynMat &dynmat, double t[3], double m[3]) {
/*	cout << dynmat.getNions() << " " << n << "\n";
	Matrix rot;
	Matrix rotT;
	getAcousticRots(dynmat, rot, rotT);
*/
	PolarIntegrator pi(dynmat, t, m, &one_on_k_2_total_rot, &i_on_k_total_rot, &k0_total_rot);

	//cerr << "k^-2 expansion nmax=" << L_MAX_N << endl;
	pi.calcGnExp(-2, L_MAX_N, gn_m2);
	/*
	for(unsigned int n=0u; n<L_MAX_N/2; ++n) {
		cerr << "n=" << 2*n << "\n";
		cerr << gn_m2[n];
	}
	cerr.flush();
	*/
	cerr << "k^-1 expansion nmax=" << L_MAX_N << endl;
	pi.calcGnExp(-1, L_MAX_N, gn_m1);
	for(unsigned int n=0u; n<L_MAX_N/2; ++n) {
		cerr << "n=" << 2*n+1 << "\n";
		cerr << gn_m1[n];
	}
	cerr.flush();
	cerr << "k^0 expansion nmax=" << L_MAX_N << endl;
	pi.calcGnExp(0, L_MAX_N, gn_0);
	for(unsigned int n=0u; n<L_MAX_N/2; ++n) {
		cerr << "n=" << 2*n << "\n";
		cerr << gn_0[n];
	}
	cerr.flush();
}

/* Semi-continuum matrix */
ZMatrix *Gsc;

/* calculate the semi-continuum correction
 * dynmat: Dynamical matrix
 * uc: crystal unit cell
 * inkpoints: kpoints in the dislocation slab plane perp to t
 * ninkpt: number of kpoints
 */
void doSemiCont(DynMat &dynmat, UnitCell &uc, double *inkpoints, unsigned int ninkpt) {
	InvExpPtrs fcptrs;

	fcptrs.one_on_k2 = &one_on_k_2_total_rot;
	fcptrs.i_on_k = &i_on_k_total_rot;
	fcptrs.k0 = &k0_total_rot;

	fcptrs.XiInvAA = &XiInvAA;
	fcptrs.XiInvAO = &XiInvAO;
	fcptrs.XiInvOO = &XiInvOO;

	fcptrs.i_on_k_sub = &i_on_k_sub;
	fcptrs.k0_sub = &k0_sub;

	SemiCont sc(dynmat, fcptrs, inkpoints, ninkpt, uc);
	sc.calcSemiCont(Gsc);
	std::cerr << "SemiCont done..." << std::endl;
	/*
	for(int i = 0; i < ninkpt; ++i) {
		std::cerr << Gsc[i];
	}
	std::cerr.flush();
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
Matrix getGbr(DynMat &dynmat, double t[3], double m[3], double R[3], UnitCell &uc, double *inkpoints, unsigned int ninkpt, double *outkpoints, unsigned int noutkpt) {
	PolarIntegrator pi(dynmat, t, m, &one_on_k_2_total_rot, &i_on_k_total_rot, &k0_total_rot);

	Matrix gbrm2 = pi.calcG_b_R(gn_m2, -2, L_MAX_N, uc.getKmax(), R, uc.getVolume());
	Matrix gbrm1 = pi.calcG_b_R(gn_m1, -1, L_MAX_N, uc.getKmax(), R, uc.getVolume());
	Matrix gbr0 = pi.calcG_b_R(gn_0, 0, L_MAX_N, uc.getKmax(), R, uc.getVolume());
	cerr << "R = (" << R[0] << "," << R[1] << "," << R[2] << ")\n";
	cerr << "g(-2)(R) =\n" << gbrm2;
	cerr << "g(-1)(R) =\n" << gbrm1;
	cerr << "g( 0)(R) =\n" << gbr0;



	/* calculate the contribution from the in plane semi-continuum piece
	 * Assume time reversal symmetry k -> -k. G_{(sc)} -> G^{*}_{(sc)}.
	 */
	Matrix gR_sc(3u*dynmat.getNions(), 3u*dynmat.getNions());
	for(unsigned int i=0u; i<ninkpt; ++i) {
		double *curr_k = inkpoints + 3u*i;
		double kdotR = vecdot(curr_k,R);
		gR_sc += std::cos(kdotR)*(Gsc[i].getReal()) + std::sin(kdotR)*(Gsc[i].getImag());
	}
	gR_sc *= 1.0/(ninkpt+noutkpt);

	std::cerr << "Semicontinuum:\n" << gR_sc;


	Matrix gR_out(3u*dynmat.getNions(), 3u*dynmat.getNions());
	/* Outer grid: just invert D(k) */
	for(unsigned int i=0u; i<noutkpt; ++i) {
		double *curr_k = outkpoints + 3u*i;
		ZMatrix Dft(3u*dynmat.getNions(), 3u*dynmat.getNions());
		dynmat.getFourierTransformRot(curr_k, Dft);
		Dft.invert();
		double kdotR = vecdot(curr_k,R);
		gR_out += std::cos(kdotR)*Dft.getReal() + std::sin(kdotR)*Dft.getImag();
	}
	gR_out *= 1.0/(ninkpt+noutkpt);

	//std::cout << "Parallel planes: " << gR_out;
	//std::cout.flush();

	return gbrm2 + gbrm1 + gbr0 + gR_sc + gR_out;
}

/* load in the dislocation description file from dislstream */
void load_disl(std::istream &dislstream, const UnitCell &uc, double t[3], double m[3]) {
	std::string s;
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
	int m_unit[3];
	iss >> m_unit[0] >> m_unit[1] >> m_unit[2];
	Matrix avec = uc.getAvec();
	double scale = uc.getScale();
	m[0] = (m_unit[0]*avec.val(0,0) + m_unit[1]*avec.val(1,0) + m_unit[2]*avec.val(2,0))*scale;
	m[1] = (m_unit[0]*avec.val(0,1) + m_unit[1]*avec.val(1,1) + m_unit[2]*avec.val(2,1))*scale;
	m[2] = (m_unit[0]*avec.val(0,2) + m_unit[1]*avec.val(1,2) + m_unit[2]*avec.val(2,2))*scale;

	t[0] = (t_unit[0]*avec.val(0,0) + t_unit[1]*avec.val(1,0) + t_unit[2]*avec.val(2,0))*scale;
	t[1] = (t_unit[0]*avec.val(0,1) + t_unit[1]*avec.val(1,1) + t_unit[2]*avec.val(2,1))*scale;
	t[2] = (t_unit[0]*avec.val(0,2) + t_unit[1]*avec.val(1,2) + t_unit[2]*avec.val(2,2))*scale;
}

int main(int argc, char *argv[]) {
	if(argc != 6) {
		cerr << "Usage: " << argv[0] << " dynamical_matrix kpoints cellvec dislocfile rvecfile"
		     << endl;
		return(-1);
	}
	double *kpoints, *weights;
	unsigned int num_kpoints;
	DynMat dynmat;
	ifstream dmstream(argv[1]);
	ifstream kpstream(argv[2]);
	ifstream cellstream(argv[3]);
	ifstream dislstream(argv[4]);
	ifstream rvecstream(argv[5]);

	Matrix avec(3u,3u);

	cout.precision(16);
	cout << scientific;
	cout.setf(cout.showpos);
	cerr.precision(16);
	cerr << scientific;
	cerr.setf(cerr.showpos);

	/* load our input data */

	dynmat.load(dmstream);

	UnitCell uc(cellstream);

	load_kpoints(kpstream, num_kpoints, kpoints, weights, uc);

	double *Rpoints;
	unsigned int num_rpoints;
	load_Rpoints(rvecstream, num_rpoints, Rpoints, uc);


	double t[3];
	double m[3];


	load_disl(dislstream, uc, t, m);

	rvecstream.close();
	dislstream.close();
	cellstream.close();
	dmstream.close();
	kpstream.close();

	unsigned int n_inplane_kpt, n_outplane_kpt;
	double *inplane_kpt, *outplane_kpt, *inplane_weights, *outplane_weights;

	// split kpoints into in threading plane, and out of threading plane
	inplane_kpoints(kpoints, weights, num_kpoints, inplane_kpt, inplane_weights, n_inplane_kpt, t);
	outplane_kpoints(kpoints, weights, num_kpoints, outplane_kpt, outplane_weights, n_outplane_kpt, t);

	//calculate polar expansion of the lattice Green function divergent
	//terms
	calcGns(dynmat, t, m);
	Gsc = new ZMatrix[n_inplane_kpt];

	//calculate semi-continuum piece
	doSemiCont(dynmat, uc, inplane_kpt, n_inplane_kpt);

	Matrix rot(dynmat.getRot());
	Matrix rotT(dynmat.getRotT());

	//rotate back from acoustic/optical to cartesian basis
	for(double *R = Rpoints; R < Rpoints + 3u*num_rpoints; R += 3) {
		Matrix gbr = rotT*getGbr(dynmat, t, m, R, uc, inplane_kpt, n_inplane_kpt, outplane_kpt, n_outplane_kpt)*rot;

		std::cout << R[0] << " " << R[1] << " " << R[2] << "\n" << gbr;
	}
	std::cout.flush();

	delete [] Gsc;
	delete [] Rpoints;
	delete [] kpoints;
	delete [] weights;
	delete [] inplane_kpt;
	delete [] outplane_kpt;
	delete [] inplane_weights;
	delete [] outplane_weights;

	return 0;
}
