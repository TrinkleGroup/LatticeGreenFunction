#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include "gsl_complex_op.hh"
#include "UnitCell.hh"
#include "Matrix.hh"
#include "VecMath.hh"

using namespace std;

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
		cerr << "Usage: " << argv[0] << " cellvec dislocfile T M N"
		     << endl;
		return(-1);
	}
	int num_t;
	int num_m;
	int num_n;
	ifstream cellstream(argv[1]);
	ifstream dislstream(argv[2]);
	num_t = atoi(argv[3]);
	num_m = atoi(argv[4]);
	num_n = atoi(argv[5]);

	cout.precision(16);
	cout << scientific;
	cout.setf(cout.showpos);
	cerr.precision(16);
	cerr << scientific;
	cerr.setf(cerr.showpos);

	UnitCell uc(cellstream);

	double t[3];
	double m[3];
	double n[3];


	load_disl(dislstream, uc, t, m);

	dislstream.close();
	cellstream.close();

	crossprod(t,m,n);

	//grid generation:
	//\vec{k} = 2\pi/|t|^2 \vec{t} t/T + 2\pi/|m|^2 \vec{m} m/M + 2\pi/|n|^2 \vec{n} n/N
	//t = -T+1..T-1, m = -M+1..M-1, n=-N+1..N-1
	//num kpoints = (2T-1)(2M-1)(2N-1)

	double *kpoints;
	int num_kpt = (2*num_t - 1)*(2*num_m - 1)*(2*num_n - 1);
	kpoints = new double[3*num_kpt];

	double ts = 2.0*M_PI/vecmag2(t)/num_t;
	double ms = 2.0*M_PI/vecmag2(m)/num_m;
	double ns = 2.0*M_PI/vecmag2(n)/num_n;

	std::cout << num_kpt << "\n";
	int q = 0;
	for(int tx = -num_t+1; tx<num_t; ++tx) {
		for(int mx = -num_m+1; mx<num_m; ++mx) {
			for(int nx = -num_n+1; nx<num_n; ++nx) {
				kpoints[3*q] = tx*ts*t[0] + mx*ms*m[0] + nx*ns*n[0];
				kpoints[3*q+1] = tx*ts*t[1] + mx*ms*m[1] + nx*ns*n[1];
				kpoints[3*q+2] = tx*ts*t[2] + mx*ms*m[2] + nx*ns*n[2];
				std::cout << kpoints[3*q] << " "
					  << kpoints[3*q+1] << " "
					  << kpoints[3*q+2] << "\n";
				++q;
			}
		}
	}
	std::cout.flush();

	return 0;
}
