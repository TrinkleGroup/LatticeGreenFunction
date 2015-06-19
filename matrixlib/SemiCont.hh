#ifndef _SEMI_CONT_H_
#define _SEMI_CONT_H_

#include <cmath>
#include <fstream>
#include "DynMat.hh"
#include "Matrix.hh"
#include "ZMatrix.hh"
#include "UnitCell.hh"

/* function pointers to:
 * 1/k^2 correction
 * i/k   correction
 * k^0   correction
 * Leading order pieces
 * next leading order pieces for subtraction
 */
struct InvExpPtrs {
	Matrix (*one_on_k2)(DynMat &, double[3]);
	Matrix (*i_on_k)(DynMat &, double[3]);
	Matrix (*k0)(DynMat &, double[3]);
	Matrix (*XiInvAA)(DynMat &, double[3]);
	Matrix (*XiInvAO)(DynMat &, double[3]);
	Matrix (*XiInvOO)(DynMat &, double[3]);
	Matrix (*i_on_k_sub)(DynMat &, double[3]);
	Matrix (*k0_sub)(DynMat &, double[3]);
};


class SemiCont {
private:
	static const double SMALLK_TOL = 1e-12;
	/* log(SMALLK_TOL): */
	static const double LOG_SMALLK_TOL = -32.236191301916641;
	static const int SMALLK_NMAX = 10;
	double *kpoints;  /* Nx3 list of kpoints for inversion */
	double kmax;       /* maximum k (sphere inscribed inside BZ */
	double ksmall;     /* switch to small k expansion below this k */
	int nkpt;

	UnitCell ucell;    /* Unit cell: lattice vectors, reciprocal lattice
			    * vectors, cell volume, etc.
			    */

	InvExpPtrs fcptrs;

	DynMat *dynmat;

	SemiCont();
	void load_kpoints(std::ifstream &ifs);
public:
	SemiCont(DynMat &dynmat, InvExpPtrs fcptrs, double *kpoints, int nkpt, UnitCell &ucell);
	SemiCont(DynMat &dynmat, InvExpPtrs fcptrs, std::ifstream &kpt_file, UnitCell &ucell);
	~SemiCont();

	void calcSemiCont(ZMatrix *gsc);
	ZMatrix directSub(double *kpt);
	ZMatrix smallKtaylor(double *kpt);
	ZMatrix smallKexp(double *kpt);
};

#endif
