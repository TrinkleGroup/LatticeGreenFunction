/* File: DynMat.cc
 * Dyanmical matrix manipulation functions
 */
#include <iostream>
#include <cmath>
#include <cstring>
#include <cassert>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf_gamma.h>
#include "gsl_complex_op.hh"
#include "DynMat.hh"
#include "Matrix.hh"
#include "ZMatrix.hh"
#include "VecMath.hh"

#define A_TOL 1.0e-8 /* Acoustic mode tolerance */

/* calculate the phonon projectors for k=0. */
void DynMat::calcProjectors() {
	if(ap != NULL) delete ap;
	if(op != NULL) delete op;
	if(apt != NULL) delete apt;
	if(opt != NULL) delete opt;
	if(apf != NULL) delete ap;
	if(opf != NULL) delete op;
	if(aptf != NULL) delete apt;
	if(optf != NULL) delete opt;
	Matrix m(*(mat[0]));
	double **evecs;
	double *evals;
	for(int i=1; i < numr; ++i) {
		m += *(mat[i]);
	}
	m.symm_eigen(evecs, evals);
	/* eigenvalues should be sorted by absolute value
	 * Acoustic modes should be the first 3 modes
	 * But we'll count */
	unsigned int acoustic_modes = 0u;
	for(unsigned int i=0u; i<3u*nions; ++i) {
		if(fabs(evals[i]) < A_TOL) {
			acoustic_modes++;
		}
	}
	std::cerr << "Acoustic Modes: " << acoustic_modes << std::endl;
	std::cerr << "Eigenvalue, Eigenvectors: " << std::endl;
	for(unsigned int i=0u; i<3u*nions; ++i) {
		std::cerr << evals[i] << ",";
		for(unsigned int j=0u; j<3u*nions; ++j) {
			std::cerr << " " << evecs[i][j];
		}
		std::cerr << std::endl;
	}
	std::cerr << std::endl;
	if(acoustic_modes != 3u) {
		std::cerr << "ERROR: There are not 3 acoustic modes."
		          << std::endl;
	}
	assert(acoustic_modes == 3u);
	/* acoustic projectors to 3x3 matrix */
	ap = new Matrix(evecs, acoustic_modes, 3u*nions);
	/* pointer arithmetic, start from the first optical eigenvector */
	/* optical projectors  to (3N-3)x(3N-3) matrix */
	op = new Matrix(evecs+acoustic_modes, 3u*nions-acoustic_modes, 3u*nions);
	/* transpose (inverse) of acoustic projectors */
	apt = new Matrix(ap->transpose());
	/* transpose (inverse) of optical projectors */
	opt = new Matrix(op->transpose());
	delete [] evals;
	for(unsigned int i=0u; i<3u*nions; ++i) {
		delete [] evecs[i];
	}
	delete [] evecs;

	/* acoustical and optical projectors for the full
	 * 3Nx3N geometry */
	apf = new Matrix(3u*nions,3u);
	opf = new Matrix(3u*nions,3u*(nions-1));

	apf->val(0,0) = 1.0;
	apf->val(1,1) = 1.0;
	apf->val(2,2) = 1.0;

	for(int i = 0; i < 3*(nions-1); ++i) {
		opf->val(i+3,i) = 1.0;
	}
	aptf = new Matrix(apf->transpose());
	optf = new Matrix(opf->transpose());

}

/* build rotation matrix for the acoustic mode/ optical mode quadrants */
void DynMat::calcAcousticRots() {
	Matrix ap = getAcousticProjector();
	Matrix op = getOpticalProjector();
	rot = new Matrix(ap.getNumColumns(), ap.getNumColumns());
	for(unsigned int i = 0u; i < 3u; ++i) {
		for(unsigned int j = 0u; j < ap.getNumColumns(); ++j) {
			(*rot)[i][j] = double(ap.val(i,j));
		}
	}
	for(unsigned int i = 0u; i < op.getNumRows(); ++i) {
		for(unsigned int j = 0u; j < op.getNumColumns(); ++j) {
			(*rot)[i+3u][j] = double(op.val(i,j));
		}
	}

	rotT = new Matrix(rot->transpose());
}

/* default constructor */
DynMat::DynMat() : orderCache(5,(Matrix *)NULL),
                   ap(NULL), op(NULL), apt(NULL), opt(NULL),
                   apf(NULL), opf(NULL), aptf(NULL), optf(NULL),
                   rot(NULL), rotT(NULL),
                   nions(0), numr(0), R(NULL), mat(NULL)
{
	cached_k[0] = 0.0;
	cached_k[1] = 0.0;
	cached_k[2] = 0.0;
}

/* Construct from standard C++ double arrays
 * matrix: The R x 3N x 3N dynamical matrix,
 * pointer indexes are: R, i, j
 * where R labels the lattice vector in R
 * i is a cartesian index and atom index x1,y1,z1,x2,y2,x2,...
 * j is a cartesian index and atom index x1,y1,z1,x2,y2,x2,...
 * nions: number of ions N
 * numr: number of lattice vectors R
 */
DynMat::DynMat(double ***matrix, double *R, unsigned int nions, unsigned int numr)
	: orderCache(5,(Matrix *)NULL), nions(nions),
	  numr(numr), ap(NULL), op(NULL), apt(NULL), opt(NULL),
	  apf(NULL), opf(NULL), aptf(NULL), optf(NULL),
          rot(NULL), rotT(NULL)
{
	cached_k[0] = 0.0;
	cached_k[1] = 0.0;
	cached_k[2] = 0.0;
	this->R = new double[3u*numr];
	this->mat = new Matrix*[numr];
	memcpy(this->R, R, sizeof(double)*3u*numr);
	for(unsigned int i=0u; i<numr; ++i) {
		this->mat[i] = new Matrix(matrix[i], 3u*nions, 3u*nions);
	}
	calcProjectors();
	calcAcousticRots();
}

/* construct from Matrix wrapper class
 * matrix: The R x 3N x 3N dynamical matrix,
 * pointer index is for: R vectors
 * Each Matrix is 3N x 3N with indexes:
 * i is a cartesian index and atom index x1,y1,z1,x2,y2,x2,...
 * j is a cartesian index and atom index x1,y1,z1,x2,y2,x2,...
 * nions: number of ions N
 * numr: number of lattice vectors R
 */
DynMat::DynMat(Matrix *matrix, double *R, unsigned int nions, unsigned int numr)
	: orderCache(5,(Matrix *)NULL), nions(nions),
	  numr(numr), ap(NULL), op(NULL), apt(NULL), opt(NULL),
	  apf(NULL), opf(NULL), aptf(NULL), optf(NULL),
          rot(NULL), rotT(NULL)
{
	cached_k[0] = 0.0;
	cached_k[1] = 0.0;
	cached_k[2] = 0.0;
	this->R = new double[3u*numr];
	this->mat = new Matrix*[numr];
	memcpy(this->R, R, sizeof(double)*3u*numr);
	for(unsigned int i=0u; i<numr; ++i) {
		this->mat[i] = new Matrix(matrix[i]);
	}
	calcProjectors();
	calcAcousticRots();
}

/* Memory clean-up */
void DynMat::clean_up() {
	if(mat != NULL) {
		for(unsigned int i=0u; i<numr; ++i) {
			if(mat[i] != NULL) delete mat[i];
		}
		delete [] mat;
		mat = NULL;
	}
	if(R != NULL) {
		delete [] R;
		R = NULL;
	}
	if(ap != NULL)
		delete ap;
	if(op != NULL)
		delete op;
	if(apt != NULL)
		delete apt;
	if(opt != NULL)
		delete opt;
	if(apf != NULL)
		delete apf;
	if(opf != NULL)
		delete opf;
	if(aptf != NULL)
		delete aptf;
	if(optf != NULL)
		delete optf;
	if(rot != NULL)
		delete rot;
	if(rotT != NULL)
		delete rotT;
	for(std::vector<Matrix *>::iterator it = orderCache.begin();
	    it < orderCache.end(); ++it)
	{
		if(*it != NULL) {
			delete *it;
			*it = NULL;
		}
	}
}

/* destructor */
DynMat::~DynMat() {
	clean_up();
}

//Load Dynamical matrix from an input stream is:
//Format:
//Description
//num_r_vectors(M) num_ions(N)
//R1_1 R1_2 R1_3 D(R1)_1,1 D_1,2 .. D(R1)_1,3N .. D(R1)_3N,3N
//R2_1 R2_2 R2_3 D(R2)_1,1 D_1,2 .. D(R2)_1,3N .. D(R2)_3N,3N
//  ...
//RM_1 RM_2 RM_3 D(RM)_1,1 D_1,2 .. D(RM)_1,3N .. D(RM)_3N,3N
void DynMat::load(std::istream &is) {
	clean_up();
	char buf[256];
	is.getline(buf, 256);
	description = buf;
	is >> numr >> nions;
	mat = new Matrix*[numr];
	R = new double[3u*numr];
	double *temp = new double[9u*nions*nions];
	for(unsigned int i=0u;i<numr;++i) {
		is >> R[i*3u] >> R[i*3u + 1u] >> R[i*3u + 2u];
		for(unsigned int j=0u;j<9u*nions*nions;++j) {
			is >> temp[j];
		}
		mat[i] = new Matrix(temp, 3u*nions, 3u*nions);
	}
	delete [] temp;
	calcProjectors();
	calcAcousticRots();
}

/* Calculate the k^n Taylor expansion term in the
 * Fourier space Dynamical matrix
 */
const Matrix &DynMat::getOrder(int n, double k[3]) {
	// check to see if this k vector is cached
	if(cached_k[0] != k[0] || cached_k[1] != k[1] || cached_k[2] != k[2]) {
		// not cached, clear the cache and add the new vector.
		for(std::vector<Matrix *>::iterator it = orderCache.begin();
		    it < orderCache.end(); ++it)
		{
			if(*it != NULL) {
				delete *it;
				*it = NULL;
			}
		}
		cached_k[0] = k[0];
		cached_k[1] = k[1];
		cached_k[2] = k[2];
	}
	//Is the cache large enough for the k^n term?
	if(orderCache.size() > n) {
		//yes, and found entry in the cache
		if(orderCache[n] != NULL) {
			return *orderCache[n];
		}
	} else {
		//grow the cache to accomodate
		orderCache.resize(n+1, (Matrix *)NULL);
	}
	// Calculate the k^n term and add to the cache
	Matrix *ret = new Matrix(3u*nions, 3u*nions);
	for(unsigned int i=0u; i<numr; ++i) {
		double coeff = gsl_pow_int(vecdot(R+3u*i,k),n)/gsl_sf_fact(n);
		*ret += (*(mat[i]))*coeff;
	}
	orderCache[n] = ret;
	return *orderCache[n];
}

const Matrix &DynMat::getAcousticProjector() const {
	return *ap;
}

const Matrix &DynMat::getOpticalProjector() const {
	return *op;
}

const Matrix &DynMat::getAcousticProjectorT() const {
	return *apt;
}

const Matrix &DynMat::getOpticalProjectorT() const {
	return *opt;
}

const Matrix &DynMat::getAcousticFullP() const {
	return *apf;
}

const Matrix &DynMat::getOpticalFullP() const {
	return *opf;
}

const Matrix &DynMat::getAcousticFullPT() const {
	return *aptf;
}

const Matrix &DynMat::getOpticalFullPT() const {
	return *optf;
}

const Matrix &DynMat::getRot() const {
	return *rot;
}

const Matrix &DynMat::getRotT() const {
	return *rotT;
}

/* calculate the full real and imaginary pieces of the Fourier transform
 * of the dynamical matrix at k
 */
void DynMat::getFourierTransform(double k[3], Matrix &realFt, Matrix &imFt) const {
	/* Zero out result matrices */
	realFt.setZero();
	imFt.setZero();

	for(unsigned int i = 0; i<numr; ++i) {
		double kr = vecdot(k,R + 3u*i);
		realFt += (*(mat[i])) * std::cos(kr);
		imFt += (*(mat[i])) * std::sin(kr);
	}
}

/* calculate the full real and imaginary pieces of the Fourier transform
 * of the dynamical matrix at k, complex matrix output
 */
void DynMat::getFourierTransform(double k[3], ZMatrix &ft) const {
	/* Zero out result matrices */
	ft.setZero();

	for(unsigned int i = 0; i<numr; ++i) {
		double kr = vecdot(k,R + 3u*i);
		gsl_complex ikr;
		GSL_SET_COMPLEX(&ikr, 0.0, kr);

		ft += (*(mat[i])) * gsl_complex_exp(ikr);
	}
}

/* calculate the full real and imaginary pieces of the Fourier transform
 * of the dynamical matrix at k in the rotated Acoustic/Optical basis
 */
void DynMat::getFourierTransformRot(double k[3], Matrix &realFt, Matrix &imFt) const {
	/* Zero out result matrices */
	realFt.setZero();
	imFt.setZero();

	for(unsigned int i = 0; i<numr; ++i) {
		double kr = vecdot(k,R + 3u*i);
		//realFt += *rot*(*(mat[i]))* *rotT * std::cos(kr);
		//imFt += *rot*(*(mat[i])) * *rotT * std::sin(kr);
		realFt += *(mat[i]) * std::cos(kr);
		imFt += *(mat[i]) * std::sin(kr);
	}
	realFt = *rot * realFt * *rotT;
	imFt = *rot * imFt * *rotT;
}

/* calculate the full real and imaginary pieces of the Fourier transform
 * of the dynamical matrix at k in the rotated Acoustic/Optical basis
 * complex matrix output
 */
void DynMat::getFourierTransformRot(double k[3], ZMatrix &ft) const {
	/* Zero out result matrices */
	ft.setZero();
	Matrix zero(3u*getNions(), 3u*getNions());

	for(unsigned int i = 0; i<numr; ++i) {
		double kr = vecdot(k,R + 3u*i);
		gsl_complex ikr;
		GSL_SET_COMPLEX(&ikr, 0.0, kr);

		//ft += *rot * (*(mat[i])) * *rotT * gsl_complex_exp(ikr);
		ft += *(mat[i]) * gsl_complex_exp(ikr);
	}
	ft = ZMatrix(*rot,zero) * ft * ZMatrix(*rotT,zero);
}

/* calculate the small k expansion real and imaginary pieces of the
 * Fourier transform of the dynamical matrix at k in the
 * rotated Acoustic/Optical basis. complex matrix output.
 */
void DynMat::getSmallFourierTransform(double k[3], ZMatrix &ft) const {

	Matrix ap = getAcousticProjector();
	Matrix apt = getAcousticProjectorT();
	Matrix op = getOpticalProjector();
	Matrix opt = getOpticalProjectorT();

	ZMatrix acac(3u,3u);
	ZMatrix acop(3u,3u*(getNions()-1u));
	ZMatrix opac(3u*(getNions()-1u),3u);
	ZMatrix opop(3u*(getNions()-1u),3u*(getNions()-1u));

	for(unsigned int i = 0; i<numr; ++i) {
		double kr = vecdot(k,R + 3u*i);
		gsl_complex ikr;
		GSL_SET_COMPLEX(&ikr, 0.0, kr);

		gsl_complex expk2 = 1.0 + ikr - 0.5*kr*kr - gsl_complex_exp(ikr);
		gsl_complex expk1 = 1.0 + ikr - gsl_complex_exp(ikr);
		gsl_complex expk0 = 1.0 - gsl_complex_exp(ikr);

		acac += ap*(*(mat[i]))*apt * expk2;
		acop += ap*(*(mat[i]))*opt * expk1;
		opac += op*(*(mat[i]))*apt * expk1;
		opop += op*(*(mat[i]))*opt * expk0;
	}
	for(int i = 0; i < 3; ++i) {
		for(int j=0; j<3; ++j) {
			ft.val(i,j) = acac.val(i,j);
		}
	}
	for(int i = 0; i < 3; ++i) {
		for(int j=3; j<3*getNions(); ++j) {
			ft.val(i,j) = acop.val(i,j-3);
		}
	}
	for(int i = 3; i < 3*getNions(); ++i) {
		for(int j=0; j<3; ++j) {
			ft.val(i,j) = opac.val(i-3,j);
		}
	}
	for(int i = 3; i < 3*getNions(); ++i) {
		for(int j=3; j<3*getNions(); ++j) {
			ft.val(i,j) = opop.val(i-3,j-3);
		}
	}
}

unsigned int DynMat::getNions() const {
	return nions;
}
unsigned int DynMat::getNumr() const {
	return numr;
}
