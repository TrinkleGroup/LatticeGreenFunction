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

#define A_TOL 1.0e-7 /* Acoustic mode tolerance */

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
	 * Acoustic modes should be the first CARTDIM modes
	 * But we'll count */
	unsigned int acoustic_modes = 0u;
	for(unsigned int i=0u; i<CARTDIM*nions; ++i) {
		if(fabs(evals[i]) < A_TOL) {
			acoustic_modes++;
		}
	}
	std::cerr << "D(k=0): " << m << std::endl;
	std::cerr << "Acoustic Modes: " << acoustic_modes << std::endl;
	std::cerr << "Eigenvalue, Eigenvectors: " << std::endl;
	for(unsigned int i=0u; i<CARTDIM*nions; ++i) {
		std::cerr << evals[i] << ",";
		for(unsigned int j=0u; j<CARTDIM*nions; ++j) {
			std::cerr << " " << evecs[i][j];
		}
		std::cerr << std::endl;
	}
	std::cerr << std::endl;
	if(acoustic_modes != CARTDIM) {
		std::cerr << "ERROR: There are not " << CARTDIM << " acoustic modes."
		          << std::endl;
	}
	assert(acoustic_modes == CARTDIM);

	if(nions == 1) { /* no optical space if only one ion */
		ap = new Matrix(Matrix::eye(CARTDIM));
		apt = new Matrix(Matrix::eye(CARTDIM));
		apf = new Matrix(Matrix::eye(CARTDIM));
		aptf = new Matrix(Matrix::eye(CARTDIM));
		op = NULL;
		opt = NULL;
		opf = NULL;
		optf = NULL;
	} else {
		/* acoustic projectors to CARTDIM x CARTDIM matrix */
		ap = new Matrix(evecs, acoustic_modes, CARTDIM*nions);
		/* pointer arithmetic, start from the first optical eigenvector */
		/* optical projectors  to (CARTDIM N-CARTDIM)x(CARTDIM N-CARTDIM) matrix */
		op = new Matrix(evecs+acoustic_modes, CARTDIM*nions-acoustic_modes, CARTDIM*nions);
		/* transpose (inverse) of acoustic projectors */
		apt = new Matrix(ap->transpose());
		/* transpose (inverse) of optical projectors */
		opt = new Matrix(op->transpose());

		/* acoustical and optical projectors for the full
		 * CARTDIM NxCARTDIM N geometry */
		apf = new Matrix(CARTDIM*nions,CARTDIM);
		opf = new Matrix(CARTDIM*nions,CARTDIM*(nions-1));

		for(unsigned int i=0ul; i<CARTDIM; ++i) {
			apf->val(i,i) = 1.0;
		}

		for(int i = 0; i < CARTDIM*(nions-1); ++i) {
			opf->val(i+CARTDIM,i) = 1.0;
		}
		aptf = new Matrix(apf->transpose());
		optf = new Matrix(opf->transpose());
	}
	delete [] evals;
	for(unsigned int i=0u; i<CARTDIM*nions; ++i) {
		delete [] evecs[i];
	}
	delete [] evecs;


}

/* build rotation matrix for the acoustic mode/ optical mode quadrants */
void DynMat::calcAcousticRots() {
	if(nions == 1) {
		rot = new Matrix(Matrix::eye(CARTDIM));
		rotT = new Matrix(Matrix::eye(CARTDIM));
		return;
	}
	Matrix ap = getAcousticProjector();
	Matrix op = getOpticalProjector();
	rot = new Matrix(ap.getNumColumns(), ap.getNumColumns());
	for(unsigned int i = 0u; i < CARTDIM; ++i) {
		for(unsigned int j = 0u; j < ap.getNumColumns(); ++j) {
			(*rot)[i][j] = double(ap.val(i,j));
		}
	}
	for(unsigned int i = 0u; i < op.getNumRows(); ++i) {
		for(unsigned int j = 0u; j < op.getNumColumns(); ++j) {
			(*rot)[i+CARTDIM][j] = double(op.val(i,j));
		}
	}

	rotT = new Matrix(rot->transpose());
}

/* default constructor */
DynMat::DynMat() : orderCache(5,(Matrix *)NULL),
                   ap(NULL), op(NULL), apt(NULL), opt(NULL),
                   apf(NULL), opf(NULL), aptf(NULL), optf(NULL),
                   rot(NULL), rotT(NULL),
                   nions(0), numr(0), R(NULL), mat(NULL), uc(NULL)
{
	cached_k = new double[CARTDIM];
	for(int idx = 0; idx < CARTDIM; ++idx) {
		cached_k[idx] = 0.0;
	}
}

/* Construct from standard C++ double arrays
 * matrix: The R x CARTDIM N x CARTDIM N dynamical matrix,
 * pointer indexes are: R, i, j
 * where R labels the lattice vector in R
 * i is a cartesian index and atom index x1,y1,z1,x2,y2,x2,...
 * j is a cartesian index and atom index x1,y1,z1,x2,y2,x2,...
 * nions: number of ions N
 * numr: number of lattice vectors R
 */
DynMat::DynMat(double ***matrix, double *R, unsigned int nions, unsigned int numr, UnitCell *uc)
	: orderCache(5,(Matrix *)NULL), nions(nions),
	  numr(numr), ap(NULL), op(NULL), apt(NULL), opt(NULL),
	  apf(NULL), opf(NULL), aptf(NULL), optf(NULL),
          rot(NULL), rotT(NULL)
{
	cached_k = new double[CARTDIM];
	for(unsigned int i=0u; i<CARTDIM; ++i) {
		cached_k[i] = 0.0;
	}
	this->R = new double[CARTDIM*numr];
	this->mat = new Matrix*[numr];
	memcpy(this->R, R, sizeof(double)*CARTDIM*numr);
	for(unsigned int i=0u; i<numr; ++i) {
		this->mat[i] = new Matrix(matrix[i], CARTDIM*nions, CARTDIM*nions);
	}
	calcProjectors();
	calcAcousticRots();
	this->uc = uc;
}

/* construct from Matrix wrapper class
 * matrix: The R x CARTDIM N x CARTDIM N dynamical matrix,
 * pointer index is for: R vectors
 * Each Matrix is CARTDIM N x CARTDIM N with indexes:
 * i is a cartesian index and atom index x1,y1,z1,x2,y2,x2,...
 * j is a cartesian index and atom index x1,y1,z1,x2,y2,x2,...
 * nions: number of ions N
 * numr: number of lattice vectors R
 */
DynMat::DynMat(Matrix *matrix, double *R, unsigned int nions, unsigned int numr, UnitCell *uc)
	: orderCache(5,(Matrix *)NULL), nions(nions),
	  numr(numr), ap(NULL), op(NULL), apt(NULL), opt(NULL),
	  apf(NULL), opf(NULL), aptf(NULL), optf(NULL),
          rot(NULL), rotT(NULL)
{
	cached_k = new double[CARTDIM];
	for(unsigned int i=0u; i<CARTDIM; ++i) {
		cached_k[i] = 0.0;
	}
	this->R = new double[CARTDIM*numr];
	this->mat = new Matrix*[numr];
	memcpy(this->R, R, sizeof(double)*CARTDIM*numr);
	for(unsigned int i=0u; i<numr; ++i) {
		this->mat[i] = new Matrix(matrix[i]);
	}
	calcProjectors();
	calcAcousticRots();

	this->uc = uc;
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
	if(cached_k != NULL) {
		delete [] cached_k;
		cached_k = NULL;
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
//R1_1 R1_2 .. R1_CARTDIM  D(R1)_1,1 D_1,2 .. D(R1)_1,CARTDIM*N .. D(R1)_CARTDIM*N,CARTDIM*N
//R2_1 R2_2 .. R2_CARTDIM  D(R2)_1,1 D_1,2 .. D(R2)_1,CARTDIM*N .. D(R2)_CARTDIM*N,CARTDIM*N
//  ...
//RM_1 RM_2 .. RM_CARTDIM D(RM)_1,1 D_1,2 .. D(RM)_1,CARTDIM*N .. D(RM)_CARTDIM*N,CARTDIM*N
void DynMat::load(std::istream &is, UnitCell &uc) {
	clean_up();

	cached_k = new double[CARTDIM];
	for(unsigned int i=0u; i<CARTDIM; ++i) {
		cached_k[i] = 0.0;
	}

	char buf[256];
	is.getline(buf, 256);
	description = buf;
	is >> numr >> nions;
	mat = new Matrix*[numr];
	R = new double[CARTDIM*numr];
	double *temp = new double[CARTDIM*CARTDIM*nions*nions];
	for(unsigned int i=0u;i<numr;++i) {
		for(unsigned int d = 0u; d < CARTDIM; d++) {
			is >> R[i*CARTDIM + d];
		}
		for(unsigned int j=0u;j<CARTDIM*CARTDIM*nions*nions;++j) {
			is >> temp[j];
		}
		mat[i] = new Matrix(temp, CARTDIM*nions, CARTDIM*nions);
	}
	delete [] temp;
	calcProjectors();
	calcAcousticRots();

	this->uc = &uc;
}

/* Calculate the k^n Taylor expansion term in the
 * Fourier space Dynamical matrix
 */
const Matrix &DynMat::getOrder(int n, double *k) {
	// check to see if this k vector is cached
	if(!vecequal(cached_k, k)) {
		// not cached, clear the cache and add the new vector.
		for(std::vector<Matrix *>::iterator it = orderCache.begin();
		    it < orderCache.end(); ++it)
		{
			if(*it != NULL) {
				delete *it;
				*it = NULL;
			}
		}
		veccpy(cached_k,k);
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
	Matrix *ret = new Matrix(CARTDIM*nions, CARTDIM*nions);
	for(unsigned int i=0u; i<numr; ++i) {
		double coeff = gsl_pow_int(vecdot(R+CARTDIM*i,k),n)/gsl_sf_fact(n);
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
void DynMat::getFourierTransform(double *k, Matrix &realFt, Matrix &imFt) const {
	/* Zero out result matrices */
	realFt.setZero();
	imFt.setZero();

	for(unsigned int i = 0; i<numr; ++i) {
		double kr = vecdot(k,R + CARTDIM*i);
		realFt += (*(mat[i])) * std::cos(kr);
		imFt += (*(mat[i])) * std::sin(kr);
	}
}

/* calculate the full real and imaginary pieces of the Fourier transform
 * of the dynamical matrix at k, complex matrix output
 */
void DynMat::getFourierTransform(double *k, ZMatrix &ft) const {
	/* Zero out result matrices */
	ft.setZero();

	for(unsigned int i = 0; i<numr; ++i) {
		double kr = vecdot(k,R + CARTDIM*i);

		gsl_complex ikr;
		GSL_SET_COMPLEX(&ikr, 0.0, kr);

		ft += (*(mat[i])) * gsl_complex_exp(ikr);
	}
}

/* calculate the full real and imaginary pieces of the Fourier transform
 * of the dynamical matrix at k in the rotated Acoustic/Optical basis
 */
void DynMat::getFourierTransformRot(double *k, Matrix &realFt, Matrix &imFt) const {
	getFourierTransform(k, realFt, imFt);
	realFt = *rot * realFt * *rotT;
	imFt = *rot * imFt * *rotT;
}

/* calculate the full real and imaginary pieces of the Fourier transform
 * of the dynamical matrix at k in the rotated Acoustic/Optical basis
 * complex matrix output
 */
void DynMat::getFourierTransformRot(double *k, ZMatrix &ft) const {
	getFourierTransform(k, ft);
	ft = ZMatrix(*rot) * ft * ZMatrix(*rotT);
}

/* calculate the small k expansion real and imaginary pieces of the
 * Fourier transform of the dynamical matrix at k in the
 * rotated Acoustic/Optical basis. complex matrix output.
 */
void DynMat::getSmallFourierTransform(double *k, ZMatrix &ft) const {
	if(nions == 1) {
		getSmallFourierTransform1a(k, ft);
	} else {
		getSmallFourierTransformMa(k, ft);
	}
}

void DynMat::getSmallFourierTransform1a(double *k, ZMatrix &ft) const {
	ft.setZero();
	for(unsigned int i = 0; i<numr; ++i) {
		double kr = vecdot(k,R + CARTDIM*i);
		gsl_complex ikr;
		GSL_SET_COMPLEX(&ikr, 0.0, kr);
		gsl_complex expk2 = 1.0 + ikr - 0.5*kr*kr - gsl_complex_exp(ikr);

		ft += (*(mat[i]))*expk2;
	}
}

void DynMat::getSmallFourierTransformMa(double *k, ZMatrix &ft) const {
	Matrix ap = getAcousticProjector();
	Matrix apt = getAcousticProjectorT();
	Matrix op = getOpticalProjector();
	Matrix opt = getOpticalProjectorT();

	ZMatrix acac(CARTDIM,CARTDIM);
	ZMatrix acop(CARTDIM,CARTDIM*(getNions()-1u));
	ZMatrix opac(CARTDIM*(getNions()-1u),CARTDIM);
	ZMatrix opop(CARTDIM*(getNions()-1u),CARTDIM*(getNions()-1u));

	for(unsigned int i = 0; i<numr; ++i) {
		double kr = vecdot(k,R + CARTDIM*i);
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
	for(int i = 0; i < CARTDIM; ++i) {
		for(int j=0; j<CARTDIM; ++j) {
			ft.val(i,j) = acac.val(i,j);
		}
	}
	for(int i = 0; i < CARTDIM; ++i) {
		for(int j=CARTDIM; j<CARTDIM*getNions(); ++j) {
			ft.val(i,j) = acop.val(i,j-CARTDIM);
		}
	}
	for(int i = CARTDIM; i < CARTDIM*getNions(); ++i) {
		for(int j=0; j<CARTDIM; ++j) {
			ft.val(i,j) = opac.val(i-CARTDIM,j);
		}
	}
	for(int i = CARTDIM; i < CARTDIM*getNions(); ++i) {
		for(int j=CARTDIM; j<CARTDIM*getNions(); ++j) {
			ft.val(i,j) = opop.val(i-CARTDIM,j-CARTDIM);
		}
	}
}

const UnitCell *DynMat::getUC() const {
	return uc;
}
unsigned int DynMat::getNions() const {
	return nions;
}
unsigned int DynMat::getNumr() const {
	return numr;
}
