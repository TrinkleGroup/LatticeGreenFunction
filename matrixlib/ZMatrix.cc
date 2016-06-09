#include <cassert>
#include <cstring>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sort_vector.h>
#include <iostream>
#include "ZMatrix.hh"

ZMatrix::ZMatrix() :  mat(NULL) { }

ZMatrix::ZMatrix(unsigned int m, unsigned int n) : mat(gsl_matrix_complex_calloc(m,n)) { }

ZMatrix::ZMatrix(gsl_complex **matrix, unsigned int m, unsigned int n) :
	mat(gsl_matrix_complex_alloc(m,n))
{
	for(unsigned int i=0u; i<m; ++i)
		for(unsigned int j=0u; j<n; ++j)
			gsl_matrix_complex_set(mat,i,j,matrix[i][j]);
}

ZMatrix::ZMatrix(gsl_complex *matrix, unsigned int m, unsigned int n) :
	mat(gsl_matrix_complex_alloc(m,n))
{
	for(unsigned int i=0u; i<m; ++i)
		for(unsigned int j=0u; j<n; ++j)
			gsl_matrix_complex_set(mat,i,j,matrix[i*n + j]);
}

ZMatrix::ZMatrix(const ZMatrix &matrix) {
	this->mat = gsl_matrix_complex_alloc(matrix.getNumRows(), matrix.getNumColumns());
	gsl_matrix_complex_memcpy(this->mat, matrix.mat);
}

ZMatrix::ZMatrix(const Matrix &real, const Matrix &imag) {
	assert(real.getNumRows() == imag.getNumRows()
	    && real.getNumColumns() == imag.getNumColumns());
	mat = gsl_matrix_complex_alloc(real.getNumRows(),
				       real.getNumColumns());
	for(int i=0; i<getNumRows(); ++i) {
		for(int j=0; j<getNumColumns(); ++j) {
			GSL_SET_COMPLEX(gsl_matrix_complex_ptr(mat, i, j),real.val(i,j), imag.val(i,j));
		}
	}
}

ZMatrix::ZMatrix(const Matrix &real) {
	mat = gsl_matrix_complex_alloc(real.getNumRows(),
				       real.getNumColumns());
	for(int i=0; i<getNumRows(); ++i) {
		for(int j=0; j<getNumColumns(); ++j) {
			GSL_SET_COMPLEX(gsl_matrix_complex_ptr(mat, i, j),real.val(i,j), 0.0);
		}
	}
}

ZMatrix::~ZMatrix() {
	if(mat != NULL) {
		gsl_matrix_complex_free(mat);
		mat = NULL;
	}
}

gsl_complex &ZMatrix::val(unsigned int i, unsigned int j) {
	return *gsl_matrix_complex_ptr(mat, i, j);
}

gsl_complex ZMatrix::val(unsigned int i, unsigned int j) const {
	return gsl_matrix_complex_get(mat, i, j);
}


//matrix determinant
gsl_complex ZMatrix::det() const {
	//must be a square matrix
	assert(getNumRows() == getNumColumns());
	int s;
	gsl_permutation *p = gsl_permutation_alloc(getNumRows());
	gsl_matrix_complex *lu = gsl_matrix_complex_alloc(getNumRows(),getNumColumns());
	gsl_matrix_complex_memcpy(lu,mat);
			
	int err;
	if((err = gsl_linalg_complex_LU_decomp(lu, p, &s))!=0) {
		std::cerr << "LU decomp problem! " << "err=" << err <<std::endl;
	}

	gsl_complex determ = gsl_linalg_complex_LU_det(lu, s);

	gsl_matrix_complex_free(lu);
	gsl_permutation_free(p);

	return determ;
}

//matrix inversion
void ZMatrix::invert() {
	//must be a square matrix
	assert(getNumRows() == getNumColumns());
	int s;
	gsl_permutation *p = gsl_permutation_alloc(getNumRows());
	gsl_matrix_complex *lu = gsl_matrix_complex_alloc(getNumRows(),getNumColumns());
	gsl_matrix_complex_memcpy(lu,mat);
			
	int err;
	if((err = gsl_linalg_complex_LU_decomp(lu, p, &s))!=0) {
		std::cerr << "LU decomp problem! " << "err=" << err <<std::endl;
	}

	if((err=gsl_linalg_complex_LU_invert(lu, p, mat))!=0) {
		std::cerr << "Inversion problem! " << "err=" << err <<std::endl;
	}

	gsl_matrix_complex_free(lu);
	gsl_permutation_free(p);
}

//matrix transposition
void ZMatrix::transp() {
	//square matrices can do an inplace transpose
	if(getNumRows() == getNumColumns()) {
		gsl_matrix_complex_transpose(mat);
	} else {
		gsl_matrix_complex *m = gsl_matrix_complex_alloc(getNumColumns(), getNumRows());
		gsl_matrix_complex_transpose_memcpy(m, mat);
		gsl_matrix_complex_free(mat);
		mat = m;
	}
}

void ZMatrix::adjInPlace() {
	for(int i=0; i<getNumRows(); ++i) {
		for(int j=0; j<i; ++j) {
			gsl_complex tmp = gsl_complex_conjugate(gsl_matrix_complex_get(mat,i,j));
			gsl_matrix_complex_set(mat, i, j, gsl_complex_conjugate(gsl_matrix_complex_get(mat,j,i)));
			gsl_matrix_complex_set(mat, j, i, tmp);
		}
		gsl_matrix_complex_set(mat, i, i, gsl_complex_conjugate(gsl_matrix_complex_get(mat,i,i)));
	}
}

//matrix adjoint
void ZMatrix::adj() {
	if(isSquare()) {
		adjInPlace();
	} else {
		gsl_matrix_complex *oldmat = mat;
		mat = gsl_matrix_complex_alloc(getNumColumns(), getNumRows());
		for(int i=0; i<getNumRows(); ++i) {
			for(int j=0; j<getNumColumns(); ++j) {
				gsl_matrix_complex_set(mat, i, j, gsl_complex_conjugate(gsl_matrix_complex_get(mat,j,i)));
			}
		}
		gsl_matrix_complex_free(oldmat);
	}
}

/* hermitian matrix eigenvalues: */
void ZMatrix::herm_evals(double *evals) const {
	assert(getNumRows() == getNumColumns());
	
	gsl_vector *ev = gsl_vector_alloc(getNumRows());

	gsl_matrix_complex *zmat = gsl_matrix_complex_alloc(getNumRows(),getNumColumns());
	gsl_matrix_complex_memcpy(zmat, mat);
	gsl_eigen_herm_workspace *w = gsl_eigen_herm_alloc(getNumRows());

	gsl_eigen_herm(zmat, ev, w);

	gsl_sort_vector(ev);

	for(unsigned int i = 0u; i < getNumRows(); ++i) {
		evals[i] = gsl_vector_get(ev,i);
	}

	gsl_eigen_herm_free(w);
	gsl_matrix_complex_free(zmat);
	gsl_vector_free(ev);
}

void ZMatrix::herm_eigen(gsl_complex **&evecs, double *&evals) const {
	assert(getNumRows() == getNumColumns());
	//gsl_eigen_hermv destroys the lower triangle and diagonal,
	//so we'll make a copy
	gsl_matrix_complex *tmpmat = gsl_matrix_complex_alloc(getNumRows(), getNumColumns());
	gsl_matrix_complex_memcpy(tmpmat, mat);

	gsl_vector *eval = gsl_vector_alloc(getNumRows());
	gsl_matrix_complex *evec = gsl_matrix_complex_alloc(getNumRows(),getNumColumns());


	gsl_eigen_hermv_workspace *w = gsl_eigen_hermv_alloc(getNumRows());
	gsl_eigen_hermv(tmpmat, eval, evec, w);

	gsl_eigen_hermv_free(w);

	gsl_matrix_complex_free(tmpmat);


	//sort eigenvalues in ascending order of magnitude
	gsl_eigen_hermv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);

	evecs = new gsl_complex*[getNumColumns()];
	evals = new double[getNumRows()];
	for(unsigned int i = 0u; i < getNumColumns(); ++i) {
		evals[i] = gsl_vector_get(eval,i);
		evecs[i] = new gsl_complex[getNumRows()];
		for(unsigned int j=0u; j<getNumRows(); ++j) {
			evecs[i][j] = gsl_matrix_complex_get(evec,j,i);
		}
	}
	gsl_vector_free(eval);
	gsl_matrix_complex_free(evec);
}

/* condition number estimate */
/* currently slow (computes inverse),
 * should implement a quicker estimate
 */
double ZMatrix::cond1est() const {
	ZMatrix inv = this->inverse();
	double normA = 0.0;
	double normAinv = 0.0;
	/* 1 norm = max_j sum_i A_ij */
	for(int j=0; j<this->getNumColumns(); ++j) {
		double norm1 = 0.0;
		for(int i=0; i<this->getNumRows(); ++i) {
			norm1 += gsl_complex_abs(this->val(i,j));
		}
		if(norm1 > normA) {
			normA = norm1;
		}
	}
	for(int j=0; j<inv.getNumColumns(); ++j) {
		double norm1 = 0.0;
		for(int i=0; i<inv.getNumRows(); ++i) {
			norm1 += gsl_complex_abs(inv.val(i,j));
		}
		if(norm1 > normAinv) {
			normAinv = norm1;
		}
	}

	return normA*normAinv;
}

bool ZMatrix::isSquare() const {
	return getNumRows() == getNumColumns();
}

ZMatrix ZMatrix::inverse() const {
	ZMatrix ret(*this);
	ret.invert();
	return ret;
}

ZMatrix ZMatrix::transpose() const {
	ZMatrix ret(*this);
	ret.transp();
	return ret;
}

ZMatrix ZMatrix::adjoint() const {
	ZMatrix ret(*this);
	ret.adj();
	return ret;
}

ZMatrix ZMatrix::eye(unsigned int n) {
	ZMatrix ret(n,n);
	gsl_matrix_complex_set_identity(ret.mat);
	return ret;
}

unsigned int ZMatrix::getNumRows() const {
	return mat->size1;
}

unsigned int ZMatrix::getNumColumns() const {
	return mat->size2;
}

void ZMatrix::setZero() {
	gsl_matrix_complex_set_zero(mat);
}

void ZMatrix::setReal(const Matrix &m) {
	for(unsigned int i=0u; i<getNumRows(); ++i) {
		for(unsigned int j=0u; j<getNumColumns(); ++j) {
			GSL_SET_REAL(gsl_matrix_complex_ptr(mat,i,j),m.val(i,j));
		}
	}
}

void ZMatrix::setImag(const Matrix &m) {
	for(unsigned int i=0u; i<getNumRows(); ++i) {
		for(unsigned int j=0u; j<getNumColumns(); ++j) {
			GSL_SET_IMAG(gsl_matrix_complex_ptr(mat,i,j),m.val(i,j));
		}
	}
}

Matrix ZMatrix::getReal() const {
	Matrix rp(getNumRows(), getNumColumns());
	for(int i = 0; i < getNumRows(); ++i) {
		for(int j = 0; j < getNumColumns(); ++j) {
			rp.val(i,j) = GSL_REAL(this->val(i,j));
		}
	}
	return rp;
}

Matrix ZMatrix::getImag() const {
	Matrix ip(getNumRows(), getNumColumns());
	for(int i = 0; i < getNumRows(); ++i) {
		for(int j = 0; j < getNumColumns(); ++j) {
			ip.val(i,j) = GSL_IMAG(this->val(i,j));
		}
	}
	return ip;
}

/*
gsl_complex *ZMatrix::operator[](int i) {
	return &(mat->data[i*mat->tda]);
}
*/

ZMatrix ZMatrix::operator+(const ZMatrix& m) const {
	assert(getNumRows() == m.getNumRows() && getNumColumns() == m.getNumColumns());

	ZMatrix ret(*this);

	gsl_matrix_complex_add(ret.mat, m.mat);

	return ret;
}

ZMatrix ZMatrix::operator-(const ZMatrix& m) const {
	assert(getNumRows() == m.getNumRows() && getNumColumns() == m.getNumColumns());

	ZMatrix ret(*this);

	gsl_matrix_complex_sub(ret.mat, m.mat);

	return ret;
}

ZMatrix ZMatrix::operator-() const {
	ZMatrix ret(*this);
	gsl_matrix_complex_scale(ret.mat, GSL_COMPLEX_NEGONE);
	return ret;
}

ZMatrix ZMatrix::operator*(const ZMatrix& m) const {
	assert(getNumColumns() == m.getNumRows());

	ZMatrix ret(getNumRows(),m.getNumColumns());

	/* ret = this * m */

	gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE,
	               this->mat, m.mat, GSL_COMPLEX_ZERO, ret.mat);

	return ret;

}


ZMatrix ZMatrix::operator*(const gsl_complex s) const {
	ZMatrix ret(*this);
	gsl_matrix_complex_scale(ret.mat, s);
	return ret;
}

ZMatrix ZMatrix::operator*(const double s) const {
	gsl_complex s2;
	GSL_SET_COMPLEX(&s2,s,0.0);
	return (*this)*s2;
}

ZMatrix& ZMatrix::operator=(const ZMatrix& m) {
	if(mat == NULL) {
		mat = gsl_matrix_complex_alloc(m.getNumRows(), m.getNumColumns());
	} else if(m.getNumRows() != getNumRows() || m.getNumColumns() != getNumColumns()) {
		gsl_matrix_complex_free(mat);
		mat = gsl_matrix_complex_alloc(m.getNumRows(), m.getNumColumns());
	}
	gsl_matrix_complex_memcpy(mat, m.mat);
	return *this;
}

ZMatrix& ZMatrix::operator+=(const ZMatrix& m) {
	assert(getNumRows() == m.getNumRows() && getNumColumns() == m.getNumColumns());

	gsl_matrix_complex_add(mat,m.mat);

	return *this;
}

ZMatrix& ZMatrix::operator-=(const ZMatrix& m) {
	assert(getNumRows() == m.getNumRows() && getNumColumns() == m.getNumColumns());
	gsl_matrix_complex_sub(mat,m.mat);

	return *this;
}

ZMatrix& ZMatrix::operator*=(const gsl_complex s) {
	gsl_matrix_complex_scale(mat,s);
	return *this;
}

ZMatrix& ZMatrix::operator*=(const double s) {
	gsl_complex s2;
	GSL_SET_COMPLEX(&s2,s,0.0);
	return (*this *= s2);
}

std::ostream& operator<<(std::ostream &os, const ZMatrix &m) {
	std::ios_base::fmtflags oldflag = os.setf(os.showpos);
	for(int i=0; i<m.getNumRows(); ++i) {
		gsl_complex c = gsl_matrix_complex_get(m.mat, i, 0);
		os << GSL_REAL(c) << GSL_IMAG(c) << "i";
		for(int j=1; j<m.getNumColumns(); ++j) {
			c = gsl_matrix_complex_get(m.mat, i, j);
			os << " " << GSL_REAL(c) << GSL_IMAG(c) << "i";
		}
		os << "\n";
	}
	os.setf(oldflag);
	return os;
}

ZMatrix operator*(const gsl_complex s, const ZMatrix &m) {
	return m * s;
}

ZMatrix operator*(const double s, const ZMatrix &m) {
	return m * s;
}

ZMatrix operator*(const gsl_complex s, const Matrix &m) {
	return m * s;
}

ZMatrix operator*(const Matrix &m, const gsl_complex s) {
	ZMatrix ret(m.getNumRows(), m.getNumColumns());
	ret.setReal(m*GSL_REAL(s));
	ret.setImag(m*GSL_IMAG(s));
	return ret;
}
