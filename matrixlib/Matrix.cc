/* File: Matrix.cc
 * wrapper for gsl_matrix
 */
#include <cassert>
#include <cstring>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_sort_vector.h>
#include "Matrix.hh"

Matrix::Matrix() :  mat(NULL) { }

/* create zero m x n matrix */
Matrix::Matrix(unsigned int m, unsigned int n) : mat(gsl_matrix_calloc(m,n)) { }

/* create m x n matrix from the double array matrix */
Matrix::Matrix(double **matrix, unsigned int m, unsigned int n) :
	mat(gsl_matrix_alloc(m,n))
{
	for(unsigned int i=0u; i<m; ++i)
		for(unsigned int j=0u; j<n; ++j)
			gsl_matrix_set(mat,i,j, matrix[i][j]);
}

/* create m x n matrix from the signle array matrix */
Matrix::Matrix(double *matrix, unsigned int m, unsigned int n) :
	mat(gsl_matrix_alloc(m,n))
{
	for(unsigned int i=0u; i<m; ++i)
		for(unsigned int j=0u; j<n; ++j)
			gsl_matrix_set(mat,i,j,matrix[i*n + j]);
}

/* copy constructor */
Matrix::Matrix(const Matrix &matrix) {
	this->mat = gsl_matrix_alloc(matrix.getNumRows(), matrix.getNumColumns());
	gsl_matrix_memcpy(this->mat, matrix.mat);
}

/* clean up */
Matrix::~Matrix() {
	if(mat != NULL) {
		gsl_matrix_free(mat);
		mat = NULL;
	}
}

/* get a value at i,k in the matrix, reference to the value */
double &Matrix::val(unsigned int i, unsigned int j) {
	return *gsl_matrix_ptr(mat, i, j);
}

/* get a value at i,k in the matrix */
double Matrix::val(unsigned int i, unsigned int j) const {
	return gsl_matrix_get(mat, i, j);
}


//matrix inversion
void Matrix::invert() {
	//must be a square matrix
	assert(getNumRows() == getNumColumns());
	int s;
	gsl_permutation *p = gsl_permutation_alloc(getNumRows());
	gsl_matrix *lu = gsl_matrix_alloc(getNumRows(), getNumColumns());
	gsl_matrix_memcpy(lu,mat);
	gsl_linalg_LU_decomp(lu, p, &s);
	gsl_linalg_LU_invert(lu, p, mat);
	gsl_matrix_free(lu);
	gsl_permutation_free(p);
}

//matrix determinant
double Matrix::det() const {
	//must be a square matrix
	assert(getNumRows() == getNumColumns());
	int s;
	gsl_permutation *p = gsl_permutation_alloc(getNumRows());
	gsl_matrix *lu = gsl_matrix_alloc(getNumRows(), getNumColumns());
	gsl_matrix_memcpy(lu,mat);
	gsl_linalg_LU_decomp(lu, p, &s);
	double determ = gsl_linalg_LU_det(lu, s);
	gsl_matrix_free(lu);
	gsl_permutation_free(p);
	return determ;
}


/* helper function to build a complex matrix from real and imaginary parts */
static gsl_matrix_complex * build_z_matrix(const Matrix &real, const Matrix imag) {
	assert(real.getNumRows() == imag.getNumRows()
	    && real.getNumColumns() == imag.getNumColumns());
	gsl_matrix_complex *lu = gsl_matrix_complex_alloc(real.getNumRows(),
	                                                  real.getNumColumns());
	for(int i=0; i<lu->size1; ++i) {
		for(int j=0; j<lu->size2; ++j) {
			gsl_complex q;
			GSL_SET_COMPLEX(&q,real.val(i,j), imag.val(i,j));
			gsl_matrix_complex_set(lu, i, j, q);
		}
	}

	return lu;

}

//complex inversion from real and imaginary parts */
void invertMatrixZ(const Matrix &real, const Matrix &imag, Matrix &realOut,
                   Matrix &imagOut)
{
	assert(real.getNumRows() == imag.getNumRows()
	    && real.getNumColumns() == imag.getNumColumns()
	    && real.getNumRows() == realOut.getNumRows()
	    && real.getNumColumns() == realOut.getNumColumns()
	    && real.getNumRows() == imagOut.getNumRows()
	    && real.getNumColumns() == imagOut.getNumColumns());

	int s;
	gsl_permutation *p = gsl_permutation_alloc(real.getNumRows());
	gsl_matrix_complex *lu = build_z_matrix(real,imag);

	int err;
	if((err = gsl_linalg_complex_LU_decomp(lu, p, &s))!=0) {
		std::cerr << "LU decomp problem! " << "err=" << err <<std::endl;
	}

	gsl_matrix_complex *inv = gsl_matrix_complex_alloc(real.getNumRows(),
	                                                  real.getNumColumns());

	if((err=gsl_linalg_complex_LU_invert(lu, p, inv))!=0) {
		std::cerr << "Inversion problem! " << "err=" << err <<std::endl;
	}

	gsl_matrix_complex_free(lu);
	gsl_permutation_free(p);
	for(int i=0; i<inv->size1; ++i) {
		for(int j=0; j<inv->size2; ++j) {
			realOut[i][j] = GSL_REAL(gsl_matrix_complex_get(inv, i, j));
			imagOut[i][j] = GSL_IMAG(gsl_matrix_complex_get(inv, i, j));
		}
	}
	gsl_matrix_complex_free(inv);
}

//matrix transposition
void Matrix::transp() {
	//square matrices can do an inplace transpose
	if(getNumRows() == getNumColumns()) {
		gsl_matrix_transpose(mat);
	} else {
		gsl_matrix *m = gsl_matrix_alloc(getNumColumns(), getNumRows());
		gsl_matrix_transpose_memcpy(m, mat);
		gsl_matrix_free(mat);
		mat = m;
	}
}

/* hermitian matrix eigenvalues: */
void herm_evals(const Matrix &real, const Matrix &imag, double *evals) {
	assert(real.getNumRows() == imag.getNumRows()
	    && real.getNumColumns() == imag.getNumColumns()
	    && real.getNumRows() == real.getNumColumns());
	
	gsl_vector *ev = gsl_vector_alloc(real.getNumRows());

	gsl_matrix_complex *zmat = build_z_matrix(real,imag);
	gsl_eigen_herm_workspace *w = gsl_eigen_herm_alloc(real.getNumRows());

	gsl_eigen_herm(zmat, ev, w);

	gsl_sort_vector(ev);

	for(unsigned int i = 0u; i < real.getNumRows(); ++i) {
		evals[i] = gsl_vector_get(ev,i);
	}

	gsl_eigen_herm_free(w);
	gsl_matrix_complex_free(zmat);
	gsl_vector_free(ev);
}

/* Real symmetric matric eigenvectors and eigenvalues */
void Matrix::symm_eigen(double **&evecs, double *&evals) const {
	assert(getNumRows() == getNumColumns());
	//gsl_eigen_symmv destroys the lower triangle and diagonal,
	//so we'll make a copy
	gsl_matrix *tmpmat = gsl_matrix_alloc(getNumRows(), getNumColumns());
	gsl_matrix_memcpy(tmpmat, mat);

	gsl_vector *eval = gsl_vector_alloc(getNumRows());
	gsl_matrix *evec = gsl_matrix_alloc(getNumRows(),getNumColumns());


	gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc(getNumRows());
	gsl_eigen_symmv(tmpmat, eval, evec, w);

	gsl_eigen_symmv_free(w);

	gsl_matrix_free(tmpmat);


	//sort eigenvalues in ascending order of magnitude
	gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);

	evecs = new double*[getNumColumns()];
	evals = new double[getNumRows()];
	for(unsigned int i = 0u; i < getNumColumns(); ++i) {
		evals[i] = gsl_vector_get(eval,i);
		gsl_vector_view gv = gsl_matrix_column(evec, i);
		evecs[i] = new double[getNumRows()];
		for(unsigned int j=0u; j<getNumRows(); ++j) {
			evecs[i][j] = gsl_vector_get(&gv.vector,j);
		}
	}
	gsl_vector_free(eval);
	gsl_matrix_free(evec);
}

/* Real symmetric matrix eigenvalues */
void Matrix::symm_evals(double *evals) const {
	assert(getNumRows() == getNumColumns());

	gsl_vector *ev = gsl_vector_alloc(getNumRows());

	gsl_matrix *zmat = gsl_matrix_alloc(getNumRows(),getNumColumns());
	gsl_matrix_memcpy(zmat, mat);
	gsl_eigen_symm_workspace *w = gsl_eigen_symm_alloc(getNumRows());

	gsl_eigen_symm(zmat, ev, w);

	gsl_sort_vector(ev);

	for(unsigned int i = 0u; i < getNumRows(); ++i) {
		evals[i] = gsl_vector_get(ev,i);
	}

	gsl_eigen_symm_free(w);
	gsl_matrix_free(zmat);
	gsl_vector_free(ev);
}

/* is the matrix a square matrix? */
bool Matrix::isSquare() const {
	return getNumRows() == getNumColumns();
}

/* matrix inversion */
Matrix Matrix::inverse() const {
	Matrix ret(*this);
	ret.invert();
	return ret;
}

/* matrix transpose */
Matrix Matrix::transpose() const {
	Matrix ret(*this);
	ret.transp();
	return ret;
}

/* real matrix, adjoint == transpose */
Matrix Matrix::adjoint() const {
	Matrix ret(*this);
	ret.transp();
	return ret;
}

/* return the n x n identity matrix */
Matrix Matrix::eye(unsigned int n) {
	Matrix ret(n,n);
	gsl_matrix_set_identity(ret.mat);
	return ret;
}

unsigned int Matrix::getNumRows() const {
	return mat->size1;
}

unsigned int Matrix::getNumColumns() const {
	return mat->size2;
}

void Matrix::setZero() {
	gsl_matrix_set_zero(mat);
}

/* index of matrix row */
double *Matrix::operator[](int i) {
	return &(mat->data[i*mat->tda]);
}

/* add two matrices */
Matrix Matrix::operator+(const Matrix& m) const {
	assert(getNumRows() == m.getNumRows() && getNumColumns() == m.getNumColumns());

	Matrix ret(*this);

	gsl_matrix_add(ret.mat, m.mat);

	return ret;
}

/* subtract two matrices */
Matrix Matrix::operator-(const Matrix& m) const {
	assert(getNumRows() == m.getNumRows() && getNumColumns() == m.getNumColumns());

	Matrix ret(*this);

	gsl_matrix_sub(ret.mat, m.mat);

	return ret;
}

/* negate a matrix */
Matrix Matrix::operator-() const {
	Matrix ret(*this);
	gsl_matrix_scale(ret.mat, -1.0);
	return ret;
}

/* matrix multiplication */
Matrix Matrix::operator*(const Matrix& m) const {
	assert(getNumColumns() == m.getNumRows());

	Matrix ret(getNumRows(),m.getNumColumns());

	/* ret = this * m */

	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0,
	               this->mat, m.mat, 0.0, ret.mat);

	return ret;

}


/* multiple matrix by a scalar */
Matrix Matrix::operator*(const double s) const {
	Matrix ret(*this);
	gsl_matrix_scale(ret.mat, s);
	return ret;
}

/* matrix assignment */
Matrix& Matrix::operator=(const Matrix& m) {
	if(mat == NULL) {
		mat = gsl_matrix_alloc(m.getNumRows(), m.getNumColumns());
	} else if(m.getNumRows() != getNumRows() || m.getNumColumns() != getNumColumns()) {
		gsl_matrix_free(mat);
		mat = gsl_matrix_alloc(m.getNumRows(), m.getNumColumns());
	}
	gsl_matrix_memcpy(mat, m.mat);
	return *this;
}

Matrix& Matrix::operator+=(const Matrix& m) {
	assert(getNumRows() == m.getNumRows() && getNumColumns() == m.getNumColumns());

	gsl_matrix_add(mat,m.mat);

	return *this;
}

Matrix& Matrix::operator-=(const Matrix& m) {
	assert(getNumRows() == m.getNumRows() && getNumColumns() == m.getNumColumns());
	gsl_matrix_sub(mat,m.mat);

	return *this;
}

Matrix& Matrix::operator*=(const double s) {
	gsl_matrix_scale(mat,s);
	return *this;
}

/* output a formatted matrix to an output stream *//
std::ostream& operator<<(std::ostream &os, const Matrix &m) {
	if(m.getNumRows() > 0 && m.getNumColumns() > 0) {
		for(int i=0; i<m.getNumRows(); ++i) {
			os << gsl_matrix_get(m.mat, i, 0);
			for(int j=1; j<m.getNumColumns(); ++j) {
				os << " " << gsl_matrix_get(m.mat, i, j);
			}
			os << "\n";
		}
	}
	return os;
}

Matrix operator*(const double s, const Matrix &m) {
	return m * s;
}
