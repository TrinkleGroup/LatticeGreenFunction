#ifndef _Z_MATRIX_H_
#define _Z_MATRIX_H_

#include <iostream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>

#include "Matrix.hh"

class ZMatrix {
private:
	gsl_matrix_complex *mat;
public:
	ZMatrix();
	ZMatrix(unsigned int m, unsigned int n);
	ZMatrix(gsl_complex **matrix, unsigned int m, unsigned int n);
	ZMatrix(gsl_complex *matrix, unsigned int m, unsigned int n);
	ZMatrix(const ZMatrix &matrix);
	ZMatrix(const Matrix &real, const Matrix &imag);
	ZMatrix(const Matrix &real); //implied zero imaginary part
	~ZMatrix();

	virtual gsl_complex &val(unsigned int i, unsigned int j);
	virtual gsl_complex val(unsigned int i, unsigned int j) const;

	virtual void herm_eigen(gsl_complex **&evecs, double *&evals) const;
	virtual void herm_evals(double *evals) const;

	virtual double cond1est() const;

	virtual bool isSquare() const;

	// These act in place, modifying the current ZMatrix
	void invert();
	void transp();
	void adj();
	void adjInPlace();

	virtual gsl_complex det() const;
	//These don't modify the matrix
	virtual ZMatrix inverse() const;
	virtual ZMatrix transpose() const;
	virtual ZMatrix adjoint() const;

	static ZMatrix eye(unsigned int n);

	virtual void setZero();
	virtual void setReal(const Matrix &m);
	virtual void setImag(const Matrix &m);
	virtual Matrix getReal() const;
	virtual Matrix getImag() const;

	virtual unsigned int getNumRows() const;
	virtual unsigned int getNumColumns() const;

	/*
	virtual gsl_complex *operator[](int i);
	*/
	virtual ZMatrix operator+(const ZMatrix& m) const;
	virtual ZMatrix operator-(const ZMatrix& m) const;
	virtual ZMatrix operator-() const;
	virtual ZMatrix operator*(const ZMatrix& m) const;
	virtual ZMatrix operator*(const gsl_complex s) const;
	virtual ZMatrix operator*(const double s) const;
	virtual ZMatrix& operator=(const ZMatrix& m);
	virtual ZMatrix& operator+=(const ZMatrix& m);
	virtual ZMatrix& operator-=(const ZMatrix& m);
	virtual ZMatrix& operator*=(const gsl_complex s);
	virtual ZMatrix& operator*=(const double s);

	friend std::ostream& operator<<(std::ostream &os, const ZMatrix &m);
};

ZMatrix operator*(const double s, const ZMatrix& m);
ZMatrix operator*(const gsl_complex s, const ZMatrix& m);
ZMatrix operator*(const gsl_complex s, const Matrix& m);
ZMatrix operator*(const Matrix& m, const gsl_complex s);
#endif
