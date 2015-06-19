#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <iostream>
#include <cmath>
#include <gsl/gsl_matrix.h>

class Matrix {
private:
	gsl_matrix *mat;

public:
	Matrix();
	Matrix(unsigned int m, unsigned int n);
	Matrix(double **matrix, unsigned int m, unsigned int n);
	Matrix(double *matrix, unsigned int m, unsigned int n);
	Matrix(const Matrix &matrix);
	~Matrix();

	virtual double &val(unsigned int i, unsigned int j);
	virtual double val(unsigned int i, unsigned int j) const;

	virtual void symm_eigen(double **&evecs, double *&evals) const;
	virtual void symm_evals(double *evals) const;

	virtual bool isSquare() const;

	// These modify the current matrix by acting inline
	void invert();
	void transp();

	virtual double det() const;
	// These leave the current matrix untouched
	virtual Matrix inverse() const;
	virtual Matrix transpose() const;
	virtual Matrix adjoint() const;

	static Matrix eye(unsigned int n);

	virtual void setZero();

	virtual unsigned int getNumRows() const;
	virtual unsigned int getNumColumns() const;

	virtual double *operator[](int i);
	virtual Matrix operator+(const Matrix& m) const;
	virtual Matrix operator-(const Matrix& m) const;
	virtual Matrix operator-() const;
	virtual Matrix operator*(const Matrix& m) const;
	virtual Matrix operator*(const double s) const;
	virtual Matrix& operator=(const Matrix& m);
	virtual Matrix& operator+=(const Matrix& m);
	virtual Matrix& operator-=(const Matrix& m);
	virtual Matrix& operator*=(const double s);

	friend std::ostream& operator<<(std::ostream &os, const Matrix &m);
	friend void invertMatrixZ(const Matrix &real, const Matrix &imag, 
	                          Matrix &realOut, Matrix &imagOut);
	friend void herm_evals(const Matrix &real, const Matrix &imag, 
	                       double *evals);
};

Matrix operator*(const double s, const Matrix& m);

#endif
