#ifndef _DYN_MAT_H_
#define _DYN_MAT_H_

#include <iostream>
#include <string>
#include <vector>
#include "Matrix.hh"
#include "ZMatrix.hh"
#include "SystemDimension.hh"
#include "UnitCell.hh"

class DynMat {
protected:
	unsigned int nions; //Number of ions in the unit cell
	unsigned int numr; //Number of R vectors in the Dynamical matrix
	std::string description; // simple description of the dynamical matrix
	double *R; // CARTDIM*numr (CARTDIM*r_i + j)
	Matrix **mat; // CARTDIMnions x CARTDIMnions x numr

	UnitCell *uc;

	std::vector<Matrix *> orderCache;  //k point cache of the k^n terms of the dynamical matrix Fourier space Taylor expansion
	double *cached_k; //the currently cached k vector

	Matrix *ap;  //projector onto the acoustic modes CARTDIMxCARTDIM
	Matrix *op;  //projector onto the optical modes CARTDIM(N-1)xCARTDIM(N-1)
	Matrix *apt; //transpose of ap
	Matrix *opt; //transpose of op

	Matrix *apf;  //projector onto the acoustic modes full CARTDIMNxCARTDIMN
	Matrix *opf;  //projector onto the optical modes full CARTDIMNxCARTDIMN
	Matrix *aptf; //transpose of apf
	Matrix *optf; //transpose of opf

	Matrix *rot;  //Rotation matrix into the acoustic/optical basis
	Matrix *rotT; //transpose of rot

	virtual void clean_up(); //clean up memory
	virtual void calcProjectors(); //calculate the projection matrices
	virtual void calcAcousticRots(); //calculate the rotation matrices

public:
	DynMat();

	/* Construct from standard C++ double arrays
	 * matrix: The R x CARTDIMN x CARTDIMN dynamical matrix,
	 * pointer indexes are: R, i, j
	 * where R labels the lattice vector in R
	 * i is a cartesian index and atom index x1,y1,z1,x2,y2,x2,...
	 * j is a cartesian index and atom index x1,y1,z1,x2,y2,x2,...
	 * nions: number of ions N
	 * numr: number of lattice vectors R
	 */
	DynMat(double ***matrix, double *R, unsigned int nions, unsigned int numr, UnitCell *uc);

	/* construct from Matrix wrapper class
	 * matrix: The R x CARTDIMN x CARTDIMN dynamical matrix,
	 * pointer index is for: R vectors
	 * Each Matrix is CARTDIMN x CARTDIMN with indexes:
	 * i is a cartesian index and atom index x1,y1,z1,x2,y2,x2,...
	 * j is a cartesian index and atom index x1,y1,z1,x2,y2,x2,...
	 * nions: number of ions N
	 * numr: number of lattice vectors R
	 */
	DynMat(Matrix *matrix, double *R, unsigned int nions, unsigned int numr, UnitCell *uc);
	~DynMat();

	/* load Dynamical matrix from an input stream */
	virtual void load(std::istream &is, UnitCell &uc);

	
	/* Calculate the k^n Taylor expansion term in the
	 * Fourier space Dynamical matrix
	 */
	virtual const Matrix &getOrder(int n, double *k);

	virtual const Matrix &getAcousticProjector() const;
	virtual const Matrix &getOpticalProjector() const;
	virtual const Matrix &getAcousticProjectorT() const;
	virtual const Matrix &getOpticalProjectorT() const;

	virtual const Matrix &getAcousticFullP() const;
	virtual const Matrix &getAcousticFullPT() const;
	virtual const Matrix &getOpticalFullP() const;
	virtual const Matrix &getOpticalFullPT() const;

	virtual const Matrix &getRot() const;
	virtual const Matrix &getRotT() const;

	/* calculate the full real and imaginary pieces of the Fourier transform
	 * of the dynamical matrix at k
	 */
	virtual void getFourierTransform(double *k, Matrix &realFt, Matrix &imFt) const;

	/* calculate the full real and imaginary pieces of the Fourier transform
	 * of the dynamical matrix at k, complex matrix output
	 */
	virtual void getFourierTransform(double *k, ZMatrix &ft) const;

	/* calculate the full real and imaginary pieces of the Fourier transform
	 * of the dynamical matrix at k in the rotated Acoustic/Optical basis
	 */
	virtual void getFourierTransformRot(double *k, Matrix &realFt, Matrix &imFt) const;

	/* calculate the full real and imaginary pieces of the Fourier transform
	 * of the dynamical matrix at k in the rotated Acoustic/Optical basis
	 * complex matrix output
	 */
	virtual void getFourierTransformRot(double *k, ZMatrix &ft) const;

	/* calculate the small k expansion real and imaginary pieces of the
	 * Fourier transform of the dynamical matrix at k in the
	 * rotated Acoustic/Optical basis. complex matrix output.
	 */
	virtual void getSmallFourierTransform(double *k, ZMatrix &ft) const;
	virtual void getSmallFourierTransform1a(double *k, ZMatrix &ft) const;
	virtual void getSmallFourierTransformMa(double *k, ZMatrix &ft) const;

	virtual const UnitCell *getUC() const;
	virtual unsigned int getNions() const;
	virtual unsigned int getNumr() const;
};
#endif
