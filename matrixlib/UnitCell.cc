/* File: UnitCell.cc
 * Loads the unit cell for the crystal
 * and stores it
 */
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <iostream>
#include "UnitCell.hh"
#include "VecMath.hh"
#include "SystemDimension.hh"


UnitCell::UnitCell()
	: avec(CARTDIM,CARTDIM), bvec(CARTDIM,CARTDIM), elasticConsts(NULL), natoms(0),
	  atomPos(NULL), atomicMass(NULL), volume(0.0), kmax(0.0)
{
}

static void readNonComment(std::ifstream &file, std::istringstream &is) {
	static int i = 0;
	std::string buf;
	do {
		if(file.eof()) {
			break;
		}
		getline(file, buf);
	} while(buf[0] == '#' && !file.eof());

	/* clear the old line */
	is.clear();
	is.str(buf);
}

UnitCell::UnitCell(std::ifstream &cellfile)
	: avec(CARTDIM,CARTDIM), bvec(CARTDIM,CARTDIM)
{
	load_file(cellfile);
}

void UnitCell::load_file(std::ifstream &cellfile)
{
	std::istringstream is;

	// read lattice scale factor
	readNonComment(cellfile, is);
	is >> scale;

	rscale = 2.0*M_PI/scale;

	// read unit cell in Cartesian coordinates
	for(int idx = 0; idx < CARTDIM; ++idx) {
		readNonComment(cellfile, is);
		for(int jdx = 0; jdx < CARTDIM; ++jdx) {
			is >> avec.val(idx,jdx);
		}
	}

	// calculate reciprocal lattice vector
	bvec = (avec.inverse()).transpose();
	//calculate volume
	
	volume = pow(scale, CARTDIM)*avec.det();
	//calculate kmax
	kmax = calcKmax();

	// read crystal class
	readNonComment(cellfile, is);
	is >> c_class;
	// ignore elastic constants for now, implement when needed
	elasticConsts = NULL;

	readNonComment(cellfile, is);
	is >> natoms;

	atomPos = new double*[natoms];
	atomicMass = new double[natoms];
	for(int i=0; i<natoms; ++i) {
		atomPos[i] = new double[CARTDIM];
		readNonComment(cellfile, is);
		double attmp[CARTDIM];
		for(int j = 0; j < CARTDIM; ++j) {
			is >> attmp[j];
		}
		// Convert from lattice to Cartesian
		for(int j = 0; j < CARTDIM; ++j) {
			atomPos[i][j] = 0.0;
			for(int k = 0; k < CARTDIM; ++k) {
				atomPos[i][j] += avec.val(k,j) * attmp[k];
			}
			atomPos[i][j] *= scale;
		}
	}
	for(int i=0; i<natoms; ++i) {
		readNonComment(cellfile, is);
		is >> atomicMass[i];
	}
}

UnitCell::UnitCell(UnitCell &c)
	: avec(c.avec), bvec(c.bvec), natoms(c.natoms), scale(c.scale),
	  rscale(c.rscale), c_class(c.c_class), volume(c.volume), kmax(c.kmax)
{
	if(natoms == 0) {
		atomPos = NULL;
		atomicMass = NULL;
		elasticConsts = NULL;
	} else {
		atomPos = new double*[natoms];
		atomicMass = new double[natoms];
		for(int i=0; i<natoms; ++i) {
			atomPos[i] = new double[CARTDIM];
			for(int j=0; j<CARTDIM; ++j) {
				atomPos[i][j] = c.atomPos[i][j];
			}
		}
		for(int i=0; i<natoms; ++i) {
			atomicMass[i] = c.atomicMass[i];
		}
		elasticConsts = NULL;
	}
}

UnitCell::~UnitCell() {
	if(atomicMass != NULL)
		delete [] atomicMass;

	for(int i=0; i<natoms; ++i) {
		delete [] atomPos[i];
	}

	if(atomPos != NULL)
		delete [] atomPos;
	
	if(elasticConsts != NULL)
		delete [] elasticConsts;
}

Matrix UnitCell::getAvec() const {
	return avec;
}

Matrix UnitCell::getBvec() const {
	return bvec;
}

double UnitCell::getScale() const {
	return scale;
}

double UnitCell::getRscale() const {
	return rscale;
}

double UnitCell::getVolume() const {
	return volume;
}

double UnitCell::calcKmax() {
	int row[CARTDIM];

	double kmin = vecmag2(bvec[0]);
	double vec[CARTDIM];

	/* initialize the nn vector of integers */
	for(int idx = 0; idx < CARTDIM; ++idx) {
		row[idx] = -KMAX_NN;
	}

	while(row[0] <= KMAX_NN) {
		for(int idx = 0; idx < CARTDIM; ++idx) {
			vec[idx] = 0.0;
			for(int jdx = 0; jdx < CARTDIM; ++jdx) {
				vec[idx] += bvec.val(jdx,idx) * row[jdx];
			}
		}
		double vmag = vecmag2(vec);
		if(vmag < kmin) kmin = vmag;
		int counter = CARTDIM - 1;

		while(++row[counter] > KMAX_NN && counter > 0) {
			row[counter] = -KMAX_NN;
			counter--;
		}
		bool iszero = true;
		for(int idx = 0; idx < CARTDIM; ++idx) {
			if(row[idx] != 0) {
				iszero = false;
				break;
			}
		}
		counter = CARTDIM - 1;
		if(iszero) {

			while(++row[counter] > KMAX_NN && counter > 0) {
				row[counter] = -KMAX_NN;
				counter--;
			}
		}

	}

	return 0.5*getRscale()*sqrt(kmin);
}

double UnitCell::getKmax() const {
	return kmax;
}

void UnitCell::getAtomDiff(int i, int j, double *rr) const {
	for(int idx=0; idx<CARTDIM; ++idx) {
		rr[idx] = atomPos[i][idx] - atomPos[j][idx];
	}
}
