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


UnitCell::UnitCell()
	: avec(3u,3u), bvec(3u,3u), elasticConsts(NULL), natoms(0),
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
	: avec(3u,3u), bvec(3u,3u)
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
	readNonComment(cellfile, is);
	is >> avec.val(0,0) >> avec.val(0,1) >> avec.val(0,2);
	readNonComment(cellfile, is);
	is >> avec.val(1,0) >> avec.val(1,1) >> avec.val(1,2);
	readNonComment(cellfile, is);
	is >> avec.val(2,0) >> avec.val(2,1) >> avec.val(2,2);

	// calculate reciprocal lattice vector
	bvec = (avec.inverse()).transpose();
	//calculate volume
	volume = scale*scale*scale*avec.det();
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
		atomPos[i] = new double[3];
		readNonComment(cellfile, is);
		is >> atomPos[i][0] >> atomPos[i][1] >> atomPos[i][2];
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
			atomPos[i] = new double[3];
			atomPos[i][0] = c.atomPos[i][0];
			atomPos[i][1] = c.atomPos[i][1];
			atomPos[i][2] = c.atomPos[i][2];
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
	int row[3];

	double kmin = vecmag2(bvec[0]);
	double vec[3];

	for(row[0] = -KMAX_NN; row[0] <= KMAX_NN; ++row[0]) {
		for(row[1] = -KMAX_NN; row[1] <= KMAX_NN; ++row[1]) {
			for(row[2] = -KMAX_NN; row[2] <= KMAX_NN; ++row[2]) {
				if(row[0] == 0 && row[1] == 0 && row[2] == 0)
					continue;
				vec[0] = bvec.val(0,0) * row[0] +
				         bvec.val(1,0) * row[1] +
				         bvec.val(2,0) * row[2];
				vec[1] = bvec.val(0,1) * row[0] +
				         bvec.val(1,1) * row[1] +
				         bvec.val(2,1) * row[2];
				vec[2] = bvec.val(0,2) * row[0] +
				         bvec.val(1,2) * row[1] +
				         bvec.val(2,2) * row[2];

				double vmag = vecmag2(vec);
				if(vmag < kmin) kmin = vmag;
			}
		}
	}

	return 0.5*getRscale()*sqrt(kmin);
}

double UnitCell::getKmax() const {
	return kmax;
}
