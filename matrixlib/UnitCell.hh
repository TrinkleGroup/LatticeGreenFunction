#ifndef _UNIT_CELL_H_
#define _UNIT_CELL_H_

#include <fstream>
#include "Matrix.hh"
#include "SystemDimension.hh"

#define KMAX_NN 3


class UnitCell {
private:
	double scale;
	double rscale; // 2*pi/scale
	double volume; //UnitCell volume: scale^3 det(avec)
	double kmax;   //radius of sphere inscribed in the BZ
	Matrix avec;
	Matrix bvec;
	int c_class;
	double *elasticConsts;
	int natoms;
	double **atomPos;
	double *atomicMass;
	double calcKmax();
public:
	UnitCell();
	UnitCell(UnitCell &);
	UnitCell(std::ifstream &cellfile);
	~UnitCell();
	virtual void load_file(std::ifstream &cellfile);

	virtual Matrix getAvec() const;
	virtual Matrix getBvec() const;
	virtual double getScale() const;
	virtual double getRscale() const;
	virtual double getVolume() const;
	virtual double getKmax() const;

	virtual void getAtomDiff(int, int, double *) const;
};

#endif
