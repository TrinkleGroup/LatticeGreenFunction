#ifndef _UNIT_CELL_H_
#define _UNIT_CELL_H_

#include <fstream>
#include "Matrix.hh"

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
	void load_file(std::ifstream &cellfile);

	Matrix getAvec() const;
	Matrix getBvec() const;
	double getScale() const;
	double getRscale() const;
	double getVolume() const;
	double getKmax() const;
};

#endif
