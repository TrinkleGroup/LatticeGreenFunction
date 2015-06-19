#ifndef _ALLMATS_H_
#define _ALLMATS_H_

#include "Matrix.hh"
#include "DynMat.hh"

Matrix one_on_k_2_total(DynMat &dynmat, double curr_k[3]);
Matrix i_on_k_total(DynMat &dynmat, double curr_k[3]);
Matrix k0_total(DynMat &dynmat, double curr_k[3]);

Matrix one_on_k_2_total_rot(DynMat &dynmat, double curr_k[3]);
Matrix i_on_k_total_rot(DynMat &dynmat, double curr_k[3]);
Matrix k0_total_rot(DynMat &dynmat, double curr_k[3]);

Matrix XiInvAA(DynMat &dynmat, double k[3]);
Matrix XiInvAO(DynMat &dynmat, double k[3]);
Matrix XiInvOO(DynMat &dynmat, double k[3]);

Matrix i_on_k_sub(DynMat &dynmat, double curr_k[3]);
Matrix k0_sub(DynMat &dynmat, double curr_k[3]);

#endif
