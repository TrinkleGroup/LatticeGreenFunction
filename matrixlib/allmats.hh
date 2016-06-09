#ifndef _ALLMATS_H_
#define _ALLMATS_H_

#include "SystemDimension.hh"
#include "Matrix.hh"
#include "DynMat.hh"

Matrix one_on_k_2_total(DynMat &dynmat, double *curr_k);
Matrix i_on_k_total(DynMat &dynmat, double *curr_k);
Matrix k0_total(DynMat &dynmat, double *curr_k);

Matrix one_on_k_2_total_rot(DynMat &dynmat, double *curr_k);
Matrix i_on_k_total_rot(DynMat &dynmat, double *curr_k);
Matrix k0_total_rot(DynMat &dynmat, double *curr_k);

Matrix XiInvAA(DynMat &dynmat, double *k);
Matrix XiInvAO(DynMat &dynmat, double *k);
Matrix XiInvOO(DynMat &dynmat, double *k);

Matrix i_on_k_sub(DynMat &dynmat, double *curr_k);
Matrix k0_sub(DynMat &dynmat, double *curr_k);

Matrix one_on_k_2_1a(DynMat &dynmat, double *curr_k);
Matrix k0_1a(DynMat &dynmat, double *curr_k);
Matrix XiInv_1a(DynMat &dynmat, double *k);
Matrix k0_sub_1a(DynMat &dynmat, double *curr_k);
#endif
