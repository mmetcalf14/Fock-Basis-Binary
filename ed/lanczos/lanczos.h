#ifndef __LANCZOS_H__
#define __LANCZOS_H__
#include "EDType.h"
#include "matrix/EigenMatrix.h"

bool LanczosEV(const size_t N, const RealSparseMatrixType A,
  RealVectorType &Vec, RealType &Val,
  size_t &max_iter, double err_tol = 1.0E-7);

bool LanczosEV(const size_t N, const ComplexSparseMatrixType A,
  ComplexVectorType &Vec, RealType &Val,
  size_t &max_iter, double err_tol = 1.0E-7);

#endif//__LANCZOS_H__