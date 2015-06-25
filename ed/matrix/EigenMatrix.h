#ifndef __EIGEN_MATRIX_HPP__
#define __EIGEN_MATRIX_HPP__

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "EDType.h"

/** Dense complex matrix. */
typedef Eigen::Matrix<ComplexType, Eigen::Dynamic, Eigen::Dynamic, Eigen::AutoAlign|Eigen::RowMajor> ComplexMatrixType;
/** Dense real matrix. */
typedef Eigen::Matrix<RealType, Eigen::Dynamic, Eigen::Dynamic, Eigen::AutoAlign|Eigen::RowMajor> RealMatrixType;
/** Default Matrix Type comes from MatrixElemType. */
// typedef Eigen::Matrix<MatrixElemType, Eigen::Dynamic, Eigen::Dynamic, Eigen::AutoAlign|Eigen::RowMajor> MatrixType;

/** Dense complex vector. */
typedef Eigen::Matrix<ComplexType, Eigen::Dynamic, 1, Eigen::AutoAlign> ComplexVectorType;
/** Dense real vector. */
typedef Eigen::Matrix<RealType, Eigen::Dynamic, 1, Eigen::AutoAlign> RealVectorType;
/** Dense vector of integers. */
typedef Eigen::Matrix<int, Eigen::Dynamic, 1, Eigen::AutoAlign> IntVectorType;
/** Default vector type comes from MelemType. */
// typedef Eigen::Matrix<MatrixElemType, Eigen::Dynamic, 1, Eigen::AutoAlign> VectorType;

/** Solvers */
// typedef Eigen::Map<MatrixType> MapMatrix;
// typedef Eigen::SelfAdjointEigenSolver<MatrixType> Dia;
// typedef Eigen::JacobiSVD<MatrixType> SVD;

/** Sparse complex matrix. */
typedef Eigen::SparseMatrix<ComplexType, Eigen::AutoAlign|Eigen::RowMajor> ComplexSparseMatrixType;
/** Sparse real matrix. */
typedef Eigen::SparseMatrix<RealType, Eigen::AutoAlign|Eigen::RowMajor> RealSparseMatrixType;
/** Default Matrix Type comes from MatrixElemType. */
// typedef Eigen::SparseMatrix<MatrixElemType, Eigen::AutoAlign|Eigen::RowMajor> SparseMatrixType;

/** Use to fill sparse matrix*/
typedef Eigen::Triplet<RealType> RealTriplet;
typedef Eigen::Triplet<ComplexType> ComplexTriplet;
// typedef Eigen::Triplet<MatrixElemType> MatrixElemTriplet;
#endif
