#ifndef __HAMILTONIAN_H__
#define __HAMILTONIAN_H__
#include <vector>
#include "EDType.h"
#include "fermion/basis.h"
#include "matrix/EigenMatrix.h"

template<typename Tnum>
class FermiHubbard: public FermionBasis{
public:
  typedef Eigen::Matrix<Tnum, Eigen::Dynamic, 1, Eigen::AutoAlign> VectorType;
  typedef Eigen::Matrix<Tnum, Eigen::Dynamic, Eigen::Dynamic, Eigen::AutoAlign|Eigen::RowMajor> MatrixType;
  typedef Eigen::SparseMatrix<Tnum, Eigen::AutoAlign|Eigen::RowMajor> SparseMatrixType;
  typedef Eigen::Triplet<Tnum> MatrixElemT;
  typedef Eigen::Map<MatrixType> MapMatrix;
  typedef Eigen::SelfAdjointEigenSolver<MatrixType> Dia;
  typedef Eigen::JacobiSVD<MatrixType> SVD;
  FermiHubbard (size_t L, size_t Nup, size_t Ndn, Tnum t_up = 1.0, Tnum t_dn = 1.0);
  virtual ~FermiHubbard ();
  void Build1DHoppingTerms();
  void BuildOneBodyTerms();
  void BuildTwoBodyTerms( const Tnum U );
  void ConstructTotalHamiltonian();
  void eigh()const;

private:
  Tnum t = 1.0;
  uint64_t dim;
  std::vector<Tnum> tList;
  std::vector<MatrixElemT> tripletList;
  SparseMatrixType H0;
  SparseMatrixType HOne;
  SparseMatrixType HTwo;
  SparseMatrixType Htot;
};

#endif// __HAMILTONIAN_H__
