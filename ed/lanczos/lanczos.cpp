#include <cstring>
#include "lanczos/lanczos.h"
#ifdef MKL
    #include "mkl.h"
#else
    #include "numeric/lapack_wrapper.h"
#endif

#ifndef DEBUG
#define DEBUG 5
#endif


bool LanczosEV(const size_t N, const ComplexSparseMatrixType A,
  ComplexVectorType &Vec, RealType &Val,
  size_t &max_iter, double err_tol)
{
  RUNTIME_ERROR("Lanczos on complex<double> is not supported yet!");
}

bool LanczosEV(const size_t N, const RealSparseMatrixType A,
  RealVectorType &Vec, RealType &Val, size_t &max_iter, double err_tol)
{
  const int min_iter = 2;
  const RealType beta_err = 1E-14;
  if( max_iter < min_iter ){
    RUNTIME_ERROR("Maximum iteration number should be set greater than 2.");
  }
  RealType alpha;
  RealType beta = 1.0;
  size_t M = max_iter;
  std::vector<RealVectorType> Vm;
  RealType *Alphas = (RealType*)malloc( M * sizeof(RealType) );
  RealType *Betas = (RealType*)malloc( M * sizeof(RealType) );
  RealType *d = (RealType*)malloc(M * sizeof(RealType));
  RealType *e = (RealType*)malloc(M * sizeof(RealType));
  memset(Alphas, 0, M * sizeof(RealType));
  memset(Betas, 0, M * sizeof(RealType));
  //NOTE: normalized Vec
  Vec.normalize();
  Vm.push_back(Vec);
  RealType e_diff = 1.0;
  RealType e0_old = 0.0;
  bool converged = false;
  int it = 0;
  while ( ( ((e_diff > err_tol) && it < max_iter) || it < min_iter ) &&
          beta > beta_err ) {
    RealVectorType work = A * Vm[it];
    if( it > 0 ) work -= beta * Vm[it-1];
    Vm.push_back(work);
    alpha = Vm[it+1].dot( Vm[it] );
    Vm[it+1] -= alpha * Vm[it];
    beta = Vm[it+1].norm();
    Alphas[it] = alpha;
    if( DEBUG > 5 ){
      INFO("alpha @ " << it << " is " << alpha);
      INFO(" beta @ " << it << " is " << beta);
    }
    if( beta > beta_err ){
      Vm[it+1].normalize();
      if( it < max_iter - 1) Betas[it] = beta;
    }
    else{
      converged = true;
    }
    it++;
    if( it > 1 ){
      RealType* z = (RealType*)malloc(it * it * sizeof(RealType));
      RealType* work = (RealType*)malloc(4 * it * sizeof(RealType));
      int info;
      memcpy(d, Alphas, it * sizeof(RealType));
      memcpy(e, Betas, it * sizeof(RealType));
      //dstev - LAPACK
      //      DSTEV computes all eigenvalues and, optionally, eigenvectors of a
      //            real symmetric tridiagonal matrix A.
      //      'N':  Compute eigenvalues only;
      //       d is diagonal.
      //       e is sub-diagonal.
      dstev((char*)"N", &it, d, e, z, &it, work, &info);
      if(info != 0){
        INFO("Lapack INFO = " << info);
        RUNTIME_ERROR("Error in Lapack function 'dstev'");
      }
      RealType base = fabs(d[0]) > 1 ? fabs(d[0]) : 1;
      e_diff = fabs(d[0] - e0_old) / base;
      e0_old = d[0];
      if(e_diff <= err_tol)
        converged = true;
      free(z), free(work);
    }
  }
  if( it > 1 ){
    memcpy(d, Alphas, it * sizeof(RealType));
    memcpy(e, Betas, it * sizeof(RealType));
    RealType* z = (RealType*)malloc(it * it * sizeof(RealType));
    RealType* work = (RealType*)malloc(4 * it * sizeof(RealType));
    int info;
    //dstev - LAPACK
    dstev((char*)"V", &it, d, e, z, &it, work, &info);
    if(info != 0){
      INFO("Lapack INFO = " << info);
      RUNTIME_ERROR("Error in Lapack function 'dstev'");
    }
    Vec.setZero();
    for(size_t k = 0; k < it; k++){
      Vec += z[k] * Vm[k];
    }
    max_iter = it;
    Val = d[0];
    free(z), free(work);
  }
  else{
    max_iter = 1;
    Val = 0;
  }
  free(Alphas), free(Betas), free(d), free(e);
  return converged;
}
