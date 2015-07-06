#include "fermion/hamiltonian.h"
#include "fermion/bitwise.h"
// #include "fermion/search.h"
#include "lanczos/lanczos.h"

#ifndef DEBUG
#define DEBUG 3
#endif

template<typename Tnum>
FermiHubbard<Tnum>::FermiHubbard(size_t L, size_t Nup, size_t Ndn,
  Tnum t_up, Tnum t_dn):
FermionBasis( L, Nup, Ndn )
{
  dim = getTotalHilbertSpace();
  // INFO(dim);
  H0.resize(dim, dim);
  H0.reserve(2*dim);
  HOne.resize(dim, dim);
  HOne.reserve(dim);
  HTwo.resize(dim, dim);
  HTwo.reserve(dim);
  Htot.resize(dim, dim);
  Htot.reserve(3*dim);

  tList.push_back(t_up);
  tList.push_back(t_dn);

  Build1DHoppingTerms();
}

template<typename Tnum>
FermiHubbard<Tnum>::~FermiHubbard()
{
  INFO("FermiHubbard Destructor called.");
  tripletList.clear();
  H0.setZero();
  HOne.setZero();
  HTwo.setZero();
  Htot.setZero();
}

template<typename Tnum>
void FermiHubbard<Tnum>::Build1DHoppingTerms()
{
  tripletList.clear();
  size_t rid, cid;
  uint64_t count_up = getEachHilbertSpace(0);
  uint64_t count_dn = getEachHilbertSpace(1);
  size_t lid = 0, pid = 0;//l and p's index
  /* For up-spin */
  for( auto l : Basis[0] ){
    for (size_t i = 0; i < getL() - 1; i++) {
      size_t j = i + 1;// hopping between j - i
      /* c^\dagger_ic_j */
      if ( btest(l, j) && !(btest(l, i)) ) {
        /* if yes, no particle in i and one particle in j. */
        size_t p = ibset(l,i);
        p = ibclr(p,j);
        // auto val = binary_locate(Basis[0].cbegin(), Basis[0].cend(), p);
        // pid = std::distance(Basis[0].cbegin(), val);
        lid = Index[0][l];
        pid = Index[0][p];
        // INFO(lid << " " << pid);
        for (size_t cnt = 0; cnt < count_dn; cnt++) {
          rid = getIndices2(lid, cnt);//cnt * count_up + lid;
          cid = getIndices2(pid, cnt);//cnt * count_up + pid;
          //NOTE: there is a convention here.
          Tnum value = -t;
          tripletList.push_back(MatrixElemT(rid, cid, value));
          tripletList.push_back(MatrixElemT(cid, rid, value));
          if( DEBUG > 4) INFO(lid << " " << pid << " " << rid << " " << cid << " " << value);
        }
      }
    }
  }
  for( auto l : Basis[1] ){
    for (size_t i = 0; i < getL(); i++) {
      size_t j = i + 1;// hopping between j - i
      /* c^\dagger_ic_j */
      if ( btest(l, j) && !(btest(l, i)) ) {
        /* if yes, no particle in i and one particle in j. */
        size_t p = ibset(l,i);
        p = ibclr(p,j);
        // auto val = binary_locate(Basis[1].cbegin(), Basis[1].cend(), p);
        // pid = std::distance(Basis[1].cbegin(), val);
        lid = Index[1][l];
        pid = Index[1][p];
        // INFO(lid << " " << pid);
        for (size_t cnt = 0; cnt < count_up; cnt++) {
          rid = getIndices2(cnt, lid);//lid * count_up + cnt;
          cid = getIndices2(cnt, pid);//pid * count_up + cnt;
          Tnum value = -t;
          tripletList.push_back(MatrixElemT(rid, cid, value));
          tripletList.push_back(MatrixElemT(cid, rid, value));
          if( DEBUG > 4) INFO(lid << " " << pid << " " << rid << " " << cid << " " << value);
        }
      }
    }
  }
  INFO("Build H0 - ");
  H0.setFromTriplets(tripletList.begin(), tripletList.end());
  if( DEBUG > 5) INFO(H0);
}

template<typename Tnum>
void FermiHubbard<Tnum>::BuildOneBodyTerms(){}

template<typename Tnum>
void FermiHubbard<Tnum>::BuildTwoBodyTerms( const Tnum U )
{
  tripletList.clear();
  std::vector< std::vector<size_t> > IndexUup;
  std::vector< std::vector<size_t> > IndexUdn;
  size_t point_up, point_dn;
  Tnum g;
  for (size_t i = 0; i < getL(); i++) {
    std::vector<size_t> work;
    point_up = 0;
    for (size_t cnt_up = 0; cnt_up < getEachHilbertSpace(0); cnt_up++) {
      if ( getTotalN() <= getL() ) {
        if ( btest(Basis.at(0).at(cnt_up), i) ) {
          point_up++;
          work.push_back(cnt_up);
        }
      }
      else{
        if ( !(btest(Basis.at(0).at(cnt_up), i)) ) {
          point_up++;
          work.push_back(cnt_up);
        }
      }
    }
    IndexUup.push_back(work);
    work.clear();
  }
  if ( getN(0) == getN(1) ) {
    point_dn = point_up;
    IndexUdn = IndexUup;
  }
  else{
    for (size_t i = 0; i < getL(); i++) {
      std::vector<size_t> work;
      point_dn = 0;
      for (size_t cnt_dn = 0; cnt_dn < getEachHilbertSpace(1); cnt_dn++) {
        if ( getTotalN() <= getL() ) {
          if ( btest(Basis.at(1).at(cnt_dn), i) ) {
            point_dn++;
            work.push_back(cnt_dn);
          }
        }
        else{
          if ( !(btest(Basis.at(1).at(cnt_dn), i)) ) {
            point_dn++;
            work.push_back(cnt_dn);
          }
        }
      }
      IndexUdn.push_back(work);
      work.clear();
    }
  }
  // INFO("point " << point_up << " " << point_dn);
  if ( getTotalN() > getL() ) {
    g = U * (Tnum)( getTotalN() - getL() );
    for (size_t cnt = 0; cnt < getTotalHilbertSpace(); cnt++) {
      tripletList.push_back(MatrixElemT(cnt, cnt, g));
      if (DEBUG > 5) INFO("1 " << cnt << " " << cnt << " " << g);
    }
  }
  if ( getN(0) == getN(1) ) {
    if ( getTotalN() <= getL() ) {
      g = U * (Tnum)getN(0);
    }
    else{
      g = U * (Tnum)( getL() - getN(0) );
    }
    for (size_t cnt_up = 0; cnt_up < getEachHilbertSpace(0); cnt_up++) {
      size_t id = getIndices2(cnt_up, cnt_up);
      tripletList.push_back(MatrixElemT(id, id, g));
      if (DEBUG > 5) INFO("2 " << id << " " << id << " " << g);
    }
  }
  for (size_t i = 0; i < getL(); i++) {
    for (size_t up = 0; up < point_up; up++) {
      for (size_t dn = 0; dn < point_dn; dn++) {
        if ( getN(0) == getN(1) ) {
          if ( up != dn ) {
            size_t id = getIndices2(IndexUup.at(i).at(up), IndexUdn.at(i).at(dn) );
            tripletList.push_back(MatrixElemT(id, id, U));
            if (DEBUG > 5) INFO("3 " << id << " " << id << " " << U);
          }
        }
        else {
          size_t id = getIndices2(IndexUup.at(i).at(up), IndexUdn.at(i).at(dn) );
          tripletList.push_back(MatrixElemT(id, id, U));
          if (DEBUG > 5){
            INFO(i << " " << IndexUup.at(i).at(up) << " " << IndexUdn.at(i).at(dn));
            INFO(Basis[0].at(IndexUup.at(i).at(up)) << " " << Basis[1].at(IndexUdn.at(i).at(dn)));
            INFO("4 " << id << " " << id << " " << U);
          }
        }
      }
    }
  }
  INFO("Build HTwo - ");
  HTwo.setFromTriplets(tripletList.begin(), tripletList.end());
  if( DEBUG > 5) INFO(HTwo);
}

template<typename Tnum>
void FermiHubbard<Tnum>::ConstructTotalHamiltonian()
{
  Htot = H0 + HOne + HTwo;
  if (DEBUG > 5) {
    INFO(Htot);
  }
}

template<typename Tnum>
void FermiHubbard<Tnum>::eigh()const
{
  VectorType Vec = RealVectorType::Random( getTotalHilbertSpace() );
  RealType Val = 0.0e0;
  size_t max_iter = 200;
  double err_tol = 1.0E-7;
  if ( LanczosEV( getTotalHilbertSpace(), Htot, Vec, Val, max_iter, err_tol) ){
    INFO("Eigenvalue = " << Val << " takes " << max_iter << " iterations.");
  }
  else{
    INFO("Lanczos is not converged!");
  }
}

template<typename Tnum>
void FermiHubbard<Tnum>::LanczosExpH(const size_t Order)const
{}

template class FermiHubbard<RealType>;
template class FermiHubbard<ComplexType>;
