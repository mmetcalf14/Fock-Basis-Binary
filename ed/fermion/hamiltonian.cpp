#include "fermion/hamiltonian.h"
#include "fermion/bitwise.h"
// #include "fermion/search.h"

template<typename Tnum>
FermiHubbard<Tnum>::FermiHubbard(size_t L, size_t Nup, size_t Ndn,
  Tnum t_up, Tnum t_dn):
FermionBasis( L, Nup, Ndn )
{
  dim = TotalHilbertSpace;
  // INFO(dim);
  H0.resize(dim, dim);
  H0.reserve(3*dim);
  HI.resize(dim, dim);
  HI.reserve(3*dim);
  Htot.resize(dim, dim);
  Htot.reserve(3*dim);
  // INFO(Hamiltonian.outerSize());
  tripletList.clear();

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
  HI.setZero();
  Htot.setZero();
}

template<typename Tnum>
void FermiHubbard<Tnum>::Build1DHoppingTerms()
{
  size_t rid, cid;
  uint64_t count_up = HilbertSpace[0];
  uint64_t count_dn = HilbertSpace[1];
  size_t lid = 0, pid = 0;//l and p's index
  /* For up-spin */
  for( auto l : Basis[0] ){
    for (size_t i = 0; i < lattice_points - 1; i++) {
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
          auto value = -t * pow(-1.0, btest(Basis[1][cnt], i) );
          tripletList.push_back(MatrixElemT(rid, cid, value));
          tripletList.push_back(MatrixElemT(cid, rid, value));
          INFO(lid << " " << pid << " " << rid << " " << cid << " " << value);
        }
      }
    }
  }
  for( auto l : Basis[1] ){
    for (size_t i = 0; i < lattice_points; i++) {
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
        for (size_t cnt = 0; cnt < count_dn; cnt++) {
          rid = getIndices2(cnt, lid);//lid * count_up + cnt;
          cid = getIndices2(cnt, pid);//pid * count_up + cnt;
          Tnum value = -t * pow(-1.0, btest(Basis[0][cnt], j) );//NOTE: there is a convention here.
          tripletList.push_back(MatrixElemT(rid, cid, value));
          tripletList.push_back(MatrixElemT(cid, rid, value));
          INFO(lid << " " << pid << " " << rid << " " << cid << " " << value);
        }
      }
    }
  }
  H0.setFromTriplets(tripletList.begin(), tripletList.end());
  INFO(H0);
}

template<typename Tnum>
void FermiHubbard<Tnum>::BuildOneBodyTerms(){}

template<typename Tnum>
void FermiHubbard<Tnum>::BuildTwoBodyTerms(){}

template class FermiHubbard<RealType>;
template class FermiHubbard<ComplexType>;
