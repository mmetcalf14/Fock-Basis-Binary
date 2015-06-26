#include <iomanip>//fixed, setw
#include <cmath>//pow
#include <cassert>//assert
#include "EDType.h"
#include "fermion/basis.h"
#include "fermion/bitwise.h"

uint64_t choose( int n, int k);
int gcd(int n,int m){return m==0?n:gcd(m,n%m);}

uint64_t choose( int n, int k)
{
  uint64_t result = 1;
  for (uint64_t d = 1; d <= k; d++)
  {
    uint64_t num = n - k + d;
    uint64_t denom = d;
    uint64_t common_factor = gcd(result, denom);
    result /= common_factor;
    denom /= common_factor;
    assert(num % denom == 0);
    num /= denom;
    if (result > UINT64_MAX / num) {
      OVERFLOW_ERROR("Overflow");
    }
    result *= num;
  }
  return result;
}

FermionBasis::FermionBasis(size_t L, size_t Nup, size_t Ndn)
{
    lattice_points = L;
    number_species = 2;
    NSpecies.push_back(Nup);
    NSpecies.push_back(Ndn);
    CalculateHilbertSpace();
    BuildBasis();
    BuildIndices2();
    PrintIndices2();
}

FermionBasis::~FermionBasis()
{
  INFO("FermionBasis Destructor called.");
  Basis.clear();
  Index.clear();
}

void FermionBasis::CalculateHilbertSpace()
{
  TotalHilbertSpace = 1;
  for (size_t cnt = 0; cnt < number_species; cnt++) {
    uint64_t hs = choose(lattice_points, NSpecies[cnt]);
    HilbertSpace.push_back(hs);
    TotalHilbertSpace *= hs;
    INFO("Hibert Space of species-" << cnt << " is " << HilbertSpace[cnt] << ".");
  }
  INFO("Total Hibert Space is " << TotalHilbertSpace << ".");
}

void FermionBasis::BuildBasis()
{
  std::vector<size_t> wbasis;
  for (size_t idx = 0; idx < number_species; idx++) {
    size_t minrange = 0;
    size_t maxrange = 0;
    for (size_t cnt = 1; cnt < NSpecies[idx] + 1; cnt++) {
      minrange += pow(2, cnt - 1);
      maxrange += pow(2, lattice_points - cnt);
    }
    // INFO(minrange << " " << maxrange);
    size_t sub_cnt = 0;
    std::vector<size_t> windex (maxrange+1, 0);
    for (size_t cnt1 = minrange; cnt1 <= maxrange; cnt1++) {
      int nbit = 0;
      for (size_t cnt2 = 0; cnt2 < lattice_points; cnt2++) {
        if ( btest(cnt1, cnt2) ){
          nbit += 1;
        }
      }
      if (nbit == NSpecies[idx]){
        sub_cnt += 1;
        wbasis.push_back(cnt1);
        windex.at(cnt1) = sub_cnt - 1;//NOTE: I start from 0!
        // INFO(cnt1 << " " << sub_cnt << " " << windex.at(cnt1));
      }
    }
    assert( wbasis.size() == HilbertSpace.at(idx) );
    Basis.push_back(wbasis);
    Index.push_back(windex);
    wbasis.clear();
  }
  /*NOTE: check basis consistancy */
  size_t nBasis = 1;
  for (size_t cnt = 0; cnt < number_species; cnt++) {
    nBasis *= Basis[cnt].size();
  }
  if ( nBasis != TotalHilbertSpace ) {
    RUNTIME_ERROR("Hilbert space is not consistence.");
  }
}

void FermionBasis::BuildIndices2()
{
  for ( auto &b1 : Basis[1] ){
    for ( auto &b0 : Basis[0] ){
      Indices.push_back( getIndices2fromBasis( b0, b1 ) );
    }
  }
}

void FermionBasis::PrintIndices2()const
{
  size_t cnt = 0;
  for ( auto &b1 : Basis[1] ){
    for ( auto &b0 : Basis[0] ){
      INFO(std::fixed << std::setw(8) <<
        "   Up: " << std::setw(4) << b0 << " " << std::setw(4) << Index[0][b0] << std::setw(8) <<
        " Down: " << std::setw(4) << b1 << " " << std::setw(4) << Index[1][b1] << std::setw(8) <<
        "Index: " << std::setw(6) << Indices[cnt]);
      cnt += 1;
    }
  }
}
