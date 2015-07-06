#include <iostream>
#include "hamiltonian.h"

int main (int argc, char const *argv[])
{
  size_t L = 12;
  size_t Nup = 6;
  size_t Ndn = 6;
  FermiHubbard<RealType> FHM(L, Nup, Ndn);
  // FermiHubbard<ComplexType> FHM(L, Nup, Ndn);
  RealType U = 0.0;
  FHM.BuildTwoBodyTerms(U);
  FHM.ConstructTotalHamiltonian();
  FHM.eigh();
  return 0;
}
