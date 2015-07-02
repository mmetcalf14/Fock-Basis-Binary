#include <iostream>
#include "hamiltonian.h"

int main (int argc, char const *argv[])
{
  size_t L = 4;
  size_t Nup = 2;
  size_t Ndn = 3;
  FermiHubbard<RealType> FHM(L, Nup, Ndn);
  // FermiHubbard<ComplexType> FHM(L, Nup, Ndn);
  RealType U = 1.0;
  FHM.BuildTwoBodyTerms(U);
  return 0;
}
