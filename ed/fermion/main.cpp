#include <iostream>
#include "hamiltonian.h"

int main (int argc, char const *argv[])
{
  size_t L = 4;
  size_t Nup = 3;
  size_t Ndn = 2;
  FermiHubbard<RealType> FHM(L, Nup, Ndn);
  // FermiHubbard<ComplexType> FHM(L, Nup, Ndn);

    return 0;
}
