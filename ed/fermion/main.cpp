#include <iostream>
#include "hamiltonian.h"

int main (int argc, char const *argv[])
{
    FermiHubbard<RealType> FHM(3, 2, 1);
    // FermiHubbard<ComplexType> FHM(3, 2, 1);

    return 0;
}
