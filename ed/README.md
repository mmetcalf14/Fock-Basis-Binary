## Eigen3 Sparse Matrix Example
 Inside `ed` folder, execute

  ```bash
  make eigen.app
  ```

This will compile `matrix/EigenMatrix.cpp`, in which initialize a sparse matrix and convert it into a dense one.

## Fermion Basis and Hamiltonian
This is a combined example for inheritance class and template class. `fermion/basis.h` defined the base class for `fermion/hamiltonian.h`. Within `class FermiHubbard`, we can access all public and protected members, but not private. The `class FermiHubbard` itself is a template class, which can be `RealType` or `ComplexType`.

Just type `make` inside `ed` folder, which compile `fermion/main.cpp`. In that main program, I only call constructor of `FermiHubbard` of `RealType` (or `ComplexType`). Those types are defined in `EDType.h`. This constructor does two things as following.

1. Build the Fermion basis from `class FermionBasis`, which are defined in `fermion/basis.h` and `fermion/basis.cpp`.

2. Using the build basis to construct the hopping Hamiltonian. Please see `fermion/hamiltonian.h` and `fermion/hamiltonian.cpp`.
