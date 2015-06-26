#ifndef __BASIS_HPP__
#define __BASIS_HPP__
#include <iostream>
#include <vector>

class FermionBasis {
protected:
  size_t lattice_points;
  size_t number_species;
  std::vector<size_t> NSpecies;
  std::vector<uint64_t> HilbertSpace;
  uint64_t TotalHilbertSpace;
  std::vector< std::vector<size_t> > Basis;
  std::vector< std::vector<size_t> > Index;
  std::vector<size_t> Indices;

public:
  FermionBasis (size_t L, size_t Nup, size_t Ndn);
  virtual ~FermionBasis ();
  void CalculateHilbertSpace();
  void BuildBasis();
  void BuildIndices2();
  void PrintIndices2();
  inline size_t getIndices2(const size_t &idx0, const size_t &idx1){
    return (idx1 * HilbertSpace[0] + idx0);
  };
  inline size_t getIndices2fromBasis(const size_t &basis0, const size_t &basis1){
    return (Index[1][basis1] * HilbertSpace[0] + Index[0][basis0]);
  };
  inline uint64_t getTotalHilbertSpace(){return TotalHilbertSpace;};
  inline uint64_t getEachHilbertSpace(const size_t id){return HilbertSpace[id];};

};
#endif// __BASIS_HPP__
