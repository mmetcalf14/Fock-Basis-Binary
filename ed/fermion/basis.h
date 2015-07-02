#ifndef __BASIS_H__
#define __BASIS_H__
#include <vector>
#include <numeric>

class FermionBasis {
private:
  size_t lattice_points;
  size_t number_species;
  uint64_t TotalHilbertSpace;
  std::vector<size_t> NSpecies;
  std::vector<uint64_t> HilbertSpace;

protected:
  std::vector< std::vector<size_t> > Basis;
  std::vector< std::vector<size_t> > Index;
  std::vector<size_t> Indices;

public:
  FermionBasis (size_t L, size_t Nup, size_t Ndn);
  virtual ~FermionBasis ();
  void CalculateHilbertSpace();
  void BuildBasis();
  void BuildIndices2();
  void PrintIndices2()const;
  inline size_t getL()const{return lattice_points;};
  inline size_t getNSpecies()const{return number_species;};
  inline size_t getTotalN()const{return std::accumulate(NSpecies.begin(), NSpecies.end(), 0);};
  inline size_t getN( const size_t species_id )const{return NSpecies.at(species_id);};
  inline uint64_t getTotalHilbertSpace()const{return TotalHilbertSpace;};
  inline uint64_t getEachHilbertSpace(const size_t id)const{return HilbertSpace[id];};
  inline size_t getIndices2(const size_t &idx0, const size_t &idx1)const{
    return (idx1 * HilbertSpace[0] + idx0);
  };
  inline size_t getIndices2fromBasis(const size_t &basis0, const size_t &basis1)const{
    return (Index[1][basis1] * HilbertSpace[0] + Index[0][basis0]);
  };
};
#endif// __BASIS_H__
