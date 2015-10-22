//
//  Basis.h
//  Fock Basis Binary
//
//  Created by mekena McGrew on 10/22/15.
//  Copyright © 2015 Mekena Metcalf. All rights reserved.
//

#ifndef Basis_h
#define Basis_h

#include <vector>
#include "Bitwise_Functions.h"

class Basis //declare class for basis creation
{
protected:

    std::vector<size_t> basis_up;
    std::vector<size_t> basis_down;
    std::vector<size_t> index_up;
    std::vector<size_t> index_dn;
    size_t count_up, count_dn;
    size_t L, Nup, Ndn;

public:
    Basis();
    Basis(size_t _L, size_t Nup, size_t Ndn);
    void BuildBasis();
    void CreateBasis(size_t N, std::vector<size_t> &basis, std::vector<size_t> &index);
    inline size_t getNsite()const{return L;};
    inline size_t TotalIndex(size_t id1, size_t id2)const{
      return id2 * count_up + id1;
    };
};
#endif /* Basis_h */
