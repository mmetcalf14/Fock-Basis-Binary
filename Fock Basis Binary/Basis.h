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
private:  //have to talk to values through constructor function
    
    
protected:
    std::vector<size_t> basis_up;
    std::vector<size_t> basis_down;
    std::vector<size_t> index_up;
    std::vector<size_t> index_dn;
    size_t count_up, count_dn;
    size_t L, Nup, Ndn;
    //how to dynamically allocate. How to talk to private data?
public:
    Basis();
    Basis(size_t _L, size_t Nup, size_t Ndn);
    void BuildBasis();
    inline size_t getNsite()const{return L;};
    inline void changeNsite(size_t New_L){L = New_L;};
    
};
#endif /* Basis_h */
