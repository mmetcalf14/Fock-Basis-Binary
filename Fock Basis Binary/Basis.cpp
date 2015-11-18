//
//  Basis.cpp
//  ED Hubbard Hamiltonian
//
//  Created by mekena McGrew on 11/4/15.
//  Copyright © 2015 Mekena Metcalf. All rights reserved.
//

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <cassert>
#include "Basis.h"


Basis::Basis(size_t _L, size_t _Nup, size_t _Ndn)
{
    L = _L;
    Nup = _Nup;
    Ndn = _Ndn;
    
    CreateBasis(Nup, basis_up, index_up);
    count_up = basis_up.size();
    CreateBasis(Ndn, basis_down, index_dn);
    count_dn = basis_down.size();
    
    
}



void Basis::CreateBasis(size_t N, std::vector<size_t> &basis,
                        std::vector<size_t> &index)
{
    size_t count = 0;
    size_t minrange = 0;
    size_t maxrange = 0;
    
    for (int i = 1; i <= N; i++)
    {
        minrange += pow(2,(i-1));
        maxrange += pow(2,(L-i));
    }
    
    index.reserve(maxrange+1);
    std::vector<size_t> work (maxrange+1, 0);
    for (size_t i = minrange; i <= maxrange; i++) //create basis and vectors
    {
        int nbit = 0;
        for(size_t j = 0; j < L; j++)
        {
            if (MY_bittest(i,j))
            {
                nbit++;
            }
        }
        if (nbit == N)
        {
            basis.push_back(i);
            work.at(i) = count;
            count++;
        }
    }
    index = work;
    assert( basis.size() == count );
}
