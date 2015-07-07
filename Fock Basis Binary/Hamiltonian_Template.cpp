//
//  Hamiltonian_Template.cpp
//  Fock Basis Binary
//
//  Created by mekena McGrew on 6/30/15.
//  Copyright (c) 2015 Mekena Metcalf. All rights reserved.
//

#include <stdio.h>
#include "Hamiltonian_Template.h"

//void Hamiltonian::Set_tbar(int _tbar)
//{ tbar = _tbar;}

Basis::Basis(size_t _L, size_t _Nup, size_t _Ndn)
//:L(_L), Nup(_Nup), Ndn(_Ndn)
{
    L = _L;
    Nup = _Nup;
    Ndn = _Ndn;
    BuildBasis();
}

void Basis::BuildBasis()
{
//    std::cout << L << " " << Nup << " " << Ndn << std::endl;
    size_t minrange_up, maxrange_up;
    size_t minrange_down, maxrange_down;
    for (int i = 1; i <= Nup; i++)
    {  minrange_up += pow(2,(i-1));
        maxrange_up += pow(2,(L-i));
    }
    
    for (int i = 1; i <= Ndn; i++)
    {  minrange_down += pow(2,(i-1));
        maxrange_down += pow(2,(L-i));
    }
    index_up.reserve(maxrange_up);
    count_up = 0;
    for (size_t i = minrange_up; i <= maxrange_up; i++) //create spin up basis and vectors
    {
        int nbit = 0;
        for(size_t j = 0; j < L; j++)
        {std::cout << "When does it break? \n";
            if (MY_bittest(i,j))
            {
                nbit++;
            }
        }
        
        if (nbit == Nup)
        {
            count_up++;
            basis_up.push_back(i);
//            index_up.at(i) = count_up;
        }
        
    }
    
    
    count_dn = 0;
    for (size_t i = minrange_down; i <= maxrange_down; i++) //create spin down basis and index
    {
        int nbit = 0;
        for(size_t j = 0; j < L; j++)
        {
            if (MY_bittest(i,j))
            {
                nbit++;
            }
        }
        
        if (nbit == Ndn && Ndn != 0)
        {
            count_dn++;
            std::cout << i<< " " << count_dn << std::endl;
            basis_down.push_back(i);
//            index_down.at(i) = count_down;
        }
        
    }

}


void Hamiltonian::BuildHopHam_up()
{
        for(size_t bs = 0; bs < count_up; bs++)
        {
            size_t p_bas = basis_up[bs];
            size_t p_ind = index_up[p_bas];
            //cout << "We are acting on basis, " << p_bas << " with index, "<< p_ind << endl;
            for(size_t i = 0; i < (L-1); i++)
            {
                size_t l_bas = MY_bitset(MY_bitclr(p_bas,i),i+1);
                size_t l_ind = index_up[l_bas];
    
                if(l_bas != p_bas && l_ind != 0)
                { //cout << "The hopping term gives the new state: " << l_bas << " With index: " << l_ind <<endl;
    
    
                if(count_dn > 0 )
                {
                    for(size_t k = 1; k <= count_dn; k++)
                    {
                        size_t r = ((k-1)*count_up) + p_ind;
                        size_t s = ((k-1)*count_up) + l_ind;
                        size_t val = -tbar ;//* pow(-1,testbit(k,i));
                        TL_up.push_back(Tp((r-1),(s-1),val));
    
                    }
                }
                else
                {   int val = -tbar;
                    TL_up.push_back(Tp(p_ind-1,l_ind-1, val ));
                   
                }
                }
            }
        }
 
}

void Hamiltonian::Set_Mat_Dim()
{
    int Tot_base = count_up*count_dn; //should this be size t or int for matrix dim
    SpMat HopHam_down(Tot_base, Tot_base);
    SpMat HopHam_up(Tot_base, Tot_base);
}

void Hamiltonian::Matrix_Build()
{
    HopHam_up.setFromTriplets(TL_up.begin(), TL_up.end());
    HopHam_down.setFromTriplets(TL_down.begin(), TL_down.end());
    
    std::cout << "The Hamiltonian for up spin is:" << HopHam_up << std::endl;
}







