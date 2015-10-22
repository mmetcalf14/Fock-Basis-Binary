//
//  Basis.cpp
//  Fock Basis Binary
//
//  Created by mekena McGrew on 10/22/15.
//  Copyright © 2015 Mekena Metcalf. All rights reserved.
//

#include <iostream>
#include <stdio.h>
#include <cmath>
#include "Basis.h"


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
    size_t minrange_up = 0, maxrange_up = 0;
    size_t minrange_down = 0, maxrange_down = 0;
    
    //cout << "min up preloop " << minrange_up << " max up " << maxrange_up << endl;
    
    for (int i = 1; i <= Nup; i++)
    {
        
        minrange_up += pow(2,(i-1));
        maxrange_up += pow(2,(L-i));
        
    }
    // cout << "min up " << minrange_up << " max up " << maxrange_up << endl;
    
    for (int i = 1; i <= Ndn; i++)
    {  minrange_down += pow(2,(i-1));
        maxrange_down += pow(2,(L-i));
    }
    //cout << "Ranges down" << minrange_down << " " << maxrange_down <<endl;
    index_up.reserve(maxrange_up+1);
    std::vector<size_t> work (maxrange_up+1, 0);
    count_up = 0;
    for (size_t i = minrange_up; i <= maxrange_up; i++) //create spin up basis and vectors
    {
        int nbit = 0;
        for(size_t j = 0; j < L; j++)
        {
            
            if (MY_bittest(i,j))
            {
                nbit++;
            }
        }
        
        if (nbit == Nup)
        {
            count_up++;
            basis_up.push_back(i);
            work.at(i) = count_up;
            //cout << i << " " << work[i] << endl;
        }
        
    }
    index_up = work;
    
    //    cout << "Basis up \n";
    //        for(int i = 0; i < maxrange_down; i++)
    //        { cout << basis_up[i]<<endl;}
    //    cout << "Index up \n";
    //    for(int i = 0; i < maxrange_up; i++)
    //    { cout << index_up[i]<<endl;}
    
    
    std::cout << "End of basis up allocation \n";
    
    index_dn.reserve(maxrange_down+1);
    std::vector<size_t> nowork (maxrange_down+1, 0);
    count_dn = 0;
    //program quitting right here
    // cout << "Testing basis down vector\n";
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
            //std::cout << i<< " " << count_dn << std::endl;//this is right
            basis_down.push_back(i);
            //cout << i << " " << basis_down[] << endl;//below are all wrong?? but it worked before?
            nowork.at(i) = count_dn;
            //cout << i << " " << nowork[i] << endl;
        }
        
    }
    index_dn = nowork;
    //    cout << "Basis down \n";
    //    for(int i = 0; i < maxrange_down; i++)
    //    { cout << basis_down[i]<<endl;}
    //    cout << "Index down \n";
    //    for(int i = 0; i <= maxrange_down; i++)
    //    { cout << index_dn[i]<<endl;}
    std::cout << "End of basis dn allocation \n";
    
}

