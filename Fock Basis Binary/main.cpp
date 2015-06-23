//
//  main.cpp
//  Fock Basis Binary
//
//  Created by mekena McGrew on 6/15/15.
//  Copyright (c) 2015 Mekena Metcalf. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <vector>
#include "/usr/local/include/c++/4.9.2/Eigen/Eigen"
#include "/usr/local/include/c++/4.9.2/Eigen/Dense"
#include "/usr/local/include/c++/4.9.2/Eigen/Eigenvalues"
#include "/usr/local/include/c++/4.9.2/Eigen/Sparse"

int bittest(int m,  int n);
int bitclr(int m,  int n);
int ibtset(int m,  int n);

int main(int argc, const char * argv[])
{
    using namespace std;
    using namespace Eigen;
    
    
    int maxrange = 0;
    int minrange = 0;
    int nbit = 0;
    int Nsig = 2;
    int Nsite = 3;
    int count_sig = 0;
    
    for (int i = 1; i <= Nsig; i++)
    {  minrange += pow(2,(i-1));
        maxrange += pow(2,(Nsite-i));
    }
    
    std::vector<int> basis_sig;
//    SparseVector<int> basis_sig(maxrange);
    SparseVector<int> index_sig(maxrange+1);//need to choose if i=3 this is element 4 in vector since
    //starting index is 0
    
    
    
    for (int i = minrange; i <= maxrange; i++)
    {
        nbit = 0;
        for(int j = 0; j < Nsite; j++)
        {
            //cout << i << " " << j << " " << bittest(i,j) << endl;
            if (bittest(i,j))
            {
                nbit++;            
            }
        }
        
        if (nbit == Nsig)
        {
            count_sig++;
           // cout << count_sig << endl;
            basis_sig.push_back(i);
//            basis_sig.coeffRef(count_sig) = i;
            index_sig.coeffRef(i) = count_sig;
        }
    }
    for( int cnt=0; cnt < basis_sig.size(); cnt++){
        cout << basis_sig[cnt] << endl;
    }
    
    int Clr = bitclr(5,0);
    cout << Clr << endl;
    int New_Int = ibtset(Clr, 1);
    cout << "The Basis 101 -> 110 should be 6: " << New_Int << endl;
//    cout << "Basis vector: " << basis_sig <<endl; //tells you where non-zero entries are
    cout << "Index vector: " << index_sig << endl;
    
    return 0;
}

// Function returns TRUE is jth bit of int i is 1. FALSE if 0
int bittest(int m, int n)// m -> basis integer, n -> site
{
    int Eval;
    Eval = (m & (1 << n));
    return Eval;
}

int bitclr(int m,  int n) // set nth bit to zero
{
    int Clr_bit;
    Clr_bit = m & ~(1 << n);
    return Clr_bit;
}

int ibtset(int m,  int n)
{
    int New_State;
    New_State = m | (1 << n);
    return New_State;
}



