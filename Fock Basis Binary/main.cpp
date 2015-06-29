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
#include <algorithm>
#include "/usr/local/include/c++/4.9.2/Eigen/Eigen"
#include "/usr/local/include/c++/4.9.2/Eigen/Dense"
#include "/usr/local/include/c++/4.9.2/Eigen/Eigenvalues"
#include "/usr/local/include/c++/4.9.2/Eigen/Sparse"
#include "/usr/local/include/c++/4.9.2/Eigen/StdVector"

int bittest(int m,  int n);
int bitclr(int m,  int n);
int ibtset(int m,  int n);


int main(int argc, const char * argv[])
{
    using namespace std;
    using namespace Eigen;
    
    
    int maxrange_up = 0;
    int minrange_up = 0;
    int maxrange_down = 0;
    int minrange_down = 0;
    int nbit = 0;
    int Nup = 2;
    int Ndown = 1;
    int Nsite = 3;
    int count_up = 0;
    int count_down = 0;
    
    double tbar = 1.0;
    
    typedef SparseMatrix<double> SpMat;
    
    typedef Triplet<double> T;
    
    vector<T> TL_up;
    vector<T> TL_down;
    
    TL_up.reserve(3);
    TL_down.reserve(3);
    
    for (int i = 1; i <= Nup; i++)
    {  minrange_up += pow(2,(i-1));
        maxrange_up += pow(2,(Nsite-i));
    }
    
    for (int i = 1; i <= Ndown; i++)
    {  minrange_down += pow(2,(i-1));
        maxrange_down += pow(2,(Nsite-i));
    }
    cout << minrange_down << " " << maxrange_down << endl;
    vector<int> basis_up;
    vector<int> basis_down;
//    SparseVector<int> basis_sig(maxrange);
   // SparseVector<int> index_sig(maxrange+1);//need to choose if i=3 this is element 4 in vector since
    //starting index is 0
    vector<int> index_up (maxrange_up+1);
    vector<int> index_down (maxrange_down+1);
    
    
    for (int i = minrange_up; i <= maxrange_up; i++) //create spin up basis and vectors
    {
        nbit = 0;
        for(int j = 0; j < Nsite; j++)
        {
            if (bittest(i,j))
            {
                nbit++;            
            }
        }
        
        if (nbit == Nup)
        {
            count_up++;
            basis_up.push_back(i);
            index_up.at(i) = count_up;
        }
        
    }
    

    
    for (int i = minrange_down; i <= maxrange_down; i++) //create spin down basis and index
    {
        nbit = 0;
        for(int j = 0; j < Nsite; j++)
        {
            if (bittest(i,j))
            {
                nbit++;
            }
        }
        
        if (nbit == Ndown && Ndown != 0)
        {
            count_down++;
            cout << i<< " " << count_down<< endl;
            basis_down.push_back(i);
            index_down.at(i) = count_down;
        }
        
    }
    
    int Tot_base = count_up*count_down;
    
    
    //Declaring Matrices
    SpMat H_up(Tot_base,Tot_base);
    SpMat H_down(Tot_base,Tot_base);

    
    //Display Basis and Index vectors
    for( int cnt=0; cnt < basis_up.size(); cnt++)
    {
        cout << basis_up[cnt] << endl;
    }
    cout << endl;
    for( int cnt=0; cnt < index_up.size(); cnt++)
    {
        cout << index_up[cnt] << endl; //yet it does work here
    }
    
    
    //Creation of the off diagonal elements through hopping
    
    for(int bs = 0; bs < count_up; bs++)
    {
        int p_bas = basis_up[bs];
        int p_ind = index_up[p_bas];
        //cout << "We are acting on basis, " << p_bas << " with index, "<< p_ind << endl;
        for(int i = 0; i < (Nsite-1); i++)
        {
            int l_bas = ibtset(bitclr(p_bas,i),i+1);
            int l_ind = index_up[l_bas];
            
            if(l_bas != p_bas && l_ind != 0)
            { //cout << "The hopping term gives the new state: " << l_bas << " With index: " << l_ind <<endl;
            

            if(count_down > 0 )
            {
                for(int k = 1; k <= count_down; k++)
                {
                    cout << "CHecking loop \n";
                    int r = ((k-1)*count_up) + p_ind;
                    int s = ((k-1)*count_up) + l_ind;
                    int val = -tbar * pow(-1,bittest(k,i));
                    TL_up.push_back(T((r-1),(s-1),val));
                  
                }
            }
            else
            {   int val = -tbar;
                TL_up.push_back(T(p_ind-1,l_ind-1, val ));
               
            }
            }
        }
    }
    
    H_up.setFromTriplets(TL_up.begin(), TL_up.end());

    cout << "The Hamiltonian for up spin is:" << H_up << endl;
    
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



