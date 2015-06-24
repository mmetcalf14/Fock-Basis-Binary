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

int bittest(int m,  int n);
int bitclr(int m,  int n);
int ibtset(int m,  int n);
int esign(int m, int n);

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
    int Ndown = 0;
    int Nsite = 3;
    int count_up = 0;
    int count_down = 0;
    
    double tbar = 1.0;
    MatrixXd H_up;
    MatrixXd H_down;
    
    for (int i = 1; i <= Nup; i++)
    {  minrange_up += pow(2,(i-1));
        maxrange_up += pow(2,(Nsite-i));
    }
    
    for (int i = 1; i <= Ndown; i++)
    {  minrange_down += pow(2,(i-1));
        maxrange_down += pow(2,(Nsite-i));
    }
    
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
            //cout << i << " " << j << " " << bittest(i,j) << endl;
            if (bittest(i,j))
            {
                nbit++;            
            }
        }
        
        if (nbit == Nup)
        {
            count_up++;
            cout << i<< " " << count_up<< endl;
            basis_up.push_back(i);
            index_up.at(i) = count_up;
        }
        
    }
    
    for (int i = minrange_down; i <= maxrange_down; i++) //create spin down basis and index
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
        
        if (nbit == Ndown)
        {
            count_down++;
            cout << i<< " " << count_down<< endl;
            basis_down.push_back(i);
            index_down.at(i) = count_down;
        }
        
    }
    
    //Display Basis and Index vectors
    for( int cnt=0; cnt < basis_up.size(); cnt++){
        cout << basis_up[cnt] << endl;
    }
    cout << endl;
    for( int cnt=0; cnt < index_up.size(); cnt++){
        cout << index_up[cnt] << endl; //yet it does work here
    }
    
    
    //Testing mapping process of Fock states
    int Clr = bitclr(3,0);
    //cout << Clr << endl;
    int New_Int = ibtset(Clr, 2);
    cout << "This is the New Int: " << New_Int << endl;
    
    cout << index_up[6] << endl;
    int Index_Sig = index_up.at(New_Int);
    cout << "This is the index: " << Index_Sig << endl;
    int New_State = basis_up.at(Index_Sig-1);
    cout << "This is the new state: " << New_State << endl;
    
    //statement for esign
    // if ((i-j)%2 == 0)
    // {esign(i,j)};

    
    
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

int esign(int m, int n)
{
    if((m-n)%2 == 0)
    {return -1;}
    else
        return 0;
}

