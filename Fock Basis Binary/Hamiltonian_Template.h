//
//  Hamiltonian_Template.h
//  Fock Basis Binary
//
//  Created by mekena McGrew on 6/29/15.
//  Copyright (c) 2015 Mekena Metcalf. All rights reserved.
//

#include <iostream>
#include <vector>
#include "/usr/local/include/c++/4.9.2/Eigen/Eigen"
#include "/usr/local/include/c++/4.9.2/Eigen/Sparse"

#ifndef Fock_Basis_Binary_Hamiltonian_Template_h
#define Fock_Basis_Binary_Hamiltonian_Template_h
 class Basis //declare class for basis creation
{
private:  //have to talk to values through constructor function
    

protected:
    std::vector<size_t> basis_up;
    std::vector<size_t> basis_down;
    std::vector<size_t> index_up;
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

//template <class T>
class Hamiltonian :public Basis //declare class for Hamiltonian matrices
{
private:
    
    int tbar;
    
    typedef Eigen::SparseMatrix<size_t> SpMat;

    typedef Eigen::Triplet<size_t> Tp;
    
    std::vector<Tp> TL_up;
    std::vector<Tp> TL_down;
//    TL_up.reserve(3); //put this in the function
//    TL_down.reserve(3);
    
    //Declaring Matrices
    SpMat HopHam_up;//declare dimension in function
    SpMat HopHam_down;
public:
    Hamiltonian( size_t _L, size_t _Nup, size_t _Ndn ):Basis(_L, _Nup, _Ndn){};
    void Set_tbar(int _tbar);
    void Set_Mat_Dim();
    void BuildHopHam_up();
    void BuildHopHam_dn();
    
    void Matrix_Build();
    
};


//Bitewise functions used on Fock states
size_t testbit(size_t m, int n)// m -> basis integer, n -> site
{
    size_t Eval;
    Eval = (m & (1 << n));
    return Eval;
}

size_t bitclr(size_t m,  int n) // set nth bit to zero
{
    size_t Clr_bit;
    Clr_bit = m & ~(1 << n);
    return Clr_bit;
}

size_t ibtset(size_t m,  int n)
{
    size_t New_State;
    New_State = m | (1 << n);
    return New_State;
}
#endif
