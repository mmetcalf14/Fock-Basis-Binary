//
//  Hamiltonian_Template.h
//  Fock Basis Binary
//
//  Created by mekena McGrew on 6/29/15.
//  Copyright (c) 2015 Mekena Metcalf. All rights reserved.
//

#include <iostream>
#include <vector>
#include "/usr/include/Eigen/Eigen"
#include "/usr/include/Eigen/Sparse"
#include "/usr/include/Eigen/Core"

#ifndef Fock_Basis_Binary_Hamiltonian_Template_h
#define Fock_Basis_Binary_Hamiltonian_Template_h


//using namespace Eigen;

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

//template <class T>
class Hamiltonian :public Basis //declare class for Hamiltonian matrices
{
    friend class Lanczos_Diag;
private:
    
    double tbar;
    double U;
    
    
    
    typedef Eigen::SparseMatrix<double> SpMat;

    typedef Eigen::Triplet<double> Tp;
    
    std::vector<Tp> TL_up;
    std::vector<Tp> TL_down;
    std::vector<Tp> TL_Ubase;
    //    TL_up.reserve(3); //put this in the function
//    TL_down.reserve(3);
    
    //Interaction Index Matrices
    Eigen::MatrixXd IndexU_up;
    Eigen::MatrixXd IndexU_dn;
    int point_up;
    int point_dn;
    

protected:
    int Tot_base;
    //Declaring Matrices
    SpMat HopHam_up;//declare dimension in function
    SpMat HopHam_down;
    SpMat Ham_Interact;
    SpMat Ham_Tot;
    
public:
    Hamiltonian( size_t _L, size_t _Nup, size_t _Ndn ):Basis(_L, _Nup, _Ndn){};//is this costructor or
    //Hamiltonian Functions
    void Set_Mat_Dim();
    void BuildHopHam_up();
    void BuildHopHam_dn();
    void Interaction_Index();
    void Build_Interactions();
    void BaseInteraction();
    
    void HopMatrix_Build();
    void IntMatrix_Build();
    void Set_Const(double _tbar, double _U);
    void Save_Ham();//input can be filename from main cpp
    void Total_Ham();
};

class Lanczos_Diag //:public Hamiltonian
{
private:
    
    //typedef Eigen::SparseMatrix<double> SpMat;
    
    Eigen::MatrixXd TriDiag;
    Eigen::MatrixXd K_Mat;
    Eigen::VectorXd Lanczos_Vec;
    Eigen::VectorXd Lanczos_Vec_Temp;
    Eigen::VectorXd r_vec;
    
    double alpha;
    double beta;
    
public:
    
    Lanczos_Diag(const Hamiltonian){};//Program not accepting this constructor::SEE ERROR
    void Set_Mat_Dim_LA(int Tot_base);
    void Random_Vector();
    
    template <typename Derived>
    void Diagonalize(const Eigen::SparseMatrixBase<Derived> &Ham_Tot);
    //why isn't it recognizing the template?
    void Get_Gstate();
    
    
};


//Bitewise functions definitions used on Fock states
//Binary Function algorithm
inline size_t MY_bittest(size_t m, size_t n)// m -> basis integer, n -> site
{
    size_t Eval; //if Eval is size_t I get a totally wrong number compared to int
    //seg fault occurring regardless of whether return value is correct or incorrect
    Eval = (m & (1 << n));
    return Eval;
}

inline size_t MY_bitclr(size_t m,  size_t n) // set nth bit to zero
{
    size_t Clr_bit;
    Clr_bit = m & ~(1 << n);
    return Clr_bit;
}

inline size_t MY_bitset(size_t m,  size_t n)
{
    size_t New_State;
    New_State = m | (1 << n);
    return New_State;
}

#endif
