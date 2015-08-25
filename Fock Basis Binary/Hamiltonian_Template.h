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
    
    double J1; //A-B hopping
    double J2; //B-A hopping
    
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
    int Tot_base;
    //Declaring Matrices
    SpMat HopHam_up;//declare dimension in function
    SpMat HopHam_down;
    SpMat Ham_Interact;
    SpMat Ham_Tot;
    

protected:
//    int Tot_base;
//    //Declaring Matrices
//    SpMat HopHam_up;//declare dimension in function
//    SpMat HopHam_down;
//    SpMat Ham_Interact;
//    SpMat Ham_Tot;
    
public:
    
    Hamiltonian( size_t _L, size_t _Nup, size_t _Ndn ):Basis(_L, _Nup, _Ndn){};//is this costructor or
    
    //making public because too difficult to pass as friend object
    //Hamiltonian Functions
    void Set_Mat_Dim();
    void BuildHopHam_up();
    void BuildHopHam_dn();
    void Interaction_Index();
    void Build_Interactions();
    void BaseInteraction();
    
    void HopMatrix_Build();
    void IntMatrix_Build();
    void Set_Const(double t_1, double t_2, double _U);
    void Save_Ham();//input can be filename from main cpp
    void Total_Ham();
};

class Lanczos_Diag //:public Hamiltonian
{
private:
    
    //typedef Eigen::SparseMatrix<double> SpMat;
    
    Eigen::MatrixXd TriDiag;
    std::vector<Eigen::VectorXd> K_Mat;
    Eigen::VectorXd Lanczos_Vec;
    Eigen::VectorXd Lanczos_Vec_Temp;
    Eigen::VectorXd r_vec;
    Eigen::Matrix4d Test_Ham;
    Eigen::Vector4d Test_Lanczos;
    Eigen::MatrixXd Evec_Mat;
    Eigen::VectorXd G_state;
    Eigen::VectorXd Evec;
    //Eigen::VectorXd G_state_realspace;
    

    
    double alpha;
    double beta;
    int cnt;
    
public:
    
    std::vector<double> n_up;//public so they can be used in main program to write the file
    std::vector<double> n_dn;
    
    Lanczos_Diag(const Hamiltonian){};//Program not accepting this constructor::SEE ERROR
    //construct new,simple matrix to test algorithm and eigen values
    //and set Lanz vec to be one from analytical example
    void Lanczos_TestM(const Eigen::Matrix4d& _Test_Ham, const Eigen::Vector4d& _Test_Lanczos);
    void Set_Mat_Dim_LA(Hamiltonian& );//int Tot_base
   // void Random_Vector();
    
   // template <typename Derived>
    void Diagonalize(const Hamiltonian &Ham, Hamiltonian&);
    //why isn't it recognizing the template?
    void Get_Gstate();
    //void Test_Tri();
    void Gstate_RealSpace(Hamiltonian& ct_up, Hamiltonian& ct_dn, Hamiltonian& Nsite,const Hamiltonian& basis_up,const Hamiltonian& basis_dn);
    
    
};


//Bitewise functions definitions used on Fock states
//Binary Function algorithm
inline size_t MY_bittest(size_t m, size_t n)// m -> basis integer, n -> site
{
    size_t Eval;//if Eval is size_t I get a totally wrong number compared to int
    //seg fault occurring regardless of whether return value is correct or incorrect
    //std::cout << m << " " << n << std::endl;
    Eval = (m & (1 << n));//I haven't changed anything why not working all of a sudden?
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
