//
//  Hamiltonian_Template.h
//  Fock Basis Binary
//
//  Created by mekena McGrew on 6/29/15.
//  Copyright (c) 2015 Mekena Metcalf. All rights reserved.
//

#include <iostream>

#include <cmath>
#include <complex>
#include </usr/include/Eigen/Eigen>
#include </usr/include/Eigen/Sparse>
#include </usr/include/Eigen/Core>
#include "Basis.h"

#ifndef Fock_Basis_Binary_Hamiltonian_Template_h
#define Fock_Basis_Binary_Hamiltonian_Template_h


//using namespace Eigen;



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
    void BuildHopHam(size_t count, size_t count_opp, std::vector<size_t> basis, std::vector<size_t> index, SpMat &HopHam);
    void BuildHopHam_up();
    void BuildHopHam_dn();
    void Interaction_Index();
    void Build_Interactions();
    void BaseInteraction();
    
    void ClearTriplet();
    void QuenchU(double _Uquench);
    
    void HopMatrix_Build();
    void IntMatrix_Build();
    void Set_Const(double t_1, double t_2);
    void Save_Ham();//input can be filename from main cpp
    void Total_Ham();
};




#endif
