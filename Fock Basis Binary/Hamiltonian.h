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
//#include <Eigen/Eigen>
//#include <Eigen/Sparse>
//#include <Eigen/Core>
 #include </usr/include/Eigen/Eigen>
 #include </usr/include/Eigen/Sparse>
 #include </usr/include/Eigen/Core>
#include "Basis.h"

#ifndef Fock_Basis_Binary_Hamiltonian_Template_h
#define Fock_Basis_Binary_Hamiltonian_Template_h


//using namespace Eigen;



template <typename Tnum>
class Hamiltonian :public Basis //declare class for Hamiltonian matrices
{
    template<typename T> //any type of Lanczos is friend of Hamiltonian class
    friend class Lanczos_Diag;
private:

    typedef Eigen::SparseMatrix<Tnum> SpMat;
    typedef Eigen::Triplet<Tnum> Tp;

    Tnum J1; //A-B hopping
    Tnum J2; //B-A hopping
    Tnum U;
    int Tot_base;
    //Declaring Matrices
    SpMat HopHam_up;//declare dimension in function
    SpMat HopHam_down;
    SpMat Ham_Interact;
    SpMat Ham_Tot;

public:

    Hamiltonian( size_t _L, size_t _Nup, size_t _Ndn ):Basis(_L, _Nup, _Ndn){
      Tot_base = count_up * count_dn;
      std::cout << Tot_base;
      HopHam_down.resize(Tot_base, Tot_base);
      HopHam_up.resize(Tot_base, Tot_base);
      Ham_Interact.resize(Tot_base,Tot_base);
    };//is this costructor or

    //making public because too difficult to pass as friend object
    //Hamiltonian Functions
    // void Set_Mat_Dim();
    // void BuildHopHam_up();
    // void BuildHopHam_dn();
    // void Interaction_Index();
    // void Build_Interactions();
    // void BaseInteraction();

    void ClearTriplet();
    void QuenchU(Tnum _Uquench);

    void Set_Const(Tnum t_1, Tnum t_2, Tnum _U = (Tnum)0.0e0);
    void HopMatrix_Build();
    void BuildHopHam(int species, size_t count, size_t count_opp, std::vector<size_t> basis, std::vector<size_t> index, SpMat &HopHam);
    void IntMatrix_Build();
    void Save_Ham();//input can be filename from main cpp
    void Total_Ham();
};




#endif
