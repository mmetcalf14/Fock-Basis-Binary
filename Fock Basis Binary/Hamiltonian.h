//
//  Hamiltonian.h
//  ED Hubbard Hamiltonian
//
//  Created by mekena McGrew on 11/4/15.
//  Copyright © 2015 Mekena Metcalf. All rights reserved.
//




#ifndef Hamiltonian_h
#define Hamiltonian_h
//using namespace Eigen;

#include <iostream>

#include <cmath>
#include <complex>
#include </usr/local/include/Eigen/Eigen>
#include </usr/local/include/Eigen/Sparse>
#include </usr/local/include/Eigen/Core>
//#include </home/mmcgrew/Library/Eigen/Eigen>
//#include </home/mmcgrew/Library/Eigen/Sparse>
//#include </home/mmcgrew/Library/Eigen/Core>
#include "Basis.h"

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
    Tnum h;
    
    //Declaring Matrices
    SpMat HopHam_up;//declare dimension in function
    SpMat HopHam_down;
    SpMat SOCHam_up;
    SpMat SOCHam_dn;
    SpMat Zeeman_Ham;
    SpMat Ham_Interact;
    SpMat Ham_Tot;
    
    double Phi_t;
    
    std::vector<double> Harm_Trap;
    
public:
    int Tot_base;
    
    Hamiltonian( size_t _L, size_t _Nup, size_t _Ndn ):Basis(_L, _Nup, _Ndn){
        Tot_base = count_up * count_dn;
        
        HopHam_down.resize(Tot_base, Tot_base);
        HopHam_up.resize(Tot_base, Tot_base);
        Ham_Interact.resize(Tot_base,Tot_base);
        SOCHam_up.resize(Tot_base,Tot_base);
        SOCHam_dn.resize(Tot_base,Tot_base);
        Zeeman_Ham.resize(Tot_base,Tot_base);
    };//is this costructor or
    

    //making public because too difficult to pass as friend object
    //Hamiltonian Functions
    // void Set_Mat_Dim();
    // void BuildHopHam_up();
    // void BuildHopHam_dn();
    // void Interaction_Index();
    // void Build_Interactions();
    // void BaseInteraction();
    
    void ClearHopTriplet();
    void ClearInteractTriplet();
    void QuenchU(Tnum _Uquench);
    void GetHarmTrap(std::vector<double> HT);
    
    void Set_Const(Tnum t_1, Tnum t_2, Tnum _U = (Tnum)0.0e0);
    void HopMatrix_Build();
    void HopMatrix_Build_Peierls();
    void HopMatrix_Build_QPump();
    void HopMatrix_Build_Periodic();
    void HopMatrix_Build_Fibonacci();
    void HopMatrix_Build_PeriodicWithSOC(int site1, int site2, double gamma, double phase);
    void HopMatrix_Build_PeriodicNNNHop(int link_num, double gamma, double phase);
    void HopMatrix_Build_NNNHop_WReal(int link_num, double t2, double gamma, double phase);
    void HopMatrix_Build_HaldHam_NoGauge(int link_num, double gamma, double phase1, double phase2);
    void HopMatrix_Build_Zeeman_WSOC(int link_num, double gamma, double  phase_SOC, double phase_B, double B);
    void MakeCut(int cut);
    
    void BuildHopHam(int species, size_t count, size_t count_opp, std::vector<size_t> basis, std::vector<size_t> index, SpMat &HopHam, std::vector<double> HT = std::vector<double>(100,0.0e0));
    void BuildHopHam_Peierls(int species, size_t count, size_t count_opp, std::vector<size_t> basis, std::vector<size_t> index, SpMat &HopHam, std::vector<double> HT = std::vector<double>(100,0.0e0));
     void BuildHopHam_QPump(int species, size_t count, size_t count_opp, std::vector<size_t> basis, std::vector<size_t> index, SpMat &HopHam);
    void BuildHopHam_Periodic(int species, size_t count, size_t count_opp, std::vector<size_t> basis, std::vector<size_t> index, SpMat &HopHam);
    void BuildHopHam_Fibonacci(int species, size_t count, size_t count_opp, std::vector<size_t> basis, std::vector<size_t> index, SpMat &HopHam);
    void BuildSOCHam(int species, size_t count, size_t count_opp, std::vector<size_t> basis, std::vector<size_t> index, SpMat &HopHam, double gamma, double phase, int site1, int site2);
    void Build_NNNphase_periodic(int species, size_t count, size_t count_opp, std::vector<size_t> basis, std::vector<size_t> index, SpMat &HopHam, double t2, double gamma, double phase1, double phase2, int link_num);
    void Build_NNPhase_periodic(int species, size_t count, size_t count_opp, std::vector<size_t> basis, std::vector<size_t> index, SpMat &HopHam, double phase);
    void Build_Zeeman_Onsite(SpMat &HopHam, double B);
    
    void GetPhi(double _Phi_t);
    void GetOnsite(Tnum _h);
    int AntiSym(size_t c, size_t d, size_t bas);
    void IntMatrix_Build();
    void Save_Ham();//input can be filename from main cpp
    void Total_Ham();
    void Total_Ham_WSOC();
    void Total_Ham_Zeeman();
    void OutHam();
    void CheckHermicity(SpMat &Mat);
    
};

//constants
const double hbar = 1.;
const std::complex<double> I(0.0,1.0);
//const std::complex<double> NegI(-0.0,-1.0);
const double Pi = 4*atan(1);


#endif /* Hamiltonian_h */
