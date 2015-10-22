//
//  Lanczos.h
//  Fock Basis Binary
//
//  Created by mekena McGrew on 10/22/15.
//  Copyright © 2015 Mekena Metcalf. All rights reserved.
//

#ifndef Lanczos_h
#define Lanczos_h
#include "Hamiltonian.h"

class Lanczos_Diag //:public Hamiltonian
{
private:
    int itmax = 200;
    //typedef Eigen::SparseMatrix<double> SpMat;
    
    Eigen::MatrixXd TriDiag;
    
    Eigen::VectorXd Lanczos_Vec;
    Eigen::VectorXd Lanczos_Vec_Temp;
    Eigen::VectorXcd rc_vec;
    Eigen::VectorXd r_vec;
    
    Eigen::MatrixXcd D_Mat;
    Eigen::MatrixXcd Q_Mat;
    Eigen::VectorXcd G_state;
    //Eigen::VectorXcd Temp_G_state;
    
    //Eigen::VectorXd G_state_realspace;
    //Test matrices
    Eigen::Matrix4d Test_Ham;
    Eigen::Vector4d Test_Lanczos;
    
    //time evolution constants
    std::complex<double> I;
    double dt;//construct?
    double hbar;
    
    
    std::complex<double> alpha;//alpha can be complex
    double beta; //beta is real because a modulus is always real
    int cnt;
    
public:
    
    std::vector<double> n_up;//public so they can be used in main program to write the file
    std::vector<double> n_dn;
    
    Lanczos_Diag(const Hamiltonian){};//Program not accepting this constructor::SEE ERROR
    void TimeEvoCoeff(const double &_dt);
    //construct new,simple matrix to test algorithm and eigen values
    //and set Lanz vec to be one from analytical example
    void Lanczos_TestM(const Eigen::Matrix4d& _Test_Ham, const Eigen::Vector4d& _Test_Lanczos);
    void Set_Mat_Dim_LA(Hamiltonian& );//int Tot_base
    // void Random_Vector();
    
    // template <typename Derived>
    void Diagonalize(const Hamiltonian &Ham, Hamiltonian&);
    //why isn't it recognizing the template?
    
    
    //void Test_Tri();
    void Density(const Hamiltonian& ct_up, const Hamiltonian& ct_dn, Hamiltonian& Nsite,const Hamiltonian& basis_up,const Hamiltonian& basis_dn);
    void ResetLanczos();
    void GetExponential(const Eigen::VectorXd& vec, int max_it);
    void Dynamics(Hamiltonian &ham, Hamiltonian &tb);
    
    
    
};


#endif /* Lanczos_h */
