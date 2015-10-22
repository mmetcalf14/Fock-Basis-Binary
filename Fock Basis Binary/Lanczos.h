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

template<typename Tnum>
class Lanczos_Diag //:public Hamiltonian
{
private:
    int itmax = 200;
    //typedef Eigen::SparseMatrix<double> SpMat;
    typedef Eigen::Matrix<Tnum, Eigen::Dynamic, Eigen::Dynamic> MatrixType;
    typedef Eigen::Matrix<Tnum, Eigen::Dynamic, 1> VectorType;
    MatrixType TriDiag;
    //Eigen::MatrixXd TriDiag;
    
    VectorType Lanczos_Vec;
    VectorType Lanczos_Vec_Temp;
    VectorType rc_vec;
    VectorType r_vec;
    
    MatrixType D_Mat;
    MatrixType Q_Mat;
    VectorType G_state;
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
    void Lanczos_TestM(const MatrixType& _Test_Ham, const VectorType& _Test_Lanczos);
    void Set_Mat_Dim_LA(Hamiltonian& );//int Tot_base
    // void Random_Vector();
    
    // template <typename Derived>
    void Diagonalize(const Hamiltonian &Ham);//, Hamiltonian&);
    //why isn't it recognizing the template?
    
    
    //void Test_Tri();
    void Density(const Hamiltonian& Ham);
    void ResetLanczos();
    void GetExponential(const Eigen::VectorXd& vec, int max_it);
    void Dynamics(Hamiltonian &ham);
    
    
    
};


#endif /* Lanczos_h */
