//
//  Lanczos.h
//  ED Hubbard Hamiltonian
//
//  Created by mekena McGrew on 11/4/15.
//  Copyright © 2015 Mekena Metcalf. All rights reserved.
//

#ifndef Lanczos_h
#define Lanczos_h

#include "Hamiltonian.h"

template<typename Tnum>
class Lanczos_Diag //:public Hamiltonian
{
    //friend class ObservableFunctions;
    
private:
    static int itmax;
    //typedef Eigen::SparseMatrix<double> SpMat;
    typedef Eigen::Matrix<Tnum, Eigen::Dynamic, Eigen::Dynamic> MatrixType;
    typedef Eigen::Matrix<Tnum, Eigen::Dynamic, 1> VectorType;
    Eigen::MatrixXd TriDiag;
    //Eigen::MatrixXd TriDiag;
    
    VectorType Lanczos_Vec;
    VectorType Lanczos_Vec_Temp;
    Eigen::VectorXcd rc_vec;
    VectorType r_vec;
    
    Eigen::MatrixXcd D_Mat;//for it to work D_mat and Q_mat need to be Matrixtype and Tnum complex
    Eigen::MatrixXcd Q_Mat;
//    MatrixType D_Mat;
//    MatrixType Q_Mat;
    


    //Eigen::VectorXcd G_state;
    
    //Eigen::VectorXcd Temp_G_state;
    
    //Eigen::VectorXd G_state_realspace;
    //Test matrices
    Eigen::Matrix4d Test_Ham;
    Eigen::Vector4d Test_Lanczos;
    std::complex<double> Jup;
    std::complex<double> Jdn;    
    
    std::complex<double> alpha;//alpha can be complex
    double beta; //beta is real because a modulus is always real
    int cnt;
    
    
    
protected:


    
public:
    
    std::vector<double> n_up;//public so they can be used in main program to write the file
    std::vector<double> n_dn;
    
    Eigen::MatrixXcd EVMat;
    VectorType G_state;
    
    
    double Correlation;
    
    Lanczos_Diag(const Hamiltonian<Tnum>&){};
    //void TimeEvoCoeff(const double &_dt);

    void Lanczos_TestM(const Eigen::Matrix4d& _Test_Ham, const Eigen::Vector4d& _Test_Lanczos);
    void Set_Mat_Dim_LA(const Hamiltonian<Tnum> &tb);//int Tot_base

    void Diagonalize(const Hamiltonian<Tnum> &Ham);

    Eigen::VectorXd FullDiagonalization(const Hamiltonian<Tnum> &Ham);
    

    //void Test_Tri();
    void Density(const Hamiltonian<Tnum> &Ham);
    
    double DensityWCorr(const Hamiltonian<Tnum> &Ham, int cut);
    double DensityWCorr_O2(const Hamiltonian<Tnum> &Ham, int cut);
    void SpinCorr(const Hamiltonian<Tnum> &Ham, std::ofstream &output, double t, int cut);
    std::complex<double> Expect_Cij(Hamiltonian<Tnum> &Ham, int spinspec, size_t count, size_t count_opp, std::vector<size_t> basis, std::vector<size_t> index, Eigen::VectorXcd EigenState, size_t s1, size_t s2);
    std::complex<double> Number(const Hamiltonian<Tnum> &Ham, int spinspec, size_t count, size_t count_opp, std::vector<size_t> basis, std::vector<size_t> index, size_t s);
    std::complex<double> NumberNumber(const Hamiltonian<Tnum> &Ham, int spinspec, size_t count, size_t count_opp, std::vector<size_t> basis, std::vector<size_t> index, size_t s, size_t q);
    void TotalCurrents(Hamiltonian<Tnum> &Ham, size_t s1, size_t s2);
    double Current_UPspin(Hamiltonian<Tnum> &Ham, int EigenNum, size_t s1, size_t s2);
    double Current_DNspin(Hamiltonian<Tnum> &Ham, int EigenNum, size_t s1, size_t s2);
    std::complex<double> CurrentVariance(const Hamiltonian<Tnum> &Ham, int spec, size_t s1, size_t s2);
    std::complex<double> CurrentSquare(const Hamiltonian<Tnum> &Ham, int spec, size_t s1, size_t s2);
    
    double DensityCorrelation(double bu, double bd, std::complex<double> cf, size_t site1, size_t site2);
    double DensityCorr_O2(double bu, double bd, std::complex<double> cf, size_t site1, size_t site2);
    double OnsiteDensity_O2(double bu, double bd, std::complex<double> cf, size_t site);
    double Calc_SC(double b1, double b2, std::complex<double> cf, size_t site1, size_t site2);
    double Calc_SameSpin(double bs, std::complex<double> cf, size_t site1, size_t site2);
    
    double SpinCurrent();
    double ChargeCurrent();
    inline double UpCurrent(){return Jup.real();};
    inline double DownCurrent(){return Jdn.real();};
    
    double Occupation_AnyLevel(int spin, const Hamiltonian<Tnum> &Ham, Eigen::VectorXcd Evec);
    void OutputOccupation(std::ofstream &output, const Hamiltonian<Tnum> &Ham);
    
    void ResetLanczos();

    void GetExponential(const Eigen::VectorXd& vec, int max_it, double dt);
    void Dynamics(Hamiltonian<Tnum> &Ham, double dt);
    void DebugDynamics(Hamiltonian<Tnum> &ham, double dt);
    void CHECK();
    inline VectorType SendGstate(){return G_state;};
    inline Eigen::MatrixXcd SendEigenMat(){return EVMat;};
    
    
    
};




#endif /* Lanczos_h */
