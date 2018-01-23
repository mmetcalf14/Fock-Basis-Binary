//
//  Observables.h
//  Fock Basis Binary
//
//  Created by mekena McGrew on 10/17/17.
//  Copyright © 2017 Mekena Metcalf. All rights reserved.
//

#ifndef Observables_h
#define Observables_h


#include "Lanczos.h"
#include "Hamiltonian.h"


template <typename Tnum>
class ObservableFunctions //: public Lanczos_Diag<Tnum>
{
    
private:
    typedef Eigen::Matrix<Tnum, Eigen::Dynamic, Eigen::Dynamic> MatrixType;
    typedef Eigen::Matrix<Tnum, Eigen::Dynamic, 1> VectorType;
    
    Eigen::MatrixXcd EVMat;
    VectorType G_state;
    
    std::complex<double> Jup;
    std::complex<double> Jdn;

    
public:
    
    std::vector<double> n_up;
    std::vector<double> n_dn;
    
    ObservableFunctions(Lanczos_Diag<Tnum> &D, Hamiltonian<Tnum> &H) {
        Jup = 0.0;
        Jdn = 0.0;
        
    };
    void Density(const Hamiltonian<Tnum> &Ham, Lanczos_Diag<Tnum> &D);
    inline void Initialize(Lanczos_Diag<Tnum> &D){
        G_state = D.SendGstate();
        EVMat = D.SendEigenMat();}//this is a silly way to do this for large system size. Just ref from LD
    
    double DensityWCorr(const Hamiltonian<Tnum> &Ham, Lanczos_Diag<Tnum> &D, int cut);
    double DensityWCorr_O2(const Hamiltonian<Tnum> &Ham, Lanczos_Diag<Tnum> &D, int cut);
    void SpinCorr(const Hamiltonian<Tnum> &Ham, Lanczos_Diag<Tnum> &D, std::ofstream &output, double t, int cut);
    std::complex<long double> Expect_Cij(Hamiltonian<Tnum> &Ham, Lanczos_Diag<Tnum> &D, int spinspec, int TwoTunnel, size_t count, size_t count_opp, std::vector<size_t> basis, std::vector<size_t> index, Eigen::VectorXcd EigenState, size_t s1, size_t s2);
    std::complex<long double> Number(const Hamiltonian<Tnum> &Ham, Lanczos_Diag<Tnum> &D, int spinspec, size_t count, size_t count_opp, std::vector<size_t> basis, std::vector<size_t> index, Eigen::VectorXcd EigenState, size_t s);
    std::complex<long double> NumberNumber(const Hamiltonian<Tnum> &Ham, Lanczos_Diag<Tnum> &D, int spinspec, size_t count, size_t count_opp, std::vector<size_t> basis, std::vector<size_t> index, Eigen::VectorXcd EigenState, size_t i, size_t j);
    void TotalCurrents(Hamiltonian<Tnum> &Ham, Lanczos_Diag<Tnum> &D, size_t s1, size_t s2);
    double Current_UPspin(Hamiltonian<Tnum> &Ham, Lanczos_Diag<Tnum> &D, int TwoTunnel, int EigenNum, double phi, size_t s1, size_t s2);
    double Current_DNspin(Hamiltonian<Tnum> &Ham, Lanczos_Diag<Tnum> &D, int TwoTunnel, int EigenNum, double phi, size_t s1, size_t s2);
    std::complex<double> CurrentVariance(Hamiltonian<Tnum> &Ham, Lanczos_Diag<Tnum> &D, int EigenNum, int spec, size_t s1, size_t s2);
    std::complex<long double> CurrentSquare(Hamiltonian<Tnum> &Ham, Lanczos_Diag<Tnum> &D, int EigenNum, int spec, size_t s1, size_t s2);
    
    double DensityCorrelation(double bu, double bd, std::complex<double> cf, size_t site1, size_t site2);
    
    double DensityCorr_O2(double bu, double bd, std::complex<double> cf, size_t site1, size_t site2);
    double OnsiteDensity_O2(double bu, double bd, std::complex<double> cf, size_t site);
    double Calc_SC(double b1, double b2, std::complex<double> cf, size_t site1, size_t site2);
    double Calc_SameSpin(double bs, std::complex<double> cf, size_t site1, size_t site2);
    
    double SpinCurrent();
    double ChargeCurrent();
    inline double UpCurrent(){return Jup.real();};
    inline double DownCurrent(){return Jdn.real();};
    
    inline double Area(double dx, double y){return dx*y;};
    inline double Trapezoidal(double dx, double y1, double y2){return dx*((y1+y2)/2.);}
    
    double Occupation_AnyLevel(int spin, const Hamiltonian<Tnum> &Ham, Eigen::VectorXcd Evec);
    void OutputOccupation(std::ofstream &output, const Hamiltonian<Tnum> &Ham, Lanczos_Diag<Tnum> &D);
    
};

#endif /* Observables_h */
