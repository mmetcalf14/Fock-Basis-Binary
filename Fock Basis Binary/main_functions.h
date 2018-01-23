//
//  main_functions.h
//  Fock Basis Binary
//
//  Created by mekena McGrew on 1/19/18.
//  Copyright © 2018 Mekena Metcalf. All rights reserved.
//

#ifndef main_functions_h
#define main_functions_h

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
//#include </usr/local/include/Eigen/Sparse>
//#include </usr/local/include/Eigen/Eigen>
//#include </usr/local/include/Eigen/Core>
#include "Hamiltonian.h"
#include "Lanczos.h"
#include "Observables.h"

using namespace std;
using namespace Eigen;
void Write_Density(ofstream &fout, vector<double> &n_up, vector<double> &n_dn, int L );
void DensityDiff(ofstream &fout, vector<double> &n_up, vector<double> &n_dn, int L );
void DensityInTime(ofstream &fout, vector<double> &n, int L, double Phi );
void Density(ofstream &fout, vector<double> &n, int L);
void EdgeDensity(ofstream &fout, vector<double> &n, double time);//can I move these functions to a different cpp file?


void Harmonic_Trap(vector<double> &HT, int L, double y);

template<typename Tnum>
void Fidelity(ofstream &output, Hamiltonian<Tnum> &h, Lanczos_Diag<Tnum> &d , double U, double J1, double J2, double Um, double Jm, int L);

template<typename Tnum>
void U_Fid(ofstream &output, Hamiltonian<Tnum> &h, Lanczos_Diag<Tnum> &d , double U, double J1, double J2, double Um, int L);

template<typename Tnum>
void U_Fid_V2(ofstream &output, Hamiltonian<Tnum> &h, Lanczos_Diag<Tnum> &d , double U, double J1, double J2, double Um, int L);

template<typename Tnum>
void Fidelity_HT(ofstream &output, Hamiltonian<Tnum> &h, Lanczos_Diag<Tnum> &d , double U, double J1, double J2, double Um, double Ym, vector<double> HT, int L);

template<typename Tnum>
void PeierlsTD(ofstream &out_up, ofstream &out_dn, Hamiltonian<Tnum> &h, Lanczos_Diag<Tnum> &d, int T_it, double dt, double tp, int L);

template<typename Tnum>
void QPump_TD(ofstream &out_up, ofstream &out_dn, Hamiltonian<Tnum> &ham, Lanczos_Diag<Tnum> &diag, const int T, const double dt, const double h_0, const double J_0, const double d_0, double U, int L);

template<typename Tnum>
void CorrTD(ofstream &output, ofstream &out, Hamiltonian<Tnum> &h, Lanczos_Diag<Tnum> &d, int cut, int T_it, double dt, int L);

void EVProbeGamma(ofstream &output1, ofstream &output2, ofstream &output3, Hamiltonian<complex<double>> &h, Lanczos_Diag<complex<double>> &d, int s1, int s2, int L, double phi1, double phi2, double link);
void EVProbeGammaVInt(ofstream &output1,  Hamiltonian<complex<double>> &h, Lanczos_Diag<complex<double>> &d,  int L, double phi1, double link);
void EVProbePhaseVInt(ofstream &output1, ofstream &output2,  Hamiltonian<complex<double>> &h, Lanczos_Diag<complex<double>> &d,ObservableFunctions<complex<double>> &O, double t1, double t2,  int L,  double link);
void HysteresisArea(ofstream &output1, ofstream &output2,  Hamiltonian<complex<double>> &h, Lanczos_Diag<complex<double>> &d,ObservableFunctions<complex<double>> &O, double U, double t1, double t2, int L,  double link);
void EVProbePhaseVNNNHop(ofstream &output1, ofstream &output2,  Hamiltonian<complex<double>> &h, Lanczos_Diag<complex<double>> &d,ObservableFunctions<complex<double>> &O, double t1, double t2, double U, int L,  double link);
void EVProbeInt(ofstream &output1, ofstream &output2, ofstream &output3, Hamiltonian<complex<double>> &h, Lanczos_Diag<complex<double>> &d, int s1, int s2, int L, double gSO, double phi1, int link);
void EVProbePhase(ofstream &output1, ofstream &output2, ofstream &output3, Hamiltonian<complex<double>> &h, Lanczos_Diag<complex<double>> &d, ObservableFunctions<complex<double>> &O, int s1, int s2, int L, double gSO, int link, double a);
void ProbePhase_wRealNNN(ofstream &output1, ofstream &output2, ofstream &output3, Hamiltonian<complex<double>> &h, Lanczos_Diag<complex<double>> &d, int s1, int s2, int L, double gSO, int link);
void CurrentVsGamma(ofstream &out1, ofstream &out2, ofstream &out3, Hamiltonian<complex<double>> &h, Lanczos_Diag<complex<double>> &d, double a, int links);
void ZeemanObservables(ofstream &output1, ofstream &output2, Hamiltonian<complex<double>> &h, Lanczos_Diag<complex<double>> &d, double gamma);
void NNNRealHop_Obs(ofstream &output1, ofstream &output2, ofstream &output3, Hamiltonian<complex<double>> &h, Lanczos_Diag<complex<double>> &d, double gamma, int L, double phi);
void TF_ForceCalc(ofstream &outForce,  Hamiltonian<complex<double>> &h, Lanczos_Diag<complex<double>> &d , int link, double Phi );

template<typename Tnum>
void LocalCurrent(Hamiltonian<Tnum> &h, Lanczos_Diag<Tnum> &d, int s1, int s2, double Phi);


double GdotG( const VectorXd &G1, const VectorXd &G2);
double GdotG( const VectorXcd &Gc1, const VectorXcd &Gc2);

#endif /* main_functions_h */
