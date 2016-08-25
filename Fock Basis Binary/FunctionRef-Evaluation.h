//
//  FunctionRef-Evaluation.h
//  Fock Basis Binary
//
//  Created by mekena McGrew on 8/25/16.
//  Copyright © 2016 Mekena Metcalf. All rights reserved.
//

#ifndef FunctionRef_Evaluation_h
#define FunctionRef_Evaluation_h

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "Hamiltonian.h"
#include "Lanczos.h"
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
void Fidelity_HT(ofstream &output, Hamiltonian<Tnum> &h, Lanczos_Diag<Tnum> &d , double U, double J1, double J2, double Um, double Ym, vector<double> HT, int L);

template<typename Tnum>
void PeierlsTD(ofstream &out_up, ofstream &out_dn, Hamiltonian<Tnum> &h, Lanczos_Diag<Tnum> &d, int T_it, double dt, double tp, int L);

template<typename Tnum>
void QPump_TD(ofstream &out_up, ofstream &out_dn, Hamiltonian<Tnum> &ham, Lanczos_Diag<Tnum> &diag, const int T, const double dt, const double h_0, const double J_0, const double d_0, double U, int L);

double GdotG( const VectorXd &G1, const VectorXd &G2);
double GdotG( const VectorXcd &Gc1, const VectorXcd &Gc2);

#endif /* FunctionRef_Evaluation_h */
