//
//  main_dynamic.cpp
//  ED Hubbard Hamiltonian
//
//  Created by mekena McGrew on 11/4/15.
//  Copyright © 2015 Mekena Metcalf. All rights reserved.
//


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
#include "main_functions.h"

using namespace std;
using namespace Eigen;

int main(int argc, const char * argv[])
{

    int Nup;
    int Ndown;
    int Nsite;
    double Y = 0.4;
    double Ymax = 1.;
    double t_1 = 1.;
    double t_2;
    double hop_rat;
    double U;
    double Umax = 5.;
    double Jmax = 2.;
    int T_f;
    double dt = .01;
    double t_p;
    int links;

    const double h0 = 0.5;
    const double d0 = 0.5;
    const double J0 = 1.0;
    double cut;
    double gamma;
    double B;
    double Phi_B;
    double Phi_SOC = Pi/2.;
    double PhiNN;
    double PhiNNN = Pi/2.;//this is the NNN phi and it determine the NN phi
    double alpha;
    int Site1;
    int Site2;



    char output[80];
//    char HTFout[80];
    char IntFidOut[80];
    char output2[80];
    char NumOut[80];
//    char PD[80];
//    char PeierlsDensity_up[100];
//    char PeierlsDensity_down[100];
//    char QPD_up[100];
//    char QPD_dn[100];
//    char CorrOut[100];
    

    Matrix4d Test_Ham;
    Vector4d Test_Lanczos;
    Test_Ham << 0, -1, 0, 0,
    -1,0,-1,0,
    0,-1,0,-1,
    0,0,-1,0;
    Test_Lanczos << 0.5,0.5,0.5,0.5;

    vector<double> Harm_Trap (Nsite, 0.0);

    ifstream ReadFile("ED_ActAro_J1J2_DataInput.cfg");
    assert(ReadFile.is_open());
    if (!ReadFile.is_open())
    {
        cout<<"NOT OPEN"<<endl;
        exit (1);
    }

    ReadFile >> Nup;
    ReadFile >> Ndown;
    ReadFile >> Nsite;
    ReadFile >> hop_rat;
    ReadFile >> U;
    ReadFile >> gamma;
    ReadFile >> B;
    ReadFile >> alpha;
    ReadFile >> links;
//    ReadFile >> Site1;
//    ReadFile >> Site2;
    ReadFile >> output;
    ReadFile >> output2;
    ReadFile >> IntFidOut;
    ReadFile >> NumOut;
    

    cout << Nup << endl;
    cout << Ndown << endl;
    cout << Nsite << endl;
    cout << t_1 << endl;
    //cout << t_2 << endl;
    cout << U << endl;
    cout << output << endl;

    

    int T_tot = T_f/dt;
    if(hop_rat != 0)
    {
        t_2 = t_1/hop_rat;
    }
    else{t_2 = 0.0;}
    
    cout << "Hopping rat: " << t_2 << endl;
    
    PhiNN = (4./(6.*alpha -2.))*PhiNNN;
    Phi_B = (sqrt(3.)/4.)*B;
    
    //cout << "Phi_B: " << Phi_B << " Phi_SOC: " << Phi_SOC << endl;
    

    ofstream fout(output);
    assert(fout.is_open());
    fout.setf(ios::scientific);
    fout.precision(11);
    
    ofstream OutFile(output2);
    assert(OutFile.is_open());
    OutFile.setf(ios::scientific);
    OutFile.precision(11);
    
    ofstream FidOut(IntFidOut);
    assert(FidOut.is_open());
    FidOut.setf(ios::scientific);
    FidOut.precision(11);
    
    ofstream OccOut(NumOut);
    assert(OccOut.is_open());
    OccOut.setf(ios::scientific);
    OccOut.precision(11);
//
//    ofstream HTOut(HTFout);
//    assert(HTOut.is_open());
//    HTOut.setf(ios::scientific);
//    HTOut.precision(11);
//
//    ofstream PDout(PD);
//    assert(PDout.is_open());
//    PDout.setf(ios::scientific);
//    PDout.precision(11);

    
//    ofstream QPout_up(QPD_up);
//    assert(QPout_up.is_open());
//    QPout_up.setf(ios::scientific);
//    QPout_up.precision(11);
//    
//    ofstream QPout_dn(QPD_dn);
//    assert(QPout_dn.is_open());
//    QPout_dn.setf(ios::scientific);
//    QPout_dn.precision(11);
    

    

//    ofstream PDout_up(PeierlsDensity_up);
//    assert(PDout_up.is_open());
//    PDout_up.setf(ios::scientific);
//    PDout_up.precision(11);
//
//    ofstream PDout_dn(PeierlsDensity_down);
//    assert(PDout_dn.is_open());
//    PDout_dn.setf(ios::scientific);
//    PDout_dn.precision(11);
//    
//    ofstream OutCorr(CorrOut);
//    assert(OutCorr.is_open());
//    OutCorr.setf(ios::scientific);
//    OutCorr.precision(11);

   //creating vector with harmoinc trap values per site
    //Harmonic_Trap(Harm_Trap, Nsite, Y);

    //Build basis and pass to Hamiltonian class through inheritance
    //COMPLEX
    Hamiltonian<complex<double>> ham(Nsite, Nup, Ndown);
    Lanczos_Diag<complex<double>> Diag(ham);
    ObservableFunctions<complex<double>> OF(Diag, ham);
    
    
    //DOUBLE
    //Hamiltonian<double> ham(Nsite, Nup, Ndown);
//    Lanczos_Diag<double> Diag(ham);


    //ham.GetHarmTrap(Harm_Trap);
    //Fidelity_HT(HTOut, ham, Diag, U, t_1, t_2, Umax, Ymax, Harm_Trap, Nsite);
    //Fidelity(FidOut, ham, Diag, U, t_1, t_2, Umax, Jmax, Nsite);
    
    //U_Fid_V2(FidOut, ham, Diag, U, t_1, t_2, Umax, Nsite);

    //set hopping and interaction coefficients
    ham.Set_Const(t_1, t_2, U);//U=0 until |G> is found for t=0


    //set hamiltonian from triplets
    //ham.HopMatrix_Build();
    //ham.HopMatrix_Build_Periodic();
    //ham.HopMatrix_Build_PeriodicWithSOC(Site1, Site2, gamma, Phi);
    //ham.HopMatrix_Build_PeriodicNNNHop(links, gamma, Phi_SOC);
    //ham.HopMatrix_Build_HaldHam_NoGauge(links, gamma, PhiNN, PhiNNN);
    //ham.HopMatrix_Build_Zeeman_WSOC(6, gamma, Phi_SOC, Phi_B, B);
    
//    long double Phi_Difficult = (5./6.-.05)*Pi;
//    ham.HopMatrix_Build_PeriodicBfield(links, Phi_Difficult);//Benzene with uniform flux
    
    //ham.GetPhi(Phi);
    //ham.HopMatrix_Build_Peierls();

    
    //build interaction matrix
    //ham.IntMatrix_Build();
//
//    //add together all three matrices for total Ham
    //ham.Total_Ham();
    //ham.Total_Ham_WSOC();
    //ham.Total_Ham_Zeeman();
//
//    //create object for diag class
//
//    //Diag.Lanczos_TestM(Test_Ham, Test_Lanczos);
//
    
//    //set Lanczos vector dimensions
//    cout << "Setting LA Dim \n";
    //Diag.Set_Mat_Dim_LA(ham);
//
    //cout << "Diagonalizing \n";
//    //Diagonalization of t=0 Hamiltonian
    //Diag.Diagonalize(ham);
    //Diag.CHECK();
    //VectorXd V = Diag.FullDiagonalization(ham);
//    cout << "EigenVals: " << V << endl;
//
    //Output occupation number on each energy level
    //Diag.OutputOccupation(OccOut, ham);
//
//    //convert |G> from Fock basis to onsite basis
//    //seperate |G> states for nup and ndn
    //cout << "Getting Density\n";
    //Diag.Density(ham);//before interaction turned on
    //Write_Density(fout, Diag.n_up, Diag.n_dn, Nsite);
    
    
    //Correlation Functions
    //Use Cut 0 to
    //double C;
   // C = Diag.DensityWCorr(ham, 0);
    
	//Diag.SpinCorr(ham,fout, 0, 0);//spin correlation
    
   // cout << "Correlation:  " << C << endl;
    
    
    //Current
    //Diag.TotalCurrents(ham, 0, 1);
    
//    complex<double> Ju = Diag.Current_UPspin(ham, 0, 1);
//    complex<double> Jd = Diag.Current_DNspin(ham, 0, 1);
//    cout << "Current up: " << Ju << " Current down: " << Jd << endl;
    
    //cout << "SpinCUrr: " << Diag.SpinCurrent() << " ChargeCurr: " << Diag.ChargeCurrent() << endl;
    
//    double Jup_NN_1 = OF.Current_UPspin(ham, Diag, 0, 0, Phi_Difficult, 0, 1);
//    double Jdn_NN_1 = OF.Current_UPspin(ham, Diag, 1, 0, Phi_Difficult, 0, 1);
//    complex<double> var_NN_1 = sqrt(OF.CurrentVariance(ham, Diag, 1, 0, 1));
////
//    cout << "Current NN: " << Jup_NN_1 << " Fluctuations NN: " << var_NN_1 << endl;

        //Triplets removed to redo Interaction matrix after quenching
        // and all non-zero elemenst of Total Ham and Ham_U are set to zero
        //ham.ClearInteractTriplet();


        //Interactions turned on after intial ground state found
        //if U=0 there is no need to build interaction ham until after ground state is found
//        ham.QuenchU(U);
//        ham.IntMatrix_Build();
//        ham.Total_Ham();

        //Lanczos_Diag<complex<double> > Diag_Comp(ham);

        //Time Evolve

    
    //QPump_TD(PDout_up, PDout, ham, Diag, T_f, dt, h0, J0, d0, U, Nsite);//QPout_up, QPout_dn
    
    //PeierlsTD(PDout_up, PDout, ham, Diag, T_tot, dt, t_p, Nsite);
    
    //CorrTD(OutCorr, fout, ham, Diag, cut, T_tot, dt, Nsite);
    
    //EVProbeGamma(fout, OutFile, FidOut, ham, Diag, Site1, Site2, Nsite, PhiNN, PhiNNN, links);
    //EVProbeInt(fout, OutFile, FidOut, ham, Diag, Site1, Site2, Nsite, gamma, PhiNNN, links);
    //EVProbeGammaVInt(fout, ham, Diag, Nsite, Phi_SOC, links);
    //EVProbePhase(fout, OutFile, FidOut, ham, Diag, OF, Site1, Site2, Nsite, gamma, links, alpha);
    HysteresisArea(OutFile, FidOut, ham, Diag, OF, U, t_1, t_2, Nsite, links);
    //EVProbePhaseVInt(OutFile, FidOut, ham, Diag, OF, t_1, t_2, Nsite, links);
    //ProbePhase_wRealNNN(fout, OutFile, FidOut, ham, Diag, Site1, Site2, Nsite, gamma, links);
    //CurrentVsGamma(fout, OutFile, FidOut, ham, Diag, alpha, links);
    //ZeemanObservables(fout, OutFile, ham, Diag, gamma);
    //NNNRealHop_Obs(fout, OutFile, FidOut, ham, Diag, gamma, Nsite, Phi_SOC);
    
    //TF_ForceCalc(FidOut, ham, Diag, links, Phi_SOC);
    
        //int NN = T_tot/10;

    FidOut.close();
    OutFile.close();
    OccOut.close();
//    HTOut.close();
//    PDout.close();
//    PDout_dn.close();
//    PDout_up.close();
//    OutCorr.close();

//    QPout_up.close();
//    QPout_dn.close();
    
    fout.close();
    cout << "Code is Done! \n";

    return 0;
}


template class Hamiltonian<double>;
template class Hamiltonian<complex<double>>;


