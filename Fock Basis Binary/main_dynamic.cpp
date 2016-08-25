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
#include "FunctionRef-Evaluation.h"

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
    double Umax = 20;
    double Jmax = 2;
    int T_f;
    double dt = .01;
    double t_p;

    const double h0 = 0.5;
    const double d0 = 0.5;
    const double J0 = 1.0;




    char output[80];
    char HTFout[80];
    char IntFidOut[80];
    char PD[80];
    char PeierlsDensity_up[100];
    char PeierlsDensity_down[100];
    char QPD_up[100];
    char QPD_dn[100];
    

    Matrix4d Test_Ham;
    Vector4d Test_Lanczos;
    Test_Ham << 0, -1, 0, 0,
    -1,0,-1,0,
    0,-1,0,-1,
    0,0,-1,0;
    Test_Lanczos << 0.5,0.5,0.5,0.5;

    vector<double> Harm_Trap (Nsite, 0.0);

    ifstream ReadFile("ED_J1J2_DataInput.cfg");
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
    ReadFile >> T_f;
    ReadFile >> t_p;
    ReadFile >> output;
    ReadFile >> HTFout;
    ReadFile >> IntFidOut;
    ReadFile >> PD;
    ReadFile >> PeierlsDensity_up;
    ReadFile >> PeierlsDensity_down;
    

    cout << Nup << endl;
    cout << Ndown << endl;
    cout << Nsite << endl;
    cout << t_1 << endl;
    //cout << t_2 << endl;
    cout << U << endl;
    cout << output << endl;

    

    int T_tot = T_f/dt;
    t_2 = t_1/hop_rat;
    cout << "Hopping rat: " << t_2 << endl;

    ofstream fout(output);
    assert(fout.is_open());
    fout.setf(ios::scientific);
    fout.precision(11);

    ofstream FidOut(IntFidOut);
    assert(FidOut.is_open());
    FidOut.setf(ios::scientific);
    FidOut.precision(11);

    ofstream HTOut(HTFout);
    assert(HTOut.is_open());
    HTOut.setf(ios::scientific);
    HTOut.precision(11);

    ofstream PDout(PD);
    assert(PDout.is_open());
    PDout.setf(ios::scientific);
    PDout.precision(11);

    
//    ofstream QPout_up(QPD_up);
//    assert(QPout_up.is_open());
//    QPout_up.setf(ios::scientific);
//    QPout_up.precision(11);
//    
//    ofstream QPout_dn(QPD_dn);
//    assert(QPout_dn.is_open());
//    QPout_dn.setf(ios::scientific);
//    QPout_dn.precision(11);
    

    

    ofstream PDout_up(PeierlsDensity_up);
    assert(PDout_up.is_open());
    PDout_up.setf(ios::scientific);
    PDout_up.precision(11);

    ofstream PDout_dn(PeierlsDensity_down);
    assert(PDout_dn.is_open());
    PDout_dn.setf(ios::scientific);
    PDout_dn.precision(11);

   //creating vector with harmoinc trap values per site
    //Harmonic_Trap(Harm_Trap, Nsite, Y);

    //Build basis and pass to Hamiltonian class through inheritance
    Hamiltonian<complex<double> > ham(Nsite, Nup, Ndown);
    Lanczos_Diag<complex<double> > Diag(ham);


    //ham.GetHarmTrap(Harm_Trap);
    //Fidelity_HT(HTOut, ham, Diag, U, t_1, t_2, Umax, Ymax, Harm_Trap, Nsite);
    //Fidelity(FidOut, ham, Diag, U, t_1, t_2, Umax, Jmax, Nsite);
    //U_Fid(FidOut, ham, Diag, U, t_1, t_2, Umax, Nsite);

    //set hopping and interaction coefficients
    ham.Set_Const(t_1, t_2, U);//U=0 until |G> is found for t=0



    //set hamiltonian from triplets
//    ham.HopMatrix_Build();
//
//
//    //build interaction matrix
    ham.IntMatrix_Build();
//
//    //add together all three matrices for total Ham
//    ham.Total_Ham();
//
//    //create object for diag class
//
//    //Diag.Lanczos_TestM(Test_Ham, Test_Lanczos);
//
//    //set Lanczos vector dimensions
//    cout << "Setting LA Dim \n";
//    Diag.Set_Mat_Dim_LA(ham);
//
//    cout << "Diagonalizing \n";
//    //Diagonalization of t=0 Hamiltonian
//    Diag.Diagonalize(ham);
    //Diag.CHECK();
//
//
//    //convert |G> from Fock basis to onsite basis
//    //seperate |G> states for nup and ndn
//    cout << "Getting Density\n";
//    Diag.Density(ham);//before interaction turned on
//    Write_Density(fout, Diag.n_up, Diag.n_dn, Nsite);

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
    
    PeierlsTD(PDout_up, PDout, ham, Diag, T_tot, dt, t_p, Nsite);
    

        //int NN = T_tot/10;

    FidOut.close();
    HTOut.close();
    fout.close();
    PDout.close();
    PDout_dn.close();
    PDout_up.close();

//    QPout_up.close();
//    QPout_dn.close();
    

    cout << "Code is Done! \n";

    return 0;
}

