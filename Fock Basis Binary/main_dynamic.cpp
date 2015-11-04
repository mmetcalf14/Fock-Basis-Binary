//
//  main.cpp
//  Fock Basis Binary
//
//  Created by mekena McGrew on 6/15/15.
//  Copyright (c) 2015 Mekena Metcalf. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include </usr/local/include/Eigen/Sparse>
#include "Hamiltonian.h"
#include "Lanczos.h"

using namespace std;
using namespace Eigen;

void Write_Density(ofstream &fout, vector<double> &n_up, vector<double> &n_dn, int L );

int main(int argc, const char * argv[])
{

    
    

    int Nup;
    int Ndown;
    int Nsite;
    double t_1;
    double t_2;
    double U;
    int T_f;
    double dt = 0.01;
    char output[60];
    Matrix4d Test_Ham;
    Vector4d Test_Lanczos;
    Test_Ham << 0, -1, 0, 0,
               -1,0,-1,0,
                0,-1,0,-1,
                0,0,-1,0;
    Test_Lanczos << 0.5,0.5,0.5,0.5;
    
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
    ReadFile >> t_1;
    ReadFile >> t_2;
    ReadFile >> U;
    ReadFile >> T_f;
    ReadFile >> output;
    
    cout << Nup << endl;
    cout << Ndown << endl;
    cout << Nsite << endl;
    cout << t_1 << endl;
    cout << t_2 << endl;
    cout << U << endl;
    cout << output << endl;
    
    int T_tot = T_f/dt;
    
    
    ofstream fout(output);
    assert(fout.is_open());
    
    fout.setf(ios::scientific);
    fout.precision(11);
    
    //Build basis and pass to Hamiltonian class through inheritance
    Hamiltonian<complex<double>> ham(Nsite, Nup, Ndown);
    
 
    //set hopping and interaction coefficients
    ham.Set_Const(t_1, t_2);//U=0 until |G> is found for t=0



    //set hamiltonian from triplets
    ham.HopMatrix_Build();
    

    //build interaction matrix
    ham.IntMatrix_Build();
    
    //add together all three matrices for total Ham
    ham.Total_Ham();
    
    //create object for diag class
    Lanczos_Diag<complex<double>> Diag(ham);//how to I do this constructor
    
    Diag.Lanczos_TestM(Test_Ham, Test_Lanczos);
    
    //set Lanczos vector dimensions
    //cout << "Setting LA Dim \n";
    Diag.Set_Mat_Dim_LA(ham);

    //cout << "Diagonalizing \n";
    //Diagonalization of t=0 Hamiltonian
    Diag.Diagonalize(ham);

    
    //convert |G> from Fock basis to onsite basis
    //seperate |G> states for nup and ndn
    //cout << "Getting Density\n";
    Diag.Density(ham);//before interaction turned on
    Write_Density(fout, Diag.n_up, Diag.n_dn, Nsite);

//    //Triplets removed to redo Interaction matrix after quenching
//    // and all non-zero elemenst of Total Ham and Ham_U are set to zero
//    ham.ClearTriplet();
//   
//
//    //Interactions turned on after intial ground state found
//    //if U=0 there is no need to build interaction ham until after ground state is found
//    ham.QuenchU(U);
//    ham.IntMatrix_Build();
//    ham.Total_Ham();
//    
//    //Lanczos_Diag<complex<double> > Diag_Comp(ham);
//    
//    //Time Evolve
//    Diag.TimeEvoCoeff(dt);
//   // Diag.Dynamics(ham, ham);
//    
//    int NN = T_tot/10;
//    int Nflag = 0;
//    for(int t = 0; t < T_tot; t++)
//    {
//        //cout << "iteration: "<< t << endl;
//        Diag.Dynamics(ham);
//        
//        if(Nflag == NN-1)
//        {
//           Diag.Density(ham);
//           Write_Density(fout, Diag.n_up, Diag.n_dn, Nsite);
//            Nflag = 0;
//        }
//        
//        Nflag++;
//    }
//    
    
    fout.close();
    cout << "Code is Done! \n";

    return 0;
}

void Write_Density(ofstream &fout, vector<double> &n_up, vector<double> &n_dn, int L )
{
    for(int i = 0; i < L; i++)
    {
        fout << i << " " << n_up[i] << " " << n_dn[i] << endl;
        
    }
    fout << endl;
    cout << endl;
    
    
}




