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
//#include "/usr/local/include/c++/4.9.2/Eigen/Eigen"
//#include "/usr/local/include/c++/4.9.2/Eigen/Dense"
//#include "/usr/local/include/c++/4.9.2/Eigen/Eigenvalues"
#include "/usr/include/Eigen/Sparse"
//#include "/usr/local/include/c++/4.9.2/Eigen/StdVector"
#include "Hamiltonian_Template.h"
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
    char output[60];
    Matrix4d Test_Ham;
    Vector4d Test_Lanczos;
    Test_Ham << 0, -1, 0, 0,
               -1,0,-1,0,
                0,-1,0,-1,
                0,0,-1,0;
    Test_Lanczos << 0.5,0.5,0.5,0.5;
    
    ifstream ReadFile("ED_J1J2_Static_DataInput.cfg");
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
    ReadFile >> output;
    
    cout << Nup << endl;
    cout << Ndown << endl;
    cout << Nsite << endl;
    cout << t_1 << endl;
    cout << t_2 << endl;
    cout << U << endl;
    cout << output << endl;
    
    ofstream fout(output);
    assert(fout.is_open());
    
//    ofstream EVout("ED_J1-1_J2-2_U100_Nu1_Nd2_L3_EV_101315.dat");
//    assert(fout.is_open());
    
    fout.setf(ios::scientific);
    fout.precision(11);
    
    //Build basis and pass to Hamiltonian class through inheritance
    Hamiltonian ham(Nsite, Nup, Ndown);
    
 
    //set hopping and interaction coefficients
    ham.Set_Const(t_1, t_2, U);

    //Set Dimensions for all matrices in Ham class
    ham.Set_Mat_Dim();
    
    //building seperate hopping hamiltonian for up and down spin
    ham.BuildHopHam_up();
    ham.BuildHopHam_dn();
    //set hamiltonian from triplets
    ham.HopMatrix_Build();
    
    //create indices in Fock basis
    ham.Interaction_Index();
    //build interaction matrix
    ham.BaseInteraction();
    ham.IntMatrix_Build();
    ham.Build_Interactions();
    
    //add together all three matrices for total Ham
    ham.Total_Ham();
    //ham.Check_Degeneracy(EVout);
    
    //create object for diag class
    Lanczos_Diag Diag(ham);//how to I do this constructor
    
    //Diag.Lanczos_TestM(Test_Ham, Test_Lanczos);
    
    //set Lanczos vector dimensions
    Diag.Set_Mat_Dim_LA(ham);
    //create random matrix
    //Diagonalize.Random_Vector(); done in dimension algorithm
    Diag.Diagonalize(ham, ham);
    //Diag.Test_Tri();
    Diag.Get_Gstate();
    //convert |G> from Fock basis to onsite basis
    //seperate |G> states for nup and ndn
    Diag.Gstate_RealSpace(ham, ham, ham, ham, ham);
    
    Write_Density(fout, Diag.n_up, Diag.n_dn, Nsite);
    
    fout.close();
    //EVout.close();
    cout << "Code is Done! \n";

    return 0;
}

void Write_Density(ofstream &fout, vector<double> &n_up, vector<double> &n_dn, int L )
{
    for(int i = 0; i < L; i++)
    {
        fout << i << " " << n_up[i] << " " << n_dn[i] << endl;
    }
}




