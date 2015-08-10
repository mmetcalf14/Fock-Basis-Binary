//
//  main.cpp
//  Fock Basis Binary
//
//  Created by mekena McGrew on 6/15/15.
//  Copyright (c) 2015 Mekena Metcalf. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
//#include "/usr/local/include/c++/4.9.2/Eigen/Eigen"
//#include "/usr/local/include/c++/4.9.2/Eigen/Dense"
//#include "/usr/local/include/c++/4.9.2/Eigen/Eigenvalues"
#include "/usr/include/Eigen/Sparse"
//#include "/usr/local/include/c++/4.9.2/Eigen/StdVector"
#include "Hamiltonian_Template.h"


int main(int argc, const char * argv[])
{
    using namespace std;
    using namespace Eigen;
    
    

    int Nup = 2;
    int Ndown = 3;
    int Nsite = 4;
    double tbar = 1.0;
    double U = 1.0;
    Matrix4d Test_Ham;
    Vector4d Test_Lanczos;
    Test_Ham << 0, -1, 0, 0,
               -1,0,-1,0,
                0,-1,0,-1,
                0,0,-1,0;
    Test_Lanczos << 0.5,0.5,0.5,0.5;
    
    
    //Build basis and pass to Hamiltonian class through inheritance
    Hamiltonian ham(Nsite, Nup, Ndown);
    
 
    //set tbar
    ham.Set_Const(tbar, U);

    ham.Set_Mat_Dim();
    
    ham.BuildHopHam_up();
    ham.BuildHopHam_dn();
    ham.HopMatrix_Build();
    
    ham.Interaction_Index();
    ham.BaseInteraction();
    ham.IntMatrix_Build();
    ham.Build_Interactions();
    ham.Total_Ham();
    
    Lanczos_Diag Diag(ham);//how to I do this constructor
    
    Diag.Lanczos_TestM(Test_Ham, Test_Lanczos);
    //set Lanczos vector dimensions
    Diag.Set_Mat_Dim_LA(ham);
    //create random matrix
    //Diagonalize.Random_Vector(); done in dimension algorithm
    Diag.Diagonalize(ham, ham);
    //Diag.Test_Tri();
    Diag.Get_Gstate();
    Diag.Gstate_RealSpace(ham, ham, ham, ham, ham, ham);
    
    
    cout << "Code is Done! \n";

    return 0;
}





