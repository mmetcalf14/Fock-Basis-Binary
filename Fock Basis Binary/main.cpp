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
#include "/usr/local/include/c++/4.9.2/Eigen/Eigen"
#include "/usr/local/include/c++/4.9.2/Eigen/Dense"
#include "/usr/local/include/c++/4.9.2/Eigen/Eigenvalues"
#include "/usr/local/include/c++/4.9.2/Eigen/Sparse"
#include "/usr/local/include/c++/4.9.2/Eigen/StdVector"
#include "Hamiltonian_Template.h"


int main(int argc, const char * argv[])
{
    using namespace std;
    using namespace Eigen;
    
    

    int Nup = 2;
    int Ndown = 1;
    int Nsite = 3;
    double tbar = 1.0;
    
    
    
    //Build basis and pass to Hamiltonian class through inheritance
    Hamiltonian ham(Nsite, Nup, Ndown);
    
    //set tbar
    ham.Set_tbar(tbar);
    ham.Set_Mat_Dim();
    ham.BuildHopHam_up();
    
    

    

    return 0;
}





