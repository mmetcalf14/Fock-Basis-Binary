//
//  Lanczos_Algorithm.cpp
//  Fock Basis Binary
//
//  Created by mekena McGrew on 7/28/15.
//  Copyright (c) 2015 Mekena Metcalf. All rights reserved.
//

#include <stdio.h>
#include "Hamiltonian_Template.h"

using namespace std;




void Lanczos_Diag::Random_Vector()
{
    Lanczos_Vec_Temp.setRandom();//do I need another vector for
    //First Col of K_Mat
    cout << Lanczos_Vec_Temp << endl;
    
}

void Lanczos_Diag::Set_Mat_Dim_LA(int Tot_base)
{
    TriDiag.resize(Tot_base, Tot_base);//set a max iteration dim limit
    K_Mat.resize(Tot_base, Tot_base);//vector class with Eigen inside
    Lanczos_Vec.resize(Tot_base);
    Lanczos_Vec_Temp.resize(Tot_base);
    r_vec.resize(Tot_base);
}

template <typename Derived>//why did I have to put this template twice?
void Lanczos_Diag::Diagonalize(const Eigen::SparseMatrixBase<Derived> &Ham_Tot)
{
    Lanczos_Vec = Lanczos_Vec_Temp;//should I use dot for inner prod?
    Lanczos_Vec.normalize();
    K_Mat.col(0)=Lanczos_Vec;
    
    r_vec = Ham_Tot*Lanczos_Vec;
    K_Mat.col(1) = r_vec;//second col of K_Mat
    alpha = Lanczos_Vec.adjoint()*r_vec;
    TriDiag(0,0) = alpha;
    r_vec = r_vec-(alpha*Lanczos_Vec);
    beta = r_vec.norm();
    TriDiag(1,0)=beta;
    
    
    
}