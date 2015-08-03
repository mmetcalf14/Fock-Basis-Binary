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



//void Lanczos_Diag::Random_Vector()
//{
//
//    //Lanczos_Vec_Temp.setRandom();//do I need another vector for
//    //First Col of K_Mat
//    cout << Lanczos_Vec_Temp << endl;//is it ok if decimal point?
//    
//}

void Lanczos_Diag::Set_Mat_Dim_LA(Hamiltonian& tb)//int Tot_base
{
    TriDiag.resize(tb.Tot_base, tb.Tot_base);//set a max iteration dim limit
    //K_Mat.resize(tb.Tot_base, tb.Tot_base);//vector class with Eigen inside
    Lanczos_Vec = Eigen::VectorXd::Random(tb.Tot_base);
    cout <<"Before normalization: \n"<< Lanczos_Vec << endl;
    r_vec.resize(tb.Tot_base);
}

//template <typename Derived>//why did I have to put this template twice?
void Lanczos_Diag::Diagonalize(const Hamiltonian &Ham)//const Eigen::SparseMatrixBase<Derived>
{
    bool Converged = false;//used to exit while loop should I set it to false here?
    
//    Lanczos_Vec = Lanczos_Vec_Temp;//should I use dot for inner prod?
    Lanczos_Vec.normalize(); //check if the new vector is normalized
    cout <<"A normalization: \n"<< Lanczos_Vec << endl;
    K_Mat.push_back(Lanczos_Vec);
    
    int it = 0;
    int itmax = 200;
    

 do
    {
        if(it == 0)
        {
            r_vec = Ham.Ham_Tot*K_Mat[it];
        }
        else{
            r_vec = Ham.Ham_Tot*K_Mat[it]-(alpha*K_Mat[it-1]);
        }
        
        alpha = K_Mat[it].adjoint()*r_vec;
        
        r_vec = r_vec-(alpha*K_Mat[it]);
        beta = r_vec.norm();
        TriDiag(it,it) = alpha;
        TriDiag(it+1,it)=beta;
        r_vec.normalize(); //this is the new Lanczos Vector
        K_Mat.push_back(r_vec);
        
        if( (beta < .00000000000001) || it == itmax)
           {
               Converged = true;
           }
        
        it++;
        
        
        
    }while(!(Converged));
 
    
    
}