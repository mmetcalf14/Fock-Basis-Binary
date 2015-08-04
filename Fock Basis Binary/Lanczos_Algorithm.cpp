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


void Lanczos_Diag::Set_Mat_Dim_LA(Hamiltonian& tb)//int Tot_base
{
    TriDiag.resize(tb.Tot_base, tb.Tot_base);//set a max iteration dim limit
    //K_Mat.resize(tb.Tot_base, tb.Tot_base);//vector class with Eigen inside
    Lanczos_Vec = Eigen::VectorXd::Random(tb.Tot_base);
    cout <<"Before normalization: \n"<< Lanczos_Vec <<endl << Lanczos_Vec.norm()<< endl;
    r_vec.resize(tb.Tot_base);
}

//template <typename Derived>//why did I have to put this template twice?
void Lanczos_Diag::Diagonalize(const Hamiltonian &Ham, Hamiltonian &tb)//const Eigen::SparseMatrixBase<Derived>
{
    bool Converged = false; //used to exit while loop should I set it to false here?
    
    Lanczos_Vec.normalize(); //check if the new vector is normalized...it is

    
    cout <<"After normalization: \n"<< Lanczos_Vec <<" " << Lanczos_Vec.norm()<< endl;
    K_Mat.push_back(Lanczos_Vec);
    cout << "First vector placed in K Mat " << K_Mat[0] << endl;
    
    int it = 0;
    int itmax = 200;

    do
    {
        if(it == 0)
        {
            r_vec = Ham.Ham_Tot*K_Mat[it];
//            cout << "Hamiltonian: \n:"<< Ham.Ham_Tot << endl;
//            cout << "K vector: \n "<< K_Mat[it] << endl;
//            cout << "r vector: \n "<< r_vec << endl;
            
        }
        else{
            r_vec = Ham.Ham_Tot*K_Mat[it]-(beta*K_Mat[it-1]);
            //cout << "Second Hamiltonian No Fail on " << it << "iteration \n";
            cout << "Beta is: " << beta << endl;
        }
       // cout << K_Mat[it] << " " << K_Mat[it].adjoint() << endl;
        alpha = K_Mat[it].adjoint()*r_vec;//this is correct
        
        r_vec -= (alpha*K_Mat[it]); //= r_vec changed to -=
        beta = r_vec.norm();
        TriDiag(it,it) = alpha;
        TriDiag(it+1,it)=beta;
        r_vec.normalize(); //this is the new Lanczos Vector
        K_Mat.push_back(r_vec);
        
        if( (beta < .00000000000001) || it == itmax)
           {
               Converged = true;
           }
        if(it == (tb.Tot_base-2))
        {
            cout << "Not Converged \n";
            Converged = true;
        }
        
        it++;
        
        
        
    }while(!(Converged));
 
    
    
}



//void Lanczos_Diag::Random_Vector()
//{
//
//    //Lanczos_Vec_Temp.setRandom();//do I need another vector for
//    //First Col of K_Mat
//    cout << Lanczos_Vec_Temp << endl;//is it ok if decimal point?
//
//}