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
#define TESTMAT

void Lanczos_Diag::Lanczos_TestM(const Eigen::Matrix4d& _Test_Ham, const Eigen::Vector4d& _Test_Lanczos)
{
    Test_Ham = _Test_Ham;
    Test_Lanczos = _Test_Lanczos;
    cout << "Test Ham: \n" << Test_Ham << endl;
    cout << "Test Lanczos: \n" << Test_Lanczos << endl;
}

void Lanczos_Diag::Set_Mat_Dim_LA(Hamiltonian& tb)//int Tot_base
{
    TriDiag = Eigen::MatrixXd::Zero(tb.Tot_base, tb.Tot_base);//set a max iteration dim limit
    
    Lanczos_Vec = Eigen::VectorXd::Random(tb.Tot_base);
    G_state = Eigen::MatrixXd::Zero(tb.Tot_base);
    
    r_vec.resize(tb.Tot_base);
}

void Lanczos_Diag::Diagonalize(const Hamiltonian &Ham, Hamiltonian &tb)
{
    bool Converged = false; //used to exit while loop should I set it to false here?
    Eigen::MatrixXd Work_Mat;
    Eigen::VectorXd Eval;
    Eigen::VectorXd Eval_P;
    
    
//#ifdef TESTMAT
//    Lanczos_Vec = Test_Lanczos;
//    Eigen::Tridiagonalization<Eigen::MatrixXd> triOfHam(Test_Ham);
//#endif
    
    Lanczos_Vec.normalize(); //check if the new vector is normalized...it is

    
   // cout <<"After normalization: \n"<< Lanczos_Vec <<" " << Lanczos_Vec.norm()<< endl;
    K_Mat.push_back(Lanczos_Vec);
    //cout << "First vector placed in K Mat " << K_Mat[0] << endl;
    
    int it = 0;
    int itmax = 200;

    do
    {
        if(it == 0)
        {

            r_vec = Ham.Ham_Tot*K_Mat[it]; //how can I use ifdef here?

//#ifdef TESTMAT
//           r_vec = Test_Ham*K_Mat[it];
//#endif
            
        }
        else
        {
            
            r_vec = Ham.Ham_Tot*K_Mat[it]-(beta*K_Mat[it-1]);
            
//            #ifdef TESTMAT
//            r_vec = Test_Ham*K_Mat[it]-(beta*K_Mat[it-1]);
//            #endif
            
            //cout << "Beta is: " << beta << endl;
        }
       // cout << K_Mat[it] << " " << K_Mat[it].adjoint() << endl;
        alpha = K_Mat[it].adjoint().dot( r_vec );//this is correct
        
        r_vec -= (alpha*K_Mat[it]); //= r_vec changed to -=
        beta = r_vec.norm();
        TriDiag(it,it) = alpha;
        TriDiag(it+1,it)=beta; //self adjoint eigensolver only uses lower triangle
        r_vec.normalize(); //this is the new Lanczos Vector
        K_Mat.push_back(r_vec);
        
        
        
        it++;
        
        if (it > 1)
        {
            Work_Mat = TriDiag.block(0,0,it,it);
            //cout <<"Work mat: \n" << Work_Mat << endl;
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> DiagMe(Work_Mat);
            Eval = DiagMe.eigenvalues();
            Evec = DiagMe.eigenvectors();
            //cout << "Eigenvalues on " << it << " iteration /n" << Eval << endl;
            if(it > 2)
            {
                double Eval_Dif = Eval_P(0) - Eval(0);
                
                if(Eval_Dif < .0000000000001)
                {
                 cout << "Dif of Eval is zero \n";
                    Converged = true;
                }
            }
            Eval_P = Eval;
        }
        if( (beta < .00000000000001) )
        {
            cout << "Beta is zero \n";
            Converged = true;
        }
        if(it == (tb.Tot_base-2)|| it == itmax)
        {
            cout << "Not Converged \n";
            Converged = true;
        }
        
        
        
    }while(!(Converged));
 
    
    
}

void Lanczos_Diag::Get_Gstate(Hamiltonian& tb)
{
    //how to multiply k_mat by eigen vec?
    Eigen::VectorXd Temp = K_Mat[0];
    cout << Temp.row(1) << " Now matrix: "<< K_Mat[0] << endl; //this doesn't work
    
//    for(int i = 0; i < K_Mat.size(); i++)
//    {
//        G_state(i) = 0.;
//        for(int j = 0; j < tb.Tot_base; i++)
//        {
//            Eigen::VectorXd Temp = K_Mat[j];
//            G_state(i) += Temp.row(i)*Evec(j,1);//help on this algorithm with mixed classes
//        //ALSO ASK HOW TO CONSTRUCT MAKEFILE
//        }
//    }
    
}

//void Lanczos_Diag::Test_Tri() //Not the same as the Lanzcos algorithm
//{
// Eigen::Tridiagonalization<Eigen::MatrixXd> triOfHam(Test_Ham);
//    TriDiag = triOfHam.matrixT();
//    Eigen::MatrixXd Q = triOfHam.matrixQ();
//    cout << "This is the tridiagonal matrix: \n" << TriDiag << endl;
//    cout << "This is the corresponding K matrix: \n" << Q << endl;
//}



//void Lanczos_Diag::Random_Vector()
//{
//
//    //Lanczos_Vec_Temp.setRandom();//do I need another vector for
//    //First Col of K_Mat
//    cout << Lanczos_Vec_Temp << endl;//is it ok if decimal point?
//
//}