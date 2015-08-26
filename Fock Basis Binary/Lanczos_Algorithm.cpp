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
//    cout << "Test Ham: \n" << Test_Ham << endl;
//    cout << "Test Lanczos: \n" << Test_Lanczos << endl;
}

void Lanczos_Diag::Set_Mat_Dim_LA(Hamiltonian& tb)//int Tot_base
{
    TriDiag = Eigen::MatrixXd::Zero(itmax, itmax);//set a max iteration dim limit
    
    Lanczos_Vec = Eigen::VectorXd::Random(tb.Tot_base);
    G_state = Eigen::VectorXd::Zero(tb.Tot_base);
    
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
            Evec_Mat = DiagMe.eigenvectors();
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
 
    cnt = it;
    Evec.resize(cnt);
    Evec = Evec_Mat.col(0);
    //cout << "We have Eigenvectors \n" << Evec_Mat << endl;
}

void Lanczos_Diag::Get_Gstate()
{

    for(int i = 0; i < cnt; i++)
    {
        G_state += K_Mat[i]*Evec.row(i);//Evec.row(i)**K_Mat[i]
    }
    
    
//    cout << "This is the eigen vector: \n" << Evec <<endl;
//    cout << "K_MAT \n";
//    for(int i = 0; i < cnt; i++)
//    {
//        cout << K_Mat[i] << "\t";
//    }
    //cout << endl << "This is the ground state in Fock Basis: \n"<< G_state << endl;
}


void Lanczos_Diag::Gstate_RealSpace(Hamiltonian& ct_up, Hamiltonian& ct_dn, Hamiltonian& Nsite, const Hamiltonian& basis_up, const Hamiltonian& basis_dn)
{
    n_up.resize(Nsite.L,0.0);
    n_dn.resize(Nsite.L,0.0);
    cout << "count up is " << ct_up.count_up <<endl;
    cout << "count down is " << ct_dn.count_dn <<endl;
    
    for(size_t i = 1; i <= ct_up.count_up; i++)//have to start at 1 change all appropriatly
    {
        for(size_t j= 1; j <= ct_dn.count_dn; j++)
        {
            size_t ind = ((j-1)*ct_up.count_up)+i; //finding appropriate Fock state
            //cout << "Index: "<< ind << endl;
            
            double cf = G_state(ind-1)*G_state(ind-1);
            
            for(int n = 0; n < Nsite.L; n++)
            {
                size_t bas_up = basis_up.basis_up[i-1];
                size_t bas_dn = basis_dn.basis_down[j-1];
                
                if(MY_bittest(bas_up, n))//testing if up particle in Fock state on site n
                {
                    n_up.at(n) += cf;
                }
                if(MY_bittest(bas_dn, n))//testing if down particle in Fock state on site n
                {
                    n_dn.at(n) += cf;
                }
            }
        }
    }
    
    cout << "Here is the ground state for up spin: \n";
    for(int i = 0; i < Nsite.L; i++)
    {
        cout << n_up[i] << endl;
    }
    cout << "Here is the ground state for down spin: \n";
    for(int i = 0; i < Nsite.L; i++)
    {
        cout << n_dn[i] << endl;
    }
    
    
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