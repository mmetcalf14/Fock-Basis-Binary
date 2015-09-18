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

void Lanczos_Diag::TimeEvoCoeff(const double &_dt)
{
    I.real(0.0);
    I.imag(1.0);
    cout << "I: " << I << endl;
    dt = _dt;
    hbar = 1.;
}

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

    G_state = Eigen::VectorXcd::Zero(tb.Tot_base);

    //Temp_G_state = Eigen::VectorXcd::Zero(tb.Tot_base);
    Q_Mat.resize(tb.Tot_base, 10);

    r_vec.resize(tb.Tot_base);
    rc_vec.resize(tb.Tot_base);
}

void Lanczos_Diag::Diagonalize(const Hamiltonian &Ham, Hamiltonian &tb)
{
    bool Converged = false; //used to exit while loop should I set it to false here?
    Eigen::MatrixXd Work_Mat;
    //Eigen::VectorXd Eval;
    Eigen::VectorXd Eval_P;
    std::vector<Eigen::VectorXd> K_Mat;

    Eigen::MatrixXd Evec_Mat;
    Eigen::VectorXd Eval;

    Eigen::VectorXd Evec;


//#ifdef TESTMAT
//    Lanczos_Vec = Test_Lanczos;
//    Eigen::Tridiagonalization<Eigen::MatrixXd> triOfHam(Test_Ham);
//#endif

    Lanczos_Vec.normalize();

    K_Mat.push_back(Lanczos_Vec);//error here
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
        alpha = K_Mat[it].dot( r_vec );//this is correct

        r_vec -= (alpha.real()*K_Mat[it]); //= r_vec changed to -=
        beta = r_vec.norm();
        TriDiag(it,it) = alpha.real();
        TriDiag(it+1,it)=beta; //self adjoint eigensolver only uses lower triangle
        TriDiag(it,it+1)=beta;
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
    TriDiag.resize(10,10);
    //cout << "We have Eigenvectors \n" << Evec_Mat << endl;

    for(int i = 0; i < cnt; i++)
    {
        G_state.real() += K_Mat[i]*Evec.row(i);//Evec.row(i)**K_Mat[i]
    }
}




void Lanczos_Diag::Density(const Hamiltonian& ct_up, const Hamiltonian& ct_dn, Hamiltonian& Nsite, const Hamiltonian& basis_up, const Hamiltonian& basis_dn)
{
    n_up.resize(Nsite.L,0.0);
    n_dn.resize(Nsite.L,0.0);

    // cout << "G_state:\n" << G_state << endl;//a couple values over 1 with time
                                            //this can't be true since the sum can only be one

    for(size_t i = 1; i <= ct_up.count_up; i++)//have to start at 1 change all appropriatly
    {
        for(size_t j= 1; j <= ct_dn.count_dn; j++)
        {
            size_t ind = ((j-1)*ct_up.count_up)+i; //finding appropriate Fock state
            //cout << "Index: "<< ind << endl;

            complex<double> cf = conj(G_state(ind-1))*G_state(ind-1);


            for(int n = 0; n < Nsite.L; n++)
            {
                size_t bas_up = basis_up.basis_up[i-1];
                size_t bas_dn = basis_dn.basis_down[j-1];

                if(MY_bittest(bas_up, n))//testing if up particle in Fock state on site n
                {
                    //cout << ind << " " << n << " " << bas_up << endl;
                    //the above are consistent with each iteration
                    //no large values for G_state
                    //cout << "cf: " << cf << endl; cf is fully real
                    n_up.at(n) += cf.real();
                }
                if(MY_bittest(bas_dn, n))//testing if down particle in Fock state on site n
                {
                    n_dn.at(n) += cf.real();
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


void Lanczos_Diag::Dynamics(Hamiltonian &ham, Hamiltonian &tb)
{
    //cout << "Beginning Dynamics\n";

    int imax = 9;
    int it = 0;
    Eigen::MatrixXd Work_Tri;
    Eigen::MatrixXcd Work_Q;
    Eigen::MatrixXd Evec_Mat;
    Eigen::VectorXd Eval;
    //D_Mat = Eigen::MatrixXcd::Zero(imax, imax);

    //cout << "G_state before norm\n" << G_state << endl;
    cout << G_state(0) << endl;
    G_state.normalize();
    //cout << "G_state after norm\n" << G_state << endl;
    Q_Mat.col(0) = G_state;//G_state input correctly
    cout << G_state(0) << endl;


    TriDiag.setZero();
    do
    {
        if(it == 0)
        {

            rc_vec = ham.Ham_Tot*Q_Mat.col(it);

        }
        else
        {

            rc_vec = ham.Ham_Tot*Q_Mat.col(it)-(beta*Q_Mat.col(it-1));

        }

        alpha = Q_Mat.col(it).dot( rc_vec );//this shouldn't be giving complex
        //cout << "alpha: "<<alpha << endl;
        TriDiag(it,it) = alpha.real();

        rc_vec -= (alpha*Q_Mat.col(it));

        beta = rc_vec.norm();//beta converges to zero but the iteration keeps going
        //cout << "beta: "<<beta << endl;

        TriDiag(it+1,it)=beta; //self adjoint eigensolver only uses lower triangle

        TriDiag(it,it+1)=beta;
        rc_vec.normalize(); //this is the new Lanczos Vector

        Q_Mat.col(it+1) = rc_vec;//it+1 since we aren't using pushback

        it++;
    }while((it < imax) && (beta > 1.0e-6));

    //testing by adding another alpha and taking full matrix. Theorem may not prove true without
    //beta=0 elements
    rc_vec = ham.Ham_Tot*Q_Mat.col(it)-(beta*Q_Mat.col(it-1));
    alpha = Q_Mat.col(it).dot( rc_vec );//this shouldn't be giving complex
    TriDiag(it,it) = alpha.real();
    it++;


    Work_Tri = TriDiag.block(0,0,it,it);//do I need to include zero for beta
    Work_Q = Q_Mat.block(0,0,tb.Tot_base,it);

//    cout <<"Tri Matrix: \n" << TriDiag << endl;
    cout <<"Block Tri Matrix: \n" << Work_Tri << endl;
    //cout << "Q_Mat:\n" << Work_Q << endl;
    //Work Tri looks good now when beta =0 and when it =imax

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> DiagMe(Work_Tri);
    Eval = DiagMe.eigenvalues(); //set Eval and Evec to real
    Evec_Mat = DiagMe.eigenvectors();
    cout << "New Eigenvalues: " << Eval << endl;

    GetExponential(Eval, it);

    Eigen::VectorXcd Temp_Gstate;
    cout << Work_Q.adjoint()*Work_Q << endl;

    Temp_Gstate = Work_Q*Evec_Mat*D_Mat*Evec_Mat.adjoint()*Work_Q.adjoint()*G_state;
    //cout << "Got the new ground state\n";
    // cout << G_state.dot(Temp_Gstate) << endl;
    cout << Temp_Gstate.dot(Temp_Gstate) << endl;
    cout << G_state(0) << " " << Temp_Gstate(0) << endl;
    G_state = Temp_Gstate;
    cin.get();


}


void Lanczos_Diag::GetExponential(const Eigen::VectorXd& vec, int max_it)
{
    //int it = max_it-1;
    Eigen::VectorXcd D(max_it);//max_it is the number of iteration done in Dynamics (it)

    for(int i = 0; i < max_it; i++)
    {

        D(i) = exp(-1.*(I*dt*vec(i))/hbar);//cos((dt*Eval(i))/hbar)
    }

    D_Mat = D.asDiagonal();

    //cout << "itmax: " << max_it << endl;
    //cout << Eval.size()  << endl;
    //cout << "D_Mat: \n" << D_Mat << endl;//these values are not changing with each iteration

    //return D_Mat;//this may cause problems since D_Mat is defined in class
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
