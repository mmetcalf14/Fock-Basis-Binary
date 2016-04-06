//
//  Lanczos.cpp
//  ED Hubbard Hamiltonian
//
//  Created by mekena McGrew on 11/4/15.
//  Copyright © 2015 Mekena Metcalf. All rights reserved.
//

#include <stdio.h>
#include "Hamiltonian.h"
#include "Lanczos.h"
#include <fstream>

using namespace std;
#define TESTMAT
template<typename Tnum>
int Lanczos_Diag<Tnum>::itmax = 200;

template<typename Tnum>
void Lanczos_Diag<Tnum>::TimeEvoCoeff(const double &_dt)
{
    I.real(0.0);
    I.imag(1.0);
    dt = _dt;
    hbar = 1.;
}

template<typename Tnum>
void Lanczos_Diag<Tnum>::Lanczos_TestM(const Eigen::Matrix4d& _Test_Ham, const Eigen::Vector4d& _Test_Lanczos)
{
    Test_Ham = _Test_Ham;
    Test_Lanczos = _Test_Lanczos;

}
template<typename Tnum>
void Lanczos_Diag<Tnum>::Set_Mat_Dim_LA(const Hamiltonian<Tnum> &tb)//int Tot_base
{
    TriDiag = Eigen::MatrixXd::Zero(itmax, itmax);//set a max iteration dim limit
    
    Lanczos_Vec = VectorType::Random(tb.Tot_base);
    
    G_state = VectorType::Zero(tb.Tot_base);
    
    //Temp_G_state = Eigen::VectorXcd::Zero(tb.Tot_base);
    Q_Mat.resize(tb.Tot_base, 10);
    
    r_vec.resize(tb.Tot_base);
    rc_vec.resize(tb.Tot_base);
}
template<typename Tnum>
void Lanczos_Diag<Tnum>::Diagonalize(const Hamiltonian<Tnum> &Ham)//, Hamiltonian &tb)
{
    bool Converged = false; //used to exit while loop should I set it to false here?
    Eigen::MatrixXd Work_Mat;
    //Eigen::VectorXd Eval;
    Eigen::VectorXd Eval_P;
    std::vector<VectorType> K_Mat;
    
    Eigen::MatrixXd Evec_Mat;
    Eigen::VectorXd Eval;
    
    Eigen::VectorXd Evec;
    
    
    //#ifdef TESTMAT
    //Lanczos_Vec = Test_Lanczos;
    //    Eigen::Tridiagonalization<Eigen::MatrixXd> triOfHam(Test_Ham);
    //#endif
    //cout << "Lanczos not norm: " << Lanczos_Vec << endl;
    Lanczos_Vec.normalize();
    //cout << "Lanczos norm: " << Lanczos_Vec << endl;
    
    K_Mat.push_back(Lanczos_Vec);//error here
    int it = 0;
    
    
    do
    {
        if(it == 0)
        {
            
            r_vec = Ham.Ham_Tot*K_Mat[it]; //how can I use ifdef here?
            
            //#ifdef TESTMAT
            //r_vec = Test_Ham*K_Mat[it];
            //#endif
            
        }
        else
        {
            
            r_vec = Ham.Ham_Tot*K_Mat[it]-(beta*K_Mat[it-1]);
            
            //            #ifdef TESTMAT
            //          r_vec = Test_Ham*K_Mat[it]-(beta*K_Mat[it-1]);
            //            #endif
            
            //cout << "Beta is: " << beta << endl;
        }
        // cout << K_Mat[it] << " " << K_Mat[it].adjoint() << endl;
        alpha = K_Mat[it].dot( r_vec );//this is correct
        //cout << "alpha: " << alpha << endl;
        r_vec -= (alpha.real()*K_Mat[it]); //= r_vec changed to -=
        beta = r_vec.norm();
        //cout << "Beta: " << beta << endl;
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
        if(it == (Ham.Tot_base-2)|| it == itmax)
        {
            cout << "Not Converged \n";
            Converged = true;
        }
        
        
        
    }while(!(Converged));
    
    //cout <<" First-GS state: " <<Eval(1)-Eval(0)<< endl;
    //cout << "T mat: " << TriDiag << endl;
//    cout << "Eigenvalues: " << Eval << endl;
//    
//    ofstream output;
//    output.open("AlternatingChain_J1-2_J2-1_Delta0_y0.4_L10_Nup6_Ndn6_EnergyEV_MLM_032316.dat");
//    assert(output.is_open());
//    output.setf(ios::scientific);
//    output.precision(11);
//    
//    output << Eval << endl;
    
    cnt = it;
    Evec.resize(cnt);
    Evec = Evec_Mat.col(0);
    TriDiag.resize(10,10);
    //cout << "We have Eigenvectors \n" << Evec_Mat << endl;
    
    for(int i = 0; i < cnt; i++)
    {
        G_state += K_Mat[i]*Evec.row(i);//Evec.row(i)**K_Mat[i]
    }
    //cout << "Test: " << G_state(1)-G_state(0) << endl;
    //cout <<"G_state: " << G_state << endl;
    //output.close();
}



template<typename Tnum>
void Lanczos_Diag<Tnum>::Density(const Hamiltonian<Tnum> &Ham)
{
    n_up.resize(Ham.L);
    n_dn.resize(Ham.L);
    for(auto &j : n_up){
        j = 0.0;
    }
    for(auto &j : n_dn){
        j = 0.0;
    }
    //cout << "G_state:\n" << G_state << endl;//a couple values over 1 with time
    //this can't be true since the sum can only be one
    double sum = 0;
    
    
    for(size_t i = 0; i < Ham.count_up; i++)//have to start at 1 change all appropriatly
    {
        for(size_t j= 0; j < Ham.count_dn; j++)
        {
            size_t ind = (j*Ham.count_up)+i; //finding appropriate Fock state
            //cout << "Index: "<< ind << endl;
            // cout << "Normalized? " << G_state.dot(G_state) << endl;
            complex<double> cf = conj(G_state(ind))*G_state(ind);
            sum += cf.real();
            
            for(int n = 0; n < Ham.L; n++)
            {
                size_t bas_up = Ham.basis_up[i];
                size_t bas_dn = Ham.basis_down[j];
                
                if(MY_bittest(bas_up, n))//testing if up particle in Fock state on site n
                {
                    //cout << ind << " " << n << " " << bas_up << " " << cf << endl;
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
    
    cout << "Sum: " << sum << endl;
    
    cout << "Here is the ground state for up spin: \n";
    for(int i = 0; i < Ham.L; i++)
    {
        cout << n_up[i] << endl;
    }
    cout << "Here is the ground state for down spin: \n";
    for(int i = 0; i < Ham.L; i++)
    {
        cout << n_dn[i] << endl;
    }
    
    
}

template<typename Tnum>
void Lanczos_Diag<Tnum>::Dynamics(Hamiltonian<Tnum> &ham)//, Eigen::VectorXd GS)
{
    cout << "Beginning Dynamics\n";

    int imax = 9;
    int it = 0;
    Eigen::MatrixXd Work_Tri;
    Eigen::MatrixXcd Work_Q;
    Eigen::MatrixXd Evec_Mat;//this is real because tri_diag mat is real
    Eigen::VectorXd Eval;//real for same reason
    
    //G_state = Eigen::VectorXcd::Zero(tb.Tot_base);
    //D_Mat = Eigen::MatrixXcd::Zero(imax, imax);

    //DEGUG
    cout << "Assigning Gstate\n";
    //G_state.real() = Test_Lanczos;

    //G_state.real() = GS;


    Q_Mat.col(0) = G_state;//G_state input correctly



    do
    {

        if(it == 0)
        {

            rc_vec = ham.Ham_Tot*Q_Mat.col(it);
            //rc_vec = Test_Ham*Q_Mat.col(it);
        }
        else
        {

            rc_vec = ham.Ham_Tot*Q_Mat.col(it)-(beta*Q_Mat.col(it-1));
            //rc_vec = Test_Ham*Q_Mat.col(it)-(beta*Q_Mat.col(it-1));

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
    }while((it < imax) && (beta > .0000000001));

    //testing by adding another alpha and taking full matrix. Theorem may not prove true without
    //beta=0 elements
    rc_vec = ham.Ham_Tot*Q_Mat.col(it)-(beta*Q_Mat.col(it-1));
    alpha = Q_Mat.col(it).dot( rc_vec );//this shouldn't be giving complex
    TriDiag(it,it) = alpha.real();
    it++;


    Work_Tri = TriDiag.block(0,0,it,it);//do I need to include zero for beta
    Work_Q = Q_Mat.block(0,0,ham.Tot_base,it);

//    cout <<"Tri Matrix: \n" << TriDiag << endl;
    //cout <<"Block Tri Matrix: \n" << Work_Tri << endl;
    //cout << "Q_Mat:\n" << Work_Q << endl;
    //Work Tri looks good now when beta =0 and when it =imax

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> DiagMe(Work_Tri);
    Eval = DiagMe.eigenvalues(); //set Eval and Evec to real
    Evec_Mat = DiagMe.eigenvectors();
    //cout << "New Eigenvalues: " << Eval << endl;

    GetExponential(Eval, it);

    VectorType Temp_Gstate;//this could be set as real if Tnum is double
    //Eigen::VectorXcd Temp_Gstate;
    //do I want this redefined for every time iteration?


    Temp_Gstate = Work_Q*Evec_Mat*D_Mat*Evec_Mat.adjoint()*Work_Q.adjoint()*G_state;//this should be
    //complex because D_Mat and Work_Q are complex
    
    G_state = Temp_Gstate;

    
    G_state.normalize();
    



}

template<typename Tnum>
void Lanczos_Diag<Tnum>::GetExponential(const Eigen::VectorXd& vec, int max_it)
{
    //int it = max_it-1;
    Eigen::VectorXcd D(max_it);//max_it is the number of iteration done in Dynamics (it)(complex vector

    for(int i = 0; i < max_it; i++)
    {

        D(i) = exp(-1.*(I*dt*vec(i))/hbar);//vec is real but the I makes D complex
    }

    D_Mat = D.asDiagonal();//Dmat is complex and so is D

    //cout << "itmax: " << max_it << endl;
    //cout << Eval.size()  << endl;
    //cout << "D_Mat: \n" << D_Mat << endl;//these values are not changing with each iteration

    //return D_Mat;//this may cause problems since D_Mat is defined in class
}

template class Lanczos_Diag<double>;
template class Lanczos_Diag<complex<double> >;



