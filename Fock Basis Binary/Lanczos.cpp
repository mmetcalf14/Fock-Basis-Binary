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

    //G_state = VectorType::Zero(tb.Tot_base);
    G_state = VectorType::Zero(tb.Tot_base);

    //Temp_G_state = Eigen::VectorXcd::Zero(tb.Tot_base);
    Q_Mat.resize(tb.Tot_base, 10);//change back when change imax

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
    //cout << "Lanczos norm: \n" << Lanczos_Vec << endl;

    K_Mat.push_back(Lanczos_Vec);//error here
    //cout << "Kmat: \n" << K_Mat[0] << endl;
    int it = 0;


    do
    {
        //cout << "it: " << it << endl;
        if(it == 0)
        {

            r_vec = Ham.Ham_Tot*K_Mat[it]; //how can I use ifdef here?
            //cout << "rvec worked\n";
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
        if(it == (Ham.Tot_base-2)|| it == (itmax-1))
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
    TriDiag.resize(20,20);
    //cout << "We have Eigenvectors \n" << Evec_Mat << endl;

    for(int i = 0; i < cnt; i++)
    {
        G_state += K_Mat[i]*Evec.row(i);//Evec.row(i)**K_Mat[i]
    }
    cout << "Have the Groundstate\n";
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

    //cout << "Density |G>: \n" << G_state << endl;
    //cout << "Sum: " << sum << endl;

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

template<>
void Lanczos_Diag<complex<double>>::DebugDynamics(Hamiltonian<complex<double> > &ham)
{
    Eigen::VectorXcd work_Gstate(ham.Tot_base);
    Eigen::VectorXcd Eval_w(ham.Tot_base);
    Eigen::VectorXcd Exp_Vec(ham.Tot_base);
    
    Eigen::MatrixXcd Evec_w(ham.Tot_base,ham.Tot_base);
    Eigen::MatrixXcd Exp_Mat = Eigen::MatrixXcd::Zero(ham.Tot_base,ham.Tot_base);
    
    //Diagonalize Ham
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> DiagMe(ham.Ham_Tot);
    Eval_w = DiagMe.eigenvalues();
    Evec_w = DiagMe.eigenvectors();
    
    //Get Time Evolution operator
    for(int i = 0; i < ham.Tot_base; i++)
    {
        Exp_Vec(i) = exp(-1.*(I*dt*Eval_w(i))/hbar);
        
    }
    Exp_Mat = Exp_Vec.asDiagonal();
    
    //Get time-evolved Gstate
    work_Gstate = (Evec_w*Exp_Mat)*(Evec_w.adjoint()*G_state);
    //work_Gstate = Eval_w;
    
    work_Gstate.normalize();
//        cout << "Hamiltonian: \n" << ham.Ham_Tot << endl;
//    cout << "Exact time-evolved Gstate: \n" << work_Gstate << endl;
//    cout << "WIth the Evolved Desnity: \n";

    G_state = work_Gstate;
    //Density(ham);
    
}

template<>
void Lanczos_Diag<complex<double>>::Dynamics(Hamiltonian<complex<double> > &ham)//, Eigen::VectorXd GS)
{
    //cout << "Beginning Dynamics\n";

    int imax = 19;
    int it = 0;
    
    TriDiag = Eigen::MatrixXd::Zero(imax+1, imax+1);
    Eigen::MatrixXd Work_Tri;
    Eigen::MatrixXcd Work_Q;
    Eigen::MatrixXcd QMat(ham.Tot_base, imax+1);
    Eigen::MatrixXd Evec_Mat(imax+1,imax+1);//this is real because tri_diag mat is real
    Eigen::VectorXd Eval(imax+1,1);//real for same reason

//    DebugDynamics(ham);
    //cout << "G_state init: \n" << G_state << endl;
    QMat.col(0) = G_state;//G_state input correctly
    
    //cout<< "ham in dynam: " <<ham.Ham_Tot;
    
    //cout << "Begin Loop\n"

    do
    {

        if(it == 0)
        {
            rc_vec = ham.Ham_Tot * QMat.col(it);
            //rc_vec = Test_Ham*Q_Mat.col(it);
        }
        else
        {
            rc_vec = ham.Ham_Tot * QMat.col(it)-( beta * QMat.col(it-1) );
            //rc_vec = Test_Ham*Q_Mat.col(it)-(beta*Q_Mat.col(it-1));
        }

//        cout << "R vec:\n" << rc_vec << endl;
//        cout << "Q vec:\n" << QMat.col(it) << endl;

        alpha = QMat.col(it).dot( rc_vec );//this shouldn't be giving complex
        //cout << "alpha: "<<alpha << endl;
        TriDiag(it,it) = alpha.real();

        rc_vec -= (alpha * QMat.col(it));

        beta = rc_vec.norm();//beta converges to zero but the iteration keeps going
        //cout << "beta: "<<beta << endl;

        TriDiag(it+1,it) = beta; //self adjoint eigensolver only uses lower triangle

        TriDiag(it,it+1) = beta;
        rc_vec.normalize(); //this is the new Lanczos Vector

        QMat.col(it+1) = rc_vec;//it+1 since we aren't using pushback

        it++;
    }while((it < (imax)) && (beta > .0000000001));

    //testing by adding another alpha and taking full matrix. Theorem may not prove true without
    //beta=0 elements
    rc_vec = ham.Ham_Tot * QMat.col(it) - ( beta * QMat.col(it-1) );
    alpha = QMat.col(it).dot( rc_vec );//this shouldn't be giving complex
    TriDiag(it,it) = alpha.real();
    //it++;
    //cout << "TriaDiag \n" << TriDiag << endl;

    Work_Tri = TriDiag.block(0,0,it,it);//do I need to include zero for beta
    //cout << "Work Tri: " << Work_Tri << endl;
    Work_Q = QMat.block(0,0,ham.Tot_base,it);
    //cout << "Diagonalize\n";
    // std::cout << "DiagMe" << std::endl;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> DiagMe(Work_Tri);
    Eval = DiagMe.eigenvalues(); //set Eval and Evec to real
    Evec_Mat = DiagMe.eigenvectors();

    //cout << "Exponential\n";
    GetExponential(Eval, it);

    Eigen::VectorXcd Temp_Gstate;
    Eigen::MatrixXcd work = Work_Q * Evec_Mat;
   // cout << "work:\n" << work.adjoint() << "\nOther way:\n" << Evec_Mat.adjoint()*Work_Q.adjoint() << endl;
    
    Temp_Gstate = ( work * D_Mat ) * ( work.adjoint() * G_state );//this should be
    // Temp_Gstate = Work_Q*Evec_Mat*D_Mat*Evec_Mat.adjoint()*Work_Q.adjoint()*G_state;//this should be

    G_state = Temp_Gstate;//if I make temp_Gstate complex then G_state is always complex

    G_state.normalize();
    //cout << "G_state fin: \n" << G_state << endl;
    // std::cout << "DONE!" << std::endl;



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

}

template<typename Tnum>
void Lanczos_Diag<Tnum>::CHECK()
{
   
    //MY NUMERICAL VS ANALYTICAL CHECK YEILDED EXACT RESULTS
    Eigen::Matrix4cd M;
    M << 0, -1, 0, 0,
    -1,0,-1,0,
    0,-1,0,-1,
    0,0,-1,0;
    
    Eigen::Vector4cd V;
    V << 0.5,0.5,0.5,0.5;
    
    int imax = 6;
    int it = 0;
    Eigen::MatrixXd Work_Tri;
    Eigen::MatrixXcd Work_Q;
    Eigen::MatrixXcd QMat(4, imax+1);
    Eigen::MatrixXd Evec_Mat(imax+1,imax+1);
    Eigen::VectorXd Eval(imax+1,1);
    
    TriDiag.resize(imax+1,imax+1);
    
    QMat.col(0) = V;
    
    do
    {
        
        if(it == 0)
        {
            rc_vec = M * QMat.col(it);
            //rc_vec = Test_Ham*Q_Mat.col(it);
        }
        else
        {
            rc_vec = M * QMat.col(it)-( beta * QMat.col(it-1) );
            //rc_vec = Test_Ham*Q_Mat.col(it)-(beta*Q_Mat.col(it-1));
        }
        
        alpha = QMat.col(it).dot( rc_vec );//this shouldn't be giving complex
        //cout << "alpha: "<<alpha << endl;
        TriDiag(it,it) = alpha.real();
        
        rc_vec -= (alpha * QMat.col(it));
        
        beta = rc_vec.norm();//beta converges to zero but the iteration keeps going
        //cout << "beta: "<<beta << endl;
        
        TriDiag(it+1,it) = beta; //self adjoint eigensolver only uses lower triangle
        
        TriDiag(it,it+1) = beta;
        rc_vec.normalize(); //this is the new Lanczos Vector
        
        QMat.col(it+1) = rc_vec;//it+1 since we aren't using pushback
        
        it++;
    }while((it < (imax)) && (beta > .0000000001));
    
    //testing by adding another alpha and taking full matrix. Theorem may not prove true without
    //beta=0 elements
    rc_vec = M * QMat.col(it) - ( beta * QMat.col(it-1) );
    alpha = QMat.col(it).dot( rc_vec );//this shouldn't be giving complex
    TriDiag(it,it) = alpha.real();
    //it++;
    cout << "TriaDiag \n" << TriDiag << endl;
    Work_Tri = TriDiag.block(0,0,it,it);
    cout << "Work-Tri \n" << Work_Tri << endl;
    
    Work_Q = QMat.block(0,0,4,it);
    //cout << "Diagonalize\n";
    // std::cout << "DiagMe" << std::endl;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> DiagMe(Work_Tri);
    Eval = DiagMe.eigenvalues(); //set Eval and Evec to real
    Evec_Mat = DiagMe.eigenvectors();
    cout << "Eigenvalues: \n" << Eval << endl;
    
    

    
    
}

template class Lanczos_Diag<double>;
template class Lanczos_Diag<complex<double> >;
