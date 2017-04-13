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

//template<typename Tnum>
//void Lanczos_Diag<Tnum>::TimeEvoCoeff(const double &_dt)
//{
//    I.real(0.0);
//    I.imag(1.0);
//    dt = _dt;
//    hbar = 1.;
//}

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
   // cout << "Have the Groundstate\n" << G_state.dot(G_state) << endl;
    //cout << "Test: " << G_state(1)-G_state(0) << endl;
    //cout <<"G_state: " << G_state.adjoint()*G_state << endl;
    //output.close();

    G_state.normalize();
    cout <<"G_state: " << G_state << endl;

}



template<typename Tnum>
void Lanczos_Diag<Tnum>::Density(const Hamiltonian<Tnum> &Ham)
{
    n_up.resize(Ham.L);
    n_dn.resize(Ham.L);
    double C = 0.;
    
    for(auto &j : n_up)
    {
        j = 0.0;
    }
    for(auto &j : n_dn)
    {
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
            //switched i and j above 08/30/16 now have correct index
            size_t bas_up = Ham.basis_up[i];
            size_t bas_dn = Ham.basis_down[j];
            //cout << "Index: "<< ind << endl;
            // cout << "Normalized? " << G_state.dot(G_state) << endl;
            complex<double> cf = conj(G_state(ind))*G_state(ind);
            sum += cf.real();

            
            for(int n = 0; n < Ham.L; n++)
            {
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
double Lanczos_Diag<Tnum>::DensityWCorr(const Hamiltonian<Tnum> &Ham, int cut)
{
    n_up.resize(Ham.L);
    n_dn.resize(Ham.L);
    double C = 0.;
    double Correlation = 0;
    
    size_t site_1;
    size_t site_2;
    if(cut == 0)
    {
        site_1 = Ham.L-1;
        site_2 = 0;
    }
    else if(cut == 1)
    {
        site_1 = 0;
        site_2 = 1;
    }
    else{
        site_1 = 0;
        site_2 = 0;
    }
    
    for(auto &j : n_up)
    {
        j = 0.0;
    }
    for(auto &j : n_dn)
    {
        j = 0.0;
    }
    
    double sum = 0;
    
    
    for(size_t i = 0; i < Ham.count_up; i++)//have to start at 1 change all appropriatly
    {
        for(size_t j= 0; j < Ham.count_dn; j++)
        {
            size_t ind = (j*Ham.count_up)+i; //finding appropriate Fock state
            //switched i and j above 08/30/16 now have correct index
            size_t bas_up = Ham.basis_up[i];
            size_t bas_dn = Ham.basis_down[j];
            //cout << "Index: "<< ind << endl;
            // cout << "Normalized? " << G_state.dot(G_state) << endl;
            complex<double> cf = conj(G_state(ind))*G_state(ind);
            sum += cf.real();
            //cout << "ind: " << ind << " coeff: " << cf << endl;
            //Correlation
            C += DensityCorrelation(bas_up, bas_dn, cf, site_1, site_2);
            //cout << "C: " << C << endl;
            
            for(int n = 0; n < Ham.L; n++)
            {
                if(MY_bittest(bas_up, n))//testing if up particle in Fock state on site n
                {

                    n_up.at(n) += cf.real();
                }
                if(MY_bittest(bas_dn, n))//testing if down particle in Fock state on site n
                {
                    n_dn.at(n) += cf.real();
                }
            }
        }
        
    }

    
    
    
    //Final Calculated Correlation
    Correlation = C - ((n_up[site_1]+n_dn[site_1])*(n_up[site_2]+n_dn[site_2])); //with spin correlation
    //Correlation = C - (n_up[LastSite]*n_up[0]); //no spin correlation

    return Correlation;
}

template<typename Tnum>
double Lanczos_Diag<Tnum>::DensityWCorr_O2(const Hamiltonian<Tnum> &Ham, int cut)
{

    double C = 0.;
    double D1 = 0.;
    double D2 = 0.;
    double Correlation = 0.;
    
    size_t site_1;
    size_t site_2;
    if(cut == 0)
    {
        site_1 = Ham.L-1;
        site_2 = 0;
    }
    else if(cut == 1)
    {
        site_1 = 0;
        site_2 = 1;
    }
    else{
        site_1 = 0;
        site_2 = 0;
    }
    
    
    for(size_t i = 0; i < Ham.count_up; i++)//have to start at 1 change all appropriatly
    {
        for(size_t j= 0; j < Ham.count_dn; j++)
        {
            size_t ind = (j*Ham.count_up)+i;
            
            size_t bas_up = Ham.basis_up[i];
            size_t bas_dn = Ham.basis_down[j];
            
            complex<double> cf = conj(G_state(ind))*G_state(ind);
           
            //Correlation
            C += DensityCorr_O2(bas_up, bas_dn, cf, site_1, site_2);
            
            //Density Site L-1
            D1 += OnsiteDensity_O2(bas_up, bas_dn, cf, site_1);
            //Density Site 0
            D2 += OnsiteDensity_O2(bas_up, bas_dn, cf, site_2);

            
        }
        
    }
    
    //Final Calculated Correlation
    Correlation = C - (D1*D2);
    
    return Correlation;
}



template<typename Tnum>
double Lanczos_Diag<Tnum>::DensityCorrelation(double bu, double bd, complex<double> cf, size_t site1, size_t site2)
{
    double Corr;
    double n1_up =0.;
    double n2_up=0.;
    double n1_dn=0.;
    double n2_dn=0.;
    
   

    if(MY_bittest(bu, site1))
    {
        n1_up = 1.;
    }
    if(MY_bittest(bu, site2))
    {
     n2_up = 1.;
    }
    if(MY_bittest(bd, site1))
    {
        n1_dn = 1.;
        
    }
    if(MY_bittest(bd, site2))
    {
        n2_dn = 1.;
    }
    
    Corr = cf.real()*(n1_up + n1_dn)*(n2_up + n2_dn);//with spin correlation
    //Corr = cf.real()*n1_up*n2_up;//no spin correlation
    
    return Corr;
}

template<typename Tnum>
double Lanczos_Diag<Tnum>::DensityCorr_O2(double bu, double bd, complex<double> cf, size_t site1, size_t site2)
{
    double Corr;
    double n1_up =0.;
    double n2_up=0.;
    double n1_dn=0.;
    double n2_dn=0.;
    
    
    
    if(MY_bittest(bu, site1))
    {
        n1_up = 1.;
    }
    if(MY_bittest(bu, site2))
    {
        n2_up = 1.;
    }
    if(MY_bittest(bd, site1))
    {
        n1_dn = 1.;
        
    }
    if(MY_bittest(bd, site2))
    {
        n2_dn = 1.;
    }
    
    Corr = cf.real()*(n1_up + n1_dn)*(n1_up + n1_dn)*(n2_up + n2_dn)*(n2_up + n2_dn);
    
    return Corr;
}

template<typename Tnum>
double Lanczos_Diag<Tnum>::OnsiteDensity_O2(double bu, double bd, complex<double> cf, size_t site)
{
    double D;
    double n_up =0.;
    double n_dn=0.;
    
    if(MY_bittest(bu, site))
    {
        n_up = 1.;
    }

    if(MY_bittest(bd, site))
    {
        n_dn = 1.;
        
    }

    
    D = cf.real()*((n_up + n_dn)*(n_up + n_dn));
    
    return D;
}

template<typename Tnum>
void Lanczos_Diag<Tnum>::SpinCorr(const Hamiltonian<Tnum> &Ham, std::ofstream &output, double t, int cut)
{
    double SC1 = 0.;
    double SC2 = 0.;
    
    size_t site_1;
    size_t site_2;
    if(cut == 0)
    {
        site_1 = Ham.L-1;
        site_2 = 0;
    }
    else if(cut == 1)
    {
        site_1 = 0;
        site_2 = 1;
    }
    else{
        site_1 = 0;
        site_2 = 0;
    }
    
    double n1_up =0.;
    double n1_dn = 0.;
    double n2_up = 0.;
    double n2_dn = 0.;
    
    for(size_t i = 0; i < Ham.count_up; i++)//have to start at 1 change all appropriatly
    {
        for(size_t j= 0; j < Ham.count_dn; j++)
        {
            size_t ind = (j*Ham.count_up)+i;
            
            size_t bas_up = Ham.basis_up[i];
            size_t bas_dn = Ham.basis_down[j];
            
            complex<double> cf = conj(G_state(ind))*G_state(ind);
            
            //SpinCorrelation
//            SC1 += Calc_SC(bas_dn, bas_up, cf, 0, Ham.L-1);//<n_1,dn*n_L,up>
//            SC2 += Calc_SC(bas_up, bas_dn, cf, 0, Ham.L-1);//<n_1,up*n_L,dn>
            
            SC1 += Calc_SameSpin(bas_up, cf, site_1, site_2);//<n_1,up*n_L,up>
            SC2 += Calc_SameSpin(bas_dn, cf, site_1, site_2);//<n_1,dn*n_L,dn>
            
            if(MY_bittest(bas_up, site_1))
            {
                n1_up += cf.real();
            }
            if(MY_bittest(bas_dn, site_1))
            {
                n1_dn += cf.real();
            }
            if(MY_bittest(bas_up, site_2))
            {
                n2_up += cf.real();
            }
            if(MY_bittest(bas_dn, site_2))
            {
                n2_dn += cf.real();
            }
    
        }
        
    }
    
    output << t << " " << SC1-(n1_up*n2_up) << " " << SC2-(n1_dn*n2_dn) << endl;
    //output << t << " " << SC1 << " " << SC2 << endl;
	cout << "Spin Corr: " << SC1-(n1_up*n2_up) << endl;

    
}

template<typename Tnum>
double Lanczos_Diag<Tnum>::Calc_SC(double b1, double b2, std::complex<double> cf, size_t site1, size_t site2)
{
    double SC;
    
    double n1 = 0;
    double n2 = 0;
    
    if(MY_bittest(b1, site1))
    {
        n1 = 1.;
    }
    if(MY_bittest(b2, site2))
    {
        n2 = 1.;
    }
    
    SC = cf.real()*n1*n2;
    
    return SC;
}

template<typename Tnum>
double Lanczos_Diag<Tnum>::Calc_SameSpin(double bs, std::complex<double> cf, size_t site1, size_t site2)
{
    double SC;
    
    double n1 = 0;
    double n2 = 0;
    
    if(MY_bittest(bs, site1))
    {
        n1 = 1.;
    }
    if(MY_bittest(bs, site2))
    {
        n2 = 1.;
    }
    
    SC = cf.real()*n1*n2;
    
    return SC;
}

template<typename Tnum>
complex<double> Lanczos_Diag<Tnum>::Expect_Cij(const Hamiltonian<Tnum> &Ham, int spinspec, size_t count, size_t count_opp,vector<size_t> basis, vector<size_t> index, size_t s1, size_t s2)
{
    
    complex<double> cf = 0.0;
    //complex<double> cfh;
    
    
    for(size_t bs = 0; bs < count; bs++)
    {
        size_t pb = basis[bs];
        
        size_t pi = index[pb];
        
        if ( MY_bittest(pb, s1) && !(MY_bittest(pb, s2)) )
        {
            size_t qb = MY_bitset(MY_bitclr(pb,s1),s2);
            size_t qi = index[qb];
            assert( qb != pb );
            
            //cout << "pb: " << pb << " qb: " << qb << endl;
            for(size_t k = 0; k < count_opp; k++)
            {
                int r,s;
                if ( spinspec == 0 ){
                    r = Ham.TotalIndex(pi, k);//(k*count) + p_ind;
                    s = Ham.TotalIndex(qi, k);//(k*count) + l_ind;
                    
                }else if ( spinspec == 1 ){
                    r = Ham.TotalIndex(k, pi);
                    s = Ham.TotalIndex(k, qi);
                    
                }else{
                    //cout << "More than 2 species fermion!!" << endl;
                }
                //cf += conj(G_state(s))*G_state(r)*Ham.Ham_Tot.coeffRef(r,s); //sum_nm <conj(c_n)*c_m*Hnm>
                //Is above the correct way to handle the Hamiltonian in Fock basis?
                //cout << "r: " << r << " s: "<< s << endl;
                cf += conj(G_state(s))*G_state(r);
                //cout << "coefficient: " << cf << endl;
            }
        }
    }
    
    //cfh = Ham.J1*exp(I*phi)*cf;
    
    //J = -2.*cfh.imag();
    
    return cf;
}

template<typename Tnum>
complex<double> Lanczos_Diag<Tnum>::Number(const Hamiltonian<Tnum> &Ham, int spinspec, size_t count, size_t count_opp,vector<size_t> basis, vector<size_t> index, size_t s)
{
    complex<double> cf;
    
    for(size_t bs = 0; bs < count; bs++)
    {
        size_t pb = basis[bs];
        
        size_t pi = index[pb];
        
        
        if (MY_bittest(pb, s))
        {
            
            for(size_t k = 0; k < count_opp; k++)
            {
                int r,s;
                if ( spinspec == 0 ){
                    r = Ham.TotalIndex(pi, k);//(k*count) + p_ind;
                    //(k*count) + l_ind;
                    
                }else if ( spinspec == 1 ){
                    r = Ham.TotalIndex(k, pi);
                    
                    
                }else{
                    //cout << "More than 2 species fermion!!" << endl;
                }

                cf += conj(G_state(r))*G_state(r);//I don't know which one is the correct one
            }
        }
    }

    
    
    return cf;
    
}

template<typename Tnum>
complex<double> Lanczos_Diag<Tnum>::NumberNumber(const Hamiltonian<Tnum> &Ham, int spinspec, size_t count, size_t count_opp, std::vector<size_t> basis, std::vector<size_t> index, size_t s, size_t q)
{
    
    complex<double> cf;
    double nq;
    double ns;
    
    for(size_t bs = 0; bs < count; bs++)
    {
        size_t pb = basis[bs];
        
        size_t pi = index[pb];
        
        
        if (MY_bittest(pb, s))
        {ns = 1.;}
        else{ns = 0.;}
        if (MY_bittest(pb, q))
        {nq = 1.;}
        else{nq = 0.;}
        
            for(size_t k = 0; k < count_opp; k++)
            {
                int r,s;
                if ( spinspec == 0 ){
                    r = Ham.TotalIndex(pi, k);//(k*count) + p_ind;
                    //(k*count) + l_ind;
                    
                }else if ( spinspec == 1 ){
                    r = Ham.TotalIndex(k, pi);
                    
                    
                }else{
                    //cout << "More than 2 species fermion!!" << endl;
                }
                
                cf += conj(G_state(r))*G_state(r)*ns*nq;
            }
        
    }

    
    return cf;
}

template<>
void Lanczos_Diag<complex<double>>::DebugDynamics(Hamiltonian<complex<double> > &ham, double dt)
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
void Lanczos_Diag<complex<double>>::Dynamics(Hamiltonian<complex<double> > &ham, double dt)//, Eigen::VectorXd GS)
{
    //cout << "Beginning Dynamics\n";

    int imax = 9;
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
    GetExponential(Eval, it, dt);

    Eigen::VectorXcd Temp_Gstate;
    Eigen::MatrixXcd work = Work_Q * Evec_Mat;
   // cout << "work:\n" << work.adjoint() << "\nOther way:\n" << Evec_Mat.adjoint()*Work_Q.adjoint() << endl;
    
    Temp_Gstate = ( work * D_Mat ) * ( work.adjoint() * G_state );//this should be
    // Temp_Gstate = Work_Q*Evec_Mat*D_Mat*Evec_Mat.adjoint()*Work_Q.adjoint()*G_state;//this should be
   // Temp_Gstate.normalize();
    //cout << "Inner Product: " << Temp_Gstate.dot(G_state) << endl;
    G_state = Temp_Gstate;//if I make temp_Gstate complex then G_state is always complex

    G_state.normalize();
    //cout << "G_state fin: \n" << G_state << endl;
    // std::cout << "DONE!" << std::endl;



}

template<typename Tnum>
void Lanczos_Diag<Tnum>::GetExponential(const Eigen::VectorXd& vec, int max_it, double dt)
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

template<>//only for spin orbit Hamiltonian with 6 sites
Eigen::VectorXd Lanczos_Diag<complex<double>>::FullDiagonalization(const Hamiltonian<std::complex<double> > &Ham)
{
    Eigen::MatrixXcd H = Eigen::MatrixXcd(Ham.Ham_Tot);
    Eigen::MatrixXcd EVMat;
    Eigen::VectorXd Ev;
    
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> FullDiag(H);
    Ev = FullDiag.eigenvalues();
    EVMat = FullDiag.eigenvectors();
    
    G_state = EVMat.col(0);//the first coloum in the lowest energy eigenstate
    G_state.normalize();
    
    //cout << "ground state: " << G_state << endl;
    
    return Ev;
}

template<typename Tnum>
void Lanczos_Diag<Tnum>::TotalCurrents(const Hamiltonian<Tnum> &Ham, size_t s1, size_t s2)
{
    complex<double> Cij_up;
    complex<double> Cij_dn;
    //Phi = 0 for SOC case
    
    cout <<"Phi: "<< Ham.Phi_t << endl;
    
    if(Ham.Phi_t == 0.0)
    {
        //cout << "UP\n";
    Cij_up = -Ham.J1*Expect_Cij(Ham, 0, Ham.count_up, Ham.count_dn, Ham.basis_up, Ham.index_up, s1, s2);
        //cout << "Down\n";
    Cij_dn = -Ham.J1*Expect_Cij(Ham, 1, Ham.count_dn, Ham.count_up, Ham.basis_down, Ham.index_dn, s1, s2);
        //cout << "In the correct loop\n";
    }
    else{
        Cij_up = -Ham.J1*exp(I*Ham.Phi_t)*Expect_Cij(Ham, 0, Ham.count_up, Ham.count_dn, Ham.basis_up, Ham.index_up, s1, s2);
        Cij_dn = -Ham.J1*exp(I*Ham.Phi_t)*Expect_Cij(Ham, 1, Ham.count_dn, Ham.count_up, Ham.basis_down, Ham.index_dn, s1, s2);
    }
    
    Jup = -2*Cij_up.imag();
    Jdn = -2*Cij_dn.imag();
    
    cout << "Jup: " << Jup << " Jdn: " << Jdn << endl;
    
    double JQ = Jup.real()+Jdn.real();
    double JS = Jup.real()-Jdn.real();
    cout << "JQ: " << JQ << " JS: " << JS << endl;
    
}

template<typename Tnum>
double Lanczos_Diag<Tnum>::SpinCurrent()
{
    return Jup.real()-Jdn.real();
}

template<typename Tnum>
double Lanczos_Diag<Tnum>::ChargeCurrent()
{
    return Jup.real()+Jdn.real();
}

template<typename Tnum>
complex<double> Lanczos_Diag<Tnum>::CurrentVariance(const Hamiltonian<Tnum> &Ham, int spec, size_t s1, size_t s2)
{

    complex<double> Jsquare;
    complex<double> Var;
    
    if(spec == 0)
    {
        Jsquare = CurrentSquare(Ham, 0, s1, s2);

        Var = Jsquare + (Jup*Jup);
    }
    else if (spec == 1)
    {
        Jsquare = CurrentSquare(Ham, 1, s1, s2);
        Var = Jsquare + (Jdn*Jdn);
    }
    else{
        cout << "Not a spin species\n";
    }
    
    
    return Var;
}

template<typename Tnum>
complex<double> Lanczos_Diag<Tnum>::CurrentSquare(const Hamiltonian<Tnum> &Ham, int spec, size_t s1, size_t s2)
{
    complex<double> Jsq;

    complex<double> N1;
    complex<double> N2;
    complex<double> N12;
    
    if(spec == 0)
    {
        N1 = Number(Ham, 0, Ham.count_up, Ham.count_dn, Ham.basis_up, Ham.index_up, s1);
        N2 = Number(Ham, 0, Ham.count_up, Ham.count_dn, Ham.basis_up, Ham.index_up, s2);
        N12 = NumberNumber(Ham, 0, Ham.count_up, Ham.count_dn, Ham.basis_up, Ham.index_up, s1, s2);
    }
    else if (spec == 1)
    {
        N1 = Number(Ham, 1, Ham.count_dn, Ham.count_up, Ham.basis_down, Ham.index_dn, s1);
        N2 = Number(Ham, 1, Ham.count_dn, Ham.count_up, Ham.basis_down, Ham.index_dn, s2);
        N12 = NumberNumber(Ham, 1, Ham.count_dn, Ham.count_up, Ham.basis_down, Ham.index_dn, s1, s2);
        
    }
    else{
        cout << "Not a spin species\n";
    }
    
    Jsq = Ham.J1*Ham.J1*(N1 + N2 - (2.*N12));
    
    return Jsq;
}

template class Lanczos_Diag<double>;
template class Lanczos_Diag<complex<double> >;
