//
//  Observables.cpp
//  Fock Basis Binary
//
//  Created by mekena McGrew on 10/17/17.
//  Copyright © 2017 Mekena Metcalf. All rights reserved.
//

#include <stdio.h>
#include "Observables.h"
#include <fstream>

using namespace std;

template<typename Tnum>
void ObservableFunctions<Tnum>::Density(const Hamiltonian<Tnum> &Ham, Lanczos_Diag<Tnum> &D)
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

    double sum = 0;
    
    
    for(size_t i = 0; i < Ham.count_up; i++)
    {
        for(size_t j= 0; j < Ham.count_dn; j++)
        {
            size_t ind = (j*Ham.count_up)+i; //finding appropriate Fock state
            //switched i and j above 08/30/16 now have correct index
            size_t bas_up = Ham.basis_up[i];
            size_t bas_dn = Ham.basis_down[j];
            
            complex<double> cf = conj(D.G_state(ind))*D.G_state(ind);
            sum += cf.real();
            
            
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
double ObservableFunctions<Tnum>::Occupation_AnyLevel(int spin, const Hamiltonian<Tnum> &Ham, Eigen::VectorXcd Evec)
{
    double Ntotal = 0.0;
    cout << "Evec: " << Evec << endl;
    
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
            complex<double> cf = conj(Evec(ind))*Evec(ind);
            sum += cf.real();
            
            
            for(int n = 0; n < Ham.L; n++)
            {
                if((spin == 0) && MY_bittest(bas_up, n))//testing if up particle in Fock state on site n
                {
                    //n_up.at(n) += cf.real();
                    Ntotal += cf.real();
                }
                if((spin == 1) && MY_bittest(bas_dn, n))//testing if down particle in Fock state on site n
                {
                    //n_dn.at(n) += cf.real();
                    Ntotal += cf.real();
                }
            }
        }
        
    }
    
    
    cout << "Sum: " << sum << endl;
    
    
    
    return Ntotal;
}
template<typename Tnum>
double ObservableFunctions<Tnum>::DensityWCorr(const Hamiltonian<Tnum> &Ham, Lanczos_Diag<Tnum> &D, int cut)
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
            complex<double> cf = conj(D.G_state(ind))*D.G_state(ind);
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
double ObservableFunctions<Tnum>::DensityWCorr_O2(const Hamiltonian<Tnum> &Ham, Lanczos_Diag<Tnum> &D, int cut)
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
            
            complex<double> cf = conj(D.G_state(ind))*D.G_state(ind);
            
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
double ObservableFunctions<Tnum>::DensityCorrelation(double bu, double bd, complex<double> cf, size_t site1, size_t site2)
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
double ObservableFunctions<Tnum>::DensityCorr_O2(double bu, double bd, complex<double> cf, size_t site1, size_t site2)
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
double ObservableFunctions<Tnum>::OnsiteDensity_O2(double bu, double bd, complex<double> cf, size_t site)
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
void ObservableFunctions<Tnum>::SpinCorr(const Hamiltonian<Tnum> &Ham, Lanczos_Diag<Tnum> &D, std::ofstream &output, double t, int cut)
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
            
            complex<double> cf = conj(D.G_state(ind))*D.G_state(ind);
            
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
double ObservableFunctions<Tnum>::Calc_SC(double b1, double b2, std::complex<double> cf, size_t site1, size_t site2)
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
double ObservableFunctions<Tnum>::Calc_SameSpin(double bs, std::complex<double> cf, size_t site1, size_t site2)
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
complex<long double> ObservableFunctions<Tnum>::Expect_Cij(Hamiltonian<Tnum> &Ham, Lanczos_Diag<Tnum> &D, int spinspec, int TwoTunnel, size_t count, size_t count_opp,vector<size_t> basis, vector<size_t> index, Eigen::VectorXcd EigenState, size_t s1, size_t s2)
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
                size_t r,s;
                
                if ( spinspec == 0 ){
                    r = Ham.TotalIndex(pi, k);//(k*count) + p_ind;
                    s = Ham.TotalIndex(qi, k);//(k*count) + l_ind;
                    
                }else if ( spinspec == 1 ){
                    r = Ham.TotalIndex(k, pi);
                    s = Ham.TotalIndex(k, qi);
                    
                }else{
                    //cout << "More than 2 species fermion!!" << endl;
                }
                
                //cf += conj(G_state(s))*G_state(r)*Ham.Ham_Tot.coeffRef(s,r); //sum_nm <conj(c_n)*c_m*Hnm>
                if(TwoTunnel == 0)
                {
                    cf += conj(EigenState(s))*EigenState(r)*Ham.Ham_Tot.coeffRef(s,r);}
                else if(TwoTunnel == 1)
                {
                    cf += conj(EigenState(s))*EigenState(r);
                }
                else{
                    cout << "Not a valid option\n";
                }
                //Is above the correct way to handle the Hamiltonian in Fock basis?
                //cout << "r: " << r << " s: "<< s << endl;
                //cf += conj(G_state(s))*G_state(r);
                //cout << "coefficient: " << cf << endl;
            }
        }
    }
    
    //cfh = Ham.J1*exp(I*phi)*cf;
    
    //J = -2.*cfh.imag();
    
    return cf;
}

template<typename Tnum>
complex<long double> ObservableFunctions<Tnum>::Number(const Hamiltonian<Tnum> &Ham, Lanczos_Diag<Tnum> &D, int spinspec, size_t count, size_t count_opp,vector<size_t> basis, vector<size_t> index, Eigen::VectorXcd EigenState, size_t s)
{
    complex<long double> cf;
    
    for(size_t bs = 0; bs < count; bs++)
    {
        size_t pb = basis[bs];
        
        size_t pi = index[pb];
        
        
        if (MY_bittest(pb, s))
        {
            
            for(size_t k = 0; k < count_opp; k++)
            {
                size_t r;
                if ( spinspec == 0 ){
                    r = Ham.TotalIndex(pi, k);//(k*count) + p_ind;
                    //(k*count) + l_ind;
                    
                }else if ( spinspec == 1 ){
                    r = Ham.TotalIndex(k, pi);
                    
                    
                }else{
                    //cout << "More than 2 species fermion!!" << endl;
                }
                
                cf += conj(EigenState(r))*EigenState(r);//I don't know which one is the correct one
            }
        }
    }
    
    
    
    return cf;
    
}

template<typename Tnum>
complex<long double> ObservableFunctions<Tnum>::NumberNumber(const Hamiltonian<Tnum> &Ham, Lanczos_Diag<Tnum> &D,int spinspec, size_t count, size_t count_opp, std::vector<size_t> basis, std::vector<size_t> index, Eigen::VectorXcd EigenState, size_t i, size_t j)
{
    
    complex<long double> cf;
    double ni;
    double nj;
    
    for(size_t bs = 0; bs < count; bs++)
    {
        size_t pb = basis[bs];
        
        size_t pi = index[pb];
        
        
        if (MY_bittest(pb, i))
        {ni = 1.;}
        else{ni = 0.;}
        if (MY_bittest(pb, j))
        {nj = 1.;}
        else{nj = 0.;}
        
        for(size_t k = 0; k < count_opp; k++)
        {
            size_t r;
            if ( spinspec == 0 ){
                r = Ham.TotalIndex(pi, k);//(k*count) + p_ind;
                //(k*count) + l_ind;
                
            }else if ( spinspec == 1 ){
                r = Ham.TotalIndex(k, pi);
                
                
            }else{
                //cout << "More than 2 species fermion!!" << endl;
            }
            
            cf += conj(EigenState(r))*EigenState(r)*ni*nj;//only added if both ni,nj = 1
        }
        
    }
    
    
    return cf;
}

template<typename Tnum>
void ObservableFunctions<Tnum>::TotalCurrents(Hamiltonian<Tnum> &Ham, Lanczos_Diag<Tnum> &D, size_t s1, size_t s2)
{
    complex<double> Cij_up;
    complex<double> Cij_dn;
    //Phi = 0 for SOC case
    Jup = 0.0;
    Jdn = 0.0;
    
    //cout <<"Phi: "<< Ham.Phi_t << endl;
    
    
    Cij_up = Expect_Cij(Ham, D, 0, 0,Ham.count_up, Ham.count_dn, Ham.basis_up, Ham.index_up, D.EVMat.col(0), s1, s2);
    //cout << "Down\n";
    Cij_dn = Expect_Cij(Ham, D, 1, 0, Ham.count_dn, Ham.count_up, Ham.basis_down, Ham.index_dn, D.EVMat.col(0), s1, s2);
    //cout << "In the correct loop\n";
    
    
    Jup = -2.*Cij_up.imag();
    Jdn = -2.*Cij_dn.imag();
    
    cout << "Jup: " << Jup << " Jdn: " << Jdn << endl;
    
    double JQ = Jup.real()+Jdn.real();
    double JS = Jup.real()-Jdn.real();
    cout << "JQ: " << JQ << " JS: " << JS << endl;
    
}

template<typename Tnum>
double ObservableFunctions<Tnum>::Current_UPspin(Hamiltonian<Tnum> &Ham, Lanczos_Diag<Tnum> &D, int TwoTunnel,int EigenNum, double phi, size_t s1, size_t s2)
{
    Jup = 0.0;
    
    complex<long double> NN_H_Element = Ham.J1*Ham.J1*exp(2.*I*phi);
    //cout << "phi: " << phi << " HamElement: " << NN_H_Element << endl;
    complex<long double> Cij =  Expect_Cij(Ham, D, 0, TwoTunnel, Ham.count_up, Ham.count_dn, Ham.basis_up, Ham.index_up, D.EVMat.col(EigenNum), s1, s2);
    
    //cout << "Cij: " << Cij << endl;
    
    complex<double> J;
    
    if(TwoTunnel == 0){
        //calculate current from correlation matrix
        J = -2.*Cij.imag();//Hamiltonian element called in Expect_Cij Algorithm. This is my general current function
    }
    else{
        complex<long double> PL = NN_H_Element*Cij;
        J = -2.*PL.imag();//special case
    }
    
    Jup = J;
    
    return J.real();
}

template<typename Tnum>
double ObservableFunctions<Tnum>::Current_DNspin(Hamiltonian<Tnum> &Ham, Lanczos_Diag<Tnum> &D, int TwoTunnel, int EigenNum, double phi, size_t s1, size_t s2)
{
    Jdn = 0.0;
    
    complex<long double> NN_H_Element = Ham.J1*Ham.J1*exp(2.*I*phi);
    
    complex<long double> Cij= Expect_Cij(Ham, D, 1, TwoTunnel, Ham.count_dn, Ham.count_up, Ham.basis_down, Ham.index_dn, D.EVMat.col(EigenNum), s1, s2);
    complex<double> J;
    if(TwoTunnel == 0){
    
        J = -2.*Cij.imag();//this expression is called except for a special case
    }
    else{
        J = -2.*NN_H_Element.imag()*Cij.imag();
    }
    
    Jdn = J;
    
    return J.real();
}

template<typename Tnum>
double ObservableFunctions<Tnum>::SpinCurrent()
{
    return Jup.real()-Jdn.real();
}

template<typename Tnum>
double ObservableFunctions<Tnum>::ChargeCurrent()
{
    return Jup.real()+Jdn.real();
}

template<typename Tnum>
complex<double> ObservableFunctions<Tnum>::CurrentVariance(Hamiltonian<Tnum> &Ham, Lanczos_Diag<Tnum> &D, int EigenNum, int spec, size_t s1, size_t s2)
{
    
    complex<double> Jsquare;
    complex<double> Var;
    
    //TotalCurrents(Ham, D, s1, s2);//generate Jup and Jdn, You might want to get rid of this function for speed up
    
    if(spec == 0)//spin 1/2
    {
        //calculation of <J^2>
        Jsquare = CurrentSquare(Ham, D, EigenNum, 0, s1, s2);
        
        //calculation of Variance \deltaJ^2 = <J^2> - <J>^2
        Var = Jsquare - (Jup*Jup);
    }
    else if (spec == 1)//spin -1/2
    {
        Jsquare = CurrentSquare(Ham, D, EigenNum, 1, s1, s2);
        Var = Jsquare - (Jdn*Jdn);
    }
    else{
        cout << "Not a spin species\n";
    }
    
    
    return Var;//this gets a square root in the main file
}

template<typename Tnum>
complex<long double> ObservableFunctions<Tnum>::CurrentSquare(Hamiltonian<Tnum> &Ham, Lanczos_Diag<Tnum> &D, int EigenNum, int spec, size_t s1, size_t s2)
{
    complex<long double> Jsq;
    
    complex<long double> N1;
    complex<long double> N2;
    complex<long double> N12;
    complex<long double> C12;
    
    if(spec == 0)//spin 1/2
    {
        N1 = Number(Ham, D, 0, Ham.count_up, Ham.count_dn, Ham.basis_up, Ham.index_up, D.EVMat.col(EigenNum), s1);//this is the number on site 1
        N2 = Number(Ham, D, 0, Ham.count_up, Ham.count_dn, Ham.basis_up, Ham.index_up, D.EVMat.col(EigenNum), s2);//this is the number on site 2
        N12 = NumberNumber(Ham, D, 0, Ham.count_up, Ham.count_dn, Ham.basis_up, Ham.index_up, D.EVMat.col(EigenNum), s1, s2);// this is my number-number <n_1n_2>
        C12 = Expect_Cij(Ham, D, 0, 1, Ham.count_up, Ham.count_dn, Ham.basis_up, Ham.index_up, D.EVMat.col(EigenNum), s1, s2);//the Hamiltonian element is not called here, and this is only used for Wick
    }
    else if (spec == 1)//spin -1/2
    {
        N1 = Number(Ham, D, 1, Ham.count_dn, Ham.count_up, Ham.basis_down, Ham.index_dn, D.EVMat.col(EigenNum), s1);
        N2 = Number(Ham, D, 1, Ham.count_dn, Ham.count_up, Ham.basis_down, Ham.index_dn, D.EVMat.col(EigenNum), s2);
        N12 = NumberNumber(Ham, D, 1, Ham.count_dn, Ham.count_up, Ham.basis_down, Ham.index_dn, D.EVMat.col(EigenNum), s1, s2);
        C12 = Expect_Cij(Ham, D, 1, 1, Ham.count_dn, Ham.count_up, Ham.basis_down, Ham.index_dn, D.EVMat.col(EigenNum), s1, s2);
        
    }
    else{
        cout << "Not a spin species\n";
    }
    
    complex<long double> tbar;
    complex<long double> Two = 2.;
    if(s2-s1 == 1)//Current Fluctutation on NN links
    {
        tbar = Ham.J1;
        //cout << "tbar should be t1: " << tbar << endl;
    }
    else if(s2-s1 == 2)
    {
        tbar = Ham.J2;
        //cout << "tbar should be t2: " << tbar << endl;
    }
    
    //<J^2> = |t12|^2(<n1> + <n2> - 2<n1n2>)
    Jsq = tbar*conj(tbar)*(N1 + N2 - (Two*N12));//No wick
    
    //Jsq = tbar*conj(tbar)*(N1 + N2 -Two*N1*N2 + Two*C12*conj(C12));//Wick decomposition
//    cout << "N1: " << N1 << endl;
//    cout << "N2: " << N2 << endl;
//    cout << "C12: " << C12 << endl;
    
    return Jsq;
}

template<typename Tnum>
void ObservableFunctions<Tnum>::OutputOccupation(ofstream &output, const Hamiltonian<Tnum> &Ham, Lanczos_Diag<Tnum> &D)
{
    vector<double> n1;
    vector<double> n2;
    
    //it might be worth changing these into vectors later
    for(int it = 0; it < Ham.Tot_base; it++)
    {
        double nu = Occupation_AnyLevel(0, Ham, D.EVMat.col(it));
        double nd = Occupation_AnyLevel(1, Ham, D.EVMat.col(it));
        
        n1.push_back(nu);
        n2.push_back(nd);
        
        output << it << " " << nu << " "  << nd << endl;
        cout << "it: " << it << " nu: " << nu << " nd: "  << nd << endl;
        cout << "Sz: " << (nu-nd)/2.  << " Qz: " << (nu+nd)/2. - Ham.L/2. << endl;
        
    }
    
}



template class ObservableFunctions<double>;
template class ObservableFunctions<complex<double> >;

