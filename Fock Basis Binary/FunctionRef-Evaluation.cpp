//
//  FunctionRef-Evaluation.cpp
//  Fock Basis Binary
//
//  Created by mekena McGrew on 8/25/16.
//  Copyright © 2016 Mekena Metcalf. All rights reserved.
//

#include <stdio.h>
#include "FunctionRef-Evaluation.h"

void Write_Density(ofstream &fout, vector<double> &n_up, vector<double> &n_dn, int L )
{
    for(int i = 0; i < L; i++)
    {
        fout << i << " " << n_up[i] << " " << n_dn[i] << endl;
        
    }
    fout << endl;
    
}

void DensityDiff(ofstream &fout, vector<double> &n_up, vector<double> &n_dn, int L )
{
    for(int i = 0; i < L; i++)
    {
        fout << i << " " << n_up[i] - n_dn[i] << endl;
        
    }
    fout << endl;
}

void DensityInTime(ofstream &fout, vector<double> &n, int L, double Phi )
{
    for(int i = 0; i < L; i++)
    {
        fout << i << " " << Phi << " "  << n[i] << endl;
        
    }
    
}

void Density(ofstream &fout, vector<double> &n, int L)
{
    for(int i = 0; i < L; i++)
    {
        fout << i << " "<< n[i] << endl;
        
    }
    fout << endl;
}

void EdgeDensity(ofstream &fout, vector<double> &n, double time)
{
    fout << time << " " << n[0] << endl;
    //cout << time << " " << n[0] << endl;
}

void Harmonic_Trap(vector<double> &HT, int L, double y)
{
    int Ro;
    if(L % 2 == 0)
    {
        
        Ro = (L)/2;
    }
    else
    {
        Ro = (L +1)/2;
    }
    
    //cout << "R0: " << Ro << endl;
    for(int i = 1; i <= L; i++)
    {
        double val = y*((i)-Ro)*((i)-Ro);
        HT.push_back(val);
        //cout << HT[i-1] << endl;
    }
    
    
}

template<typename Tnum>
void Fidelity(ofstream &output, Hamiltonian<Tnum> &h, Lanczos_Diag<Tnum> &d, double U, double J1, double J2, double Um, double Jm, int L)

{
    typedef Eigen::Matrix<Tnum, Eigen::Dynamic, 1> VectorType;
    
    double dU = 1.;
    double dJ = .1;
    int Jtot = Jm/dJ;
    int Utot = Um/dU;
    double g;
    double F;
    VectorType Gstate;
    VectorType Gstate_Temp;
    // U = 0;
    
    for(int i = 0; i <= Utot ; i++)
    {
        cout << "Loop: " << i << endl;
        for(int j = 0; j <= Jtot; j++)
        {
            
            h.Set_Const(J1, J2, U);//not having U will be a problem
            h.HopMatrix_Build();
            if(j == 0)
            {
                h.IntMatrix_Build();
            }
            h.Total_Ham();
            d.Set_Mat_Dim_LA(h);
            d.Diagonalize(h);
            Gstate = d.SendGstate();//hmmm...
            if(j > 0)
            {
                F = GdotG(Gstate, Gstate_Temp);
                
                g = 2*(1-F)/(L*(dU*dU));//this is wrong
                output << U << " " << J2 << " " << g << endl;
                
            }
            Gstate_Temp = Gstate;
            J2 += dJ;
            h.ClearHopTriplet();
            
        }
        U += dU;
        J2 = 0;
        h.ClearInteractTriplet();
    }
    
    
}

template<typename Tnum>
void U_Fid(ofstream &output, Hamiltonian<Tnum> &h, Lanczos_Diag<Tnum> &d , double U, double J1, double J2, double Um, int L)
{
    typedef Eigen::Matrix<Tnum, Eigen::Dynamic, 1> VectorType;
    
    double ddU = .0001;
    double dU = 1.;
    double UU;
    int Utot = Um/dU;
    double g;
    double glog;
    double F;
    VectorType Gstate_i;
    VectorType Gstate_f;
    // U = 0;
    
    for(int i = 0; i <= Utot ; i++)
    {
        cout << "Loop: " << i << " U: "<< U<< endl;
        
        
        h.Set_Const(J1, J2, U);//change value of U
        
        if(i == 0)
        {
            h.HopMatrix_Build();//only need to build hopham once
        }
        
        //get preliminary Gstate
        h.IntMatrix_Build();//build new matrix for every iteration
        h.Total_Ham();//combine hopham and hopint
        
        d.Set_Mat_Dim_LA(h);
        d.Diagonalize(h);
        Gstate_i = d.SendGstate();
        
        //add ddu to U and get next G-state
        UU = U+ddU;//increase U
        cout << "U+ddu: "<< UU << endl;
        h.ClearInteractTriplet();//clear interaction matrix
        h.Set_Const(J1, J2, UU);//reset constant with U+ddU
        h.IntMatrix_Build();//build new matrix for every iteration
        h.Total_Ham();//combine hopham and hopint
        d.Set_Mat_Dim_LA(h);
        d.Diagonalize(h);
        Gstate_f = d.SendGstate();//Gstate of U+ddU
        
        //get fidelity between U and U+ddU
        
        F = GdotG(Gstate_f, Gstate_i);
        g = 2.*(1.-F)/(L*(ddU*ddU));//changing variable is U
        glog = (-(2./L)*log(F))/(ddU*ddU);
        output << U << " " << g << endl;
        
        if(i == 0)
        {
            cout << scientific;
            cout << "F: " << F << " g: " << g << " glog: " << glog <<"Log(F) " << log(F) <<  endl;
        }
        
        //clear interaction matrix for next iteration
        U++;
        h.ClearInteractTriplet();
    }
    
    
}

template<typename Tnum>
void Fidelity_HT(ofstream &output, Hamiltonian<Tnum> &h, Lanczos_Diag<Tnum> &d , double U, double J1, double J2, double Um, double Ym, vector<double> HT, int L)
{
    typedef Eigen::Matrix<Tnum, Eigen::Dynamic, 1> VectorType;
    
    double dU = 1.;
    double dw = .01;
    double y = 0;
    int Wtot = Ym/dw;
    int Utot = Um/dU;
    long double g;
    long double F;
    VectorType Gstate;
    VectorType Gstate_Temp;
    // U = 0;
    
    
    for(int w = 0; w <= Wtot; w++)
    {
        cout << "iteration: " << w << endl;
        Harmonic_Trap(HT, L, y); //this is correct
        //            cout << "Harmonic Trap: " << y << endl;//this is correct
        //            for(int i = 0; i < HT.size(); i++)
        //            {
        //                cout << HT[i] << endl;
        //            }
        h.GetHarmTrap(HT);
        h.Set_Const(J1, J2, U);//not having U will be a problem
        h.HopMatrix_Build();
        if(w == 0)
        {
            h.IntMatrix_Build();
        }
        h.Total_Ham();
        d.Set_Mat_Dim_LA(h);
        d.Diagonalize(h);
        Gstate = d.SendGstate();//ground state is flipped compared to master
        //cout << "Gstate: " << Gstate << endl;
        if(w > 0)
        {
            F = GdotG(Gstate, Gstate_Temp);
            cout << "F: " << F << endl;
            g = 2*(1-F)/(L*(dw*dw));//this is wrong
            output << y << " " << g << endl;
            
        }
        Gstate_Temp = Gstate;
        y += dw;
        h.ClearHopTriplet();
        HT.clear();
        
    }
    //        U += dU;
    //        J2 = 0;
    //        h.ClearInteractTriplet();
    //    }
    
}

template<typename Tnum>
void PeierlsTD(ofstream &out_up, ofstream &out_dn, Hamiltonian<Tnum> &h, Lanczos_Diag<Tnum> &d, int T_it, double dt, double tp, int L)
{
    int Nflag = 0;
    double t = 0.;
    double Pi = (4*atan(1.0));
    double Phi_max = Pi/4.;
    double Phi_t = 0.0;
    
    for(int it = 0; it <= T_it; it++)
    {
        t = it*dt;
        cout << "t: " << t << endl;
        
        
        if(t <= tp)
        {
            Phi_t = (t/tp)*Phi_max;
            
            h.GetPhi(Phi_t);
            h.HopMatrix_Build_Peierls();
            
            h.Total_Ham();
            if(it == 0)
            {
                d.Set_Mat_Dim_LA(h);
                d.Diagonalize(h);
                
                d.Density(h);
                //Write_Density(out, d.n_up, d.n_dn, L);
                Density(out_up, d.n_up, L);
                //Density(out_dn, d.n_dn, L);
                EdgeDensity(out_dn, d.n_up, t);
                
            }
            else
            {
                d.Dynamics(h);
                //d.DebugDynamics(h);
            }
            
            if(t < tp)
            {
                h.ClearHopTriplet();
            }
        }
        else
        {
            //cout << "in loop 2\n";
            //ham.OutHam();
            d.Dynamics(h);
            //d.DebugDynamics(h);
        }
        
        
        
        if(Nflag == 100)
        {
            cout << "Getting Density for t=" << t << endl;
            d.Density(h);
            Density(out_up, d.n_up, L);
            //            Density(out_dn, d.n_dn, L);
            EdgeDensity(out_dn, d.n_up, t);
            
            //DensityDiff(PDout, Diag.n_up, Diag.n_dn, Nsite);
            Nflag = 0;
        }
        
        Nflag++;
        
        
    }
    
}

template<typename Tnum>//why does this program work and not the other???????
void QPump_TD(ofstream &out_up, ofstream &out_dn, Hamiltonian<Tnum> &ham, Lanczos_Diag<Tnum> &diag, const int T, const double dt, const double h_0, const double J_0, const double d_0, double U, int L)
{
    double delta;
    double J1, J2;
    double t;
    double h;
    
    diag.TimeEvoCoeff(dt);
    ham.ClearHopTriplet();
    
    int T_it = T /dt;
    
    const double omega = (8*atan(1.0))/T;
    
    for( int it = 0; it <= T_it; it++)
    {
        t = it*dt;
        
        if( it % 100 == 0)
        {
            cout << "t: " << t << endl;
        }
        
        delta = d_0*cos(omega*t);
        J1 = J_0 + delta;
        J2 = J_0 - delta;
        //cout << "J_1: " << J1 << " J2: " << J2 << endl;
        h = h_0*sin(omega*t);
        //cout << "h: " << h << endl;
        ham.GetOnsite(h);
        
        ham.Set_Const(J1, J2, U);
        ham.HopMatrix_Build_QPump();
        ham.Total_Ham();
        //ham.OutHam();
        
        if(it == 0)
        {
            diag.Set_Mat_Dim_LA(ham);
            diag.Diagonalize(ham);
        }
        else
        {
            //cout << "Beginning Dynamics" << endl;
            diag.Dynamics(ham);
            //diag.DebugDynamics(ham);
        }
        
        ham.ClearHopTriplet();
        //cout << "Triplet Clear\n";
        if(it == 0 || it == T_it/2. || it == T_it)
        {
            cout << "Getting Density for t=" << t << endl;
            //cout << "Hamiltonian: \n" << ham.Ham_Tot << endl;
            diag.Density(ham);
            //Write_Density(out, diag.n_up, diag.n_dn, L);
            Density(out_up, diag.n_up, L);
            Density(out_dn, diag.n_dn, L);
            
        }
    }
}

double GdotG( const VectorXd &G1, const VectorXd &G2)
{
    double F;
    F = abs(G1.dot(G2));
    
    //cout << "Fidelity: " << F << endl;
    return F;
}
double GdotG( const VectorXcd &Gc1, const VectorXcd &Gc2)
{
    complex<double> F;
    complex<double> Fp;
    F = Gc1.dot(Gc2);
    Fp = F*conj(F);
    F = sqrt(Fp);
    
    return F.real();
}

template class Hamiltonian<double>;
template class Hamiltonian<complex<double> >;
