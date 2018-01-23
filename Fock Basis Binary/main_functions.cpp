//
//  main_functions.cpp
//  Fock Basis Binary
//
//  Created by mekena McGrew on 1/19/18.
//  Copyright © 2018 Mekena Metcalf. All rights reserved.
//

#include "main_functions.h"

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
    double dU = 0.1;
    double lambda = ddU/dU;
    double UU;
    int Div = 1/dU;
    int Utot = Um*Div;
    double g;
    double glog;
    double F;
    VectorType Gstate_i;
    VectorType Gstate_f;
    // U = 0;
    
    for(int i = 0; i <= Utot ; i++)
    {
        if( i % Div == 0)
        {
            cout << "Loop: " << i << " U: "<< U<< endl;
        }
        
        
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
        cout << scientific;
        cout << "Fidelity: " << F << endl;
        g = 2.*(1.-F)/(L*(lambda*lambda));//changing variable is U
        glog = (-(2./L)*log(F))/(lambda*lambda);
        output << U << " " << g << endl;
        
        //        if(i == 0)
        //        {
        //            cout << scientific;
        //            cout << "F: " << F << " g: " << g << " glog: " << glog <<"Log(F) " << log(F) <<  endl;
        //        }
        
        //clear interaction matrix for next iteration
        U += dU;
        h.ClearInteractTriplet();
    }
    
    
}

template<typename Tnum>
void U_Fid_V2(ofstream &output, Hamiltonian<Tnum> &h, Lanczos_Diag<Tnum> &d , double U, double J1, double J2, double Um, int L)
{
    typedef Eigen::Matrix<Tnum, Eigen::Dynamic, 1> VectorType;
    
    
    double dU = 0.01;
    
    
    int Div = 1/dU;
    int Utot = Um*Div;
    double g;
    double glog;
    double F;
    VectorType Gstate;
    VectorType Gstate_temp;
    VectorXd V;
    // U = 0;
    
    for(int i = 0; i <= Utot ; i++)
    {
        if( i % Div == 0)
        {
            cout << "Loop: " << i << " U: "<< U<< endl;
        }
        
        
        h.Set_Const(J1, J2, U);//change value of U
        
        if(i == 0)
        {
            //h.HopMatrix_Build();//only need to build hopham once
            h.HopMatrix_Build_Periodic();
        }
        
        //get preliminary Gstate
        h.IntMatrix_Build();//build new matrix for every iteration
        h.Total_Ham();//combine hopham and hopint
        
        d.Set_Mat_Dim_LA(h);
        //d.Diagonalize(h);
        V = d.FullDiagonalization(h);//only for small L
        Gstate = d.SendGstate();
        
        //get fidelity between U and U-dU
        if( i != 0)
        {
            F = GdotG(Gstate, Gstate_temp);
            cout << scientific;
            //cout << "Fidelity: " << F << endl;
            g = 2.*(1.-F)/(L*(dU*dU));//changing variable is U
            glog = (-(2./L)*log(F))/(dU*dU);
            output << U << " " << g << endl;
            cout << "U: " << U << endl;
        }
        //clear interaction matrix for next iteration
        Gstate_temp = Gstate;
        U += dU;
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
            g = 2.*(1.-F)/(L*(dw*dw));//this is wrong
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
        if (it%100 == 0)
        {
            cout << "t: " << t << endl;
            
        }
        
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
                d.Dynamics(h, dt);
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
            d.Dynamics(h, dt);
            //d.DebugDynamics(h);
        }
        
        
        
        if(Nflag == 10)
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

template<typename Tnum>
void CorrTD(ofstream &output, ofstream &out, Hamiltonian<Tnum> &h, Lanczos_Diag<Tnum> &d, int cut, int T_it, double dt, int L)
{
    double t = 0.;
    
    double C;
    
    for(int it = 0; it <= T_it; it++)
    {
        t = it*dt;
        
        cout << "it: " << it << endl;
        
        if(it % 100 == 0)
        {
            cout << "t: " << t << endl;
        }
        
        
        if(it == 0)
        {
            
            h.HopMatrix_Build_Periodic();
            //h.HopMatrix_Build();
            
            h.Total_Ham();
            
            d.Set_Mat_Dim_LA(h);
            d.Diagonalize(h);
            
            //C = d.DensityWCorr_O2(h);
            C = d.DensityWCorr(h, cut);
            //d.SpinCorr(h,output, t, cut);//spin correlation
            
            output << t << " " << C << endl;
            
            //h.ClearHopTriplet();
            
        }
        else if(it == 1)
        {
            
            //h.HopMatrix_Build();
            //h.Total_Ham();
            h.MakeCut(cut);
            d.Dynamics(h, dt);
            //d.DebugDynamics(h);
            
            C = d.DensityWCorr(h, cut);  //1st order
            //C = d.DensityWCorr_O2(h); //2nd order
            //d.SpinCorr(h,output, t, cut);//spin correlation
            output << t << " " << C << endl;
        }
        else{
            d.Dynamics(h, dt);
            
            if( (it % 50) == 0)
            {
                C = d.DensityWCorr(h, cut); //1st order
                //C = d.DensityWCorr_O2(h); //2nd order
                //d.SpinCorr(h,output, t, cut);//spin correlation
                output << t << " " << C << endl;
            }
        }
        
        //        if(it%300 == 0)
        //        {
        //            if(cut == 0)
        //            {Write_Density(out, d.n_up, d.n_dn, L);}
        //            else if(cut == 1)
        //            {
        //                for(int i = 1; i <= L; i++)
        //                {
        //                    if( i < L)
        //                    {
        //                    out << i << " " << d.n_up[i] << " " << d.n_dn[i] << endl;
        //                    }
        //                    else{
        //                        out << i << " " << d.n_up[0] << " " << d.n_dn[0] << endl;
        //                    }
        //
        //                }
        //                out << endl;
        //            }
        //            else{
        //                out << "Wrong Cut\n";
        //            }
        //
        //        }
    }
    
}

void EVProbeGamma(ofstream &output1, ofstream &output2, ofstream &output3, Hamiltonian<complex<double>> &h, Lanczos_Diag<complex<double>> &d, int s1, int s2, int L, double phi1, double phi2, double link)
{
    
    double dg = 0.01;
    double Gm = 10.0/dg;
    double gamma = 0.0;
    Eigen::VectorXd E;
    
    
    //cout << "TB: " << h.Tot_base << endl;
    
    
    for(int g = 0; g <= Gm; g++)
    {
        gamma = g*dg;
        cout << "gamma: " << gamma << endl;
        // h.HopMatrix_Build_PeriodicWithSOC(s1, s2, gamma, Pi/2.);
        h.HopMatrix_Build_PeriodicNNNHop(link, gamma, Pi/2.);
        //h.HopMatrix_Build_HaldHam_NoGauge(link, gamma, phi1, phi2);
        if (g == 0)
        {
            h.IntMatrix_Build();
        }
        h.Total_Ham_WSOC();
        
        E = d.FullDiagonalization(h);
        //output Eigenvalues
        output1 << gamma << " " << E(0) <<" " << E(1) << " " << E(2)<< " " << E(3) <<" " << E(4) << " " << E(5) << endl;
        //output1 << gamma << " " << E(0) << endl;
        //cout << gamma << " " << E(0) << endl;
        
        //output Currents
        //        d.TotalCurrents(h, 0, 1);//outer current
        //        output1 << gamma << " " << d.SpinCurrent() << " " << d.ChargeCurrent() << endl;
        //        d.TotalCurrents(h, 0, 2);//inner current
        //        output2 << gamma << " " << d.SpinCurrent() << " " << d.ChargeCurrent() << endl;
        //        output2 << gamma << " " << d.Current_UPspin(h, 1, 0, 1) + d.Current_DNspin(h,1, 0, 1) << " " << d.Current_UPspin(h,1, 0, 2)+d.Current_DNspin(h,1, 0, 2) << " " << d.Current_UPspin(h, 1, 2, 3)+d.Current_DNspin(h, 1, 2, 3) << endl;
        //        output3 << gamma << " " << d.Current_UPspin(h, 1, 0, 1) - d.Current_DNspin(h, 1, 0, 1) << " " << d.Current_UPspin(h, 1, 0, 2)-d.Current_DNspin(h, 1, 0, 2) << " " << d.Current_UPspin(h, 1, 2, 3)-d.Current_DNspin(h, 1, 2, 3) << endl;
        output2 << gamma << " " << d.Current_UPspin(h, 0, 0, 1) << " " << d.Current_UPspin(h,0, 0, 2) << " " << d.Current_UPspin(h, 0, 2, 3)<< endl;
        output3 << gamma << " " << d.Current_DNspin(h, 0, 0, 1) << " " << d.Current_DNspin(h, 0, 0, 2) << " " << d.Current_DNspin(h, 0, 2, 3) << endl;
        
        h.ClearHopTriplet();
        
    }
    
}

void EVProbeInt(ofstream &output1, ofstream &output2, ofstream &output3, Hamiltonian<complex<double>> &h, Lanczos_Diag<complex<double>> &d, int s1, int s2, int L, double gSO, double phi1, int link)
{
    double dU = 0.1;
    double Um = 10.0/dU;
    double U;
    //cout << "g: " <<  gSO << endl;
    
    Eigen::VectorXd E;
    //cout << "TB: " << h.Tot_base << endl;
    
    for(int u = 0; u <= Um; u++)
    {
        U = u*dU;
        h.Set_Const(1.0, 1.0, U);
        cout << "U: " << U <<endl;
        
        if (u == 0)
        {
            //h.HopMatrix_Build_PeriodicWithSOC(s1, s2, gSO, Pi/4.);
            h.HopMatrix_Build_PeriodicNNNHop(link, gSO, Pi/2.);
            //h.HopMatrix_Build_HaldHam_NoGauge(6, gSO, phi1, phi2);
        }
        h.IntMatrix_Build();
        h.Total_Ham_WSOC();
        
        //Eigenvalues
        E = d.FullDiagonalization(h);
        //output1 << U << " " << E(0) <<" " << E(1) << " " << E(2)<< " " << E(3) <<" " << E(4) << " " << E(5)  << endl;
        output1 << U << " " << E.transpose() << endl;
        //        for(int i = 0; i < h.Tot_base; i++)
        //        {
        //            if(i == 0)
        //            {
        //                output << U << " " << E(i) << " ";
        //            }
        //            else{
        //                output << E(i) << " ";
        //            }
        //        }
        //        output << endl;
        
        //Current
        output2 << U << " " << d.Current_UPspin(h, 0, 0, 1) + d.Current_DNspin(h, 0, 0, 1) << " " << d.Current_UPspin(h, 0, 0, 2)+d.Current_DNspin(h, 0, 0, 2) << endl;
        output3 << U << " " << d.Current_UPspin(h, 0, 0, 1) - d.Current_DNspin(h, 0, 0, 1) << " " << d.Current_UPspin(h, 0, 0, 2)-d.Current_DNspin(h, 0, 0, 2) << endl;
        
        h.ClearInteractTriplet();
        
    }
    
}


void EVProbePhase(ofstream &output1, ofstream &output2, ofstream &output3, Hamiltonian<complex<double>> &h, Lanczos_Diag<complex<double>> &d, ObservableFunctions<complex<double>> &O, int s1, int s2, int L, double gSO, int link, double a)
{
    long double Phi_max = 1.5;//1.5
    long double Phi_init = 0.5;//0.5
    double dp = .001;
    double P_it = (Phi_max-Phi_init)/dp;
    double phi1;
    double phi2 = Phi_init*Pi;
    Eigen::VectorXd E;
    
    Eigen::VectorXcd Gstate;
    Eigen::VectorXcd Gstate_temp;
    double F;
    double g;
    
    
    for(int p = 0; p <= P_it+1; p++)
    {
        //cout << p<< endl;
        //phi2 = dp*p*Pi;
        
        cout << "phi: " << phi2/Pi << endl;
        //phi1 = (4./((6.*a) - 2.))* phi2;
        //cout << "Phi: " << phi << endl;
        //h.HopMatrix_Build_PeriodicWithSOC(s1, s2, gSO, phi);
        //h.HopMatrix_Build_PeriodicNNNHop(link, gSO, phi2);
        h.HopMatrix_Build_PeriodicBfield(link, phi2);
        //h.HopMatrix_Build_HaldHam_NoGauge(6, gSO, phi1, phi2);
        
        if(p == 0)
        {
            h.IntMatrix_Build();
        }
        h.Total_Ham_WSOC();
        //Eigenvalues
        E = d.FullDiagonalization(h);
        
        //Get Currents
        //d.TotalCurrents(h, 2, 3);//is a void function
        //output << phi/Pi << " " << d.SpinCurrent() << " " << d.ChargeCurrent() << endl;
        //output << phi/Pi << " " << d.UpCurrent() << " " << d.DownCurrent() << endl;
        double Jup_NN_1 = O.Current_UPspin(h, d, 0, 0, phi2, 0, 1);
        complex<double> var_NN_1 = sqrt(O.CurrentVariance(h, d, 0, 0, 0, 1));
        
        double Jup_NN_FE = O.Current_UPspin(h, d, 0, 1, phi2, 0, 1);
        complex<double> var_NN_FE = sqrt(O.CurrentVariance(h, d, 1, 0, 0, 1));
        
        double Jup_NN_2 = O.Current_UPspin(h, d, 0, 0, phi2, 1, 2);
        complex<double> var_NN_2 = sqrt(O.CurrentVariance(h, d, 0, 0, 1, 2));
        
        double Jup_NNN = O.Current_UPspin(h, d, 0, 0, phi2, 0, 2);//TT is 0 for full current, 1 for uncertainty relation
        complex<double> var_NNN = sqrt(O.CurrentVariance(h, d, 0, 0, 0, 2));
        
        complex<double> Uncty_L = var_NN_1*var_NN_1*var_NN_2*var_NN_2;//left side of uncertainty eqn
        double Uncty_R = (1./4.)*Jup_NNN*Jup_NNN;//right side of uncertainty eqn
        
        
        
        //Return Currents
        //        output1 << phi/Pi << " " << d.Current_UPspin(h, 0, 1) << " " << d.Current_UPspin(h, 0, 2) << " " << d.Current_UPspin(h, 2, 3) << endl;
        //        output2 << phi/Pi << " " << d.Current_DNspin(h, 0, 1) << " " << d.Current_DNspin(h, 0, 2) << " " << d.Current_DNspin(h, 2, 3) << endl;
        //        output2 << phi2/Pi << " " << d.Current_UPspin(h, 0, 0, 1) + d.Current_DNspin(h, 0, 0, 1) << " " << d.Current_UPspin(h, 0, 0, 2)+d.Current_DNspin(h, 0, 0, 2) << endl;
        //        output3 << phi2/Pi << " " << d.Current_UPspin(h, 0, 0, 1) - d.Current_DNspin(h, 0, 0, 1) << " " << d.Current_UPspin(h, 0, 0, 2)-d.Current_DNspin(h, 0, 0, 2) << endl;
        output1 << phi2/Pi << " " << E.transpose() << endl;
        
        output2 << phi2/Pi << " " << Jup_NN_FE << " " << var_NN_FE.real() << endl;
        //output2 << Jup_NN_1 << " " << var_NN_1.real() << endl;
        //output3 << phi2/Pi << " " << Uncty_L.real() << " " << Uncty_R << endl;
        output3 << phi2/Pi << " " << Jup_NN_1 << " " << var_NN_1.real() << endl;
        
        h.ClearHopTriplet();
        
        phi2 += dp*Pi;
        //phi2 -= dp*Pi;
        
    }
}

void ProbePhase_wRealNNN(ofstream &output1, ofstream &output2, ofstream &output3, Hamiltonian<complex<double>> &h, Lanczos_Diag<complex<double>> &d, int s1, int s2, int L, double gSO, int link)
{
    double Phi_max = 0.5;
    double dp = .01;
    double P_it = Phi_max/dp;
    
    double phi;
    Eigen::VectorXd E;
    
    double t2;
    double lam_b;
    
    for(int p = 1; p < P_it; p++)
    {
        //cout << p<< endl;
        phi = p*Pi*dp;
        
        t2 = gSO/tan(phi);
        lam_b = sqrt((t2*t2)+(gSO*gSO));
        //cout << "Phi: " << phi <<" t2: " << t2 << endl;
        
        h.HopMatrix_Build_PeriodicNNNHop(link, lam_b, phi);
        
        if(p == 1)
        {
            h.IntMatrix_Build();
        }
        h.Total_Ham_WSOC();
        //Eigenvalues
        E = d.FullDiagonalization(h);
        output1 << phi/Pi << " " << t2 << " " << E.transpose() << endl;
        
        //Get Currents
        //d.TotalCurrents(h, 2, 3);//is a void function
        //output << phi/Pi << " " << d.SpinCurrent() << " " << d.ChargeCurrent() << endl;
        //output << phi/Pi << " " << d.UpCurrent() << " " << d.DownCurrent() << endl;
        
        //Return Currents
        //        output1 << phi/Pi << " " << d.Current_UPspin(h, 0, 1) << " " << d.Current_UPspin(h, 0, 2) << " " << d.Current_UPspin(h, 2, 3) << endl;
        //        output2 << phi/Pi << " " << d.Current_DNspin(h, 0, 1) << " " << d.Current_DNspin(h, 0, 2) << " " << d.Current_DNspin(h, 2, 3) << endl;
        output2 << phi/Pi << " " << t2 <<" " << d.Current_UPspin(h, 1, 0, 1) + d.Current_DNspin(h, 1, 0, 1) << " " << d.Current_UPspin(h, 1, 0, 2)+d.Current_DNspin(h, 1, 0, 2)<< " " << d.Current_UPspin(h, 1, 2, 3)+d.Current_DNspin(h, 1, 2, 3) << endl;
        output3 << phi/Pi <<" " << t2 << " " << d.Current_UPspin(h, 1, 0, 1) - d.Current_DNspin(h, 1, 0, 1) << " " << d.Current_UPspin(h, 1, 0, 2)-d.Current_DNspin(h, 1, 0, 2) << " " << d.Current_UPspin(h, 1, 2, 3)-d.Current_DNspin(h, 1, 2, 3) << endl;
        
        
        h.ClearHopTriplet();
        
    }
}

void CurrentVsGamma(ofstream &out1, ofstream &out2, ofstream &out3, Hamiltonian<complex<double>> &h, Lanczos_Diag<complex<double>> &d, double a, int links)
{
    double Phi_max = 0.5;
    double dp = .01;
    double P_it = Phi_max/dp;
    double phi1;
    double phi2;
    
    double dg = 0.1;
    double Gm = 5.0/dg;
    double gamma = 0.0;
    double t2;
    double gamma_p;
    
    
    Eigen::VectorXd E;
    
    for(int git = 0;  git <= Gm; git++)
    {
        gamma = git*dg;
        cout << "gamma: " << gamma << endl;
        
        for(int p = 1; p < P_it; p++)
        {
            //cout << p<< endl;
            phi2 = dp*p*Pi;
            //phi1 = (4./(6.*a -2.))* phi2;
            t2 = gamma/tan(phi2);
            gamma_p = sqrt((t2*t2)+(gamma*gamma));
            
            h.HopMatrix_Build_PeriodicNNNHop(links, gamma_p, phi2);
            //h.HopMatrix_Build_HaldHam_NoGauge(6, gamma, phi1, phi2);
            if(p == 1)
            {
                h.IntMatrix_Build();
            }
            h.Total_Ham_WSOC();
            //Eigenvalues
            E = d.FullDiagonalization(h);
            
            out1 << gamma << " " << phi2/Pi << " " << E(0) << endl;
            out2 << gamma << " " << phi2/Pi << " "  << d.Current_UPspin(h, 0, 0, 1) - d.Current_DNspin(h, 0, 0, 1) << endl;//only spin current with equal pop
            out3 << gamma << " " << phi2/Pi << " "  << d.Current_UPspin(h, 0, 0, 2) - d.Current_DNspin(h, 0, 0, 2) << endl;
            
            h.ClearHopTriplet();
            
        }
        
    }
    //end function
}

void EVProbeGammaVInt(ofstream &output1,  Hamiltonian<complex<double>> &h, Lanczos_Diag<complex<double>> &d,  int L, double phi1, double link)
{
    double dg = 0.1;
    double Gm = 5.0/dg;
    double gamma = 0.0;
    
    double dU = 0.1;
    double Um = 5.0/dU;
    double U;
    
    Eigen::VectorXd E;
    
    for(int u = 0; u <= Um; u++)
    {
        U = u*dU;
        h.Set_Const(1.0, 1.0, U);
        cout << "U: " << U <<endl;
        
        for(int g = 0; g <= Gm; g++)
        {
            gamma = g*dg;
            
            //cout << "gamma: " << gamma << endl;
            
            h.HopMatrix_Build_PeriodicNNNHop(link, gamma, phi1);
            
            if (g == 0)
            {
                h.IntMatrix_Build();
            }
            h.Total_Ham_WSOC();
            
            E = d.FullDiagonalization(h);
            
            output1 << U << " " << g << " " << E(0) <<" " << E(1) << " " << E(2)<< " " << E(3) <<" " << E(4) << " " << E(5)  << endl;
            
            h.ClearHopTriplet();
            
        }
        
        h.ClearInteractTriplet();
    }
}

void EVProbePhaseVInt(ofstream &output1, ofstream &output2,  Hamiltonian<complex<double>> &h, Lanczos_Diag<complex<double>> &d, ObservableFunctions<complex<double>> &O, double t1, double t2, int L,  double link)
{
    long double Phi_max = (5.)/6.;//1.5]''
    long double Phi_init = 2./3.;
    double dp = .0001;
    //double dp1 = .0001;
    //double dp2 = .00001;
    double P_it = (Phi_max-Phi_init)/dp;
    //double P_it1 = (Phi_max-Phi_init)/dp1;
    //double P_it2 = (Phi_max-Phi_init)/dp2;
    
    double phi2 = Phi_init*Pi;
    Eigen::VectorXd E;
    
    double dU = 0.1;
    double Um = 4.0/dU;
    double U;
    
    for(int u = 1; u < Um; u++)
    {
        U = u*dU;
        h.Set_Const(t1, t2, U);
        cout << "U: " << U <<endl;
        double Area = 0.0;
        double AreaL = 0.0;
        double AreaR = 0.0;
        long double Jtemp = 0.0;
        complex<long double> var_temp = 0.0;
        complex<long double> var = 0.0;
        long double J = 0.0;
        
        double phi2 = Phi_init*Pi;
        
        
        
        for(int p = 0; p <= P_it+1; p++)
        {
            cout << "phi: " << phi2/Pi << endl;
            
            h.HopMatrix_Build_PeriodicBfield(link, phi2);
            
            if(p == 0)
            {
                h.IntMatrix_Build();
            }
            h.Total_Ham_WSOC();
            
            //Eigenvalues
            E = d.FullDiagonalization(h);
            //output1 << phi2/Pi << " " << E.transpose() << endl;
            
            J = O.Current_UPspin(h, d, 0, 0, phi2, 0, 1);
            complex<long double> var = sqrt(O.CurrentVariance(h, d, 0, 0, 0, 1));
            
            //            double Jup_NN_2 = O.Current_UPspin(h, d, 0, 0, phi2, 1, 2);
            //            complex<double> var_NN_2 = sqrt(O.CurrentVariance(h, d, 0, 1, 2));
            //
            //            double Jup_NNN = O.Current_UPspin(h, d, 0, 0, phi2, 0, 2);//TT is 0 for full current, 1 for uncertainty relation
            //            complex<double> var_NNN = sqrt(O.CurrentVariance(h, d, 0, 0, 2));
            //
            //            complex<double> Uncty_L = var_NN_1*var_NN_1*var_NN_2*var_NN_2;//left side of uncertainty eqn
            //            double Uncty_R = (1./4.)*Jup_NNN*Jup_NNN;//right side of uncertainty eqn
            
            
            
            //output2 << phi2/Pi << " " << Jup_NN_1 << " " << var_NN_1.real() << endl;
            output1 << J << " " << var.real() << endl;
            //output3 << phi2/Pi << " " << Uncty_L.real() << " " << Uncty_R << endl;
            //output3 << phi2/Pi << " " << Jup_NNN << " " << var_NNN.real() << endl;
            
            if(p > 0)
            {
                if(Jtemp < J)
                {
                    double dJ = J - Jtemp;
                    //cout << "dJ 1st loop: " << dJ << endl;
                    
                    //Atop += O.Area(dJ, var.real());
                    AreaL += O.Trapezoidal(dJ, var.real(), var_temp.real());
                }
                else{
                    double dJ = Jtemp-J;
                    //cout << "dJ 2nd loop: " << dJ << endl;
                    
                    //Abtm += O.Area(dJ, var_temp.real());
                    AreaR += O.Trapezoidal(dJ, var.real(), var_temp.real());
                }
            }
            //cout << "Aleft: " << AreaL <<endl;
            //cout << "Aright: " << AreaR << endl;
            
            h.ClearHopTriplet();
            
            phi2 += dp*Pi;
            
            //cout << "J: " << J << " JTemp: " << Jtemp  << endl;
            //cout << "var: " << var << " var temp: " << var_temp  << endl;
            
            Jtemp = J;
            
            
            var_temp = var;
            
        }
        
        cout << "Aleft total: " << AreaL << endl;
        cout << "Aright total: " << AreaR << endl;
        Area = 2.*(AreaR - AreaL);
        
        output1 << endl;
        output2 << U << " " << Area << endl;
        cout << U << " " << Area << endl;
        
        h.ClearInteractTriplet();
        
    }
    
    
}

void HysteresisArea(ofstream &output1, ofstream &output2,  Hamiltonian<complex<double>> &h, Lanczos_Diag<complex<double>> &d,ObservableFunctions<complex<double>> &O, double U, double t1, double t2, int L,  double link)
{
    long double Phi_init = 2./3.;
    double phi = Phi_init*Pi;
    double dpL = 0.0001;
    double dpR = 0.00001;
    long double Phi_max = ((5.)/6.+dpR)*Pi;
    
    int maxit1 = 5000;
    int maxit2 = 5000;
    
//    double Area = 0.0;
//    double AreaL = 0.0;
//    double AreaR = 0.0;
//    long double Jtemp = 0.0;
    long double Jcond = -1.0;
//    complex<long double> var_temp = 0.0;
//    complex<long double> var = 0.0;
    //long double J = 0.0;
    
    int it = 0;
    
    Eigen::VectorXd E;
    
    double dU = 0.1;
    double Um = 2.0/dU;
    U = 0.0;
    
    for(int u = 1; u < Um; u++)
    {
        U = u*dU;
        h.Set_Const(t1, t2, U);
        cout << "U: " << U <<endl;
        double Area = 0.0;
        double AreaL = 0.0;
        double AreaR = 0.0;
        long double Jtemp = 0.0;
        complex<long double> var_temp = 0.0;
        complex<long double> var = 0.0;
        long double J = 0.0;
        long double Jcond = -1.0;
        
        phi = Phi_init*Pi;
    
    do{
        if(it > 1)
        {
            Jcond = Jtemp;
        }
        
        cout << "phi: " << phi/Pi << endl;
        
        h.HopMatrix_Build_PeriodicBfield(link, phi);
        
        if(it == 0)
        {
            h.IntMatrix_Build();
        }
        h.Total_Ham_WSOC();
        
        //Eigenvalues
        E = d.FullDiagonalization(h);
        
        //Current
        J = O.Current_UPspin(h, d, 0, 0, phi, 0, 1);
        var = sqrt(O.CurrentVariance(h, d, 0, 0, 0, 1));
        
        //output Hysteresis Curve
        output1 << J << " " << var.real() << endl;
        
        if(it > 0)
        {
            
            double dJ = J - Jtemp;
            AreaL += O.Trapezoidal(dJ, var.real(), var_temp.real());
            cout << "AreaL: " << AreaL << endl;
        }
        
        h.ClearHopTriplet();
        
        Jtemp = J;
        var_temp = var;
        
        phi += dpL*Pi;
        it++;
        
        cout << "J: " << J << " Jcond: " << Jcond << endl;
        
    }while(Jcond < J && it < (maxit1+1));//this is failing because line JT = J in do
    
    cout << "End loop one and begin loop two\n";
    it = 0;
    output1 << endl;
    
    do{
        cout << "phi: " << phi/Pi << endl;
        
        h.HopMatrix_Build_PeriodicBfield(link, phi);
        //no need to build interaction in this loop process
        h.Total_Ham_WSOC();
        
        //Eigenvalues
        E = d.FullDiagonalization(h);
        
        //Current
        J = O.Current_UPspin(h, d, 0, 0, phi, 0, 1);
        complex<long double> var = sqrt(O.CurrentVariance(h, d, 0, 0, 0, 1));
        
        //output Hysteresis Curve
        output1 << J << " " << var.real() << endl;
        
//        if(it > 0)
//        {
            double dJ = Jtemp-J;
            AreaR += O.Trapezoidal(dJ, var.real(), var_temp.real());
            cout << "AreaR: " << AreaR << endl;
        //}
        
        h.ClearHopTriplet();
        
        Jtemp = J;
        var_temp = var;
        
        cout << "it: " << it << endl;
        it++;
        phi += dpR*Pi;
        //cout << "J 2nd loop: " << Jtemp << endl;
        
    }while( phi <= Phi_max && it < (maxit2+1));//(J >= 0)phi/Pi <= Phi_max
    cout << "AR: " << AreaR << " AL: " << AreaL << endl;
    Area = 2.*(AreaR - AreaL);
    
    output2 << U << " " << Area <<endl;
        output1 << endl;
    }
}

void EVProbePhaseVNNNHop(ofstream &output1, ofstream &output2,  Hamiltonian<complex<double>> &h, Lanczos_Diag<complex<double>> &d,ObservableFunctions<complex<double>> &O, double t1, double t2, double U,  int L,  double link)
{
    long double Phi_max = (5.)/6.;//1.5]''
    long double Phi_init = 2./3.;
    double dp = .0001;
    //double dp1 = .0001;
    //double dp2 = .00001;
    double P_it = (Phi_max-Phi_init)/dp;
    //double P_it1 = (Phi_max-Phi_init)/dp1;
    //double P_it2 = (Phi_max-Phi_init)/dp2;
    
    double phi2 = Phi_init*Pi;
    Eigen::VectorXd E;
    
    double Hoprat = 0.0;
    double HM = 5.0;
    double dH = 0.1;
    double HTot = HM/dH;
    
    
}

void ZeemanObservables(ofstream &output1, ofstream &output2, Hamiltonian<complex<double>> &h, Lanczos_Diag<complex<double>> &d, double gamma)
{
    double B;
    double db = .01;
    double Bmax = 2.;
    int Btot = Bmax/db;
    
    double Phi_Z;
    
    
    for(int b = 0; b <= Btot; b++)
    {
        B = b*db;
        Phi_Z = (sqrt(3.)/4.)*B;
        
        h.HopMatrix_Build_Zeeman_WSOC(6, gamma, Pi/2., Phi_Z, B);
        
        if(b == 0)
        {
            h.IntMatrix_Build();
        }
        
        h.Total_Ham_Zeeman();
        
        VectorXd E = d.FullDiagonalization(h);
        output1 << B << " " << E(0) << endl;
        
        h.ClearHopTriplet();
        
    }
    
}

void NNNRealHop_Obs(ofstream &output1, ofstream &output2, ofstream &output3, Hamiltonian<complex<double>> &h, Lanczos_Diag<complex<double>> &d, double gamma, int L, double phi)
{
    double J_NNN;
    double dJ = .01;
    double Jmax = 3.0;
    int Jtot = Jmax/dJ;
    
    
    for(int j = 0; j < Jtot; j++)
    {
        J_NNN = j*dJ;
        cout <<"JNNN: " << J_NNN << endl;
        
        h.HopMatrix_Build_NNNHop_WReal(6, J_NNN, gamma, phi);
        if(j == 0)
        {
            h.IntMatrix_Build();
        }
        h.Total_Ham_WSOC();
        
        //Eigenvalues
        VectorXd E = d.FullDiagonalization(h);
        output1 << J_NNN << " "  << E.transpose() << endl;
        
        //Current
        output2 << phi/Pi << " " << d.Current_UPspin(h, 0, 0, 1) + d.Current_DNspin(h, 0, 0, 1) << " " << d.Current_UPspin(h, 0, 0, 2)+d.Current_DNspin(h, 0, 0, 2) << endl;
        output3 << phi/Pi << " " << d.Current_UPspin(h, 0, 0, 1) - d.Current_DNspin(h, 0, 0, 1) << " " << d.Current_UPspin(h, 0, 0, 2)-d.Current_DNspin(h, 0, 0, 2) << endl;
        
        h.ClearHopTriplet();
    }
}

void TF_ForceCalc(ofstream &outForce,  Hamiltonian<complex<double>> &h, Lanczos_Diag<complex<double>> &d, int link, double Phi )
{
    
    //ifstream DataIn("ThomasFermi_VSOC_VNNN_rat_ovrR_090717.dat");//1NNN Link Thomas Fermi
    ifstream DataIn("ColoumbPotential_Z92_VSOC_ovrActdist_rb-lo2_121217.dat");//6NNN Link Coloumb
    assert(DataIn.is_open());
    if (!DataIn.is_open())
    {
        cout<<"NOT OPEN"<<endl;
        exit (1);
    }
    
    int row = 9440;//786;//122;
    int col = 2;//4;
    
    double E0_temp = 0.0;
    double E0;
    double dx = 0.001;
    
    double tBenz = 4.069e-19;//J of hopping
    double Angs = 0.529;//Angstrom in 1 ab
    double NewtonConversion = tBenz/(Angs*1e-10);//J/meter
    
    MatrixXd DataArray;
    VectorXd E;
    
    DataArray.resize(row,col);
    
    for(int i = 0; i < row; i++)
    {
        for(int j = 0; j < col; j++)
        {
            DataIn >> DataArray(i,j);
        }
    }
    //col 0 is R
    //col 1 is gamma
    //cout << DataArray << endl;
    
    for(int r = 0; r < row; r++)
    {
        
        double gamma = DataArray(r,1);
        double x = DataArray(r,0);//in Bohr radii Convert to Angstrom in output
        cout << "r: " << r << endl;
        
        h.HopMatrix_Build_PeriodicNNNHop(link, gamma, Phi);
        
        if (r == 0)
        {
            h.IntMatrix_Build();
        }
        h.Total_Ham_WSOC();
        
        E = d.FullDiagonalization(h);
        E0 = E(0);
        
        if(r != 0)
        {
            double dE = E0 - E0_temp;
            outForce << x*Angs << " "  << E0 << " " << -1.*NewtonConversion*dE/dx << endl;
        }
        
        E0_temp = E0;
        
        h.ClearHopTriplet();
    }
    
}

template<typename Tnum>
void LocalCurrent(Hamiltonian<Tnum> &h, Lanczos_Diag<Tnum> &d, int s1, int s2, double Phi)
{
    d.TotalCurrents(h, s1, s2, Phi);
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
template class Hamiltonian<complex<double>>;