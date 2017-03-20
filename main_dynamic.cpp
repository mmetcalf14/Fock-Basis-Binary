//
//  main_dynamic.cpp
//  ED Hubbard Hamiltonian
//
//  Created by mekena McGrew on 11/4/15.
//  Copyright © 2015 Mekena Metcalf. All rights reserved.
//


#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
//#include </usr/local/include/Eigen/Sparse>
//#include </usr/local/include/Eigen/Eigen>
//#include </usr/local/include/Eigen/Core>
#include "Hamiltonian.h"
#include "Lanczos.h"


using namespace std;
using namespace Eigen;
void Write_Density(ofstream &fout, vector<double> &n_up, vector<double> &n_dn, int L );
void DensityDiff(ofstream &fout, vector<double> &n_up, vector<double> &n_dn, int L );
void DensityInTime(ofstream &fout, vector<double> &n, int L, double Phi );
void Density(ofstream &fout, vector<double> &n, int L);
void EdgeDensity(ofstream &fout, vector<double> &n, double time);//can I move these functions to a different cpp file?


void Harmonic_Trap(vector<double> &HT, int L, double y);

template<typename Tnum>
void Fidelity(ofstream &output, Hamiltonian<Tnum> &h, Lanczos_Diag<Tnum> &d , double U, double J1, double J2, double Um, double Jm, int L);

template<typename Tnum>
void U_Fid(ofstream &output, Hamiltonian<Tnum> &h, Lanczos_Diag<Tnum> &d , double U, double J1, double J2, double Um, int L);

template<typename Tnum>
void U_Fid_V2(ofstream &output, Hamiltonian<Tnum> &h, Lanczos_Diag<Tnum> &d , double U, double J1, double J2, double Um, int L);

template<typename Tnum>
void Fidelity_HT(ofstream &output, Hamiltonian<Tnum> &h, Lanczos_Diag<Tnum> &d , double U, double J1, double J2, double Um, double Ym, vector<double> HT, int L);

template<typename Tnum>
void PeierlsTD(ofstream &out_up, ofstream &out_dn, Hamiltonian<Tnum> &h, Lanczos_Diag<Tnum> &d, int T_it, double dt, double tp, int L);

template<typename Tnum>
void QPump_TD(ofstream &out_up, ofstream &out_dn, Hamiltonian<Tnum> &ham, Lanczos_Diag<Tnum> &diag, const int T, const double dt, const double h_0, const double J_0, const double d_0, double U, int L);

template<typename Tnum>
void CorrTD(ofstream &output, ofstream &out, Hamiltonian<Tnum> &h, Lanczos_Diag<Tnum> &d, int cut, int T_it, double dt, int L);

void EVProbeGamma(ofstream &output, Hamiltonian<complex<double>> &h, Lanczos_Diag<complex<double>> &d, int s1, int s2, int L);
void EVProbeInt(ofstream &output, Hamiltonian<complex<double>> &h, Lanczos_Diag<complex<double>> &d, int s1, int s2, int L, double gSO);

template<typename Tnum>
void LocalCurrent(Hamiltonian<Tnum> &h, Lanczos_Diag<Tnum> &d, int s1, int s2, double Phi);

double GdotG( const VectorXd &G1, const VectorXd &G2);
double GdotG( const VectorXcd &Gc1, const VectorXcd &Gc2);



int main(int argc, const char * argv[])
{



    int Nup;
    int Ndown;
    int Nsite;
    double Y = 0.4;
    double Ymax = 1.;
    double t_1 = 1.;
    double t_2;
    double hop_rat;
    double U;
    double Umax = 5.;
    double Jmax = 2.;
    int T_f;
    double dt = .01;
    double t_p;

    const double h0 = 0.5;
    const double d0 = 0.5;
    const double J0 = 1.0;
    double cut;
    double gamma;
    double Phi = Pi/4.;
    int Site1;
    int Site2;



    char output[80];
//    char HTFout[80];
    char IntFidOut[80];
//    char PD[80];
//    char PeierlsDensity_up[100];
//    char PeierlsDensity_down[100];
//    char QPD_up[100];
//    char QPD_dn[100];
//    char CorrOut[100];
    

    Matrix4d Test_Ham;
    Vector4d Test_Lanczos;
    Test_Ham << 0, -1, 0, 0,
    -1,0,-1,0,
    0,-1,0,-1,
    0,0,-1,0;
    Test_Lanczos << 0.5,0.5,0.5,0.5;

    vector<double> Harm_Trap (Nsite, 0.0);

    ifstream ReadFile("ED_ActAro_J1J2_DataInput.cfg");
    assert(ReadFile.is_open());
    if (!ReadFile.is_open())
    {
        cout<<"NOT OPEN"<<endl;
        exit (1);
    }

    ReadFile >> Nup;
    ReadFile >> Ndown;
    ReadFile >> Nsite;
    ReadFile >> hop_rat;
    ReadFile >> U;
    ReadFile >> gamma;
    ReadFile >> Site1;
    ReadFile >> Site2;
    ReadFile >> output;
    ReadFile >> IntFidOut;
    

    cout << Nup << endl;
    cout << Ndown << endl;
    cout << Nsite << endl;
    cout << t_1 << endl;
    //cout << t_2 << endl;
    cout << U << endl;
    cout << output << endl;

    

    int T_tot = T_f/dt;
    t_2 = t_1/hop_rat;
    cout << "Hopping rat: " << t_2 << endl;
    

    ofstream fout(output);
    assert(fout.is_open());
    fout.setf(ios::scientific);
    fout.precision(11);
//
    ofstream FidOut(IntFidOut);
    assert(FidOut.is_open());
    FidOut.setf(ios::scientific);
    FidOut.precision(11);
//
//    ofstream HTOut(HTFout);
//    assert(HTOut.is_open());
//    HTOut.setf(ios::scientific);
//    HTOut.precision(11);
//
//    ofstream PDout(PD);
//    assert(PDout.is_open());
//    PDout.setf(ios::scientific);
//    PDout.precision(11);

    
//    ofstream QPout_up(QPD_up);
//    assert(QPout_up.is_open());
//    QPout_up.setf(ios::scientific);
//    QPout_up.precision(11);
//    
//    ofstream QPout_dn(QPD_dn);
//    assert(QPout_dn.is_open());
//    QPout_dn.setf(ios::scientific);
//    QPout_dn.precision(11);
    

    

//    ofstream PDout_up(PeierlsDensity_up);
//    assert(PDout_up.is_open());
//    PDout_up.setf(ios::scientific);
//    PDout_up.precision(11);
//
//    ofstream PDout_dn(PeierlsDensity_down);
//    assert(PDout_dn.is_open());
//    PDout_dn.setf(ios::scientific);
//    PDout_dn.precision(11);
//    
//    ofstream OutCorr(CorrOut);
//    assert(OutCorr.is_open());
//    OutCorr.setf(ios::scientific);
//    OutCorr.precision(11);

   //creating vector with harmoinc trap values per site
    //Harmonic_Trap(Harm_Trap, Nsite, Y);

    //Build basis and pass to Hamiltonian class through inheritance
    //COMPLEX
    Hamiltonian<complex<double> > ham(Nsite, Nup, Ndown);
    Lanczos_Diag<complex<double> > Diag(ham);
    
    //DOUBLE
    //Hamiltonian<double> ham(Nsite, Nup, Ndown);
//    Lanczos_Diag<double> Diag(ham);


    //ham.GetHarmTrap(Harm_Trap);
    //Fidelity_HT(HTOut, ham, Diag, U, t_1, t_2, Umax, Ymax, Harm_Trap, Nsite);
    //Fidelity(FidOut, ham, Diag, U, t_1, t_2, Umax, Jmax, Nsite);
    //U_Fid(FidOut, ham, Diag, U, t_1, t_2, Umax, Nsite);
    //U_Fid_V2(FidOut, ham, Diag, U, t_1, t_2, Umax, Nsite);

    //set hopping and interaction coefficients
    ham.Set_Const(t_1, t_2, U);//U=0 until |G> is found for t=0



    //set hamiltonian from triplets
    //ham.HopMatrix_Build();
    //ham.HopMatrix_Build_Periodic();
   //ham.HopMatrix_Build_PeriodicWithSOC(Site1, Site2, gamma);
    
    ham.GetPhi(Phi);
    ham.HopMatrix_Build_Peierls();

    
    //build interaction matrix
    ham.IntMatrix_Build();
//
//    //add together all three matrices for total Ham
    //ham.Total_Ham();
    ham.Total_Ham_WSOC();
//
//    //create object for diag class
//
//    //Diag.Lanczos_TestM(Test_Ham, Test_Lanczos);
//
//    //set Lanczos vector dimensions
//    cout << "Setting LA Dim \n";
    Diag.Set_Mat_Dim_LA(ham);
//
    //cout << "Diagonalizing \n";
//    //Diagonalization of t=0 Hamiltonian
    //Diag.Diagonalize(ham);
    //Diag.CHECK();
    VectorXd V = Diag.FullDiagonalization(ham);
    //cout << "EigenVals: " << V << endl;
//
//
//    //convert |G> from Fock basis to onsite basis
//    //seperate |G> states for nup and ndn
    //cout << "Getting Density\n";
    //Diag.Density(ham);//before interaction turned on
    //Write_Density(fout, Diag.n_up, Diag.n_dn, Nsite);
    
    
    //Correlation Functions
    //Use Cut 0 to
    //double C;
   // C = Diag.DensityWCorr(ham, 0);
    
	//Diag.SpinCorr(ham,fout, 0, 0);//spin correlation
    
   // cout << "Correlation:  " << C << endl;
    
    
    //Current
    Diag.TotalCurrents(ham, 3, 4);

        //Triplets removed to redo Interaction matrix after quenching
        // and all non-zero elemenst of Total Ham and Ham_U are set to zero
        //ham.ClearInteractTriplet();


        //Interactions turned on after intial ground state found
        //if U=0 there is no need to build interaction ham until after ground state is found
//        ham.QuenchU(U);
//        ham.IntMatrix_Build();
//        ham.Total_Ham();

        //Lanczos_Diag<complex<double> > Diag_Comp(ham);

        //Time Evolve

    
    //QPump_TD(PDout_up, PDout, ham, Diag, T_f, dt, h0, J0, d0, U, Nsite);//QPout_up, QPout_dn
    
    //PeierlsTD(PDout_up, PDout, ham, Diag, T_tot, dt, t_p, Nsite);
    
    //CorrTD(OutCorr, fout, ham, Diag, cut, T_tot, dt, Nsite);
    
    //EVProbeGamma(fout, ham, Diag, Site1, Site2, Nsite);
    //EVProbeInt(fout, ham, Diag, Site1, Site2, Nsite, gamma);
        //int NN = T_tot/10;

    FidOut.close();
//    HTOut.close();
//    PDout.close();
//    PDout_dn.close();
//    PDout_up.close();
//    OutCorr.close();

//    QPout_up.close();
//    QPout_dn.close();
    
    fout.close();
    cout << "Code is Done! \n";

    return 0;
}

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

void EVProbeGamma(ofstream &output, Hamiltonian<complex<double>> &h, Lanczos_Diag<complex<double>> &d, int s1, int s2, int L)
{
 
    double dg = 0.1;
    double Gm = 2.0/dg;
    double gamma = 0.0;
    Eigen::VectorXd E;
    cout << "TB: " << h.Tot_base << endl;
    
    for(int g = 0; g <= Gm; g++)
    {
        gamma = g*dg;
        cout << "gamma: " << gamma << endl;
        h.HopMatrix_Build_PeriodicWithSOC(s1, s2, gamma);
        if (g == 0)
        {
            h.IntMatrix_Build();
        }
        h.Total_Ham_WSOC();
        
        E = d.FullDiagonalization(h);
        for(int i = 0; i < h.Tot_base; i++)
        {
            if(i == 0)
            {
                output << gamma << " " << E(i) << " ";
            }
            else{
                output << E(i) << " ";
            }
        }
        output << endl;
        
        h.ClearHopTriplet();
        
    }
    
}

void EVProbeInt(ofstream &output, Hamiltonian<complex<double>> &h, Lanczos_Diag<complex<double>> &d, int s1, int s2, int L, double gSO)
{
    double dU = 0.01;
    double Um = 3.0/dU;
    double U;
    cout << "g: " <<  gSO << endl;

    Eigen::VectorXd E;
    cout << "TB: " << h.Tot_base << endl;
    
    for(int u = 0; u <= Um; u++)
    {
        U = u*dU;
        h.Set_Const(1.0, 1.0, U);
        cout << "U: " << U <<endl;
        
        if (u == 0)
        {
            
            h.HopMatrix_Build_PeriodicWithSOC(s1, s2, gSO);
        }
        h.IntMatrix_Build();
        h.Total_Ham_WSOC();
        
        E = d.FullDiagonalization(h);
        for(int i = 0; i < h.Tot_base; i++)
        {
            if(i == 0)
            {
                output << U << " " << E(i) << " ";
            }
            else{
                output << E(i) << " ";
            }
        }
        output << endl;
        
        h.ClearInteractTriplet();
        
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


