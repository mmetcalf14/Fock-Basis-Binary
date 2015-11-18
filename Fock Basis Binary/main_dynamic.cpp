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
void Harmonic_Trap(vector<double> &HT, int L, double y);

template<typename Tnum>
void Fidelity(ofstream &output, Hamiltonian<Tnum> &h, Lanczos_Diag<Tnum> &d , double U, double J1, double J2, double Um, double Jm, int L);
template<typename Tnum>
void Fidelity_HT(ofstream &output, Hamiltonian<Tnum> &h, Lanczos_Diag<Tnum> &d , double U, double J1, double J2, double Um, double Ym, vector<double> HT, int L);
double GdotG( const VectorXd &G1, const VectorXd &G2);
double GdotG( const VectorXcd &Gc1, const VectorXcd &Gc2);

int main(int argc, const char * argv[])
{
    

    
    int Nup;
    int Ndown;
    int Nsite;
    double Y = 0.;
    double Ymax = 1.;
    double t_1;
    double t_2;
    double U;
    double Umax = 10;
    double Jmax = 2;
    int T_f;
    double dt = 1.;
    char output[60];
    Matrix4d Test_Ham;
    Vector4d Test_Lanczos;
    Test_Ham << 0, -1, 0, 0,
    -1,0,-1,0,
    0,-1,0,-1,
    0,0,-1,0;
    Test_Lanczos << 0.5,0.5,0.5,0.5;
    
    vector<double> Harm_Trap (Nsite, 0.0);

    ifstream ReadFile("ED_J1J2_DataInput.cfg");
    assert(ReadFile.is_open());
    if (!ReadFile.is_open())
    {
        cout<<"NOT OPEN"<<endl;
        exit (1);
    }
    
    ReadFile >> Nup;
    ReadFile >> Ndown;
    ReadFile >> Nsite;
    ReadFile >> t_1;
    ReadFile >> t_2;
    ReadFile >> U;
    ReadFile >> T_f;
    ReadFile >> output;
    
    cout << Nup << endl;
    cout << Ndown << endl;
    cout << Nsite << endl;
    cout << t_1 << endl;
    cout << t_2 << endl;
    cout << U << endl;
    cout << output << endl;
    
    int T_tot = T_f/dt;
    
    
    ofstream fout(output);
    assert(fout.is_open());
    
    fout.setf(ios::scientific);
    fout.precision(11);
    
    ofstream FidOut("ED_Nu2_Nd3_L5_Fidelity_110515.dat");
    assert(FidOut.is_open());
    ofstream HTOut("ED_Nu2_Nd3_L5_Fidelity_HT_110515.dat");
    assert(HTOut.is_open());
    
    FidOut.setf(ios::scientific);
    FidOut.precision(11);
    HTOut.setf(ios::scientific);
    HTOut.precision(11);
    
   //creating vector with harmoinc trap values per site
    Harmonic_Trap(Harm_Trap, Nsite, Y);
    
    //Build basis and pass to Hamiltonian class through inheritance
    Hamiltonian<double> ham(Nsite, Nup, Ndown);
    Lanczos_Diag<double> Diag(ham);

    
    //ham.GetHarmTrap(Harm_Trap);
    Fidelity_HT(HTOut, ham, Diag, U, t_1, t_2, Umax, Ymax, Harm_Trap, Nsite);
    //Fidelity(FidOut, ham, Diag, U, t_1, t_2, Umax, Jmax, Nsite);

//    //set hopping and interaction coefficients
//    ham.Set_Const(t_1, t_2, U);//U=0 until |G> is found for t=0
//    
//    
//    
//    //set hamiltonian from triplets
//    ham.HopMatrix_Build();
//    
//    
//    //build interaction matrix
//    ham.IntMatrix_Build();
//    
//    //add together all three matrices for total Ham
//    ham.Total_Ham();
//    
//    //create object for diag class
//    
//    //Diag.Lanczos_TestM(Test_Ham, Test_Lanczos);
//    
//    //set Lanczos vector dimensions
//    //cout << "Setting LA Dim \n";
//    Diag.Set_Mat_Dim_LA(ham);
//    
//    //cout << "Diagonalizing \n";
//    //Diagonalization of t=0 Hamiltonian
//    Diag.Diagonalize(ham);
//    
    
//    //convert |G> from Fock basis to onsite basis
//    //seperate |G> states for nup and ndn
//    //cout << "Getting Density\n";
//    Diag.Density(ham);//before interaction turned on
//    Write_Density(fout, Diag.n_up, Diag.n_dn, Nsite);
//    
    //    //Triplets removed to redo Interaction matrix after quenching
    //    // and all non-zero elemenst of Total Ham and Ham_U are set to zero
    //    ham.ClearTriplet();
    //
    //
    //    //Interactions turned on after intial ground state found
    //    //if U=0 there is no need to build interaction ham until after ground state is found
    //    ham.QuenchU(U);
    //    ham.IntMatrix_Build();
    //    ham.Total_Ham();
    //
    //    //Lanczos_Diag<complex<double> > Diag_Comp(ham);
    //
    //    //Time Evolve
    //    Diag.TimeEvoCoeff(dt);
    //   // Diag.Dynamics(ham, ham);
    //
    //    int NN = T_tot/10;
    //    int Nflag = 0;
    //    for(int t = 0; t < T_tot; t++)
    //    {
    //        //cout << "iteration: "<< t << endl;
    //        Diag.Dynamics(ham);
    //
    //        if(Nflag == NN-1)
    //        {
    //           Diag.Density(ham);
    //           Write_Density(fout, Diag.n_up, Diag.n_dn, Nsite);
    //            Nflag = 0;
    //        }
    //
    //        Nflag++;
    //    }
    //
    FidOut.close();
    HTOut.close();
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
    cout << endl;
    
    
}

void Harmonic_Trap(vector<double> &HT, int L, double y)
{
    int Ro = (L +1)/2;
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

                g = (1-F)/(L*dU);//this is wrong
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
void Fidelity_HT(ofstream &output, Hamiltonian<Tnum> &h, Lanczos_Diag<Tnum> &d , double U, double J1, double J2, double Um, double Ym, vector<double> HT, int L)
{
    typedef Eigen::Matrix<Tnum, Eigen::Dynamic, 1> VectorType;
    
    double dU = 1.;
    double dw = .01;
    double y = 0;
    int Wtot = Ym/dw;
    int Utot = Um/dU;
    double g;
    double F;
    VectorType Gstate;
    VectorType Gstate_Temp;
    // U = 0;
    

        for(int w = 0; w <= Wtot; w++)
        {
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
                
                g = (1-F)/(L*dw);//this is wrong
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


