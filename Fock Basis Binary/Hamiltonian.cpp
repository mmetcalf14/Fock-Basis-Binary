//
//  Hamiltonian.cpp
//  ED Hubbard Hamiltonian
//
//  Created by mekena McGrew on 11/4/15.
//  Copyright © 2015 Mekena Metcalf. All rights reserved.
//

#include <stdio.h>
#include <cassert>
#include "Hamiltonian.h"
using namespace std;
//void Hamiltonian::Set_tbar(int _tbar)
//{ tbar = _tbar;}

template<typename Tnum>
void Hamiltonian<Tnum>::Set_Const(Tnum t_1, Tnum t_2, Tnum _U)
{
    J1 = t_1;
    J2 = t_2;
    U = _U;

}

template<typename Tnum>
void Hamiltonian<Tnum>::QuenchU(Tnum _Uquench)
{
    U = _Uquench;
}

template<typename Tnum>
void Hamiltonian<Tnum>::GetPhi(double _Phi_t)
{
    Phi_t = _Phi_t;
    //cout <<"Phi: " << Phi_t << endl;
}

template<typename Tnum>
void Hamiltonian<Tnum>::GetOnsite(Tnum _h)
{
    h = _h;

}

template<typename Tnum>
void Hamiltonian<Tnum>::GetHarmTrap(std::vector<double> HT)
{
    Harm_Trap = HT;
}

template<typename Tnum>
void Hamiltonian<Tnum>::BuildHopHam(int species, size_t count, size_t count_opp,
                                    vector<size_t> basis, vector<size_t> index, SpMat &HopHam, std::vector<double> HT)
{
    std::vector<Tp> TL;

    for(size_t bs = 0; bs < count; bs++)
    {
        size_t p_bas = basis[bs];

        size_t p_ind = index[p_bas];
        //cout << "We are acting on basis, " << p_bas << " with index, "<< p_ind << endl;
        for(size_t i = 0; i < (L-1); i++)
        {
            if ( MY_bittest(p_bas, i) && !(MY_bittest(p_bas, i+1)) )//next to nearest neigbor i+2
                //if i+1 is occupied and going to to i+2 must get minus sign
            {
                size_t l_bas = MY_bitset(MY_bitclr(p_bas,i),i+1);
                size_t l_ind = index[l_bas];
                assert( l_bas != p_bas );

                Tnum val;

                if( count_opp )
                {
                    for(size_t k = 0; k < count_opp; k++)
                    {
                        int r,s;
                        if ( species == 0 ){
                            r = TotalIndex(p_ind, k);//(k*count) + p_ind;
                            s = TotalIndex(l_ind, k);//(k*count) + l_ind;
                        }else if ( species == 1 ){
                            r = TotalIndex(k, p_ind);
                            s = TotalIndex(k, l_ind);
                        }else{
                            //cout << "More than 2 species fermion!!" << endl;
                        }
                        //cout << "pbas: "<<p_bas<<" lbas: "<<l_bas<< " " << r << " " << s << endl;

                        if((i%2) == 0)//even number sites have J1 hop to nn (A)
                        {
                            val = -J1 ;
                        }
                        else //odd number sites have J2 (B)
                        {
                            val = -J2;
                        }


                        TL.push_back(Tp(r,s,val));
                        TL.push_back(Tp(s,r,val));

                        //I need to reference the total index


                    }
                }
                else
                {
                    if((i%2) == 0)
                    {
                        val = -J1 ;
                    }
                    else
                    {
                        val = -J2;
                    }

                    TL.push_back(Tp(p_ind,l_ind, val ));
                    TL.push_back(Tp(l_ind,p_ind, val ));

                }
            }
        }
    }

   if(HT.size() > 0)
   {
    for(size_t bs = 0; bs < count; bs++)//harmonic trap
    {
        size_t bas = basis[bs];

        size_t ind = index[bas];
        for(size_t i = 0; i < L; i++)
        {
            if(MY_bittest(bas, i))
            {
                Tnum val = HT[i];
                for(size_t k = 0; k < count_opp; k++)
                {
                    int q;
                    if(species == 0)
                    {
                        q = TotalIndex(ind, k);
                    }
                    if(species == 1)
                    {
                        q = TotalIndex(k,ind);
                    }
                    TL.push_back(Tp(q,q,val));

                }
            }
        }
    }
   }

    HopHam.setFromTriplets(TL.begin(), TL.end());
    cout << "Hop Ham set \n";

}

template<>
void Hamiltonian<complex<double> >::BuildHopHam_Peierls(int species, size_t count, size_t count_opp,
                                    vector<size_t> basis, vector<size_t> index, SpMat &HopHam, std::vector<double> HT)
{
    std::vector<Tp> TL;
    complex<double> I;
    I.real(0.0);
    I.imag(1.0);
    

    for(size_t bs = 0; bs < count; bs++)
    {
        size_t p_bas = basis[bs];

        size_t p_ind = index[p_bas];
        //cout << "We are acting on basis, " << p_bas << " with index, "<< p_ind << endl;
        for(size_t i = 0; i < (L-1); i++)
        {
            if ( MY_bittest(p_bas, i) && !(MY_bittest(p_bas, i+1)) )//next to nearest neigbor i+2
                //if i+1 is occupied and going to to i+2 must get minus sign
            {
                size_t l_bas = MY_bitset(MY_bitclr(p_bas,i),i+1);
                size_t l_ind = index[l_bas];
                assert( l_bas != p_bas );

                complex<double> val;

                if( count_opp )
                {
                    for(size_t k = 0; k < count_opp; k++)
                    {
                        int r,s;
                        if ( species == 0 ){
                            r = TotalIndex(p_ind, k);//(k*count) + p_ind;
                            s = TotalIndex(l_ind, k);//(k*count) + l_ind;
                        }else if ( species == 1 ){
                            r = TotalIndex(k, p_ind);
                            s = TotalIndex(k, l_ind);
                        }else{
                            //cout << "More than 2 species fermion!!" << endl;
                        }
                        //cout << "pbas: "<<p_bas<<" lbas: "<<l_bas<< r << " " << s << endl;

                        if((i%2) == 0)//even number sites have J1 hop to nn (A)
                        {
                            val = -J1*exp(I*Phi_t) ;
                        }
                        else //odd number sites have J2 (B)
                        {
                            val = -J2*exp(I*Phi_t);
                        }

                       // cout << "Setting Triplet\n";
                        TL.push_back(Tp(r,s,val));
                        TL.push_back(Tp(s,r,conj(val)));//take conjugate of exponential for other half of matrix
                        ///cout << "Triplet Set\n";
                        //I need to reference the total index


                    }
                }
                else
                {
                    if((i%2) == 0)
                    {
                        val = -J1*exp(I*Phi_t) ;
                    }
                    else
                    {
                        val = -J2*exp(I*Phi_t);
                    }

                    TL.push_back(Tp(p_ind,l_ind, val ));
                    TL.push_back(Tp(l_ind,p_ind, conj(val) ));

                }
            }
        }
    }

    
    //cout << "exponential: " << exp(I*Phi_t) << endl;
    HopHam.setFromTriplets(TL.begin(), TL.end());
    //cout << "Hop Ham set \n";

}

template<typename Tnum>
void Hamiltonian<Tnum>::BuildHopHam_QPump(int species, size_t count, size_t count_opp, std::vector<size_t> basis, std::vector<size_t> index, SpMat &HopHam)
{
    std::vector<Tp> TL;

    for(size_t bs = 0; bs < count; bs++)
    {
        size_t p_bas = basis[bs];

        size_t p_ind = index[p_bas];
        //cout << "We are acting on basis, " << p_bas << " with index, "<< p_ind << endl;
        for(size_t i = 0; i < (L-1); i++)
        {
            if ( MY_bittest(p_bas, i) && !(MY_bittest(p_bas, i+1)) )//next to nearest neigbor i+2
                //if i+1 is occupied and going to to i+2 must get minus sign
            {
                size_t l_bas = MY_bitset(MY_bitclr(p_bas,i),i+1);
                size_t l_ind = index[l_bas];
                assert( l_bas != p_bas );

                Tnum val;

                if( count_opp )
                {
                    for(size_t k = 0; k < count_opp; k++)
                    {
                        int r,s;
                        if ( species == 0 ){
                            r = TotalIndex(p_ind, k);//(k*count) + p_ind;
                            s = TotalIndex(l_ind, k);//(k*count) + l_ind;
                        }else if ( species == 1 ){
                            r = TotalIndex(k, p_ind);
                            s = TotalIndex(k, l_ind);
                        }else{
                            //cout << "More than 2 species fermion!!" << endl;
                        }
                        //cout << "pbas: "<<p_bas<<" lbas: "<<l_bas<< " " << r << " " << s << endl;

                        if((i%2) == 0)//even number sites have J1 hop to nn (A)
                        {
                            val = -J1 ;
                        }
                        else //odd number sites have J2 (B)
                        {
                            val = -J2;
                        }


                        TL.push_back(Tp(r,s,val));
                        TL.push_back(Tp(s,r,val));

                        //I need to reference the total index


                    }
                }
                else
                {
                    if((i%2) == 0)
                    {
                        val = -J1 ;
                    }
                    else
                    {
                        val = -J2;
                    }

                    TL.push_back(Tp(p_ind,l_ind, val ));
                    TL.push_back(Tp(l_ind,p_ind, val ));

                }
            }
        }
    }

    //Checked that Hamiltonian matches SPM exactly
        for(size_t bs = 0; bs < count; bs++)//harmonic trap
        {
            size_t bas = basis[bs];

            size_t ind = index[bas];
            for(size_t i = 0; i < L; i++)
            {
                Tnum val;
                if(MY_bittest(bas, i))
                {
                   if( (i%2) == 0)
                   {
                       val = -h;
                   }
                   else
                   {
                       val = h;
                   }


                    for(size_t k = 0; k < count_opp; k++)
                    {
                        int q;
                        if(species == 0)
                        {
                            q = TotalIndex(ind, k);
                        }
                        if(species == 1)
                        {
                            q = TotalIndex(k,ind);
                        }
                        TL.push_back(Tp(q,q,val));

                    }
                }
            }
        }

    HopHam.setZero();
    HopHam.setFromTriplets(TL.begin(), TL.end());
    //cout << "Hop Ham set \n";
}

template<typename Tnum>
void Hamiltonian<Tnum>::IntMatrix_Build()
{
    //build interaction index for up spin
    size_t point_up = 0;
    size_t point_dn = 0;
    Eigen::MatrixXd IndexU_up = Eigen::MatrixXd::Zero(L, count_up);
    Eigen::MatrixXd IndexU_dn = Eigen::MatrixXd::Zero(L, count_dn);
    for(size_t i = 0; i < L ; i++)
    {
        point_up = 0;
        for(size_t j = 0; j < count_up; j++)
        {
            size_t j_bas = basis_up[j];
            if((Nup + Ndn) <= L)
            {
                if(MY_bittest(j_bas,i)) //testing if site is occupied
                {
                    //                    ++point_up;
                    //cout << i << " " << j_bas << " " << MY_bittest(j_bas,i) << endl;
                    IndexU_up(i,point_up) = j;//was j+1
                    //std::cout << point_up << std::endl;
                    if(Nup == Ndn)
                    {
                        point_dn = point_up;
                        IndexU_dn(i,point_dn) = j;//was j+1
                        point_dn++; //why is this causing errors?

                    }
                    point_up++;
                }
            }
            else //testing if site is unoccupied, Nup + Ndn > L
            {
                if(!(MY_bittest(j_bas,i))) //testing if site is unoccupied ALG NOT WORKING
                {
                    //                    ++point_up;
                    //cout << i << " " << j_bas << " " << ~MY_bittest(j_bas,i) << endl;
                    IndexU_up(i,point_up) = j;//was j+1
                    if(Nup == Ndn)
                    {
                        point_dn = point_up;
                        IndexU_dn(i,point_dn) = j;//was j+1
                        point_dn++;
                    }
                    point_up++;
                }
            }
        }
    }

    if(Nup != Ndn)
    {
        for(size_t i = 0; i < L ; i++)
        {
            point_dn = 0;
            for(size_t j = 0; j < count_dn; j++)
            {
                size_t j_bas = basis_down[j];
                if((Nup + Ndn) <= L)
                {
                    if(MY_bittest(j_bas,i)) //testing if site is occupied
                    {
                        //                        ++point_dn;
                        IndexU_dn(i,point_dn) = j;//was j+1
                        point_dn++;
                    }
                }
                else //testing if site is unoccupied, Nup + Ndn > L
                {
                    if(!(MY_bittest(j_bas,i))) //testing if site is unoccupied
                    {
                        //                        ++point_dn;
                        IndexU_dn(i,point_dn) = j;//j+1
                        point_dn++;
                    }
                }
            }
        }

    }

    //    std::cout << "The spin up Index matrix is: \n" << IndexU_up << std::endl;
    //    std::cout << "The spin down Index matrix is: \n" << IndexU_dn << std::endl;
    // }
    //
    // void Hamiltonian::BaseInteraction()
    // {
    //cout << "Was Hamiltonian::BaseInteraction" << endl;
    std::vector<Tp> TL_Ubase;
    Tnum g;
    Tnum NNup = Nup;
    Tnum NNdn = Ndn;
    Tnum LL = L;

    if((Nup + Ndn) > L)//
    {
        g = U*(NNup + NNdn -LL);

        for(size_t i = 0; i < Tot_base; i++)
        {
            TL_Ubase.push_back(Tp(i,i, g ));//do I need to have two different Hamiltonians?
        }

    }

    if(Nup == Ndn)
    {
        if( (Nup+Ndn) <= L)
        {
            g = U*NNup;

        }
        else
        {
            g = U*(LL-NNup);//should I do this on top of other build with g?

        }
        for(size_t i = 0; i < count_up; i++)
        {
            size_t k = TotalIndex(i, i);
            // size_t k = (count_up +1)*i - count_up - 1;//added a -1 because for loop starts at 1
            TL_Ubase.push_back(Tp(k,k, g ));//do I need to have two different Hamiltonians?
           // cout << k << " " << k << " " << g << endl;

        }
    }
    //end function algorithm
    // }
    //
    // void Hamiltonian::Build_Interactions()
    // {
    //cout << "Was Hamiltonian::Build_Interactions" << endl;
    //Construct Diagonal elements from Fock space mixed indices
    for(int i = 0; i < L; i++)//The Matrix class in Eigen only allows ints for index
    {
        for( int k = 0; k < point_up; k++)
        {
            for(int l = 0; l < point_dn; l++)//Currently I'm not adding in the if statement because
            {  int r;
                if(Nup == Ndn)
                { //cout << "k: " << k << " l: " << l << endl;
                    if( l != k )
                    {
                        //cout << i << " " << k << " " << l << endl;
                        //cout << IndexU_dn(i,l) << " " << IndexU_up(i,k) << endl;
                        r = (IndexU_dn(i,l) *count_up) + IndexU_up(i,k);//=TotalIndex(IndexU_up(i,k), IndexU_dn(i,l)
                        //cout << "r: " << r << endl;
                        TL_Ubase.push_back(Tp(r,r, U ));//do we need a -1? double check here
                        //cout << r << " " << r << " " << U << endl;
                    }
                }
                else
                {
                    r = (IndexU_dn(i,l)*count_up) + IndexU_up(i,k);
                    TL_Ubase.push_back(Tp(r,r, U ));
                }
            }
        }
    }
    Ham_Interact.setFromTriplets(TL_Ubase.begin(), TL_Ubase.end());
    cout << "No problem building interaction Ham \n";
    //cout << "The total interaction Hamiltonian is: /n"<< Ham_Interact << endl;
    //end function algorithm
}

// void Hamiltonian::Set_Mat_Dim()
// {
//     //std::cout << "Entering dimension alg \n";
//     // Tot_base = count_up*count_dn; //should this be size t or int for matrix dim
//     // HopHam_down.resize(Tot_base, Tot_base);
//     // HopHam_up.resize(Tot_base, Tot_base);
//     // Ham_Interact.resize(Tot_base,Tot_base);
//
//     IndexU_up.resize(L,count_up);
//     IndexU_dn.resize(L,count_dn);
// }

template<typename Tnum>
void Hamiltonian<Tnum>::HopMatrix_Build()
{
    std::cout << "No problem before setting Triplet\n";
    BuildHopHam(0, count_up, count_dn, basis_up, index_up, HopHam_up, Harm_Trap);
    std::cout << "No problem after up spin\n";
    //cout << "HopHam_up: " << HopHam_up << endl;
    BuildHopHam(1, count_dn, count_up, basis_down, index_dn, HopHam_down, Harm_Trap);
    std::cout << "No problem after down spin\n";
}

template<>
void Hamiltonian<complex<double> >::HopMatrix_Build_Peierls()
{
    //std::cout << "No problem before setting Triplet\n";
    BuildHopHam_Peierls(0, count_up, count_dn, basis_up, index_up, HopHam_up, Harm_Trap);
    
    //std::cout << "No problem after up spin\n";
    BuildHopHam_Peierls(1, count_dn, count_up, basis_down, index_dn, HopHam_down, Harm_Trap);
    
    //std::cout << "No problem after down spin\n";
}

template<typename Tnum>
void Hamiltonian<Tnum>::HopMatrix_Build_QPump()
{
    //std::cout << "No problem before setting Triplet\n";
    BuildHopHam_QPump(0, count_up, count_dn, basis_up, index_up, HopHam_up);
    //std::cout << "No problem after up spin\n";
    BuildHopHam_QPump(1, count_dn, count_up, basis_down, index_dn, HopHam_down);
    //std::cout << "No problem after down spin\n";
}

// void Hamiltonian::IntMatrix_Build()
// {
//     Ham_Interact.setFromTriplets(TL_Ubase.begin(), TL_Ubase.end());
//     cout << "No problem building interaction Ham \n";
//     //std::cout << "The Interaction Hamiltonian is: \n" << Ham_Interact  << std::endl;
// }
template<typename Tnum>
void Hamiltonian<Tnum>::Total_Ham()
{
    Ham_Tot = HopHam_up + HopHam_down + Ham_Interact;
   // cout << "Total Hamiltonian: \n" << Ham_Tot << endl;
}

template<typename Tnum>
void Hamiltonian<Tnum>::ClearInteractTriplet()
{
    // TL_Ubase.clear();
    //cout << "Is triplet set to 0? " << TL_Ubase.size() << endl;
    Ham_Interact.setZero();
    Ham_Tot.setZero();
}

template<typename Tnum>
void Hamiltonian<Tnum>::ClearHopTriplet()
{
    // TL_Ubase.clear();
    //cout << "Is triplet set to 0? " << TL_Ubase.size() << endl;
    HopHam_up.setZero();
    HopHam_down.setZero();
    Ham_Tot.setZero();
}

template<typename Tnum>
void Hamiltonian<Tnum>::OutHam()
{
   cout << "Total Hamiltonian: \n" << Ham_Tot << endl;
}


//template class Hamiltonian<int>;
template class Hamiltonian<double>;
template class Hamiltonian<complex<double> >;
