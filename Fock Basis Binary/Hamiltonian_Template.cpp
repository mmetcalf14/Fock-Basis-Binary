//
//  Hamiltonian_Template.cpp
//  Fock Basis Binary
//
//  Created by mekena McGrew on 6/30/15.
//  Copyright (c) 2015 Mekena Metcalf. All rights reserved.
//

#include <stdio.h>
#include "Hamiltonian_Template.h"
using namespace std;
//void Hamiltonian::Set_tbar(int _tbar)
//{ tbar = _tbar;}

Basis::Basis(size_t _L, size_t _Nup, size_t _Ndn)
//:L(_L), Nup(_Nup), Ndn(_Ndn)
{
    L = _L;
    Nup = _Nup;
    Ndn = _Ndn;
    BuildBasis();
}

void Hamiltonian::Set_Const(double t_1, double t_2, double _U)
{
    J1 = t_1;
    J2 = t_2;
    //cout << J1 << " " << J2 << endl;
    U = _U;
}

void Basis::BuildBasis()
{
//    std::cout << L << " " << Nup << " " << Ndn << std::endl;
    size_t minrange_up = 0, maxrange_up = 0;
    size_t minrange_down = 0, maxrange_down = 0;
    
    //cout << "min up preloop " << minrange_up << " max up " << maxrange_up << endl;
    
    for (int i = 1; i <= Nup; i++)
    {
        
        minrange_up += pow(2,(i-1));
        maxrange_up += pow(2,(L-i));
       
    }
   // cout << "min up " << minrange_up << " max up " << maxrange_up << endl;
    
    for (int i = 1; i <= Ndn; i++)
    {  minrange_down += pow(2,(i-1));
        maxrange_down += pow(2,(L-i));
    }
    //cout << "Ranges down" << minrange_down << " " << maxrange_down <<endl;
    index_up.reserve(maxrange_up+1);
    std::vector<size_t> work (maxrange_up+1, 0);
    count_up = 0;
    for (size_t i = minrange_up; i <= maxrange_up; i++) //create spin up basis and vectors
    {
        int nbit = 0;
        for(size_t j = 0; j < L; j++)
        {
            
            if (MY_bittest(i,j))
            {
                nbit++;
            }
        }
        
        if (nbit == Nup)
        {
            count_up++;
            basis_up.push_back(i);
            work.at(i) = count_up;
            //cout << i << " " << work[i] << endl;
        }
        
    }
    index_up = work;
    
//    cout << "Basis up \n";
//        for(int i = 0; i < maxrange_down; i++)
//        { cout << basis_up[i]<<endl;}
//    cout << "Index up \n";
//    for(int i = 0; i < maxrange_up; i++)
//    { cout << index_up[i]<<endl;}
    
    
    std::cout << "End of basis up allocation \n";
    
    index_dn.reserve(maxrange_down+1);
    std::vector<size_t> nowork (maxrange_down+1, 0);
    count_dn = 0;
    //program quitting right here
   // cout << "Testing basis down vector\n";
    for (size_t i = minrange_down; i <= maxrange_down; i++) //create spin down basis and index
    {
        int nbit = 0;
        for(size_t j = 0; j < L; j++)
        {
            if (MY_bittest(i,j))
            {
                nbit++;
            }
        }
        
        if (nbit == Ndn && Ndn != 0)
        {
            count_dn++;
            //std::cout << i<< " " << count_dn << std::endl;//this is right
            basis_down.push_back(i);
            //cout << i << " " << basis_down[] << endl;//below are all wrong?? but it worked before?
            nowork.at(i) = count_dn;
            //cout << i << " " << nowork[i] << endl;
        }
        
    }
    index_dn = nowork;
//    cout << "Basis down \n";
//    for(int i = 0; i < maxrange_down; i++)
//    { cout << basis_down[i]<<endl;}
//    cout << "Index down \n";
//    for(int i = 0; i <= maxrange_down; i++)
//    { cout << index_dn[i]<<endl;}
    std::cout << "End of basis dn allocation \n";

}


void Hamiltonian::BuildHopHam_up()
{
        for(size_t bs = 0; bs < count_up; bs++)
        {
            size_t p_bas = basis_up[bs];
            
            size_t p_ind = index_up[p_bas];
            //cout << "We are acting on basis, " << p_bas << " with index, "<< p_ind << endl;
            for(size_t i = 0; i < (L-1); i++)
            {
                size_t l_bas = MY_bitset(MY_bitclr(p_bas,i),i+1);
                size_t l_ind = index_up[l_bas]; //going outside array dimension for bs
                
                if( l_bas >= index_up.size())
                {
                   l_ind = 0;
                }
                if(p_ind == (count_up) && l_bas > p_bas)// if l_ind is outside memory of index_up: CRASH!
                {
                    l_ind = 0;
                }
                //cout << p_bas << " " << p_ind <<" " << l_bas << " " << l_ind << endl;
                if(l_bas != p_bas && l_ind != 0 )
                {
                    //cout << p_bas << " " << p_ind <<" " << l_bas << " " << l_ind << endl;
                    
                double val;
                    
                if(count_dn > 0 )
                {
                    for(size_t k = 1; k <= count_dn; k++)
                    {
                        
                        
                        int r = ((k-1)*count_up) + p_ind-1;
                        int s = ((k-1)*count_up) + l_ind-1;
                        //cout << r << " " << s << endl;
                        
                        if((i%2) == 0)//even number sites have J1 hop to nn (A)
                        {
                            val = -J1 ;
                        }
                        else //odd number sites have J2 (B)
                        {
                            val = -J2;
                        }
                        
                        TL_up.push_back(Tp((r),(s),val));
                        TL_up.push_back(Tp((s),(r),val));
    
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
                    
                    TL_up.push_back(Tp(p_ind-1,l_ind-1, val ));
                    TL_up.push_back(Tp(l_ind-1,p_ind-1, val ));
                }
                }
            }
        }
 
}

void Hamiltonian::BuildHopHam_dn()
{cout << "HopHam down \n";
    for(size_t bs = 0; bs < count_dn; bs++)
    {
        size_t p_bas = basis_down[bs];
        size_t p_ind = index_dn[p_bas];
       // cout << "We are acting on basis, " << p_bas << " with index, "<< p_ind << endl;
        for(size_t i = 0; i < (L-1); i++)
        {
            size_t l_bas = MY_bitset(MY_bitclr(p_bas,i),i+1);
            size_t l_ind = index_dn[l_bas];
            
            if( l_bas >= index_dn.size())
            {
                l_ind = 0;
            }
            if(p_ind == (count_dn) && l_bas > p_bas) // if l_ind is outside memory of index_dn: CRASH!
            {
                l_ind = 0;
            }
            
            
            //cout << p_bas << " " << p_ind <<" " << l_bas << " " << l_ind << endl;
            if(l_bas != p_bas && l_ind != 0)
            {
                //cout << p_bas << " " << p_ind <<" " << l_bas << " " << l_ind << endl;
                //cout << "The hopping term gives the new state: " << l_bas << " With index: " << l_ind <<endl;
                double val;
                
                if(count_up > 0 )
                {
                    for(size_t k = 1; k <= count_up; k++)
                    {
                        int r = ((p_ind-1)*count_up) + k;
                        int s = ((l_ind-1)*count_up) + k;
                        
                        if((i%2) == 0)
                        {
                            val = -J1 ;
                        }
                        else
                        {
                            val = -J2;
                        }
                        TL_down.push_back(Tp((r-1),(s-1),val));
                        TL_down.push_back(Tp((s-1),(r-1),val));
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
                    
                    TL_down.push_back(Tp(p_ind-1,l_ind-1, val ));
                    TL_down.push_back(Tp(l_ind-1,p_ind-1, val ));
                }
            }
        }
    }
    
}

void Hamiltonian::Interaction_Index()
{
    //build interaction index for up spin
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
                    ++point_up;
                    //cout << i << " " << j_bas << " " << MY_bittest(j_bas,i) << endl;
                    IndexU_up(i,point_up-1) = j+1;
                    //std::cout << point_up << std::endl;
                    
                    if(Nup == Ndn)
                    {
                        point_dn = point_up;
                        IndexU_dn(i,point_dn-1) = j+1;
                    }
                }
            }
            else //testing if site is unoccupied, Nup + Ndn > L
            {
                if(!(MY_bittest(j_bas,i))) //testing if site is unoccupied ALG NOT WORKING
                {
                    ++point_up;
                    //cout << i << " " << j_bas << " " << ~MY_bittest(j_bas,i) << endl;
                    IndexU_up(i,point_up-1) = j+1;
                    
                    if(Nup == Ndn)
                    {
                        point_dn = point_up;
                        IndexU_dn(i,point_dn-1) = j+1;
                    }
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
                        ++point_dn;
                        
                        IndexU_dn(i,point_dn-1) = j+1;
                        
                    }
                }
                else //testing if site is unoccupied, Nup + Ndn > L
                {
                    if(!(MY_bittest(j_bas,i))) //testing if site is unoccupied
                    {
                        ++point_dn;
                        IndexU_dn(i,point_dn-1) = j+1;
                        
                    }
                }
            }
        }
   
    }
    
//    std::cout << "The spin up Index matrix is: \n" << IndexU_up << std::endl;
//    std::cout << "The spin down Index matrix is: \n" << IndexU_dn << std::endl;
}

void Hamiltonian::BaseInteraction()
{
    int g;
    if((Nup + Ndn) > L)//
    {
         g = U*(Nup + Ndn -L);
        
        for(size_t i = 0; i < Tot_base; i++)
        {
            TL_Ubase.push_back(Tp(i,i, g ));//do I need to have two different Hamiltonians?
        }
        
    }
    
    if(Nup == Ndn)
    {
        if( (Nup+Ndn) <= L)
        {
            g = U*Nup;
        }
        else
        {
            g = U*(L-Nup);//should I do this on top of other build with g?
           
        }
        for(size_t i = 1; i <= count_up; i++)
        {
            size_t k = (count_up +1)*i - count_up;//added a -1 because for loop starts at 1
            TL_Ubase.push_back(Tp(k-1,k-1, g ));//do I need to have two different Hamiltonians?
           
        }
    }
    //end function algorithm
}

void Hamiltonian::Build_Interactions()
{
    //Construct Diagonal elements from Fock space mixed indices
    for(int i = 0; i < L; i++)//The Matrix class in Eigen only allows ints for index
    {
        for( int k = 0; k < point_up; k++)
        {
            for(int l = 0; l < point_dn; l++)//Currently I'm not adding in the if statement because
            {  int r;
                if(Nup == Ndn)
                { //cout << "k: " << k << " l: " << l << endl;
                    if( l != k)
                    {
                        
                    r = ((IndexU_dn(i,l)-1)*count_up) + IndexU_up(i,k);
                        //cout << "r: " << r << endl;
                    Ham_Interact.coeffRef((r-1), (r-1)) += U;//do we need a -1? double check here
                    }
                }
                else
                {
                r = ((IndexU_dn(i,l)-1)*count_up) + IndexU_up(i,k);
                    Ham_Interact.coeffRef((r-1), (r-1)) += U;
                }
            }
        }
    }
    //cout << "The total interaction Hamiltonian is: /n"<< Ham_Interact << endl;
    //end function algorithm
}

void Hamiltonian::Set_Mat_Dim()
{
    //std::cout << "Entering dimension alg \n";
    Tot_base = count_up*count_dn; //should this be size t or int for matrix dim
    HopHam_down.resize(Tot_base, Tot_base);
    HopHam_up.resize(Tot_base, Tot_base);
    Ham_Interact.resize(Tot_base,Tot_base);
    
    IndexU_up.resize(L,count_up);
    IndexU_dn.resize(L,count_dn);
    //std::cout << "dimesnion set \n";
    
}

void Hamiltonian::HopMatrix_Build()
{
    std::cout << "No problem before setting Triplet\n";
    HopHam_up.setFromTriplets(TL_up.begin(), TL_up.end());
    std::cout << "No problem after up spin\n";
    HopHam_down.setFromTriplets(TL_down.begin(), TL_down.end());
    std::cout << "No problem after down spin\n";
   // std::cout << "The Hopping Hamiltonian up is: \n" << HopHam_up  << std::endl;
}

void Hamiltonian::IntMatrix_Build()
{
    Ham_Interact.setFromTriplets(TL_Ubase.begin(), TL_Ubase.end());
    cout << "No problem building interaction Ham \n";
//    std::cout << "The Base Interaction Hamiltonian is: \n" << Ham_Interact  << std::endl;
}

void Hamiltonian::Total_Ham()
{
    Ham_Tot = HopHam_up + HopHam_down + Ham_Interact;
    
    //cout << "Total Hamiltonian: \n" << Ham_Tot << endl;
}

void Hamiltonian::Check_Degeneracy(ofstream &EVout)
{
    Eigen::VectorXd Eval;
    Eigen::MatrixXd Evec_q;
    Eigen:: MatrixXd Temp_Ham;
    Temp_Ham = Eigen::MatrixXd(Ham_Tot);
    
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> DiagMe(Temp_Ham);
    Eval = DiagMe.eigenvalues();
    //Evec_q = DiagMe.eigenvectors();
    cout << "Full Hamiltonian Eigenvalues: \n" << Eval << endl;
    for(int i = 0; i < Eval.size(); i++)
    {
        EVout << i << " " << Eval(i) << endl;
    }
}





