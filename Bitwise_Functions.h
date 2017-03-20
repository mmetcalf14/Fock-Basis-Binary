//
//  Bitwise_Functions.h
//  ED Hubbard Hamiltonian
//
//  Created by mekena McGrew on 11/4/15.
//  Copyright © 2015 Mekena Metcalf. All rights reserved.
//

#ifndef Bitwise_Functions_h
#define Bitwise_Functions_h

//Bitewise functions definitions used on Fock states
//Binary Function algorithm
inline size_t MY_bittest(size_t m, size_t n)// m -> basis integer, n -> site
{
    size_t Eval;//if Eval is size_t I get a totally wrong number compared to int
    //seg fault occurring regardless of whether return value is correct or incorrect
    //std::cout << m << " " << n << std::endl;
    Eval = (m & (1 << n));//I haven't changed anything why not working all of a sudden?
    return Eval;
}

inline size_t MY_bitclr(size_t m,  size_t n) // set nth bit to zero
{
    size_t Clr_bit;
    Clr_bit = m & ~(1 << n);
    return Clr_bit;
}

inline size_t MY_bitset(size_t m,  size_t n)
{
    size_t New_State;
    New_State = m | (1 << n);
    return New_State;
}

#endif /* Bitwise_Functions_h */
