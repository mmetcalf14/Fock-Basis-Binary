#ifndef __EDTYPE_H__
#define __EDTYPE_H__
#include <iostream>
#include <stdexcept>
#include <complex>

#define INFO(MSG)             std::cout << MSG << std::endl
#define INFO_NONEWLINE(MSG)   std::cout << MSG << std::flush

#define S(x) #x
#define S_(x) S(x)
#define S__LINE__ S_(__LINE__)
/* use S__LINE__ instead of __LINE__ */
#define RUNTIME_ERROR(MSG)    throw std::runtime_error(MSG __FILE__ ": " S__LINE__ ": ")
#define LOGIC_ERROR(MSG)      throw std::logic_error(MSG __FILE__ ": " S__LINE__ ": ");
#define OVERFLOW_ERROR(MSG)   throw std::overflow_error(MSG __FILE__ ": " S__LINE__ ": ");

/** Real floating point type. */
typedef double RealType;
/** Complex type. */
typedef std::complex<double> ComplexType;

/** Matrix element type. */
// #ifdef COMPLEX_MATRIX_ELEMENS
// #warning using complex matrix elements
// typedef ComplexType MatrixElemType;
// #else
// typedef RealType MatrixElemType;
// #endif

#endif//__EDTYPE_H__
