/*
 * LapackLA.cpp
 *
 *  Created on: Apr 8, 2019
 *      Author: castelletto1
 */

// Include the corresponding header file.
#include "BlasLapackLA.hpp"

// Put everything under the geosx namespace.
namespace geosx
{

//-------------------------------------------------------Mathematical methods---

real64 BlasLapackLA::vectorNorm1(array1d<real64> const & X) const
{
  return cblas_dasum( integer_conversion<int>(X.size()),
                      X.data(),
                      1);
}

real64 BlasLapackLA::vectorNorm2(array1d<real64> const & X) const
{
  return cblas_dnrm2( integer_conversion<int>(X.size()),
                      X.data(),
                      1);
}

real64 BlasLapackLA::vectorNormInf(array1d<real64> const & X) const
{
  int ind = cblas_idamax( integer_conversion<int>(X.size()),
                          X.data(),
                          1);
  return std::abs(X[ind]);
}

void BlasLapackLA::vectorVectorAdd( array1d< real64> const & X,
                                    array1d< real64> & Y,
                                    real64 const alpha)
{

  GEOS_ASSERT_MSG( X.size() == Y.size(),
                   "Vector dimensions not compatible for sum" );

  cblas_daxpy( integer_conversion<int>( X.size() ),
               alpha,
               X.data(),
               1,
               Y.data(),
               1 );
  return;
}

void BlasLapackLA::vectorScale(array1d< real64> & X,
                               real64 alpha)
{

  cblas_dscal(integer_conversion<int>(X.size()),
              alpha,
              X.data(),
              1);
  return;
}

real64 BlasLapackLA::vectorDot( array1d<real64> const & X,
                                array1d<real64> const & Y)
{
  GEOS_ASSERT_MSG( X.size() == Y.size(),
                   "Vector dimensions not compatible for dot product" );

  return cblas_ddot(integer_conversion<int>(X.size()),
                    X.data(),
                    1,
                    Y.data(),
                    1);
}

void BlasLapackLA::vectorCopy(array1d< real64> const & X,
                              array1d< real64> & Y)
{
  GEOS_ASSERT_MSG( X.size() == Y.size(),
                   "Vector dimensions not compatible for copying" );

  // Call to BLAS using CBLAS interface
  cblas_dcopy(integer_conversion<int>(X.size()),
              X.data(),
              1,
              Y.data(),
              1);
  return;
}

//----------------------------------------------------------------I/O methods---
// matrix nice output
void BlasLapackLA::print(array1d<real64> const & X)
{
  for( int i = 0 ; i < X.size() ; ++i )
  {
    printf( "%10.2e ", X[i] );
    printf( "\n" );
  }
}

} // end geosx namespace
