/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#ifndef GEOSX_TESTCOMPFLOWUTILS_HPP
#define GEOSX_TESTCOMPFLOWUTILS_HPP

#include "gtest/gtest.h"
#include "common/DataTypes.hpp"
#include "codingUtilities/UnitTestUtilities.hpp"

namespace geosx
{

template< typename T, int NDIM >
using Array = LvArray::Array< T, NDIM, localIndex >;

namespace testing
{

template< typename T >
void checkDerivative( T valueEps,
                      T value,
                      T deriv,
                      real64 eps,
                      real64 relTol,
                      string const & name,
                      string const & var )
{
  T const numDeriv = ( valueEps - value ) / eps;
  checkRelativeError( deriv, numDeriv, relTol, "d(" + name + ")/d(" + var + ")" );
}

template< typename T, typename ... Args >
void
checkDerivative( arraySlice1d< T > const & valueEps,
                 arraySlice1d< T > const & value,
                 arraySlice1d< T > const & deriv,
                 real64 eps,
                 real64 relTol,
                 string const & name,
                 string const & var,
                 string_array const & labels,
                 Args ... label_lists )
{
  localIndex const size = labels.size( 0 );

  for( localIndex i = 0; i < size; ++i )
  {
    checkDerivative( valueEps[i], value[i], deriv[i], eps, relTol,
                     name + "[" + labels[i] + "]", var, label_lists... );
  }
}

template< typename T, int DIM, typename ... Args >
typename std::enable_if< ( DIM > 1 ), void >::type
checkDerivative( array_slice< T, DIM > const & valueEps,
                 array_slice< T, DIM > const & value,
                 array_slice< T, DIM > const & deriv,
                 real64 eps,
                 real64 relTol,
                 string const & name,
                 string const & var,
                 string_array const & labels,
                 Args ... label_lists )
{
  localIndex const size = labels.size( 0 );

  for( localIndex i = 0; i < size; ++i )
  {
    checkDerivative( valueEps[i], value[i], deriv[i], eps, relTol,
                     name + "[" + labels[i] + "]", var, label_lists... );
  }
}

// invert compositional derivative array layout to move innermost slice on the top
// (this is needed so we can use checkDerivative() to check derivative w.r.t. for each compositional var)
array1d< real64 > invertLayout( arraySlice1d< real64 const > const & input,
                                localIndex N )
{
  Array< real64, 1 > output( N );
  for( int i = 0; i < N; ++i )
  {
    output[i] = input[i];
  }

  return output;
}

array2d< real64 > invertLayout( arraySlice2d< real64 const > const & input,
                                localIndex N1,
                                localIndex N2 )
{
  Array< real64, 2 > output( N2, N1 );

  for( localIndex i = 0; i < N1; ++i )
  {
    for( localIndex j = 0; j < N2; ++j )
    {
      output[j][i] = input[i][j];
    }
  }

  return output;
}

array3d< real64 > invertLayout( arraySlice3d< real64 const > const & input,
                                localIndex N1,
                                localIndex N2,
                                localIndex N3 )
{
  Array< real64, 3 > output( N3, N1, N2 );

  for( localIndex i = 0; i < N1; ++i )
  {
    for( localIndex j = 0; j < N2; ++j )
    {
      for( localIndex k = 0; k < N3; ++k )
      {
        output[k][i][j] = input[i][j][k];
      }
    }
  }

  return output;
}

} // namespace testing

} // namespace geosx

#endif //GEOSX_TESTCOMPFLOWUTILS_HPP
