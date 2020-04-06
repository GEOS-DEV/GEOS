/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOSX_TESTCOMPFLOWUTILS_HPP
#define GEOSX_TESTCOMPFLOWUTILS_HPP

#include "gtest/gtest.h"
#include "common/DataTypes.hpp"
#include "codingUtilities/UnitTestUtilities.hpp"

namespace geosx
{

namespace testing
{

void checkDerivative( real64 const valueEps,
                      real64 const value,
                      real64 const deriv,
                      real64 const eps,
                      real64 const relTol,
                      real64 const absTol,
                      string const & name,
                      string const & var )
{
  real64 const numDeriv = (valueEps - value) / eps;
  checkRelativeError( deriv, numDeriv, relTol, absTol, "d(" + name + ")/d(" + var + ")" );
}

void checkDerivative( real64 const valueEps,
                      real64 const value,
                      real64 const deriv,
                      real64 const eps,
                      real64 const relTol,
                      string const & name,
                      string const & var )
{ return checkDerivative( valueEps, value, deriv, eps, relTol, DEFAULT_ABS_TOL, name, var ); }

void checkDerivative( arraySlice1d< real64 const > const & valueEps,
                      arraySlice1d< real64 const > const & value,
                      arraySlice1d< real64 const > const & deriv,
                      real64 const eps,
                      real64 const relTol,
                      real64 const absTol,
                      string const & name,
                      string const & var,
                      arrayView1d< string const > const & labels )
{
  localIndex const size = labels.size( 0 );

  for( localIndex i = 0; i < size; ++i )
  {
    checkDerivative( valueEps[i], value[i], deriv[i], eps, relTol, absTol,
                     name + "[" + labels[i] + "]", var );
  }
}

template< int DIM, typename ... Args >
void checkDerivative( ArraySlice< real64 const, DIM > const & valueEps,
                      ArraySlice< real64 const, DIM > const & value,
                      ArraySlice< real64 const, DIM > const & deriv,
                      real64 const eps,
                      real64 const relTol,
                      real64 const absTol,
                      string const & name,
                      string const & var,
                      arrayView1d< string const > const & labels,
                      Args ... label_lists )
{
  localIndex const size = labels.size( 0 );

  for( localIndex i = 0; i < size; ++i )
  {
    checkDerivative( valueEps[i], value[i], deriv[i], eps, relTol, absTol,
                     name + "[" + labels[i] + "]", var, label_lists ... );
  }
}

template< int DIM, typename ... Args >
void checkDerivative( ArraySlice< real64 const, DIM > const & valueEps,
                      ArraySlice< real64 const, DIM > const & value,
                      ArraySlice< real64 const, DIM > const & deriv,
                      real64 const eps,
                      real64 const relTol,
                      string const & name,
                      string const & var,
                      arrayView1d< string const > const & labels,
                      Args ... label_lists )
{ return checkDerivative( valueEps, value, deriv, eps, relTol, DEFAULT_ABS_TOL, name, var, labels, label_lists ... ); }

// invert compositional derivative array layout to move innermost slice on the top
// (this is needed so we can use checkDerivative() to check derivative w.r.t. for each compositional var)
array1d< real64 > invertLayout( arraySlice1d< real64 const > const & input,
                                localIndex N )
{
  array1d< real64 > output( N );
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
  array2d< real64 > output( N2, N1 );

  for( localIndex i = 0; i < N1; ++i )
  {
    for( localIndex j = 0; j < N2; ++j )
    {
      output( j, i ) = input( i, j );
    }
  }

  return output;
}

array3d< real64 > invertLayout( arraySlice3d< real64 const > const & input,
                                localIndex N1,
                                localIndex N2,
                                localIndex N3 )
{
  array3d< real64 > output( N3, N1, N2 );

  for( localIndex i = 0; i < N1; ++i )
  {
    for( localIndex j = 0; j < N2; ++j )
    {
      for( localIndex k = 0; k < N3; ++k )
      {
        output( k, i, j ) = input( i, j, k );
      }
    }
  }

  return output;
}

} // namespace testing

} // namespace geosx

#endif //GEOSX_TESTCOMPFLOWUTILS_HPP
