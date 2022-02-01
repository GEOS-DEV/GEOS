/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file BilinearFormUtilities.hpp
 */

#ifndef GEOSX_FINITEELEMENT_LINEARFORMUTILITIES_HPP_
#define GEOSX_FINITEELEMENT_LINEARFORMUTILITIES_HPP_

#include "finiteElement/PDEUtilities.hpp"

namespace geosx
{

namespace LinearFormUtilities
{

template< PDEUtilities::FunctionSpace V,
          PDEUtilities::DifferentialOperator OpV >
struct Helper
{};

template<>
struct Helper< PDEUtilities::FunctionSpace::P0,
               PDEUtilities::DifferentialOperator::Identity >
{
  template< int numTestDOF >
  GEOSX_HOST_DEVICE
  void static compute( real64 (& vec)[numTestDOF],
                       real64 const (&Nv)[numTestDOF],
                       real64 const A,
                       real64 const weight )
  {
    for( int a = 0; a < numTestDOF; ++a )
    {
      vec[a] = vec[a] + Nv[a] * A * weight;
    }
  }
};

template<>
struct Helper< PDEUtilities::FunctionSpace::H1vector,
               PDEUtilities::DifferentialOperator::Identity >
{
  template< int numTestDOF >
  GEOSX_HOST_DEVICE
  void static compute( double (& vec)[numTestDOF],
                       double const (&Nv)[numTestDOF/3],
                       double const (&A)[3],
                       real64 const weight )
  {
    for( int a = 0; a < numTestDOF/3; ++a )
    {
      vec[a*3+0] = vec[a*3+0] + Nv[a] * A[0] * weight;
      vec[a*3+1] = vec[a*3+1] + Nv[a] * A[1] * weight;
      vec[a*3+2] = vec[a*3+2] + Nv[a] * A[2] * weight;
    }
  }
};

template<>
struct Helper< PDEUtilities::FunctionSpace::H1vector,
               PDEUtilities::DifferentialOperator::SymmetricGradient >
{
  // symmetric second-order tensor A
  template< int numTestDOF >
  GEOSX_HOST_DEVICE
  void static compute( real64 (& vec)[numTestDOF],
                       real64 const (&dNdX)[numTestDOF/3][3],
                       real64 const (&A)[6],
                       real64 const weight )
  {
    for( int a = 0; a < numTestDOF/3; ++a )
    {
      vec[a*3+0] = vec[a*3+0] + ( dNdX[a][0] * A[0] + dNdX[a][1] * A[5] + dNdX[a][2] * A[4] ) * weight;
      vec[a*3+1] = vec[a*3+1] + ( dNdX[a][0] * A[5] + dNdX[a][1] * A[1] + dNdX[a][2] * A[3] ) * weight;
      vec[a*3+2] = vec[a*3+2] + ( dNdX[a][0] * A[4] + dNdX[a][1] * A[3] + dNdX[a][2] * A[2] ) * weight;
    }
  }

  // diagonal second-order tensor
  template< int numTestDOF >
  GEOSX_HOST_DEVICE
  void static compute( real64 (& vec)[numTestDOF],
                       real64 const (&dNdX)[numTestDOF/3][3],
                       real64 const (&A)[3],
                       real64 const weight )
  {
    for( int a = 0; a < numTestDOF/3; ++a )
    {
      vec[a*3+0] = vec[a*3+0] + dNdX[a][0] * A[0] * weight;
      vec[a*3+1] = vec[a*3+1] + dNdX[a][1] * A[1] * weight;
      vec[a*3+2] = vec[a*3+2] + dNdX[a][2] * A[2] * weight;
    }
  }

  // scalar*identity second-order tensor
  template< int numTestDOF >
  GEOSX_HOST_DEVICE
  void static compute( real64 (& vec)[numTestDOF],
                       real64 const (&dNdX)[numTestDOF/3][3],
                       real64 const A,
                       real64 const weight )
  {
    for( int a = 0; a < numTestDOF/3; ++a )
    {
      vec[a*3+0] = vec[a*3+0] + dNdX[a][0] * A * weight;
      vec[a*3+1] = vec[a*3+1] + dNdX[a][1] * A * weight;
      vec[a*3+2] = vec[a*3+2] + dNdX[a][2] * A * weight;
    }
  }
};

// Generic linear form template f(v)  = op1(V)^T * A
// V = matrix storing test space basis
template< PDEUtilities::FunctionSpace V,
          PDEUtilities::DifferentialOperator OpV,
          typename VECTOR,
          typename V_SPACE_OPV_BASIS_VALUES,
          typename TENSOR >
GEOSX_HOST_DEVICE
static void compute( VECTOR && vec,
                     V_SPACE_OPV_BASIS_VALUES const & v,
                     TENSOR const & A,
                     real64 const weight )
{
  Helper< V, OpV >::compute( vec, v, A, weight );
}

} // namespace LinearFormUtilities

} // namespace geosx

#endif //GEOSX_FINITEELEMENT_LINEARFORMUTILITIES_HPP_
