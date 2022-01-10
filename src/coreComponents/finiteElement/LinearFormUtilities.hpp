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

template< PDEUtilities::Space V,
          PDEUtilities::DifferentialOperator OpV >
struct Helper
{};

template<>
struct Helper< PDEUtilities::Space::L2,
               PDEUtilities::DifferentialOperator::Identity >
{
  template< int numTestDOF >
  void static compute( double (& vec)[numTestDOF],
                       double const (&Nv)[numTestDOF],
                       double const A )
  {
    for( int a = 0; a < numTestDOF; ++a )
    {
      vec[a] = vec[a] + Nv[a] * A;
    }
  }
};

template<>
struct Helper< PDEUtilities::Space::H1vector,
               PDEUtilities::DifferentialOperator::Identity >
{
  template< int numTestDOF >
  void static compute( double (& vec)[numTestDOF],
                       double const (&Nv)[numTestDOF/3],
                       double const (&A)[3] )
  {
    for( int a = 0; a < numTestDOF/3; ++a )
    {
      vec[a*3+0] = vec[a*3+0] + Nv[a] * A[0];
      vec[a*3+1] = vec[a*3+1] + Nv[a] * A[1];
      vec[a*3+2] = vec[a*3+2] + Nv[a] * A[2];
    }
  }
};

template<>
struct Helper< PDEUtilities::Space::H1vector,
               PDEUtilities::DifferentialOperator::SymmetricGradient >
{
  // symmetric second-order tensor A
  template< int numTestDOF >
  GEOSX_HOST_DEVICE
  void static compute( real64 (& vec)[numTestDOF],
                       real64 const (&dNdX)[numTestDOF/3][3],
                       real64 const (&A)[6] )
  {
    for( int a = 0; a < numTestDOF/3; ++a )
    {
      vec[a*3+0] = vec[a*3+0] + dNdX[a][0] * A[0] + dNdX[a][1] * A[5] + dNdX[a][2] * A[4];
      vec[a*3+1] = vec[a*3+1] + dNdX[a][0] * A[5] + dNdX[a][1] * A[1] + dNdX[a][2] * A[3];
      vec[a*3+2] = vec[a*3+2] + dNdX[a][0] * A[4] + dNdX[a][1] * A[3] + dNdX[a][2] * A[2];
    }
  }
//
//  // diagonal second-order tensor
//  template< int numTrialDOF, int numTestDOF >
//  GEOSX_HOST_DEVICE
//  void static compute( real64 (& mat)[numTrialDOF][numTestDOF],
//                       real64 const (&dNdX)[numTrialDOF/3][3],
//                       real64 const (&A)[3],
//                       real64 const (&Np)[numTestDOF] )
//  {
//    for( int a = 0; a < numTrialDOF/3; ++a )
//    {
//      for( int b = 0; b < numTestDOF; ++b )
//      {
//        mat[a*3+0][b] = mat[a*3+0][b] + dNdX[a][0] * A[0] * Np[b];
//        mat[a*3+1][b] = mat[a*3+1][b] + dNdX[a][1] * A[1] * Np[b];
//        mat[a*3+2][b] = mat[a*3+2][b] + dNdX[a][2] * A[2] * Np[b];
//      }
//    }
//  }
//
//  // scalar*identity second-order tensor
//  template< int numTrialDOF, int numTestDOF >
//  GEOSX_HOST_DEVICE
//  void static compute( real64 (& mat)[numTrialDOF][numTestDOF],
//                       real64 const (&dNdX)[numTrialDOF/3][3],
//                       real64 const A,
//                       real64 const (&Np)[numTestDOF] )
//  {
//    for( int a = 0; a < numTrialDOF/3; ++a )
//    {
//      for( int b = 0; b < numTestDOF; ++b )
//      {
//        mat[a*3+0][b] = mat[a*3+0][b] + dNdX[a][0] * A * Np[b];
//        mat[a*3+1][b] = mat[a*3+1][b] + dNdX[a][1] * A * Np[b];
//        mat[a*3+2][b] = mat[a*3+2][b] + dNdX[a][2] * A * Np[b];
//      }
//    }
//  }
};

// Generic linear form template f(v)  = op1(V)^T * A
// V = matrix storing test space basis
template< PDEUtilities::Space V,
          PDEUtilities::DifferentialOperator OpV,
          typename VECTOR,
          typename V_SPACE_OPV_BASIS_VALUES,
          typename TENSOR >
GEOSX_HOST_DEVICE
static void compute( VECTOR && vec,
                     V_SPACE_OPV_BASIS_VALUES const & v,
                     TENSOR const & A )
{
  Helper< V, OpV >::compute( vec, v, A );
}

} // namespace LinearFormUtilities

} // namespace geosx

#endif //GEOSX_FINITEELEMENT_LINEARFORMUTILITIES_HPP_
