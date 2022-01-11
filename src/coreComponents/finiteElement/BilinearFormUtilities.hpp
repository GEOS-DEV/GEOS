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

#ifndef GEOSX_FINITEELEMENT_BILINEARFORMUTILITIES_HPP_
#define GEOSX_FINITEELEMENT_BILINEARFORMUTILITIES_HPP_

#include "finiteElement/PDEUtilities.hpp"

namespace geosx
{

namespace BilinearFormUtilities
{

template< PDEUtilities::Space V,
          PDEUtilities::Space U,
          PDEUtilities::DifferentialOperator OpV,
          PDEUtilities::DifferentialOperator OpU >
struct Helper
{};

template<>
struct Helper< PDEUtilities::Space::L2,
               PDEUtilities::Space::H1vector,
               PDEUtilities::DifferentialOperator::Identity,
               PDEUtilities::DifferentialOperator::Divergence >
{
  template< int numTestDOF, int numTrialDOF >
  GEOSX_HOST_DEVICE
  void static compute( real64 (& mat)[numTestDOF][numTrialDOF],
                       real64 const (&Nv)[numTestDOF],
                       real64 const A,
                       real64 const (&dNudX)[numTrialDOF/3][3] )
  {
    for( int a = 0; a < numTestDOF; ++a )
    {
      for( int b = 0; b < numTrialDOF/3; ++b )
      {
        mat[a][b*3+0] = mat[a][b*3+0] + Nv[a] * A * dNudX[b][0];
        mat[a][b*3+1] = mat[a][b*3+1] + Nv[a] * A * dNudX[b][1];
        mat[a][b*3+2] = mat[a][b*3+2] + Nv[a] * A * dNudX[b][2];
      }
    }
  }

};

template<>
struct Helper< PDEUtilities::Space::L2,
               PDEUtilities::Space::L2,
               PDEUtilities::DifferentialOperator::Identity,
               PDEUtilities::DifferentialOperator::Identity >
{
  template< int numTestDOF, int numTrialDOF >
  GEOSX_HOST_DEVICE
  void static compute( double (& mat)[numTestDOF][numTrialDOF],
                       double const (&Nv)[numTestDOF],
                       double const A,
                       double const (&Nu)[numTrialDOF] )
  {
    for( int a = 0; a < numTestDOF; ++a )
    {
      for( int b = 0; b < numTrialDOF; ++b )
      {
        mat[a][b] = mat[a][b] + Nv[a] * A * Nu[b];
      }
    }
  }
};

template<>
struct Helper< PDEUtilities::Space::H1vector,
               PDEUtilities::Space::L2,
               PDEUtilities::DifferentialOperator::Identity,
               PDEUtilities::DifferentialOperator::Identity >
{
  template< int numTestDOF, int numTrialDOF >
  GEOSX_HOST_DEVICE
  void static compute( double (& mat)[numTestDOF][numTrialDOF],
                       double const (&N)[numTestDOF/3],
                       double const (&A)[3],
                       double const (&Np)[numTrialDOF] )
  {
    for( int a = 0; a < numTestDOF/3; ++a )
    {
      for( int b = 0; b < numTrialDOF; ++b )
      {
        mat[a*3+0][b] = mat[a*3+0][b] + N[a] * A[0] * Np[b];
        mat[a*3+1][b] = mat[a*3+1][b] + N[a] * A[1] * Np[b];
        mat[a*3+2][b] = mat[a*3+2][b] + N[a] * A[2] * Np[b];
      }
    }
  }
};

template<>
struct Helper< PDEUtilities::Space::H1vector,
               PDEUtilities::Space::L2,
               PDEUtilities::DifferentialOperator::SymmetricGradient,
               PDEUtilities::DifferentialOperator::Identity >
{
  // symmetric second-order tensor A
  template< int numTestDOF, int numTrialDOF >
  GEOSX_HOST_DEVICE
  void static compute( real64 (& mat)[numTestDOF][numTrialDOF],
                       real64 const (&dNdX)[numTestDOF/3][3],
                       real64 const (&A)[6],
                       real64 const (&Np)[numTrialDOF] )
  {
    for( int a = 0; a < numTestDOF/3; ++a )
    {
      for( int b = 0; b < numTrialDOF; ++b )
      {
        mat[a*3+0][b] = mat[a*3+0][b] + dNdX[a][0] * A[0] * Np[b] + dNdX[a][1] * A[5] * Np[b] + dNdX[a][2] * A[4] * Np[b];
        mat[a*3+1][b] = mat[a*3+1][b] + dNdX[a][0] * A[5] * Np[b] + dNdX[a][1] * A[1] * Np[b] + dNdX[a][2] * A[3] * Np[b];
        mat[a*3+2][b] = mat[a*3+2][b] + dNdX[a][0] * A[4] * Np[b] + dNdX[a][1] * A[3] * Np[b] + dNdX[a][2] * A[2] * Np[b];
      }
    }
  }

  // diagonal second-order tensor
  template< int numTestDOF, int numTrialDOF >
  GEOSX_HOST_DEVICE
  void static compute( real64 (& mat)[numTestDOF][numTrialDOF],
                       real64 const (&dNdX)[numTestDOF/3][3],
                       real64 const (&A)[3],
                       real64 const (&Np)[numTrialDOF] )
  {
    for( int a = 0; a < numTestDOF/3; ++a )
    {
      for( int b = 0; b < numTrialDOF; ++b )
      {
        mat[a*3+0][b] = mat[a*3+0][b] + dNdX[a][0] * A[0] * Np[b];
        mat[a*3+1][b] = mat[a*3+1][b] + dNdX[a][1] * A[1] * Np[b];
        mat[a*3+2][b] = mat[a*3+2][b] + dNdX[a][2] * A[2] * Np[b];
      }
    }
  }

  // scalar*identity second-order tensor
  template< int numTestDOF, int numTrialDOF >
  GEOSX_HOST_DEVICE
  void static compute( real64 (& mat)[numTestDOF][numTrialDOF],
                       real64 const (&dNdX)[numTestDOF/3][3],
                       real64 const A,
                       real64 const (&Np)[numTrialDOF] )
  {
    for( int a = 0; a < numTestDOF/3; ++a )
    {
      for( int b = 0; b < numTrialDOF; ++b )
      {
        mat[a*3+0][b] = mat[a*3+0][b] + dNdX[a][0] * A * Np[b];
        mat[a*3+1][b] = mat[a*3+1][b] + dNdX[a][1] * A * Np[b];
        mat[a*3+2][b] = mat[a*3+2][b] + dNdX[a][2] * A * Np[b];
      }
    }
  }
};

// Generic bilinear form template a(v,u)  = op1(V)^T * A * op2( U )
// V = matrix storing test space basis
// U = matrix storing trial space basis
template< PDEUtilities::Space V,
          PDEUtilities::Space U,
          PDEUtilities::DifferentialOperator OpV,
          PDEUtilities::DifferentialOperator OpU,
          typename MATRIX,
          typename V_SPACE_OPV_BASIS_VALUES,
          typename BILINEAR_FORM_MATRIX,
          typename U_SPACE_OPU_BASIS_VALUES >
GEOSX_HOST_DEVICE
static void compute( MATRIX && mat,
                     V_SPACE_OPV_BASIS_VALUES const & v,
                     BILINEAR_FORM_MATRIX const & A,
                     U_SPACE_OPU_BASIS_VALUES const & u )
{
  Helper< V, U, OpV, OpU >::compute( mat, v, A, u );
}

} // namespace BilinearFormUtilities

} // namespace geosx

#endif //GEOSX_FINITEELEMENT_BILINEARFORMUTILITIES_HPP_
