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

template< PDEUtilities::FunctionSpace V,
          PDEUtilities::FunctionSpace U,
          PDEUtilities::DifferentialOperator OpV,
          PDEUtilities::DifferentialOperator OpU >
struct Helper
{};

template<>
struct Helper< PDEUtilities::FunctionSpace::P0,
               PDEUtilities::FunctionSpace::H1vector,
               PDEUtilities::DifferentialOperator::Identity,
               PDEUtilities::DifferentialOperator::Divergence >
{
  template< int numTestDOF, int numTrialDOF >
  GEOSX_HOST_DEVICE
  void static compute( real64 (& mat)[numTestDOF][numTrialDOF],
                       real64 const & Nv,
                       real64 const A,
                       real64 const (&dNudX)[numTrialDOF/3][3],
                       real64 const weight )
  {
    GEOSX_UNUSED_VAR( Nv );
    for( int a = 0; a < numTestDOF; ++a )
    {
      for( int b = 0; b < numTrialDOF/3; ++b )
      {
        mat[a][b*3+0] = mat[a][b*3+0] + A * dNudX[b][0] * weight;
        mat[a][b*3+1] = mat[a][b*3+1] + A * dNudX[b][1] * weight;
        mat[a][b*3+2] = mat[a][b*3+2] + A * dNudX[b][2] * weight;
      }
    }
  }

  template< int numTestDOF, int numTrialDOF >
  GEOSX_HOST_DEVICE
  void static compute( real64 (& mat)[numTestDOF][numTrialDOF],
                       real64 const & Nv,
                       real64 const (&A)[numTestDOF],
                       real64 const (&dNudX)[numTrialDOF/3][3],
                       real64 const weight )
  {
    GEOSX_UNUSED_VAR( Nv );
    for( int a = 0; a < numTestDOF; ++a )
    {
      for( int b = 0; b < numTrialDOF/3; ++b )
      {
        mat[a][b*3+0] = mat[a][b*3+0] + A[a] * dNudX[b][0] * weight;
        mat[a][b*3+1] = mat[a][b*3+1] + A[a] * dNudX[b][1] * weight;
        mat[a][b*3+2] = mat[a][b*3+2] + A[a] * dNudX[b][2] * weight;
      }
    }
  }

};

template<>
struct Helper< PDEUtilities::FunctionSpace::P0,
               PDEUtilities::FunctionSpace::P0,
               PDEUtilities::DifferentialOperator::Identity,
               PDEUtilities::DifferentialOperator::Identity >
{
  template< int numTestDOF, int numTrialDOF >
  GEOSX_HOST_DEVICE
  void static compute( double (& mat)[numTestDOF][numTrialDOF],
                       double const & Nv,
                       double const A,
                       double const & Nu,
                       real64 const weight )
  {
    GEOSX_UNUSED_VAR( Nv );
    GEOSX_UNUSED_VAR( Nu );
    for( int a = 0; a < numTestDOF; ++a )
    {
      for( int b = 0; b < numTrialDOF; ++b )
      {
        mat[a][b] = mat[a][b] + A * weight;
      }
    }
  }

  template< int numTestDOF, int numTrialDOF >
  GEOSX_HOST_DEVICE
  void static compute( double (& mat)[numTestDOF][numTrialDOF],
                       double const & Nv,
                       double const A[numTestDOF],
                       double const & Nu,
                       real64 const weight )
  {
    GEOSX_UNUSED_VAR( Nv );
    GEOSX_UNUSED_VAR( Nu );
    for( int a = 0; a < numTestDOF; ++a )
    {
      for( int b = 0; b < numTrialDOF; ++b )
      {
        mat[a][b] = mat[a][b] + A[a] * weight;
      }
    }
  }

  template< int numTestDOF, int numTrialDOF >
  GEOSX_HOST_DEVICE
  void static compute( double (& mat)[numTestDOF][numTrialDOF],
                       double const & Nv,
                       double const A[numTestDOF][numTrialDOF],
                       double const & Nu,
                       real64 const weight )
  {
    GEOSX_UNUSED_VAR( Nv );
    GEOSX_UNUSED_VAR( Nu );
    for( int a = 0; a < numTestDOF; ++a )
    {
      for( int b = 0; b < numTrialDOF; ++b )
      {
        mat[a][b] = mat[a][b] + A[a][b] * weight;
      }
    }
  }
};

template<>
struct Helper< PDEUtilities::FunctionSpace::H1vector,
               PDEUtilities::FunctionSpace::P0,
               PDEUtilities::DifferentialOperator::Identity,
               PDEUtilities::DifferentialOperator::Identity >
{
  template< int numTestDOF, int numTrialDOF >
  GEOSX_HOST_DEVICE
  void static compute( real64 (& mat)[numTestDOF][numTrialDOF],
                       real64 const (&N)[numTestDOF/3],
                       real64 const (&A)[3],
                       real64 const & Np,
                       real64 const weight )
  {
    GEOSX_UNUSED_VAR( Np );
    for( int a = 0; a < numTestDOF/3; ++a )
    {
      for( int b = 0; b < numTrialDOF; ++b )
      {
        mat[a*3+0][b] = mat[a*3+0][b] + N[a] * A[0] * weight;
        mat[a*3+1][b] = mat[a*3+1][b] + N[a] * A[1] * weight;
        mat[a*3+2][b] = mat[a*3+2][b] + N[a] * A[2] * weight;
      }
    }
  }

  template< int numTestDOF, int numTrialDOF >
  GEOSX_HOST_DEVICE
  void static compute( real64 (& mat)[numTestDOF][numTrialDOF],
                       real64 const (&N)[numTestDOF/3],
                       real64 const (&A)[3][numTrialDOF],
                       real64 const & Np,
                       real64 const weight )
  {
    GEOSX_UNUSED_VAR( Np );
    for( int a = 0; a < numTestDOF/3; ++a )
    {
      for( int b = 0; b < numTrialDOF; ++b )
      {
        mat[a*3+0][b] = mat[a*3+0][b] + N[a] * A[0][b] * weight;
        mat[a*3+1][b] = mat[a*3+1][b] + N[a] * A[1][b] * weight;
        mat[a*3+2][b] = mat[a*3+2][b] + N[a] * A[2][b] * weight;
      }
    }
  }
};

template<>
struct Helper< PDEUtilities::FunctionSpace::H1vector,
               PDEUtilities::FunctionSpace::P0,
               PDEUtilities::DifferentialOperator::SymmetricGradient,
               PDEUtilities::DifferentialOperator::Identity >
{
  // symmetric second-order tensor A
  template< int numTestDOF, int numTrialDOF >
  GEOSX_HOST_DEVICE
  void static compute( real64 (& mat)[numTestDOF][numTrialDOF],
                       real64 const (&dNdX)[numTestDOF/3][3],
                       real64 const (&A)[6],
                       real64 const & Np,
                       real64 const weight )
  {
    GEOSX_UNUSED_VAR( Np );
    for( int a = 0; a < numTestDOF/3; ++a )
    {
      for( int b = 0; b < numTrialDOF; ++b )
      {
        mat[a*3+0][b] = mat[a*3+0][b] + ( dNdX[a][0] * A[0] + dNdX[a][1] * A[5] + dNdX[a][2] * A[4] ) * weight;
        mat[a*3+1][b] = mat[a*3+1][b] + ( dNdX[a][0] * A[5] + dNdX[a][1] * A[1] + dNdX[a][2] * A[3] ) * weight;
        mat[a*3+2][b] = mat[a*3+2][b] + ( dNdX[a][0] * A[4] + dNdX[a][1] * A[3] + dNdX[a][2] * A[2] ) * weight;
      }
    }
  }

  // diagonal second-order tensor
  template< int numTestDOF, int numTrialDOF >
  GEOSX_HOST_DEVICE
  void static compute( real64 (& mat)[numTestDOF][numTrialDOF],
                       real64 const (&dNdX)[numTestDOF/3][3],
                       real64 const (&A)[3],
                       real64 const & Np,
                       real64 const weight )
  {
    GEOSX_UNUSED_VAR( Np );
    for( int a = 0; a < numTestDOF/3; ++a )
    {
      for( int b = 0; b < numTrialDOF; ++b )
      {
        mat[a*3+0][b] = mat[a*3+0][b] + dNdX[a][0] * A[0] * weight;
        mat[a*3+1][b] = mat[a*3+1][b] + dNdX[a][1] * A[1] * weight;
        mat[a*3+2][b] = mat[a*3+2][b] + dNdX[a][2] * A[2] * weight;
      }
    }
  }

  // scalar*identity second-order tensor
  template< int numTestDOF, int numTrialDOF >
  GEOSX_HOST_DEVICE
  void static compute( real64 (& mat)[numTestDOF][numTrialDOF],
                       real64 const (&dNdX)[numTestDOF/3][3],
                       real64 const A,
                       real64 const & Np,
                       real64 const weight )
  {
    GEOSX_UNUSED_VAR( Np );
    for( int a = 0; a < numTestDOF/3; ++a )
    {
      for( int b = 0; b < numTrialDOF; ++b )
      {
        mat[a*3+0][b] = mat[a*3+0][b] + dNdX[a][0] * A * weight;
        mat[a*3+1][b] = mat[a*3+1][b] + dNdX[a][1] * A * weight;
        mat[a*3+2][b] = mat[a*3+2][b] + dNdX[a][2] * A * weight;
      }
    }
  }
};

template<>
struct Helper< PDEUtilities::FunctionSpace::H1vector,
               PDEUtilities::FunctionSpace::H1vector,
               PDEUtilities::DifferentialOperator::SymmetricGradient,
               PDEUtilities::DifferentialOperator::SymmetricGradient >
{
  // Fourth-order tensor A
  template< int numTestDOF, int numTrialDOF, typename DISCRETIZATION_OPS >
  GEOSX_HOST_DEVICE
  void static compute( real64 (& mat)[numTestDOF][numTrialDOF],
                       real64 const (&dNvdX)[numTestDOF/3][3],
                       DISCRETIZATION_OPS const & A,
                       real64 const (&dNudX)[numTrialDOF/3][3],
                       real64 const weight )
  {
    GEOSX_UNUSED_VAR( dNvdX );
    const_cast< DISCRETIZATION_OPS & >(A).template BTDB< numTestDOF/3 >( dNudX, weight, mat );
  }
};

template<>
struct Helper< PDEUtilities::FunctionSpace::H1vector,
               PDEUtilities::FunctionSpace::H1vector,
               PDEUtilities::DifferentialOperator::Identity,
               PDEUtilities::DifferentialOperator::Divergence >
{
  template< int numTestDOF, int numTrialDOF >
  GEOSX_HOST_DEVICE
  void static compute( real64 (& mat)[numTestDOF][numTrialDOF],
                       real64 const (&Nv)[numTestDOF/3],
                       real64 const (&A)[3],
                       real64 const (&dNudX)[numTestDOF/3][3],
                       real64 const weight )
  {
    for( int a = 0; a < numTestDOF/3; ++a )
    {
      for( int b = 0; b < numTrialDOF/3; ++b )
      {
        mat[a*3+0][b*3+0] = mat[a*3+0][b*3+0] + Nv[a] * A[0] * dNudX[b][0] * weight;
        mat[a*3+0][b*3+1] = mat[a*3+0][b*3+1] + Nv[a] * A[0] * dNudX[b][1] * weight;
        mat[a*3+0][b*3+2] = mat[a*3+0][b*3+2] + Nv[a] * A[0] * dNudX[b][2] * weight;
        mat[a*3+1][b*3+0] = mat[a*3+1][b*3+0] + Nv[a] * A[1] * dNudX[b][0] * weight;
        mat[a*3+1][b*3+1] = mat[a*3+1][b*3+1] + Nv[a] * A[1] * dNudX[b][1] * weight;
        mat[a*3+1][b*3+2] = mat[a*3+1][b*3+2] + Nv[a] * A[1] * dNudX[b][2] * weight;
        mat[a*3+2][b*3+0] = mat[a*3+2][b*3+0] + Nv[a] * A[2] * dNudX[b][0] * weight;
        mat[a*3+2][b*3+1] = mat[a*3+2][b*3+1] + Nv[a] * A[2] * dNudX[b][1] * weight;
        mat[a*3+2][b*3+2] = mat[a*3+2][b*3+2] + Nv[a] * A[2] * dNudX[b][2] * weight;
      }
    }
  }

};

// Generic bilinear form template a(v,u)  = op1(V)^T * A * op2( U )
// V = matrix storing test space basis
// U = matrix storing trial space basis
template< PDEUtilities::FunctionSpace V,
          PDEUtilities::FunctionSpace U,
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
                     U_SPACE_OPU_BASIS_VALUES const & u,
                     real64 const weight )
{
  Helper< V, U, OpV, OpU >::compute( mat, v, A, u, weight );
}

} // namespace BilinearFormUtilities

} // namespace geosx

#endif //GEOSX_FINITEELEMENT_BILINEARFORMUTILITIES_HPP_
