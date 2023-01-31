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
#include "LvArray/src/tensorOps.hpp"

namespace geosx
{

namespace BilinearFormUtilities
{

template< PDEUtilities::FunctionSpace V,
          PDEUtilities::FunctionSpace U,
          PDEUtilities::DifferentialOperator OP_V,
          PDEUtilities::DifferentialOperator OP_U >
struct Helper
{};

////////////////P0xP0 space BilinearForm functions////////////////

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

////////////////End P0xP0 space BilinearForm functions////////////////

////////////////Begin H1vector x P0 space BilinearForm functions////////////////


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

////////////////End H1vector x P0 space BilinearForm functions////////////////

////////////////Begin H1 x H1 space BilinearForm functions////////////////////

template<>
struct Helper< PDEUtilities::FunctionSpace::H1,
               PDEUtilities::FunctionSpace::H1,
               PDEUtilities::DifferentialOperator::Gradient,
               PDEUtilities::DifferentialOperator::Gradient >
{
  // Second-order tensor A
  template< int numTestDOF, int numTrialDOF >
  GEOSX_HOST_DEVICE
  void static compute( real64 (& mat)[numTestDOF][numTrialDOF],
                       real64 const (&dNvdX)[numTestDOF][3],
                       real64 const (&A)[6],
                       real64 const (&dNudX)[numTrialDOF][3],
                       real64 const weight )
  {
    for( int a = 0; a < numTestDOF; ++a )
    {
      for( int b = 0; b < numTrialDOF; ++b )
      {
        //this is working in the trivial case where A is diagonal, but I haven't tested it in general
        mat[a][b] = mat[a][b] + dNvdX[a][0] * A[0] * dNudX[b][0] * weight; //A11
        mat[a][b] = mat[a][b] + dNvdX[a][1] * A[1] * dNudX[b][1] * weight; //A22
        mat[a][b] = mat[a][b] + dNvdX[a][2] * A[2] * dNudX[b][2] * weight; //A33
        mat[a][b] = mat[a][b] + dNvdX[a][0] * A[5] * dNudX[b][1] * weight; //A12
        mat[a][b] = mat[a][b] + dNvdX[a][1] * A[5] * dNudX[b][0] * weight; //A21
        mat[a][b] = mat[a][b] + dNvdX[a][0] * A[4] * dNudX[b][2] * weight; //A13
        mat[a][b] = mat[a][b] + dNvdX[a][2] * A[4] * dNudX[b][0] * weight; //A31
        mat[a][b] = mat[a][b] + dNvdX[a][1] * A[3] * dNudX[b][2] * weight; //A23
        mat[a][b] = mat[a][b] + dNvdX[a][2] * A[3] * dNudX[b][1] * weight; //A32
      }
    }
  }

  // Scalar A
  template< int numTestDOF, int numTrialDOF >
  GEOSX_HOST_DEVICE
  void static compute( real64 (& mat)[numTestDOF][numTrialDOF],
                       real64 const (&dNvdX)[numTestDOF][3],
                       real64 const & A,
                       real64 const (&dNudX)[numTrialDOF][3],
                       real64 const weight )
  {
    //real[numTestDOF][3] dummy = LVARRAY_TENSOROPS_INIT_LOCAL_3(dNudX);
    for( int a = 0; a < numTestDOF; ++a )
    {
      for( int b = 0; b < numTrialDOF; ++b )
      {
        mat[a][b] = mat[a][b] + A * LvArray::tensorOps::AiBi< 3 >( dNvdX[a], dNudX[b] ) * weight;
      }
    }
  }

};


template<>
struct Helper< PDEUtilities::FunctionSpace::H1,
               PDEUtilities::FunctionSpace::H1,
               PDEUtilities::DifferentialOperator::Identity,
               PDEUtilities::DifferentialOperator::Identity >
{
  //scalar Trial function x scalar x Test function
  template< int numTestDOF, int numTrialDOF >
  GEOSX_HOST_DEVICE
  void static compute( real64 (& mat)[numTestDOF][numTrialDOF],
                       real64 const (&Nv)[numTestDOF],
                       real64 const (&A),
                       real64 const (&Nu)[numTrialDOF],
                       real64 const weight )
  {
    for( int a = 0; a < numTestDOF; ++a )
    {
      for( int b = 0; b < numTrialDOF; ++b )
      {
        mat[a][b] = mat[a][b] + Nv[a] * A * Nu[b] * weight;
      }
    }
  }

};

/////////////////End H1 x H1 space BilinearForm functions/////////////////////


////////////Begin H1vector x H1vector space BilinearForm functions///////////


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

////////////End H1vector x H1vector space BilinearForm functions//////////////


/**
 * @brief Generic bilinear form template to assemble elemental matrices.
 *
 * This function compute the elemental matrix for a grid cell \f$ \Omega^e \f$
 * with entries defined by the bilinear form
 *
 * \f[
 *   \int_{\Omega^e}
 *     \texttt{OP}_V ( v^h ) \; \cdot \;
 *     \mathsf{T}  \; \cdot \;
 *     \texttt{OP}_U ( u^h ) \;
 *     \mathrm{d}\Omega
 * \f]
 *
 * where \f$ v^h \f$ and \f$ u^h \f$ are finite element functions belonging to
 * the test (@p V) and  trial (@p U) space, respectively, \f$ \texttt{OP}_V \f$
 * and \f$ \texttt{OP}_U \f$ denote differential operators (zero- or
 * first-order) applied to the test and trial function, respectively, and
 * \f$ \mathsf{T} \f$ indicates a tensor.
 *
 * @tparam V Test function space.
 * @tparam U Trial function space.
 * @tparam OP_V Differential operator applied to functions in @p V.
 * @tparam OP_U Differential operator applied to functions in @p U.
 * @tparam MATRIX Derived matrix type.
 * @tparam FE_VALUES_V Derived test shape function values or derivatives type.
 * @tparam TENSOR Derived tensor type.
 * @tparam FE_VALUES_U Derived trial shape function values or derivatives type.
 * @param matrix The elemental matrix.
 * @param feValuesV Test shape function values or derivatives.
 * @param tensor The tensor defining the bilinear form.
 * @param feValuesU Trial shape function values or derivatives.
 * @param weight Quadrature weight.
 *
 */
template< PDEUtilities::FunctionSpace V,
          PDEUtilities::FunctionSpace U,
          PDEUtilities::DifferentialOperator OP_V,
          PDEUtilities::DifferentialOperator OP_U,
          typename MATRIX,
          typename FE_VALUES_V,
          typename TENSOR,
          typename FE_VALUES_U >
GEOSX_HOST_DEVICE
static void compute( MATRIX && matrix,
                     FE_VALUES_V const & feValuesV,
                     TENSOR const & tensor,
                     FE_VALUES_U const & feValuesU,
                     real64 const weight )
{
  Helper< V, U, OP_V, OP_U >::compute( matrix, feValuesV, tensor, feValuesU, weight );
}

} // namespace BilinearFormUtilities

} // namespace geosx

#endif //GEOSX_FINITEELEMENT_BILINEARFORMUTILITIES_HPP_
