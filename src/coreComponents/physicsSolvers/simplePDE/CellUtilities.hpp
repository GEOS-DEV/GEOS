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
 * @file CellUtilities.hpp
 */

#ifndef GEOS_PHYSICSSOLVER_SIMPLEPDE_CELLUTILITIES_HPP_
#define GEOS_PHYSICSSOLVER_SIMPLEPDE_CELLUTILITIES_HPP_

#include "common/DataTypes.hpp"

namespace geos
{

//namespace cells TODO?

struct Dense3x3Tensor
{
  real64 data[3][3]{{}};

  void leftMultiplyTranspose( real64 const (&src)[3],
                              real64 (&dst)[3]  ) const
  {
      dst[0] = data[0][0] * src[0] + data[1][0] * src[1] + data[2][0] * src[2];
      dst[1] = data[0][1] * src[0] + data[1][1] * src[1] + data[2][1] * src[2];
      dst[2] = data[0][2] * src[0] + data[1][2] * src[1] + data[2][2] * src[2];
  };

  real64 inPlaceInvert()
  {
    real64 const temp[3][3] = 
    { { data[1][1]*data[2][2] - data[1][2]*data[2][1], data[0][2]*data[2][1] - data[0][1]*data[2][2], data[0][1]*data[1][2] - data[0][2]*data[1][1] },
      { data[1][2]*data[2][0] - data[1][0]*data[2][2], data[0][0]*data[2][2] - data[0][2]*data[2][0], data[0][2]*data[1][0] - data[0][0]*data[1][2] },
      { data[1][0]*data[2][1] - data[1][1]*data[2][0], data[0][1]*data[2][0] - data[0][0]*data[2][1], data[0][0]*data[1][1] - data[0][1]*data[1][0] } };   

    real64 const det = data[0][0] * temp[0][0] + data[1][0] * temp[0][1] + data[2][0] * temp[0][2];
    real64 const invDet = 1.0 / det;

    for( int i=0; i<3; ++i )
    {
      for( int j=0; j<3; ++j )
      {
        data[i][j] = temp[i][j] * invDet;
      }
    }
    return det;
  };
  
  void add_XiYj( real64 const (&X)[3],
                 real64 const (&Y)[3] )
  {
    data[0][0] += X[0] * Y[0]; 
    data[0][1] += X[0] * Y[1]; 
    data[0][2] += X[0] * Y[2]; 
    data[1][0] += X[1] * Y[0]; 
    data[1][1] += X[1] * Y[1]; 
    data[1][2] += X[1] * Y[2]; 
    data[2][0] += X[2] * Y[0]; 
    data[2][1] += X[2] * Y[1]; 
    data[2][2] += X[2] * Y[2]; 
  };


};

// namespace Polynomials
// {

// /**
//  * @brief Polynomial types
//  */
// enum class PolynomialName : integer
// {
//   Lagrange
// };

// template< PolynomialName POLYNOMIAL,
//           int NUM_NODES >
// struct Helper
// {};

// template<>
// struct Helper< PolynomialName::Lagrange,
//                2 >
// {
//   GEOS_HOST_DEVICE
//   // constexpr
//   static stackArray1d< real64, 2 > getDerivatives( real64 const Xiq[3] )
//   { 
//     stackArray1d< real64, 2 > deriv;
//     deriv[0] = -0.5;
//     deriv[1] =  0.5;
//     return deriv;
//   }

//   GEOS_HOST_DEVICE
//   static void getDerivatives( real64 const (& Xiq)[3], real64 ( & deriv )[2] )
//   { 
//     deriv[0] = -0.5;
//     deriv[1] =  0.5;
//   }

// };

// template< PolynomialName POLYNOMIAL,
//           int NUM_NODES >
// GEOS_HOST_DEVICE
// // constexpr
// static stackArray1d< real64, NUM_NODES > getDerivatives( real64 const Xiq[3] )
// {
//   return Helper< POLYNOMIAL, NUM_NODES >::getDerivatives( Xiq );
// }

// template< PolynomialName POLYNOMIAL,
//           int NUM_NODES >
// GEOS_HOST_DEVICE
// // constexpr
// static void getDerivatives( real64 const (& Xiq)[3], real64 ( & derivatives )[NUM_NODES] )
// {
//   return Helper< POLYNOMIAL, NUM_NODES >::getDerivatives( Xiq, derivatives );
// }


// /// Declare strings associated with enumeration values.
// ENUM_STRINGS( PolynomialName,
//               "Lagrange" );

// } // namespace Polynomials



namespace CellUtilities
{

template< typename CELL_TYPE >
GEOS_HOST_DEVICE
static camp::tuple< real64, typename CELL_TYPE::JacobianType >
getJacobianDeterminantAndJacobianInverse( CELL_TYPE cell,
                                          real64 const refPointCoords[3] )
{
  // Compute Jacobian
  typename CELL_TYPE::JacobianType J = cell.getJacobian( refPointCoords );

  // Compute determinant and invert Jacobian in place
  real64 const detJ = LvArray::tensorOps::invert< 3 >( J.data );
  
  return camp::make_tuple( detJ, J );
}

} // namespace CellUtilities








class HexahedronCell
{
public:
  constexpr static int numVertex = 8;

  using JacobianType = Dense3x3Tensor;

  GEOS_HOST_DEVICE
  HexahedronCell( real64 const nodeCoords[numVertex][3] )
  {
    for( int i = 0; i < numVertex; ++i )
    {
      for( int j = 0; j < 3; ++j )
      {
        m_nodeCoords[i][j] = nodeCoords[i][j];
      }
    }
  }

  GEOS_HOST_DEVICE
  JacobianType getJacobian( real64 const refPointCoords[3] ) const
  {
    JacobianType J;
    
    // Compute Jacobian
    real64 dPhiLin[2] = { -1.0, 1.0 };

    for( int k = 0; k < 2; ++k )
    {
      for( int j = 0; j < 2; ++j )
      {
        for( int i = 0; i < 2; ++i )
        {
          real64 gradPhi[3]{ 0.125 * (       dPhiLin[i]                     ) * ( 1.0 + dPhiLin[j] * refPointCoords[1] ) * ( 1.0 + dPhiLin[k] * refPointCoords[2] ),
                             0.125 * ( 1.0 + dPhiLin[i] * refPointCoords[0] ) * (       dPhiLin[j]                     ) * ( 1.0 + dPhiLin[k] * refPointCoords[2] ),
                             0.125 * ( 1.0 + dPhiLin[i] * refPointCoords[0] ) * ( 1.0 + dPhiLin[j] * refPointCoords[1] ) * (       dPhiLin[k]                     ) };

          int vertexInd = 4 * k + 2 * j + i;

          J.add_XiYj( m_nodeCoords[vertexInd], gradPhi );
        }
      }
    }

    return J;
  }

  // GEOS_HOST_DEVICE real64[3] mapping( real64 refPointCoords[3] ) const
  // {

  // }

  GEOS_HOST_DEVICE
  void getLocalCoordinates( real64 (& xLocal)[numVertex][3] ) const
  {
    for( int i = 0; i < numVertex; ++i )
    {
      xLocal[i][0] = m_nodeCoords[i][0]; 
      xLocal[i][1] = m_nodeCoords[i][1]; 
      xLocal[i][2] = m_nodeCoords[i][2]; 
    }
  }

private:
  real64 m_nodeCoords[numVertex][3];
};


// class HexadronIJKCell
// {
// public:
  // constexpr static int numVertex = 8;
// 
  // using JacobianType = real64[3][3];
// 
  // GEOS_HOST_DEVICE HexadronIJKCell( real64[3] h )
    // :
    // m_h( h )
  // {
// 
  // }
// 
  // GEOS_HOST_DEVICE JacobianType getJacobian( real64[3] refPointCoords )
  // {
    // return 0.5 * m_h;
  // }
// 
  // GEOS_HOST_DEVICE real64[3] mapping( real64[3] refPointCoords ) const
  // {
// 
  // }
// 
// private:
  // real64[3] m_h;
// };

// ***************************************************

class WedgeCell
{
public:
  constexpr static int numVertex = 6;

  using JacobianType = Dense3x3Tensor;
 
  // using IndexType = tripleIndex // to be added for IJK hex meshes

  GEOS_HOST_DEVICE
  WedgeCell( real64 const nodeCoords[numVertex][3] )
  {
    for( int i = 0; i < numVertex; ++i )
    {
      for( int j = 0; j < 3; ++j )
      {
        m_nodeCoords[i][j] = nodeCoords[i][j];
      }
    }
  }
  
  GEOS_HOST_DEVICE
  inline
  constexpr static int linearMap( int const indexT, int const indexL )
  {
    return 2 * indexT + indexL;
  }

  GEOS_HOST_DEVICE
  JacobianType getJacobian( real64 const refPointCoords[3] ) const
  {
    JacobianType J;
    
    // Compute Jacobian
    real64 const psiTRI[3] = { 1.0 - refPointCoords[0] - refPointCoords[1], refPointCoords[0], refPointCoords[1] };
    real64 const psiLIN[2] = { 0.5 - 0.5*refPointCoords[2], 0.5 + 0.5*refPointCoords[2] };
    constexpr real64 dpsiTRI[2][3] = { { -1.0, 1.0, 0.0 }, { -1.0, 0.0, 1.0 } };
    constexpr real64 dpsiLIN[2] = { -0.5, 0.5 };

    

    for( int a=0; a<3; ++a )
    {
      for( int b=0; b<2; ++b )
      {
        real64 gradPhi[3]{ dpsiTRI[0][a] * psiLIN[b],
                           dpsiTRI[1][a] * psiLIN[b],
                           psiTRI[a] * dpsiLIN[b] };

        int vertexInd = linearMap( a, b );

        J.add_XiYj( m_nodeCoords[vertexInd], gradPhi );
      }
    }

    return J;
  }

  // GEOS_HOST_DEVICE real64[3] mapping( real64 refPointCoords[3] ) const
  // {

  // }

  GEOS_HOST_DEVICE
  void getLocalCoordinates( real64 (& xLocal)[numVertex][3] ) const
  {
    for( int i = 0; i < numVertex; ++i )
    {
      xLocal[i][0] = m_nodeCoords[i][0]; 
      xLocal[i][1] = m_nodeCoords[i][1]; 
      xLocal[i][2] = m_nodeCoords[i][2]; 
    }
  }

private:
  real64 m_nodeCoords[numVertex][3];
};

// ***************************************************

class TetrahedronCell
{
public:
  constexpr static int numVertex = 4;

  using JacobianType = Dense3x3Tensor;
 
  GEOS_HOST_DEVICE
  TetrahedronCell( real64 const nodeCoords[numVertex][3] )
  {
    for( int i = 0; i < numVertex; ++i )
    {
      for( int j = 0; j < 3; ++j )
      {
        m_nodeCoords[i][j] = nodeCoords[i][j];
      }
    }
  }
  
  GEOS_HOST_DEVICE
  JacobianType getJacobian( real64 const refPointCoords[3] ) const
  {
    GEOS_UNUSED_VAR( refPointCoords );
    JacobianType J;
    J.data[0][0] = - m_nodeCoords[0][0] + m_nodeCoords[1][0];
    J.data[0][1] = - m_nodeCoords[0][0] + m_nodeCoords[2][0];
    J.data[0][2] = - m_nodeCoords[0][0] + m_nodeCoords[3][0];

    J.data[1][0] = - m_nodeCoords[0][1] + m_nodeCoords[1][1];
    J.data[1][1] = - m_nodeCoords[0][1] + m_nodeCoords[2][1];
    J.data[1][2] = - m_nodeCoords[0][1] + m_nodeCoords[3][1];

    J.data[2][0] = - m_nodeCoords[0][2] + m_nodeCoords[1][2];
    J.data[2][1] = - m_nodeCoords[0][2] + m_nodeCoords[2][2];
    J.data[2][2] = - m_nodeCoords[0][2] + m_nodeCoords[3][2];

    return J;
  }

  // GEOS_HOST_DEVICE real64[3] mapping( real64 refPointCoords[3] ) const
  // {

  // }

  GEOS_HOST_DEVICE
  void getLocalCoordinates( real64 (& xLocal)[numVertex][3] ) const
  {
    for( int i = 0; i < numVertex; ++i )
    {
      xLocal[i][0] = m_nodeCoords[i][0]; 
      xLocal[i][1] = m_nodeCoords[i][1]; 
      xLocal[i][2] = m_nodeCoords[i][2]; 
    }
  }

private:
  real64 m_nodeCoords[numVertex][3];
};

// ***************************************************

class PyramidCell
{
public:
  constexpr static int numVertex = 5;

  using JacobianType = Dense3x3Tensor;
 
  // using IndexType = tripleIndex // to be added for IJK hex meshes

  GEOS_HOST_DEVICE
  PyramidCell( real64 const nodeCoords[numVertex][3] )
  {
    for( int i = 0; i < numVertex; ++i )
    {
      for( int j = 0; j < 3; ++j )
      {
        m_nodeCoords[i][j] = nodeCoords[i][j];
      }
    }
  }
  
  GEOS_HOST_DEVICE
  inline
  constexpr static int linearMap( int const i, int const j )
  {
    return i + 2 * j;
  }

  GEOS_HOST_DEVICE
  JacobianType getJacobian( real64 const refPointCoords[3] ) const
  {
    JacobianType J;
    
    // Compute Jacobian
    real64 dPhiLin[2] = { -1.0, 1.0 };

    for( int j = 0; j < 2; ++j )
    {
      for( int i = 0; i < 2; ++i )
      {
        real64 gradPhi[3]{ 0.125 * (       dPhiLin[i]                     ) * ( 1.0 + dPhiLin[j] * refPointCoords[1] ) * ( 1.0 - refPointCoords[2] ),
                           0.125 * ( 1.0 + dPhiLin[i] * refPointCoords[0] ) * (       dPhiLin[j]                     ) * ( 1.0 - refPointCoords[2] ),
                           - 0.125 * ( 1.0 + dPhiLin[i] * refPointCoords[0] ) * ( 1.0 + dPhiLin[j] * refPointCoords[1] ) };

        int vertexInd = linearMap( i, j );

        J.add_XiYj( m_nodeCoords[vertexInd], gradPhi );
      }
    }
    
    // Contribution from the basis function paired with the apex nodes
    J.data[0][2] += m_nodeCoords[4][0] * 0.5;
    J.data[1][2] += m_nodeCoords[4][1] * 0.5;
    J.data[2][2] += m_nodeCoords[4][2] * 0.5;

    return J;
  }

  // GEOS_HOST_DEVICE real64[3] mapping( real64 refPointCoords[3] ) const
  // {

  // }

  GEOS_HOST_DEVICE
  void getLocalCoordinates( real64 (& xLocal)[numVertex][3] ) const
  {
    for( int i = 0; i < numVertex; ++i )
    {
      xLocal[i][0] = m_nodeCoords[i][0]; 
      xLocal[i][1] = m_nodeCoords[i][1]; 
      xLocal[i][2] = m_nodeCoords[i][2]; 
    }
  }

private:
  real64 m_nodeCoords[numVertex][3];
};

} // namespace geos

#endif // GEOS_PHYSICSSOLVER_SIMPLEPDE_CELLUTILITIES_HPP_