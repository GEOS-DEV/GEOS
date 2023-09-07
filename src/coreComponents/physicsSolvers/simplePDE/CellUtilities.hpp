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
  real64 val[3][3]{{}};
  
  void inPlaceInvert()
  {
    LvArray::tensorOps::invert< 3 >( val );
  }
};

class HexahedronCell
{
public:
  constexpr static int numVertex = 8;

  using JacobianType = Dense3x3Tensor;
  using IndexType = localIndex; // 
  // using IndexType = tripleIndex // to be added for IJK hex meshes

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

          J.val[0][0] += m_nodeCoords[vertexInd][0] * gradPhi[0]; 
          J.val[0][1] += m_nodeCoords[vertexInd][0] * gradPhi[1]; 
          J.val[0][2] += m_nodeCoords[vertexInd][0] * gradPhi[2]; 
          J.val[1][0] += m_nodeCoords[vertexInd][1] * gradPhi[0]; 
          J.val[1][1] += m_nodeCoords[vertexInd][1] * gradPhi[1]; 
          J.val[1][2] += m_nodeCoords[vertexInd][1] * gradPhi[2]; 
          J.val[2][0] += m_nodeCoords[vertexInd][2] * gradPhi[0]; 
          J.val[2][1] += m_nodeCoords[vertexInd][2] * gradPhi[1]; 
          J.val[2][2] += m_nodeCoords[vertexInd][2] * gradPhi[2]; 
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
  using IndexType = localIndex; // 
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

        J.val[0][0] += m_nodeCoords[vertexInd][0] * gradPhi[0]; 
        J.val[0][1] += m_nodeCoords[vertexInd][0] * gradPhi[1]; 
        J.val[0][2] += m_nodeCoords[vertexInd][0] * gradPhi[2]; 
        J.val[1][0] += m_nodeCoords[vertexInd][1] * gradPhi[0]; 
        J.val[1][1] += m_nodeCoords[vertexInd][1] * gradPhi[1]; 
        J.val[1][2] += m_nodeCoords[vertexInd][1] * gradPhi[2]; 
        J.val[2][0] += m_nodeCoords[vertexInd][2] * gradPhi[0]; 
        J.val[2][1] += m_nodeCoords[vertexInd][2] * gradPhi[1]; 
        J.val[2][2] += m_nodeCoords[vertexInd][2] * gradPhi[2]; 
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

} // namespace geos

#endif // GEOS_PHYSICSSOLVER_SIMPLEPDE_CELLUTILITIES_HPP_