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
 * @file IsoparametricMeshUtilities.hpp
 */

#ifndef GEOS_PHYSICSSOLVER_SIMPLEPDE_ISOPARAMETRICMESHUTILITIES_HPP_
#define GEOS_PHYSICSSOLVER_SIMPLEPDE_ISOPARAMETRICMESHUTILITIES_HPP_

#include "common/DataTypes.hpp"

namespace geos
{

//namespace cells TODO?
struct Dense3x3Tensor
{
  real64 val[3][3]{{}};
};

class HexadronCell
{
public:
  constexpr static int numVertex = 8;

  using JacobianType = Dense3x3Tensor;

  GEOS_HOST_DEVICE
  HexadronCell( real64 nodeCoords[numVertex][3] )
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
  JacobianType getJacobian( real64 refPointCoords[3] ) const
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
          real64 gradPhi[3]{ 0.125 * dPhiLin[i] * ( 1.0 - refPointCoords[1]) * ( 1.0 - refPointCoords[2]),
                             0.125 * ( 1.0 - refPointCoords[0]) * dPhiLin[j]  * ( 1.0 - refPointCoords[2]),
                             0.125 * ( 1.0 - refPointCoords[0]) * ( 1.0 - refPointCoords[1]) * dPhiLin[k] };

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


template< typename ARRAY_VIEW_TYPE  >
class isoparametricMesh
{
public:
  isoparametricMesh( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X,
                     ARRAY_VIEW_TYPE const elementToNodes )
   : m_X( X ), m_elementToNodes( elementToNodes )
  {

  }
  
  localIndex numCells() const
  {
    return m_elementToNodes.size( 0 );
  }

protected:
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;
  // traits::ViewTypeConst< typename SUBREGION_TYPE::NodeMapType::base_type > const m_elementToNodes;
  ARRAY_VIEW_TYPE const m_elementToNodes;
};

template< typename ARRAY_VIEW_TYPE  >
class isoparametricHexahedronMesh : public isoparametricMesh< ARRAY_VIEW_TYPE >
{
public :
  // using const_iterator = CellIterator< HexadronCell >;
  //...

 using CellType = HexadronCell;
  // ... 
  constexpr static int numCellVertex = HexadronCell::numVertex;

  // ...
  // constructor
  isoparametricHexahedronMesh( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X,
                               ARRAY_VIEW_TYPE const elementToNodes ) : isoparametricMesh< ARRAY_VIEW_TYPE >( X, elementToNodes )
  {

  }

  // 
  HexadronCell getCell( localIndex k ) const
  {
    real64 xLocal[numCellVertex][3]{};

    for( int i=0; i<numCellVertex; ++i )
    {
      localIndex const localNodeIndex = this->m_elementToNodes( k, i );
      xLocal[ i ][ 0 ] = this->m_X[ localNodeIndex ][ 0 ];
      xLocal[ i ][ 1 ] = this->m_X[ localNodeIndex ][ 1 ];
      xLocal[ i ][ 2 ] = this->m_X[ localNodeIndex ][ 2 ];
    }

    return HexadronCell( xLocal );
  }

  // localIndex numCells() const;
};

// class ijkMesh
// {
  // ...
  // constexpr static nVertex = 8;
  // ...
  // constructor
// 
  
  // HexadronIJKCell getCell( localIndex k )
  // {
    // return ...;
  // }
// private:
  // ArrayOfArrays< real64 const > const m_h;
// 
// };

// class isoparametricWedgeMesh : public isoparametricMesh
// {
//   // ... 
//   constexpr static nVertex = 6;
//   // ...
//   WedgeCell getCell( localIndex k )
//   {
//     return ...;
//   }
// };

// Factory selectign dynamically the appropriate isopametric mesh
template< typename ARRAY_VIEW_TYPE  >
isoparametricMesh< ARRAY_VIEW_TYPE > selectIsoparametricMesh( ElementType elemType,
                                           arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const nodePositions,
                                           ARRAY_VIEW_TYPE const elementToNodes )
{
  if( elemType == ElementType::Hexahedron)
  {
    return isoparametricHexahedronMesh< ARRAY_VIEW_TYPE >( nodePositions, elementToNodes );
  }
  // else if( elemType == ElementType::Wedge) )
  // {
    // return isoparametricWedgenMesh( nodePositions, elementToNodes );
  // }
  else
  {
    GEOS_ERROR( "finiteElement::dispatchlowOrder3D() is not implemented for input of " << elemType );
    return isoparametricHexahedronMesh< ARRAY_VIEW_TYPE >( nodePositions, elementToNodes );
  }
}

// ***************************************************

template< typename ARRAY_VIEW_TYPE , int NUM_VERTEX >
struct NumVertexToSubregionMesh
{
  using type = isoparametricHexahedronMesh< ARRAY_VIEW_TYPE >;
};

// template< typename ARRAY_VIEW_TYPE >
// struct NumVertexToSubregionMesh< ARRAY_VIEW_TYPE, 6 >
// {
//   using type = isoparametricWedgeMesh< ARRAY_VIEW_TYPE >;
// };

template< typename ARRAY_VIEW_TYPE >
struct NumVertexToSubregionMesh< ARRAY_VIEW_TYPE, 8 >
{
  using type = isoparametricHexahedronMesh< ARRAY_VIEW_TYPE >;
};


} // namespace geos

#endif // GEOS_PHYSICSSOLVER_SIMPLEPDE_ISOPARAMETRICMESHUTILITIES_HPP_
