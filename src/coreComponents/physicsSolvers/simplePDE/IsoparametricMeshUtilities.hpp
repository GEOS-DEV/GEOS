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

#include "CellUtilities.hpp"

namespace geos
{

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

 using CellType = HexahedronCell;
  // ... 
  constexpr static int numCellVertex = HexahedronCell::numVertex;

  // ...
  // constructor
  isoparametricHexahedronMesh( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X,
                               ARRAY_VIEW_TYPE const elementToNodes ) : isoparametricMesh< ARRAY_VIEW_TYPE >( X, elementToNodes )
  {

  }

  // 
  HexahedronCell getCell( localIndex k ) const
  {
    real64 xLocal[numCellVertex][3]{};

    for( int i=0; i<numCellVertex; ++i )
    {
      localIndex const localNodeIndex = this->m_elementToNodes( k, i );
      xLocal[ i ][ 0 ] = this->m_X[ localNodeIndex ][ 0 ];
      xLocal[ i ][ 1 ] = this->m_X[ localNodeIndex ][ 1 ];
      xLocal[ i ][ 2 ] = this->m_X[ localNodeIndex ][ 2 ];
    }

    return HexahedronCell( xLocal );
  }

  // localInde{x numCells() const;
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
