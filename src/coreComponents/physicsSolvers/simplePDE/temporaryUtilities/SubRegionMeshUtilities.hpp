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
 * @file SubRegionMeshUtilities.hpp
 */

#ifndef GEOS_PHYSICSSOLVER_SIMPLEPDE_SUBREGIONMESHUTILITIES_HPP_
#define GEOS_PHYSICSSOLVER_SIMPLEPDE_SUBREGIONMESHUTILITIES_HPP_

#include "common/DataTypes.hpp"

#include "CellUtilities.hpp"

namespace geos
{

template< typename ARRAY_VIEW_TYPE >
class isoparametricMesh
{
public:
  using CellIndexType = localIndex;

  isoparametricMesh( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X,
                     ARRAY_VIEW_TYPE const elementToNodes )
    : m_X( X ), m_elementToNodes( elementToNodes )
  {}

  localIndex numCells() const
  {
    return m_elementToNodes.size( 0 );
  }

protected:
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;
  // traits::ViewTypeConst< typename SUBREGION_TYPE::NodeMapType::base_type > const m_elementToNodes;
  ARRAY_VIEW_TYPE const m_elementToNodes;
};

template< typename ARRAY_VIEW_TYPE >
class IJKMesh
{
public:
  using CellIndexType = localIndex; // Should be a triple index
  using CellType = CubeCell;

  IJKMesh( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X,
           ARRAY_VIEW_TYPE const elementToNodes )
    : m_h( 0.1 )
  {
    GEOS_UNUSED_VAR( X );
    GEOS_UNUSED_VAR( elementToNodes );

  }

  //
  CubeCell getCell( CellIndexType k ) const
  {
    GEOS_UNUSED_VAR( k );
    return CubeCell( this->m_h );
  }

  // localIndex numCells() const
  // {
  //   return m_elementToNodes.size( 0 );
  // }

protected:
  real64 m_h;
  // traits::ViewTypeConst< typename SUBREGION_TYPE::NodeMapType::base_type > const m_elementToNodes;
  ARRAY_VIEW_TYPE const m_elementToNodes;
};


template< typename ARRAY_VIEW_TYPE >
class isoparametricCuboidMesh : public isoparametricMesh< ARRAY_VIEW_TYPE >
{
public:
  // using const_iterator = CellIterator< HexadronCell >;
  //...
  using CellIndexType = typename isoparametricMesh< ARRAY_VIEW_TYPE >::CellIndexType;

  using CellType = CuboidCell;
  constexpr static int numCellVertex = CuboidCell::numVertex;

  // ...

  // ...
  // constructor
  isoparametricCuboidMesh( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X,
                           ARRAY_VIEW_TYPE const elementToNodes ): isoparametricMesh< ARRAY_VIEW_TYPE >( X, elementToNodes )
  {}

  //
  CuboidCell getCell( CellIndexType k ) const
  {
    real64 xLocal[numCellVertex][3]{};

    for( int i=0; i<numCellVertex; ++i )
    {
      localIndex const localNodeIndex = this->m_elementToNodes( k, i );
      xLocal[ i ][ 0 ] = this->m_X[ localNodeIndex ][ 0 ];
      xLocal[ i ][ 1 ] = this->m_X[ localNodeIndex ][ 1 ];
      xLocal[ i ][ 2 ] = this->m_X[ localNodeIndex ][ 2 ];
    }

    return CuboidCell( xLocal );
  }

  // localInde{x numCells() const;
};


template< typename ARRAY_VIEW_TYPE >
class isoparametricWedgeMesh : public isoparametricMesh< ARRAY_VIEW_TYPE >
{
public:

  using CellIndexType = typename isoparametricMesh< ARRAY_VIEW_TYPE >::CellIndexType;
  using CellType = WedgeCell;
  // ...
  constexpr static int numCellVertex = WedgeCell::numVertex;

  // ...
  // constructor
  isoparametricWedgeMesh( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X,
                          ARRAY_VIEW_TYPE const elementToNodes ): isoparametricMesh< ARRAY_VIEW_TYPE >( X, elementToNodes )
  {}

  //
  WedgeCell getCell( CellIndexType k ) const
  {
    real64 xLocal[numCellVertex][3]{};

    for( int i=0; i<numCellVertex; ++i )
    {
      localIndex const localNodeIndex = this->m_elementToNodes( k, i );
      xLocal[ i ][ 0 ] = this->m_X[ localNodeIndex ][ 0 ];
      xLocal[ i ][ 1 ] = this->m_X[ localNodeIndex ][ 1 ];
      xLocal[ i ][ 2 ] = this->m_X[ localNodeIndex ][ 2 ];
    }

    return WedgeCell( xLocal );
  }

  // localInde{x numCells() const;
};


template< typename ARRAY_VIEW_TYPE >
class isoparametricTetrahedronMesh : public isoparametricMesh< ARRAY_VIEW_TYPE >
{
public:

  using CellIndexType = typename isoparametricMesh< ARRAY_VIEW_TYPE >::CellIndexType;
  using CellType = TetrahedronCell;
  // ...
  constexpr static int numCellVertex = TetrahedronCell::numVertex;

  // ...
  // constructor
  isoparametricTetrahedronMesh( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X,
                                ARRAY_VIEW_TYPE const elementToNodes ): isoparametricMesh< ARRAY_VIEW_TYPE >( X, elementToNodes )
  {}

  //
  TetrahedronCell getCell( CellIndexType k ) const
  {
    real64 xLocal[numCellVertex][3]{};

    for( int i=0; i<numCellVertex; ++i )
    {
      localIndex const localNodeIndex = this->m_elementToNodes( k, i );
      xLocal[ i ][ 0 ] = this->m_X[ localNodeIndex ][ 0 ];
      xLocal[ i ][ 1 ] = this->m_X[ localNodeIndex ][ 1 ];
      xLocal[ i ][ 2 ] = this->m_X[ localNodeIndex ][ 2 ];
    }

    return TetrahedronCell( xLocal );
  }

  // localInde{x numCells() const;
};


template< typename ARRAY_VIEW_TYPE >
class isoparametricPyramidMesh : public isoparametricMesh< ARRAY_VIEW_TYPE >
{
public:

  using CellIndexType = typename isoparametricMesh< ARRAY_VIEW_TYPE >::CellIndexType;
  using CellType = PyramidCell;
  // ...
  constexpr static int numCellVertex = PyramidCell::numVertex;

  // ...
  // constructor
  isoparametricPyramidMesh( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X,
                            ARRAY_VIEW_TYPE const elementToNodes ): isoparametricMesh< ARRAY_VIEW_TYPE >( X, elementToNodes )
  {}

  //
  PyramidCell getCell( CellIndexType k ) const
  {
    real64 xLocal[numCellVertex][3]{};

    for( int i=0; i<numCellVertex; ++i )
    {
      localIndex const localNodeIndex = this->m_elementToNodes( k, i );
      xLocal[ i ][ 0 ] = this->m_X[ localNodeIndex ][ 0 ];
      xLocal[ i ][ 1 ] = this->m_X[ localNodeIndex ][ 1 ];
      xLocal[ i ][ 2 ] = this->m_X[ localNodeIndex ][ 2 ];
    }

    return PyramidCell( xLocal );
  }

  // localInde{x numCells() const;
};


// // Factory selectign dynamically the appropriate isopametric mesh
// template< typename ARRAY_VIEW_TYPE >
// isoparametricMesh< ARRAY_VIEW_TYPE > selectIsoparametricMesh( ElementType elemType,
//                                                               arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const
// nodePositions,
//                                                               ARRAY_VIEW_TYPE const elementToNodes )
// {
//   if( elemType == ElementType::Hexahedron )
//   {
//     return isoparametricHexahedronMesh< ARRAY_VIEW_TYPE >( nodePositions, elementToNodes );
//   }
//   // else if( elemType == ElementType::Wedge) )
//   // {
//   // return isoparametricWedgenMesh( nodePositions, elementToNodes );
//   // }
//   else
//   {
//     GEOS_ERROR( "finiteElement::dispatchlowOrder3D() is not implemented for input of " << elemType );
//     return isoparametricHexahedronMesh< ARRAY_VIEW_TYPE >( nodePositions, elementToNodes );
//   }
// }

// ***************************************************

template< typename ARRAY_VIEW_TYPE, int NUM_VERTEX >
struct NumVertexToSubregionMesh
{
  using type = isoparametricCuboidMesh< ARRAY_VIEW_TYPE >;
};

// template< typename ARRAY_VIEW_TYPE >
// struct NumVertexToSubregionMesh< ARRAY_VIEW_TYPE, 6 >
// {
//   using type = isoparametricWedgeMesh< ARRAY_VIEW_TYPE >;

template< typename ARRAY_VIEW_TYPE >
struct NumVertexToSubregionMesh< ARRAY_VIEW_TYPE, 4 >
{
  using type = isoparametricTetrahedronMesh< ARRAY_VIEW_TYPE >;
};

template< typename ARRAY_VIEW_TYPE >
struct NumVertexToSubregionMesh< ARRAY_VIEW_TYPE, 5 >
{
  using type = isoparametricPyramidMesh< ARRAY_VIEW_TYPE >;
};

template< typename ARRAY_VIEW_TYPE >
struct NumVertexToSubregionMesh< ARRAY_VIEW_TYPE, 6 >
{
  using type = isoparametricWedgeMesh< ARRAY_VIEW_TYPE >;
};
// };

// #define USE_IJK_MESH

#ifdef USE_IJK_MESH
template< typename ARRAY_VIEW_TYPE >
struct NumVertexToSubregionMesh< ARRAY_VIEW_TYPE, 8 >
{
  using type = IJKMesh< ARRAY_VIEW_TYPE >;
};
#else
template< typename ARRAY_VIEW_TYPE >
struct NumVertexToSubregionMesh< ARRAY_VIEW_TYPE, 8 >
{
  using type = isoparametricCuboidMesh< ARRAY_VIEW_TYPE >;
};
#endif

} // namespace geos

#endif // GEOS_PHYSICSSOLVER_SIMPLEPDE_SUBREGIONMESHUTILITIES_HPP_
