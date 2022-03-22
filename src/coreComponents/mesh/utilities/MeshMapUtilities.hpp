/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file MeshMapUtilities.hpp
 */
#ifndef GEOSX_MESH_UTILITIES_MESHMAPUTILITIES_HPP
#define GEOSX_MESH_UTILITIES_MESHMAPUTILITIES_HPP

#include "common/DataTypes.hpp"
#include "mesh/ElementRegionManager.hpp"


namespace geosx
{

/**
 * @brief This namespace contains helper functions that facilitate access
 *        into the assortment of maps used by GEOSX mesh object managers
 *        (e.g. array2d/array1d(array1d)/ArrayOfArrays/ArrayOfSets, etc.)
 */
namespace meshMapUtilities
{

//////////////////////////////////////////////////////////////////////////////////

/**
 * @brief @return the size of the map along first dimension
 * @tparam T type of map element
 * @tparam USD unit-stride dimension of the map
 * @param map reference to the map
 */
template< typename T, int USD >
GEOSX_HOST_DEVICE
inline localIndex size0( arrayView2d< T, USD > const & map )
{
  return map.size( 0 );
}

/**
 * @brief @return the size of the map along first dimension
 * @tparam T type of map element
 * @param map reference to the map
 */
template< typename T >
GEOSX_HOST_DEVICE
inline localIndex size0( ArrayOfArraysView< T > const & map )
{
  return map.size();
}

/**
 * @brief @return the size of the map along first dimension
 * @tparam T type of map element
 * @param map reference to the map
 */
template< typename T >
GEOSX_HOST_DEVICE
inline localIndex size0( ArrayOfSetsView< T > const & map )
{
  return map.size();
}

//////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Transposes an input map (array2d, ArrayOfArrays or ArrayOfSets)
 * @tparam POLICY execution policy to use
 * @tparam VIEW_TYPE type of view of the source map
 * @param srcMap the source map
 * @param dstSize number of target objects ("cols" of @p srcMap)
 * @param overAlloc overallocation (extra capacity per row) for resulting map
 * @return the transpose of @p srcMap stored as ArrayOfArrays (most general form)
 */
template< typename POLICY, typename VIEW_TYPE >
ArrayOfArrays< std::remove_const_t< typename VIEW_TYPE::ValueType > >
transposeIndexMap( VIEW_TYPE const & srcMap,
                   localIndex const dstSize,
                   localIndex const overAlloc = 0 )
{
  // Count the number of elements in each set
  array1d< localIndex > counts( dstSize );
  counts.setValues< POLICY >( overAlloc );
  forAll< POLICY >( size0( srcMap ), [srcMap, counts = counts.toView()] ( localIndex const srcIndex )
  {
    for( localIndex const dstIndex : srcMap[ srcIndex ] )
    {
      RAJA::atomicInc< AtomicPolicy< POLICY > >( &counts[ dstIndex ] );
    }
  } );

  // Allocate storage for the transpose map
  ArrayOfArrays< localIndex > dstMap;
  dstMap.resizeFromCapacities< parallelHostPolicy >( dstSize, counts.data() );

  // Fill the sub-arrays with unsorted entries
  forAll< POLICY >( size0( srcMap ), [srcMap, dstMap = dstMap.toView()] ( localIndex const srcIndex )
  {
    for( localIndex const dstIndex : srcMap[ srcIndex ] )
    {
      dstMap.emplaceBackAtomic< AtomicPolicy< POLICY > >( dstIndex, srcIndex );
    }
  } );

  return dstMap;
}

/**
 * @brief Convert ToCellRelation into ToElementRelation.
 * @tparam POLICY execution policy
 * @param blockToSubRegion a map from cell blocks to region/subregion pairs
 * @param srcMap source map (object-to-cells)
 * @param dstMap target map (object-to-elements)
 */
template< typename POLICY >
void transformCellBlockToRegionMap( arrayView2d< localIndex const > const & blockToSubRegion,
                                    ToCellRelation< ArrayOfArrays< localIndex > > const & srcMap,
                                    ToElementRelation< ArrayOfArrays< localIndex > > & dstMap )
{
  GEOSX_ASSERT_EQ( blockToSubRegion.size(), 2 );
  localIndex const numObjects = srcMap.toCellIndex.size();

  localIndex const * offsets = srcMap.toCellIndex.toViewConst().getOffsets();
  dstMap.m_toElementRegion.resizeFromOffsets( numObjects, offsets );
  dstMap.m_toElementSubRegion.resizeFromOffsets( numObjects, offsets );
  dstMap.m_toElementIndex.resizeFromOffsets( numObjects, offsets );

  forAll< POLICY >( numObjects, [toCell = srcMap.toCellIndex.toViewConst(),
                                 toBlock = srcMap.toBlockIndex.toViewConst(),
                                 toRegion = dstMap.m_toElementRegion.toView(),
                                 toSubRegion = dstMap.m_toElementSubRegion.toView(),
                                 toElement = dstMap.m_toElementIndex.toView(),
                                 blockToSubRegion]( localIndex const objIndex )
  {
    arraySlice1d< localIndex const > const cells = toCell[objIndex];
    for( localIndex i = 0; i < cells.size(); ++i )
    {
      localIndex const blockIndex = toBlock( objIndex, i );
      localIndex const er = blockToSubRegion( blockIndex, 0 );
      localIndex const esr = blockToSubRegion( blockIndex, 1 );

      // Check needed because some blocks may remain unused
      if( er >= 0 && esr >= 0 )
      {
        toRegion.emplaceBack( objIndex, er );
        toSubRegion.emplaceBack( objIndex, esr );
        toElement.emplaceBack( objIndex, cells[i] );
      }
    }
  } );
}

/**
 * @brief Convert ToCellRelation into ToElementRelation.
 * @tparam POLICY execution policy
 * @param blockToSubRegion a map from cell blocks to region/subregion pairs
 * @param srcMap source map (object-to-cells)
 * @param dstMap target map (object-to-elements)
 */
template< typename POLICY, typename PERM >
void transformCellBlockToRegionMap( arrayView2d< localIndex const > const & blockToSubRegion,
                                    ToCellRelation< array2d< localIndex, PERM > > const & srcMap,
                                    ToElementRelation< array2d< localIndex, PERM > > & dstMap )
{
  GEOSX_ASSERT_EQ( blockToSubRegion.size(), 2 );
  localIndex const numObjects = srcMap.toCellIndex.size( 0 );
  localIndex const maxCellsPerObject = srcMap.toCellIndex.size( 1 );

  dstMap.m_toElementRegion.resize( numObjects, maxCellsPerObject );
  dstMap.m_toElementSubRegion.resize( numObjects, maxCellsPerObject );
  dstMap.m_toElementIndex.resize( numObjects, maxCellsPerObject );

  forAll< POLICY >( numObjects, [toCell = srcMap.toCellIndex.toViewConst(),
                                 toBlock = srcMap.toBlockIndex.toViewConst(),
                                 toRegion = dstMap.m_toElementRegion.toView(),
                                 toSubRegion = dstMap.m_toElementSubRegion.toView(),
                                 toElement = dstMap.m_toElementIndex.toView(),
                                 blockToSubRegion,
                                 maxCellsPerObject]( localIndex const objIndex )
  {
    arraySlice1d< localIndex const > const cells = toCell[objIndex];
    localIndex cellCount = 0;
    for( localIndex i = 0; i < maxCellsPerObject && cells[i] >= 0; ++i )
    {
      localIndex const blockIndex = toBlock( objIndex, i );
      localIndex const er = blockToSubRegion( blockIndex, 0 );
      localIndex const esr = blockToSubRegion( blockIndex, 1 );

      // Check needed because some blocks may remain unused
      if( er >= 0 && esr >= 0 )
      {
        toRegion( objIndex, cellCount ) = er;
        toSubRegion( objIndex, cellCount ) = esr;
        toElement( objIndex, cellCount ) = cells[i];
        ++cellCount;
      }
    }
  } );
}

} // namespace meshMapUtilities

} // namespace geosx

#endif //GEOSX_MESH_UTILITIES_MESHMAPUTILITIES_HPP
