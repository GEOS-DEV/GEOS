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
 * @param map reference to the map
 */
template< typename T >
GEOSX_HOST_DEVICE
inline localIndex size0( arrayView1d< arrayView1d< T const > const > const & map )
{
  return map.size();
}

/**
 * @copydoc size0(arrayView1d< arrayView1d< T const > const > const &)
 * @tparam USD unit-stride dimension of the map
 */
template< typename T, int USD >
GEOSX_HOST_DEVICE
inline localIndex size0( arrayView2d< T, USD > const & map )
{
  return map.size( 0 );
}

/**
 * @copydoc size0(arrayView1d< arrayView1d< T const > const > const &)
 */
template< typename T >
GEOSX_HOST_DEVICE
inline localIndex size0( ArrayOfArraysView< T const > const & map )
{
  return map.size();
}

/**
 * @copydoc size0(arrayView1d< arrayView1d< T const > const > const &)
 */
template< typename T >
GEOSX_HOST_DEVICE
inline localIndex size0( ArrayOfSetsView< T const > const & map )
{
  return map.size();
}

//////////////////////////////////////////////////////////////////////////////////

/**
 * @brief @return the size of the map along second dimension
 * @tparam T type of map element
 * @param map reference to the map
 * @param i0 first-dimension index into the map
 */
template< typename T >
GEOSX_HOST_DEVICE
inline localIndex size1( arrayView1d< arrayView1d< T const > const > const & map, localIndex const i0 )
{
  return map[i0].size();
}

/**
 * @copydoc size1(arrayView1d< arrayView1d< T const > const > const &, localIndex const)
 * @tparam USD unit-stride dimension of the map
 */
template< typename T, int USD >
GEOSX_HOST_DEVICE
inline localIndex size1( arrayView2d< T, USD > const & map, localIndex const i0 )
{
  GEOSX_UNUSED_VAR( i0 );
  return map.size( 1 );
}

/**
 * @copydoc size1(arrayView1d< arrayView1d< T const > const > const &, localIndex const)
 */
template< typename T >
GEOSX_HOST_DEVICE
inline localIndex size1( ArrayOfArraysView< T const > const & map, localIndex const i0 )
{
  return map.sizeOfArray( i0 );
}

/**
 * @copydoc size1(arrayView1d< arrayView1d< T const > const > const &, localIndex const)
 */
template< typename T >
GEOSX_HOST_DEVICE
inline localIndex size1( ArrayOfSetsView< T const > const & map, localIndex const i0 )
{
  return map.sizeOfSet( i0 );
}

//////////////////////////////////////////////////////////////////////////////////

/**
 * @brief @return the value of the map
 * @tparam T type of map element
 * @param map reference to the map
 * @param i0 first-dimension index into the map
 * @param i1 second-dimension index into the map
 */
template< typename T >
GEOSX_HOST_DEVICE
inline T const & value( arrayView1d< arrayView1d< T const > const > const & map, localIndex const i0, localIndex const i1 )
{
  return map[i0][i1];
}

/**
 * @copydoc value(arrayView1d< arrayView1d< T const > const > const &, localIndex const, localIndex const)
 * @tparam USD unit-stride dimension of the map
 */
template< typename T, int USD >
GEOSX_HOST_DEVICE
inline T const & value( arrayView2d< T, USD > const & map, localIndex const i0, localIndex const i1 )
{
  return map( i0, i1 );
}

/**
 * @copydoc value(arrayView1d< arrayView1d< T const > const > const &, localIndex const, localIndex const)
 */
template< typename T >
GEOSX_HOST_DEVICE
inline T const & value( ArrayOfArraysView< T const > const & map, localIndex const i0, localIndex const i1 )
{
  return map( i0, i1 );
}

/**
 * @copydoc value(arrayView1d< arrayView1d< T const > const > const &, localIndex const, localIndex const)
 */
template< typename T >
GEOSX_HOST_DEVICE
inline T const & value( ArrayOfSetsView< T const > const & map, localIndex const i0, localIndex const i1 )
{
  return map( i0, i1 );
}

//////////////////////////////////////////////////////////////////////////////////

/**
 * @brief @return the total number of elements in the map
 * @tparam T type of map element
 * @param map reference to the map
 */
template< typename T >
GEOSX_HOST_DEVICE
inline localIndex numElements( arrayView1d< arrayView1d< T const > const > const & map )
{
  return map.size();
}

/**
 * @copydoc numElements(arrayView1d< arrayView1d< T const > const > const &)
 * @tparam USD unit-stride dimension of the map
 */
template< typename T, int USD >
GEOSX_HOST_DEVICE
inline localIndex numElements( arrayView2d< T, USD > const & map )
{
  return map.size();
}

/**
 * @copydoc numElements(arrayView1d< arrayView1d< T const > const > const &)
 */
template< typename T >
GEOSX_HOST_DEVICE
inline localIndex numElements( ArrayOfArraysView< T const > const & map )
{
  localIndex size = 0;
  for( localIndex i = 0; i < map.size(); ++i )
  {
    size += map.sizeOfArray( i );
  }
  return size;
}

/**
 * @copydoc numElements(arrayView1d< arrayView1d< T const > const > const &)
 */
template< typename T >
GEOSX_HOST_DEVICE
inline localIndex numElements( ArrayOfSetsView< T const > const & map )
{
  localIndex size = 0;
  for( localIndex i = 0; i < map.size(); ++i )
  {
    size += map.sizeOfSet( i );
  }
  return size;
}

} // namespace meshMapUtilities

} // namespace geosx

#endif //GEOSX_MESH_UTILITIES_MESHMAPUTILITIES_HPP
