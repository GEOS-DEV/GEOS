/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOSX_DATAREPOSITORY_BUFFEROPSDEVICE_H_
#define GEOSX_DATAREPOSITORY_BUFFEROPSDEVICE_H_

#include "common/DataTypes.hpp"
#include "codingUtilities/Utilities.hpp"
#include "codingUtilities/static_if.hpp"
#include "codingUtilities/traits.hpp"
#include "cxx-utilities/src/IntegerConversion.hpp"
#include "BufferOps.hpp"

#include "rajaInterface/GEOS_RAJA_Interface.hpp"

#include <type_traits>

namespace geosx
{

namespace bufferOps
{

//------------------------------------------------------------------------------
template< typename POLICY, bool DO_PACKING, typename T, int NDIM, int UNIT_STRIDE_DIM >
typename std::enable_if< can_memcpy< T >, localIndex >::type
PackDevice( buffer_unit_type * & buffer,
            ArrayView< T, NDIM, UNIT_STRIDE_DIM > const & var );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, int NDIM, int UNIT_STRIDE_DIM >
GEOSX_HOST_DEVICE
typename std::enable_if< (!can_memcpy< T > || (NDIM > 1)), localIndex >::type
PackDevice( buffer_unit_type * & GEOSX_UNUSED_ARG( buffer ),
            ArraySlice< T, NDIM, UNIT_STRIDE_DIM > const & GEOSX_UNUSED_ARG( var ))
{
  GEOSX_ERROR( "Trying to pack data type ("<<typeid(T).name()<<") on device but type is not memcpy-able." );
  return 0;
}

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T >
GEOSX_HOST_DEVICE
typename std::enable_if< !can_memcpy< T >, localIndex >::type
PackDevice( buffer_unit_type * & GEOSX_UNUSED_ARG( buffer ),
            T & GEOSX_UNUSED_ARG( var ) )
{
  GEOSX_ERROR( "Trying to pack data type ("<<typeid(T).name()<<") on device but type is not packable." );
  return 0;
}

//------------------------------------------------------------------------------
template< typename POLICY, bool DO_PACKING, typename T, int NDIM, typename T_INDICES >
typename std::enable_if< can_memcpy< T >, localIndex >::type
PackByIndexDevice( buffer_unit_type * & buffer,
                   ArrayView< T, NDIM > const & var,
                   T_INDICES const & indices );

//------------------------------------------------------------------------------
template< typename POLICY, bool DO_PACKING, typename T, typename T_INDICES >
typename std::enable_if< (!is_packable_by_index< T > || !can_memcpy< T >), localIndex >::type
PackByIndexDevice( buffer_unit_type * & GEOSX_UNUSED_ARG( buffer ),
                   T const & GEOSX_UNUSED_ARG( var ),
                   T_INDICES const & GEOSX_UNUSED_ARG( indices ) )
{
  GEOSX_ERROR( "Trying to pack data type ("<<typeid(T).name()<<") on device but type is not packable by index." );
  return 0;
}

//------------------------------------------------------------------------------
template< typename POLICY, typename T, int NDIM, int UNIT_STRIDE_DIM >
typename std::enable_if< can_memcpy< T >, localIndex >::type
UnpackDevice( buffer_unit_type const * & buffer,
              ArrayView< T, NDIM, UNIT_STRIDE_DIM > const & var );

//------------------------------------------------------------------------------
template< typename POLICY, typename T >
typename std::enable_if< !can_memcpy< T >, localIndex >::type
UnpackDevice( buffer_unit_type const * & GEOSX_UNUSED_ARG( buffer ),
              T & GEOSX_UNUSED_ARG( var ) )
{
  GEOSX_ERROR( "Trying to unpack data type ("<<typeid(T).name()<<") on device but type is not packable." );
  return 0;
}

//------------------------------------------------------------------------------
template< typename POLICY, typename T, int NDIM, typename T_INDICES >
typename std::enable_if< can_memcpy< T >, localIndex >::type
UnpackByIndexDevice ( buffer_unit_type const * & buffer,
                      ArrayView< T, NDIM > const & var,
                      T_INDICES const & indices );

//------------------------------------------------------------------------------
template< typename POLICY, typename T, typename T_INDICES >
typename std::enable_if< (!is_packable_by_index< T > || !can_memcpy< T >), localIndex >::type
UnpackByIndexDevice( buffer_unit_type const * & GEOSX_UNUSED_ARG( buffer ),
                     T & GEOSX_UNUSED_ARG( var ),
                     T_INDICES const & GEOSX_UNUSED_ARG( indices ) )
{
  GEOSX_ERROR( "Trying to unpack data type ("<<typeid(T).name()<<") but type is not packable by index." );
  return 0;
}
}
}
#endif
