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

#ifndef GEOSX_DATAREPOSITORY_DYNAMICCASTS_HPP_
#define GEOSX_DATAREPOSITORY_DYNAMICCASTS_HPP_

#include "common/Logger.hpp"

#include <type_traits>

namespace geosx
{
namespace dataRepository
{

// The DynamicCast functions rely on `GroupDownCastHelper` to down cast instances with private inheritance to `Group`.
// And `GroupDownCastHelper` depends on `Group` which uses the DynamicCast functions...
template< class T >
class GroupDownCastHelper;

/**
 * A technical namespace dedicated to the dataRepository namespace components.
 */
namespace details
{

/**
 * @brief Perform a type cast of base to derived pointer.
 * @tparam NEW_TYPE      derived pointer type
 * @tparam EXISTING_TYPE base type
 * @param val            base pointer to cast
 * @return               pointer cast to derived type or @p nullptr
 *
 * With respect to standard dynamic_cast, this function adds some treatment to dataRepository instances
 * since it takes into account the potential private inheritance to `Group`.
 */
template< typename NEW_TYPE, typename EXISTING_TYPE >
NEW_TYPE DynamicCast( EXISTING_TYPE * const val )
{
  static_assert( std::is_pointer< NEW_TYPE >::value, "NEW_TYPE must be a pointer." );

  typedef const GroupDownCastHelper< std::remove_cv_t< std::remove_pointer_t< NEW_TYPE > > > GDCH;

  if( dynamic_cast< GDCH * >( val ) )
  {
    return reinterpret_cast< NEW_TYPE >( val );
  }
  else
  {
    return dynamic_cast< NEW_TYPE >( val );
  }
}

/**
 * @brief Perform a type cast of base to derived reference.
 * @tparam NEW_TYPE      derived reference type
 * @tparam EXISTING_TYPE base type
 * @param val            base reference to cast
 * @return               reference cast to derived type or @p nullptr
 *
 * With respect to standard dynamic_cast, this function adds some treatment to dataRepository instances
 * since it takes into account the potential private inheritance to `Group`.
 */
template< typename NEW_TYPE, typename EXISTING_TYPE >
NEW_TYPE DynamicCast( EXISTING_TYPE & val )
{
  static_assert( std::is_reference< NEW_TYPE >::value, "NEW_TYPE must be a reference." );

  using POINTER_TO_NEW_TYPE = std::remove_reference_t< NEW_TYPE > *;
  POINTER_TO_NEW_TYPE ptr = DynamicCast< POINTER_TO_NEW_TYPE >( &val );
  GEOSX_ERROR_IF( ptr == nullptr, "Cast failed." );

  return *ptr;
}

} /* end of namespace details */
} /* end namespace dataRepository */
} /* end namespace geosx */

#endif /* GEOSX_DATAREPOSITORY_DYNAMICCASTS_HPP_ */
