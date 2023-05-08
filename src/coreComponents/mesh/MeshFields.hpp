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
 * @file MeshFields.hpp
 */

#ifndef GEOS_MESH_FIELDS_HPP_
#define GEOS_MESH_FIELDS_HPP_

#include "codingUtilities/traits.hpp"
#include "dataRepository/RestartFlags.hpp"
#include "common/DataTypes.hpp"

/**
 * @brief Generates a traits struct.
 * @param NAME Name of the traits struct.
 * @param KEY The string literal that will be used as the key to register and
 *   lookup the data in the repository.
 * @param TYPE The type of data that will be registered.
 * @param DEFAULT The default value for the data.
 * @param PLOTLEVEL The default plot level for the wrapper that will contain the data.
 * @param RESTARTFLAG The default restart flag for the wrapper that contains the data.
 * @param DESCRIPTION A string literal that contains a description of the data for
 *   use in sphinx documentation.
 */
#define DECLARE_FIELD( NAME, \
                       KEY, \
                       TYPE, \
                       DEFAULT, \
                       PLOTLEVEL, \
                       RESTARTFLAG, \
                       DESCRIPTION ) \
/** @struct NAME */ \
/** @brief Trait struct for NAME data */ \
  struct NAME \
  { \
    /** @brief @return The key for registration with the data repository. */ \
    static constexpr char const * key() \
    { return KEY; } \
    /** The actual type to be registered. */ \
    using type = TYPE; \
    /** The template type T for registration of a container<T>. */ \
    using dataType = ::geos::fields::internal::typeHelper_t< TYPE >; \
    /** @brief @return The default data value for NAME. */ \
    static constexpr dataType defaultValue() \
    { return DEFAULT; } \
    /** The default dataRepository::PlotLevel for NAME. */ \
    static constexpr dataRepository::PlotLevel plotLevel = dataRepository::PlotLevel::PLOTLEVEL; \
    /** The default dataRepository::RestartFlags for NAME. */ \
    static constexpr dataRepository::RestartFlags restartFlag = dataRepository::RestartFlags::RESTARTFLAG; \
    /** Description of the NAME data for use in sphinx documentation */ \
    static constexpr char const * description = DESCRIPTION; \
  }

namespace geos
{
/**
 * A scope for field traits.
 */
namespace fields
{

namespace internal
{
template< typename TYPE, bool HAS_TYPE = std::enable_if_t< traits::HasAlias_value_type< TYPE >, std::true_type >::value >
struct typeHelper
{
  using type = typename TYPE::value_type;
};

template< typename TYPE >
struct typeHelper< TYPE, false >
{
  using type = TYPE;  //typename std::enable_if< std::is_fundamental<TYPE>::value, TYPE>::type;
};

template< typename T >
using typeHelper_t = typename typeHelper< T >::type;
}

DECLARE_FIELD( StructuredIndex,
               "structuredIndex",
               array2d< integer >,
               -1,
               NOPLOT,
               NO_WRITE,
               "Structured cell index provided by mesh generator" );

DECLARE_FIELD( ghostRank,
               "ghostRank",
               array1d< integer >,
               -2,
               LEVEL_0,
               WRITE_AND_READ,
               "Ghost rank." );

DECLARE_FIELD( elementVolume,
               "elementVolume",
               array1d< real64 >,
               0,
               LEVEL_1,
               WRITE_AND_READ,
               "Element volume." );

DECLARE_FIELD( elementAperture,
               "elementAperture",
               array1d< real64 >,
               1e-5,
               LEVEL_0,
               WRITE_AND_READ,
               "Element aperture." );

DECLARE_FIELD( parentIndex,
               "parentIndex",
               array1d< localIndex >,
               -1,
               LEVEL_2,
               WRITE_AND_READ,
               "Index of parent within the mesh object it is registered on." );

DECLARE_FIELD( parentEdgeIndex,
               "parentEdgeIndex",
               array1d< localIndex >,
               -1,
               LEVEL_2,
               WRITE_AND_READ,
               "Index of parent edge within the mesh object it is registered on." );

DECLARE_FIELD( childIndex,
               "childIndex",
               array1d< localIndex >,
               -1,
               LEVEL_2,
               WRITE_AND_READ,
               "Index of child within the mesh object it is registered on." );

DECLARE_FIELD( ruptureTime,
               "ruptureTime",
               array1d< real64 >,
               1.0e9,
               LEVEL_0,
               WRITE_AND_READ,
               "Time that the object was ruptured/split." );


} // namespace fields
} // namespace geos

#endif /* GEOS_MESH_FIELDS_HPP_ */
