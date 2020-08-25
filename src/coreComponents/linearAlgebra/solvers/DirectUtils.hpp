/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file HypreUtils.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_DIRECTUTILS_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_DIRECTUTILS_HPP_

#include "common/DataTypes.hpp"

#ifdef GEOSX_USE_SUPERLU_DIST

#include <superlu_ddefs.h>

namespace geosx
{

/**
 * @brief Convert GEOSX global index value to SLUD int_t
 * @param index the input value
 * @return the converted value
 */
inline int_t toSuperlu_intT( globalIndex const index )
{
  return LvArray::integerConversion< int_t >( index );
}

/**
 * @brief Converts a non-const array from GEOSX globalIndex type to int_t
 * @param[in] index the input array
 * @return the converted array
 */
inline int_t * toSuperlu_intT( globalIndex * const index )
{
  return reinterpret_cast< int_t * >(index);
}

/**
 * @brief Converts a const array from GEOSX globalIndex type to int_t
 * @param[in] index the input array
 * @return the converted array
 */
inline int_t const * toSuperlu_intT( globalIndex const * const index )
{
  return reinterpret_cast< int_t const * >(index);
}

}

#endif

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_DIRECTUTILS_HPP_*/
