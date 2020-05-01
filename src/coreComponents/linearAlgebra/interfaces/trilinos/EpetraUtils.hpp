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

/**
 * @file EpetraUtils.hpp
 */
#ifndef GEOSX_LINEARALGEBRA_INTERFACES_EPETRAUTILS_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_EPETRAUTILS_HPP_

#include "common/DataTypes.hpp"

namespace geosx
{

inline long long * toEpetraLongLong( globalIndex * const index )
{
  return reinterpret_cast< long long * >(index);
}

inline long long const * toEpetraLongLong( globalIndex const * const index )
{
  return reinterpret_cast< long long const * >(index);
}

}

#endif //GEOSX_LINEARALGEBRA_INTERFACES_EPETRAUTILS_HPP_
