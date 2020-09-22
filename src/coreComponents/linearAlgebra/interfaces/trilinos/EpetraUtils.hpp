/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
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

#ifdef GEOSX_USE_MPI
#include <Epetra_MpiComm.h>

/// Alias for specific EpetraComm implementation used.
using EpetraComm = Epetra_MpiComm;
#else
#include <Epetra_SerialComm.h>

/// Alias for specific EpetraComm implementation used.
using EpetraComm = Epetra_SerialComm;
#endif

namespace geosx
{

/**
 * @brief Converts a non-const array from GEOSX globalIndex type to Epetra long long
 * @param[in] index the input array
 * @return the converted array
 */
inline long long * toEpetraLongLong( globalIndex * const index )
{
  return reinterpret_cast< long long * >(index);
}

/**
 * @brief Converts a const array from GEOSX globalIndex type to Epetra long long
 * @param[in] index the input array
 * @return the converted array
 */
inline long long const * toEpetraLongLong( globalIndex const * const index )
{
  return reinterpret_cast< long long const * >(index);
}

}

#endif //GEOSX_LINEARALGEBRA_INTERFACES_EPETRAUTILS_HPP_
