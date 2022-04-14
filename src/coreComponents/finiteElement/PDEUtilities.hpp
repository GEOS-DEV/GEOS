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
 * @file PDEUtilities.hpp
 */

#ifndef GEOSX_FINITEELEMENT_PDEUTILITIES_HPP_
#define GEOSX_FINITEELEMENT_PDEUTILITIES_HPP_

#include "codingUtilities/EnumStrings.hpp"

namespace geosx
{

namespace PDEUtilities
{

/**
 * @brief Differential operator type.
 */
enum class DifferentialOperator : integer
{
  Divergence,
  Gradient,
  Identity,
  SymmetricGradient
};

/**
 * @brief Function space type.
 */
enum class FunctionSpace : integer
{
  invalid,
  P0,
  H1,
  H1vector
};

/// Declare strings associated with enumeration values.
ENUM_STRINGS( FunctionSpace,
              "invalid",
              "P0",
              "H1",
              "H1vector" );

} // namespace PDEUtilities

} // namespace geosx

#endif //GEOSX_FINITEELEMENT_PDEUTILITIES_HPP_
