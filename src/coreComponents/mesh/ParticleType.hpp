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
 * @file ParticleType.hpp
 */

#ifndef GEOSX_MESH_PARTICLETYPE_HPP
#define GEOSX_MESH_PARTICLETYPE_HPP

#include "codingUtilities/EnumStrings.hpp"

namespace geos
{

/**
 * @brief Denotes type of particle shape/interpolation scheme
 */
enum class ParticleType : integer
{
  SinglePoint,          ///< Single-point (delta dirac characteristic function)
  CPDI,                 ///< Convected particle domain interpolation, parallelepiped domain
  CPTI,                 ///< Convected particle tetrahedral-domain interpolation, tet domain
  CPDI2                 ///< "2nd-order" CPDI, hexahedral domain
};

/// Strings for ParticleType
ENUM_STRINGS( ParticleType,
              "SinglePoint",
              "CPDI",
              "CPTI",
              "CPDI2" );

} // namespace geos

#endif //GEOSX_MESH_PARTICLETYPE_HPP
