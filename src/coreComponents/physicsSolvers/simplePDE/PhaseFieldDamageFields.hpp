/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file PhaseFieldDamageFields.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_SIMPLEPDE_PHASEFIELDDAMAGEFIELDS_HPP_
#define GEOS_PHYSICSSOLVERS_SIMPLEPDE_PHASEFIELDDAMAGEFIELDS_HPP_

#include "common/DataLayouts.hpp"
#include "mesh/MeshFields.hpp"

namespace geos
{
/**
 * A scope for field traits.
 */
namespace fields
{

namespace phaseField
{

DECLARE_FIELD( damage,
               "damage",
               array1d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Phase-field damage at the nodes" );

}

}

}

#endif // GEOS_PHYSICSSOLVERS_SIMPLEPDE_PHASEFIELDDAMAGEFIELDS_HPP_
