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
 * @file DispersionFields.hpp
 */

#ifndef GEOS_CONSTITUTIVE_DISPERSION_DISPERSIONFIELDS_HPP_
#define GEOS_CONSTITUTIVE_DISPERSION_DISPERSIONFIELDS_HPP_

#include "mesh/MeshFields.hpp"
#include "constitutive/dispersion/Layout.hpp"

namespace geos
{

namespace fields
{

namespace dispersion
{

using array3dLayoutPhase = array3d< real64, constitutive::dispersion::LAYOUT_PHASE_VELOCITY >;

DECLARE_FIELD( dispersivity,
               "dispersivity",
               array4d< real64 >,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Dispersivity" );

DECLARE_FIELD( phaseVelocity,
               "cellCenterPhaseVelocity",
               array3dLayoutPhase,
               1,
               LEVEL_0,
               NO_WRITE,
               "Phase velocities reconstructed at cell center" );


}

}

}

#endif // GEOS_CONSTITUTIVE_DISPERSION_DISPERSIONFIELDS_HPP_
