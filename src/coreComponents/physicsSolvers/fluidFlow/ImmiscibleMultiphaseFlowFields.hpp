/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CompositionalMultiphaseBaseFields.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_IMMISCIBLEMULTIPHASEFLOWFIELDS_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_IMMISCIBLEMULTIPHASEFLOWFIELDS_HPP_

#include "common/DataLayouts.hpp"
#include "mesh/MeshFields.hpp"

namespace geos
{
/**
 * A scope for field traits.
 */
namespace fields
{

namespace immiscibleMultiphaseFlow
{

using array2dLayoutPhase = array2d< real64, immiscibleFlow::LAYOUT_PHASE >;
using array3dLayoutPhase_dS = array3d< real64, immiscibleFlow::LAYOUT_PHASE_DS >;

DECLARE_FIELD( phaseVolumeFraction,
               "phaseVolumeFraction",
               array2dLayoutPhase,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Phase volume fraction" );

DECLARE_FIELD( phaseVolumeFraction_n,
               "phaseVolumeFraction_n",
               array2dLayoutPhase,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Phase volume fraction at the previous converged time step" );

DECLARE_FIELD( phaseMass,
               "phaseMass",
               array2dLayoutPhase,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Phase mass" );

DECLARE_FIELD( phaseMass_n,
               "phaseMass_n",
               array2dLayoutPhase,
               0,
               NOPLOT,
               WRITE_AND_READ,
               "Phase mass at the previous converged time step" );

DECLARE_FIELD( phaseMobility,
               "phaseMobility",
               array2dLayoutPhase,
               0,
               LEVEL_0,
               WRITE_AND_READ,
               "Phase mobility" );

DECLARE_FIELD( dPhaseVolumeFraction,
               "dPhaseVolumeFraction",
               array3dLayoutPhase_dS,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivative of phase volume fraction with respect to pressure, temperature, global component density" );

DECLARE_FIELD( dPhaseMobility,
               "dPhaseMobility",
               array3dLayoutPhase_dS,
               0,
               NOPLOT,
               NO_WRITE,
               "Derivative of phase mobility with respect to pressure, temperature, global component density" );


DECLARE_FIELD( phaseOutflux,
               "phaseOutflux",
               array2dLayoutPhase,
               0,
               NOPLOT,
               NO_WRITE,
               "Phase outflux" );

DECLARE_FIELD( phaseCFLNumber,
               "phaseCFLNumber",
               array1d< real64 >,
               0,
               LEVEL_0,
               NO_WRITE,
               "Phase CFL number" );

DECLARE_FIELD( phaseDensity,
               "phaseDensity",
               array2dLayoutPhase,
               1000.0,
               NOPLOT,
               WRITE_AND_READ,
               "Phase density" );


DECLARE_FIELD( dPhaseDensity,
               "dPhaseDensity",
               array2dLayoutPhase,
               0.0,
               NOPLOT,
               WRITE_AND_READ,
               "Phase density derivative with respect to P" );

DECLARE_FIELD( phaseViscosity,
               "phaseViscosity",
               array2dLayoutPhase,
               1.0e-4,
               NOPLOT,
               WRITE_AND_READ,
               "Phase density" );


DECLARE_FIELD( dPhaseViscosity,
               "dPhaseViscosity",
               array2dLayoutPhase,
               0.0,
               NOPLOT,
               WRITE_AND_READ,
               "Phase density derivative with respect to P" );


}
}

}



#endif // GEOS_PHYSICSSOLVERS_FLUIDFLOW_IMMISCIBLEMULTIPHASEFLOWFIELDS_HPP_
