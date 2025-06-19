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
 * @file CompositionalMultiphaseBaseFields.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEBASEFIELDS_HPP_
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEBASEFIELDS_HPP_

#include "common/DataLayouts.hpp"
#include "mesh/MeshFields.hpp"

namespace geos
{
/**
 * A scope for field traits.
 */
namespace fields
{

namespace flow
{

using array2dLayoutPhase = array2d< real64, compflow::LAYOUT_PHASE >;
using array3dLayoutPhase_dC = array3d< real64, compflow::LAYOUT_PHASE_DC >;
using array2dLayoutComp = array2d< real64, compflow::LAYOUT_COMP >;
using array3dLayoutComp_dC = array3d< real64, compflow::LAYOUT_COMP_DC >;
using array3dLayoutPhaseComp = array3d< real64, compflow::LAYOUT_PHASE_COMP >;

DECLARE_FIELD_WITH_NAMESPACE( globalCompDensity,
                              "globalCompDensity",
                              array2dLayoutComp,
                              0,
                              LEVEL_0,
                              WRITE_AND_READ,
                              "Global component density",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( globalCompDensity_n,
                              "globalCompDensity_n",
                              array2dLayoutComp,
                              0,
                              NOPLOT,
                              NO_WRITE,
                              "Global component density updates at the previous converged time step",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( globalCompDensity_k,
                              "globalCompDensity_k",
                              array2dLayoutComp,
                              0,
                              NOPLOT,
                              NO_WRITE,
                              "Global component density updates at the previous sequential iteration",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( globalCompFraction,
                              "globalCompFraction",
                              array2dLayoutComp,
                              0,
                              LEVEL_0,
                              WRITE_AND_READ,
                              "Global component fraction",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( globalCompFraction_n,
                              "globalCompFraction_n",
                              array2dLayoutComp,
                              0,
                              NOPLOT,
                              NO_WRITE,
                              "Global component fraction at the previous converged time step",
                              "flow" );

// may be needed later for sequential poromechanics implementation
//DECLARE_FIELD_WITH_NAMESPACE( globalCompFraction_k,
//               "globalCompFraction_k",
//               array2dLayoutComp,
//               0,
//               NOPLOT,
//               NO_WRITE,
//               "Global component fraction updates at the previous sequential iteration" );

DECLARE_FIELD_WITH_NAMESPACE( faceGlobalCompFraction,
                              "faceGlobalCompFraction",
                              array2dLayoutComp,
                              0,
                              LEVEL_0,
                              WRITE_AND_READ,
                              "Face global component fraction",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( dGlobalCompFraction_dGlobalCompDensity,
                              "dGlobalCompFraction_dGlobalCompDensity",
                              array3dLayoutComp_dC,
                              0,
                              NOPLOT,
                              NO_WRITE,
                              "Derivative of global component fraction with respect to component density",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( phaseVolumeFraction,
                              "phaseVolumeFraction",
                              array2dLayoutPhase,
                              0,
                              LEVEL_0,
                              WRITE_AND_READ,
                              "Phase volume fraction",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( dPhaseVolumeFraction,
                              "dPhaseVolumeFraction",
                              array3dLayoutPhase_dC,
                              0,
                              NOPLOT,
                              NO_WRITE,
                              "Derivative of phase volume fraction with respect to pressure, temperature, global component density",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( phaseMobility,
                              "phaseMobility",
                              array2dLayoutPhase,
                              0,
                              LEVEL_0,
                              WRITE_AND_READ,
                              "Phase mobility",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( dPhaseMobility,
                              "dPhaseMobility",
                              array3dLayoutPhase_dC,
                              0,
                              NOPLOT,
                              NO_WRITE,
                              "Derivative of phase volume fraction with respect to pressure, temperature, global component density",
                              "flow" );

// this is needed for time step selector
DECLARE_FIELD_WITH_NAMESPACE( phaseVolumeFraction_n,
                              "phaseVolumeFraction_n",
                              array2dLayoutPhase,
                              0,
                              NOPLOT,
                              WRITE_AND_READ,
                              "Phase volume fraction at the previous converged time step",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( phaseOutflux,
                              "phaseOutflux",
                              array2dLayoutPhase,
                              0,
                              NOPLOT,
                              NO_WRITE,
                              "Phase outflux",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( componentOutflux,
                              "componentOutflux",
                              array2dLayoutComp,
                              0,
                              NOPLOT,
                              NO_WRITE,
                              "Component outflux",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( phaseCFLNumber,
                              "phaseCFLNumber",
                              array1d< real64 >,
                              0,
                              LEVEL_0,
                              NO_WRITE,
                              "Phase CFL number",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( componentCFLNumber,
                              "componentCFLNumber",
                              array1d< real64 >,
                              0,
                              LEVEL_0,
                              NO_WRITE,
                              "Component CFL number",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( globalCompDensityScalingFactor,
                              "globalCompDensityScalingFactor",
                              array1d< real64 >,
                              1,
                              NOPLOT,
                              NO_WRITE,
                              "Scaling factors for global component densities",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( globalCompFractionScalingFactor,
                              "globalCompFractionScalingFactor",
                              array1d< real64 >,
                              1,
                              NOPLOT,
                              NO_WRITE,
                              "Scaling factors for global component fractions",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( compAmount,
                              "compAmount",
                              array2dLayoutComp,
                              0,
                              LEVEL_0,
                              WRITE_AND_READ,
                              "Component amount",
                              "flow" );

DECLARE_FIELD_WITH_NAMESPACE( compAmount_n,
                              "compAmount_n",
                              array2dLayoutComp,
                              0,
                              LEVEL_0,
                              WRITE_AND_READ,
                              "Component amount at the previous converged time step",
                              "flow" );

}

}

}

#endif // GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEBASEFIELDS_HPP_
