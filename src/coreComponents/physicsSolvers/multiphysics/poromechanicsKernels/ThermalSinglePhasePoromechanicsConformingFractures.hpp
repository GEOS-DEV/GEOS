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
 * @file SinglePhasePoromechanicsConformingFractures.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_THERMALSINGLEPHASEPOROMECHANICSCONFORMINGFRACTURES_HPP
#define GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_THERMALSINGLEPHASEPOROMECHANICSCONFORMINGFRACTURES_HPP

#include "physicsSolvers/multiphysics/poromechanicsKernels/SinglePhasePoromechanicsConformingFractures.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/ThermalSinglePhasePoromechanicsConformingFracturesKernelBase.hpp"

namespace geos
{

namespace thermalSinglePhasePoromechanicsConformingFracturesKernels
{

template< integer NUM_EQN, integer NUM_DOF >
using ConnectorBasedAssemblyKernel =
  thermalConformingFracturesKernels::ThermalConformingFracturesConnectorBasedAssemblyKernel< NUM_EQN, NUM_DOF,
                                                                                             singlePhasePoromechanicsConformingFracturesKernels::ConnectorBasedAssemblyKernel >;

using ConnectorBasedAssemblyKernelFactory =
  thermalConformingFracturesKernels::ThermalConformingFracturesConnectorBasedAssemblyKernelFactory<
    singlePhasePoromechanicsConformingFracturesKernels::ConnectorBasedAssemblyKernel >;

} // namespace thermalSinglePhasePoromechanicsConformingFracturesKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_MULTIPHYSICS_POROMECHANICSKERNELS_THERMALSINGLEPHASEPOROMECHANICSCONFORMINGFRACTURES_HPP
