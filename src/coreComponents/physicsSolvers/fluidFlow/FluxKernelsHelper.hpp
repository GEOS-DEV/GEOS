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
 * @file FluxKernelsHelper.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_FLUXKERNELSHELPER_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_FLUXKERNELSHELPER_HPP

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"


namespace geosx
{

namespace FluxKernelsHelper
{

template< typename VIEWTYPE >
 using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

GEOSX_HOST_DEVICE
void computeSinglePhaseFlux();


} // namespace FluxKernelsHelper

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_FLUXKERNELSHELPER_HPP
