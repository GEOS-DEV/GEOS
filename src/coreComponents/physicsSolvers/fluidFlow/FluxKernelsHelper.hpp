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
#include "mesh/ElementRegionManager.hpp"

namespace geosx
{

namespace FluxKernelsHelper
{

template< typename VIEWTYPE >
 using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

GEOSX_HOST_DEVICE
void computeSinglePhaseFlux( arraySlice1d< localIndex const > const & seri,
                             arraySlice1d< localIndex const > const & sesri,
                             arraySlice1d< localIndex const > const & sei,
                             real64 const ( &transmissibility )[2],
                             real64 const ( &dTrans_dPres )[2],
                             real64 const ( &dTrans_dAper )[2],
                             ElementViewConst< arrayView1d< real64 const > > const & pres,
                             ElementViewConst< arrayView1d< real64 const > > const & dPres,
                             ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                             ElementViewConst< arrayView2d< real64 const > > const & dens,
                             ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
                             ElementViewConst< arrayView1d< real64 const > > const & mob,
                             ElementViewConst< arrayView1d< real64 const > > const & dMob_dPres,
                             real64 & fluxVal,
                             real64 (&dFlux_dP)[2],
                             real64 (&dFlux_dTrans)[2] );


} // namespace FluxKernelsHelper

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_FLUXKERNELSHELPER_HPP
