/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file TwoPhaseHybridFVMKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_TWOPHASEBASEKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_TWOPHASEBASEKERNELS_HPP

#include "common/DataTypes.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "physicsSolvers/fluidFlow/TwoPhaseBase.hpp"

namespace geosx
{

namespace TwoPhaseBaseKernels
{

/******************************** PhaseMobilityKernel ********************************/

struct PhaseMobilityKernel
{

  static void
  Compute( arraySlice1d< real64 const > phaseDens,
           arraySlice1d< real64 const > dPhaseDens_dPres,
           arraySlice1d< real64 const > phaseVisc,
           arraySlice1d< real64 const > dPhaseVisc_dPres,
           arraySlice1d< real64 const > phaseRelPerm,
           arraySlice2d< real64 const > dPhaseRelPerm_dSat,
           arraySlice1d< real64 > phaseMob,
           arraySlice1d< real64 > dPhaseMob_dPres,
           arraySlice1d< real64 > dPhaseMob_dSat );

  static void
  Launch( localIndex const size,
          arrayView3d< real64 const > const & phaseDens,
          arrayView3d< real64 const > const & dPhaseDens_dPres,
          arrayView3d< real64 const > const & phaseVisc,
          arrayView3d< real64 const > const & dPhaseVisc_dPres,
          arrayView3d< real64 const > const & phaseRelPerm,
          arrayView4d< real64 const > const & dPhaseRelPerm_dSat,
          arrayView2d< real64 > const & phaseMob,
          arrayView2d< real64 > const & dPhaseMob_dPres,
          arrayView2d< real64 > const & dPhaseMob_dSat );

};

/******************************** AccumulationKernel ********************************/

struct AccumulationKernel
{

  static void
  Compute( real64 const & volume,
           real64 const & porosityOld,
           real64 const & porosityRef,
           real64 const & pvMult,
           real64 const & dPvMult_dPres,
           arraySlice1d< real64 const > const phaseSat,
           arraySlice1d< real64 const > const dPhaseSat,
           arraySlice1d< real64 const > const phaseDensOld,
           arraySlice1d< real64 const > const phaseDens,
           arraySlice1d< real64 const > const dPhaseDens_dPres,
           arraySlice1d< real64 > const & localAccum,
           arraySlice2d< real64 > const & localAccumJacobian );

};


} // namespace TwoPhaseBaseKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_TWOPHASEBASEKERNELS_HPP
