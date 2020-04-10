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
 * @file TwoPhaseBaseKernels.cpp
 */

#include "TwoPhaseBaseKernels.hpp"
#include "TwoPhaseBase.hpp"

namespace geosx
{

namespace TwoPhaseBaseKernels
{

void PhaseMobilityKernel::Compute( arraySlice1d< real64 const > phaseDens,
                                   arraySlice1d< real64 const > dPhaseDens_dPres,
                                   arraySlice1d< real64 const > phaseVisc,
                                   arraySlice1d< real64 const > dPhaseVisc_dPres,
                                   arraySlice1d< real64 const > phaseRelPerm,
                                   arraySlice2d< real64 const > dPhaseRelPerm_dSat,
                                   arraySlice1d< real64 > phaseMob,
                                   arraySlice1d< real64 > dPhaseMob_dPres,
                                   arraySlice1d< real64 > dPhaseMob_dSat )
{
  localIndex constexpr numPhases  = TwoPhaseBase::NUM_PHASES;
  for( localIndex ip = 0; ip < numPhases; ++ip )
  {
    real64 const relPermOverVisc        = phaseRelPerm[ip] / phaseVisc[ip];
    real64 const dRelPermOverVisc_dPres = -phaseRelPerm[ip] * dPhaseVisc_dPres[ip]
                                          / (phaseVisc[ip] * phaseVisc[ip]);
    real64 const dRelPermOverVisc_dSat  = dPhaseRelPerm_dSat[ip][ip] / phaseVisc[ip];

    phaseMob[ip]        = phaseDens[ip] * relPermOverVisc;
    dPhaseMob_dPres[ip] = phaseDens[ip] * dRelPermOverVisc_dPres
                          + dPhaseDens_dPres[ip] * relPermOverVisc;
    dPhaseMob_dSat[ip]  = phaseDens[ip] * dRelPermOverVisc_dSat;
  }
  dPhaseMob_dSat[1] *= -1; // we assume that the first saturation is the primary variable
}

void PhaseMobilityKernel::Launch( localIndex const size,
                                  arrayView3d< real64 const > const & phaseDens,
                                  arrayView3d< real64 const > const & dPhaseDens_dPres,
                                  arrayView3d< real64 const > const & phaseVisc,
                                  arrayView3d< real64 const > const & dPhaseVisc_dPres,
                                  arrayView3d< real64 const > const & phaseRelPerm,
                                  arrayView4d< real64 const > const & dPhaseRelPerm_dSat,
                                  arrayView2d< real64 > const & phaseMob,
                                  arrayView2d< real64 > const & dPhaseMob_dPres,
                                  arrayView2d< real64 > const & dPhaseMob_dSat )
{
  forAll< serialPolicy >( size, [=] ( localIndex const a )
  {
    Compute( phaseDens[a][0],
             dPhaseDens_dPres[a][0],
             phaseVisc[a][0],
             dPhaseVisc_dPres[a][0],
             phaseRelPerm[a][0],
             dPhaseRelPerm_dSat[a][0],
             phaseMob[a],
             dPhaseMob_dPres[a],
             dPhaseMob_dSat[a] );
  } );
}

void AccumulationKernel::Compute( real64 const & volume,
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
                                  arraySlice2d< real64 > const & localAccumJacobian )
{
  localIndex constexpr numPhases  = TwoPhaseBase::NUM_PHASES;

  localIndex constexpr dp = TwoPhaseBase::ColOffset::DPRES;
  localIndex constexpr dS = TwoPhaseBase::ColOffset::DSAT;

  real64 const poreVolNew        = volume * porosityRef * pvMult;
  real64 const dPoreVolNew_dPres = volume * porosityRef * dPvMult_dPres;
  real64 const poreVolOld        = volume * porosityOld;

  stackArray1d< real64, numPhases > phaseSatNew( numPhases );
  stackArray1d< real64, numPhases > dPhaseSatNew_dS( numPhases );

  for( localIndex ip = 0; ip < numPhases; ++ip )
  {
    phaseSatNew[ip] = phaseSat[ip] + dPhaseSat[ip];
  }
  dPhaseSatNew_dS[0] = 1;
  dPhaseSatNew_dS[1] = -1;

  for( localIndex ip = 0; ip < numPhases; ++ip )
  {
    // residual
    localAccum[ip]  = poreVolNew * phaseDens[ip] * phaseSatNew[ip];
    localAccum[ip] -= poreVolOld * phaseDensOld[ip] * phaseSat[ip];

    // jacobian
    localAccumJacobian[ip][dp] = ( dPoreVolNew_dPres * phaseDens[ip]
                                   + poreVolNew * dPhaseDens_dPres[ip] ) * phaseSatNew[ip];
    localAccumJacobian[ip][dS] = poreVolNew * phaseDens[ip] * dPhaseSatNew_dS[ip];
  }
}


} // namespace TwoPhaseBaseKernels

} // namespace geosx
