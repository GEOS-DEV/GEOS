/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file CompositionalMultiphaseFlowKernels.hpp
 */

#ifndef GEOSX_COMPOSITIONALMULTIPHASEFLOWKERNELS_HPP
#define GEOSX_COMPOSITIONALMULTIPHASEFLOWKERNELS_HPP

#include "common/DataTypes.hpp"
#include "constitutive/Fluid/MultiFluidBase.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

namespace geosx
{

namespace CompositionalMultiphaseFlowKernels
{

/******************************** ComponentFractionKernel ********************************/

/**
 * @brief Functions to compute component fractions from global component densities (mass or molar)
 */
struct ComponentFractionKernel
{

  template<localIndex NC>
  static inline RAJA_HOST_DEVICE void
  Compute( arraySlice1d<real64 const> compDens,
           arraySlice1d<real64 const> dCompDens,
           arraySlice1d<real64> compFrac,
           arraySlice2d<real64> dCompFrac_dCompDens );

  static inline RAJA_HOST_DEVICE void
  Compute( localIndex NC,
           arraySlice1d<real64 const> compDens,
           arraySlice1d<real64 const> dCompDens,
           arraySlice1d<real64> compFrac,
           arraySlice2d<real64> dCompFrac_dCompDens );

  template<localIndex NC>
  static void Launch( localIndex begin, localIndex end,
                      arrayView2d<real64 const> const & compDens,
                      arrayView2d<real64 const> const & dCompDens,
                      arrayView2d<real64> const & compFrac,
                      arrayView3d<real64> const & dCompFrac_dCompDens );

  static void Launch( localIndex NC,
                      localIndex begin, localIndex end,
                      arrayView2d<real64 const> const & compDens,
                      arrayView2d<real64 const> const & dCompDens,
                      arrayView2d<real64> const & compFrac,
                      arrayView3d<real64> const & dCompFrac_dCompDens );

  template<localIndex NC>
  static void Launch( set<localIndex> const & targetSet,
                      arrayView2d<real64 const> const & compDens,
                      arrayView2d<real64 const> const & dCompDens,
                      arrayView2d<real64> const & compFrac,
                      arrayView3d<real64> const & dCompFrac_dCompDens );

  static void Launch( localIndex NC,
                      set<localIndex> const & targetSet,
                      arrayView2d<real64 const> const & compDens,
                      arrayView2d<real64 const> const & dCompDens,
                      arrayView2d<real64> const & compFrac,
                      arrayView3d<real64> const & dCompFrac_dCompDens );

};

#define INST_ComponentFractionKernel(NC) \
extern template \
void ComponentFractionKernel::Launch<NC>( localIndex begin, localIndex end, \
                                          arrayView2d<real64 const> const & compDens, \
                                          arrayView2d<real64 const> const & dCompDens, \
                                          arrayView2d<real64> const & compFrac, \
                                          arrayView3d<real64> const & dCompFrac_dCompDens ); \
extern template \
void ComponentFractionKernel::Launch<NC>( set<localIndex> const & targetSet, \
                                          arrayView2d<real64 const> const & compDens, \
                                          arrayView2d<real64 const> const & dCompDens, \
                                          arrayView2d<real64> const & compFrac, \
                                          arrayView3d<real64> const & dCompFrac_dCompDens )

INST_ComponentFractionKernel(1);
INST_ComponentFractionKernel(2);
INST_ComponentFractionKernel(3);
INST_ComponentFractionKernel(4);
INST_ComponentFractionKernel(5);

#undef INST_ComponentFractionKernel

/******************************** PhaseVolumeFractionKernel ********************************/

/**
 * @brief Functions to compute phase volume fractions (saturations) and derivatives
 */
struct PhaseVolumeFractionKernel
{

  template<localIndex NC, localIndex NP>
  static inline RAJA_HOST_DEVICE void
  Compute( arraySlice1d<real64 const> compDens,
           arraySlice1d<real64 const> dCompDens,
           arraySlice2d<real64 const> dCompFrac_dCompDens,
           arraySlice1d<real64 const> phaseDens,
           arraySlice1d<real64 const> dPhaseDens_dPres,
           arraySlice2d<real64 const> dPhaseDens_dComp,
           arraySlice1d<real64 const> phaseFrac,
           arraySlice1d<real64 const> dPhaseFrac_dPres,
           arraySlice2d<real64 const> dPhaseFrac_dComp,
           arraySlice1d<real64> phaseVolFrac,
           arraySlice1d<real64> dPhaseVolFrac_dPres,
           arraySlice2d<real64> dPhaseVolFrac_dComp );

  static inline RAJA_HOST_DEVICE void
  Compute( localIndex NC, localIndex NP,
           arraySlice1d<real64 const> compDens,
           arraySlice1d<real64 const> dCompDens,
           arraySlice2d<real64 const> dCompFrac_dCompDens,
           arraySlice1d<real64 const> phaseDens,
           arraySlice1d<real64 const> dPhaseDens_dPres,
           arraySlice2d<real64 const> dPhaseDens_dComp,
           arraySlice1d<real64 const> phaseFrac,
           arraySlice1d<real64 const> dPhaseFrac_dPres,
           arraySlice2d<real64 const> dPhaseFrac_dComp,
           arraySlice1d<real64> phaseVolFrac,
           arraySlice1d<real64> dPhaseVolFrac_dPres,
           arraySlice2d<real64> dPhaseVolFrac_dComp );

  template<localIndex NC, localIndex NP>
  static void Launch( localIndex begin, localIndex end,
                      arrayView2d<real64 const> const & compDens,
                      arrayView2d<real64 const> const & dCompDens,
                      arrayView3d<real64 const> const & dCompFrac_dCompDens,
                      arrayView3d<real64 const> const & phaseDens,
                      arrayView3d<real64 const> const & dPhaseDens_dPres,
                      arrayView4d<real64 const> const & dPhaseDens_dComp,
                      arrayView3d<real64 const> const & phaseFrac,
                      arrayView3d<real64 const> const & dPhaseFrac_dPres,
                      arrayView4d<real64 const> const & dPhaseFrac_dComp,
                      arrayView2d<real64> const & phaseVolFrac,
                      arrayView2d<real64> const & dPhaseVolFrac_dPres,
                      arrayView3d<real64> const & dPhaseVolFrac_dComp );

  static void Launch( localIndex NC, localIndex NP,
                      localIndex begin, localIndex end,
                      arrayView2d<real64 const> const & compDens,
                      arrayView2d<real64 const> const & dCompDens,
                      arrayView3d<real64 const> const & dCompFrac_dCompDens,
                      arrayView3d<real64 const> const & phaseDens,
                      arrayView3d<real64 const> const & dPhaseDens_dPres,
                      arrayView4d<real64 const> const & dPhaseDens_dComp,
                      arrayView3d<real64 const> const & phaseFrac,
                      arrayView3d<real64 const> const & dPhaseFrac_dPres,
                      arrayView4d<real64 const> const & dPhaseFrac_dComp,
                      arrayView2d<real64> const & phaseVolFrac,
                      arrayView2d<real64> const & dPhaseVolFrac_dPres,
                      arrayView3d<real64> const & dPhaseVolFrac_dComp );

  template<localIndex NC, localIndex NP>
  static void Launch( set<localIndex> const & targetSet,
                      arrayView2d<real64 const> const & compDens,
                      arrayView2d<real64 const> const & dCompDens,
                      arrayView3d<real64 const> const & dCompFrac_dCompDens,
                      arrayView3d<real64 const> const & phaseDens,
                      arrayView3d<real64 const> const & dPhaseDens_dPres,
                      arrayView4d<real64 const> const & dPhaseDens_dComp,
                      arrayView3d<real64 const> const & phaseFrac,
                      arrayView3d<real64 const> const & dPhaseFrac_dPres,
                      arrayView4d<real64 const> const & dPhaseFrac_dComp,
                      arrayView2d<real64> const & phaseVolFrac,
                      arrayView2d<real64> const & dPhaseVolFrac_dPres,
                      arrayView3d<real64> const & dPhaseVolFrac_dComp );

  static void Launch( localIndex NC, localIndex NP,
                      set<localIndex> const & targetSet,
                      arrayView2d<real64 const> const & compDens,
                      arrayView2d<real64 const> const & dCompDens,
                      arrayView3d<real64 const> const & dCompFrac_dCompDens,
                      arrayView3d<real64 const> const & phaseDens,
                      arrayView3d<real64 const> const & dPhaseDens_dPres,
                      arrayView4d<real64 const> const & dPhaseDens_dComp,
                      arrayView3d<real64 const> const & phaseFrac,
                      arrayView3d<real64 const> const & dPhaseFrac_dPres,
                      arrayView4d<real64 const> const & dPhaseFrac_dComp,
                      arrayView2d<real64> const & phaseVolFrac,
                      arrayView2d<real64> const & dPhaseVolFrac_dPres,
                      arrayView3d<real64> const & dPhaseVolFrac_dComp );

};

#define INST_PhaseVolumeFractionKernel(NC,NP) \
extern template \
void PhaseVolumeFractionKernel::Launch<NC,NP>( localIndex begin, localIndex end, \
                                               arrayView2d<real64 const> const & compDens, \
                                               arrayView2d<real64 const> const & dCompDens, \
                                               arrayView3d<real64 const> const & dCompFrac_dCompDens, \
                                               arrayView3d<real64 const> const & phaseDens, \
                                               arrayView3d<real64 const> const & dPhaseDens_dPres, \
                                               arrayView4d<real64 const> const & dPhaseDens_dComp, \
                                               arrayView3d<real64 const> const & phaseFrac, \
                                               arrayView3d<real64 const> const & dPhaseFrac_dPres, \
                                               arrayView4d<real64 const> const & dPhaseFrac_dComp, \
                                               arrayView2d<real64> const & phaseVolFrac, \
                                               arrayView2d<real64> const & dPhaseVolFrac_dPres, \
                                               arrayView3d<real64> const & dPhaseVolFrac_dComp ); \
extern template \
void PhaseVolumeFractionKernel::Launch<NC,NP>( set<localIndex> const & targetSet, \
                                               arrayView2d<real64 const> const & compDens, \
                                               arrayView2d<real64 const> const & dCompDens, \
                                               arrayView3d<real64 const> const & dCompFrac_dCompDens, \
                                               arrayView3d<real64 const> const & phaseDens, \
                                               arrayView3d<real64 const> const & dPhaseDens_dPres, \
                                               arrayView4d<real64 const> const & dPhaseDens_dComp, \
                                               arrayView3d<real64 const> const & phaseFrac, \
                                               arrayView3d<real64 const> const & dPhaseFrac_dPres, \
                                               arrayView4d<real64 const> const & dPhaseFrac_dComp, \
                                               arrayView2d<real64> const & phaseVolFrac, \
                                               arrayView2d<real64> const & dPhaseVolFrac_dPres, \
                                               arrayView3d<real64> const & dPhaseVolFrac_dComp )

INST_PhaseVolumeFractionKernel(1,1);
INST_PhaseVolumeFractionKernel(2,1);
INST_PhaseVolumeFractionKernel(3,1);
INST_PhaseVolumeFractionKernel(4,1);
INST_PhaseVolumeFractionKernel(5,1);

INST_PhaseVolumeFractionKernel(1,2);
INST_PhaseVolumeFractionKernel(2,2);
INST_PhaseVolumeFractionKernel(3,2);
INST_PhaseVolumeFractionKernel(4,2);
INST_PhaseVolumeFractionKernel(5,2);

INST_PhaseVolumeFractionKernel(1,3);
INST_PhaseVolumeFractionKernel(2,3);
INST_PhaseVolumeFractionKernel(3,3);
INST_PhaseVolumeFractionKernel(4,3);
INST_PhaseVolumeFractionKernel(5,3);

#undef INST_PhaseVolumeFractionKernel

/******************************** PhaseMobilityKernel ********************************/

/**
 * @brief Functions to compute phase mobilities and derivatives from density, viscosity and relperm
 */
struct PhaseMobilityKernel
{

  template<localIndex NC, localIndex NP>
  static inline RAJA_HOST_DEVICE void
  Compute( arraySlice2d<real64 const> dCompFrac_dCompDens,
           arraySlice1d<real64 const> phaseDens,
           arraySlice1d<real64 const> dPhaseDens_dPres,
           arraySlice2d<real64 const> dPhaseDens_dComp,
           arraySlice1d<real64 const> phaseVisc,
           arraySlice1d<real64 const> dPhaseVisc_dPres,
           arraySlice2d<real64 const> dPhaseVisc_dComp,
           arraySlice1d<real64 const> phaseRelPerm,
           arraySlice2d<real64 const> dPhaseRelPerm_dPhaseVolFrac,
           arraySlice1d<real64 const> dPhaseVolFrac_dPres,
           arraySlice2d<real64 const> dPhaseVolFrac_dComp,
           arraySlice1d<real64> phaseMob,
           arraySlice1d<real64> dPhaseMob_dPres,
           arraySlice2d<real64> dPhaseMob_dComp );

  static inline RAJA_HOST_DEVICE void
  Compute( localIndex NC, localIndex NP,
           arraySlice2d<real64 const> dCompFrac_dCompDens,
           arraySlice1d<real64 const> phaseDens,
           arraySlice1d<real64 const> dPhaseDens_dPres,
           arraySlice2d<real64 const> dPhaseDens_dComp,
           arraySlice1d<real64 const> phaseVisc,
           arraySlice1d<real64 const> dPhaseVisc_dPres,
           arraySlice2d<real64 const> dPhaseVisc_dComp,
           arraySlice1d<real64 const> phaseRelPerm,
           arraySlice2d<real64 const> dPhaseRelPerm_dPhaseVolFrac,
           arraySlice1d<real64 const> dPhaseVolFrac_dPres,
           arraySlice2d<real64 const> dPhaseVolFrac_dComp,
           arraySlice1d<real64> phaseMob,
           arraySlice1d<real64> dPhaseMob_dPres,
           arraySlice2d<real64> dPhaseMob_dComp );

  template<localIndex NC, localIndex NP>
  static void Launch( localIndex begin, localIndex end,
                      arrayView3d<real64 const> const & dCompFrac_dCompDens,
                      arrayView3d<real64 const> const & phaseDens,
                      arrayView3d<real64 const> const & dPhaseDens_dPres,
                      arrayView4d<real64 const> const & dPhaseDens_dComp,
                      arrayView3d<real64 const> const & phaseVisc,
                      arrayView3d<real64 const> const & dPhaseVisc_dPres,
                      arrayView4d<real64 const> const & dPhaseVisc_dComp,
                      arrayView3d<real64 const> const & phaseRelPerm,
                      arrayView4d<real64 const> const & dPhaseRelPerm_dPhaseVolFrac,
                      arrayView2d<real64 const> const & dPhaseVolFrac_dPres,
                      arrayView3d<real64 const> const & dPhaseVolFrac_dComp,
                      arrayView2d<real64> const & phaseMob,
                      arrayView2d<real64> const & dPhaseMob_dPres,
                      arrayView3d<real64> const & dPhaseMob_dComp );

  static void Launch( localIndex NC, localIndex NP,
                      localIndex begin, localIndex end,
                      arrayView3d<real64 const> const & dCompFrac_dCompDens,
                      arrayView3d<real64 const> const & phaseDens,
                      arrayView3d<real64 const> const & dPhaseDens_dPres,
                      arrayView4d<real64 const> const & dPhaseDens_dComp,
                      arrayView3d<real64 const> const & phaseVisc,
                      arrayView3d<real64 const> const & dPhaseVisc_dPres,
                      arrayView4d<real64 const> const & dPhaseVisc_dComp,
                      arrayView3d<real64 const> const & phaseRelPerm,
                      arrayView4d<real64 const> const & dPhaseRelPerm_dPhaseVolFrac,
                      arrayView2d<real64 const> const & dPhaseVolFrac_dPres,
                      arrayView3d<real64 const> const & dPhaseVolFrac_dComp,
                      arrayView2d<real64> const & phaseMob,
                      arrayView2d<real64> const & dPhaseMob_dPres,
                      arrayView3d<real64> const & dPhaseMob_dComp );

  template<localIndex NC, localIndex NP>
  static void Launch( set<localIndex> const & targetSet,
                      arrayView3d<real64 const> const & dCompFrac_dCompDens,
                      arrayView3d<real64 const> const & phaseDens,
                      arrayView3d<real64 const> const & dPhaseDens_dPres,
                      arrayView4d<real64 const> const & dPhaseDens_dComp,
                      arrayView3d<real64 const> const & phaseVisc,
                      arrayView3d<real64 const> const & dPhaseVisc_dPres,
                      arrayView4d<real64 const> const & dPhaseVisc_dComp,
                      arrayView3d<real64 const> const & phaseRelPerm,
                      arrayView4d<real64 const> const & dPhaseRelPerm_dPhaseVolFrac,
                      arrayView2d<real64 const> const & dPhaseVolFrac_dPres,
                      arrayView3d<real64 const> const & dPhaseVolFrac_dComp,
                      arrayView2d<real64> const & phaseMob,
                      arrayView2d<real64> const & dPhaseMob_dPres,
                      arrayView3d<real64> const & dPhaseMob_dComp );

  static void Launch( localIndex NC, localIndex NP,
                      set<localIndex> const & targetSet,
                      arrayView3d<real64 const> const & dCompFrac_dCompDens,
                      arrayView3d<real64 const> const & phaseDens,
                      arrayView3d<real64 const> const & dPhaseDens_dPres,
                      arrayView4d<real64 const> const & dPhaseDens_dComp,
                      arrayView3d<real64 const> const & phaseVisc,
                      arrayView3d<real64 const> const & dPhaseVisc_dPres,
                      arrayView4d<real64 const> const & dPhaseVisc_dComp,
                      arrayView3d<real64 const> const & phaseRelPerm,
                      arrayView4d<real64 const> const & dPhaseRelPerm_dPhaseVolFrac,
                      arrayView2d<real64 const> const & dPhaseVolFrac_dPres,
                      arrayView3d<real64 const> const & dPhaseVolFrac_dComp,
                      arrayView2d<real64> const & phaseMob,
                      arrayView2d<real64> const & dPhaseMob_dPres,
                      arrayView3d<real64> const & dPhaseMob_dComp );

};

#define INST_PhaseMobilityKernel(NC,NP) \
extern template \
void PhaseMobilityKernel::Launch<NC,NP>( localIndex begin, localIndex end, \
                                         arrayView3d<real64 const> const & dCompFrac_dCompDens, \
                                         arrayView3d<real64 const> const & phaseDens, \
                                         arrayView3d<real64 const> const & dPhaseDens_dPres, \
                                         arrayView4d<real64 const> const & dPhaseDens_dComp, \
                                         arrayView3d<real64 const> const & phaseVisc, \
                                         arrayView3d<real64 const> const & dPhaseVisc_dPres, \
                                         arrayView4d<real64 const> const & dPhaseVisc_dComp, \
                                         arrayView3d<real64 const> const & phaseRelPerm, \
                                         arrayView4d<real64 const> const & dPhaseRelPerm_dPhaseVolFrac, \
                                         arrayView2d<real64 const> const & dPhaseVolFrac_dPres, \
                                         arrayView3d<real64 const> const & dPhaseVolFrac_dComp, \
                                         arrayView2d<real64> const & phaseMob, \
                                         arrayView2d<real64> const & dPhaseMob_dPres, \
                                         arrayView3d<real64> const & dPhaseMob_dComp ); \
extern template \
void PhaseMobilityKernel::Launch<NC,NP>( set<localIndex> const & targetSet, \
                                         arrayView3d<real64 const> const & dCompFrac_dCompDens, \
                                         arrayView3d<real64 const> const & phaseDens, \
                                         arrayView3d<real64 const> const & dPhaseDens_dPres, \
                                         arrayView4d<real64 const> const & dPhaseDens_dComp, \
                                         arrayView3d<real64 const> const & phaseVisc, \
                                         arrayView3d<real64 const> const & dPhaseVisc_dPres, \
                                         arrayView4d<real64 const> const & dPhaseVisc_dComp, \
                                         arrayView3d<real64 const> const & phaseRelPerm, \
                                         arrayView4d<real64 const> const & dPhaseRelPerm_dPhaseVolFrac, \
                                         arrayView2d<real64 const> const & dPhaseVolFrac_dPres, \
                                         arrayView3d<real64 const> const & dPhaseVolFrac_dComp, \
                                         arrayView2d<real64> const & phaseMob, \
                                         arrayView2d<real64> const & dPhaseMob_dPres, \
                                         arrayView3d<real64> const & dPhaseMob_dComp )

INST_PhaseMobilityKernel(1,1);
INST_PhaseMobilityKernel(2,1);
INST_PhaseMobilityKernel(3,1);
INST_PhaseMobilityKernel(4,1);
INST_PhaseMobilityKernel(5,1);

INST_PhaseMobilityKernel(1,2);
INST_PhaseMobilityKernel(2,2);
INST_PhaseMobilityKernel(3,2);
INST_PhaseMobilityKernel(4,2);
INST_PhaseMobilityKernel(5,2);

INST_PhaseMobilityKernel(1,3);
INST_PhaseMobilityKernel(2,3);
INST_PhaseMobilityKernel(3,3);
INST_PhaseMobilityKernel(4,3);
INST_PhaseMobilityKernel(5,3);

#undef INST_PhaseMobilityKernel

/******************************** AccumulationKernel ********************************/

/**
 * @brief Functions to assemble accumulation term contributions to residual and Jacobian
 */
struct AccumulationKernel
{

  static inline RAJA_HOST_DEVICE void
  Compute( localIndex const NC, localIndex const NP,
           real64 const & volume,
           real64 const & porosityOld,
           real64 const & porosityRef,
           real64 const & pvMult,
           real64 const & dPvMult_dPres,
           arraySlice2d<real64 const> const dCompFrac_dCompDens,
           arraySlice1d<real64 const> const phaseVolFracOld,
           arraySlice1d<real64 const> const phaseVolFrac,
           arraySlice1d<real64 const> const dPhaseVolFrac_dPres,
           arraySlice2d<real64 const> const dPhaseVolFrac_dCompDens,
           arraySlice1d<real64 const> const phaseDensOld,
           arraySlice1d<real64 const> const phaseDens,
           arraySlice1d<real64 const> const dPhaseDens_dPres,
           arraySlice2d<real64 const> const dPhaseDens_dComp,
           arraySlice2d<real64 const> const phaseCompFracOld,
           arraySlice2d<real64 const> const phaseCompFrac,
           arraySlice2d<real64 const> const dPhaseCompFrac_dPres,
           arraySlice3d<real64 const> const dPhaseCompFrac_dComp,
           arraySlice1d<real64> const & localAccum,
           arraySlice2d<real64> const & localAccumJacobian )
  {
    localIndex constexpr maxNumComp = constitutive::MultiFluidBase::MAX_NUM_COMPONENTS;

    // temporary work arrays
    stackArray1d<real64, maxNumComp> const dPhaseAmount_dC( NC );
    stackArray1d<real64, maxNumComp> const dPhaseCompFrac_dC( NC );

    // reset the local values
    for (localIndex ic = 0; ic < NC; ++ic)
    {
      localAccum[ic] = 0.0;
      for (localIndex jc = 0; jc < NC; ++jc)
      {
        localAccumJacobian[ic][jc] = 0.0;
      }
    }

    // compute fluid-independent (pore volume) part
    real64 const volNew = volume;
    real64 const volOld = volume;
    real64 const dVol_dP = 0.0; // used in poroelastic solver

    real64 const poroNew = porosityRef * pvMult;
    real64 const poroOld = porosityOld;
    real64 const dPoro_dP = porosityRef * dPvMult_dPres;

    real64 const poreVolNew = volNew * poroNew;
    real64 const poreVolOld = volOld * poroOld;
    real64 const dPoreVol_dP = dVol_dP * poroNew + volNew * dPoro_dP;

    // sum contributions to component accumulation from each phase
    for (localIndex ip = 0; ip < NP; ++ip)
    {
      real64 const phaseAmountNew = poreVolNew * phaseVolFrac[ip] * phaseDens[ip];
      real64 const phaseAmountOld = poreVolOld * phaseVolFracOld[ip] * phaseDensOld[ip];

      real64 const dPhaseAmount_dP = dPoreVol_dP * phaseVolFrac[ip] * phaseDens[ip]
                                    + poreVolNew * (dPhaseVolFrac_dPres[ip] * phaseDens[ip]
                                                         + phaseVolFrac[ip] * dPhaseDens_dPres[ip]);

      // assemble density dependence
      applyChainRule(NC, dCompFrac_dCompDens, dPhaseDens_dComp[ip], dPhaseAmount_dC);
      for (localIndex jc = 0; jc < NC; ++jc)
      {
        dPhaseAmount_dC[jc] = dPhaseAmount_dC[jc] * phaseVolFrac[ip]
                              + phaseDens[ip] * dPhaseVolFrac_dCompDens[ip][jc];
        dPhaseAmount_dC[jc] *= poreVolNew;
      }

      // ic - index of component whose conservation equation is assembled
      // (i.e. row number in local matrix)
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        real64 const phaseCompAmountNew = phaseAmountNew * phaseCompFrac[ip][ic];
        real64 const phaseCompAmountOld = phaseAmountOld * phaseCompFracOld[ip][ic];

        real64 const dPhaseCompAmount_dP = dPhaseAmount_dP * phaseCompFrac[ip][ic]
                                          + phaseAmountNew * dPhaseCompFrac_dPres[ip][ic];

        localAccum[ic] += phaseCompAmountNew - phaseCompAmountOld;
        localAccumJacobian[ic][0] += dPhaseCompAmount_dP;

        // jc - index of component w.r.t. whose compositional var the derivative is being taken
        // (i.e. col number in local matrix)

        // assemble phase composition dependence
        applyChainRule(NC, dCompFrac_dCompDens, dPhaseCompFrac_dComp[ip][ic], dPhaseCompFrac_dC);
        for (localIndex jc = 0; jc < NC; ++jc)
        {
          real64 const dPhaseCompAmount_dC = dPhaseCompFrac_dC[jc] * phaseAmountNew
                                           + phaseCompFrac[ip][ic] * dPhaseAmount_dC[jc];
          localAccumJacobian[ic][jc + 1] += dPhaseCompAmount_dC;
        }
      }
    }
  }

  template<localIndex NC, localIndex NP>
  static inline RAJA_HOST_DEVICE void
  Compute( real64 const & volume,
           real64 const & porosityOld,
           real64 const & porosityRef,
           real64 const & pvMult,
           real64 const & dPvMult_dPres,
           arraySlice2d<real64 const> const dCompFrac_dCompDens,
           arraySlice1d<real64 const> const phaseVolFracOld,
           arraySlice1d<real64 const> const phaseVolFrac,
           arraySlice1d<real64 const> const dPhaseVolFrac_dPres,
           arraySlice2d<real64 const> const dPhaseVolFrac_dCompDens,
           arraySlice1d<real64 const> const phaseDensOld,
           arraySlice1d<real64 const> const phaseDens,
           arraySlice1d<real64 const> const dPhaseDens_dPres,
           arraySlice2d<real64 const> const dPhaseDens_dComp,
           arraySlice2d<real64 const> const phaseCompFracOld,
           arraySlice2d<real64 const> const phaseCompFrac,
           arraySlice2d<real64 const> const dPhaseCompFrac_dPres,
           arraySlice3d<real64 const> const dPhaseCompFrac_dComp,
           arraySlice1d<real64> const & localAccum,
           arraySlice2d<real64> const & localAccumJacobian )
  {
    // temporary work arrays
    stackArray1d<real64, NC> dPhaseAmount_dC( NC );
    stackArray1d<real64, NC> dPhaseCompFrac_dC( NC );

    // reset the local values
    for (localIndex ic = 0; ic < NC; ++ic)
    {
      localAccum[ic] = 0.0;
      for (localIndex jc = 0; jc < NC; ++jc)
      {
        localAccumJacobian[ic][jc] = 0.0;
      }
    }

    // compute fluid-independent (pore volume) part
    real64 const volNew = volume;
    real64 const volOld = volume;
    real64 const dVol_dP = 0.0; // used in poroelastic solver

    real64 const poroNew = porosityRef * pvMult;
    real64 const poroOld = porosityOld;
    real64 const dPoro_dP = porosityRef * dPvMult_dPres;

    real64 const poreVolNew = volNew * poroNew;
    real64 const poreVolOld = volOld * poroOld;
    real64 const dPoreVol_dP = dVol_dP * poroNew + volNew * dPoro_dP;

    // sum contributions to component accumulation from each phase
    for (localIndex ip = 0; ip < NP; ++ip)
    {
      real64 const phaseAmountNew = poreVolNew * phaseVolFrac[ip] * phaseDens[ip];
      real64 const phaseAmountOld = poreVolOld * phaseVolFracOld[ip] * phaseDensOld[ip];

      real64 const dPhaseAmount_dP = dPoreVol_dP * phaseVolFrac[ip] * phaseDens[ip]
                                    + poreVolNew * (dPhaseVolFrac_dPres[ip] * phaseDens[ip]
                                                         + phaseVolFrac[ip] * dPhaseDens_dPres[ip]);

      // assemble density dependence
      applyChainRule(NC, dCompFrac_dCompDens, dPhaseDens_dComp[ip], dPhaseAmount_dC);
      for (localIndex jc = 0; jc < NC; ++jc)
      {
        dPhaseAmount_dC[jc] = dPhaseAmount_dC[jc] * phaseVolFrac[ip]
                                  + phaseDens[ip] * dPhaseVolFrac_dCompDens[ip][jc];
        dPhaseAmount_dC[jc] *= poreVolNew;
      }

      // ic - index of component whose conservation equation is assembled
      // (i.e. row number in local matrix)
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        real64 const phaseCompAmountNew = phaseAmountNew * phaseCompFrac[ip][ic];
        real64 const phaseCompAmountOld = phaseAmountOld * phaseCompFracOld[ip][ic];

        real64 const dPhaseCompAmount_dP = dPhaseAmount_dP * phaseCompFrac[ip][ic]
                                          + phaseAmountNew * dPhaseCompFrac_dPres[ip][ic];

        localAccum[ic] += phaseCompAmountNew - phaseCompAmountOld;
        localAccumJacobian[ic][0] += dPhaseCompAmount_dP;

        // jc - index of component w.r.t. whose compositional var the derivative is being taken
        // (i.e. col number in local matrix)

        // assemble phase composition dependence
        applyChainRule(NC, dCompFrac_dCompDens, dPhaseCompFrac_dComp[ip][ic], dPhaseCompFrac_dC);
        for (localIndex jc = 0; jc < NC; ++jc)
        {
          real64 const dPhaseCompAmount_dC = dPhaseCompFrac_dC[jc] * phaseAmountNew
                                           + phaseCompFrac[ip][ic] * dPhaseAmount_dC[jc];
          localAccumJacobian[ic][jc + 1] += dPhaseCompAmount_dC;
        }
      }
    }
  }

};

/******************************** FluxKernel ********************************/

/**
 * @brief Functions to assemble flux term contributions to residual and Jacobian
 */
struct FluxKernel
{

  static inline RAJA_HOST_DEVICE void
  Compute( localIndex const NC, localIndex const NP,
           localIndex const stencilSize,
           FluxApproximationBase::CellStencil::Entry const * const stencil,
           ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const & pres,
           ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const & dPres,
           ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const & gravDepth,
           ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> const & phaseMob,
           ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> const & dPhaseMob_dPres,
           ElementRegionManager::ElementViewAccessor<arrayView3d<real64>> const & dPhaseMob_dComp,
           ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> const & dPhaseVolFrac_dPres,
           ElementRegionManager::ElementViewAccessor<arrayView3d<real64>> const & dPhaseVolFrac_dComp,
           ElementRegionManager::ElementViewAccessor<arrayView3d<real64>> const & dCompFrac_dCompDens,
           ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & phaseDens ,
           ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & dPhaseDens_dPres,
           ElementRegionManager::MaterialViewAccessor<arrayView4d<real64>> const & dPhaseDens_dComp,
           ElementRegionManager::MaterialViewAccessor<arrayView4d<real64>> const & phaseCompFrac,
           ElementRegionManager::MaterialViewAccessor<arrayView4d<real64>> const & dPhaseCompFrac_dPres,
           ElementRegionManager::MaterialViewAccessor<arrayView5d<real64>> const & dPhaseCompFrac_dComp,
           ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & phaseCapPressure,
           ElementRegionManager::MaterialViewAccessor<arrayView4d<real64>> const & dPhaseCapPressure_dPhaseVolFrac,
           localIndex const fluidIndex,
           localIndex const capPressureIndex,
           integer const gravityFlag,
           integer const capPressureFlag,
           real64 const dt,
           arraySlice1d<real64> const & localFlux,
           arraySlice2d<real64> const & localFluxJacobian )
  {
    localIndex constexpr numElems   = FluxApproximationBase::CellStencil::NUM_POINT_IN_FLUX;
    localIndex constexpr maxStencil = FluxApproximationBase::CellStencil::MAX_STENCIL_SIZE;
    localIndex constexpr maxNumComp = constitutive::MultiFluidBase::MAX_NUM_COMPONENTS;

    localIndex const NDOF = NC + 1;

    // create local work arrays

    stackArray1d<real64, maxNumComp> dPhaseCompFrac_dCompDens( NC );

    stackArray1d<real64, maxStencil>              dPhaseFlux_dP( stencilSize );
    stackArray2d<real64, maxStencil * maxNumComp> dPhaseFlux_dC( stencilSize, NC );

    stackArray1d<real64, maxNumComp>                           compFlux( NC );
    stackArray2d<real64, maxStencil * maxNumComp>              dCompFlux_dP( stencilSize, NC );
    stackArray3d<real64, maxStencil * maxNumComp * maxNumComp> dCompFlux_dC( stencilSize, NC, NC );

    stackArray1d<real64, maxNumComp> dCapPressure_dC( NC );
    stackArray1d<real64, maxNumComp> dDens_dC( NC );

    stackArray1d<real64, numElems>              dDensMean_dP( numElems );
    stackArray2d<real64, numElems * maxNumComp> dDensMean_dC( numElems, NC );

    stackArray1d<real64, maxStencil>              dPresGrad_dP( stencilSize );
    stackArray2d<real64, maxStencil * maxNumComp> dPresGrad_dC( stencilSize, NC );

    stackArray1d<real64, numElems>                dGravHead_dP( numElems );
    stackArray2d<real64, numElems * maxNumComp>   dGravHead_dC( numElems, NC );

    // reset the local values to zero
    compFlux = 0.0;
    dCompFlux_dP = 0.0;
    dCompFlux_dC = 0.0;

    for (localIndex i = 0; i < numElems * NC; ++i)
    {
      localFlux[i] = 0.0;
      for (localIndex j = 0; j < stencilSize * NDOF; ++j)
      {
        localFluxJacobian[i][j] = 0.0;
      }
    }

    // loop over phases, compute and upwind phase flux and sum contributions to each component's flux
    for (localIndex ip = 0; ip < NP; ++ip)
    {
      // clear working arrays
      real64 densMean = 0.0;
      dDensMean_dP = 0.0;
      dDensMean_dC = 0.0;

      real64 presGrad = 0.0;
      dPresGrad_dP = 0.0;
      dPresGrad_dC = 0.0;

      real64 gravHead = 0.0;
      dGravHead_dP = 0.0;
      dGravHead_dC = 0.0;

      real64 phaseFlux;
      dPhaseFlux_dP = 0.0;
      dPhaseFlux_dC = 0.0;

      // calculate quantities on primary connected cells
      for (localIndex i = 0; i < numElems; ++i)
      {
        CellDescriptor const & cell = stencil[i].index;

        localIndex const er  = cell.region;
        localIndex const esr = cell.subRegion;
        localIndex const ei  = cell.index;

        // density
        real64 const density  = phaseDens[er][esr][fluidIndex][ei][0][ip];
        real64 const dDens_dP = dPhaseDens_dPres[er][esr][fluidIndex][ei][0][ip];

        applyChainRule( NC,
                        dCompFrac_dCompDens[er][esr][ei],
                        dPhaseDens_dComp[er][esr][fluidIndex][ei][0][ip],
                        dDens_dC );

        // average density and pressure derivative
        densMean += 0.5 * density;
        dDensMean_dP[i] = 0.5 * dDens_dP;

        // compositional derivatives
        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dDensMean_dC[i][jc] = 0.5 * dDens_dC[jc];
        }
      }

      //***** calculation of flux *****

      // compute potential difference MPFA-style
      for (localIndex i = 0; i < stencilSize; ++i)
      {
        FluxApproximationBase::CellStencil::Entry const & entry = stencil[i];

        localIndex const er  = entry.index.region;
        localIndex const esr = entry.index.subRegion;
        localIndex const ei  = entry.index.index;
        real64 const weight  = entry.weight;

        //capillary pressure
        real64 capPressure     = 0.0;
        real64 dCapPressure_dP = 0.0;
        dCapPressure_dC = 0.0;

        if (capPressureFlag)
        {
          capPressure = phaseCapPressure[er][esr][capPressureIndex][ei][0][ip];

          for (localIndex jp = 0; jp < NP; ++jp)
          {
            real64 const dCapPressure_dS = dPhaseCapPressure_dPhaseVolFrac[er][esr][capPressureIndex][ei][0][ip][jp];
            dCapPressure_dP += dCapPressure_dS * dPhaseVolFrac_dPres[er][esr][ei][jp];

            for (localIndex jc = 0; jc < NC; ++jc)
            {
              dCapPressure_dC[jc] += dCapPressure_dS * dPhaseVolFrac_dComp[er][esr][ei][jp][jc];
            }
          }
        }

        presGrad += weight * (pres[er][esr][ei] + dPres[er][esr][ei] - capPressure);
        dPresGrad_dP[i] += weight * (1 - dCapPressure_dP);
        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dPresGrad_dC[i][jc] += - weight * dCapPressure_dC[jc];
        }

        if (gravityFlag)
        {
          real64 const gravD = weight * gravDepth[er][esr][ei];
          gravHead += densMean * gravD;

          // need to add contributions from both cells the mean density depends on
          for (localIndex j = 0; j < numElems; ++j)
          {
            dGravHead_dP[j] += dDensMean_dP[j] * gravD;
            for (localIndex jc = 0; jc < NC; ++jc)
            {
              dGravHead_dC[j][jc] += dDensMean_dC[j][jc] * gravD;
            }
          }
        }
      }

      // *** upwinding ***

      // use PPU currently; advanced stuff like IHU would go here
      // TODO isolate into a kernel?

      // compute phase potential gradient
      real64 potGrad = presGrad + gravHead;

      // choose upstream cell
      localIndex const k_up = (potGrad >= 0) ? 0 : 1;

      CellDescriptor const & cell_up = stencil[k_up].index;
      localIndex er_up  = cell_up.region;
      localIndex esr_up = cell_up.subRegion;
      localIndex ei_up  = cell_up.index;

      real64 const mobility = phaseMob[er_up][esr_up][ei_up][ip];

      // skip the phase flux if phase not present or immobile upstream
      if (std::fabs(mobility) < 1e-20) // TODO better constant
      {
        continue;
      }

      // pressure gradient depends on all points in the stencil
      for (localIndex ke = 0; ke < stencilSize; ++ke)
      {
        dPhaseFlux_dP[ke] += dPresGrad_dP[ke];
        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dPhaseFlux_dC[ke][jc] += dPresGrad_dC[ke][jc];
        }

      }

      // gravitational head depends only on the two cells connected (same as mean density)
      for (localIndex ke = 0; ke < numElems; ++ke)
      {
        dPhaseFlux_dP[ke] += dGravHead_dP[ke];
        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dPhaseFlux_dC[ke][jc] += dGravHead_dC[ke][jc];
        }
      }

      // compute the phase flux and derivatives using upstream cell mobility
      phaseFlux = mobility * potGrad;
      for (localIndex ke = 0; ke < stencilSize; ++ke)
      {
        dPhaseFlux_dP[ke] *= mobility;
        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dPhaseFlux_dC[ke][jc] *= mobility;
        }
      }

      real64 const dMob_dP  = dPhaseMob_dPres[er_up][esr_up][ei_up][ip];
      arraySlice1d<real64 const> dPhaseMob_dCompSub = dPhaseMob_dComp[er_up][esr_up][ei_up][ip];

      // add contribution from upstream cell mobility derivatives
      dPhaseFlux_dP[k_up] += dMob_dP * potGrad;
      for (localIndex jc = 0; jc < NC; ++jc)
      {
        dPhaseFlux_dC[k_up][jc] += dPhaseMob_dCompSub[jc] * potGrad;
      }

      // slice some constitutive arrays to avoid too much indexing in component loop
      arraySlice1d<real64 const> phaseCompFracSub = phaseCompFrac[er_up][esr_up][fluidIndex][ei_up][0][ip];
      arraySlice1d<real64 const> dPhaseCompFrac_dPresSub = dPhaseCompFrac_dPres[er_up][esr_up][fluidIndex][ei_up][0][ip];
      arraySlice2d<real64 const> dPhaseCompFrac_dCompSub = dPhaseCompFrac_dComp[er_up][esr_up][fluidIndex][ei_up][0][ip];

      // compute component fluxes and derivatives using upstream cell composition
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        real64 const ycp = phaseCompFracSub[ic];
        compFlux[ic] += phaseFlux * ycp;

        // derivatives stemming from phase flux
        for (localIndex ke = 0; ke < stencilSize; ++ke)
        {
          dCompFlux_dP[ke][ic] += dPhaseFlux_dP[ke] * ycp;
          for (localIndex jc = 0; jc < NC; ++jc)
          {
            dCompFlux_dC[ke][ic][jc] += dPhaseFlux_dC[ke][jc] * ycp;
          }
        }

        // additional derivatives stemming from upstream cell phase composition
        dCompFlux_dP[k_up][ic] += phaseFlux * dPhaseCompFrac_dPresSub[ic];

        // convert derivatives of component fraction w.r.t. component fractions to derivatives w.r.t. component densities
        applyChainRule( NC, dCompFrac_dCompDens[er_up][esr_up][ei_up], dPhaseCompFrac_dCompSub[ic], dPhaseCompFrac_dCompDens );
        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dCompFlux_dC[k_up][ic][jc] += phaseFlux * dPhaseCompFrac_dCompDens[jc];
        }
      }
    }

    // *** end of upwinding

    // populate local flux vector and derivatives
    for (localIndex ic = 0; ic < NC; ++ic)
    {
      localFlux[ic]      =  dt * compFlux[ic];
      localFlux[NC + ic] = -dt * compFlux[ic];

      for (localIndex ke = 0; ke < stencilSize; ++ke)
      {
        localIndex const localDofIndexPres = ke * NDOF;
        localFluxJacobian[ic][localDofIndexPres] = dt * dCompFlux_dP[ke][ic];
        localFluxJacobian[NC + ic][localDofIndexPres] = -dt * dCompFlux_dP[ke][ic];

        for (localIndex jc = 0; jc < NC; ++jc)
        {
          localIndex const localDofIndexComp = localDofIndexPres + jc + 1;
          localFluxJacobian[ic][localDofIndexComp] = dt * dCompFlux_dC[ke][ic][jc];
          localFluxJacobian[NC + ic][localDofIndexComp] = -dt * dCompFlux_dC[ke][ic][jc];
        }
      }
    }
  }

};

/******************************** VolumeBalanceKernel ********************************/

/**
 * @brief Functions to assemble volume balance contributions to residual and Jacobian
 */
struct VolumeBalanceKernel
{

  static inline RAJA_HOST_DEVICE void
  Compute( localIndex const NC, localIndex const NP,
           real64 const & volume,
           real64 const & porosityRef,
           real64 const & pvMult,
           real64 const & dPvMult_dPres,
           arraySlice1d<real64 const> const phaseVolFrac,
           arraySlice1d<real64 const> const dPhaseVolFrac_dPres,
           arraySlice2d<real64 const> const dPhaseVolFrac_dCompDens,
           real64 & localVolBalance,
           arraySlice1d<real64> const & localVolBalanceJacobian )
  {
    localIndex const NDOF = NC + 1;

    real64 const poro     = porosityRef * pvMult;
    real64 const dPoro_dP = porosityRef * dPvMult_dPres;

    real64 const poreVol     = volume * poro;
    real64 const dPoreVol_dP = volume * dPoro_dP;

    localVolBalance = 1.0;
    for (localIndex i = 0; i < NDOF; ++i)
    {
      localVolBalanceJacobian[i] = 0.0;
    }

    // sum contributions to component accumulation from each phase
    for (localIndex ip = 0; ip < NP; ++ip)
    {
      localVolBalance -= phaseVolFrac[ip];
      localVolBalanceJacobian[0] -= dPhaseVolFrac_dPres[ip];

      for (localIndex jc = 0; jc < NC; ++jc)
      {
        localVolBalanceJacobian[jc+1] -= dPhaseVolFrac_dCompDens[ip][jc];
      }
    }

    // scale saturation-based volume balance by pore volume (for better scaling w.r.t. other equations)
    for (localIndex idof = 0; idof < NDOF; ++idof)
    {
      localVolBalanceJacobian[idof] *= poreVol;
    }
    localVolBalanceJacobian[0] += dPoreVol_dP * localVolBalance;
    localVolBalance *= poreVol;
  }

  template<localIndex NC, localIndex NP>
  static inline RAJA_HOST_DEVICE void
  Compute( real64 const & volume,
           real64 const & porosityRef,
           real64 const & pvMult,
           real64 const & dPvMult_dPres,
           arraySlice1d<real64 const> const phaseVolFrac,
           arraySlice1d<real64 const> const dPhaseVolFrac_dPres,
           arraySlice2d<real64 const> const dPhaseVolFrac_dCompDens,
           real64 & localVolBalance,
           arraySlice1d<real64> const & localVolBalanceJacobian )
  {
    localIndex constexpr NDOF = NC + 1;

    real64 const poro     = porosityRef * pvMult;
    real64 const dPoro_dP = porosityRef * dPvMult_dPres;

    real64 const poreVol     = volume * poro;
    real64 const dPoreVol_dP = volume * dPoro_dP;

    localVolBalance = 1.0;
    for (localIndex i = 0; i < NDOF; ++i)
    {
      localVolBalanceJacobian[i] = 0.0;
    }

    // sum contributions to component accumulation from each phase
    for (localIndex ip = 0; ip < NP; ++ip)
    {
      localVolBalance -= phaseVolFrac[ip];
      localVolBalanceJacobian[0] -= dPhaseVolFrac_dPres[ip];

      for (localIndex jc = 0; jc < NC; ++jc)
      {
        localVolBalanceJacobian[jc+1] -= dPhaseVolFrac_dCompDens[ip][jc];
      }
    }

    // scale saturation-based volume balance by pore volume (for better scaling w.r.t. other equations)
    for (localIndex idof = 0; idof < NDOF; ++idof)
    {
      localVolBalanceJacobian[idof] *= poreVol;
    }
    localVolBalanceJacobian[0] += dPoreVol_dP * localVolBalance;
    localVolBalance *= poreVol;
  }

};

/******************************** Kernel launch machinery ********************************/

namespace helpers
{

template<typename T, typename LAMBDA>
bool KernelLaunchSelectorCompSwitch( T value, LAMBDA && lambda )
{
  static_assert( std::is_integral<T>::value, "KernelLaunchSelectorCompSwitch: type should be integral" );

  switch (value)
  {
    case 1:  lambda( std::integral_constant<T, 1>() );  return true;
    case 2:  lambda( std::integral_constant<T, 2>() );  return true;
    case 3:  lambda( std::integral_constant<T, 3>() );  return true;
    case 4:  lambda( std::integral_constant<T, 4>() );  return true;
    case 5:  lambda( std::integral_constant<T, 5>() );  return true;
    default: return false;
  }
}

template<typename T, typename LAMBDA>
bool KernelLaunchSelectorPhaseSwitch( T value, LAMBDA && lambda )
{
  static_assert( std::is_integral<T>::value, "KernelLaunchSelectorPhaseSwitch: type should be integral" );

  switch (value)
  {
    case 1:  lambda( std::integral_constant<T, 1>() ); return true;
    case 2:  lambda( std::integral_constant<T, 2>() ); return true;
    case 3:  lambda( std::integral_constant<T, 3>() ); return true;
    default: return false;
  }
}

} // namespace helpers

template<typename KERNELWRAPPER, typename... ARGS>
void KernelLaunchSelector1( localIndex numComp, ARGS && ... args )
{
  bool const run =
  helpers::KernelLaunchSelectorCompSwitch( numComp, [&] (auto NC)
  {
    KERNELWRAPPER::template Launch<NC()>( std::forward<ARGS>(args)... );
  });

  if (!run)
  {
    KERNELWRAPPER::Launch( numComp, std::forward<ARGS>(args)... );
  }
}

template<typename KERNELWRAPPER, typename... ARGS>
void KernelLaunchSelector2( localIndex numComp, localIndex numPhase, ARGS && ... args )
{
  bool run2 = false;
  // gcc-7 produces bugged code without explicit capture list here...
  bool const run =
  helpers::KernelLaunchSelectorCompSwitch( numComp, [ &numPhase, &run2, &args... ] ( auto NC )
  {
    run2 =
    helpers::KernelLaunchSelectorPhaseSwitch( numPhase, [&] ( auto NP )
    {
      // damn you stupid C++ rules (https://stackoverflow.com/questions/43665610)
      auto constexpr NC_ = decltype(NC)::value;
      KERNELWRAPPER::template Launch<NC_, NP()>( std::forward<ARGS>(args)... );
    });
  });

  if (!run || !run2)
  {
    KERNELWRAPPER::Launch( numComp, numPhase, std::forward<ARGS>(args)... );
  }
}

} // namespace CompositionalMultiphaseFlowKernels

} // namespace geosx


#endif //GEOSX_COMPOSITIONALMULTIPHASEFLOWKERNELS_HPP
