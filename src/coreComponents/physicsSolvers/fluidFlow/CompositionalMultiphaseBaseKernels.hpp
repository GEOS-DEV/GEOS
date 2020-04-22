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
 * @file CompositionalMultiphaseBaseKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEBASEKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEBASEKERNELS_HPP

#include "common/DataTypes.hpp"
#include "constitutive/fluid/MultiFluidBase.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

namespace geosx
{

namespace CompositionalMultiphaseBaseKernels
{

/******************************** ComponentFractionKernel ********************************/

/**
 * @brief Functions to compute component fractions from global component densities (mass or molar)
 */
struct ComponentFractionKernel
{

  template< localIndex NC >
  static inline void
  Compute( arraySlice1d< real64 const > compDens,
           arraySlice1d< real64 const > dCompDens,
           arraySlice1d< real64 > compFrac,
           arraySlice2d< real64 > dCompFrac_dCompDens );

  static inline void
  Compute( localIndex NC,
           arraySlice1d< real64 const > compDens,
           arraySlice1d< real64 const > dCompDens,
           arraySlice1d< real64 > compFrac,
           arraySlice2d< real64 > dCompFrac_dCompDens );

  template< localIndex NC >
  static void Launch( localIndex const size,
                      arrayView2d< real64 const > const & compDens,
                      arrayView2d< real64 const > const & dCompDens,
                      arrayView2d< real64 > const & compFrac,
                      arrayView3d< real64 > const & dCompFrac_dCompDens );

  static void Launch( localIndex const NC,
                      localIndex const size,
                      arrayView2d< real64 const > const & compDens,
                      arrayView2d< real64 const > const & dCompDens,
                      arrayView2d< real64 > const & compFrac,
                      arrayView3d< real64 > const & dCompFrac_dCompDens );

  template< localIndex NC >
  static void Launch( SortedArrayView< localIndex const > const & targetSet,
                      arrayView2d< real64 const > const & compDens,
                      arrayView2d< real64 const > const & dCompDens,
                      arrayView2d< real64 > const & compFrac,
                      arrayView3d< real64 > const & dCompFrac_dCompDens );

  static void Launch( localIndex NC,
                      SortedArrayView< localIndex const > const & targetSet,
                      arrayView2d< real64 const > const & compDens,
                      arrayView2d< real64 const > const & dCompDens,
                      arrayView2d< real64 > const & compFrac,
                      arrayView3d< real64 > const & dCompFrac_dCompDens );

};

#define INST_ComponentFractionKernel( NC ) \
  extern template \
  void ComponentFractionKernel::Launch< NC >( localIndex const size, \
                                              arrayView2d< real64 const > const & compDens, \
                                              arrayView2d< real64 const > const & dCompDens, \
                                              arrayView2d< real64 > const & compFrac, \
                                              arrayView3d< real64 > const & dCompFrac_dCompDens ); \
  extern template \
  void ComponentFractionKernel::Launch< NC >( SortedArrayView< localIndex const > const & targetSet, \
                                              arrayView2d< real64 const > const & compDens, \
                                              arrayView2d< real64 const > const & dCompDens, \
                                              arrayView2d< real64 > const & compFrac, \
                                              arrayView3d< real64 > const & dCompFrac_dCompDens )

INST_ComponentFractionKernel( 1 );
INST_ComponentFractionKernel( 2 );
INST_ComponentFractionKernel( 3 );
INST_ComponentFractionKernel( 4 );
INST_ComponentFractionKernel( 5 );

#undef INST_ComponentFractionKernel

/******************************** PhaseVolumeFractionKernel ********************************/

/**
 * @brief Functions to compute phase volume fractions (saturations) and derivatives
 */
struct PhaseVolumeFractionKernel
{

  template< localIndex NC, localIndex NP >
  static inline void
  Compute( arraySlice1d< real64 const > compDens,
           arraySlice1d< real64 const > dCompDens,
           arraySlice2d< real64 const > dCompFrac_dCompDens,
           arraySlice1d< real64 const > phaseDens,
           arraySlice1d< real64 const > dPhaseDens_dPres,
           arraySlice2d< real64 const > dPhaseDens_dComp,
           arraySlice1d< real64 const > phaseFrac,
           arraySlice1d< real64 const > dPhaseFrac_dPres,
           arraySlice2d< real64 const > dPhaseFrac_dComp,
           arraySlice1d< real64 > phaseVolFrac,
           arraySlice1d< real64 > dPhaseVolFrac_dPres,
           arraySlice2d< real64 > dPhaseVolFrac_dComp );

  static inline void
  Compute( localIndex NC, localIndex NP,
           arraySlice1d< real64 const > compDens,
           arraySlice1d< real64 const > dCompDens,
           arraySlice2d< real64 const > dCompFrac_dCompDens,
           arraySlice1d< real64 const > phaseDens,
           arraySlice1d< real64 const > dPhaseDens_dPres,
           arraySlice2d< real64 const > dPhaseDens_dComp,
           arraySlice1d< real64 const > phaseFrac,
           arraySlice1d< real64 const > dPhaseFrac_dPres,
           arraySlice2d< real64 const > dPhaseFrac_dComp,
           arraySlice1d< real64 > phaseVolFrac,
           arraySlice1d< real64 > dPhaseVolFrac_dPres,
           arraySlice2d< real64 > dPhaseVolFrac_dComp );

  template< localIndex NC, localIndex NP >
  static void Launch( localIndex const size,
                      arrayView2d< real64 const > const & compDens,
                      arrayView2d< real64 const > const & dCompDens,
                      arrayView3d< real64 const > const & dCompFrac_dCompDens,
                      arrayView3d< real64 const > const & phaseDens,
                      arrayView3d< real64 const > const & dPhaseDens_dPres,
                      arrayView4d< real64 const > const & dPhaseDens_dComp,
                      arrayView3d< real64 const > const & phaseFrac,
                      arrayView3d< real64 const > const & dPhaseFrac_dPres,
                      arrayView4d< real64 const > const & dPhaseFrac_dComp,
                      arrayView2d< real64 > const & phaseVolFrac,
                      arrayView2d< real64 > const & dPhaseVolFrac_dPres,
                      arrayView3d< real64 > const & dPhaseVolFrac_dComp );

  static void Launch( localIndex const NC, localIndex const NP,
                      localIndex const size,
                      arrayView2d< real64 const > const & compDens,
                      arrayView2d< real64 const > const & dCompDens,
                      arrayView3d< real64 const > const & dCompFrac_dCompDens,
                      arrayView3d< real64 const > const & phaseDens,
                      arrayView3d< real64 const > const & dPhaseDens_dPres,
                      arrayView4d< real64 const > const & dPhaseDens_dComp,
                      arrayView3d< real64 const > const & phaseFrac,
                      arrayView3d< real64 const > const & dPhaseFrac_dPres,
                      arrayView4d< real64 const > const & dPhaseFrac_dComp,
                      arrayView2d< real64 > const & phaseVolFrac,
                      arrayView2d< real64 > const & dPhaseVolFrac_dPres,
                      arrayView3d< real64 > const & dPhaseVolFrac_dComp );

  template< localIndex NC, localIndex NP >
  static void Launch( SortedArrayView< localIndex const > const & targetSet,
                      arrayView2d< real64 const > const & compDens,
                      arrayView2d< real64 const > const & dCompDens,
                      arrayView3d< real64 const > const & dCompFrac_dCompDens,
                      arrayView3d< real64 const > const & phaseDens,
                      arrayView3d< real64 const > const & dPhaseDens_dPres,
                      arrayView4d< real64 const > const & dPhaseDens_dComp,
                      arrayView3d< real64 const > const & phaseFrac,
                      arrayView3d< real64 const > const & dPhaseFrac_dPres,
                      arrayView4d< real64 const > const & dPhaseFrac_dComp,
                      arrayView2d< real64 > const & phaseVolFrac,
                      arrayView2d< real64 > const & dPhaseVolFrac_dPres,
                      arrayView3d< real64 > const & dPhaseVolFrac_dComp );

  static void Launch( localIndex NC, localIndex NP,
                      SortedArrayView< localIndex const > const & targetSet,
                      arrayView2d< real64 const > const & compDens,
                      arrayView2d< real64 const > const & dCompDens,
                      arrayView3d< real64 const > const & dCompFrac_dCompDens,
                      arrayView3d< real64 const > const & phaseDens,
                      arrayView3d< real64 const > const & dPhaseDens_dPres,
                      arrayView4d< real64 const > const & dPhaseDens_dComp,
                      arrayView3d< real64 const > const & phaseFrac,
                      arrayView3d< real64 const > const & dPhaseFrac_dPres,
                      arrayView4d< real64 const > const & dPhaseFrac_dComp,
                      arrayView2d< real64 > const & phaseVolFrac,
                      arrayView2d< real64 > const & dPhaseVolFrac_dPres,
                      arrayView3d< real64 > const & dPhaseVolFrac_dComp );

};

#define INST_PhaseVolumeFractionKernel( NC, NP ) \
  extern template \
  void PhaseVolumeFractionKernel::Launch< NC, NP >( localIndex const size, \
                                                    arrayView2d< real64 const > const & compDens, \
                                                    arrayView2d< real64 const > const & dCompDens, \
                                                    arrayView3d< real64 const > const & dCompFrac_dCompDens, \
                                                    arrayView3d< real64 const > const & phaseDens, \
                                                    arrayView3d< real64 const > const & dPhaseDens_dPres, \
                                                    arrayView4d< real64 const > const & dPhaseDens_dComp, \
                                                    arrayView3d< real64 const > const & phaseFrac, \
                                                    arrayView3d< real64 const > const & dPhaseFrac_dPres, \
                                                    arrayView4d< real64 const > const & dPhaseFrac_dComp, \
                                                    arrayView2d< real64 > const & phaseVolFrac, \
                                                    arrayView2d< real64 > const & dPhaseVolFrac_dPres, \
                                                    arrayView3d< real64 > const & dPhaseVolFrac_dComp ); \
  extern template \
  void PhaseVolumeFractionKernel::Launch< NC, NP >( SortedArrayView< localIndex const > const & targetSet, \
                                                    arrayView2d< real64 const > const & compDens, \
                                                    arrayView2d< real64 const > const & dCompDens, \
                                                    arrayView3d< real64 const > const & dCompFrac_dCompDens, \
                                                    arrayView3d< real64 const > const & phaseDens, \
                                                    arrayView3d< real64 const > const & dPhaseDens_dPres, \
                                                    arrayView4d< real64 const > const & dPhaseDens_dComp, \
                                                    arrayView3d< real64 const > const & phaseFrac, \
                                                    arrayView3d< real64 const > const & dPhaseFrac_dPres, \
                                                    arrayView4d< real64 const > const & dPhaseFrac_dComp, \
                                                    arrayView2d< real64 > const & phaseVolFrac, \
                                                    arrayView2d< real64 > const & dPhaseVolFrac_dPres, \
                                                    arrayView3d< real64 > const & dPhaseVolFrac_dComp )

INST_PhaseVolumeFractionKernel( 1, 1 );
INST_PhaseVolumeFractionKernel( 2, 1 );
INST_PhaseVolumeFractionKernel( 3, 1 );
INST_PhaseVolumeFractionKernel( 4, 1 );
INST_PhaseVolumeFractionKernel( 5, 1 );

INST_PhaseVolumeFractionKernel( 1, 2 );
INST_PhaseVolumeFractionKernel( 2, 2 );
INST_PhaseVolumeFractionKernel( 3, 2 );
INST_PhaseVolumeFractionKernel( 4, 2 );
INST_PhaseVolumeFractionKernel( 5, 2 );

INST_PhaseVolumeFractionKernel( 1, 3 );
INST_PhaseVolumeFractionKernel( 2, 3 );
INST_PhaseVolumeFractionKernel( 3, 3 );
INST_PhaseVolumeFractionKernel( 4, 3 );
INST_PhaseVolumeFractionKernel( 5, 3 );

#undef INST_PhaseVolumeFractionKernel

/******************************** PhaseMobilityKernel ********************************/

/**
 * @brief Functions to compute phase mobilities and derivatives from density, viscosity and relperm
 */
struct PhaseMobilityKernel
{

  template< localIndex NC, localIndex NP >
  static inline void
  Compute( arraySlice2d< real64 const > dCompFrac_dCompDens,
           arraySlice1d< real64 const > phaseDens,
           arraySlice1d< real64 const > dPhaseDens_dPres,
           arraySlice2d< real64 const > dPhaseDens_dComp,
           arraySlice1d< real64 const > phaseVisc,
           arraySlice1d< real64 const > dPhaseVisc_dPres,
           arraySlice2d< real64 const > dPhaseVisc_dComp,
           arraySlice1d< real64 const > phaseRelPerm,
           arraySlice2d< real64 const > dPhaseRelPerm_dPhaseVolFrac,
           arraySlice1d< real64 const > dPhaseVolFrac_dPres,
           arraySlice2d< real64 const > dPhaseVolFrac_dComp,
           arraySlice1d< real64 > phaseMob,
           arraySlice1d< real64 > dPhaseMob_dPres,
           arraySlice2d< real64 > dPhaseMob_dComp );

  static inline void
  Compute( localIndex NC, localIndex NP,
           arraySlice2d< real64 const > dCompFrac_dCompDens,
           arraySlice1d< real64 const > phaseDens,
           arraySlice1d< real64 const > dPhaseDens_dPres,
           arraySlice2d< real64 const > dPhaseDens_dComp,
           arraySlice1d< real64 const > phaseVisc,
           arraySlice1d< real64 const > dPhaseVisc_dPres,
           arraySlice2d< real64 const > dPhaseVisc_dComp,
           arraySlice1d< real64 const > phaseRelPerm,
           arraySlice2d< real64 const > dPhaseRelPerm_dPhaseVolFrac,
           arraySlice1d< real64 const > dPhaseVolFrac_dPres,
           arraySlice2d< real64 const > dPhaseVolFrac_dComp,
           arraySlice1d< real64 > phaseMob,
           arraySlice1d< real64 > dPhaseMob_dPres,
           arraySlice2d< real64 > dPhaseMob_dComp );

  template< localIndex NC, localIndex NP >
  static void Launch( localIndex const size,
                      arrayView3d< real64 const > const & dCompFrac_dCompDens,
                      arrayView3d< real64 const > const & phaseDens,
                      arrayView3d< real64 const > const & dPhaseDens_dPres,
                      arrayView4d< real64 const > const & dPhaseDens_dComp,
                      arrayView3d< real64 const > const & phaseVisc,
                      arrayView3d< real64 const > const & dPhaseVisc_dPres,
                      arrayView4d< real64 const > const & dPhaseVisc_dComp,
                      arrayView3d< real64 const > const & phaseRelPerm,
                      arrayView4d< real64 const > const & dPhaseRelPerm_dPhaseVolFrac,
                      arrayView2d< real64 const > const & dPhaseVolFrac_dPres,
                      arrayView3d< real64 const > const & dPhaseVolFrac_dComp,
                      arrayView2d< real64 > const & phaseMob,
                      arrayView2d< real64 > const & dPhaseMob_dPres,
                      arrayView3d< real64 > const & dPhaseMob_dComp );

  static void Launch( localIndex const NC, localIndex const NP,
                      localIndex const size,
                      arrayView3d< real64 const > const & dCompFrac_dCompDens,
                      arrayView3d< real64 const > const & phaseDens,
                      arrayView3d< real64 const > const & dPhaseDens_dPres,
                      arrayView4d< real64 const > const & dPhaseDens_dComp,
                      arrayView3d< real64 const > const & phaseVisc,
                      arrayView3d< real64 const > const & dPhaseVisc_dPres,
                      arrayView4d< real64 const > const & dPhaseVisc_dComp,
                      arrayView3d< real64 const > const & phaseRelPerm,
                      arrayView4d< real64 const > const & dPhaseRelPerm_dPhaseVolFrac,
                      arrayView2d< real64 const > const & dPhaseVolFrac_dPres,
                      arrayView3d< real64 const > const & dPhaseVolFrac_dComp,
                      arrayView2d< real64 > const & phaseMob,
                      arrayView2d< real64 > const & dPhaseMob_dPres,
                      arrayView3d< real64 > const & dPhaseMob_dComp );

  template< localIndex NC, localIndex NP >
  static void Launch( SortedArrayView< localIndex const > const & targetSet,
                      arrayView3d< real64 const > const & dCompFrac_dCompDens,
                      arrayView3d< real64 const > const & phaseDens,
                      arrayView3d< real64 const > const & dPhaseDens_dPres,
                      arrayView4d< real64 const > const & dPhaseDens_dComp,
                      arrayView3d< real64 const > const & phaseVisc,
                      arrayView3d< real64 const > const & dPhaseVisc_dPres,
                      arrayView4d< real64 const > const & dPhaseVisc_dComp,
                      arrayView3d< real64 const > const & phaseRelPerm,
                      arrayView4d< real64 const > const & dPhaseRelPerm_dPhaseVolFrac,
                      arrayView2d< real64 const > const & dPhaseVolFrac_dPres,
                      arrayView3d< real64 const > const & dPhaseVolFrac_dComp,
                      arrayView2d< real64 > const & phaseMob,
                      arrayView2d< real64 > const & dPhaseMob_dPres,
                      arrayView3d< real64 > const & dPhaseMob_dComp );

  static void Launch( localIndex NC, localIndex NP,
                      SortedArrayView< localIndex const > const & targetSet,
                      arrayView3d< real64 const > const & dCompFrac_dCompDens,
                      arrayView3d< real64 const > const & phaseDens,
                      arrayView3d< real64 const > const & dPhaseDens_dPres,
                      arrayView4d< real64 const > const & dPhaseDens_dComp,
                      arrayView3d< real64 const > const & phaseVisc,
                      arrayView3d< real64 const > const & dPhaseVisc_dPres,
                      arrayView4d< real64 const > const & dPhaseVisc_dComp,
                      arrayView3d< real64 const > const & phaseRelPerm,
                      arrayView4d< real64 const > const & dPhaseRelPerm_dPhaseVolFrac,
                      arrayView2d< real64 const > const & dPhaseVolFrac_dPres,
                      arrayView3d< real64 const > const & dPhaseVolFrac_dComp,
                      arrayView2d< real64 > const & phaseMob,
                      arrayView2d< real64 > const & dPhaseMob_dPres,
                      arrayView3d< real64 > const & dPhaseMob_dComp );

};

#define INST_PhaseMobilityKernel( NC, NP ) \
  extern template \
  void PhaseMobilityKernel::Launch< NC, NP >( localIndex const size, \
                                              arrayView3d< real64 const > const & dCompFrac_dCompDens, \
                                              arrayView3d< real64 const > const & phaseDens, \
                                              arrayView3d< real64 const > const & dPhaseDens_dPres, \
                                              arrayView4d< real64 const > const & dPhaseDens_dComp, \
                                              arrayView3d< real64 const > const & phaseVisc, \
                                              arrayView3d< real64 const > const & dPhaseVisc_dPres, \
                                              arrayView4d< real64 const > const & dPhaseVisc_dComp, \
                                              arrayView3d< real64 const > const & phaseRelPerm, \
                                              arrayView4d< real64 const > const & dPhaseRelPerm_dPhaseVolFrac, \
                                              arrayView2d< real64 const > const & dPhaseVolFrac_dPres, \
                                              arrayView3d< real64 const > const & dPhaseVolFrac_dComp, \
                                              arrayView2d< real64 > const & phaseMob, \
                                              arrayView2d< real64 > const & dPhaseMob_dPres, \
                                              arrayView3d< real64 > const & dPhaseMob_dComp ); \
  extern template \
  void PhaseMobilityKernel::Launch< NC, NP >( SortedArrayView< localIndex const > const & targetSet, \
                                              arrayView3d< real64 const > const & dCompFrac_dCompDens, \
                                              arrayView3d< real64 const > const & phaseDens, \
                                              arrayView3d< real64 const > const & dPhaseDens_dPres, \
                                              arrayView4d< real64 const > const & dPhaseDens_dComp, \
                                              arrayView3d< real64 const > const & phaseVisc, \
                                              arrayView3d< real64 const > const & dPhaseVisc_dPres, \
                                              arrayView4d< real64 const > const & dPhaseVisc_dComp, \
                                              arrayView3d< real64 const > const & phaseRelPerm, \
                                              arrayView4d< real64 const > const & dPhaseRelPerm_dPhaseVolFrac, \
                                              arrayView2d< real64 const > const & dPhaseVolFrac_dPres, \
                                              arrayView3d< real64 const > const & dPhaseVolFrac_dComp, \
                                              arrayView2d< real64 > const & phaseMob, \
                                              arrayView2d< real64 > const & dPhaseMob_dPres, \
                                              arrayView3d< real64 > const & dPhaseMob_dComp )

INST_PhaseMobilityKernel( 1, 1 );
INST_PhaseMobilityKernel( 2, 1 );
INST_PhaseMobilityKernel( 3, 1 );
INST_PhaseMobilityKernel( 4, 1 );
INST_PhaseMobilityKernel( 5, 1 );

INST_PhaseMobilityKernel( 1, 2 );
INST_PhaseMobilityKernel( 2, 2 );
INST_PhaseMobilityKernel( 3, 2 );
INST_PhaseMobilityKernel( 4, 2 );
INST_PhaseMobilityKernel( 5, 2 );

INST_PhaseMobilityKernel( 1, 3 );
INST_PhaseMobilityKernel( 2, 3 );
INST_PhaseMobilityKernel( 3, 3 );
INST_PhaseMobilityKernel( 4, 3 );
INST_PhaseMobilityKernel( 5, 3 );

#undef INST_PhaseMobilityKernel

/******************************** AccumulationKernel ********************************/

/**
 * @brief Functions to assemble accumulation term contributions to residual and Jacobian
 */
struct AccumulationKernel
{

  static inline void
  Compute( localIndex const NC, localIndex const NP,
           real64 const & volume,
           real64 const & porosityOld,
           real64 const & porosityRef,
           real64 const & pvMult,
           real64 const & dPvMult_dPres,
           arraySlice2d< real64 const > const dCompFrac_dCompDens,
           arraySlice1d< real64 const > const phaseVolFracOld,
           arraySlice1d< real64 const > const phaseVolFrac,
           arraySlice1d< real64 const > const dPhaseVolFrac_dPres,
           arraySlice2d< real64 const > const dPhaseVolFrac_dCompDens,
           arraySlice1d< real64 const > const phaseDensOld,
           arraySlice1d< real64 const > const phaseDens,
           arraySlice1d< real64 const > const dPhaseDens_dPres,
           arraySlice2d< real64 const > const dPhaseDens_dComp,
           arraySlice2d< real64 const > const phaseCompFracOld,
           arraySlice2d< real64 const > const phaseCompFrac,
           arraySlice2d< real64 const > const dPhaseCompFrac_dPres,
           arraySlice3d< real64 const > const dPhaseCompFrac_dComp,
           arraySlice1d< real64 > const & localAccum,
           arraySlice2d< real64 > const & localAccumJacobian )
  {
    localIndex constexpr maxNumComp = constitutive::MultiFluidBase::MAX_NUM_COMPONENTS;

    // temporary work arrays
    stackArray1d< real64, maxNumComp > const dPhaseAmount_dC( NC );
    stackArray1d< real64, maxNumComp > const dPhaseCompFrac_dC( NC );

    // reset the local values
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      localAccum[ic] = 0.0;
      for( localIndex jc = 0; jc < NC; ++jc )
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
    for( localIndex ip = 0; ip < NP; ++ip )
    {
      real64 const phaseAmountNew = poreVolNew * phaseVolFrac[ip] * phaseDens[ip];
      real64 const phaseAmountOld = poreVolOld * phaseVolFracOld[ip] * phaseDensOld[ip];

      real64 const dPhaseAmount_dP = dPoreVol_dP * phaseVolFrac[ip] * phaseDens[ip]
                                     + poreVolNew * (dPhaseVolFrac_dPres[ip] * phaseDens[ip]
                                                     + phaseVolFrac[ip] * dPhaseDens_dPres[ip]);

      // assemble density dependence
      applyChainRule( NC, dCompFrac_dCompDens, dPhaseDens_dComp[ip], dPhaseAmount_dC );
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dPhaseAmount_dC[jc] = dPhaseAmount_dC[jc] * phaseVolFrac[ip]
                              + phaseDens[ip] * dPhaseVolFrac_dCompDens[ip][jc];
        dPhaseAmount_dC[jc] *= poreVolNew;
      }

      // ic - index of component whose conservation equation is assembled
      // (i.e. row number in local matrix)
      for( localIndex ic = 0; ic < NC; ++ic )
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
        applyChainRule( NC, dCompFrac_dCompDens, dPhaseCompFrac_dComp[ip][ic], dPhaseCompFrac_dC );
        for( localIndex jc = 0; jc < NC; ++jc )
        {
          real64 const dPhaseCompAmount_dC = dPhaseCompFrac_dC[jc] * phaseAmountNew
                                             + phaseCompFrac[ip][ic] * dPhaseAmount_dC[jc];
          localAccumJacobian[ic][jc + 1] += dPhaseCompAmount_dC;
        }
      }
    }
  }

  template< localIndex NC, localIndex NP >
  static inline void
  Compute( real64 const & volume,
           real64 const & porosityOld,
           real64 const & porosityRef,
           real64 const & pvMult,
           real64 const & dPvMult_dPres,
           arraySlice2d< real64 const > const dCompFrac_dCompDens,
           arraySlice1d< real64 const > const phaseVolFracOld,
           arraySlice1d< real64 const > const phaseVolFrac,
           arraySlice1d< real64 const > const dPhaseVolFrac_dPres,
           arraySlice2d< real64 const > const dPhaseVolFrac_dCompDens,
           arraySlice1d< real64 const > const phaseDensOld,
           arraySlice1d< real64 const > const phaseDens,
           arraySlice1d< real64 const > const dPhaseDens_dPres,
           arraySlice2d< real64 const > const dPhaseDens_dComp,
           arraySlice2d< real64 const > const phaseCompFracOld,
           arraySlice2d< real64 const > const phaseCompFrac,
           arraySlice2d< real64 const > const dPhaseCompFrac_dPres,
           arraySlice3d< real64 const > const dPhaseCompFrac_dComp,
           arraySlice1d< real64 > const & localAccum,
           arraySlice2d< real64 > const & localAccumJacobian )
  {
    // temporary work arrays
    stackArray1d< real64, NC > dPhaseAmount_dC( NC );
    stackArray1d< real64, NC > dPhaseCompFrac_dC( NC );

    // reset the local values
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      localAccum[ic] = 0.0;
      for( localIndex jc = 0; jc < NC; ++jc )
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
    for( localIndex ip = 0; ip < NP; ++ip )
    {
      real64 const phaseAmountNew = poreVolNew * phaseVolFrac[ip] * phaseDens[ip];
      real64 const phaseAmountOld = poreVolOld * phaseVolFracOld[ip] * phaseDensOld[ip];

      real64 const dPhaseAmount_dP = dPoreVol_dP * phaseVolFrac[ip] * phaseDens[ip]
                                     + poreVolNew * (dPhaseVolFrac_dPres[ip] * phaseDens[ip]
                                                     + phaseVolFrac[ip] * dPhaseDens_dPres[ip]);

      // assemble density dependence
      applyChainRule( NC, dCompFrac_dCompDens, dPhaseDens_dComp[ip], dPhaseAmount_dC );
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dPhaseAmount_dC[jc] = dPhaseAmount_dC[jc] * phaseVolFrac[ip]
                              + phaseDens[ip] * dPhaseVolFrac_dCompDens[ip][jc];
        dPhaseAmount_dC[jc] *= poreVolNew;
      }

      // ic - index of component whose conservation equation is assembled
      // (i.e. row number in local matrix)
      for( localIndex ic = 0; ic < NC; ++ic )
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
        applyChainRule( NC, dCompFrac_dCompDens, dPhaseCompFrac_dComp[ip][ic], dPhaseCompFrac_dC );
        for( localIndex jc = 0; jc < NC; ++jc )
        {
          real64 const dPhaseCompAmount_dC = dPhaseCompFrac_dC[jc] * phaseAmountNew
                                             + phaseCompFrac[ip][ic] * dPhaseAmount_dC[jc];
          localAccumJacobian[ic][jc + 1] += dPhaseCompAmount_dC;
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

  static inline void
  Compute( localIndex const NC, localIndex const NP,
           real64 const & volume,
           real64 const & porosityRef,
           real64 const & pvMult,
           real64 const & dPvMult_dPres,
           arraySlice1d< real64 const > const phaseVolFrac,
           arraySlice1d< real64 const > const dPhaseVolFrac_dPres,
           arraySlice2d< real64 const > const dPhaseVolFrac_dCompDens,
           real64 & localVolBalance,
           arraySlice1d< real64 > const & localVolBalanceJacobian )
  {
    localIndex const NDOF = NC + 1;

    real64 const poro     = porosityRef * pvMult;
    real64 const dPoro_dP = porosityRef * dPvMult_dPres;

    real64 const poreVol     = volume * poro;
    real64 const dPoreVol_dP = volume * dPoro_dP;

    localVolBalance = 1.0;
    for( localIndex i = 0; i < NDOF; ++i )
    {
      localVolBalanceJacobian[i] = 0.0;
    }

    // sum contributions to component accumulation from each phase
    for( localIndex ip = 0; ip < NP; ++ip )
    {
      localVolBalance -= phaseVolFrac[ip];
      localVolBalanceJacobian[0] -= dPhaseVolFrac_dPres[ip];

      for( localIndex jc = 0; jc < NC; ++jc )
      {
        localVolBalanceJacobian[jc+1] -= dPhaseVolFrac_dCompDens[ip][jc];
      }
    }

    // scale saturation-based volume balance by pore volume (for better scaling w.r.t. other equations)
    for( localIndex idof = 0; idof < NDOF; ++idof )
    {
      localVolBalanceJacobian[idof] *= poreVol;
    }
    localVolBalanceJacobian[0] += dPoreVol_dP * localVolBalance;
    localVolBalance *= poreVol;
  }

  template< localIndex NC, localIndex NP >
  static inline void
  Compute( real64 const & volume,
           real64 const & porosityRef,
           real64 const & pvMult,
           real64 const & dPvMult_dPres,
           arraySlice1d< real64 const > const phaseVolFrac,
           arraySlice1d< real64 const > const dPhaseVolFrac_dPres,
           arraySlice2d< real64 const > const dPhaseVolFrac_dCompDens,
           real64 & localVolBalance,
           arraySlice1d< real64 > const & localVolBalanceJacobian )
  {
    localIndex constexpr NDOF = NC + 1;

    real64 const poro     = porosityRef * pvMult;
    real64 const dPoro_dP = porosityRef * dPvMult_dPres;

    real64 const poreVol     = volume * poro;
    real64 const dPoreVol_dP = volume * dPoro_dP;

    localVolBalance = 1.0;
    for( localIndex i = 0; i < NDOF; ++i )
    {
      localVolBalanceJacobian[i] = 0.0;
    }

    // sum contributions to component accumulation from each phase
    for( localIndex ip = 0; ip < NP; ++ip )
    {
      localVolBalance -= phaseVolFrac[ip];
      localVolBalanceJacobian[0] -= dPhaseVolFrac_dPres[ip];

      for( localIndex jc = 0; jc < NC; ++jc )
      {
        localVolBalanceJacobian[jc+1] -= dPhaseVolFrac_dCompDens[ip][jc];
      }
    }

    // scale saturation-based volume balance by pore volume (for better scaling w.r.t. other equations)
    for( localIndex idof = 0; idof < NDOF; ++idof )
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

template< typename T, typename LAMBDA >
void KernelLaunchSelectorCompSwitch( T value, LAMBDA && lambda )
{
  static_assert( std::is_integral< T >::value, "KernelLaunchSelectorCompSwitch: type should be integral" );

  switch( value )
  {
    case 1:
    { lambda( std::integral_constant< T, 1 >() ); return;}
    case 2:
    { lambda( std::integral_constant< T, 2 >() ); return;}
    case 3:
    { lambda( std::integral_constant< T, 3 >() ); return;}
    case 4:
    { lambda( std::integral_constant< T, 4 >() ); return;}
    case 5:
    { lambda( std::integral_constant< T, 5 >() ); return;}
    default: GEOSX_ERROR( "Unknown numComp value: " << value );
  }
}

} // namespace helpers

template< typename KERNELWRAPPER, typename ... ARGS >
void KernelLaunchSelector1( localIndex numComp, ARGS && ... args )
{
  helpers::KernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
  {
    KERNELWRAPPER::template Launch< NC() >( std::forward< ARGS >( args )... );
  } );
}

template< typename KERNELWRAPPER, typename ... ARGS >
void KernelLaunchSelector2( localIndex numComp, localIndex numPhase, ARGS && ... args )
{
  helpers::KernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
  {
    switch( numPhase )
    {
      case 1:
        { KERNELWRAPPER::template Launch< NC(), 1 >( std::forward< ARGS >( args )... ); return;}
      case 2:
        { KERNELWRAPPER::template Launch< NC(), 2 >( std::forward< ARGS >( args )... ); return;}
      case 3:
        { KERNELWRAPPER::template Launch< NC(), 3 >( std::forward< ARGS >( args )... ); return;}
      default: GEOSX_ERROR( "Unknown numPhase value: " << numPhase );
    }
  } );
}

} // namespace CompositionalMultiphaseBaseKernels

} // namespace geosx


#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEBASEKERNELS_HPP
