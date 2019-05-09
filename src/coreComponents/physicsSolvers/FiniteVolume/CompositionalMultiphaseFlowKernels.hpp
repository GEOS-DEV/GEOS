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
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

namespace geosx
{

namespace CompositionalMultiphaseFlowKernels
{

/******************************** UpdateComponentFractionKernel ********************************/

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
  static void Launch( set<localIndex> targetSet,
                      arrayView2d<real64 const> const & compDens,
                      arrayView2d<real64 const> const & dCompDens,
                      arrayView2d<real64> const & compFrac,
                      arrayView3d<real64> const & dCompFrac_dCompDens );

  static void Launch( localIndex NC,
                      set<localIndex> targetSet,
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
void ComponentFractionKernel::Launch<NC>( set<localIndex> targetSet, \
                                          arrayView2d<real64 const> const & compDens, \
                                          arrayView2d<real64 const> const & dCompDens, \
                                          arrayView2d<real64> const & compFrac, \
                                          arrayView3d<real64> const & dCompFrac_dCompDens )

INST_ComponentFractionKernel(1);
INST_ComponentFractionKernel(2);
INST_ComponentFractionKernel(3);
INST_ComponentFractionKernel(4);
INST_ComponentFractionKernel(5);
INST_ComponentFractionKernel(6);
INST_ComponentFractionKernel(7);
INST_ComponentFractionKernel(8);
INST_ComponentFractionKernel(9);
INST_ComponentFractionKernel(10);

#undef INST_ComponentFractionKernel

/******************************** UpdatePhaseVolumeFractionKernel ********************************/

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
  static void Launch( set<localIndex> targetSet,
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
                      set<localIndex> targetSet,
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
void PhaseVolumeFractionKernel::Launch<NC,NP>( set<localIndex> targetSet, \
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
INST_PhaseVolumeFractionKernel(6,1);
INST_PhaseVolumeFractionKernel(7,1);
INST_PhaseVolumeFractionKernel(8,1);
INST_PhaseVolumeFractionKernel(9,1);
INST_PhaseVolumeFractionKernel(10,1);

INST_PhaseVolumeFractionKernel(1,2);
INST_PhaseVolumeFractionKernel(2,2);
INST_PhaseVolumeFractionKernel(3,2);
INST_PhaseVolumeFractionKernel(4,2);
INST_PhaseVolumeFractionKernel(5,2);
INST_PhaseVolumeFractionKernel(6,2);
INST_PhaseVolumeFractionKernel(7,2);
INST_PhaseVolumeFractionKernel(8,2);
INST_PhaseVolumeFractionKernel(9,2);
INST_PhaseVolumeFractionKernel(10,2);

INST_PhaseVolumeFractionKernel(1,3);
INST_PhaseVolumeFractionKernel(2,3);
INST_PhaseVolumeFractionKernel(3,3);
INST_PhaseVolumeFractionKernel(4,3);
INST_PhaseVolumeFractionKernel(5,3);
INST_PhaseVolumeFractionKernel(6,3);
INST_PhaseVolumeFractionKernel(7,3);
INST_PhaseVolumeFractionKernel(8,3);
INST_PhaseVolumeFractionKernel(9,3);
INST_PhaseVolumeFractionKernel(10,3);

#undef INST_PhaseVolumeFractionKernel

/******************************** UpdatePhaseMobilityKernel ********************************/

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
  static void Launch( set<localIndex> targetSet,
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
                      set<localIndex> targetSet,
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
void PhaseMobilityKernel::Launch<NC,NP>( set<localIndex> targetSet, \
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
INST_PhaseMobilityKernel(6,1);
INST_PhaseMobilityKernel(7,1);
INST_PhaseMobilityKernel(8,1);
INST_PhaseMobilityKernel(9,1);
INST_PhaseMobilityKernel(10,1);

INST_PhaseMobilityKernel(1,2);
INST_PhaseMobilityKernel(2,2);
INST_PhaseMobilityKernel(3,2);
INST_PhaseMobilityKernel(4,2);
INST_PhaseMobilityKernel(5,2);
INST_PhaseMobilityKernel(6,2);
INST_PhaseMobilityKernel(7,2);
INST_PhaseMobilityKernel(8,2);
INST_PhaseMobilityKernel(9,2);
INST_PhaseMobilityKernel(10,2);

INST_PhaseMobilityKernel(1,3);
INST_PhaseMobilityKernel(2,3);
INST_PhaseMobilityKernel(3,3);
INST_PhaseMobilityKernel(4,3);
INST_PhaseMobilityKernel(5,3);
INST_PhaseMobilityKernel(6,3);
INST_PhaseMobilityKernel(7,3);
INST_PhaseMobilityKernel(8,3);
INST_PhaseMobilityKernel(9,3);
INST_PhaseMobilityKernel(10,3);

#undef INST_PhaseMobilityKernel

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
    case 6:  lambda( std::integral_constant<T, 6>() );  return true;
    case 7:  lambda( std::integral_constant<T, 7>() );  return true;
    case 8:  lambda( std::integral_constant<T, 8>() );  return true;
    case 9:  lambda( std::integral_constant<T, 9>() );  return true;
    case 10: lambda( std::integral_constant<T, 10>() ); return true;
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
