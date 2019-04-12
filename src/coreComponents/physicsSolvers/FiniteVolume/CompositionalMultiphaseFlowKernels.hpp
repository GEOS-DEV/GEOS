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

struct UpdateComponentFractionKernel
{

  template<localIndex NC>
  static inline RAJA_HOST_DEVICE void
  Compute( arraySlice1d<real64 const> const & compDens,
           arraySlice1d<real64 const> const & dCompDens,
           arraySlice1d<real64> const & compFrac,
           arraySlice2d<real64> const & dCompFrac_dCompDens );

  static inline RAJA_HOST_DEVICE void
  Compute( localIndex NC,
           arraySlice1d<real64 const> const & compDens,
           arraySlice1d<real64 const> const & dCompDens,
           arraySlice1d<real64> const & compFrac,
           arraySlice2d<real64> const & dCompFrac_dCompDens );

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

};

#define LAUNCH_UpdateComponentFractionKernel(NC) \
extern template \
void UpdateComponentFractionKernel::Launch<NC>( localIndex begin, localIndex end, \
                                                arrayView2d<real64 const> const & compDens, \
                                                arrayView2d<real64 const> const & dCompDens, \
                                                arrayView2d<real64> const & compFrac, \
                                                arrayView3d<real64> const & dCompFrac_dCompDens )

LAUNCH_UpdateComponentFractionKernel(1);
LAUNCH_UpdateComponentFractionKernel(2);
LAUNCH_UpdateComponentFractionKernel(3);
LAUNCH_UpdateComponentFractionKernel(4);
LAUNCH_UpdateComponentFractionKernel(5);
LAUNCH_UpdateComponentFractionKernel(6);
LAUNCH_UpdateComponentFractionKernel(7);
LAUNCH_UpdateComponentFractionKernel(8);
LAUNCH_UpdateComponentFractionKernel(9);
LAUNCH_UpdateComponentFractionKernel(10);

#undef LAUNCH_UpdateComponentFractionKernel

/******************************** UpdatePhaseVolumeFractionKernel ********************************/

struct UpdatePhaseVolumeFractionKernel
{

  template<localIndex NC, localIndex NP>
  static inline RAJA_HOST_DEVICE void
  Compute( arraySlice1d<real64 const> const & compDens,
           arraySlice1d<real64 const> const & dCompDens,
           arraySlice2d<real64 const> const & dCompFrac_dCompDens,
           arraySlice1d<real64 const> const & phaseDens,
           arraySlice1d<real64 const> const & dPhaseDens_dPres,
           arraySlice2d<real64 const> const & dPhaseDens_dComp,
           arraySlice1d<real64 const> const & phaseFrac,
           arraySlice1d<real64 const> const & dPhaseFrac_dPres,
           arraySlice2d<real64 const> const & dPhaseFrac_dComp,
           arraySlice1d<real64> const & phaseVolFrac,
           arraySlice1d<real64> const & dPhaseVolFrac_dPres,
           arraySlice2d<real64> const & dPhaseVolFrac_dComp );

  static inline RAJA_HOST_DEVICE void
  Compute( localIndex NC, localIndex NP,
           arraySlice1d<real64 const> const & compDens,
           arraySlice1d<real64 const> const & dCompDens,
           arraySlice2d<real64 const> const & dCompFrac_dCompDens,
           arraySlice1d<real64 const> const & phaseDens,
           arraySlice1d<real64 const> const & dPhaseDens_dPres,
           arraySlice2d<real64 const> const & dPhaseDens_dComp,
           arraySlice1d<real64 const> const & phaseFrac,
           arraySlice1d<real64 const> const & dPhaseFrac_dPres,
           arraySlice2d<real64 const> const & dPhaseFrac_dComp,
           arraySlice1d<real64> const & phaseVolFrac,
           arraySlice1d<real64> const & dPhaseVolFrac_dPres,
           arraySlice2d<real64> const & dPhaseVolFrac_dComp );

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

};

#define LAUNCH_UpdatePhaseVolumeFractionKernel(NC,NP) \
extern template \
void UpdatePhaseVolumeFractionKernel::Launch<NC,NP>( localIndex begin, localIndex end, \
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

LAUNCH_UpdatePhaseVolumeFractionKernel(1,1);
LAUNCH_UpdatePhaseVolumeFractionKernel(2,1);
LAUNCH_UpdatePhaseVolumeFractionKernel(3,1);
LAUNCH_UpdatePhaseVolumeFractionKernel(4,1);
LAUNCH_UpdatePhaseVolumeFractionKernel(5,1);
LAUNCH_UpdatePhaseVolumeFractionKernel(6,1);
LAUNCH_UpdatePhaseVolumeFractionKernel(7,1);
LAUNCH_UpdatePhaseVolumeFractionKernel(8,1);
LAUNCH_UpdatePhaseVolumeFractionKernel(9,1);
LAUNCH_UpdatePhaseVolumeFractionKernel(10,1);

LAUNCH_UpdatePhaseVolumeFractionKernel(1,2);
LAUNCH_UpdatePhaseVolumeFractionKernel(2,2);
LAUNCH_UpdatePhaseVolumeFractionKernel(3,2);
LAUNCH_UpdatePhaseVolumeFractionKernel(4,2);
LAUNCH_UpdatePhaseVolumeFractionKernel(5,2);
LAUNCH_UpdatePhaseVolumeFractionKernel(6,2);
LAUNCH_UpdatePhaseVolumeFractionKernel(7,2);
LAUNCH_UpdatePhaseVolumeFractionKernel(8,2);
LAUNCH_UpdatePhaseVolumeFractionKernel(9,2);
LAUNCH_UpdatePhaseVolumeFractionKernel(10,2);

LAUNCH_UpdatePhaseVolumeFractionKernel(1,3);
LAUNCH_UpdatePhaseVolumeFractionKernel(2,3);
LAUNCH_UpdatePhaseVolumeFractionKernel(3,3);
LAUNCH_UpdatePhaseVolumeFractionKernel(4,3);
LAUNCH_UpdatePhaseVolumeFractionKernel(5,3);
LAUNCH_UpdatePhaseVolumeFractionKernel(6,3);
LAUNCH_UpdatePhaseVolumeFractionKernel(7,3);
LAUNCH_UpdatePhaseVolumeFractionKernel(8,3);
LAUNCH_UpdatePhaseVolumeFractionKernel(9,3);
LAUNCH_UpdatePhaseVolumeFractionKernel(10,3);

#undef LAUNCH_UpdatePhaseVolumeFractionKernel

/******************************** UpdatePhaseMobilityKernel ********************************/

struct UpdatePhaseMobilityKernel
{

  template<localIndex NC, localIndex NP>
  static inline RAJA_HOST_DEVICE void
  Compute( arraySlice2d<real64 const> const & dCompFrac_dCompDens,
           arraySlice1d<real64 const> const & phaseDens,
           arraySlice1d<real64 const> const & dPhaseDens_dPres,
           arraySlice2d<real64 const> const & dPhaseDens_dComp,
           arraySlice1d<real64 const> const & phaseVisc,
           arraySlice1d<real64 const> const & dPhaseVisc_dPres,
           arraySlice2d<real64 const> const & dPhaseVisc_dComp,
           arraySlice1d<real64 const> const & phaseRelPerm,
           arraySlice2d<real64 const> const & dPhaseRelPerm_dPhaseVolFrac,
           arraySlice1d<real64 const> const & dPhaseVolFrac_dPres,
           arraySlice2d<real64 const> const & dPhaseVolFrac_dComp,
           arraySlice1d<real64> const & phaseMob,
           arraySlice1d<real64> const & dPhaseMob_dPres,
           arraySlice2d<real64> const & dPhaseMob_dComp );

  static inline RAJA_HOST_DEVICE void
  Compute( localIndex NC, localIndex NP,
           arraySlice2d<real64 const> const & dCompFrac_dCompDens,
           arraySlice1d<real64 const> const & phaseDens,
           arraySlice1d<real64 const> const & dPhaseDens_dPres,
           arraySlice2d<real64 const> const & dPhaseDens_dComp,
           arraySlice1d<real64 const> const & phaseVisc,
           arraySlice1d<real64 const> const & dPhaseVisc_dPres,
           arraySlice2d<real64 const> const & dPhaseVisc_dComp,
           arraySlice1d<real64 const> const & phaseRelPerm,
           arraySlice2d<real64 const> const & dPhaseRelPerm_dPhaseVolFrac,
           arraySlice1d<real64 const> const & dPhaseVolFrac_dPres,
           arraySlice2d<real64 const> const & dPhaseVolFrac_dComp,
           arraySlice1d<real64> const & phaseMob,
           arraySlice1d<real64> const & dPhaseMob_dPres,
           arraySlice2d<real64> const & dPhaseMob_dComp );

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

};

#define LAUNCH_UpdatePhaseMobilityKernel(NC,NP) \
extern template \
void UpdatePhaseMobilityKernel::Launch<NC,NP>( localIndex begin, localIndex end, \
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

LAUNCH_UpdatePhaseMobilityKernel(1,1);
LAUNCH_UpdatePhaseMobilityKernel(2,1);
LAUNCH_UpdatePhaseMobilityKernel(3,1);
LAUNCH_UpdatePhaseMobilityKernel(4,1);
LAUNCH_UpdatePhaseMobilityKernel(5,1);
LAUNCH_UpdatePhaseMobilityKernel(6,1);
LAUNCH_UpdatePhaseMobilityKernel(7,1);
LAUNCH_UpdatePhaseMobilityKernel(8,1);
LAUNCH_UpdatePhaseMobilityKernel(9,1);
LAUNCH_UpdatePhaseMobilityKernel(10,1);

LAUNCH_UpdatePhaseMobilityKernel(1,2);
LAUNCH_UpdatePhaseMobilityKernel(2,2);
LAUNCH_UpdatePhaseMobilityKernel(3,2);
LAUNCH_UpdatePhaseMobilityKernel(4,2);
LAUNCH_UpdatePhaseMobilityKernel(5,2);
LAUNCH_UpdatePhaseMobilityKernel(6,2);
LAUNCH_UpdatePhaseMobilityKernel(7,2);
LAUNCH_UpdatePhaseMobilityKernel(8,2);
LAUNCH_UpdatePhaseMobilityKernel(9,2);
LAUNCH_UpdatePhaseMobilityKernel(10,2);

LAUNCH_UpdatePhaseMobilityKernel(1,3);
LAUNCH_UpdatePhaseMobilityKernel(2,3);
LAUNCH_UpdatePhaseMobilityKernel(3,3);
LAUNCH_UpdatePhaseMobilityKernel(4,3);
LAUNCH_UpdatePhaseMobilityKernel(5,3);
LAUNCH_UpdatePhaseMobilityKernel(6,3);
LAUNCH_UpdatePhaseMobilityKernel(7,3);
LAUNCH_UpdatePhaseMobilityKernel(8,3);
LAUNCH_UpdatePhaseMobilityKernel(9,3);
LAUNCH_UpdatePhaseMobilityKernel(10,3);

#undef LAUNCH_UpdatePhaseMobilityKernel

/******************************** Kernel launch machinery ********************************/

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

template<typename KERNELWRAPPER, typename... ARGS>
void KernelLaunchSelector1( localIndex numComp, ARGS && ... args )
{
  bool const run =
  KernelLaunchSelectorCompSwitch( numComp, [&] (auto NC)
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
  KernelLaunchSelectorCompSwitch( numComp, [ &numPhase, &run2, &args... ] ( auto NC )
  {
    run2 =
    KernelLaunchSelectorPhaseSwitch( numPhase, [&] ( auto NP )
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
