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
 * @file HydrostaticPressureKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_HYDROSTATICPRESSUREKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_HYDROSTATICPRESSUREKERNEL_HPP

#include "common/DataTypes.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "functions/TableFunction.hpp"

namespace geos
{

namespace isothermalCompositionalMultiphaseBaseKernels
{

/******************************** HydrostaticPressureKernel ********************************/

template< typename FLUID_WRAPPER >
struct HydrostaticPressureKernel
{
  enum class ReturnType : integer
  {
    FAILED_TO_CONVERGE = 0,
    DETECTED_MULTIPHASE_FLOW = 1,
    SUCCESS = 2,
    DETECTED_SINGLEPHASE_FLOW = 3,
    PHASE_CORRECTION_NOT_NEEDED = 4
  };

  static ReturnType
  computeHydrostaticPressure( integer const & numComps,
                              integer const & numPhases,
                              integer const & ipGas,
                              integer const & ipOil,
                              integer const & ipWater,
                              integer const & ipInit,
                              arrayView1d< real64 const > const & phaseContacts,
                              arrayView1d< real64 const > const & phaseMinVolumeFraction,
                              integer const maxNumEquilIterations,
                              real64 const & equilTolerance,
                              real64 const (&gravVector)[ 3 ],
                              FLUID_WRAPPER fluidWrapper,
                              arrayView1d< TableFunction::KernelWrapper const > compFracTableWrappers,
                              TableFunction::KernelWrapper tempTableWrapper,
                              real64 const & refElevation,
                              arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > const & refPres,
                              arraySlice1d< real64 const > const & refPhaseMassDens,
                              real64 const & newElevation,
                              arraySlice1d< real64, constitutive::multifluid::USD_PHASE - 2 > const & newPres,
                              arraySlice1d< real64 > const & newPhaseMassDens,
                              arraySlice1d< real64, constitutive::multifluid::USD_PHASE - 2 > const & newPhaseDens,
                              arraySlice2d< real64, constitutive::multifluid::USD_PHASE_COMP - 2 > const & newPhaseCompFrac );

  static ReturnType
  computeHydrostaticPressureAtMultipleElevations( localIndex const & startElevationIndex,
                                                  localIndex const & endElevationIndex,
                                                  integer const & numComps,
                                                  integer const & numPhases,
                                                  integer const & ipGas,
                                                  integer const & ipOil,
                                                  integer const & ipWater,
                                                  integer const & ipInit,
                                                  arrayView1d< real64 const > const & phaseContacts,
                                                  arrayView1d< real64 const > const & phaseMinVolumeFraction,
                                                  integer const & maxNumEquilIterations,
                                                  real64 const & equilTolerance,
                                                  real64 const (&gravVector)[ 3 ],
                                                  FLUID_WRAPPER fluidWrapper,
                                                  arrayView1d< TableFunction::KernelWrapper const > compFracTableWrappers,
                                                  TableFunction::KernelWrapper tempTableWrapper,
                                                  arrayView1d< arrayView1d< real64 > const > elevationValues,
                                                  arrayView3d< real64, constitutive::multifluid::USD_PHASE > const & pressureValues,
                                                  arrayView2d< real64 > const & phaseMassDens,
                                                  arrayView3d< real64, constitutive::multifluid::USD_PHASE > const & phaseDens,
                                                  arrayView4d< real64, constitutive::multifluid::USD_PHASE_COMP > const & phaseCompFrac );

  static ReturnType
  marchBetweenTwoElevations( real64 const & startElevation,
                             real64 const & endElevation,
                             integer const & numComps,
                             integer const & numPhases,
                             integer const & ipGas,
                             integer const & ipOil,
                             integer const & ipWater,
                             integer const & ipInit,
                             arrayView1d< real64 const > const & phaseContacts,
                             arrayView1d< real64 const > const & phaseMinVolumeFraction,
                             integer const maxNumEquilIterations,
                             real64 const & equilTolerance,
                             real64 const (&gravVector)[ 3 ],
                             FLUID_WRAPPER fluidWrapper,
                             arrayView1d< TableFunction::KernelWrapper const > compFracTableWrappers,
                             TableFunction::KernelWrapper tempTableWrapper,
                             arrayView1d< arrayView1d< real64 > const > elevationValues,
                             arrayView3d< real64, constitutive::multifluid::USD_PHASE > const & pressureValues,
                             arrayView2d< real64 > const & phaseMassDens,
                             arrayView3d< real64, constitutive::multifluid::USD_PHASE > const & phaseDens,
                             arrayView4d< real64, constitutive::multifluid::USD_PHASE_COMP > const & phaseCompFrac );

  static ReturnType
  launch( localIndex const & size,
          integer const & numComps,
          integer const & numPhases,
          integer const & ipGas,
          integer const & ipOil,
          integer const & ipWater,
          integer const & ipInit,
          integer const & maxNumEquilIterations,
          arrayView1d< real64 const > const & phaseContacts,
          arrayView1d< real64 const > const & phaseMinVolumeFraction,
          real64 const equilTolerance,
          real64 const (&gravVector)[ 3 ],
          real64 const & datumElevation,
          real64 const & datumPres,
          FLUID_WRAPPER fluidWrapper,
          arrayView1d< TableFunction::KernelWrapper const > compFracTableWrappers,
          TableFunction::KernelWrapper tempTableWrapper,
          arrayView1d< arrayView1d< real64 > const > elevationValues,
          arrayView3d< real64, constitutive::multifluid::USD_PHASE > const & pressureValues,
          arrayView3d< real64, constitutive::multifluid::USD_PHASE > const & phaseDens,
          arrayView4d< real64, constitutive::multifluid::USD_PHASE_COMP > const & phaseCompFrac );

  static ReturnType
  phaseCorrection( integer const & numComps,
                   integer const & numPhases,
                   integer const & ipGas,
                   integer const & ipWater,
                   arrayView1d< real64 const > const & phaseMinVolumeFraction,
                   real64 const & pres,
                   real64 const & temp,
                   arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & inputComposition,
                   arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > const & phaseFrac,
                   arraySlice1d< real64, compflow::USD_COMP - 1 > const & compFrac,
                   FLUID_WRAPPER fluidWrapper );

  static ReturnType
  applyPhaseCorrection( integer const & numComps,
                        integer const & numPhases,
                        integer const & ip_phase,
                        integer const & ip_otherPhase,
                        real64 const & pres,
                        real64 const & temp,
                        arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & inputComposition,
                        arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & addedCompFrac,
                        arraySlice1d< real64, compflow::USD_COMP - 1 > const & compFrac,
                        FLUID_WRAPPER fluidWrapper );

  static void mixingStep( integer const & numComps,
                          real64 const & a,
                          arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & inputComposition,
                          arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & addedCompFrac,
                          arraySlice1d< real64, compflow::USD_COMP - 1 > const & compFrac );

  static integer evaluateFlashPhaseIndex( integer const & numPhases,
                                          integer const & ipGas,
                                          integer const & ipOil,
                                          integer const & ipWater,
                                          integer const & ipInit,
                                          real64 const & elevation,
                                          arrayView1d< real64 const > const & phaseContacts );

  static void evaluatePrimaryAndContactPhaseIndices( integer const & numPhases,
                                                     integer const & ipGas,
                                                     integer const & ipOil,
                                                     integer const & ipWater,
                                                     real64 const & startElevation,
                                                     real64 const & endElevation,
                                                     arrayView1d< real64 const > const & phaseContacts,
                                                     integer & ipPP,
                                                     integer & ipCP );

  static void
  computeDatumPhaseMassDens( integer const & numComps,
                             integer const & numPhases,
                             integer const & ipGas,
                             integer const & ipWater,
                             arrayView1d< real64 const > const & phaseMinVolumeFraction,
                             real64 const & datumElevation,
                             real64 const & datumPres,
                             real64 const & datumTemp,
                             arrayView3d< real64, constitutive::multifluid::USD_PHASE > const & datumPhaseMassDens,
                             arrayView1d< TableFunction::KernelWrapper const > compFracTableWrappers,
                             bool singlePhase,
                             FLUID_WRAPPER fluidWrapper );
};

} // namespace isothermalCompositionalMultiphaseBaseKernels
} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_HYDROSTATICPRESSUREKERNEL_HPP
