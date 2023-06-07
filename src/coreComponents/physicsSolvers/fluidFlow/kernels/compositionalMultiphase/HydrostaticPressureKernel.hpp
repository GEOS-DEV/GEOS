/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file HydrostaticPressureKernel.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_HYDROSTATICPRESSUREKERNEL_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_HYDROSTATICPRESSUREKERNEL_HPP


namespace geos
{

namespace isothermalCompositionalMultiphaseBaseKernels
{

/******************************** HydrostaticPressureKernel ********************************/

struct HydrostaticPressureKernel
{

  // TODO: this type of constants should be centralized somewhere or provided by fluid model
  static geos::real64 constexpr MIN_FOR_PHASE_PRESENCE = 1e-12;

  enum class ReturnType : geos::integer
  {
    FAILED_TO_CONVERGE = 0,
    DETECTED_MULTIPHASE_FLOW = 1,
    SUCCESS = 2
  };

  template< typename FLUID_WRAPPER >
  static ReturnType
  computeHydrostaticPressure( geos::integer const numComps,
                              geos::integer const numPhases,
                              geos::integer const ipInit,
                              geos::integer const maxNumEquilIterations,
                              geos::real64 const & equilTolerance,
                              geos::real64 const (&gravVector)[3],
                              FLUID_WRAPPER fluidWrapper,
                              geos::arrayView1d< geos::TableFunction::KernelWrapper const > compFracTableWrappers,
                              geos::TableFunction::KernelWrapper tempTableWrapper,
                              geos::real64 const & refElevation,
                              geos::real64 const & refPres,
                              geos::arraySlice1d< geos::real64 const > const & refPhaseMassDens,
                              geos::real64 const & newElevation,
                              geos::real64 & newPres,
                              geos::arraySlice1d< geos::real64 > const & newPhaseMassDens )
  {
    // fluid properties at this elevation
    geos::StackArray< geos::real64, 2, geos::constitutive::MultiFluidBase::MAX_NUM_COMPONENTS, geos::compflow::LAYOUT_COMP > compFrac( 1, numComps );
    geos::StackArray< geos::real64, 3, geos::constitutive::MultiFluidBase::MAX_NUM_PHASES, geos::constitutive::multifluid::LAYOUT_PHASE > phaseFrac( 1, 1, numPhases );
    geos::StackArray< geos::real64, 3, geos::constitutive::MultiFluidBase::MAX_NUM_PHASES, geos::constitutive::multifluid::LAYOUT_PHASE > phaseDens( 1, 1, numPhases );
    geos::StackArray< geos::real64, 3, geos::constitutive::MultiFluidBase::MAX_NUM_PHASES, geos::constitutive::multifluid::LAYOUT_PHASE > phaseMassDens( 1, 1, numPhases );
    geos::StackArray< geos::real64, 3, geos::constitutive::MultiFluidBase::MAX_NUM_PHASES, geos::constitutive::multifluid::LAYOUT_PHASE > phaseVisc( 1, 1, numPhases );
    geos::StackArray< geos::real64, 3, geos::constitutive::MultiFluidBase::MAX_NUM_PHASES, geos::constitutive::multifluid::LAYOUT_PHASE > phaseEnthalpy( 1, 1, numPhases );
    geos::StackArray< geos::real64, 3, geos::constitutive::MultiFluidBase::MAX_NUM_PHASES, geos::constitutive::multifluid::LAYOUT_PHASE > phaseInternalEnergy( 1, 1, numPhases );
    geos::StackArray< geos::real64, 4, geos::constitutive::MultiFluidBase::MAX_NUM_PHASES * geos::constitutive::MultiFluidBase::MAX_NUM_COMPONENTS,
                      geos::constitutive::multifluid::LAYOUT_PHASE_COMP > phaseCompFrac( 1, 1, numPhases, numComps );
    geos::real64 totalDens = 0.0;

    bool isSinglePhaseFlow = true;

    // Step 1: compute the hydrostatic pressure at the current elevation

    geos::real64 const gravCoef = gravVector[2] * ( refElevation - newElevation );
    geos::real64 const temp = tempTableWrapper.compute( &newElevation );
    for( geos::integer ic = 0; ic < numComps; ++ic )
    {
      compFrac[0][ic] = compFracTableWrappers[ic].compute( &newElevation );
    }

    // Step 2: guess the pressure with the refPhaseMassDensity

    geos::real64 pres0 = refPres - refPhaseMassDens[ipInit] * gravCoef;
    geos::real64 pres1 = 0.0;

    // Step 3: compute the mass density at this elevation using the guess, and update pressure

    fluidWrapper.compute( pres0,
                          temp,
                          compFrac[0],
                          phaseFrac[0][0],
                          phaseDens[0][0],
                          phaseMassDens[0][0],
                          phaseVisc[0][0],
                          phaseEnthalpy[0][0],
                          phaseInternalEnergy[0][0],
                          phaseCompFrac[0][0],
                          totalDens );
    pres1 = refPres - 0.5 * ( refPhaseMassDens[ipInit] + phaseMassDens[0][0][ipInit] ) * gravCoef;

    // Step 4: fixed-point iteration until convergence

    bool equilHasConverged = false;
    for( geos::integer eqIter = 0; eqIter < maxNumEquilIterations; ++eqIter )
    {

      // check convergence
      equilHasConverged = ( LvArray::math::abs( pres0 - pres1 ) < equilTolerance );
      pres0 = pres1;

      // if converged, check number of phases and move on
      if( equilHasConverged )
      {
        // make sure that the fluid is single-phase, other we have to issue a warning (for now)
        // if only one phase is mobile, we are in good shape (unfortunately it is hard to access relperm from here)
        geos::localIndex numberOfPhases = 0;
        for( geos::integer ip = 0; ip < numPhases; ++ip )
        {
          if( phaseFrac[0][0][ip] > MIN_FOR_PHASE_PRESENCE )
          {
            numberOfPhases++;
          }
        }
        if( numberOfPhases > 1 )
        {
          isSinglePhaseFlow = false;
        }

        break;
      }

      // compute the mass density at this elevation using the previous pressure, and compute the new pressure
      fluidWrapper.compute( pres0,
                            temp,
                            compFrac[0],
                            phaseFrac[0][0],
                            phaseDens[0][0],
                            phaseMassDens[0][0],
                            phaseVisc[0][0],
                            phaseEnthalpy[0][0],
                            phaseInternalEnergy[0][0],
                            phaseCompFrac[0][0],
                            totalDens );
      pres1 = refPres - 0.5 * ( refPhaseMassDens[ipInit] + phaseMassDens[0][0][ipInit] ) * gravCoef;
    }

    // Step 5: save the hydrostatic pressure and the corresponding density

    newPres = pres1;
    for( geos::integer ip = 0; ip < numPhases; ++ip )
    {
      newPhaseMassDens[ip] = phaseMassDens[0][0][ip];
    }

    if( !equilHasConverged )
    {
      return ReturnType::FAILED_TO_CONVERGE;
    }
    else if( !isSinglePhaseFlow )
    {
      return ReturnType::DETECTED_MULTIPHASE_FLOW;
    }
    else
    {
      return ReturnType::SUCCESS;
    }
  }

  template< typename FLUID_WRAPPER >
  static ReturnType
  launch( geos::localIndex const size,
          geos::integer const numComps,
          geos::integer const numPhases,
          geos::integer const ipInit,
          geos::integer const maxNumEquilIterations,
          geos::real64 const equilTolerance,
          geos::real64 const (&gravVector)[3],
          geos::real64 const & minElevation,
          geos::real64 const & elevationIncrement,
          geos::real64 const & datumElevation,
          geos::real64 const & datumPres,
          FLUID_WRAPPER fluidWrapper,
          geos::arrayView1d< geos::TableFunction::KernelWrapper const > compFracTableWrappers,
          geos::TableFunction::KernelWrapper tempTableWrapper,
          geos::arrayView1d< geos::arrayView1d< geos::real64 > const > elevationValues,
          geos::arrayView1d< geos::real64 > pressureValues )
  {

    ReturnType returnVal = ReturnType::SUCCESS;

    // Step 1: compute the phase mass densities at datum

    // datum fluid properties
    geos::array2d< geos::real64, geos::compflow::LAYOUT_COMP > datumCompFrac( 1, numComps );
    geos::array3d< geos::real64, geos::constitutive::multifluid::LAYOUT_PHASE > datumPhaseFrac( 1, 1, numPhases );
    geos::array3d< geos::real64, geos::constitutive::multifluid::LAYOUT_PHASE > datumPhaseDens( 1, 1, numPhases );
    geos::array3d< geos::real64, geos::constitutive::multifluid::LAYOUT_PHASE > datumPhaseMassDens( 1, 1, numPhases );
    geos::array3d< geos::real64, geos::constitutive::multifluid::LAYOUT_PHASE > datumPhaseVisc( 1, 1, numPhases );
    geos::array3d< geos::real64, geos::constitutive::multifluid::LAYOUT_PHASE > datumPhaseEnthalpy( 1, 1, numPhases );
    geos::array3d< geos::real64, geos::constitutive::multifluid::LAYOUT_PHASE > datumPhaseInternalEnergy( 1, 1, numPhases );
    geos::array4d< geos::real64, geos::constitutive::multifluid::LAYOUT_PHASE_COMP > datumPhaseCompFrac( 1, 1, numPhases, numComps );
    geos::real64 datumTotalDens = 0.0;

    geos::real64 const datumTemp = tempTableWrapper.compute( &datumElevation );
    for( geos::integer ic = 0; ic < numComps; ++ic )
    {
      datumCompFrac[0][ic] = compFracTableWrappers[ic].compute( &datumElevation );
    }
    fluidWrapper.compute( datumPres,
                          datumTemp,
                          datumCompFrac[0],
                          datumPhaseFrac[0][0],
                          datumPhaseDens[0][0],
                          datumPhaseMassDens[0][0],
                          datumPhaseVisc[0][0],
                          datumPhaseEnthalpy[0][0],
                          datumPhaseInternalEnergy[0][0],
                          datumPhaseCompFrac[0][0],
                          datumTotalDens );

    // Step 2: find the closest elevation to datumElevation

    geos::forAll< geos::parallelHostPolicy >( size, [=]( geos::localIndex const i )
    {
      geos::real64 const elevation = minElevation + i * elevationIncrement;
      elevationValues[0][i] = elevation;
    } );
    geos::integer const iRef = LvArray::sortedArrayManipulation::find( elevationValues[0].begin(),
                                                                       elevationValues[0].size(),
                                                                       datumElevation );

    // Step 3: compute the mass density and pressure at the reference elevation

    geos::array2d< geos::real64 > phaseMassDens( pressureValues.size(), numPhases );
    // temporary array without permutation to compile on Lassen
    geos::array1d< geos::real64 > datumPhaseMassDensTmp( numPhases );
    for( geos::integer ip = 0; ip < numPhases; ++ip )
    {
      datumPhaseMassDensTmp[ip] = datumPhaseMassDens[0][0][ip];
    }

    ReturnType const refReturnVal =
      computeHydrostaticPressure( numComps,
                                  numPhases,
                                  ipInit,
                                  maxNumEquilIterations,
                                  equilTolerance,
                                  gravVector,
                                  fluidWrapper,
                                  compFracTableWrappers,
                                  tempTableWrapper,
                                  datumElevation,
                                  datumPres,
                                  datumPhaseMassDensTmp,
                                  elevationValues[0][iRef],
                                  pressureValues[iRef],
                                  phaseMassDens[iRef] );
    if( refReturnVal == ReturnType::FAILED_TO_CONVERGE )
    {
      return ReturnType::FAILED_TO_CONVERGE;
    }
    else if( refReturnVal == ReturnType::DETECTED_MULTIPHASE_FLOW )
    {
      returnVal = ReturnType::DETECTED_MULTIPHASE_FLOW;
    }

    // Step 4: for each elevation above the reference elevation, compute the pressure

    geos::localIndex const numEntriesAboveRef = size - iRef - 1;
    geos::forAll< geos::serialPolicy >( numEntriesAboveRef, [=, &returnVal]( geos::localIndex const i )
    {
      ReturnType const returnValAboveRef =
        computeHydrostaticPressure( numComps,
                                    numPhases,
                                    ipInit,
                                    maxNumEquilIterations,
                                    equilTolerance,
                                    gravVector,
                                    fluidWrapper,
                                    compFracTableWrappers,
                                    tempTableWrapper,
                                    elevationValues[0][iRef + i],
                                    pressureValues[iRef + i],
                                    phaseMassDens[iRef + i],
                                    elevationValues[0][iRef + i + 1],
                                    pressureValues[iRef + i + 1],
                                    phaseMassDens[iRef + i + 1] );
      if( returnValAboveRef == ReturnType::FAILED_TO_CONVERGE )
      {
        returnVal = ReturnType::FAILED_TO_CONVERGE;
      }
      else if( ( returnValAboveRef == ReturnType::DETECTED_MULTIPHASE_FLOW ) &&
               ( returnVal != ReturnType::FAILED_TO_CONVERGE ) )
      {
        returnVal = ReturnType::DETECTED_MULTIPHASE_FLOW;
      }

    } );

    // Step 5: for each elevation below the reference elevation, compute the pressure

    geos::localIndex const numEntriesBelowRef = iRef;
    geos::forAll< geos::serialPolicy >( numEntriesBelowRef, [=, &returnVal]( geos::localIndex const i )
    {
      ReturnType const returnValBelowRef =
        computeHydrostaticPressure( numComps,
                                    numPhases,
                                    ipInit,
                                    maxNumEquilIterations,
                                    equilTolerance,
                                    gravVector,
                                    fluidWrapper,
                                    compFracTableWrappers,
                                    tempTableWrapper,
                                    elevationValues[0][iRef - i],
                                    pressureValues[iRef - i],
                                    phaseMassDens[iRef - i],
                                    elevationValues[0][iRef - i - 1],
                                    pressureValues[iRef - i - 1],
                                    phaseMassDens[iRef - i - 1] );
      if( returnValBelowRef == ReturnType::FAILED_TO_CONVERGE )
      {
        returnVal = ReturnType::FAILED_TO_CONVERGE;
      }
      else if( ( returnValBelowRef == ReturnType::DETECTED_MULTIPHASE_FLOW ) &&
               ( returnVal != ReturnType::FAILED_TO_CONVERGE ) )
      {
        returnVal = ReturnType::DETECTED_MULTIPHASE_FLOW;
      }

    } );

    return returnVal;
  }

};

}

}

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_HYDROSTATICPRESSUREKERNEL_HPP
