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

#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "functions/TableFunction.hpp"


#include <fstream>
#include <iostream>


namespace geos
{

namespace isothermalCompositionalMultiphaseBaseKernels
{

/******************************** HydrostaticPressureKernel ********************************/

struct HydrostaticPressureKernel
{

  // TODO: this type of constants should be centralized somewhere or provided by fluid model
  static real64 constexpr MIN_FOR_PHASE_PRESENCE = 1e-12;
  static real64 constexpr MIN_FOR_COMP_PRESENCE = 1e-12;

  enum class ReturnType : integer
  {
    FAILED_TO_CONVERGE = 0,
    DETECTED_MULTIPHASE_FLOW = 1,
    SUCCESS = 2,
    DETECTED_SINGLEPHASE_FLOW = 3,
    PHASE_CORRECTION_NOT_NEEDED = 4
  };

  template< typename FLUID_WRAPPER >
  static ReturnType
  computeHydrostaticPressure( integer const & numComps,
                              integer const & numPhases,
                              integer const & ip_gas,
                              integer const & ip_oil,
                              integer const & ip_water,
                              integer const & ip_pp,
                              arrayView1d< real64 const > const & phaseContacts,
                              arrayView1d< real64 const > const & phaseMinVolumeFraction,
                              integer const maxNumEquilIterations,
                              real64 const & equilTolerance,
                              real64 const (&gravVector)[ 3 ],
                              FLUID_WRAPPER fluidWrapper,
                              arrayView1d< TableFunction::KernelWrapper const > compFracTableWrappers,
                              TableFunction::KernelWrapper tempTableWrapper,
                              real64 const & refElevation,
                              arraySlice1d< real64 const , constitutive::multifluid::USD_PHASE - 2 > const & refPres,
                              arraySlice1d< real64 const > const & refPhaseMassDens,
                              real64 const & newElevation,
                              arraySlice1d< real64, constitutive::multifluid::USD_PHASE - 2 > const & newPres,
                              arraySlice1d< real64 > const & newPhaseMassDens,
                              arraySlice1d< real64, constitutive::multifluid::USD_PHASE - 2 > const & newPhaseDens,
                              arraySlice2d< real64, constitutive::multifluid::USD_PHASE_COMP - 2 > const & newPhaseCompFrac,
                              localIndex const index )
  {
    // fluid properties at this elevation
    StackArray< real64, 2, constitutive::MultiFluidBase::MAX_NUM_COMPONENTS, compflow::LAYOUT_COMP > uncorrCompFrac( 1, numComps );
    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, constitutive::multifluid::LAYOUT_PHASE > phaseFrac( 1, 1, numPhases );
    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, constitutive::multifluid::LAYOUT_PHASE > phaseDens( 1, 1, numPhases );
    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, constitutive::multifluid::LAYOUT_PHASE > phaseMassDens( 1, 1, numPhases );
    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, constitutive::multifluid::LAYOUT_PHASE > phaseVisc( 1, 1, numPhases );
    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, constitutive::multifluid::LAYOUT_PHASE > phaseEnthalpy( 1, 1, numPhases );
    StackArray< real64, 3, constitutive::MultiFluidBase::MAX_NUM_PHASES, constitutive::multifluid::LAYOUT_PHASE > phaseInternalEnergy( 1, 1, numPhases );
    StackArray< real64, 4, constitutive::MultiFluidBase::MAX_NUM_PHASES *constitutive::MultiFluidBase::MAX_NUM_COMPONENTS,
                constitutive::multifluid::LAYOUT_PHASE_COMP > phaseCompFrac( 1, 1, numPhases, numComps );
    real64 totalDens = 0.0;

    bool isSinglePhaseFlow = true;

    // Step 1: compute the hydrostatic pressure at the current elevation

    real64 const gravCoef = gravVector[2] * ( refElevation - newElevation );
    real64 const temp = tempTableWrapper.compute( &newElevation );
    for( integer ic = 0; ic < numComps; ++ic )
    {
      uncorrCompFrac[0][ic] = compFracTableWrappers[ic].compute( &newElevation );
    }

    // Step 2: guess the pressure with the refPhaseMassDensity
    array1d< real64 > pres0( numPhases );
    array1d< real64 > pres1( numPhases );
    for ( localIndex ip = 0; ip < numPhases; ++ip )
    {
      pres0[ip] = refPres[ip] - refPhaseMassDens[ip] * gravCoef;
      pres1[ip] = 0.0;
    }

    // Step 3: determine which phase pressure to use as the flash pressure
    // integer const ip_pres = primaryPressure( numPhases,
    //                                          ip_gas,
    //                                          ip_oil,
    //                                          ip_water,
    //                                          newElevation,
    //                                          phaseContacts );

    // Step 3: compute the mass density at this elevation using the guess, and update pressure

    constitutive::MultiFluidBase::KernelWrapper::computeValues( fluidWrapper,
                                                                pres0[ip_pp], // flash pressure
                                                                temp,
                                                                uncorrCompFrac[0],
                                                                phaseFrac[0][0],
                                                                phaseDens[0][0],
                                                                phaseMassDens[0][0],
                                                                phaseVisc[0][0],
                                                                phaseEnthalpy[0][0],
                                                                phaseInternalEnergy[0][0],
                                                                phaseCompFrac[0][0],
                                                                totalDens );
    
    // Step 4: Ensure the correct phases exist. If not, apply phase correction.
    StackArray< real64, 2, constitutive::MultiFluidBase::MAX_NUM_COMPONENTS, compflow::LAYOUT_COMP > compFrac( 1, numComps );
    ReturnType phaseCorr = phaseCorrection( numComps,
                                            numPhases,
                                            ip_gas,
                                            ip_water,
                                            phaseMinVolumeFraction,
                                            pres0[ip_pp],
                                            temp,
                                            uncorrCompFrac[0],
                                            phaseFrac[0][0],
                                            compFrac[0],
                                            fluidWrapper );
    if ( phaseCorr == ReturnType::SUCCESS )
      std::cout << "phase correction function completed" << std::endl;


    for ( localIndex ip = 0; ip < numPhases; ++ip )
    {
      pres1[ip] = refPres[ip] - 0.5 * ( refPhaseMassDens[ip] + phaseMassDens[0][0][ip] ) * gravCoef;
    }

    // Step 4: fixed-point iteration until convergence

    // bool equilHasConverged = false;
    bool equilHasConverged;
    int iters = 0;
    for( integer eqIter = 0; eqIter < maxNumEquilIterations; ++eqIter )
    {
      iters += 1;
      // check convergence
      // equilHasConverged = ( LvArray::math::abs( pres0 - pres1 ) < equilTolerance );
      // pres0 = pres1;
      equilHasConverged = true;
      for ( localIndex ip = 0; ip < numPhases; ++ip )
      {
        equilHasConverged = equilHasConverged && ( LvArray::math::abs( pres0[ip] - pres1[ip] ) < equilTolerance );
        pres0[ip] = pres1[ip];
      }

      // if converged, check number of phases and move on
      if( equilHasConverged )
      {
        // make sure that the fluid is single-phase, other we have to issue a warning (for now)
        // if only one phase is mobile, we are in good shape (unfortunately it is hard to access relperm from here)
        localIndex numberOfPhases = 0;
        for( integer ip = 0; ip < numPhases; ++ip )
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
      constitutive::MultiFluidBase::KernelWrapper::computeValues( fluidWrapper,
                                                                  pres0[ip_pp], // flash pressure
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
      for ( localIndex ip = 0; ip < numPhases; ++ip )
      {
        pres1[ip] = refPres[ip] - 0.5 * ( refPhaseMassDens[ip] + phaseMassDens[0][0][ip] ) * gravCoef;
      }
    }

    // Step 5: save the hydrostatic pressure and the corresponding density

    // newPres[1] = pres1;
    for ( localIndex ip = 0; ip < numPhases; ++ip )
    {
      newPres[ip] = pres1[ip];
    }

    if( equilHasConverged )
    {
      std::cout << "At index = " << index << ", newPres[0] = " << newPres[0] << ", newPres[1] = " 
                << newPres[1] << " at elevation = " << newElevation << ", eqIters = " << iters 
                << ", ip_oil = " << ip_oil << ", phaseContacts = " << phaseContacts << std::endl;
    }
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      newPhaseMassDens[ip] = phaseMassDens[0][0][ip];
      newPhaseDens[ip] = phaseDens[0][0][ip];
      for ( integer ic = 0; ic < numComps; ++ic )
      {
        newPhaseCompFrac[ip][ic] = phaseCompFrac[0][0][ip][ic];
      }
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
  computeHydrostaticPressureAtMultipleElevations( localIndex const & startElevationIndex,
                                                  localIndex const & endElevationIndex,
                                                  integer const & numComps,
                                                  integer const & numPhases,
                                                  integer const & ip_gas,
                                                  integer const & ip_oil,
                                                  integer const & ip_water,
                                                  integer const & ip_pp,
                                                  arrayView1d< real64 const > const & phaseContacts,
                                                  arrayView1d< real64 const > const & phaseMinVolumeFraction,
                                                  integer const & maxNumEquilIterations,
                                                  real64 const & equilTolerance,
                                                  real64 const (&gravVector)[ 3 ],
                                                  FLUID_WRAPPER fluidWrapper,
                                                  arrayView1d< TableFunction::KernelWrapper const > compFracTableWrappers,
                                                  TableFunction::KernelWrapper tempTableWrapper,
                                                  arrayView1d< arrayView1d< real64 > const > elevationValues,
                                                  // arrayView2d< real64 > const pressureValues,
                                                  arrayView3d< real64, constitutive::multifluid::USD_PHASE > const & pressureValues,
                                                  arrayView2d< real64 > const & phaseMassDens,
                                                  arrayView3d< real64, constitutive::multifluid::USD_PHASE > const & phaseDens,
                                                  arrayView4d< real64, constitutive::multifluid::USD_PHASE_COMP > const & phaseCompFrac )
  {
    // startElevIndex is the reference point
    localIndex const numEntries = LvArray::math::abs( startElevationIndex - endElevationIndex );
    localIndex const step = ( endElevationIndex >= startElevationIndex ) ? 1 : -1;
    ReturnType returnVal = ReturnType::SUCCESS;
    forAll< serialPolicy >( numEntries, [=, &returnVal] (localIndex const i)
    {
      localIndex const ref = startElevationIndex + i * step;
      localIndex const next = ref + step;
      ReturnType const iReturnVal =
        computeHydrostaticPressure( numComps,
                                    numPhases,
                                    ip_gas,
                                    ip_oil,
                                    ip_water,
                                    ip_pp,
                                    phaseContacts,
                                    phaseMinVolumeFraction,
                                    maxNumEquilIterations,
                                    equilTolerance,
                                    gravVector,
                                    fluidWrapper,
                                    compFracTableWrappers,
                                    tempTableWrapper,
                                    elevationValues[0][ref],
                                    pressureValues[ref][0],
                                    phaseMassDens[ref],
                                    elevationValues[0][next],
                                    pressureValues[next][0],
                                    phaseMassDens[next],
                                    phaseDens[next][0],
                                    phaseCompFrac[next][0],
                                    next );
      // std::cout << "ForAll: ref_pres = " << pressureValues[ref][0] << ", newPres = " << pressureValues[next][0] << std::endl;
      if( iReturnVal == ReturnType::FAILED_TO_CONVERGE )
      {
        returnVal = ReturnType::FAILED_TO_CONVERGE;
      }
      else if( ( iReturnVal == ReturnType::DETECTED_MULTIPHASE_FLOW ) &&
               ( returnVal != ReturnType::FAILED_TO_CONVERGE ) )
      {
        returnVal = ReturnType::DETECTED_MULTIPHASE_FLOW;
      }

    } );

    return returnVal;
  }

  template< typename FLUID_WRAPPER >
  static ReturnType
  marchBetweenTwoElevations( localIndex const & size,
                             real64 const & startElevation,
                             real64 const & endElevation,
                             integer const & numComps,
                             integer const & numPhases,
                             integer const & ip_gas,
                             integer const & ip_oil,
                             integer const & ip_water,
                             integer const & ip_pp,
                             integer const & ip_cp,
                             integer const & ip_op,
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
                             arrayView4d< real64, constitutive::multifluid::USD_PHASE_COMP > const & phaseCompFrac )
  {
    std::cout << "In marchBetweenTwoElevations, ip_op = " << ip_op << std::endl;
    // Find index of the closest element in elevationValues to the start elevation
    localIndex const iStart = LvArray::sortedArrayManipulation::find( elevationValues[0].begin(),
                                                                      elevationValues[0].size(),
                                                                      startElevation );
    // Find index of the closest element in elevationValues to the end elevation
    localIndex const iEnd = LvArray::sortedArrayManipulation::find( elevationValues[0].begin(),
                                                                    elevationValues[0].size(),
                                                                    endElevation );
    // March from iStart to iEnd
    ReturnType returnVal = 
      computeHydrostaticPressureAtMultipleElevations( iStart,
                                                      iEnd,
                                                      numComps,
                                                      numPhases,
                                                      ip_gas,
                                                      ip_oil,
                                                      ip_water,
                                                      ip_pp,
                                                      phaseContacts,
                                                      phaseMinVolumeFraction,
                                                      maxNumEquilIterations,
                                                      equilTolerance,
                                                      gravVector,
                                                      fluidWrapper,
                                                      compFracTableWrappers,
                                                      tempTableWrapper,
                                                      elevationValues,
                                                      pressureValues,
                                                      phaseMassDens,
                                                      phaseDens,
                                                      phaseCompFrac );
    
    // Compute phase presssures and densities at end elevation using iEnd as the reference
    array2d< real64 > endPressure( 1, numPhases );
    array2d< real64 > endPhaseMassDens( 1, numPhases );
    array2d< real64 > endPhaseDens( 1, numPhases );
    array3d< real64 > endPhaseCompFrac( 1, numPhases, numComps );
    returnVal =
      computeHydrostaticPressure( numComps,
                                  numPhases,
                                  ip_gas,
                                  ip_oil,
                                  ip_water,
                                  ip_pp,
                                  phaseContacts,
                                  phaseMinVolumeFraction,
                                  maxNumEquilIterations,
                                  equilTolerance,
                                  gravVector,
                                  fluidWrapper,
                                  compFracTableWrappers,
                                  tempTableWrapper,
                                  elevationValues[0][iEnd],
                                  pressureValues[iEnd][0],
                                  phaseMassDens[iEnd],
                                  endElevation,
                                  endPressure[0],
                                  endPhaseMassDens[0],
                                  endPhaseDens[0],
                                  endPhaseCompFrac[0],
                                  iEnd );
    // Compute relative error defined as the relative difference between the phase pressures at end elevation                              
    real64 err = LvArray::math::abs( endPressure[0][ip_cp] - endPressure[0][ip_pp] ) / endPressure[0][ip_pp];
    std::cout << "march iter = 1: err = " << err << ", p[ip_pp] = " << endPressure[0][ip_pp] 
              << ", p[ip_cp] = " << endPressure[0][ip_cp] << std::endl;
    int maxMarchIterations = 10;

    // Marching Loop
    for ( int marchIter = 1; marchIter < maxMarchIterations; ++marchIter )
    {
      if ( err < 1e-5 ) // maybe use equilTolerance here as well and modify err from relative to absolute
      {
        std::cout << "March converged, march iter = " << marchIter << ": p[ip_pp] = " 
        << endPressure[0][ip_pp] << ", p[ip_cp] = " << endPressure[0][ip_cp] << std::endl;
        break;
      }
      // equate the phase pressure at the end elevation
      endPressure[0][ip_cp] = endPressure[0][ip_pp];
      real64 const iStartPrimaryPressure = pressureValues[iStart][0][ip_pp]; // saves the known phase pressure at iStart
      // Estimate the unknown phase (contact) pressure and density at iStart using the updated pressures
      returnVal =
        computeHydrostaticPressure( numComps,
                                    numPhases,
                                    ip_gas,
                                    ip_oil,
                                    ip_water,
                                    ip_pp,
                                    phaseContacts,
                                    phaseMinVolumeFraction,
                                    maxNumEquilIterations,
                                    equilTolerance,
                                    gravVector,
                                    fluidWrapper,
                                    compFracTableWrappers,
                                    tempTableWrapper,
                                    endElevation,
                                    endPressure[0],
                                    endPhaseMassDens[0],
                                    elevationValues[0][iStart],
                                    pressureValues[iStart][0],
                                    phaseMassDens[iStart],
                                    phaseDens[iStart][0],
                                    phaseCompFrac[iStart][0],
                                    iStart );
      pressureValues[iStart][0][ip_pp] = iStartPrimaryPressure;
      // March from iStart to iEnd
      returnVal = 
      computeHydrostaticPressureAtMultipleElevations( iStart,
                                                      iEnd,
                                                      numComps,
                                                      numPhases,
                                                      ip_gas,
                                                      ip_oil,
                                                      ip_water,
                                                      ip_pp,
                                                      phaseContacts,
                                                      phaseMinVolumeFraction,
                                                      maxNumEquilIterations,
                                                      equilTolerance,
                                                      gravVector,
                                                      fluidWrapper,
                                                      compFracTableWrappers,
                                                      tempTableWrapper,
                                                      elevationValues,
                                                      pressureValues,
                                                      phaseMassDens,
                                                      phaseDens,
                                                      phaseCompFrac );
      // Compute phase presssures and densities at the end elevation using iEnd as the reference
      returnVal =
        computeHydrostaticPressure( numComps,
                                    numPhases,
                                    ip_gas,
                                    ip_oil,
                                    ip_water,
                                    ip_pp,
                                    phaseContacts,
                                    phaseMinVolumeFraction,
                                    maxNumEquilIterations,
                                    equilTolerance,
                                    gravVector,
                                    fluidWrapper,
                                    compFracTableWrappers,
                                    tempTableWrapper,
                                    elevationValues[0][iEnd],
                                    pressureValues[iEnd][0],
                                    phaseMassDens[iEnd],
                                    endElevation,
                                    endPressure[0],
                                    endPhaseMassDens[0],
                                    endPhaseDens[0],
                                    endPhaseCompFrac[0],
                                    iEnd );
      err = LvArray::math::abs( endPressure[0][ip_cp] - endPressure[0][ip_pp] ) / endPressure[0][ip_pp];
      std::cout << "march iter = 1: err = " << err << ", p[ip_pp] = " << endPressure[0][ip_pp] 
                << ", p[ip_cp] = " << endPressure[0][ip_cp] << std::endl;
    }

    // Compute reference condition for points above iEnd if they exist in case where iEnd is above iStart
    if ( iEnd + 1 < size && iEnd > iStart )
    {
      returnVal =
        computeHydrostaticPressure( numComps,
                                    numPhases,
                                    ip_gas,
                                    ip_oil,
                                    ip_water,
                                    ip_pp,
                                    phaseContacts,
                                    phaseMinVolumeFraction,
                                    maxNumEquilIterations,
                                    equilTolerance,
                                    gravVector,
                                    fluidWrapper,
                                    compFracTableWrappers,
                                    tempTableWrapper,
                                    endElevation,
                                    endPressure[0],
                                    endPhaseMassDens[0],
                                    elevationValues[0][iEnd + 1],
                                    pressureValues[iEnd + 1][0],
                                    phaseMassDens[iEnd + 1],
                                    phaseDens[iEnd + 1][0],
                                    phaseCompFrac[iEnd + 1][0],
                                    iEnd+1 );
    }
    // Compute reference condition for points below iEnd if they exist in case where iEnd is below iStart
    else if ( iEnd + 1 > 0 && iEnd < iStart )
    {
      returnVal =
        computeHydrostaticPressure( numComps,
                                    numPhases,
                                    ip_gas,
                                    ip_oil,
                                    ip_water,
                                    ip_pp,
                                    phaseContacts,
                                    phaseMinVolumeFraction,
                                    maxNumEquilIterations,
                                    equilTolerance,
                                    gravVector,
                                    fluidWrapper,
                                    compFracTableWrappers,
                                    tempTableWrapper,
                                    endElevation,
                                    endPressure[0],
                                    endPhaseMassDens[0],
                                    elevationValues[0][iEnd - 1],
                                    pressureValues[iEnd - 1][0],
                                    phaseMassDens[iEnd - 1],
                                    phaseDens[iEnd - 1][0],
                                    phaseCompFrac[iEnd - 1][0],
                                    iEnd-1 );
    }

    return returnVal;
  }

  // template< typename FLUID_WRAPPER >
  // static ReturnType
  // marchBetweenContactAndDatum( localIndex const & size,
  //                              real64 const & contactElevation,
  //                              localIndex const & iDatum,
  //                              localIndex const & iContact,
  //                              integer const & numComps,
  //                              integer const & numPhases,
  //                              integer const & ip_gas,
  //                              integer const & ip_oil,
  //                              integer const & ip_water,
  //                              arrayView1d< real64 const > const & phaseContacts,
  //                              arrayView1d< real64 const > const & phaseMinVolumeFraction,
  //                              integer const maxNumEquilIterations,
  //                              real64 const & equilTolerance,
  //                              real64 const (&gravVector)[ 3 ],
  //                              FLUID_WRAPPER fluidWrapper,
  //                              arrayView1d< TableFunction::KernelWrapper const > compFracTableWrappers,
  //                              TableFunction::KernelWrapper tempTableWrapper,
  //                              arrayView1d< arrayView1d< real64 > const > elevationValues,
  //                             //  arrayView2d< real64 > const pressureValues,
  //                              arrayView3d< real64, constitutive::multifluid::USD_PHASE > const & pressureValues,
  //                              arrayView2d< real64 > const & phaseMassDens,
  //                              arrayView3d< real64, constitutive::multifluid::USD_PHASE > const & phaseDens,
  //                              arrayView4d< real64, constitutive::multifluid::USD_PHASE_COMP > const & phaseCompFrac )
  // {
  //   // March from iDatum to iContact
  //   ReturnType returnVal = 
  //     computeHydrostaticPressureAtMultipleElevations( iDatum,
  //                                                     iContact,
  //                                                     numComps,
  //                                                     numPhases,
  //                                                     ip_gas,
  //                                                     ip_oil,
  //                                                     ip_water,
  //                                                     phaseContacts,
  //                                                     phaseMinVolumeFraction,
  //                                                     maxNumEquilIterations,
  //                                                     equilTolerance,
  //                                                     gravVector,
  //                                                     fluidWrapper,
  //                                                     compFracTableWrappers,
  //                                                     tempTableWrapper,
  //                                                     elevationValues,
  //                                                     pressureValues,
  //                                                     phaseMassDens,
  //                                                     phaseDens,
  //                                                     phaseCompFrac );
  //   // Compute phase presssures and densities at the contact using iContact as the reference
  //   array2d< real64 > contactPressure( 1, numPhases );
  //   array2d< real64 > contactPhaseMassDens( 1, numPhases );
  //   array2d< real64 > contactPhaseDens( 1, numPhases );
  //   array3d< real64 > contactPhaseCompFrac( 1, numPhases, numComps );
  //   returnVal =
  //     computeHydrostaticPressure( numComps,
  //                                 numPhases,
  //                                 ip_gas,
  //                                 ip_oil,
  //                                 ip_water,
  //                                 phaseContacts,
  //                                 phaseMinVolumeFraction,
  //                                 maxNumEquilIterations,
  //                                 equilTolerance,
  //                                 gravVector,
  //                                 fluidWrapper,
  //                                 compFracTableWrappers,
  //                                 tempTableWrapper,
  //                                 elevationValues[0][iContact],
  //                                 pressureValues[iContact][0],
  //                                 phaseMassDens[iContact],
  //                                 contactElevation,
  //                                 contactPressure[0],
  //                                 contactPhaseMassDens[0],
  //                                 contactPhaseDens[0],
  //                                 contactPhaseCompFrac[0],
  //                                 iContact );
  //   // Compute relative error defined as the relative difference between the phase pressures at the contact                              
  //   real64 err = LvArray::math::abs( contactPressure[0][0] - contactPressure[0][1] ) / contactPressure[0][1];
  //   std::cout << "march iter = 1: err = " << err << ", gasPres = " << contactPressure[0][0] << ", liqPres = " << contactPressure[0][1] << std::endl;
  //   int maxMarchIterations = 10;

  //   // Marching Loop
  //   for ( int marchIter = 1; marchIter < maxMarchIterations; ++marchIter )
  //   {
  //     if ( err < 1e-5 ) // maybe use equilTolerance here as well and modify err from relative to absolute
  //     {
  //       std::cout << "March converged, march iter = " << marchIter << ": gasPres = " 
  //       << contactPressure[0][0] << ", liqPres = " << contactPressure[0][1] << std::endl;
  //       break;
  //     }
  //     // equate the phase pressure at the contact
  //     contactPressure[0][0] = contactPressure[0][1];
  //     // array2d< real64 > tmpRefPressure( 1, numPhases );
  //     // array2d< real64 > tmpRefPhaseMassDens( 1, numPhases );
  //     // array3d< real64 > tmpRefPhaseCompFrac( 1, numPhases, numComps );
  //     real64 const saveRefPhasePressure = pressureValues[iDatum][0][1]; // saves the known phase pressure at iDatum
  //     // Estimate the unknown phase pressure and density at iDatum using the updated contact phase pressures
  //     returnVal =
  //       computeHydrostaticPressure( numComps,
  //                                   numPhases,
  //                                   ip_gas,
  //                                   ip_oil,
  //                                   ip_water,
  //                                   phaseContacts,
  //                                   phaseMinVolumeFraction,
  //                                   maxNumEquilIterations,
  //                                   equilTolerance,
  //                                   gravVector,
  //                                   fluidWrapper,
  //                                   compFracTableWrappers,
  //                                   tempTableWrapper,
  //                                   contactElevation,
  //                                   contactPressure[0],
  //                                   contactPhaseMassDens[0],
  //                                   elevationValues[0][iDatum],
  //                                   pressureValues[iDatum][0],
  //                                   phaseMassDens[iDatum],
  //                                   phaseDens[iDatum][0],
  //                                   phaseCompFrac[iDatum][0],
  //                                   iDatum );
  //     pressureValues[iDatum][0][1] = saveRefPhasePressure;
  //     // March from iDatum to iContact 
  //     returnVal = 
  //     computeHydrostaticPressureAtMultipleElevations( iDatum,
  //                                                     iContact,
  //                                                     numComps,
  //                                                     numPhases,
  //                                                     ip_gas,
  //                                                     ip_oil,
  //                                                     ip_water,
  //                                                     phaseContacts,
  //                                                     phaseMinVolumeFraction,
  //                                                     maxNumEquilIterations,
  //                                                     equilTolerance,
  //                                                     gravVector,
  //                                                     fluidWrapper,
  //                                                     compFracTableWrappers,
  //                                                     tempTableWrapper,
  //                                                     elevationValues,
  //                                                     pressureValues,
  //                                                     phaseMassDens,
  //                                                     phaseDens,
  //                                                     phaseCompFrac );
  //   // Compute phase presssures and densities at the contact using iContact as the reference
  //   returnVal =
  //     computeHydrostaticPressure( numComps,
  //                                 numPhases,
  //                                 ip_gas,
  //                                 ip_oil,
  //                                 ip_water,
  //                                 phaseContacts,
  //                                 phaseMinVolumeFraction,
  //                                 maxNumEquilIterations,
  //                                 equilTolerance,
  //                                 gravVector,
  //                                 fluidWrapper,
  //                                 compFracTableWrappers,
  //                                 tempTableWrapper,
  //                                 elevationValues[0][iContact],
  //                                 pressureValues[iContact][0],
  //                                 phaseMassDens[iContact],
  //                                 contactElevation,
  //                                 contactPressure[0],
  //                                 contactPhaseMassDens[0],
  //                                 contactPhaseDens[0],
  //                                 contactPhaseCompFrac[0],
  //                                 iContact );
  //   err = LvArray::math::abs( contactPressure[0][0] - contactPressure[0][1] ) / contactPressure[0][1];
  //   std::cout << "march iter = " << marchIter << ": err = " << err << ", gasPres = " << contactPressure[0][0] << ", liqPres = " << contactPressure[0][1] << std::endl;

  //   }

  //   // Compute reference condition for points above iContact if they exist in case where iContact is above iDatum
  //   if ( iContact + 1 < size && iContact > iDatum )
  //   {
  //     returnVal =
  //       computeHydrostaticPressure( numComps,
  //                                   numPhases,
  //                                   ip_gas,
  //                                   ip_oil,
  //                                   ip_water,
  //                                   phaseContacts,
  //                                   phaseMinVolumeFraction,
  //                                   maxNumEquilIterations,
  //                                   equilTolerance,
  //                                   gravVector,
  //                                   fluidWrapper,
  //                                   compFracTableWrappers,
  //                                   tempTableWrapper,
  //                                   contactElevation,
  //                                   contactPressure[0],
  //                                   contactPhaseMassDens[0],
  //                                   elevationValues[0][iContact + 1],
  //                                   pressureValues[iContact + 1][0],
  //                                   phaseMassDens[iContact + 1],
  //                                   phaseDens[iContact + 1][0],
  //                                   phaseCompFrac[iContact + 1][0],
  //                                   iContact+1 );
  //   }
  //   // Compute reference condition for points below iContact if they exist in case where iContact is below iDatum
  //   else if ( iContact + 1 > 0 && iContact < iDatum )
  //   {
  //     returnVal =
  //       computeHydrostaticPressure( numComps,
  //                                   numPhases,
  //                                   ip_gas,
  //                                   ip_oil,
  //                                   ip_water,
  //                                   phaseContacts,
  //                                   phaseMinVolumeFraction,
  //                                   maxNumEquilIterations,
  //                                   equilTolerance,
  //                                   gravVector,
  //                                   fluidWrapper,
  //                                   compFracTableWrappers,
  //                                   tempTableWrapper,
  //                                   contactElevation,
  //                                   contactPressure[0],
  //                                   contactPhaseMassDens[0],
  //                                   elevationValues[0][iContact - 1],
  //                                   pressureValues[iContact - 1][0],
  //                                   phaseMassDens[iContact - 1],
  //                                   phaseDens[iContact - 1][0],
  //                                   phaseCompFrac[iContact - 1][0],
  //                                   iContact-1 );
  //   }

  //   return returnVal;
  // }


  template< typename FLUID_WRAPPER >
  static ReturnType
  launch( localIndex const & size,
          integer const & numComps,
          integer const & numPhases,
          integer const & ip_gas,
          integer const & ip_oil,
          integer const & ip_water,
          integer const & maxNumEquilIterations,
          arrayView1d< real64 const > const & phaseContacts,
          arrayView1d< real64 const > const & phaseMinVolumeFraction,
          real64 const equilTolerance,
          real64 const (&gravVector)[ 3 ],
          real64 const & minElevation,
          real64 const & elevationIncrement,
          real64 const & datumElevation,
          real64 const & datumPres,
          FLUID_WRAPPER fluidWrapper,
          arrayView1d< TableFunction::KernelWrapper const > compFracTableWrappers,
          TableFunction::KernelWrapper tempTableWrapper,
          arrayView1d< arrayView1d< real64 > const > elevationValues,
          // arrayView2d< real64 > const pressureValues,
          arrayView3d< real64, constitutive::multifluid::USD_PHASE > const & pressureValues,
          arrayView3d< real64, constitutive::multifluid::USD_PHASE > const & phaseDens,
          arrayView4d< real64, constitutive::multifluid::USD_PHASE_COMP > const & phaseCompFrac )
  {
    ReturnType returnVal = ReturnType::SUCCESS;

    std::cout << "ip_gas = " << ip_gas << ", ip_oil = " << ip_oil << ", ip_water = " << ip_water << std::endl;

    // Temporarily set all phase pressures at datum to input datum pressure.
    array3d< real64, constitutive::multifluid::LAYOUT_PHASE > datumPresInput( 1, 1, numPhases );
    for ( localIndex i = 0; i < numPhases; ++i)
    {
      datumPresInput[0][0][i] = datumPres;
    }

    std::cout << "Phase Contacts from HydrostaticPressureKernel: ";
    for ( const real64& val : phaseContacts )
    {
      std::cout << val << ", ";
    }
    std::cout << std::endl;
    std::cout << "numComps = " << numComps << ", numPhases = " << numPhases << std::endl;

    // Step 1: compute the phase mass densities at datum
    real64 const datumTemp = tempTableWrapper.compute( &datumElevation );
    array3d< real64, constitutive::multifluid::LAYOUT_PHASE > datumPhaseMassDens( 1, 1, numPhases );
    computeDatumPhaseMassDens( numComps,
                               numPhases,
                               ip_gas,
                               ip_water,
                               phaseMinVolumeFraction,
                               datumElevation,
                               datumPres,
                               datumTemp,
                               datumPhaseMassDens,
                               compFracTableWrappers,
                               fluidWrapper );

    // Step 2: populate the elevation array.
    forAll< parallelHostPolicy >( size, [&] ( localIndex const i )
    {
      real64 const elevation = minElevation + i * elevationIncrement;
      elevationValues[0][i] = elevation;
    } );

    // Step 3: find the closest elevation to datumElevation. Its index is denoted as iDatum
    localIndex const iDatum = LvArray::sortedArrayManipulation::find( elevationValues[0].begin(),
                                                                      elevationValues[0].size(),
                                                                      datumElevation );
    std::cout << "Nearest index to the datum = " << iDatum << std::endl;

    // Step 4: compute the mass density and pressure at iDatum elevation
    array2d< real64 > phaseMassDens( size, numPhases );
    // temporary array without permutation to compile on Lassen (TODO: what is this treatment?)
    array1d< real64 > datumPhaseMassDensTmp( numPhases );
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      datumPhaseMassDensTmp[ip] = datumPhaseMassDens[0][0][ip];
    }

    integer ip_pp = primaryPressure( numPhases,
                                     ip_gas,
                                     ip_oil,
                                     ip_water,
                                     datumElevation,
                                     phaseContacts );
    integer ip_cp = -1;
    integer ip_op = -1;

    ReturnType const iDatumReturnVal =
      computeHydrostaticPressure( numComps,
                                  numPhases,
                                  ip_gas,
                                  ip_oil,
                                  ip_water,
                                  ip_pp,
                                  phaseContacts,
                                  phaseMinVolumeFraction,
                                  maxNumEquilIterations,
                                  equilTolerance,
                                  gravVector,
                                  fluidWrapper,
                                  compFracTableWrappers,
                                  tempTableWrapper,
                                  datumElevation,
                                  datumPresInput[0][0],
                                  datumPhaseMassDensTmp,
                                  elevationValues[0][iDatum],
                                  pressureValues[iDatum][0],
                                  phaseMassDens[iDatum],
                                  phaseDens[iDatum][0],
                                  phaseCompFrac[iDatum][0],
                                  iDatum );
    if( iDatumReturnVal == ReturnType::FAILED_TO_CONVERGE )
    {
      return ReturnType::FAILED_TO_CONVERGE;
    }
    else if( iDatumReturnVal == ReturnType::DETECTED_MULTIPHASE_FLOW )
    {
      returnVal = ReturnType::DETECTED_MULTIPHASE_FLOW;
    }

    // Contacts classified as "close" and "far"
    int numContacts = phaseContacts.size();
    // Find the index of the contact that is close to the datum
    integer iContactClose = 0;
    integer iContactFar = -1;
    if ( numContacts > 1 )
    {
      // Assumes phaseContacts is of size 2
      if ( LvArray::math::abs( phaseContacts[0] - datumElevation ) <= 
           LvArray::math::abs( phaseContacts[1] - datumElevation ) )
      {
        iContactClose = 0;
        iContactFar = 1;
      }
      else
      {
        iContactClose = 1;
        iContactFar = 0;
      }
    }
    std::cout << "phaseContacts = " << phaseContacts << ", numContacts, iContactClose, iContactFar = " 
              << numContacts << ", " << iContactClose << ", " << iContactFar << std::endl;


    integer iContactCloseIndex = iDatum;
    integer iContactTop = iDatum;
    integer iContactBottom = iDatum;
    // March from datum to the close contact
    if ( LvArray::math::abs( datumElevation - phaseContacts[iContactClose] ) > 1e-12 )
    {
      // Find index of the closest element in elevationValues to the close contact, denoted as iContact
      iContactCloseIndex = LvArray::sortedArrayManipulation::find( elevationValues[0].begin(),
                                                                   elevationValues[0].size(),
                                                                   phaseContacts[iContactClose] );
      std::cout << "March from datum to closest phase contact: " << std::endl;
      std::cout << "Nearest index to the close contact = " << iContactCloseIndex << std::endl;

      evaluatePressureIndices( numPhases,
                               ip_gas,
                               ip_oil,
                               ip_water,
                               datumElevation,
                               phaseContacts[iContactClose],
                               phaseContacts,
                               ip_pp,
                               ip_cp,
                               ip_op );
      std::cout << "ip_pp = " << ip_pp << ", ip_cp = " << ip_cp << ", ip_op = " << ip_op << std::endl;

      // Step: for each elevation between gwc and datum
      returnVal = marchBetweenTwoElevations( size,
                                             datumElevation,
                                             phaseContacts[iContactClose],
                                             numComps,
                                             numPhases,
                                             ip_gas,
                                             ip_oil,
                                             ip_water,
                                             ip_pp,
                                             ip_cp,
                                             ip_op,
                                             phaseContacts,
                                             phaseMinVolumeFraction,
                                             maxNumEquilIterations,
                                             equilTolerance,
                                             gravVector,
                                             fluidWrapper,
                                             compFracTableWrappers,
                                             tempTableWrapper,
                                             elevationValues,
                                             pressureValues,
                                             phaseMassDens,
                                             phaseDens,
                                             phaseCompFrac );
    }

    integer iContactFarIndex = iContactCloseIndex;
    // March from close contact to far contact
    if ( LvArray::math::abs( phaseContacts[iContactClose] - phaseContacts[iContactFar] ) > 1e-12 )
    {
      // Find index of the closest element in elevationValues to the far contact, denoted as iContact
      iContactFarIndex = LvArray::sortedArrayManipulation::find( elevationValues[0].begin(),
                                                                 elevationValues[0].size(),
                                                                 phaseContacts[iContactFar] );
      std::cout << "March from the closest contact to the farthest contact: " << std::endl;
      std::cout << "Nearest index to the far contact = " << iContactFarIndex << std::endl;

      evaluatePressureIndices( numPhases,
                               ip_gas,
                               ip_oil,
                               ip_water,
                               phaseContacts[iContactClose],
                               phaseContacts[iContactFar],
                               phaseContacts,
                               ip_pp,
                               ip_cp,
                               ip_op );
      std::cout << "ip_pp = " << ip_pp << ", ip_cp = " << ip_cp << ", ip_op = " << ip_op << std::endl;

      // Step: for each elevation between gwc and datum
      returnVal = marchBetweenTwoElevations( size,
                                             phaseContacts[iContactClose],
                                             phaseContacts[iContactFar],
                                             numComps,
                                             numPhases,
                                             ip_gas,
                                             ip_oil,
                                             ip_water,
                                             ip_pp,
                                             ip_cp,
                                             ip_op,
                                             phaseContacts,
                                             phaseMinVolumeFraction,
                                             maxNumEquilIterations,
                                             equilTolerance,
                                             gravVector,
                                             fluidWrapper,
                                             compFracTableWrappers,
                                             tempTableWrapper,
                                             elevationValues,
                                             pressureValues,
                                             phaseMassDens,
                                             phaseDens,
                                             phaseCompFrac );
    }

    
    if (phaseContacts[iContactClose] > phaseContacts[iContactFar])
    {
      iContactTop = iContactCloseIndex;
      iContactBottom = iContactFarIndex;
    }
    else
    {
      iContactTop = iContactFarIndex;
      iContactBottom = iContactCloseIndex;
    }


    // Step 4: for each elevation above the max( top of in-between region, datumElevation ) 
    // Assumes phaseContacts = [ goc, owc ] and goc >= owc
    // localIndex iForUpperElevations = ( iDatum >= iContactTop ) ? iDatum : iContactTop + 1;
    localIndex iForUpperElevations = iContactTop;
    integer ip_pp_above = primaryPressure( numPhases,
                                           ip_gas,
                                           ip_oil,
                                           ip_water,
                                          //  ( datumElevation > elevationValues[0][iContactTop] ) ? datumElevation + 1e-2 : elevationValues[0][iContactTop] + 1e-2,
                                           elevationValues[0][iForUpperElevations],
                                           phaseContacts );
    std::cout << "ip_pp_above = " << ip_pp_above << ", iContactTop = " << iContactTop << ", elevationValues[0][iContactTop] = " << elevationValues[0][iContactTop] << std::endl;
    returnVal = computeHydrostaticPressureAtMultipleElevations( iForUpperElevations,
                                                                size - 1,
                                                                numComps,
                                                                numPhases,
                                                                ip_gas,
                                                                ip_oil,
                                                                ip_water,
                                                                ip_pp_above,
                                                                phaseContacts,
                                                                phaseMinVolumeFraction,
                                                                maxNumEquilIterations,
                                                                equilTolerance,
                                                                gravVector,
                                                                fluidWrapper,
                                                                compFracTableWrappers,
                                                                tempTableWrapper,
                                                                elevationValues,
                                                                pressureValues,
                                                                phaseMassDens,
                                                                phaseDens,
                                                                phaseCompFrac );

    // Step 5: for each elevation below the in-between region
    // localIndex iForLowerElevations = ( iDatum <= iContactBottom ) ? iDatum : iContactBottom - 1;
    localIndex iForLowerElevations = iContactBottom;
    integer ip_pp_below = primaryPressure( numPhases,
                                           ip_gas,
                                           ip_oil,
                                           ip_water,
                                          //  ( datumElevation < elevationValues[0][iContactBottom] ) ? datumElevation - 1.0 : elevationValues[0][iContactBottom] - 1.0,
                                          elevationValues[0][iForLowerElevations],
                                           phaseContacts );
    std::cout << "ip_pp_below = " << ip_pp_below << std::endl;
    returnVal = computeHydrostaticPressureAtMultipleElevations( iForLowerElevations,
                                                                0,
                                                                numComps,
                                                                numPhases,
                                                                ip_gas,
                                                                ip_oil,
                                                                ip_water,
                                                                ip_pp_below,
                                                                phaseContacts,
                                                                phaseMinVolumeFraction,
                                                                maxNumEquilIterations,
                                                                equilTolerance,
                                                                gravVector,
                                                                fluidWrapper,
                                                                compFracTableWrappers,
                                                                tempTableWrapper,
                                                                elevationValues,
                                                                pressureValues,
                                                                phaseMassDens,
                                                                phaseDens,
                                                                phaseCompFrac );

    const std::string & f0 = "hydro_pres.txt";
    std::ofstream outFile0(f0);
    for ( localIndex i = 0; i < size; ++i )
    {
      outFile0 << elevationValues[0][i] << "\t";
      for ( localIndex ip = 0; ip < numPhases; ++ip )
      {
        outFile0 << pressureValues[i][0][ip] << "\t";
      }
      for ( localIndex ip = 0; ip < numPhases; ++ip )
      {
        outFile0 << phaseDens[i][0][ip] << "\t";
      }
      for ( localIndex ip = 0; ip < numPhases; ++ip )
      {
        for ( localIndex ic = 0; ic < numComps; ++ic )
        {
          outFile0 << phaseCompFrac[i][0][ip][ic] << "\t";
        }
      }
      outFile0 << std::endl;
    }
    outFile0.close();

    return returnVal;
  }

  // template< typename FLUID_WRAPPER >
  // static ReturnType
  // compute( integer const numComps,
  //          integer const numPhases,
  //          integer const & ip_gas,
  //          integer const & ip_oil,
  //          integer const & ip_water,
  //          integer const maxNumEquilIterations,
  //          real64 const & equilTolerance,
  //          real64 const (&gravVector)[ 3 ],
  //          arrayView1d< real64 const > phaseContacts,
  //          FLUID_WRAPPER fluidWrapper,
  //          arrayView1d< TableFunction::KernelWrapper const > compFracTableWrappers,
  //          TableFunction::KernelWrapper tempTableWrapper,
  //          real64 const & refElevation,
  //          arraySlice1d< real64 const > const & refPres,
  //          arraySlice1d< real64 const > const & refPhaseMassDens,
  //          real64 const & newElevation,
  //          arraySlice1d< real64 > const & newPres,
  //          arraySlice1d< real64 > const & newPhaseMassDens,
  //          arraySlice1d< real64 > const & newPhaseDens,
  //          arraySlice2d< real64 > const & newPhaseCompFrac,
  //          localIndex const index )
  // {
  //   integer ip_pres = -1;
  //   // Choose the primary phase
  //   if ( ip_gas >= 0 && ip_oil >= 0 && ip_water >= 0 )
  //   {
  //     // Three phase case ( assume phaseContacts[0] = goc, phaseContacts[1] = owc )    
  //     ip_pres = ( elevation <= phaseContacts[1] ) ? ip_water :
  //               ( elevation <= phaseContacts[0] ) ? ip_oil : ip_gas;
  //   }
  //   else if ( ip_oil >= 0 && ip_water >= 0 )
  //   {
  //     // phases = oil + water ( assume phaseContacts[0] = owc )
  //     ip_pres = ( elevation <= phaseContacts[0] ) ? ip_water : ip_oil;
  //   }
  //   else if ( ip_gas >= 0 && ip_water >= 0 )
  //   {
  //     // phases = gas + water ( assume phaseContacts[0] = gwc )
  //     ip_pres = ( elevation <= phaseContacts[0] ) ? ip_water : ip_gas;
  //   }
  //   else if ( ip_gas >= 0 && ip_oil >= 0 )
  //   {
  //     // phases = gas + oil ( assume phaseContacts[0] = goc )
  //     ip_pres = ( elevation <= phaseContacts[0] ) ? ip_oil : ip_gas;
  //   }
  //   else
  //   {
  //     return ReturnType::DETECTED_SINGLEPHASE_FLOW;
  //   }
  //   std::cout << "ip_pres = " << ip_pres << std::endl;

  //   real64 const temp = tempTableWrapper.compute( &newElevation );
  //   StackArray< real64, 2, constitutive::MultiFluidBase::MAX_NUM_COMPONENTS, compflow::LAYOUT_COMP > compFrac( 1, numComps );
  //   for( integer ic = 0; ic < numComps; ++ic )
  //   {
  //     compFrac[0][ic] = compFracTableWrappers[ic].compute( &newElevation );
  //   }

  // }

  template< typename FLUID_WRAPPER >
  static ReturnType
  phaseCorrection( integer const & numComps,
                   integer const & numPhases,
                   integer const & ip_gas,
                   integer const & ip_water,
                   arrayView1d< real64 const > const & phaseMinVolumeFraction,
                   real64 const & pres,
                   real64 const & temp,
                   arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & uncorrCompFrac,
                   arraySlice1d< real64 const, constitutive::multifluid::USD_PHASE - 2 > const & phaseFrac,
                   arraySlice1d< real64, compflow::USD_COMP - 1 > const & compFrac,
                   FLUID_WRAPPER fluidWrapper )
  {
    // TODO: should the if conditions be independent or else ifs after the first if
    std::cout << "Inside phaseCorrection, phaseFrac = " << phaseFrac << std::endl;
    StackArray< real64, 2, constitutive::MultiFluidBase::MAX_NUM_COMPONENTS, compflow::LAYOUT_COMP > addedCompFrac( 1, numComps );
    integer ip_phase = -1;
    integer ip_otherPhase = -1;
    bool phaseCorrectionNeeded = false;
    if ( ip_gas >= 0 && LvArray::math::abs( phaseFrac[ip_gas] ) < MIN_FOR_PHASE_PRESENCE
                     && LvArray::math::abs( phaseMinVolumeFraction[ip_gas] ) > MIN_FOR_PHASE_PRESENCE )
    {
      addedCompFrac[0][0] = 1.0; // hard-coded (assumes co2 is at index 0)
      ip_phase = ip_gas;
      ip_otherPhase = ip_water;
      phaseCorrectionNeeded = true;
    }
    // if ( ip_oil >= 0 && LvArray::math::abs( phaseFrac[ip_oil] ) < MIN_FOR_PHASE_PRESENCE
    //                       && LvArray::math::abs( phaseMinVolumeFraction[ip_oil] ) > MIN_FOR_PHASE_PRESENCE )
    // {
    //   addedCompFrac[0][0] = 1.0; // hard-coded (assumes some hydrocarbon is at index 0)
    //   ip_phase = ip_oil;
    //   phaseCorrectionNeeded = true;
    // }

    if ( ip_water >= 0 && LvArray::math::abs( phaseFrac[ip_water] ) < MIN_FOR_PHASE_PRESENCE
                     && LvArray::math::abs( phaseMinVolumeFraction[ip_water] ) > MIN_FOR_PHASE_PRESENCE )
    {
      addedCompFrac[0][1] = 1.0; // hard-coded (assumes water is at index 1)
      ip_phase = ip_water;
      ip_otherPhase = ip_gas;
      phaseCorrectionNeeded = true;
    }

    if ( phaseCorrectionNeeded )
    {
      return applyPhaseCorrection( numComps,
                                   numPhases,
                                   ip_phase,
                                   ip_otherPhase,
                                   pres,
                                   temp,
                                   uncorrCompFrac,
                                   addedCompFrac[0],
                                   compFrac,
                                   fluidWrapper );
    }
    else
    {
      for ( localIndex ic = 0; ic < numComps; ++ic )
      {
        compFrac[ic] = uncorrCompFrac[ic];
      }
      std::cout << "Phase correction not needed: phaseFrac = " << phaseFrac << std::endl;
      return ReturnType::PHASE_CORRECTION_NOT_NEEDED;
    }

  }

  template< typename FLUID_WRAPPER >
  static ReturnType
  applyPhaseCorrection( integer const & numComps,
                        integer const & numPhases,
                        integer const & ip_phase,
                        integer const & ip_otherPhase,
                        real64 const & pres,
                        real64 const & temp,
                        arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & uncorrCompFrac,
                        arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & addedCompFrac,
                        arraySlice1d< real64, compflow::USD_COMP - 1 > const & compFrac,
                        FLUID_WRAPPER fluidWrapper )
  {
    // current implementation only works for when the added phase is gas or water
    // flash inputs
    array3d< real64, constitutive::multifluid::LAYOUT_PHASE > phaseFrac( 1, 1, numPhases );
    array3d< real64, constitutive::multifluid::LAYOUT_PHASE > phaseDens( 1, 1, numPhases );
    array3d< real64, constitutive::multifluid::LAYOUT_PHASE > phaseMassDens( 1, 1, numPhases );
    array3d< real64, constitutive::multifluid::LAYOUT_PHASE > phaseVisc( 1, 1, numPhases );
    array3d< real64, constitutive::multifluid::LAYOUT_PHASE > phaseEnthalpy( 1, 1, numPhases );
    array3d< real64, constitutive::multifluid::LAYOUT_PHASE > phaseInternalEnergy( 1, 1, numPhases );
    array4d< real64, constitutive::multifluid::LAYOUT_PHASE_COMP > phaseCompFrac( 1, 1, numPhases, numComps );
    real64 totalDens = 0.0;

    real64 a_low = 0.0;
    real64 a_high = 1.0;
    real64 a = 0.0;
    int maxIters = 100;
    real64 targetPhaseFrac = 1e-4;
    real64 err = 1e10;
    for ( int iter = 0; iter < maxIters; ++iter )
    {
      a = ( a_high + a_low ) * 0.5;
      mixingStep( numComps,
                  a,
                  uncorrCompFrac,
                  addedCompFrac,
                  compFrac );
      constitutive::MultiFluidBase::KernelWrapper::computeValues( fluidWrapper,
                                                                  pres,
                                                                  temp,
                                                                  compFrac,
                                                                  phaseFrac[0][0],
                                                                  phaseDens[0][0],
                                                                  phaseMassDens[0][0],
                                                                  phaseVisc[0][0],
                                                                  phaseEnthalpy[0][0],
                                                                  phaseInternalEnergy[0][0],
                                                                  phaseCompFrac[0][0],
                                                                  totalDens );
      err = ( phaseFrac[0][0][ip_phase] - targetPhaseFrac ) / targetPhaseFrac;
      // std::cout << "iter = " << iter << ", err = " << err << ", phase frac = " 
      //           << phaseFrac[0][0][ip_phase] << ", dens = " << phaseDens[0][0] 
      //           << ", phaseComp = " << phaseCompFrac[0][0] << std::endl;
      if ( LvArray::math::abs( phaseFrac[0][0][ip_otherPhase] - 1.0 ) < MIN_FOR_PHASE_PRESENCE )
      {
        // std::cout << "water phase encountered, phase fraction = " << phaseFrac[0][0][ip_otherPhase] << std::endl;
        a_low = a;
      }
      else if ( LvArray::math::abs( phaseFrac[0][0][ip_phase] - 1.0 ) < MIN_FOR_PHASE_PRESENCE )
      {
        // std::cout << "gas phase encountered, phase fraction = " << phaseFrac[0][0][ip_phase] << std::endl;
        a_high = a;
      }
      else
      {
        if ( err > 0 )
        {
          a_high = a;
        }
        else
        {
          a_low = a;
        }
        if ( LvArray::math::abs( err ) < 1e-5 )
        {
          break;
        }
      }
    }
    std::cout << "err = " << err << ", uncorrected compFrac = " << uncorrCompFrac 
              << ", new compFrac = " << compFrac << std::endl;

    return ReturnType::SUCCESS;
  }

  static void mixingStep( integer const & numComps,
                          real64 const & a,
                          arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & uncorrCompFrac,
                          arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & addedCompFrac,
                          arraySlice1d< real64, compflow::USD_COMP - 1 > const & compFrac )
  {
    real64 tot = 0.0;
    for ( localIndex ic = 0; ic < numComps; ++ic )
    {
       compFrac[ic] = a * addedCompFrac[ic] + ( 1 - a ) * uncorrCompFrac[ic];
       if ( compFrac[ic] < MIN_FOR_COMP_PRESENCE )
       {
         compFrac[ic] = 0.0; 
       }
       tot += compFrac[ic];
    }
    for ( localIndex ic = 0; ic < numComps; ++ic )
    {
      compFrac[ic] /= tot;
    }
  }

  static integer primaryPressure( integer const & numPhases,
                                  integer const & ip_gas,
                                  integer const & ip_oil,
                                  integer const & ip_water,
                                  real64 const & elevation,
                                  arrayView1d< real64 const > const & phaseContacts )
  {
    integer ip_pp = -1;
    // Choose the primary phase
    if ( numPhases == 1 )
    {
      ip_pp = ( ip_gas >= 0 ) ? ip_gas :
                ( ip_oil >= 0 ) ? ip_oil :
                ( ip_water >= 0 ) ? ip_water : 0;
    }
    else if ( ip_gas >= 0 && ip_oil >= 0 && ip_water >= 0 )
    {
      // Three phase case ( assume phaseContacts[0] = goc, phaseContacts[1] = owc )    
      ip_pp = ( elevation <= phaseContacts[1] ) ? ip_water :
                ( elevation <= phaseContacts[0] ) ? ip_oil : ip_gas;
      std::cout << "elevation = " << elevation << ", ip_pp = " << ip_pp << std::endl;
    }
    else if ( ip_oil >= 0 && ip_water >= 0 )
    {
      // phases = oil + water ( assume phaseContacts[0] = owc )
      ip_pp = ( elevation <= phaseContacts[0] ) ? ip_water : ip_oil;
    }
    else if ( ip_gas >= 0 && ip_water >= 0 )
    {
      // phases = gas + water ( assume phaseContacts[0] = gwc )
      ip_pp = ( elevation <= phaseContacts[0] ) ? ip_water : ip_gas;
    }
    else if ( ip_gas >= 0 && ip_oil >= 0 )
    {
      // phases = gas + oil ( assume phaseContacts[0] = goc )
      ip_pp = ( elevation <= phaseContacts[0] ) ? ip_oil : ip_gas;
    }
    // TODO: deal with the case if ip_pp is still -1
    return ip_pp;
  }

  static void evaluatePressureIndices( integer const & numPhases,
                                       integer const & ip_gas,
                                       integer const & ip_oil,
                                       integer const & ip_water,
                                       real64 const & startElevation,
                                       real64 const & endElevation,
                                       arrayView1d< real64 const > const & phaseContacts,
                                       integer & ip_pp,
                                       integer & ip_cp,
                                       integer & ip_op )
  {
    ip_pp = -1;
    ip_cp = -1;
    ip_op = -1;
    if ( numPhases > 1 )
    {
      if ( ip_gas >= 0 && ip_oil >= 0 && ip_water >= 0 )
      {
        real64 const goc = phaseContacts[0];
        real64 const owc = phaseContacts[1];
        if ( startElevation < owc && LvArray::math::abs( endElevation - owc ) < 1e-12 )
        {
          // we are below owc and will march towards it
          ip_pp = ip_water;
          ip_cp = ip_oil;
          ip_op = ip_gas;
        }
        else if ( startElevation > owc && LvArray::math::abs( endElevation - owc ) < 1e-12 )
        {
          // we are above owc and will march towards it
          ip_pp = ip_oil;
          ip_cp = ip_water;
          ip_op = ip_gas;
        }
        else if ( startElevation < goc && LvArray::math::abs( endElevation - goc ) < 1e-12 )
        {
          // we are between owc and goc and will march towards goc
          ip_pp = ip_oil;
          ip_cp = ip_gas;
          ip_op = ip_water;
        }
        else if ( startElevation > goc && LvArray::math::abs( endElevation - goc ) < 1e-12 )
        {
          // we are above owc and intend to march towards it
          ip_pp = ip_gas;
          ip_cp = ip_oil;
          ip_op = ip_water;
        }
      }
      else if ( ip_oil >= 0 && ip_water >= 0 )
      {
        real64 const owc = phaseContacts[0];
        if ( startElevation < owc )
        {
          // we below below the owc and intend to march towards it
          ip_pp = ip_water;
          ip_cp = ip_oil;
        }
        else
        {
          // we below aobve the owc and intend to march towards it
          ip_pp = ip_oil;
          ip_cp = ip_water; 
        }
      }
      else if ( ip_gas >= 0 && ip_water >= 0 )
      {
        real64 const gwc = phaseContacts[0];
        if ( startElevation < gwc )
        {
          // we below below the gwc and intend to march towards it
          ip_pp = ip_water;
          ip_cp = ip_gas;
        }
        else
        {
          // we below above the gwc and intend to march towards it
          ip_pp = ip_gas;
          ip_cp = ip_water; 
        }
      }
      else if ( ip_gas >= 0 && ip_oil >= 0 )
      {
        real64 const goc = phaseContacts[0];
        if ( startElevation < goc )
        {
          // we below below the goc and intend to march towards it
          ip_pp = ip_oil;
          ip_cp = ip_gas;
        }
        else
        {
          // we below above the goc and intend to march towards it
          ip_pp = ip_gas;
          ip_cp = ip_oil; 
        }
      }
    }
  }

  template< typename FLUID_WRAPPER >
  static void
  computeDatumPhaseMassDens( integer const & numComps,
                             integer const & numPhases,
                             integer const & ip_gas,
                             integer const & ip_water,
                             arrayView1d< real64 const > const & phaseMinVolumeFraction,
                             real64 const & datumElevation,
                             real64 const & datumPres,
                             real64 const & datumTemp,
                             arrayView3d< real64, constitutive::multifluid::USD_PHASE > const & datumPhaseMassDens,
                             arrayView1d< TableFunction::KernelWrapper const > compFracTableWrappers,
                             FLUID_WRAPPER fluidWrapper )
  {
    // datum fluid properties
    array2d< real64, compflow::LAYOUT_COMP > datumUncorrCompFrac( 1, numComps );
    array2d< real64, compflow::LAYOUT_COMP > datumCompFrac( 1, numComps );
    array3d< real64, constitutive::multifluid::LAYOUT_PHASE > datumPhaseFrac( 1, 1, numPhases );
    array3d< real64, constitutive::multifluid::LAYOUT_PHASE > datumPhaseDens( 1, 1, numPhases );
    array3d< real64, constitutive::multifluid::LAYOUT_PHASE > datumPhaseVisc( 1, 1, numPhases );
    array3d< real64, constitutive::multifluid::LAYOUT_PHASE > datumPhaseEnthalpy( 1, 1, numPhases );
    array3d< real64, constitutive::multifluid::LAYOUT_PHASE > datumPhaseInternalEnergy( 1, 1, numPhases );
    array4d< real64, constitutive::multifluid::LAYOUT_PHASE_COMP > datumPhaseCompFrac( 1, 1, numPhases, numComps );
    real64 datumTotalDens = 0.0;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      datumUncorrCompFrac[0][ic] = compFracTableWrappers[ic].compute( &datumElevation );
    }
    constitutive::MultiFluidBase::KernelWrapper::computeValues( fluidWrapper,
                                                                datumPres,
                                                                datumTemp,
                                                                datumUncorrCompFrac[0],
                                                                datumPhaseFrac[0][0],
                                                                datumPhaseDens[0][0],
                                                                datumPhaseMassDens[0][0],
                                                                datumPhaseVisc[0][0],
                                                                datumPhaseEnthalpy[0][0],
                                                                datumPhaseInternalEnergy[0][0],
                                                                datumPhaseCompFrac[0][0],
                                                                datumTotalDens );
    ReturnType datumPhaseCorr = phaseCorrection( numComps,
                                                 numPhases,
                                                 ip_gas,
                                                 ip_water,
                                                 phaseMinVolumeFraction,
                                                 datumPres,
                                                 datumTemp,
                                                 datumUncorrCompFrac[0],
                                                 datumPhaseFrac[0][0],
                                                 datumCompFrac[0],
                                                 fluidWrapper );
    if ( datumPhaseCorr == ReturnType::SUCCESS )
      std::cout << "Datum phase correction function completed!" << std::endl;
    constitutive::MultiFluidBase::KernelWrapper::computeValues( fluidWrapper,
                                                                datumPres,
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
    std::cout << "Datum flash results at p = " << datumPres << ", T = " << datumTemp
    << ", datum overall compositions:" << datumCompFrac[0] << std::endl;
    std::cout << "phase dens = " << datumPhaseDens[0][0] << ", phase mass dens = " <<
    datumPhaseMassDens[0][0] << ", Phase composition = " << datumPhaseCompFrac[0][0]
    << std::endl;                                     
  }

};

} // namespace isothermalCompositionalMultiphaseBaseKernels

} // namespace geos


#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONAL_HYDROSTATICPRESSUREKERNEL_HPP
