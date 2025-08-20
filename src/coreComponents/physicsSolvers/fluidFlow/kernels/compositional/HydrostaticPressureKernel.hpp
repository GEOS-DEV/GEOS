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
    integer const ip_pres = primaryPressure( numPhases,
                                             ip_gas,
                                             ip_oil,
                                             ip_water,
                                             newElevation,
                                             phaseContacts );

    // Step 3: compute the mass density at this elevation using the guess, and update pressure

    constitutive::MultiFluidBase::KernelWrapper::computeValues( fluidWrapper,
                                                                pres0[ip_pres], // flash pressure
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
                                            pres0[ip_pres],
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
                                                                  pres0[ip_pres], // flash pressure
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
                << std::endl;
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
  marchBetweenContactAndDatum( localIndex const & size,
                               real64 const & contactElevation,
                               localIndex const & iRef,
                               localIndex const & iContact,
                               integer const & numComps,
                               integer const & numPhases,
                               integer const & ip_gas,
                               integer const & ip_oil,
                               integer const & ip_water,
                               arrayView1d< real64 const > const & phaseContacts,
                               arrayView1d< real64 const > const & phaseMinVolumeFraction,
                               integer const maxNumEquilIterations,
                               real64 const & equilTolerance,
                               real64 const (&gravVector)[ 3 ],
                               FLUID_WRAPPER fluidWrapper,
                               arrayView1d< TableFunction::KernelWrapper const > compFracTableWrappers,
                               TableFunction::KernelWrapper tempTableWrapper,
                               arrayView1d< arrayView1d< real64 > const > elevationValues,
                              //  arrayView2d< real64 > const pressureValues,
                               arrayView3d< real64, constitutive::multifluid::USD_PHASE > const & pressureValues,
                               arrayView2d< real64 > const & phaseMassDens,
                               arrayView3d< real64, constitutive::multifluid::USD_PHASE > const & phaseDens,
                               arrayView4d< real64, constitutive::multifluid::USD_PHASE_COMP > const & phaseCompFrac )
  {
    // March from iRef to iContact
    ReturnType returnVal = 
      computeHydrostaticPressureAtMultipleElevations( iRef,
                                                      iContact,
                                                      numComps,
                                                      numPhases,
                                                      ip_gas,
                                                      ip_oil,
                                                      ip_water,
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
    // Compute phase presssures and densities at the contact using iContact as the reference
    array2d< real64 > contactPressure( 1, numPhases );
    array2d< real64 > contactPhaseMassDens( 1, numPhases );
    array2d< real64 > contactPhaseDens( 1, numPhases );
    array3d< real64 > contactPhaseCompFrac( 1, numPhases, numComps );
    returnVal =
      computeHydrostaticPressure( numComps,
                                  numPhases,
                                  ip_gas,
                                  ip_oil,
                                  ip_water,
                                  phaseContacts,
                                  phaseMinVolumeFraction,
                                  maxNumEquilIterations,
                                  equilTolerance,
                                  gravVector,
                                  fluidWrapper,
                                  compFracTableWrappers,
                                  tempTableWrapper,
                                  elevationValues[0][iContact],
                                  pressureValues[iContact][0],
                                  phaseMassDens[iContact],
                                  contactElevation,
                                  contactPressure[0],
                                  contactPhaseMassDens[0],
                                  contactPhaseDens[0],
                                  contactPhaseCompFrac[0],
                                  iContact );
    // Compute relative error defined as the relative difference between the phase pressures at the contact                              
    real64 err = LvArray::math::abs( contactPressure[0][0] - contactPressure[0][1] ) / contactPressure[0][1];
    std::cout << "march iter = 1: err = " << err << ", gasPres = " << contactPressure[0][0] << ", liqPres = " << contactPressure[0][1] << std::endl;
    int maxMarchIterations = 10;

    // Marching Loop
    for ( int marchIter = 1; marchIter < maxMarchIterations; ++marchIter )
    {
      if ( err < 1e-5 ) // maybe use equilTolerance here as well and modify err from relative to absolute
      {
        std::cout << "March converged, march iter = " << marchIter << ": gasPres = " 
        << contactPressure[0][0] << ", liqPres = " << contactPressure[0][1] << std::endl;
        break;
      }
      // equate the phase pressure at the contact
      contactPressure[0][0] = contactPressure[0][1];
      // array2d< real64 > tmpRefPressure( 1, numPhases );
      // array2d< real64 > tmpRefPhaseMassDens( 1, numPhases );
      // array3d< real64 > tmpRefPhaseCompFrac( 1, numPhases, numComps );
      real64 const saveRefPhasePressure = pressureValues[iRef][0][1]; // saves the known phase pressure at iRef
      // Estimate the unknown phase pressure and density at iRef using the updated contact phase pressures
      returnVal =
        computeHydrostaticPressure( numComps,
                                    numPhases,
                                    ip_gas,
                                    ip_oil,
                                    ip_water,
                                    phaseContacts,
                                    phaseMinVolumeFraction,
                                    maxNumEquilIterations,
                                    equilTolerance,
                                    gravVector,
                                    fluidWrapper,
                                    compFracTableWrappers,
                                    tempTableWrapper,
                                    contactElevation,
                                    contactPressure[0],
                                    contactPhaseMassDens[0],
                                    elevationValues[0][iRef],
                                    pressureValues[iRef][0],
                                    phaseMassDens[iRef],
                                    phaseDens[iRef][0],
                                    phaseCompFrac[iRef][0],
                                    iRef );
      pressureValues[iRef][0][1] = saveRefPhasePressure;
      // March from iRef to iContact 
      returnVal = 
      computeHydrostaticPressureAtMultipleElevations( iRef,
                                                      iContact,
                                                      numComps,
                                                      numPhases,
                                                      ip_gas,
                                                      ip_oil,
                                                      ip_water,
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
    // Compute phase presssures and densities at the contact using iContact as the reference
    returnVal =
      computeHydrostaticPressure( numComps,
                                  numPhases,
                                  ip_gas,
                                  ip_oil,
                                  ip_water,
                                  phaseContacts,
                                  phaseMinVolumeFraction,
                                  maxNumEquilIterations,
                                  equilTolerance,
                                  gravVector,
                                  fluidWrapper,
                                  compFracTableWrappers,
                                  tempTableWrapper,
                                  elevationValues[0][iContact],
                                  pressureValues[iContact][0],
                                  phaseMassDens[iContact],
                                  contactElevation,
                                  contactPressure[0],
                                  contactPhaseMassDens[0],
                                  contactPhaseDens[0],
                                  contactPhaseCompFrac[0],
                                  iContact );
    err = LvArray::math::abs( contactPressure[0][0] - contactPressure[0][1] ) / contactPressure[0][1];
    std::cout << "march iter = " << marchIter << ": err = " << err << ", gasPres = " << contactPressure[0][0] << ", liqPres = " << contactPressure[0][1] << std::endl;

    }

    // Compute reference condition for points above iContact if they exist in case where iContact is above iRef
    if ( iContact + 1 < size && iContact > iRef )
    {
      returnVal =
        computeHydrostaticPressure( numComps,
                                    numPhases,
                                    ip_gas,
                                    ip_oil,
                                    ip_water,
                                    phaseContacts,
                                    phaseMinVolumeFraction,
                                    maxNumEquilIterations,
                                    equilTolerance,
                                    gravVector,
                                    fluidWrapper,
                                    compFracTableWrappers,
                                    tempTableWrapper,
                                    contactElevation,
                                    contactPressure[0],
                                    contactPhaseMassDens[0],
                                    elevationValues[0][iContact + 1],
                                    pressureValues[iContact + 1][0],
                                    phaseMassDens[iContact + 1],
                                    phaseDens[iContact + 1][0],
                                    phaseCompFrac[iContact + 1][0],
                                    iContact+1 );
    }
    // Compute reference condition for points below iContact if they exist in case where iContact is below iRef
    else if ( iContact + 1 > 0 && iContact < iRef )
    {
      returnVal =
        computeHydrostaticPressure( numComps,
                                    numPhases,
                                    ip_gas,
                                    ip_oil,
                                    ip_water,
                                    phaseContacts,
                                    phaseMinVolumeFraction,
                                    maxNumEquilIterations,
                                    equilTolerance,
                                    gravVector,
                                    fluidWrapper,
                                    compFracTableWrappers,
                                    tempTableWrapper,
                                    contactElevation,
                                    contactPressure[0],
                                    contactPhaseMassDens[0],
                                    elevationValues[0][iContact - 1],
                                    pressureValues[iContact - 1][0],
                                    phaseMassDens[iContact - 1],
                                    phaseDens[iContact - 1][0],
                                    phaseCompFrac[iContact - 1][0],
                                    iContact-1 );
    }

    return returnVal;
  }


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
    std::cout << "ip_gas = " << ip_gas << std::endl;
    std::cout << "ip_water = " << ip_water << std::endl;
    std::cout << "ip_oil = " << ip_oil << std::endl;

    ReturnType returnVal = ReturnType::SUCCESS;
    // bool march_flag = false;

    // assume for now that contact = datum
    // array1d< real64 > datumPresInput( numPhases );
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

    // Step 2: find the closest elevation to datumElevation

    forAll< parallelHostPolicy >( size, [&] ( localIndex const i )
    {
      real64 const elevation = minElevation + i * elevationIncrement;
      elevationValues[0][i] = elevation;
    } );

    // localIndex const uniqueNumPoints = LvArray::sortedArrayManipulation::makeSortedUnique(elevationValues[0].begin(), elevationValues[0].end()); // sort in ascending order
    // elevationValues[0].resize( uniqueNumPoints );
    // pressureValues.resize( uniqueNumPoints );



    localIndex const iRef = LvArray::sortedArrayManipulation::find( elevationValues[0].begin(),
                                                                    elevationValues[0].size(),
                                                                    datumElevation );
    std::cout << "Nearest index to the datum = " << iRef << std::endl;

    // Step 3: compute the mass density and pressure at the reference elevation

    array2d< real64 > phaseMassDens( size, numPhases );
    // temporary array without permutation to compile on Lassen
    array1d< real64 > datumPhaseMassDensTmp( numPhases );
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      datumPhaseMassDensTmp[ip] = datumPhaseMassDens[0][0][ip];
    }

    ReturnType const refReturnVal =
      computeHydrostaticPressure( numComps,
                                  numPhases,
                                  ip_gas,
                                  ip_oil,
                                  ip_water,
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
                                  elevationValues[0][iRef],
                                  pressureValues[iRef][0],
                                  phaseMassDens[iRef],
                                  phaseDens[iRef][0],
                                  phaseCompFrac[iRef][0],
                                  iRef );
    if( refReturnVal == ReturnType::FAILED_TO_CONVERGE )
    {
      return ReturnType::FAILED_TO_CONVERGE;
    }
    else if( refReturnVal == ReturnType::DETECTED_MULTIPHASE_FLOW )
    {
      returnVal = ReturnType::DETECTED_MULTIPHASE_FLOW;
    }

    localIndex iContact = iRef;
    if ( LvArray::math::abs( datumElevation - phaseContacts[0] ) > 1e-12 )
    {
      // Find index of the closed element in elevationValues to the contact
      iContact = LvArray::sortedArrayManipulation::find( elevationValues[0].begin(),
                                                                          elevationValues[0].size(),
                                                                          phaseContacts[0] );
      std::cout << "Nearest index to the contact = " << iContact << std::endl;
      // Step: for each elevation between goc and datum
      returnVal = marchBetweenContactAndDatum( size,
                                               phaseContacts[0],
                                               iRef,
                                               iContact,
                                               numComps,
                                               numPhases,
                                               ip_gas,
                                               ip_oil,
                                               ip_water,
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

    // Step 4: for each elevation above the in-between region
    localIndex iForUpperElevations = ( iRef >= iContact ) ? iRef : iContact + 1;
    returnVal = computeHydrostaticPressureAtMultipleElevations( iForUpperElevations,
                                                                size - 1,
                                                                numComps,
                                                                numPhases,
                                                                ip_gas,
                                                                ip_oil,
                                                                ip_water,
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
    localIndex iForLowerElevations = ( iRef <= iContact ) ? iRef : iContact - 1;
    returnVal = computeHydrostaticPressureAtMultipleElevations( iForLowerElevations,
                                                                0,
                                                                numComps,
                                                                numPhases,
                                                                ip_gas,
                                                                ip_oil,
                                                                ip_water,
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

    // const std::string & f0 = "hydro_pres.txt";
    // std::ofstream outFile0(f0);
    // for ( localIndex i = 0; i < size; ++i )
    // {
    //   outFile0 << elevationValues[0][i] << "\t" << pressureValues[i][0] << "\t" << pressureValues[i][1] << "\t";
    //   for ( localIndex ip = 0; ip < numPhases; ++ip )
    //   {
    //     outFile0 << phaseDens[i][ip] << "\t";
    //   }
    //   for ( localIndex ip = 0; ip < numPhases; ++ip )
    //   {
    //     for ( localIndex ic = 0; ic < numComps; ++ic )
    //     {
    //       outFile0 << phaseCompFrac[i][ip][ic] << "\t";
    //     }
    //   }
    //   outFile0 << std::endl;
    // }
    // outFile0.close();

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
                                  integer const ip_gas,
                                  integer const ip_oil,
                                  integer const ip_water,
                                  real64 const & elevation,
                                  arrayView1d< real64 const > const & phaseContacts )
  {
    integer ip_pres = -1;
    // Choose the primary phase
    if ( numPhases == 1 )
    {
      ip_pres = ( ip_gas >= 0 ) ? ip_gas :
                ( ip_oil >= 0 ) ? ip_oil :
                ( ip_water >= 0 ) ? ip_water : 0;
    }
    else if ( ip_gas >= 0 && ip_oil >= 0 && ip_water >= 0 )
    {
      // Three phase case ( assume phaseContacts[0] = goc, phaseContacts[1] = owc )    
      ip_pres = ( elevation <= phaseContacts[1] ) ? ip_water :
                ( elevation <= phaseContacts[0] ) ? ip_oil : ip_gas;
    }
    else if ( ip_oil >= 0 && ip_water >= 0 )
    {
      // phases = oil + water ( assume phaseContacts[0] = owc )
      ip_pres = ( elevation <= phaseContacts[0] ) ? ip_water : ip_oil;
    }
    else if ( ip_gas >= 0 && ip_water >= 0 )
    {
      // phases = gas + water ( assume phaseContacts[0] = gwc )
      ip_pres = ( elevation <= phaseContacts[0] ) ? ip_water : ip_gas;
    }
    else if ( ip_gas >= 0 && ip_oil >= 0 )
    {
      // phases = gas + oil ( assume phaseContacts[0] = goc )
      ip_pres = ( elevation <= phaseContacts[0] ) ? ip_oil : ip_gas;
    }
    // TODO: deal with the case if ip_pres is still -1 
    std::cout << "ip_pres = " << ip_pres << std::endl;
    return ip_pres;
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
