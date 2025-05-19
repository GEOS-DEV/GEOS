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
 * @file HydrostaticPressureKernel.cpp
 */

#include "HydrostaticPressureKernel.hpp"
#include "constitutive/fluid/multifluid/MultiFluidSelector.hpp"

namespace geos
{
namespace isothermalCompositionalMultiphaseBaseKernels
{

static void calculateDensity( integer const numPhases,
                              constitutive::MultiFluidBase & fluid,
                              real64 const pressure,
                              real64 const temperature,
                              arraySlice1d< real64 const, compflow::USD_COMP - 1 > const composition,
                              arraySlice1d< real64 > const massDensity,
                              HydrostaticPressureStackVariables & stackVars )
{
  constitutive::constitutiveUpdatePassThru( fluid, [&] ( auto & castFluid )
  {
    using FluidType = TYPEOFREF( castFluid );
    using FluidUpdateType = typename FluidType::KernelWrapper;
    FluidUpdateType fluidWrapper = castFluid.createKernelWrapper();
    HydrostaticPressureFluidUpdate< FluidUpdateType >::getMassDensity( numPhases,
                                                                       fluidWrapper,
                                                                       pressure,
                                                                       temperature,
                                                                       composition,
                                                                       massDensity,
                                                                       stackVars );

  } );
}

static HydrostaticPressureKernel::ReturnType
computeHydrostaticPressure( integer const numComps,
                            integer const numPhases,
                            integer const ipInit,
                            integer const maxNumEquilIterations,
                            real64 const & equilTolerance,
                            real64 const (&gravVector)[ 3 ],
                            constitutive::MultiFluidBase & fluid,
                            arrayView1d< TableFunction::KernelWrapper const > compFracTableWrappers,
                            TableFunction::KernelWrapper tempTableWrapper,
                            real64 const & refElevation,
                            real64 const & refPres,
                            arraySlice1d< real64 const > const & refPhaseMassDens,
                            real64 const & newElevation,
                            real64 & newPres,
                            arraySlice1d< real64 > const & newPhaseMassDens,
                            HydrostaticPressureStackVariables & stackVars )
{
  // fluid properties at this elevation
  StackArray< real64, 2, constitutive::MultiFluidBase::MAX_NUM_COMPONENTS, compflow::LAYOUT_COMP > compFrac( 1, numComps );

  auto const phaseFraction = stackVars.phaseFraction[0][0].toSliceConst();

  bool isSinglePhaseFlow = true;

  // Step 1: compute the hydrostatic pressure at the current elevation

  real64 const gravCoef = gravVector[2] * ( refElevation - newElevation );
  real64 const temp = tempTableWrapper.compute( &newElevation );
  for( integer ic = 0; ic < numComps; ++ic )
  {
    compFrac[0][ic] = compFracTableWrappers[ic].compute( &newElevation );
  }

  // Step 2: guess the pressure with the refPhaseMassDensity

  real64 pres0 = refPres - refPhaseMassDens[ipInit] * gravCoef;
  real64 pres1 = 0.0;

  // Step 3: compute the mass density at this elevation using the guess, and update pressure
  calculateDensity( numPhases,
                    fluid,
                    pres0,
                    temp,
                    compFrac[0],
                    newPhaseMassDens,
                    stackVars );

  pres1 = refPres - 0.5 * ( refPhaseMassDens[ipInit] + newPhaseMassDens[ipInit] ) * gravCoef;

  // Step 4: fixed-point iteration until convergence

  bool equilHasConverged = false;
  for( integer eqIter = 0; eqIter < maxNumEquilIterations; ++eqIter )
  {

    // check convergence
    equilHasConverged = ( LvArray::math::abs( pres0 - pres1 ) < equilTolerance );
    pres0 = pres1;

    // if converged, check number of phases and move on
    if( equilHasConverged )
    {
      // make sure that the fluid is single-phase, other we have to issue a warning (for now)
      // if only one phase is mobile, we are in good shape (unfortunately it is hard to access relperm from here)
      localIndex numberOfPhases = 0;
      for( integer ip = 0; ip < numPhases; ++ip )
      {
        if( constitutive::MultiFluidConstants::minForSpeciesPresence < phaseFraction[ip] )
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
    calculateDensity( numPhases,
                      fluid,
                      pres0,
                      temp,
                      compFrac[0],
                      newPhaseMassDens,
                      stackVars );

    pres1 = refPres - 0.5 * ( refPhaseMassDens[ipInit] + newPhaseMassDens[ipInit] ) * gravCoef;
  }

  // Step 5: save the hydrostatic pressure and the corresponding density

  newPres = pres1;

  if( !equilHasConverged )
  {
    return HydrostaticPressureKernel::ReturnType::FAILED_TO_CONVERGE;
  }
  else if( !isSinglePhaseFlow )
  {
    return HydrostaticPressureKernel::ReturnType::DETECTED_MULTIPHASE_FLOW;
  }
  else
  {
    return HydrostaticPressureKernel::ReturnType::SUCCESS;
  }
}

HydrostaticPressureKernel::ReturnType
HydrostaticPressureKernel::launch( localIndex const size,
                                   integer const numComps,
                                   integer const numPhases,
                                   integer const ipInit,
                                   integer const maxNumEquilIterations,
                                   real64 const equilTolerance,
                                   real64 const (&gravVector)[ 3 ],
                                   real64 const & minElevation,
                                   real64 const & elevationIncrement,
                                   real64 const & datumElevation,
                                   real64 const & datumPres,
                                   constitutive::MultiFluidBase & fluid,
                                   arrayView1d< TableFunction::KernelWrapper const > compFracTableWrappers,
                                   TableFunction::KernelWrapper tempTableWrapper,
                                   arrayView1d< arrayView1d< real64 > const > elevationValues,
                                   arrayView1d< real64 > pressureValues )
{
  ReturnType returnVal = ReturnType::SUCCESS;

  // Temporary workspace
  HydrostaticPressureStackVariables stackVars( numPhases, numComps );

  // Step 1: compute the phase mass densities at datum

  // Datum fluid properties
  array2d< real64, compflow::LAYOUT_COMP > datumCompFrac( 1, numComps );
  array1d< real64 > datumPhaseMassDens( numPhases );

  real64 const datumTemp = tempTableWrapper.compute( &datumElevation );
  for( integer ic = 0; ic < numComps; ++ic )
  {
    datumCompFrac[0][ic] = compFracTableWrappers[ic].compute( &datumElevation );
  }
  calculateDensity( numPhases,
                    fluid,
                    datumPres,
                    datumTemp,
                    datumCompFrac[0],
                    datumPhaseMassDens,
                    stackVars );

  // Step 2: find the closest elevation to datumElevation

  forAll< parallelHostPolicy >( size, [=] ( localIndex const i )
  {
    real64 const elevation = minElevation + i * elevationIncrement;
    elevationValues[0][i] = elevation;
  } );
  integer const iRef = LvArray::sortedArrayManipulation::find( elevationValues[0].begin(),
                                                               elevationValues[0].size(),
                                                               datumElevation );

  // Step 3: compute the mass density and pressure at the reference elevation

  array2d< real64 > phaseMassDens( pressureValues.size(), numPhases );

  ReturnType const refReturnVal =
    computeHydrostaticPressure( numComps,
                                numPhases,
                                ipInit,
                                maxNumEquilIterations,
                                equilTolerance,
                                gravVector,
                                fluid,
                                compFracTableWrappers,
                                tempTableWrapper,
                                datumElevation,
                                datumPres,
                                datumPhaseMassDens,
                                elevationValues[0][iRef],
                                pressureValues[iRef],
                                phaseMassDens[iRef],
                                stackVars );
  if( refReturnVal == ReturnType::FAILED_TO_CONVERGE )
  {
    return ReturnType::FAILED_TO_CONVERGE;
  }
  else if( refReturnVal == ReturnType::DETECTED_MULTIPHASE_FLOW )
  {
    returnVal = ReturnType::DETECTED_MULTIPHASE_FLOW;
  }

  // Step 4: for each elevation above the reference elevation, compute the pressure

  localIndex const numEntriesAboveRef = size - iRef - 1;
  forAll< serialPolicy >( numEntriesAboveRef, [=, &fluid, &stackVars, &returnVal] ( localIndex const i )
  {
    ReturnType const returnValAboveRef =
      computeHydrostaticPressure( numComps,
                                  numPhases,
                                  ipInit,
                                  maxNumEquilIterations,
                                  equilTolerance,
                                  gravVector,
                                  fluid,
                                  compFracTableWrappers,
                                  tempTableWrapper,
                                  elevationValues[0][iRef+i],
                                  pressureValues[iRef+i],
                                  phaseMassDens[iRef+i],
                                  elevationValues[0][iRef+i+1],
                                  pressureValues[iRef+i+1],
                                  phaseMassDens[iRef+i+1],
                                  stackVars );
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

  localIndex const numEntriesBelowRef = iRef;
  forAll< serialPolicy >( numEntriesBelowRef, [=, &fluid, &stackVars, &returnVal] ( localIndex const i )
  {
    ReturnType const returnValBelowRef =
      computeHydrostaticPressure( numComps,
                                  numPhases,
                                  ipInit,
                                  maxNumEquilIterations,
                                  equilTolerance,
                                  gravVector,
                                  fluid,
                                  compFracTableWrappers,
                                  tempTableWrapper,
                                  elevationValues[0][iRef-i],
                                  pressureValues[iRef-i],
                                  phaseMassDens[iRef-i],
                                  elevationValues[0][iRef-i-1],
                                  pressureValues[iRef-i-1],
                                  phaseMassDens[iRef-i-1],
                                  stackVars );
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

} // namespace isothermalCompositionalMultiphaseBaseKernels
} // namespace geos
