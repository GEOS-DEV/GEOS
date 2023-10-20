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
 * @file PVTDriverRunTest.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_PVTDRIVERRUNTEST_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_PVTDRIVERRUNTEST_HPP_

#include "PVTDriver.hpp"

#include "constitutive/fluid/multifluid/Layouts.hpp"

namespace geos
{

template< typename FLUID_TYPE >
void PVTDriver::runTest( FLUID_TYPE & fluid, arrayView2d< real64 > const & table )
{
  // get number of phases and components

  integer const numComponents = fluid.numFluidComponents();
  integer const numPhases = fluid.numFluidPhases();

  // prefer output in mass

  fluid.setMassFlag( true );

  // create kernel wrapper

  typename FLUID_TYPE::KernelWrapper const kernelWrapper = fluid.createKernelWrapper();

  // set composition to user specified feed
  // it is more convenient to provide input in molar, so perform molar to mass conversion here

  GEOS_ASSERT_EQ( numComponents, m_feed.size() );
  array2d< real64, compflow::LAYOUT_COMP > const compositionValues( 1, numComponents );

  // extract the molecular weights
  array1d< real64 > const componentMolarWeights( numComponents );
  for( integer i = 0; i < numComponents; ++i )
  {
    componentMolarWeights[i] = fluid.componentMolarWeights()[i];
  }

  real64 sum = 0.0;
  for( integer i = 0; i < numComponents; ++i )
  {
    compositionValues[0][i] = m_feed[i] * componentMolarWeights[i];
    sum += compositionValues[0][i];
  }
  for( integer i = 0; i < numComponents; ++i )
  {
    compositionValues[0][i] /= sum;
  }

  arrayView2d< real64 const, compflow::USD_COMP > const composition = compositionValues;

  // perform fluid update using table (P,T) and save resulting total density, etc.
  // note: column indexing should be kept consistent with output file header below.

  integer const numSteps = m_numSteps;
  using ExecPolicy = typename FLUID_TYPE::exec_policy;
  forAll< ExecPolicy >( composition.size( 0 ),
                        [numPhases, numComponents, numSteps, kernelWrapper, table, composition,
                         molarWeights=componentMolarWeights.toViewConst(),
                         outputCompressibility=m_outputCompressibility,
                         outputPhaseComposition=m_outputPhaseComposition]
                        GEOS_HOST_DEVICE ( localIndex const i )
  {
    constexpr real64 epsilon = LvArray::NumericLimits< real64 >::epsilon;

    // Pressure derivative index
    constexpr integer dPIndex = constitutive::multifluid::DerivativeOffset::dP;

    // Index for start of phase properties
    integer const PHASE = outputCompressibility != 0 ? TEMP + 3 : TEMP + 2;

    for( integer n = 0; n <= numSteps; ++n )
    {
      kernelWrapper.update( i, 0, table( n, PRES ), table( n, TEMP ), composition[i] );
      real64 const totalDensity = kernelWrapper.totalDensity()( i, 0 );
      table( n, TEMP + 1 ) = totalDensity;

      if( outputCompressibility != 0 )
      {
        real64 compressibility = kernelWrapper.dTotalDensity()( i, 0, dPIndex ) / totalDensity;
        table( n, TEMP + 2 ) =  LvArray::math::max( 0.0, compressibility );
      }

      for( integer p = 0; p < numPhases; ++p )
      {
        table( n, PHASE + p ) = kernelWrapper.phaseFraction()( i, 0, p );
        table( n, PHASE + p + numPhases ) = kernelWrapper.phaseDensity()( i, 0, p );
        table( n, PHASE + p + 2 * numPhases ) = kernelWrapper.phaseViscosity()( i, 0, p );
      }
      if( outputPhaseComposition != 0 )
      {
        for( integer p = 0; p < numPhases; ++p )
        {
          integer const compStartIndex = PHASE + 3 * numPhases + p * numComponents;

          // Convert mass fractions to molar fractions
          real64 molarSum = 0.0;
          for( integer ic = 0; ic < numComponents; ++ic )
          {
            real64 const massFraction = kernelWrapper.phaseCompFraction()( i, 0, p, ic );
            table( n, compStartIndex + ic ) = massFraction / molarWeights[ic];
            molarSum += table( n, compStartIndex + ic );
          }

          real64 const invMolarSum = epsilon < molarSum ? 1.0/molarSum : 0.0;

          for( integer ic = 0; ic < numComponents; ++ic )
          {
            table( n, compStartIndex + ic ) *= invMolarSum;
          }
        }
      }
    }
  } );

}

}

#endif /* GEOS_CONSTITUTIVE_PVTDRIVERRUNTEST_HPP_ */
