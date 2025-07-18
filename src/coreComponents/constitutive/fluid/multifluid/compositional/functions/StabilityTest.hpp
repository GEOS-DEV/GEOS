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
 * @file StabilityTest.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_STABILITYTEST_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_STABILITYTEST_HPP_

#include "KValueInitialization.hpp"
#include "FugacityCalculator.hpp"
#include "Utilities.hpp"
#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/ComponentProperties.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/FlashParameters.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/FlashData.hpp"
#include "common/TimingMacros.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

struct StabilityTest
{
private:
  static constexpr integer maxNumComps = MultiFluidConstants::MAX_NUM_COMPONENTS;
  using Deriv = constitutive::multifluid::DerivativeOffset;
public:
  /**
   * @brief Perform a two-phase stability test
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the mixture
   * @param[in] componentProperties The compositional component properties
   * @param[in] equationOfState The equation of state
   * @param[in] flashData The parameters required for the flash
   * @param[in] kValues the k-values to use to create test samples
   * @param[in] continuousFlashParameters List of continuous (float) parameters for flash
   * @param[in] discreteFlashParameters List of discrete (integer) parameters for flash
   * @param[out] tangentPlaneDistance the minimum tangent plane distance (TPD)
   * @param[out] incipientComposition The composition of the incipient phase. This will be
   *             set to the incipient phase composition at the stationary point that is furthest
   *             from the trivial solution.
   * @return a flag indicating that mixture is unstable or all test points converged
   */
  template< integer USD1, integer USD2 >
  GEOS_HOST_DEVICE
  static bool compute( integer const numComps,
                       real64 const pressure,
                       real64 const temperature,
                       arraySlice1d< real64 const, USD1 > const & composition,
                       ComponentProperties::KernelWrapper const & componentProperties,
                       EquationOfStateType const & equationOfState,
                       FlashData const & flashData,
                       arraySlice1d< real64 const, USD2 > const & kValues,
                       arraySlice1d< real64 const > const & continuousFlashParameters,
                       arraySlice1d< integer const > const & discreteFlashParameters,
                       real64 & tangentPlaneDistance,
                       arraySlice1d< real64 > const & incipientComposition )
  {
    GEOS_MARK_FUNCTION;

    stackArray2d< real64, 4*maxNumComps > workSpace( 4, numComps );
    arraySlice1d< real64 > logFugacity = workSpace[0];
    arraySlice1d< real64 > normalizedComposition = workSpace[1];
    arraySlice1d< real64 > logTrialComposition = workSpace[2];
    arraySlice1d< real64 > hyperplane = workSpace[3]; // h-parameter
    stackArray1d< integer, maxNumComps > availableComponents( numComps );

    //integer const presentCount =
    calculatePresentComponents( numComps, composition, availableComponents );
    auto const presentComponents = availableComponents.toSliceConst();

    LvArray::forValuesInSlice( workSpace.toSlice(), setZero );

    // Extract flash parameters
    integer const maxIterations = discreteFlashParameters[FlashParameters::STABILITY_MAX_ITERATIONS];
    real64 const stabilityThreshold = continuousFlashParameters[FlashParameters::STABILITY_THRESHOLD];
    real64 const stabilityTolerance = continuousFlashParameters[FlashParameters::STABILITY_TOLERANCE];

    // Calculate the hyperplane parameter
    // h_i = log( z_i ) + log( phi_i )
    FugacityCalculator::computeLogFugacity( numComps,
                                            pressure,
                                            temperature,
                                            composition,
                                            componentProperties,
                                            equationOfState,
                                            flashData,
                                            logFugacity );
    for( integer const ic : presentComponents )
    {
      hyperplane[ic] = LvArray::math::log( composition[ic] ) + logFugacity[ic];
    }

    // The mimimum TPD over all trial compositions
    tangentPlaneDistance = LvArray::NumericLimits< real64 >::max;

    // Flag to indicate all trial compositions converged
    bool allConverged = true;

    // Measure of distance to trivial solution
    real64 maxDistanceToTrivialSolution = 0.0;

    for( real64 const alpha : { 1.0, -1.0 } )
    {
      // Initialise next sample
      for( integer const ic : presentComponents )
      {
        logTrialComposition[ic] = LvArray::math::log( composition[ic] ) + alpha*LvArray::math::log( kValues[ic] );
        normalizedComposition[ic] = LvArray::math::exp( logTrialComposition[ic] );
      }
      // Start iterations for this sample
      bool converged = false;
      for( localIndex iterationCount = 0; iterationCount < maxIterations; ++iterationCount )
      {
        // Normalise the composition and calculate the fugacity
        real64 const totalMoles = normalizeComposition( numComps, normalizedComposition );
        FugacityCalculator::computeLogFugacity( numComps,
                                                pressure,
                                                temperature,
                                                normalizedComposition.toSliceConst(),
                                                componentProperties,
                                                equationOfState,
                                                flashData,
                                                logFugacity );

        // Calculate the TPD and stationarity
        real64 tpd = 0.0;
        real64 error = 0.0;
        for( integer const ic : presentComponents )
        {
          real64 const dG = logTrialComposition[ic] + logFugacity[ic] - hyperplane[ic];
          tpd += composition[ic] + totalMoles * normalizedComposition[ic] * (dG - 1.0);
          error += (dG*dG);
        }
        error = LvArray::math::sqrt( error );
        if( tpd < tangentPlaneDistance )
        {
          tangentPlaneDistance = tpd;
        }
        std::cout
          << std::fixed << std::setprecision( 1 ) << std::setw( 4 ) << alpha << " "
          << std::setw( 3 ) << iterationCount << " "
          << std::fixed << std::setprecision( 6 ) << std::setw( 4 ) << normalizedComposition << " "
          << std::scientific << std::setprecision( 6 ) << logFugacity << " "
          << std::scientific << std::setprecision( 6 ) << std::setw( 13 ) << tpd << " "
          << std::scientific << std::setprecision( 6 ) << std::setw( 13 ) << error << " "
          << "\n";

        // Check stationarity
        if( error < stabilityTolerance )
        {
          converged = true;
          break;
        }

        // Update to next step
        for( integer const ic : presentComponents )
        {
          logTrialComposition[ic] = hyperplane[ic] - logFugacity[ic];
          normalizedComposition[ic] = LvArray::math::exp( logTrialComposition[ic] );
        }
      }
      allConverged = allConverged && converged;

      // Calculate the distance to the trivial solution
      real64 distance = 0.0;
      for( integer const ic : presentComponents )
      {
        real64 const dZ = normalizedComposition[ic] - composition[ic];
        distance += (dZ*dZ);
      }
      if( maxDistanceToTrivialSolution < distance )
      {
        maxDistanceToTrivialSolution = distance;
        for( integer ic = 0; ic < numComps; ++ic )
        {
          incipientComposition[ic] = normalizedComposition[ic];
        }
      }
      if( tangentPlaneDistance < stabilityThreshold )
      {
        break;
      }
    }

    // The test is successful if either we have a negative TPD or all test compositions converged to stationarity
    return (tangentPlaneDistance < stabilityThreshold) || (allConverged);
  }

  /**
   * @brief Compute derivatives of the incipient phase composition
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the mixture
   * @param[in] componentProperties The compositional component properties
   * @param[in] equationOfState The equation of state
   * @param[in] flashData The parameters required for the flash
   * @param[in] incipientComposition The composition of the incipient phase
   * @param[out] incipientCompositionDerivs Derivatives of the composition of the incipient phase
   * @param[out] compositionDerivs Workspace
   */
  template< integer USD1, integer USD2 >
  GEOS_HOST_DEVICE
  static void computeDerivatives( integer const numComps,
                                  real64 const pressure,
                                  real64 const temperature,
                                  arraySlice1d< real64 const, USD1 > const & composition,
                                  ComponentProperties::KernelWrapper const & componentProperties,
                                  EquationOfStateType const & equationOfState,
                                  FlashData const & flashData,
                                  arraySlice1d< real64 const > const & incipientComposition,
                                  arraySlice2d< real64, USD2 > const & incipientCompositionDerivs,
                                  arraySlice2d< real64, USD2 > const & compositionDerivs )
  {
    GEOS_MARK_FUNCTION;

    integer constexpr maxNumRows = MultiFluidConstants::MAX_NUM_COMPONENTS + 1;
    integer constexpr maxDofs = MultiFluidConstants::MAX_NUM_COMPONENTS + 2;

    integer const numDofs = 2 + numComps;

    StackArray< real64, 1, maxNumComps > logFugacity( numComps );
    auto const & logFugacityIncipientDerivs = incipientCompositionDerivs;
    auto const & logFugacitySampleDerivs = compositionDerivs;

    FugacityCalculator::computeLogFugacityDerivatives( numComps,
                                                       pressure,
                                                       temperature,
                                                       composition,
                                                       componentProperties,
                                                       equationOfState,
                                                       flashData,
                                                       logFugacity.toSlice(),
                                                       logFugacitySampleDerivs );
    FugacityCalculator::computeLogFugacityDerivatives( numComps,
                                                       pressure,
                                                       temperature,
                                                       incipientComposition,
                                                       componentProperties,
                                                       equationOfState,
                                                       flashData,
                                                       logFugacity.toSlice(),
                                                       logFugacityIncipientDerivs );

    StackArray< real64, 2, maxNumRows * maxNumRows, MatrixLayout::COL_MAJOR_PERM > A( numComps + 1, numComps + 1 );
    StackArray< real64, 2, maxNumRows * maxDofs, MatrixLayout::COL_MAJOR_PERM > X( numComps + 1, numDofs );

    for( integer ic = 0; ic < numComps; ++ic )
    {
      real64 const yi = incipientComposition[ic];
      for( integer jc = 0; jc < numComps; ++jc )
      {
        real64 const df_dyj = logFugacityIncipientDerivs( ic, Deriv::dC+jc );
        A( ic, jc ) = yi*df_dyj;
      }
      A( ic, ic ) += 1.0;
      A( ic, numComps ) = yi;
      A( numComps, ic ) = 1.0;
    }
    A( numComps, numComps ) = 0.0;

    for( integer ic = 0; ic < numComps; ++ic )
    {
      real64 const yi = incipientComposition[ic];
      for( integer const idof : {Deriv::dP, Deriv::dT} )
      {
        X( ic, idof ) = yi*(-logFugacityIncipientDerivs( ic, idof ) + logFugacitySampleDerivs( ic, idof ));
      }
      for( integer jc = 0; jc < numComps; ++jc )
      {
        integer const idof = Deriv::dC + jc;
        X( ic, idof ) = yi*logFugacitySampleDerivs( ic, idof );
      }
      if( MultiFluidConstants::epsilon < composition[ic] )
      {
        X( ic, Deriv::dC+ic ) += yi/composition[ic];
      }
    }
    for( integer idof = 0; idof < numDofs; ++idof )
    {
      X( numComps, idof ) = 0.0;
    }

    // Solve linear system
    solveLinearSystem( A.toSlice(), X.toSlice() );

    for( integer idof = 0; idof < numDofs; ++idof )
    {
      for( integer ic = 0; ic < numComps; ++ic )
      {
        incipientCompositionDerivs( ic, idof ) = X( ic, idof );
      }
    }
  }
};

} // namespace compositional

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_STABILITYTEST_HPP_
