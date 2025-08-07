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
   * @param[in] flashData The parameters required for the flash
   * @param[in] kValues the k-values to use to create test samples
   * @param[in] continuousFlashParameters List of continuous (float) parameters for flash
   * @param[in] discreteFlashParameters List of discrete (integer) parameters for flash
   * @param[out] unstableMixture The stability test result indicating if the mixture is unstable
   * @param[out] incipientEquationOfState The equation of state of the incipient phase
   * @param[out] incipientComposition The composition of the incipient phase. This will be
   *             set to the incipient phase composition at the stationary point that is furthest
   *             from the trivial solution.
   * @return a flag indicating that the stability test is successful i.e. mixture is unstable or all test points converged
   */
  template< integer USD1, integer USD2 >
  GEOS_HOST_DEVICE
  static bool compute( integer const numComps,
                       real64 const pressure,
                       real64 const temperature,
                       arraySlice1d< real64 const, USD1 > const & composition,
                       ComponentProperties::KernelWrapper const & componentProperties,
                       FlashData const & flashData,
                       arraySlice1d< real64 const, USD2 > const & kValues,
                       arraySlice1d< real64 const > const & continuousFlashParameters,
                       arraySlice1d< integer const > const & discreteFlashParameters,
                       bool & unstableMixture,
                       EquationOfStateType & incipientEquationOfState,
                       arraySlice1d< real64 > const & incipientComposition );

  /**
   * @brief Compute derivatives of the incipient phase composition
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the mixture
   * @param[in] componentProperties The compositional component properties
   * @param[in] sampleEquationOfState The equation of state for the sample
   * @param[in] sampleEquationOfState The equation of state for the incipient phase
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
                                  EquationOfStateType const & sampleEquationOfState,
                                  EquationOfStateType const & incipientEquationOfState,
                                  FlashData const & flashData,
                                  arraySlice1d< real64 const > const & incipientComposition,
                                  arraySlice2d< real64, USD2 > const & incipientCompositionDerivs,
                                  arraySlice2d< real64, USD2 > const & compositionDerivs )
  {
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
                                                       sampleEquationOfState,
                                                       flashData,
                                                       logFugacity.toSlice(),
                                                       logFugacitySampleDerivs );
    FugacityCalculator::computeLogFugacityDerivatives( numComps,
                                                       pressure,
                                                       temperature,
                                                       incipientComposition,
                                                       componentProperties,
                                                       incipientEquationOfState,
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
      A( ic, ic ) += 1.0;
      A( ic, numComps ) = yi;
      A( numComps, ic ) = 1.0;
    }
  }

  template< int USD1, int USD2, int USD3 >
  GEOS_HOST_DEVICE
  static real64 calculateResidual( integer const numComps,
                                   arraySlice1d< integer const > const & presentComponents,
                                   real64 const pressure,
                                   real64 const temperature,
                                   ComponentProperties::KernelWrapper const & componentProperties,
                                   EquationOfStateType const equationOfState,
                                   FlashData const & flashData,
                                   arraySlice1d< real64 const, USD1 > const & logTestComposition,
                                   real64 & tangentPlaneDistance,
                                   arraySlice1d< real64, USD2 > const & testComposition,
                                   arraySlice1d< real64 > const & testLogFugacity,
                                   arraySlice1d< real64, USD3 > const & residual );

  template< int USD1, int USD2, int USD3, int USD4 >
  GEOS_HOST_DEVICE
  static real64 calculateResidualAndJacobian( integer const numComps,
                                              arraySlice1d< integer const > const & presentComponents,
                                              real64 const pressure,
                                              real64 const temperature,
                                              ComponentProperties::KernelWrapper const & componentProperties,
                                              EquationOfStateType const equationOfState,
                                              FlashData const & flashData,
                                              arraySlice1d< real64 const, USD1 > const & logTestComposition,
                                              real64 & tangentPlaneDistance,
                                              arraySlice1d< real64, USD2 > const & testComposition,
                                              arraySlice1d< real64 > const & testFugacity,
                                              arraySlice2d< real64 > const & logTestFugacityDerivs,
                                              arraySlice1d< real64, USD3 > const & residual,
                                              arraySlice2d< real64, USD4 > const & jacobian );
};

} // namespace compositional

} // namespace constitutive

} // namespace geos

// Include the implementation
#include "StabilityTest_impl.hpp"

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_STABILITYTEST_HPP_
