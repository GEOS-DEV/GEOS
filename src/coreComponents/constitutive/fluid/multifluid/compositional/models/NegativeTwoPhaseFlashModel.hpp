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
 * @file NegativeTwoPhaseFlashModel.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_NEGATIVETWOPHASEFLASHMODEL_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_NEGATIVETWOPHASEFLASHMODEL_HPP_

#include "FunctionBase.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/EquationOfState.hpp"

#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/Utilities.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/FlashData.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/StabilityTest.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/NegativeTwoPhaseFlash.hpp"
#include "common/format/StringUtilities.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

class ModelParameters;

class NegativeTwoPhaseFlashModelUpdate final : public FunctionBaseUpdate
{
private:
  // Tolerance on temperature for phase labelling
  static constexpr real64 temperatureTolerance = 0.1;

public:
  using PhaseProp = MultiFluidVar< real64, 3, constitutive::multifluid::LAYOUT_PHASE, constitutive::multifluid::LAYOUT_PHASE_DC >;
  using PhaseComp = MultiFluidVar< real64, 4, constitutive::multifluid::LAYOUT_PHASE_COMP, constitutive::multifluid::LAYOUT_PHASE_COMP_DC >;
  using Deriv = constitutive::multifluid::DerivativeOffset;

  NegativeTwoPhaseFlashModelUpdate( integer const numComponents,
                                    integer const liquidIndex,
                                    integer const vapourIndex,
                                    EquationOfStateType const liquidEos,
                                    EquationOfStateType const vapourEos,
                                    real64 const salinity,
                                    arrayView1d< real64 const > const componentCriticalVolume,
                                    arrayView1d< real64 const > const continuousFlashParameters,
                                    arrayView1d< integer const > const discreteFlashParameters );

  // Mark as a 2-phase flash
  GEOS_HOST_DEVICE
  static constexpr integer getNumberOfPhases() { return 2; }

  template< int USD1, int USD2 >
  GEOS_HOST_DEVICE
  void compute( ComponentProperties::KernelWrapper const & componentProperties,
                real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & compFraction,
                arraySlice2d< real64, USD2 > const & kValues,
                PhaseProp::SliceType const phaseFraction,
                PhaseComp::SliceType const phaseCompFraction ) const
  {
    constexpr integer maxNumComps = MultiFluidConstants::MAX_NUM_COMPONENTS;

    integer const numDofs = 2 + m_numComponents;

    StackArray< real64, 1, maxNumComps > incipientComposition( m_numComponents );

    integer pureComponent = -1;
    integer almostPureComponent = -1;
    checkPureMixture( m_numComponents, compFraction, pureComponent, almostPureComponent );

    // Check if k-Values need to be initialised
    auto kVapourLiquid = kValues[0];
    bool const needInitialisation = (pureComponent < 0) && hasZero( m_numComponents, kVapourLiquid.toSliceConst() );
    if( needInitialisation )
    {
      if( m_flashData.liquidEos == EquationOfStateType::SoreideWhitson )
      {
        // Initialise the trial compositions using SW uniform values
        KValueInitialization::computeSoreideWhitsonKvalue( m_numComponents,
                                                           pressure,
                                                           temperature,
                                                           componentProperties,
                                                           compFraction,
                                                           kVapourLiquid );
      }
      else
      {
        // Initialise the trial compositions using Wilson k-Values
        KValueInitialization::computeWilsonGasLiquidKvalue( m_numComponents,
                                                            pressure,
                                                            temperature,
                                                            componentProperties,
                                                            kVapourLiquid );
      }
    }

    integer const stabilityIterations = m_discreteFlashParameters[FlashParameters::STABILITY_MAX_ITERATIONS];

    bool unstableMixture = false;
    bool stabilityStatus = (0 <= pureComponent);
    EquationOfStateType incipientEquationOfState = m_flashData.liquidEos;

    if( ((pureComponent < 0) && (0 < stabilityIterations)) || needInitialisation )
    {
      stabilityStatus = StabilityTest::compute( m_numComponents,
                                                pressure,
                                                temperature,
                                                compFraction,
                                                componentProperties,
                                                m_flashData,
                                                kVapourLiquid.toSliceConst(),
                                                m_continuousFlashParameters.toSliceConst(),
                                                m_discreteFlashParameters.toSliceConst(),
                                                unstableMixture,
                                                incipientEquationOfState,
                                                incipientComposition.toSlice() );
    }
    else
    {
      for( integer ic = 0; ic < m_numComponents; ++ic )
      {
        incipientComposition[ic] = compFraction[ic];
      }
    }

    // Pure mixtures are always stable
    if( 0 <= pureComponent )
    {
      unstableMixture = false;
    }

    // If the stability test failed to converge to a stationary point then we will assume the mixture is unstable
    if( unstableMixture || !stabilityStatus )
    {
      // Unstable mixture
      // Iterative solve to converge flash
      bool const flashStatus = NegativeTwoPhaseFlash::compute( m_numComponents,
                                                               pressure,
                                                               temperature,
                                                               compFraction,
                                                               componentProperties,
                                                               m_flashData,
                                                               m_continuousFlashParameters.toSliceConst(),
                                                               m_discreteFlashParameters.toSliceConst(),
                                                               kValues,
                                                               phaseFraction.value[m_vapourIndex],
                                                               phaseCompFraction.value[m_liquidIndex],
                                                               phaseCompFraction.value[m_vapourIndex] );

      GEOS_ERROR_IF( !flashStatus,
                     GEOS_FMT( "Negative two phase flash failed to converge at pressure {:.5e}, temperature {:.3f} and composition {}.",
                               pressure,
                               temperature,
                               stringutilities::concat( "", compFraction ) ) );

      // Calculate derivatives
      NegativeTwoPhaseFlash::computeDerivatives( m_numComponents,
                                                 pressure,
                                                 temperature,
                                                 compFraction,
                                                 componentProperties,
                                                 m_flashData,
                                                 phaseFraction.value[m_vapourIndex],
                                                 phaseCompFraction.value[m_liquidIndex].toSliceConst(),
                                                 phaseCompFraction.value[m_vapourIndex].toSliceConst(),
                                                 phaseFraction.derivs[m_vapourIndex],
                                                 phaseCompFraction.derivs[m_liquidIndex],
                                                 phaseCompFraction.derivs[m_vapourIndex] );

      // Handle the negative flash
      integer singlePhaseIndex = -1;
      real64 const & V = phaseFraction.value[m_vapourIndex];
      if( V < MultiFluidConstants::epsilon )
      {
        singlePhaseIndex = m_liquidIndex;
      }
      else if( 1.0 - V < MultiFluidConstants::epsilon )
      {
        singlePhaseIndex = m_vapourIndex;
      }
      if( 0 <= singlePhaseIndex )
      {
        phaseFraction.value[m_vapourIndex] = (singlePhaseIndex == m_vapourIndex) ? 1.0 : 0.0;
        LvArray::forValuesInSlice( phaseFraction.derivs[m_vapourIndex], setZero );
        LvArray::forValuesInSlice( phaseCompFraction.derivs[singlePhaseIndex], setZero );
        for( integer ic = 0; ic < m_numComponents; ++ic )
        {
          phaseCompFraction.value( singlePhaseIndex, ic ) = compFraction[ic];
          phaseCompFraction.derivs( singlePhaseIndex, ic, Deriv::dC + ic ) = 1.0;
        }
      }
    }
    else
    {
      // Stable mixture - simply label
      // The gas phase has a lower Li-temperature
      real64 const sampleTemperature = calculateLiTemperature( componentProperties, compFraction );
      real64 const incipientTemperature = calculateLiTemperature( componentProperties, incipientComposition.toSliceConst());

      int sampleIndex = m_liquidIndex;
      int incipientIndex = m_vapourIndex;
      if( sampleTemperature < incipientTemperature - temperatureTolerance )
      {
        sampleIndex = m_vapourIndex;
        incipientIndex = m_liquidIndex;
      }
      else if( incipientTemperature < sampleTemperature - temperatureTolerance )
      {
        /* Do nothing */
      }
      else if( sampleTemperature < temperature )
      {
        sampleIndex = m_vapourIndex;
        incipientIndex = m_liquidIndex;
      }

      phaseFraction.value[m_vapourIndex] = (sampleIndex == m_vapourIndex) ? 1.0 : 0.0;

      LvArray::forValuesInSlice( phaseFraction.derivs[m_vapourIndex], setZero );

      EquationOfStateType const sampleEquationOfState = (incipientEquationOfState == m_flashData.liquidEos) ? m_flashData.vapourEos : m_flashData.liquidEos;
      StabilityTest::computeDerivatives( m_numComponents,
                                         pressure,
                                         temperature,
                                         compFraction,
                                         componentProperties,
                                         sampleEquationOfState,
                                         incipientEquationOfState,
                                         m_flashData,
                                         incipientComposition.toSliceConst(),
                                         phaseCompFraction.derivs[incipientIndex],
                                         phaseCompFraction.derivs[sampleIndex] );

      LvArray::forValuesInSlice( phaseCompFraction.derivs[sampleIndex], setZero );
      for( integer ic = 0; ic < m_numComponents; ++ic )
      {
        phaseCompFraction.value( sampleIndex, ic ) = compFraction[ic];
        phaseCompFraction.derivs( sampleIndex, ic, Deriv::dC + ic ) = 1.0;

        phaseCompFraction.value( incipientIndex, ic ) = incipientComposition[ic];
      }
    }

    // Complete by calculating liquid phase fraction
    phaseFraction.value[m_liquidIndex] = 1.0 - phaseFraction.value[m_vapourIndex];
    for( integer ic = 0; ic < numDofs; ic++ )
    {
      phaseFraction.derivs[m_liquidIndex][ic] = -phaseFraction.derivs[m_vapourIndex][ic];
    }
  }

  template< int USD >
  GEOS_HOST_DEVICE
  real64 calculateLiTemperature( ComponentProperties::KernelWrapper const & componentProperties,
                                 arraySlice1d< real64 const, USD > const & composition ) const
  {
    real64 sumVz = 0.0;
    real64 sumVzt = 0.0;
    arrayView1d< real64 const > const & criticalTemperature = componentProperties.m_componentCriticalTemperature;
    for( integer ic = 0; ic < m_numComponents; ++ic )
    {
      real64 const Vz = m_componentCriticalVolume[ic] * composition[ic];
      sumVz += Vz;
      sumVzt += Vz * criticalTemperature[ic];
    }
    return sumVzt / sumVz;
  }

private:
  integer const m_numComponents;
  integer const m_liquidIndex;
  integer const m_vapourIndex;
  FlashData m_flashData;
  arrayView1d< real64 const > const m_componentCriticalVolume;
  arrayView1d< real64 const > const m_continuousFlashParameters;
  arrayView1d< integer const > const m_discreteFlashParameters;
};

class NegativeTwoPhaseFlashModel : public FunctionBase
{
public:
  NegativeTwoPhaseFlashModel( string const & name,
                              ComponentProperties const & componentProperties,
                              ModelParameters const & modelParameters,
                              arrayView1d< integer const > const phaseTypes );

  static string catalogName();

  FunctionType functionType() const override
  {
    return FunctionType::FLASH;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = NegativeTwoPhaseFlashModelUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

  // Create parameters unique to this model
  static std::unique_ptr< ModelParameters > createParameters( std::unique_ptr< ModelParameters > parameters );

  // Determine phase ordering
  static void calculatePhaseOrdering( arrayView1d< integer const > const & phaseTypes,
                                      arrayView1d< integer > const & phaseOrder );

private:
  ModelParameters const & m_parameters;
  arrayView1d< integer const > const m_phaseTypes;
};

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_NEGATIVETWOPHASEFLASHMODEL_HPP_
