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
 * @file NegativeTwoPhaseFlashModel.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_NEGATIVETWOPHASEFLASHMODEL_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_NEGATIVETWOPHASEFLASHMODEL_HPP_

#include "FunctionBase.hpp"

#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/CubicEOSPhaseModel.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/NegativeTwoPhaseFlash.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/StabilityTest.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

template< typename EOS_TYPE_LIQUID, typename EOS_TYPE_VAPOUR >
class NegativeTwoPhaseFlashModelUpdate final : public FunctionBaseUpdate
{
public:
  using PhaseProp = MultiFluidVar< real64, 3, multifluid::LAYOUT_PHASE, multifluid::LAYOUT_PHASE_DC >;
  using PhaseComp = MultiFluidVar< real64, 4, multifluid::LAYOUT_PHASE_COMP, multifluid::LAYOUT_PHASE_COMP_DC >;

  NegativeTwoPhaseFlashModelUpdate( integer const numComponents,
                                    integer const liquidIndex,
                                    integer const vapourIndex );

  // Mark as a 2-phase flash
  GEOS_HOST_DEVICE
  static constexpr integer getNumberOfPhases() { return 2; }

  template< int USD1, int USD2 >
  GEOS_HOST_DEVICE
  void compute( ComponentProperties::KernelWrapper const & componentProperties,
                real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & composition,
                arraySlice2d< real64, USD2 > const & kValues,
                PhaseProp::SliceType const phaseFraction,
                PhaseComp::SliceType const phaseComposition ) const
  {
    integer const numDofs = 2 + m_numComponents;

    // Start with a stability test to check for single phase fluid
    // The stability will initialise the k-values if successful
    real64 tangentPlaneDistance = LvArray::NumericLimits< real64 >::max;
    bool const stabilitySuccess = StabilityTest::compute< EOS_TYPE_LIQUID >(
      m_numComponents,
      pressure,
      temperature,
      composition,
      componentProperties,
      tangentPlaneDistance,
      kValues[0] );

    // If the stability test is not successful, the k-values will be set to Wilson k-values which we can then take into the flash
    bool const isSinglePhase = stabilitySuccess && (-MultiFluidConstants::fugacityTolerance < tangentPlaneDistance);

    if( isSinglePhase )
    {
      // Simply label using the Li-correlation
      applyLiCorrelationLabel(
        temperature,
        composition,
        componentProperties,
        phaseFraction.value[m_vapourIndex],
        phaseComposition.value[m_liquidIndex],
        phaseComposition.value[m_vapourIndex] );
    }
    else
    {
      // Iterative solve to converge flash
      NegativeTwoPhaseFlash::compute< EOS_TYPE_LIQUID, EOS_TYPE_VAPOUR >(
        m_numComponents,
        pressure,
        temperature,
        composition,
        componentProperties,
        kValues,
        phaseFraction.value[m_vapourIndex],
        phaseComposition.value[m_liquidIndex],
        phaseComposition.value[m_vapourIndex] );
    }

    // Calculate derivatives
    NegativeTwoPhaseFlash::computeDerivatives< EOS_TYPE_LIQUID, EOS_TYPE_VAPOUR >(
      m_numComponents,
      pressure,
      temperature,
      composition,
      componentProperties,
      phaseFraction.value[m_vapourIndex],
      phaseComposition.value[m_liquidIndex].toSliceConst(),
      phaseComposition.value[m_vapourIndex].toSliceConst(),
      phaseFraction.derivs[m_vapourIndex],
      phaseComposition.derivs[m_liquidIndex],
      phaseComposition.derivs[m_vapourIndex] );

    // Complete by calculating liquid phase fraction
    phaseFraction.value[m_liquidIndex] = 1.0 - phaseFraction.value[m_vapourIndex];
    for( integer ic = 0; ic < numDofs; ic++ )
    {
      phaseFraction.derivs[m_liquidIndex][ic] = -phaseFraction.derivs[m_vapourIndex][ic];
    }
  }

private:
  /**
   * @brief Calculate the label of a mixture using the Li correlation
   * @param[in] numComps number of components
   * @param[in] temperature temperature
   * @param[in] composition composition of the mixture
   * @param[in] componentProperties The compositional component properties
   * @param[out] vapourPhaseMoleFraction the calculated vapour (gas) mole fraction
   * @param[out] liquidComposition the calculated liquid phase composition
   * @param[out] vapourComposition the calculated vapour phase composition
   */
  template< int USD1, int USD2 >
  GEOS_HOST_DEVICE
  void applyLiCorrelationLabel( real64 const & temperature,
                                arraySlice1d< real64 const, USD1 > const & composition,
                                ComponentProperties::KernelWrapper const & componentProperties,
                                real64 & vapourPhaseMoleFraction,
                                arraySlice1d< real64, USD2 > const & liquidComposition,
                                arraySlice1d< real64, USD2 > const & vapourComposition ) const
  {
    auto const & criticalTemperature = componentProperties.m_componentCriticalTemperature;
    auto const & criticalVolume = componentProperties.m_componentCriticalVolume;

    real64 sumTV = 0.0;
    real64 sumV = 0.0;
    for( integer ic = 0; ic < m_numComponents; ++ic )
    {
      sumV += composition[ic] * criticalVolume[ic];
      sumTV += composition[ic] * criticalVolume[ic] * criticalTemperature[ic];
    }
    real64 const temperatureLi = sumTV / sumV;

    if( temperature < temperatureLi )
    {
      // Liquid
      vapourPhaseMoleFraction = 0.0;
    }
    else
    {
      // Vapour
      vapourPhaseMoleFraction = 1.0;
    }
    // In either case we set both compositions equal to the feed
    for( integer ic = 0; ic < m_numComponents; ++ic )
    {
      liquidComposition[ic] = composition[ic];
      vapourComposition[ic] = composition[ic];
    }
  }

private:
  integer const m_numComponents;
  integer const m_liquidIndex{0};
  integer const m_vapourIndex{1};
};

template< typename EOS_TYPE_LIQUID, typename EOS_TYPE_VAPOUR >
class NegativeTwoPhaseFlashModel : public FunctionBase
{
  // Currently limit to the two EOS types being the same
  static_assert( std::is_same< EOS_TYPE_LIQUID, EOS_TYPE_VAPOUR >::value );

public:
  NegativeTwoPhaseFlashModel( string const & name,
                              ComponentProperties const & componentProperties );

  static string catalogName();

  FunctionType functionType() const override
  {
    return FunctionType::FLASH;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = NegativeTwoPhaseFlashModelUpdate< EOS_TYPE_LIQUID, EOS_TYPE_VAPOUR >;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;
};

using NegativeTwoPhaseFlashPRPR = NegativeTwoPhaseFlashModel<
  CubicEOSPhaseModel< PengRobinsonEOS >,
  CubicEOSPhaseModel< PengRobinsonEOS > >;
using NegativeTwoPhaseFlashSRKSRK = NegativeTwoPhaseFlashModel<
  CubicEOSPhaseModel< SoaveRedlichKwongEOS >,
  CubicEOSPhaseModel< SoaveRedlichKwongEOS > >;

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_NEGATIVETWOPHASEFLASHMODEL_HPP_
