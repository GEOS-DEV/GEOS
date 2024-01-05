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

  template< int USD1 >
  GEOS_HOST_DEVICE
  void compute( ComponentProperties::KernelWrapper const & componentProperties,
                real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & compFraction,
                PhaseProp::SliceType const phaseFraction,
                PhaseComp::SliceType const phaseCompFraction ) const
  {
    integer const numDofs = 2 + m_numComponents;

    // Iterative solve to converge flash
    NegativeTwoPhaseFlash::compute< EOS_TYPE_LIQUID, EOS_TYPE_VAPOUR >(
      m_numComponents,
      pressure,
      temperature,
      compFraction,
      componentProperties,
      phaseFraction.value[m_vapourIndex],
      phaseCompFraction.value[m_liquidIndex],
      phaseCompFraction.value[m_vapourIndex] );

    // Calculate derivatives
    NegativeTwoPhaseFlash::computeDerivatives< EOS_TYPE_LIQUID, EOS_TYPE_VAPOUR >(
      m_numComponents,
      pressure,
      temperature,
      compFraction,
      componentProperties,
      phaseFraction.value[m_vapourIndex],
      phaseCompFraction.value[m_liquidIndex],
      phaseCompFraction.value[m_vapourIndex],
      phaseFraction.derivs[m_vapourIndex],
      phaseCompFraction.derivs[m_liquidIndex],
      phaseCompFraction.derivs[m_vapourIndex] );

    // Complete by calculating liquid phase fraction
    phaseFraction.value[m_liquidIndex] = 1.0 - phaseFraction.value[m_vapourIndex];
    for( integer ic = 0; ic < numDofs; ic++ )
    {
      phaseFraction.derivs[m_liquidIndex][ic] = -phaseFraction.derivs[m_vapourIndex][ic];
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
