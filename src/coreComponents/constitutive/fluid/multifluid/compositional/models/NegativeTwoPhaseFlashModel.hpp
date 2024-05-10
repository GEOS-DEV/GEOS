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

#include "constitutive/fluid/multifluid/compositional/functions/NegativeTwoPhaseFlash.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

class NegativeTwoPhaseFlashModelUpdate final : public FunctionBaseUpdate
{
public:
  NegativeTwoPhaseFlashModelUpdate( integer const numComponents,
                                    integer const liquidIndex,
                                    integer const vapourIndex );

  // Mark as a 2-phase flash
  GEOS_HOST_DEVICE
  static constexpr integer getNumberOfPhases() { return 2; }

  GEOS_HOST_DEVICE
  void compute( ComponentProperties::KernelWrapper const & componentProperties,
                EquationOfState::KernelWrapper const & equationOfState,
                real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const > const & composition,
                multifluid::PhaseComp::SliceType::ValueType const & kValues,
                multifluid::PhaseProp::SliceType const phaseFraction,
                multifluid::PhaseComp::SliceType const phaseComposition ) const
  {
    integer const numDofs = m_numComponents + 2;

    NegativeTwoPhaseFlash::compute( m_numComponents,
                                    pressure,
                                    temperature,
                                    composition,
                                    componentProperties,
                                    equationOfState,
                                    kValues,
                                    phaseFraction.value[m_vapourIndex],
                                    phaseComposition.value[m_liquidIndex],
                                    phaseComposition.value[m_vapourIndex] );

    NegativeTwoPhaseFlash::computeDerivatives( m_numComponents,
                                               pressure,
                                               temperature,
                                               composition,
                                               componentProperties,
                                               equationOfState,
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
  integer const m_numComponents;
  integer const m_liquidIndex{0};
  integer const m_vapourIndex{1};
};

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
  using KernelWrapper = NegativeTwoPhaseFlashModelUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;
};

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_NEGATIVETWOPHASEFLASHMODEL_HPP_
