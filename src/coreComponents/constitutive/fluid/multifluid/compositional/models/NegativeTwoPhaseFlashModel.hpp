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

  explicit NegativeTwoPhaseFlashModelUpdate( ComponentProperties const & componentProperties );

  template< int USD1, int USD2, int USD3 >
  GEOS_HOST_DEVICE
  void compute( real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & compFraction,
                arraySlice1d< real64, USD2 > const & phaseFraction,
                arraySlice2d< real64, USD3 > const & phaseCompFraction ) const
  {
    GEOS_UNUSED_VAR( pressure, temperature );
    // TODO: Constant values for now. To be linked with the static function call
    phaseFraction[m_liquidIndex] = 0.5;
    phaseFraction[m_vapourIndex] = 0.5;
    for( integer ic = 0; ic < m_numComponents; ++ic )
    {
      phaseCompFraction[m_liquidIndex][ic] = compFraction[ic];
      phaseCompFraction[m_vapourIndex][ic] = compFraction[ic];
    }
  }

  template< int USD1 >
  GEOS_HOST_DEVICE
  void compute( real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & compFraction,
                PhaseProp::SliceType const phaseFraction,
                PhaseComp::SliceType const phaseCompFraction ) const
  {
    GEOS_UNUSED_VAR( pressure, temperature );

    // TODO: Constant values for now. To be linked with the static function call
    phaseFraction.value[m_liquidIndex] = 0.5;
    phaseFraction.value[m_vapourIndex] = 0.5;
    for( integer ic = 0; ic < m_numComponents; ++ic )
    {
      phaseCompFraction.value[m_liquidIndex][ic] = compFraction[ic];
      phaseCompFraction.value[m_vapourIndex][ic] = compFraction[ic];
    }

    LvArray::forValuesInSlice( phaseFraction.derivs, setZero );
    LvArray::forValuesInSlice( phaseCompFraction.derivs, setZero );
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

  // Mark as a 2-phase flash
  static constexpr integer getNumberOfPhases() { return 2; }

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
