/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file KValueFlashModel.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_KVALUEFLASHMODEL_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_KVALUEFLASHMODEL_HPP_

#include "FunctionBase.hpp"

#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"

#include "functions/TableFunction.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

template< integer NUM_PHASE >
class KValueFlashParameters;

template< integer NUM_PHASE >
class KValueFlashModelUpdate final : public FunctionBaseUpdate
{
  static constexpr integer numPhases = NUM_PHASE;
  // Must be 2-phase or 3-phase
  static_assert( NUM_PHASE == 2 || NUM_PHASE == 3, "KValue flash must be 2-phase or 3-phase" );

public:
  using PhaseProp = MultiFluidVar< real64, 3, multifluid::LAYOUT_PHASE, multifluid::LAYOUT_PHASE_DC >;
  using PhaseComp = MultiFluidVar< real64, 4, multifluid::LAYOUT_PHASE_COMP, multifluid::LAYOUT_PHASE_COMP_DC >;

  KValueFlashModelUpdate( integer const numComponents,
                          TableFunction const & pressureTable,
                          TableFunction const & temperatureTable,
                          arrayView4d< real64 const > const & kValues );

  // Mark as a 2-phase or 3-phase flash
  GEOS_HOST_DEVICE
  static constexpr integer getNumberOfPhases() { return numPhases; }

  template< int USD1, int USD2 >
  GEOS_HOST_DEVICE
  void compute( ComponentProperties::KernelWrapper const & componentProperties,
                real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & compFraction,
                arraySlice2d< real64, USD2 > const & kValues,
                PhaseProp::SliceType const phaseFraction,
                PhaseComp::SliceType const phaseCompFraction ) const;

private:
  integer const m_numComponents{0};
  integer const m_numPressurePoints{0};
  integer const m_numTemperaturePoints{0};

  /// Table with pressure indices
  TableFunction::KernelWrapper m_pressureTable;

  /// Table with temperature indices
  TableFunction::KernelWrapper m_temperatureTable;

  /// Hypercube of k-values
  arrayView4d< real64 const > m_kValues;
};

template< integer NUM_PHASE >
class KValueFlashModel : public FunctionBase
{
public:
  KValueFlashModel( string const & name,
                    ComponentProperties const & componentProperties,
                    ModelParameters const & modelParameters );

  static string catalogName();

  FunctionType functionType() const override
  {
    return FunctionType::FLASH;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = KValueFlashModelUpdate< NUM_PHASE >;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

  // Create parameters unique to this model
  static std::unique_ptr< ModelParameters > createParameters( std::unique_ptr< ModelParameters > parameters );

private:
  KValueFlashParameters< NUM_PHASE > const * m_parameters{};
};

template< integer NUM_PHASE >
template< int USD1, int USD2 >
GEOS_HOST_DEVICE
void
KValueFlashModelUpdate< NUM_PHASE >::compute( ComponentProperties::KernelWrapper const & componentProperties,
                                              real64 const & pressure,
                                              real64 const & temperature,
                                              arraySlice1d< real64 const, USD1 > const & compFraction,
                                              arraySlice2d< real64, USD2 > const & kValues,
                                              PhaseProp::SliceType const phaseFraction,
                                              PhaseComp::SliceType const phaseCompFraction ) const
{
  GEOS_UNUSED_VAR( componentProperties );
  GEOS_UNUSED_VAR( pressure );
  GEOS_UNUSED_VAR( temperature );
  GEOS_UNUSED_VAR( compFraction );
  GEOS_UNUSED_VAR( kValues );
  GEOS_UNUSED_VAR( phaseFraction );
  GEOS_UNUSED_VAR( phaseCompFraction );
}

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_KVALUEFLASHMODEL_HPP_