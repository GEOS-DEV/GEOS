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
 * @file CompositionalEnthalpy.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_COMPOSITIONALENTHALPY_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_COMPOSITIONALENTHALPY_HPP_

#include "FunctionBase.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/EquationOfState.hpp"

#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/CompositionalProperties.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/CubicEOSPhaseModel.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/SoreideWhitsonEOSPhaseModel.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

class CompositionalEnthalpyUpdate final : public FunctionBaseUpdate
{
public:
  explicit CompositionalEnthalpyUpdate( EquationOfStateType const equationOfState );

  template< integer USD1, integer USD2 >
  GEOS_HOST_DEVICE
  void compute( ComponentProperties::KernelWrapper const & componentProperties,
                real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & phaseComposition,
                real64 & enthalpy,
                arraySlice1d< real64, USD2 > const & dEnthalpy,
                bool useMass ) const;

private:
  EquationOfStateType const m_equationOfState;
};

class CompositionalEnthalpy : public FunctionBase
{
public:
  CompositionalEnthalpy( string const & name,
                         ComponentProperties const & componentProperties,
                         integer const phaseIndex,
                         ModelParameters const & modelParameters );

  static string catalogName() { return "CompositionalEnthalpy"; }

  virtual FunctionType functionType() const override
  {
    return FunctionType::ENTHALPY;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = CompositionalEnthalpyUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

  // Create parameters unique to this model
  static std::unique_ptr< ModelParameters > createParameters( std::unique_ptr< ModelParameters > parameters );

private:
  EquationOfStateType m_equationOfState;
};

template< integer USD1, integer USD2 >
GEOS_HOST_DEVICE
void CompositionalEnthalpyUpdate::compute(
  ComponentProperties::KernelWrapper const & componentProperties,
  real64 const & pressure,
  real64 const & temperature,
  arraySlice1d< real64 const, USD1 > const & phaseComposition,
  real64 & enthalpy,
  arraySlice1d< real64, USD2 > const & dEnthalpy,
  bool useMass ) const
{
  GEOS_UNUSED_VAR( useMass );
  GEOS_UNUSED_VAR( pressure );
  GEOS_UNUSED_VAR( componentProperties );
  GEOS_UNUSED_VAR( temperature );
  GEOS_UNUSED_VAR( phaseComposition );
  GEOS_UNUSED_VAR( enthalpy );
  GEOS_UNUSED_VAR( dEnthalpy );
  GEOS_UNUSED_VAR( pressure );
}

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_COMPOSITIONALENTHALPY_HPP_
