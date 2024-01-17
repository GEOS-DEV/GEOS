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
 * @file ConstantViscosity.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_CONSTANTVISCOSITY_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_CONSTANTVISCOSITY_HPP_

#include "FunctionBase.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

class ConstantViscosityUpdate final : public FunctionBaseUpdate
{
public:
  ConstantViscosityUpdate() = default;

  GEOS_HOST_DEVICE
  void compute( ComponentProperties::KernelWrapper const & componentProperties,
                real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const > const & phaseComposition,
                arraySlice2d< real64 const > const & dPhaseComposition,
                real64 const & density,
                arraySlice1d< real64 const > const & dDensity,
                real64 & viscosity,
                arraySlice1d< real64 > const & dViscosity,
                bool useMass ) const;
};

class ConstantViscosity : public FunctionBase
{
public:
  ConstantViscosity( string const & name,
                     ComponentProperties const & componentProperties );

  static string catalogName() { return ""; }

  FunctionType functionType() const override
  {
    return FunctionType::VISCOSITY;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = ConstantViscosityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;
};

GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void ConstantViscosityUpdate::compute( ComponentProperties::KernelWrapper const & componentProperties,
                                       real64 const & pressure,
                                       real64 const & temperature,
                                       arraySlice1d< real64 const > const & phaseComposition,
                                       arraySlice2d< real64 const > const & dPhaseComposition,
                                       real64 const & density,
                                       arraySlice1d< real64 const > const & dDensity,
                                       real64 & viscosity,
                                       arraySlice1d< real64 > const & dViscosity,
                                       bool useMass ) const
{
  GEOS_UNUSED_VAR( componentProperties, pressure, temperature, useMass );
  GEOS_UNUSED_VAR( phaseComposition, dPhaseComposition );
  GEOS_UNUSED_VAR( density, dDensity );

  viscosity = 0.001;

  LvArray::forValuesInSlice( dViscosity, setZero );
}

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_CONSTANTVISCOSITY_HPP_
