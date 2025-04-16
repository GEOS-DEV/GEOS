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
 * @file PhillipsBrineViscosity.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_PHILLIPSBRINEVISCOSITY_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_PHILLIPSBRINEVISCOSITY_HPP_

#include "FunctionBase.hpp"

#include "constitutive/fluid/multifluid/Layouts.hpp"

#include "functions/TableFunction.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

class PhillipsBrineViscosityUpdate final : public FunctionBaseUpdate
{
public:
  explicit PhillipsBrineViscosityUpdate( TableFunction const & brineViscosityTable );

  template< integer USD1, integer USD2 >
  GEOS_HOST_DEVICE
  void compute( ComponentProperties::KernelWrapper const & componentProperties,
                real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & phaseComposition,
                real64 const & density,
                arraySlice1d< real64 const, USD2 > const & dDensity,
                real64 & viscosity,
                arraySlice1d< real64, USD2 > const & dViscosity,
                bool useMass ) const;

protected:
  /// Brine viscosity tabulated as a function (P,T)
  TableFunction::KernelWrapper m_brineViscosityTable;
};

class PhillipsBrineViscosity : public FunctionBase
{
public:
  PhillipsBrineViscosity( string const & name,
                          ComponentProperties const & componentProperties,
                          integer const phaseIndex,
                          ModelParameters const & modelParameters );

  static string catalogName() { return "PhillipsBrine"; }

  FunctionType functionType() const override
  {
    return FunctionType::VISCOSITY;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = PhillipsBrineViscosityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

  // Create parameters unique to this model
  static std::unique_ptr< ModelParameters > createParameters( std::unique_ptr< ModelParameters > parameters );

private:
  static TableFunction const * makeViscosityTable( string const & name,
                                                   ComponentProperties const & componentProperties,
                                                   ModelParameters const & modelParameters );

private:
  /// Brine viscosity tabulated as a function (P,T)
  TableFunction const * m_brineViscosityTable;
};

template< integer USD1, integer USD2 >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void PhillipsBrineViscosityUpdate::compute( ComponentProperties::KernelWrapper const & componentProperties,
                                            real64 const & pressure,
                                            real64 const & temperature,
                                            arraySlice1d< real64 const, USD1 > const & phaseComposition,
                                            real64 const & density,
                                            arraySlice1d< real64 const, USD2 > const & dDensity,
                                            real64 & viscosity,
                                            arraySlice1d< real64, USD2 > const & dViscosity,
                                            bool useMass ) const
{
  GEOS_UNUSED_VAR( componentProperties, pressure, useMass );
  GEOS_UNUSED_VAR( phaseComposition );
  GEOS_UNUSED_VAR( density, dDensity );

  using Deriv = constitutive::multifluid::DerivativeOffset;

  LvArray::forValuesInSlice( dViscosity, setZero );

  viscosity = m_brineViscosityTable.compute( &temperature, &dViscosity[Deriv::dT] );
}

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_PHILLIPSBRINEVISCOSITY_HPP_
