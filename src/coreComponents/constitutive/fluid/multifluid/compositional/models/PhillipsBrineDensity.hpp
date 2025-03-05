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
 * @file PhillipsBrineDensity.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_PHILLIPSBRINEDENSITY_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_PHILLIPSBRINEDENSITY_HPP_

#include "FunctionBase.hpp"

#include "constitutive/fluid/multifluid/Layouts.hpp"
#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"

#include "functions/TableFunction.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

class PhillipsBrineDensityUpdate final : public FunctionBaseUpdate
{
public:
  PhillipsBrineDensityUpdate( TableFunction const & brineDensityTable,
                              integer const waterIndex );

  template< int USD1, int USD2 >
  GEOS_HOST_DEVICE
  void compute( ComponentProperties::KernelWrapper const & componentProperties,
                real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< real64 const, USD1 > const & phaseComposition,
                real64 & molarDensity,
                arraySlice1d< real64, USD2 > const & dMolarDensity,
                real64 & massDensity,
                arraySlice1d< real64, USD2 > const & dMassDensity,
                bool useMass ) const;

  virtual void move( LvArray::MemorySpace const space, bool const touch ) override
  {
    FunctionBaseUpdate::move( space, touch );
    m_brineDensityTable.move( space, touch );
  }

protected:
  /// Table with brine density tabulated as a function (P,T,sal)
  TableFunction::KernelWrapper m_brineDensityTable;

  /// Index of the water component
  integer m_waterIndex;
};

class PhillipsBrineDensity : public FunctionBase
{
public:
  PhillipsBrineDensity( string const & name,
                        ComponentProperties const & componentProperties,
                        integer const phaseIndex,
                        ModelParameters const & modelParameters );

  static string catalogName() { return "PhillipsBrineDensity"; }

  virtual FunctionType functionType() const override
  {
    return FunctionType::DENSITY;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = PhillipsBrineDensityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const;

  // Create parameters unique to this model
  static std::unique_ptr< ModelParameters > createParameters( std::unique_ptr< ModelParameters > parameters );

  // Parameters for the Phillips model viscosity model
  class Parameters : public ModelParameters
  {
public:
    Parameters( std::unique_ptr< ModelParameters > parameters );
    ~Parameters() override = default;

    real64 m_waterCompressibility{4.5e-10};

private:
    void registerParametersImpl( MultiFluidBase * fluid ) override;
    void postInputInitializationImpl( MultiFluidBase const * fluid, ComponentProperties const & componentProperties ) override;

    struct viewKeyStruct
    {
      static constexpr char const * waterCompressibilityString() { return "waterCompressibility"; }
    };
  };

private:
  static void calculateBrineDensity( arraySlice1d< real64 const > const & pressureCoords,
                                     arraySlice1d< real64 const > const & temperatureCoords,
                                     real64 const & salinity,
                                     arraySlice1d< real64 > const & densities );

  static void calculatePureWaterDensity( arraySlice1d< real64 const > const & pressureCoords,
                                         arraySlice1d< real64 const > const & temperatureCoords,
                                         real64 const & compressibility,
                                         arraySlice1d< real64 > const & densities );

private:
  /// Table with brine density tabulated as a function of (P,T,sal)
  TableFunction const * m_brineDensityTable;

  /// Index of the water phase
  integer m_waterIndex;
};

template< int USD1, int USD2 >
GEOS_HOST_DEVICE
GEOS_FORCE_INLINE
void PhillipsBrineDensityUpdate::compute(
  ComponentProperties::KernelWrapper const & componentProperties,
  real64 const & pressure,
  real64 const & temperature,
  arraySlice1d< real64 const, USD1 > const & phaseComposition,
  real64 & molarDensity,
  arraySlice1d< real64, USD2 > const & dMolarDensity,
  real64 & massDensity,
  arraySlice1d< real64, USD2 > const & dMassDensity,
  bool useMass ) const
{
  //using Deriv = constitutive::multifluid::DerivativeOffset;
  GEOS_UNUSED_VAR( componentProperties );
  GEOS_UNUSED_VAR( pressure );
  GEOS_UNUSED_VAR( temperature );
  GEOS_UNUSED_VAR( phaseComposition );
  GEOS_UNUSED_VAR( molarDensity );
  GEOS_UNUSED_VAR( dMolarDensity );
  GEOS_UNUSED_VAR( massDensity );
  GEOS_UNUSED_VAR( dMassDensity );
  GEOS_UNUSED_VAR( useMass );
}

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_PHILLIPSBRINEDENSITY_HPP_
