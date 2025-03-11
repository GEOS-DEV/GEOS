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
#include "CompositionalDensity.hpp"

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

int constexpr USD1 = 0;
int constexpr USD2 = 0;

class PhillipsBrineDensityUpdate final : public FunctionBaseUpdate
{
public:
  PhillipsBrineDensityUpdate( TableFunction const & brineVolumeShiftTable,
                              integer const waterIndex,
                              real64 const brineMolarWeight,
                              EquationOfStateType const equationOfState );

  //template< int USD1, int USD2 >
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
    m_volumeShiftTable.move( space, touch );
  }

protected:
  /// Volume shift of the water component tabulated as a function (P,T)
  TableFunction::KernelWrapper m_volumeShiftTable;

  /// Index of the water component
  integer const m_waterIndex;

  /// The brine molecular weight
  real64 const m_brineMolarWeight;

  /// Equation of state for the density correction
  EquationOfStateType const m_equationOfState;
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
    real64 m_saltMolarWeight{58.44e-3};

private:
    void registerParametersImpl( MultiFluidBase * fluid ) override;
    void postInputInitializationImpl( MultiFluidBase const * fluid, ComponentProperties const & componentProperties ) override;

    struct viewKeyStruct
    {
      static constexpr char const * waterCompressibilityString() { return "waterCompressibility"; }
      static constexpr char const * saltMolarWeightString() { return "saltMolarWeight"; }
    };
  };

private:
  static real64 calculateBrineMolarWeight( real64 const & waterMolarWeight,
                                           real64 const & saltMolarWeight,
                                           real64 const & salinity );

  static void calculateBrineDensity( arraySlice1d< real64 const > const & pressureCoords,
                                     arraySlice1d< real64 const > const & temperatureCoords,
                                     real64 const & salinity,
                                     arraySlice1d< real64 > const & densities );

  static void calculatePureWaterDensity( arraySlice1d< real64 const > const & pressureCoords,
                                         arraySlice1d< real64 const > const & temperatureCoords,
                                         real64 const & compressibility,
                                         arraySlice1d< real64 > const & densities );

  static void calculateEosWaterMolarVolume( arraySlice1d< real64 const > const & pressureCoords,
                                            arraySlice1d< real64 const > const & temperatureCoords,
                                            ComponentProperties const & componentProperties,
                                            EquationOfStateType const equationOfState,
                                            integer const waterIndex,
                                            arraySlice1d< real64 > const & molarVolume );

  static TableFunction const * makeVolumeShiftTable( string const & name,
                                                     ComponentProperties const & componentProperties,
                                                     ModelParameters const & modelParameters,
                                                     EquationOfStateType const equationOfState,
                                                     real64 const brineMolarWeight,
                                                     integer const waterIndex );

private:
  /// Volume shift of the water component tabulated as a function (P,T)
  TableFunction const * m_volumeShiftTable;

  /// Index of the water phase
  integer m_waterIndex;

  /// Equation of state for the density correction
  EquationOfStateType m_equationOfState;

  /// The brine molecular weight
  real64 m_brineMolarWeight;
};

/**
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
   GEOS_UNUSED_VAR( molarDensity );
   GEOS_UNUSED_VAR( dMolarDensity );
   GEOS_UNUSED_VAR( massDensity );
   GEOS_UNUSED_VAR( dMassDensity );
   GEOS_UNUSED_VAR( useMass );

   integer const numComps = componentProperties.m_componentMolarWeight.size();
   integer const numDofs = 2 + numComps;

   real64 const input[2] = { pressure, temperature };
   real64 brineDensityDeriv[2]{};
   real64 const brineMassDensity = m_volumeShiftTable.compute( input, brineDensityDeriv );
   real64 const brineMolarVolume = m_brineMolarWeight / brineMassDensity;

   real64 compressibilityFactor = 0.0;
   stackArray1d< real64, 2+MultiFluidConstants::MAX_NUM_COMPONENTS > tempDerivs( numDofs );
   stackArray1d< real64, MultiFluidConstants::MAX_NUM_COMPONENTS > waterComposition( numComps );
   real64 const x_h2o = phaseComposition[m_waterIndex];
   for (integer ic = 0; ic < numComps; ++ic)
   {
    waterComposition[ic] = 0.0;
   }
   waterComposition[m_waterIndex] = 1.0;

   std::cout << "COMP " << phaseComposition << std::endl;

   CompositionalDensityUpdate::computeCompressibilityFactor( numComps,
                                                            pressure,
                                                            temperature,
                                                            phaseComposition,
                                                            componentProperties,
                                                            m_equationOfState,
                                                            compressibilityFactor,
                                                            tempDerivs.toSlice() );

   real64 mixtureMolarVolume = constants::gasConstant * temperature * compressibilityFactor / pressure;

   CompositionalDensityUpdate::computeCompressibilityFactor( numComps,
                                                            pressure,
                                                            temperature,
                                                            waterComposition.toSliceConst(),
                                                            componentProperties,
                                                            m_equationOfState,
                                                            compressibilityFactor,
                                                            tempDerivs.toSlice() );

   real64 waterMolarVolume = constants::gasConstant * temperature * compressibilityFactor / pressure;

   std::cout << "VOLUME " << brineMolarVolume << " " << waterMolarVolume  << " " << x_h2o*waterMolarVolume  << " " << mixtureMolarVolume <<
      "\n";
   }
 */

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_MODELS_PHILLIPSBRINEDENSITY_HPP_
