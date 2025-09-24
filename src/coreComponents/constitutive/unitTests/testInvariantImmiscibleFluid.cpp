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


// File: testInvariantImmiscibleFluid.cpp
// Unit tests for InvariantImmiscibleFluid constitutive model
//
// Tests a three-phase (water/oil/gas) immiscible fluid model with constant properties:
// - Phase densities: 1000/800/100 kg/m³
// - Phase viscosities: 0.001/0.005/0.00002 Pa·s
// - Component molar weights: 0.018/0.2/0.016 kg/mol
//
// Verifies:
// - Phase property calculations (density, viscosity)
// - Mass vs molar basis conversions
// - Total density derivatives
// - KernelWrapper compute function

#include "FluidModelTest.hpp"
#include "constitutive/fluid/multifluid/constant/InvariantImmiscibleFluid.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "common/initializeEnvironment.hpp"

using namespace geos::constitutive;

namespace geos
{
namespace testing
{

using CompSlice = arraySlice1d< real64 const, compflow::USD_COMP - 1 >;
using PhasePropSlice = MultiFluidBase::PhaseProp::SliceType;
using PhaseCompSlice = MultiFluidBase::PhaseComp::SliceType;
using FluidPropSlice = MultiFluidBase::FluidProp::SliceType;

// Base test fixture template for validation tests with configurable component/phase counts
template< int NUM_COMPONENTS, int NUM_PHASES >
class InvariantImmiscibleFluidValidationTestFixture : public FluidModelTest< InvariantImmiscibleFluid, NUM_COMPONENTS, NUM_PHASES >
{
public:
  using Base = FluidModelTest< InvariantImmiscibleFluid, NUM_COMPONENTS, NUM_PHASES >;
  static constexpr real64 relTol = 1.0e-10;
  static constexpr real64 absTol = 1.0e-10;

  InvariantImmiscibleFluidValidationTestFixture()
  {
    Base::createFluid( "InvariantImmiscibleFluid", [this]( InvariantImmiscibleFluid & fluid ){
      setupFluidConfiguration( fluid );

      // create a dummy discretization with one quadrature point for
      // storing constitutive data

      conduit::Node node;
      dataRepository::Group rootGroup( "root", node );
      dataRepository::Group discretization( "discretization", &rootGroup );

      discretization.resize( 1 );   // one element
      fluid.allocateConstitutiveData( discretization, 1 );   // one quadrature point
      fluid.initializeState();
    } );
  }

private:
  void setupFluidConfiguration( InvariantImmiscibleFluid & fluid )
  {
    // Set up component names
    string_array & componentNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
    componentNames.resize( NUM_COMPONENTS );
    for( int i = 0; i < NUM_COMPONENTS; ++i )
    {
      componentNames[i] = "Component" + std::to_string( i );
    }

    // Set up phase names
    string_array & phaseNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );
    phaseNames.resize( NUM_PHASES );
    for( int i = 0; i < NUM_PHASES; ++i )
    {
      phaseNames[i] = "Phase" + std::to_string( i );
    }

    // Set up component molar weights
    array1d< real64 > & componentMolarWeight = fluid.getReference< array1d< real64 > >( MultiFluidBase::viewKeyStruct::componentMolarWeightString() );
    componentMolarWeight.resize( NUM_COMPONENTS );
    for( int i = 0; i < NUM_COMPONENTS; ++i )
    {
      componentMolarWeight[i] = 0.018 + i * 0.01; // Varying weights
    }

    // Set up constant densities (size must match numPhases)
    array1d< real64 > & densities = fluid.getReference< array1d< real64 > >( "densities" );
    densities.resize( NUM_PHASES );
    for( int i = 0; i < NUM_PHASES; ++i )
    {
      densities[i] = 1000.0 - i * 100.0; // Decreasing densities
    }

    // Set up constant viscosities (size must match numPhases)
    array1d< real64 > & viscosities = fluid.getReference< array1d< real64 > >( "viscosities" );
    viscosities.resize( NUM_PHASES );
    for( int i = 0; i < NUM_PHASES; ++i )
    {
      viscosities[i] = 0.001 + i * 0.001; // Increasing viscosities
    }

    fluid.setMassFlag( true ); // Use mass basis for simplicity
  }

public:

  // Helper method to get the fluid
  InvariantImmiscibleFluid & getFluid()
  {
    return static_cast< InvariantImmiscibleFluid & >( *Base::getFluid( InvariantImmiscibleFluid::catalogName() ) );
  }

};

template< typename MODEL_TYPE >
class InvariantImmiscibleFluidTestFixture : public FluidModelTest< InvariantImmiscibleFluid, 3, 3 >
{
public:
  using Base = FluidModelTest< InvariantImmiscibleFluid, 3, 3 >;
  static constexpr bool USE_MASS = std::tuple_element_t< 0, MODEL_TYPE >::value;

  static constexpr real64 relTol = 1.0e-10;
  static constexpr real64 absTol = 1.0e-10;

public:
  InvariantImmiscibleFluidTestFixture() = default;

  void SetUp() override
  {
    Base::createFluid( "InvariantImmiscibleFluid", [this]( InvariantImmiscibleFluid & fluid ){

      // Set up component names
      string_array & componentNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
      componentNames.resize( 3 );
      componentNames[0] = "H2O";
      componentNames[1] = "C5H12";
      componentNames[2] = "CH4";

      // Set up phase names
      string_array & phaseNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );
      phaseNames.resize( 3 );
      phaseNames[0] = "water";
      phaseNames[1] = "oil";
      phaseNames[2] = "gas";

      // Set up component molar weights (water, oil, gas)
      array1d< real64 > & componentMolarWeight = fluid.getReference< array1d< real64 > >( MultiFluidBase::viewKeyStruct::componentMolarWeightString() );
      componentMolarWeight.resize( 3 );
      componentMolarWeight[0] = 0.018;  // kg/mol
      componentMolarWeight[1] = 0.2;    // kg/mol
      componentMolarWeight[2] = 0.016;  // kg/mol

      // Set up constant densities
      array1d< real64 > & densities = fluid.getReference< array1d< real64 > >( "densities" );
      densities.resize( 3 );
      densities[0] = 1000.0;  // kg/m3
      densities[1] = 800.0;   // kg/m3
      densities[2] = 100.0;   // kg/m3

      // Set up constant viscosities
      array1d< real64 > & viscosities = fluid.getReference< array1d< real64 > >( "viscosities" );
      viscosities.resize( 3 );
      viscosities[0] = 0.001;     // Pa.s
      viscosities[1] = 0.005;     // Pa.s
      viscosities[2] = 0.00002;   // Pa.s

      fluid.setMassFlag( USE_MASS );

      // Need to allocate data for the fluid properties since we're now using KernelWrapper
      // create a dummy discretization with one quadrature point for
      // storing constitutive data

      conduit::Node node;
      dataRepository::Group rootGroup( "root", node );
      dataRepository::Group discretization( "discretization", &rootGroup );

      discretization.resize( 1 );   // one element
      fluid.allocateConstitutiveData( discretization, 1 );   // one quadrature point

      // Initialize the fluid by computing properties with a test composition
      StackArray< real64, 2, 3, compflow::LAYOUT_COMP > compositionData( 1, 3 );
      auto composition = compositionData[0];
      composition[0] = 0.4;
      composition[1] = 0.3;
      composition[2] = 0.3;
      CompSlice compositionSlice = composition.toSliceConst();
      fluid.initializeState();

      // Create a test point
      localIndex k = 0;
      localIndex q = 0;
      real64 pressure = 1.0e5;    // Example pressure (Pa)
      real64 temperature = 300.0; // Example temperature (K)

      // Update the fluid properties at this point
      this->getFluid().createKernelWrapper().update( k, q, pressure, temperature, compositionSlice );
    } );
  }

  // Add a helper method to get the fluid with const access
  const InvariantImmiscibleFluid & getFluid() const
  {
    return static_cast< const InvariantImmiscibleFluid & >(
      *const_cast< InvariantImmiscibleFluidTestFixture * >(this)->Base::getFluid( InvariantImmiscibleFluid::catalogName())
      );
  }

  // Add a helper method to get KernelWrapper for testing
  typename InvariantImmiscibleFluid::KernelWrapper getKernelWrapper() const
  {
    return this->getFluid().createKernelWrapper();
  }

};

using TestTypes = ::testing::Types<
  std::tuple< std::true_type >,  // Mass-based
  std::tuple< std::false_type >  // Mole-based
  >;

TYPED_TEST_SUITE( InvariantImmiscibleFluidTestFixture, TestTypes, );

TYPED_TEST( InvariantImmiscibleFluidTestFixture, PhaseProperties )
{
  auto kernelWrapper = this->getKernelWrapper();

  bool useMassQ = this->USE_MASS;

  // Constants for StackArray sizes
  constexpr integer numPhases = 3;
  constexpr integer numComponents = 3;
  constexpr integer numDerivs = numComponents + 2; // P, T, and compositions

  // Setup a test composition
  StackArray< real64, 2, numComponents, compflow::LAYOUT_COMP > compositionData( 1, numComponents );
  auto composition = compositionData[0];
  composition[0] = 0.4;
  composition[1] = 0.3;
  composition[2] = 0.3;

  // Use StackArray for CUDA compatibility
  StackArray< real64, 3, numPhases, multifluid::LAYOUT_PHASE > phaseFractionData( 1, 1, numPhases );
  StackArray< real64, 4, numPhases * numDerivs, multifluid::LAYOUT_PHASE_DC > dPhaseFractionData( 1, 1, numPhases, numDerivs );

  StackArray< real64, 3, numPhases, multifluid::LAYOUT_PHASE > phaseDensityData( 1, 1, numPhases );
  StackArray< real64, 4, numPhases * numDerivs, multifluid::LAYOUT_PHASE_DC > dPhaseDensityData( 1, 1, numPhases, numDerivs );

  StackArray< real64, 3, numPhases, multifluid::LAYOUT_PHASE > phaseMassDensityData( 1, 1, numPhases );
  StackArray< real64, 4, numPhases * numDerivs, multifluid::LAYOUT_PHASE_DC > dPhaseMassDensityData( 1, 1, numPhases, numDerivs );

  StackArray< real64, 3, numPhases, multifluid::LAYOUT_PHASE > phaseViscosityData( 1, 1, numPhases );
  StackArray< real64, 4, numPhases * numDerivs, multifluid::LAYOUT_PHASE_DC > dPhaseViscosityData( 1, 1, numPhases, numDerivs );

  StackArray< real64, 3, numPhases, multifluid::LAYOUT_PHASE > phaseEnthalpyData( 1, 1, numPhases );
  StackArray< real64, 4, numPhases * numDerivs, multifluid::LAYOUT_PHASE_DC > dPhaseEnthalpyData( 1, 1, numPhases, numDerivs );

  StackArray< real64, 3, numPhases, multifluid::LAYOUT_PHASE > phaseInternalEnergyData( 1, 1, numPhases );
  StackArray< real64, 4, numPhases * numDerivs, multifluid::LAYOUT_PHASE_DC > dPhaseInternalEnergyData( 1, 1, numPhases, numDerivs );

  StackArray< real64, 4, numPhases * numComponents, multifluid::LAYOUT_PHASE_COMP > phaseCompFractionData( 1, 1, numPhases, numComponents );
  StackArray< real64, 5, numPhases * numComponents * numDerivs, multifluid::LAYOUT_PHASE_COMP_DC > dPhaseCompFractionData( 1, 1, numPhases, numComponents, numDerivs );

  StackArray< real64, 2, 1, multifluid::LAYOUT_FLUID > totalDensityData( 1, 1 );
  StackArray< real64, 3, numDerivs, multifluid::LAYOUT_FLUID_DC > dtotalDensityData( 1, 1, numDerivs );

  // Get views from StackArray

  auto phaseFraction = phaseFractionData[0][0];
  auto dPhaseFraction = dPhaseFractionData[0][0];
  auto phaseDensity = phaseDensityData[0][0];
  auto dPhaseDensity = dPhaseDensityData[0][0];
  auto phaseMassDensity = phaseMassDensityData[0][0];
  auto dPhaseMassDensity = dPhaseMassDensityData[0][0];
  auto phaseViscosity = phaseViscosityData[0][0];
  auto dPhaseViscosity = dPhaseViscosityData[0][0];
  auto phaseEnthalpy = phaseEnthalpyData[0][0];
  auto dPhaseEnthalpy = dPhaseEnthalpyData[0][0];
  auto phaseInternalEnergy = phaseInternalEnergyData[0][0];
  auto dPhaseInternalEnergy = dPhaseInternalEnergyData[0][0];
  auto phaseCompFraction = phaseCompFractionData[0][0];
  auto dPhaseCompFraction = dPhaseCompFractionData[0][0];
  auto totalDensity = totalDensityData[0][0];
  auto dtotalDensity = dtotalDensityData[0][0];

  // Create slices for the compute call
  CompSlice compositionSlice = composition.toSliceConst();
  PhasePropSlice phaseFractionSlice = PhasePropSlice( phaseFraction, dPhaseFraction );
  PhasePropSlice phaseDensitySlice = PhasePropSlice( phaseDensity, dPhaseDensity );
  PhasePropSlice phaseMassDensitySlice = PhasePropSlice( phaseMassDensity, dPhaseMassDensity );
  PhasePropSlice phaseViscositySlice = PhasePropSlice( phaseViscosity, dPhaseViscosity );
  PhasePropSlice phaseEnthalpySlice = PhasePropSlice( phaseEnthalpy, dPhaseEnthalpy );
  PhasePropSlice phaseInternalEnergySlice = PhasePropSlice( phaseInternalEnergy, dPhaseInternalEnergy );
  PhaseCompSlice phaseCompFractionSlice = PhaseCompSlice( phaseCompFraction, dPhaseCompFraction );
  FluidPropSlice totalDensitySlice = FluidPropSlice( totalDensity, dtotalDensity );

  // Compute fluid properties
  kernelWrapper.compute( 1.0e5, 300.0,
                         compositionSlice,
                         phaseFractionSlice,
                         phaseDensitySlice,
                         phaseMassDensitySlice,
                         phaseViscositySlice,
                         phaseEnthalpySlice,
                         phaseInternalEnergySlice,
                         phaseCompFractionSlice,
                         totalDensitySlice );

  // Test phase densities using computed values
  for( integer ip = 0; ip < 3; ip++ )
  {
    real64 expectedMassDensity = (ip == 0) ? 1000.0 : (ip == 1) ? 800.0 : 100.0;
    real64 computedMassDensity = phaseMassDensity[ip];
    EXPECT_NEAR( computedMassDensity, expectedMassDensity, this->absTol );

    real64 computedDensity = phaseDensity[ip];
    if( useMassQ )
    {
      EXPECT_NEAR( computedDensity, expectedMassDensity, this->absTol );
    }
    else
    {
      real64 expectedDensity = (ip == 0) ? 1000.0/0.018 : (ip == 1) ? 800.0/0.2 : 100.0/0.016;
      EXPECT_NEAR( computedDensity, expectedDensity, this->absTol );
    }
  }

  // Test phase viscosities using computed values
  for( integer ip = 0; ip < 3; ip++ )
  {
    real64 expectedViscosity = (ip == 0) ? 0.001 : (ip == 1) ? 0.005 : 0.00002;
    EXPECT_NEAR( phaseViscosity[ip], expectedViscosity, this->absTol );
  }

}

TYPED_TEST( InvariantImmiscibleFluidTestFixture, WaterPhaseIndex )
{
  auto & fluid = this->getFluid();
  EXPECT_EQ( fluid.getWaterPhaseIndex(), 0 );
}

TYPED_TEST( InvariantImmiscibleFluidTestFixture, ComponentPhaseCountValidation )
{
  auto & fluid = this->getFluid();
  EXPECT_EQ( fluid.numFluidComponents(), 3 );
  EXPECT_EQ( fluid.numFluidPhases(), 3 );

  // Test the basic assumption that components = phases for this model
  EXPECT_EQ( fluid.numFluidComponents(), fluid.numFluidPhases() );
}

TYPED_TEST( InvariantImmiscibleFluidTestFixture, TotalDensityDerivative )
{
  auto kernelWrapper = this->getKernelWrapper();

  // Constants for StackArray sizes
  constexpr integer numPhases = 3;
  constexpr integer numComponents = 3;
  constexpr integer numDerivs = numComponents + 2; // P, T, and compositions

  // Setup a test composition
  StackArray< real64, 2, numComponents, compflow::LAYOUT_COMP > compositionData( 1, numComponents );
  auto composition = compositionData[0];
  composition[0] = 0.4;
  composition[1] = 0.3;
  composition[2] = 0.3;

  // Use StackArray for CUDA compatibility
  StackArray< real64, 3, numPhases, multifluid::LAYOUT_PHASE > phaseFractionData( 1, 1, numPhases );
  StackArray< real64, 4, numPhases * numDerivs, multifluid::LAYOUT_PHASE_DC > dPhaseFractionData( 1, 1, numPhases, numDerivs );

  StackArray< real64, 3, numPhases, multifluid::LAYOUT_PHASE > phaseDensityData( 1, 1, numPhases );
  StackArray< real64, 4, numPhases * numDerivs, multifluid::LAYOUT_PHASE_DC > dPhaseDensityData( 1, 1, numPhases, numDerivs );

  StackArray< real64, 3, numPhases, multifluid::LAYOUT_PHASE > phaseMassDensityData( 1, 1, numPhases );
  StackArray< real64, 4, numPhases * numDerivs, multifluid::LAYOUT_PHASE_DC > dPhaseMassDensityData( 1, 1, numPhases, numDerivs );

  StackArray< real64, 3, numPhases, multifluid::LAYOUT_PHASE > phaseViscosityData( 1, 1, numPhases );
  StackArray< real64, 4, numPhases * numDerivs, multifluid::LAYOUT_PHASE_DC > dPhaseViscosityData( 1, 1, numPhases, numDerivs );

  StackArray< real64, 3, numPhases, multifluid::LAYOUT_PHASE > phaseEnthalpyData( 1, 1, numPhases );
  StackArray< real64, 4, numPhases * numDerivs, multifluid::LAYOUT_PHASE_DC > dPhaseEnthalpyData( 1, 1, numPhases, numDerivs );

  StackArray< real64, 3, numPhases, multifluid::LAYOUT_PHASE > phaseInternalEnergyData( 1, 1, numPhases );
  StackArray< real64, 4, numPhases * numDerivs, multifluid::LAYOUT_PHASE_DC > dPhaseInternalEnergyData( 1, 1, numPhases, numDerivs );

  StackArray< real64, 4, numPhases * numComponents, multifluid::LAYOUT_PHASE_COMP > phaseCompFractionData( 1, 1, numPhases, numComponents );
  StackArray< real64, 5, numPhases * numComponents * numDerivs, multifluid::LAYOUT_PHASE_COMP_DC > dPhaseCompFractionData( 1, 1, numPhases, numComponents, numDerivs );

  StackArray< real64, 2, 1, multifluid::LAYOUT_FLUID > totalDensityData( 1, 1 );
  StackArray< real64, 3, numDerivs, multifluid::LAYOUT_FLUID_DC > dtotalDensityData( 1, 1, numDerivs );

  // Get views from StackArray
  auto phaseFraction = phaseFractionData[0][0];
  auto dPhaseFraction = dPhaseFractionData[0][0];
  auto phaseDensity = phaseDensityData[0][0];
  auto dPhaseDensity = dPhaseDensityData[0][0];
  auto phaseMassDensity = phaseMassDensityData[0][0];
  auto dPhaseMassDensity = dPhaseMassDensityData[0][0];
  auto phaseViscosity = phaseViscosityData[0][0];
  auto dPhaseViscosity = dPhaseViscosityData[0][0];
  auto phaseEnthalpy = phaseEnthalpyData[0][0];
  auto dPhaseEnthalpy = dPhaseEnthalpyData[0][0];
  auto phaseInternalEnergy = phaseInternalEnergyData[0][0];
  auto dPhaseInternalEnergy = dPhaseInternalEnergyData[0][0];
  auto phaseCompFraction = phaseCompFractionData[0][0];
  auto dPhaseCompFraction = dPhaseCompFractionData[0][0];
  auto totalDensity = totalDensityData[0][0];
  auto dtotalDensity = dtotalDensityData[0][0];

  // Create slices for the compute call
  CompSlice compositionSlice = composition.toSliceConst();
  PhasePropSlice phaseFractionSlice = PhasePropSlice( phaseFraction, dPhaseFraction );
  PhasePropSlice phaseDensitySlice = PhasePropSlice( phaseDensity, dPhaseDensity );
  PhasePropSlice phaseMassDensitySlice = PhasePropSlice( phaseMassDensity, dPhaseMassDensity );
  PhasePropSlice phaseViscositySlice = PhasePropSlice( phaseViscosity, dPhaseViscosity );
  PhasePropSlice phaseEnthalpySlice = PhasePropSlice( phaseEnthalpy, dPhaseEnthalpy );
  PhasePropSlice phaseInternalEnergySlice = PhasePropSlice( phaseInternalEnergy, dPhaseInternalEnergy );
  PhaseCompSlice phaseCompFractionSlice = PhaseCompSlice( phaseCompFraction, dPhaseCompFraction );
  FluidPropSlice totalDensitySlice = FluidPropSlice( totalDensity, dtotalDensity );

  kernelWrapper.compute( 1.0e5, 300.0,
                         compositionSlice,
                         phaseFractionSlice,
                         phaseDensitySlice,
                         phaseMassDensitySlice,
                         phaseViscositySlice,
                         phaseEnthalpySlice,
                         phaseInternalEnergySlice,
                         phaseCompFractionSlice,
                         totalDensitySlice );

  for( integer ip = 0; ip < numPhases; ++ip )
  {
    EXPECT_NEAR( composition[ip], phaseFraction[ip], this->absTol );
  }

  // The derivative of total density w.r.t phase fraction for phase ip is phase density for that phase
  real64 inverseTotalDensity = 0.0;
  for( integer ip = 0; ip < numPhases; ++ip )
  {
    if( phaseFraction[ip] > 0.0 )
    {
      inverseTotalDensity += phaseFraction[ip] / phaseDensity[ip];
    }
  }
  real64 expectedTotalDensity = 1.0 / inverseTotalDensity;

  for( integer ic = 0; ic < 3; ++ic ) // 3 component derivatives
  {
    real64 dInverseDensity_dC = 0.0;
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      if( phaseFraction[ip] > 0.0 )
      {
        real64 densInv = 1.0 / phaseDensity[ip];
        real64 value = phaseFraction[ip] * densInv;
        real64 dPhi_dC = dPhaseFraction[ip][2 + ic];
        real64 dRho_dC = dPhaseDensity[ip][2 + ic];

        dInverseDensity_dC += (dPhi_dC - value * dRho_dC) * densInv;
      }
    }
    real64 expectedDeriv = -expectedTotalDensity * expectedTotalDensity * dInverseDensity_dC;
    EXPECT_NEAR( dtotalDensity[2 + ic], expectedDeriv, this->absTol );
  }
}

} // namespace testing
} // namespace geos
