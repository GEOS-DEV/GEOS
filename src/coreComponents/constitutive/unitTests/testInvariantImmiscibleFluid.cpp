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
 * @file testInvariantImmiscibleFluid.cpp
 * @brief Unit tests for InvariantImmiscibleFluid constitutive model
 *
 * Tests a three-phase (water/oil/gas) immiscible fluid model with constant properties:
 * - Phase densities: 1000/800/100 kg/m³
 * - Phase viscosities: 0.001/0.005/0.00002 Pa·s
 * - Component molar weights: 0.018/0.2/0.016 kg/mol
 *
 * Verifies:
 * - Phase property calculations (density, viscosity)
 * - Mass vs molar basis conversions
 * - Total density derivatives
 * - KernelWrapper compute function
 */

#include "FluidModelTest.hpp"
#include "constitutive/fluid/multifluid/constant/InvariantImmiscibleFluid.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "common/initializeEnvironment.hpp"

using namespace geos::constitutive;

namespace geos
{
namespace testing
{

// Test helper class to access protected members of InvariantImmiscibleFluid
class TestableInvariantImmiscibleFluid : public InvariantImmiscibleFluid
{
public:
  TestableInvariantImmiscibleFluid( string const & name, Group * parent )
    : InvariantImmiscibleFluid( name, parent ) {}

  // Expose the protected method for testing
  void testPostInputInitialization()
  {
    postInputInitialization();
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
  InvariantImmiscibleFluidTestFixture()
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
      fluid.initialize( 1, 1 );

      // Initialize the fluid by computing properties with a test composition
      array1d< real64 > composition;
      composition.resize( 3 );
      composition[0] = 0.4;
      composition[1] = 0.3;
      composition[2] = 0.3;
      fluid.initializeState();

      // Create a test point
      localIndex k = 0;
      localIndex q = 0;
      real64 pressure = 1.0e5;    // Example pressure (Pa)
      real64 temperature = 300.0; // Example temperature (K)

      // Update the fluid properties at this point
      this->getFluid().createKernelWrapper().update( k, q, pressure, temperature, composition );
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

  // Test compute directly using KernelWrapper
  void testCompute( real64 const pressure,
                    real64 const temperature,
                    array1d< real64 > const & composition ) const
  {
    // Get the wrapper
    auto kernelWrapper = getKernelWrapper();

    // Create local slices for results
    integer numPhases = 3;
    integer numComponents = 3;
    integer numDerivs = numComponents + 2; // P, T, and compositions

    // Create arrays for the values and derivatives
    array1d< real64 > phaseFractionValues( numPhases );
    array2d< real64 > phaseFractionDerivs( numPhases, numDerivs );

    array1d< real64 > phaseDensityValues( numPhases );
    array2d< real64 > phaseDensityDerivs( numPhases, numDerivs );

    array1d< real64 > phaseMassDensityValues( numPhases );
    array2d< real64 > phaseMassDensityDerivs( numPhases, numDerivs );

    array1d< real64 > phaseViscosityValues( numPhases );
    array2d< real64 > phaseViscosityDerivs( numPhases, numDerivs );

    array1d< real64 > phaseEnthalpyValues( numPhases );
    array2d< real64 > phaseEnthalpyDerivs( numPhases, numDerivs );

    array1d< real64 > phaseInternalEnergyValues( numPhases );
    array2d< real64 > phaseInternalEnergyDerivs( numPhases, numDerivs );

    array2d< real64 > phaseCompFractionValues( numPhases, numComponents );
    array3d< real64 > phaseCompFractionDerivs( numPhases, numComponents, numDerivs );

    real64 totalDensityValue;
    array1d< real64 > totalDensityDerivs( numDerivs );

    // Create slices from these arrays
    MultiFluidBase::PhaseProp::SliceType phaseFractionSlice{phaseFractionValues, phaseFractionDerivs};
    MultiFluidBase::PhaseProp::SliceType phaseDensitySlice{phaseDensityValues, phaseDensityDerivs};
    MultiFluidBase::PhaseProp::SliceType phaseMassDensitySlice{phaseMassDensityValues, phaseMassDensityDerivs};
    MultiFluidBase::PhaseProp::SliceType phaseViscositySlice{phaseViscosityValues, phaseViscosityDerivs};
    MultiFluidBase::PhaseProp::SliceType phaseEnthalpySlice{phaseEnthalpyValues, phaseEnthalpyDerivs};
    MultiFluidBase::PhaseProp::SliceType phaseInternalEnergySlice{phaseInternalEnergyValues, phaseInternalEnergyDerivs};
    MultiFluidBase::PhaseComp::SliceType phaseCompFractionSlice{phaseCompFractionValues, phaseCompFractionDerivs};
    MultiFluidBase::FluidProp::SliceType totalDensitySlice{totalDensityValue, totalDensityDerivs};

    // Call compute directly
    kernelWrapper.compute( pressure,
                           temperature,
                           composition,
                           phaseFractionSlice,
                           phaseDensitySlice,
                           phaseMassDensitySlice,
                           phaseViscositySlice,
                           phaseEnthalpySlice,
                           phaseInternalEnergySlice,
                           phaseCompFractionSlice,
                           totalDensitySlice );

    // Return the results for testing
    return;
  }
};

using TestTypes = ::testing::Types<
  std::tuple< std::true_type >,  // Mass-based
  std::tuple< std::false_type >  // Mole-based
  >;

TYPED_TEST_SUITE( InvariantImmiscibleFluidTestFixture, TestTypes, );

TYPED_TEST( InvariantImmiscibleFluidTestFixture, PhaseProperties )
{
  auto & fluid = this->getFluid();

  bool useMassQ = this->USE_MASS;

  // Test phase densities - expectations depend on mass vs molar basis
  for( integer ip = 0; ip < 3; ip++ )
  {
    real64 expectedMassDensity = (ip == 0) ? 1000.0 : (ip == 1) ? 800.0 : 100.0;
    real64 computedMassDensity = fluid.phaseMassDensity().operator()( 0, 0, ip );
    EXPECT_NEAR( computedMassDensity, expectedMassDensity, this->absTol );
    real64 computedDensity = fluid.phaseDensity().operator()( 0, 0, ip );
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

  // Test phase viscosities - same regardless of mass/molar basis
  for( integer ip = 0; ip < 3; ip++ )
  {
    real64 expectedViscosity = (ip == 0) ? 0.001 : (ip == 1) ? 0.005 : 0.00002;
    EXPECT_NEAR( fluid.phaseViscosity().operator()( 0, 0, ip ), expectedViscosity, this->absTol );
  }

  // Test KernelWrapper compute functionality
  array1d< real64 > testComp;
  testComp.resize( 3 );
  testComp[0] = 0.4;
  testComp[1] = 0.3;
  testComp[2] = 0.3;
  this->testCompute( 1.0e5, 300.0, testComp );
}

TYPED_TEST( InvariantImmiscibleFluidTestFixture, WaterPhaseIndex )
{
  auto & fluid = this->getFluid();
  EXPECT_EQ( fluid.getWaterPhaseIndex(), 0 );
}

TYPED_TEST( InvariantImmiscibleFluidTestFixture, ComponentPhaseCountValidation )
{
  // Test that the fluid correctly validates that number of components equals number of phases
  
  // This should work fine - our fixture sets up 3 components and 3 phases
  auto & fluid = this->getFluid();
  EXPECT_EQ( fluid.numFluidComponents(), 3 );
  EXPECT_EQ( fluid.numFluidPhases(), 3 );
  
  // Test the basic assumption that components = phases for this model
  EXPECT_EQ( fluid.numFluidComponents(), fluid.numFluidPhases() );
  
  // Test case 1: More components than phases (2 components, 1 phase)
  {
    TestableInvariantImmiscibleFluid testFluid1( "testFluid1", nullptr );
    
    // Set up components
    string_array & componentNames = testFluid1.getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
    componentNames.resize( 2 );
    componentNames[0] = "Component0";
    componentNames[1] = "Component1";
    
    // Set up phases
    string_array & phaseNames = testFluid1.getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );
    phaseNames.resize( 1 );
    phaseNames[0] = "Phase0";
    
    // Set up component molar weights
    array1d< real64 > & componentMolarWeight = testFluid1.getReference< array1d< real64 > >( MultiFluidBase::viewKeyStruct::componentMolarWeightString() );
    componentMolarWeight.resize( 2 );
    componentMolarWeight[0] = 0.018;
    componentMolarWeight[1] = 0.028;
    
    // Set up densities (size must match numPhases)
    array1d< real64 > & densities = testFluid1.getReference< array1d< real64 > >( "densities" );
    densities.resize( 1 );
    densities[0] = 1000.0;
    
    // Set up viscosities (size must match numPhases)
    array1d< real64 > & viscosities = testFluid1.getReference< array1d< real64 > >( "viscosities" );
    viscosities.resize( 1 );
    viscosities[0] = 0.001;
    
    EXPECT_THROW( testFluid1.testPostInputInitialization(), InputError );
  }
  
  // Test case 2: More phases than components (2 components, 3 phases)
  {
    TestableInvariantImmiscibleFluid testFluid2( "testFluid2", nullptr );
    
    // Set up components
    string_array & componentNames = testFluid2.getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
    componentNames.resize( 2 );
    componentNames[0] = "Component0";
    componentNames[1] = "Component1";
    
    // Set up phases
    string_array & phaseNames = testFluid2.getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );
    phaseNames.resize( 3 );
    phaseNames[0] = "Phase0";
    phaseNames[1] = "Phase1";
    phaseNames[2] = "Phase2";
    
    // Set up component molar weights
    array1d< real64 > & componentMolarWeight = testFluid2.getReference< array1d< real64 > >( MultiFluidBase::viewKeyStruct::componentMolarWeightString() );
    componentMolarWeight.resize( 2 );
    componentMolarWeight[0] = 0.018;
    componentMolarWeight[1] = 0.028;
    
    // Set up densities (size must match numPhases)
    array1d< real64 > & densities = testFluid2.getReference< array1d< real64 > >( "densities" );
    densities.resize( 3 );
    densities[0] = 1000.0;
    densities[1] = 900.0;
    densities[2] = 800.0;
    
    // Set up viscosities (size must match numPhases)
    array1d< real64 > & viscosities = testFluid2.getReference< array1d< real64 > >( "viscosities" );
    viscosities.resize( 3 );
    viscosities[0] = 0.001;
    viscosities[1] = 0.002;
    viscosities[2] = 0.003;
    
    EXPECT_THROW( testFluid2.testPostInputInitialization(), InputError );
  }
  
  // Test case 3: Matching counts should not throw (3 components, 3 phases)
  {
    TestableInvariantImmiscibleFluid testFluid3( "testFluid3", nullptr );
    
    // Set up components
    string_array & componentNames = testFluid3.getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
    componentNames.resize( 3 );
    componentNames[0] = "Component0";
    componentNames[1] = "Component1";
    componentNames[2] = "Component2";
    
    // Set up phases
    string_array & phaseNames = testFluid3.getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );
    phaseNames.resize( 3 );
    phaseNames[0] = "Phase0";
    phaseNames[1] = "Phase1";
    phaseNames[2] = "Phase2";
    
    // Set up component molar weights
    array1d< real64 > & componentMolarWeight = testFluid3.getReference< array1d< real64 > >( MultiFluidBase::viewKeyStruct::componentMolarWeightString() );
    componentMolarWeight.resize( 3 );
    componentMolarWeight[0] = 0.018;
    componentMolarWeight[1] = 0.028;
    componentMolarWeight[2] = 0.016;
    
    // Set up densities (size must match numPhases)
    array1d< real64 > & densities = testFluid3.getReference< array1d< real64 > >( "densities" );
    densities.resize( 3 );
    densities[0] = 1000.0;
    densities[1] = 900.0;
    densities[2] = 800.0;
    
    // Set up viscosities (size must match numPhases)
    array1d< real64 > & viscosities = testFluid3.getReference< array1d< real64 > >( "viscosities" );
    viscosities.resize( 3 );
    viscosities[0] = 0.001;
    viscosities[1] = 0.002;
    viscosities[2] = 0.003;
    
    EXPECT_NO_THROW( testFluid3.testPostInputInitialization() );
  }
}

TYPED_TEST( InvariantImmiscibleFluidTestFixture, TotalDensityDerivative )
{
  auto kernelWrapper = this->getKernelWrapper();

  // Setup a test composition
  array1d< real64 > composition( 3 );
  composition[0] = 0.4;
  composition[1] = 0.3;
  composition[2] = 0.3;

  // Prepare slices for results
  integer numPhases = 3;
  integer numDerivs = 5; // 3 components + P + T
  array1d< real64 > phaseFractionValues( numPhases );
  array2d< real64 > phaseFractionDerivs( numPhases, numDerivs );
  array1d< real64 > phaseDensityValues( numPhases );
  array2d< real64 > phaseDensityDerivs( numPhases, numDerivs );
  array1d< real64 > phaseMassDensityValues( numPhases );
  array2d< real64 > phaseMassDensityDerivs( numPhases, numDerivs );
  array1d< real64 > phaseViscosityValues( numPhases );
  array2d< real64 > phaseViscosityDerivs( numPhases, numDerivs );
  array1d< real64 > phaseEnthalpyValues( numPhases );
  array2d< real64 > phaseEnthalpyDerivs( numPhases, numDerivs );
  array1d< real64 > phaseInternalEnergyValues( numPhases );
  array2d< real64 > phaseInternalEnergyDerivs( numPhases, numDerivs );
  array2d< real64 > phaseCompFractionValues( numPhases, numPhases );
  array3d< real64 > phaseCompFractionDerivs( numPhases, numPhases, numDerivs );
  real64 totalDensityValue = 0.0;
  array1d< real64 > totalDensityDerivs( numDerivs );

  MultiFluidBase::PhaseProp::SliceType phaseFractionSlice{phaseFractionValues, phaseFractionDerivs};
  MultiFluidBase::PhaseProp::SliceType phaseDensitySlice{phaseDensityValues, phaseDensityDerivs};
  MultiFluidBase::PhaseProp::SliceType phaseMassDensitySlice{phaseMassDensityValues, phaseMassDensityDerivs};
  MultiFluidBase::PhaseProp::SliceType phaseViscositySlice{phaseViscosityValues, phaseViscosityDerivs};
  MultiFluidBase::PhaseProp::SliceType phaseEnthalpySlice{phaseEnthalpyValues, phaseEnthalpyDerivs};
  MultiFluidBase::PhaseProp::SliceType phaseInternalEnergySlice{phaseInternalEnergyValues, phaseInternalEnergyDerivs};
  MultiFluidBase::PhaseComp::SliceType phaseCompFractionSlice{phaseCompFractionValues, phaseCompFractionDerivs};
  MultiFluidBase::FluidProp::SliceType totalDensitySlice{totalDensityValue, totalDensityDerivs};

  // Call compute
  kernelWrapper.compute( 1.0e5, 300.0, composition,
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
    EXPECT_NEAR( composition[ip], phaseFractionValues[ip], this->absTol );
  }

  // The derivative of total density w.r.t phase fraction for phase ip is phase density for that phase
  real64 inverseTotalDensity = 0.0;
  for( integer ip = 0; ip < numPhases; ++ip )
  {
    if( phaseFractionValues[ip] > 0.0 )
    {
      inverseTotalDensity += phaseFractionValues[ip] / phaseDensityValues[ip];
    }
  }
  real64 expectedTotalDensity = 1.0 / inverseTotalDensity;

  for( integer ic = 0; ic < 3; ++ic ) // 3 component derivatives
  {
    real64 dInverseDensity_dC = 0.0;
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      if( phaseFractionValues[ip] > 0.0 )
      {
        real64 densInv = 1.0 / phaseDensityValues[ip];
        real64 value = phaseFractionValues[ip] * densInv;
        real64 dPhi_dC = phaseFractionDerivs[ip][2 + ic];
        real64 dRho_dC = phaseDensityDerivs[ip][2 + ic];

        dInverseDensity_dC += (dPhi_dC - value * dRho_dC) * densInv;
      }
    }
    real64 expectedDeriv = -expectedTotalDensity * expectedTotalDensity * dInverseDensity_dC;
    EXPECT_NEAR( totalDensityDerivs[2 + ic], expectedDeriv, this->absTol );
  }
}

} // namespace testing
} // namespace geos
