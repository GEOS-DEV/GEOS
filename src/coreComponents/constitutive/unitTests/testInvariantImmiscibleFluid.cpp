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

#include "FluidModelTest.hpp"
#include "constitutive/fluid/multifluid/constant/InvariantImmiscibleFluid.hpp"
#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "common/initializeEnvironment.hpp"

using namespace geos::constitutive;

namespace geos
{
namespace testing
{

template< typename MODEL_TYPE >
class InvariantImmiscibleFluidTestFixture : public FluidModelTest< InvariantImmiscibleFluid, 3, 3 >
{
public:
  using Base = FluidModelTest< InvariantImmiscibleFluid, 3, 3 >;
  static constexpr bool USE_MASS = std::tuple_element_t< 0, MODEL_TYPE >::value;

  static constexpr real64 relTol = 1.0e-15;
  static constexpr real64 absTol = 1.0e-15;

public:
  InvariantImmiscibleFluidTestFixture()
  {
    Base::createFluid( "InvariantImmiscibleFluid", []( InvariantImmiscibleFluid & fluid ){
      
      // Set up phase names
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
      fluid.update( k, q, pressure, temperature, composition );
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

  // Test using direct property access instead of getFunctions
  auto kernelWrapper = this->getKernelWrapper();

  // Test phase densities
  for( integer ip = 0; ip < 3; ip++ )
  {
    real64 expectedDensity = (ip == 0) ? 1000.0 : (ip == 1) ? 800.0 : 100.0;
    EXPECT_NEAR( fluid.phaseDensity().operator()( 0, 0, ip ), expectedDensity, this->absTol );
  }
  // Test phase viscosities
  for( integer ip = 0; ip < 3; ip++ )
  {
    real64 expectedViscosity = (ip == 0) ? 0.001 : (ip == 1) ? 0.005 : 0.00002;
    EXPECT_NEAR( fluid.phaseViscosity().operator()( 0, 0, ip ), expectedViscosity, this->absTol );
  }

  // Test direct compute with KernelWrapper - proper array initialization
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

//TYPED_TEST( InvariantImmiscibleFluidTestFixture, TotalDensityDerivative )
//{
//  auto & fluid = this->getFluid();
//  auto kernelWrapper = this->getKernelWrapper();
//
//  // Setup a test composition
//  array1d< real64 > composition( 3 );
//  composition[0] = 0.4;
//  composition[1] = 0.3;
//  composition[2] = 0.3;
//
//  // Prepare slices for results
//  integer numPhases = 3;
//  integer numDerivs = 5; // 3 compositions + P + T
//  array1d< real64 > phaseFractionValues( numPhases );
//  array2d< real64 > phaseFractionDerivs( numPhases, numDerivs );
//  array1d< real64 > phaseDensityValues( numPhases );
//  array2d< real64 > phaseDensityDerivs( numPhases, numDerivs );
//  array1d< real64 > phaseMassDensityValues( numPhases );
//  array2d< real64 > phaseMassDensityDerivs( numPhases, numDerivs );
//  array1d< real64 > phaseViscosityValues( numPhases );
//  array2d< real64 > phaseViscosityDerivs( numPhases, numDerivs );
//  array1d< real64 > phaseEnthalpyValues( numPhases );
//  array2d< real64 > phaseEnthalpyDerivs( numPhases, numDerivs );
//  array1d< real64 > phaseInternalEnergyValues( numPhases );
//  array2d< real64 > phaseInternalEnergyDerivs( numPhases, numDerivs );
//  array2d< real64 > phaseCompFractionValues( numPhases, numPhases );
//  array3d< real64 > phaseCompFractionDerivs( numPhases, numPhases, numDerivs );
//  real64 totalDensityValue = 0.0;
//  array1d< real64 > totalDensityDerivs( numDerivs );
//
//  MultiFluidBase::PhaseProp::SliceType phaseFractionSlice{phaseFractionValues, phaseFractionDerivs};
//  MultiFluidBase::PhaseProp::SliceType phaseDensitySlice{phaseDensityValues, phaseDensityDerivs};
//  MultiFluidBase::PhaseProp::SliceType phaseMassDensitySlice{phaseMassDensityValues, phaseMassDensityDerivs};
//  MultiFluidBase::PhaseProp::SliceType phaseViscositySlice{phaseViscosityValues, phaseViscosityDerivs};
//  MultiFluidBase::PhaseProp::SliceType phaseEnthalpySlice{phaseEnthalpyValues, phaseEnthalpyDerivs};
//  MultiFluidBase::PhaseProp::SliceType phaseInternalEnergySlice{phaseInternalEnergyValues, phaseInternalEnergyDerivs};
//  MultiFluidBase::PhaseComp::SliceType phaseCompFractionSlice{phaseCompFractionValues, phaseCompFractionDerivs};
//  MultiFluidBase::FluidProp::SliceType totalDensitySlice{totalDensityValue, totalDensityDerivs};
//
//  // Call compute
//  kernelWrapper.compute( 1.0e5, 300.0, composition,
//                         phaseFractionSlice,
//                         phaseDensitySlice,
//                         phaseMassDensitySlice,
//                         phaseViscositySlice,
//                         phaseEnthalpySlice,
//                         phaseInternalEnergySlice,
//                         phaseCompFractionSlice,
//                         totalDensitySlice );
//
//  // The derivative of total density w.r.t phase fraction for phase ip is phase density for that phase
//  // By convention, composition derivatives are at indices 2, 3, 4 (after P, T)
//  for( integer ip = 0; ip < numPhases; ++ip )
//  {
//    // The derivative index for phase fraction ip is ip + 2
//    EXPECT_NEAR( totalDensityDerivs[ip + 2], phaseDensityValues[ip], this->absTol );
//  }
//}

} // namespace testing
} // namespace geos
