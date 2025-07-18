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

// Source includes
#include "codingUtilities/UnitTestUtilities.hpp"
#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/EquationOfState.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/CriticalVolume.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/FlashParameters.hpp"
#include "constitutive/fluid/multifluid/compositional/models/NegativeTwoPhaseFlashModel.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/PhaseType.hpp"
#include "TestFluid.hpp"
#include "TestFluidUtilities.hpp"

using namespace geos::constitutive;
using namespace geos::constitutive::compositional;

namespace geos
{
namespace testing
{

constexpr integer numTestComps = 3;

template< integer NC >
using FlashData = std::tuple<
  real64 const,                 // pressure
  real64 const,                 // temperature
  Feed< NC > const,             // phase composition
  real64 const,                 // expected liquid fraction
  real64 const,                 // expected vapour fraction
  Feed< numTestComps > const,   // expected liquid mole fractions
  Feed< numTestComps > const    // expected vapour mole fractions
  >;

template< integer NC >
struct FluidData {};

template<>
struct FluidData< 3 >
{
  static constexpr integer testComponents[numTestComps] = {0, 1, 2};
  static std::unique_ptr< TestFluid< 3 > > createFluid()
  {
    auto fluid = TestFluid< 3 >::create( {Fluid::CH4, Fluid::C10H22, Fluid::H2O} );
    const std::array< real64 const, 3 > bics = { 0.25, 0.0, 0.0 };
    fluid->setBinaryCoefficients( bics );
    return fluid;
  }
};

template<>
struct FluidData< 4 >
{
  static constexpr integer testComponents[numTestComps] = {0, 1, 3};
  static std::unique_ptr< TestFluid< 4 > > createFluid()
  {
    auto fluid = TestFluid< 4 >::create( {Fluid::CH4, Fluid::CO2, Fluid::N2, Fluid::H2O} );
    const std::array< real64 const, 6 > bics = { 0.0, 0.0, 0.0, 0.4850, 0.1896, 0.4778 };
    fluid->setBinaryCoefficients( bics );
    return fluid;
  }
};

template<>
struct FluidData< 9 >
{
  static constexpr integer testComponents[numTestComps] = {0, 3, 8};
  static std::unique_ptr< TestFluid< 9 > > createFluid()
  {
    auto fluid = TestFluid< 9 >::create( {Fluid::H2O, Fluid::CO2, Fluid::N2, Fluid::CH4, Fluid::C2H6, Fluid::C3H8, Fluid::C4H10, Fluid::C5H12, Fluid::C8H18} );
    const std::array< real64 const, 36 > bics = {
      0.01, 0, 0.003732, 0, 0.01, 0, 0, 0.01, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0.01, 0, 0.028, 0.01, 0.01, 0, 0, 0.01, 0, 0.04532, 0.01, 0.01, 0, 0, 0
    };
    fluid->setBinaryCoefficients( bics );
    return fluid;
  }
};

template< integer NC, EquationOfStateType EOS = EquationOfStateType::PengRobinson >
class NegativeTwoPhaseFlashModelTestFixture :  public ::testing::TestWithParam< FlashData< NC > >
{
  static constexpr real64 relTol = 1.0e-5;
  static constexpr real64 absTol = 1.0e-7;
  static constexpr integer numPhases = 2;
  static constexpr integer numComps = NC;
  static constexpr integer numDofs = NC + 2;
  using Deriv = geos::constitutive::multifluid::DerivativeOffset;
  using PhasePropSlice = NegativeTwoPhaseFlashModelUpdate::PhaseProp::SliceType;
  using PhaseCompSlice = NegativeTwoPhaseFlashModelUpdate::PhaseComp::SliceType;

  using CompProp =  StackArray< real64, 3, numPhases, multifluid::LAYOUT_PHASE >;
  using PhaseProp = StackArray< real64, 4, numPhases *numDofs, multifluid::LAYOUT_PHASE_DC >;
  using PhaseComp = StackArray< real64, 5, numPhases *numComps *numDofs, multifluid::LAYOUT_PHASE_COMP_DC >;

public:
  NegativeTwoPhaseFlashModelTestFixture()
    : m_fluid( FluidData< NC >::createFluid() )
  {
    ComponentProperties const & componentProperties = this->m_fluid->getComponentProperties();

    m_parameters = NegativeTwoPhaseFlashModel::createParameters( std::move( m_parameters ) );

    auto * equationOfState = const_cast< EquationOfState * >(m_parameters->get< EquationOfState >());
    string eosName = EnumStrings< EquationOfStateType >::toString( EOS );
    equationOfState->m_equationsOfStateNames.emplace_back( eosName );
    eosName = EnumStrings< EquationOfStateType >::toString( EquationOfStateType::PengRobinson );
    equationOfState->m_equationsOfStateNames.emplace_back( eosName );

    auto * criticalVolume = const_cast< CriticalVolume * >(m_parameters->get< CriticalVolume >());
    criticalVolume->m_componentCriticalVolume.resize( numComps );
    TestFluid< NC >::populateArray( criticalVolume->m_componentCriticalVolume, this->m_fluid->criticalVolume );

    auto * flashParameters = const_cast< FlashParameters * >(m_parameters->get< FlashParameters >());
    flashParameters->m_continuousParameters[FlashParameters::STABILITY_TOLERANCE] = 1.0e-5;

    m_phaseTypes.emplace_back( static_cast< integer >(PhaseType::LIQUID));
    m_phaseTypes.emplace_back( static_cast< integer >(PhaseType::VAPOUR));

    m_flash = std::make_unique< NegativeTwoPhaseFlashModel >( "FlashModel", componentProperties, *m_parameters, m_phaseTypes );
  }

  ~NegativeTwoPhaseFlashModelTestFixture() = default;

  void testFlash( FlashData< NC > const & data )
  {
    real64 const pressure = std::get< 0 >( data );
    real64 const temperature = std::get< 1 >( data );
    stackArray1d< real64, numComps > composition;
    TestFluid< NC >::createArray( composition, std::get< 2 >( data ));

    real64 const expectedPhaseFraction[numPhases] = {std::get< 3 >( data ), std::get< 4 >( data ) };
    Feed< numTestComps > const expectedPhaseComponentFraction[numPhases] = {std::get< 5 >( data ), std::get< 6 >( data ) };

    stackArray2d< real64, (numPhases-1)*numComps > kValues( numPhases-1, numComps );
    LvArray::forValuesInSlice( kValues.toSlice(), []( real64 & v ){ v = 0.0; } );

    CompProp phaseFractionData( 1, 1, numPhases );
    PhaseProp dPhaseFractionData( 1, 1, numPhases, numDofs );
    PhaseProp phaseComponentFractionData( 1, 1, numPhases, numComps );
    PhaseComp dPhaseComponentFractionData( 1, 1, numPhases, numComps, numDofs );

    auto phaseFraction = phaseFractionData[0][0];
    auto dPhaseFraction = dPhaseFractionData[0][0];
    auto phaseComponentFraction = phaseComponentFractionData[0][0];
    auto dPhaseComponentFraction = dPhaseComponentFractionData[0][0];

    auto componentProperties = m_fluid->createKernelWrapper();
    auto flashKernelWrapper = m_flash->createKernelWrapper();

    flashKernelWrapper.compute( componentProperties,
                                pressure,
                                temperature,
                                composition.toSliceConst(),
                                kValues.toSlice(),
                                PhasePropSlice( phaseFraction, dPhaseFraction ),
                                PhaseCompSlice( phaseComponentFraction, dPhaseComponentFraction ) );

    for( integer ip = 0; ip < numPhases; ip++ )
    {
      checkRelativeError( phaseFraction[ip], expectedPhaseFraction[ip], relTol, absTol );
      for( integer i = 0; i < numTestComps; ++i )
      {
        integer const ic = FluidData< numComps >::testComponents[i];
        checkRelativeError( phaseComponentFraction[ip][ic], expectedPhaseComponentFraction[ip][i], relTol, absTol );
      }
    }
  }

  void testFlashDerivatives( FlashData< NC > const & data )
  {
    std::cout << std::fixed << std::setprecision( 8 );
    // Number of output values from each flash calculation
    constexpr integer numValues = numPhases * (1 + numComps);

    real64 const pressure = std::get< 0 >( data );
    real64 const temperature = std::get< 1 >( data );
    stackArray1d< real64, numComps > composition;
    TestFluid< NC >::createArray( composition, std::get< 2 >( data ));

    stackArray2d< real64, (numPhases-1)*numComps > kValues( numPhases-1, numComps );
    LvArray::forValuesInSlice( kValues.toSlice(), []( real64 & v ){ v = 0.0; } );

    CompProp phaseFractionData( 1, 1, numPhases );
    PhaseProp dPhaseFractionData( 1, 1, numPhases, numDofs );
    PhaseProp phaseComponentFractionData( 1, 1, numPhases, numComps );
    PhaseComp dPhaseComponentFractionData( 1, 1, numPhases, numComps, numDofs );

    auto phaseFraction = phaseFractionData[0][0];
    auto dPhaseFraction = dPhaseFractionData[0][0];
    auto phaseComponentFraction = phaseComponentFractionData[0][0];
    auto dPhaseComponentFraction = dPhaseComponentFractionData[0][0];

    stackArray1d< real64, numValues > derivatives( numValues );

    auto componentProperties = m_fluid->createKernelWrapper();
    auto flashKernelWrapper = m_flash->createKernelWrapper();
    flashKernelWrapper.compute( componentProperties,
                                pressure,
                                temperature,
                                composition.toSliceConst(),
                                kValues.toSlice(),
                                PhasePropSlice( phaseFraction, dPhaseFraction ),
                                PhaseCompSlice( phaseComponentFraction, dPhaseComponentFraction ) );

    // Combine derivatives into a single output
    auto const concatDerivatives = []( integer const kc, auto & derivs, auto const & phaseFractionDerivs,
                                       auto const & phaseComponentFractionDerivs ){
      integer j = 0;
      for( integer ip = 0; ip < numPhases; ++ip )
      {
        derivs[j++] = phaseFractionDerivs( ip, kc );
        for( integer ic = 0; ic < numComps; ++ic )
        {
          derivs[j++] = phaseComponentFractionDerivs( ip, ic, kc );
        }
      }
    };

    auto const evaluateFlash = [&]( real64 const p, real64 const t, auto const & zmf, auto & values ){
      CompProp displacedPhaseFractionData( 1, 1, numPhases );
      PhaseProp displacedPhaseFractionDerivsData( 1, 1, numPhases, numDofs );
      PhaseProp displacedPhaseComponentFractionData( 1, 1, numPhases, numComps );
      PhaseComp displacedPhaseComponentFractionDerivsData( 1, 1, numPhases, numComps, numDofs );

      auto displacedPhaseFraction = displacedPhaseFractionData[0][0];
      auto displacedPhaseFractionDerivs = displacedPhaseFractionDerivsData[0][0];
      auto displacedPhaseComponentFraction = displacedPhaseComponentFractionData[0][0];
      auto displacedPhaseComponentFractionDerivs = displacedPhaseComponentFractionDerivsData[0][0];

      LvArray::forValuesInSlice( kValues.toSlice(), []( real64 & v ){ v = 0.0; } );

      flashKernelWrapper.compute( componentProperties,
                                  p,
                                  t,
                                  zmf,
                                  kValues.toSlice(),
                                  PhasePropSlice( displacedPhaseFraction, displacedPhaseFractionDerivs ),
                                  PhaseCompSlice( displacedPhaseComponentFraction, displacedPhaseComponentFractionDerivs ) );
      integer j = 0;
      for( integer ip = 0; ip < numPhases; ++ip )
      {
        values[j++] = displacedPhaseFraction[ip];
        for( integer ic = 0; ic < numComps; ++ic )
        {
          values[j++] = displacedPhaseComponentFraction( ip, ic );
        }
      }
    };

    // Test against numerically calculated values
    // --- Pressure derivatives ---
    concatDerivatives( Deriv::dP, derivatives, dPhaseFraction, dPhaseComponentFraction );
    real64 const dp = 1.0e-4 * pressure;
    geos::testing::internal::testNumericalDerivative< numValues >(
      pressure, dp, derivatives,
      [&]( real64 const p, auto & values ) {
      evaluateFlash( p, temperature, composition.toSliceConst(), values );
    } );

    // -- Temperature derivative
    concatDerivatives( Deriv::dT, derivatives, dPhaseFraction, dPhaseComponentFraction );
    real64 const dT = 1.0e-4 * temperature;
    geos::testing::internal::testNumericalDerivative< numValues >(
      temperature, dT, derivatives,
      [&]( real64 const t, auto & values ) {
      evaluateFlash( pressure, t, composition.toSliceConst(), values );
    } );
    // -- Composition derivatives derivative
    real64 const dz = 1.0e-6;
    for( integer const ic : FluidData< numComps >::testComponents )
    {
      real64 sumZ = 0.0;
      for( integer jc = 0; jc < numComps; ++jc )
      {
        sumZ += composition[jc];
      }
      sumZ -= composition[ic];
      if( sumZ < absTol )
      {
        continue;
      }
      concatDerivatives( Deriv::dC+ic, derivatives, dPhaseFraction, dPhaseComponentFraction );
      geos::testing::internal::testNumericalDerivative< numValues >(
        0.0, dz, derivatives,
        [&]( real64 const z, auto & values ) {
        stackArray1d< real64, numComps > zmf( numComps );
        for( integer jc = 0; jc < numComps; ++jc )
        {
          zmf[jc] = composition[jc];
        }
        zmf[ic] += z;
        evaluateFlash( pressure, temperature, zmf.toSliceConst(), values );
      } );
    }
  }

protected:
  std::unique_ptr< TestFluid< NC > > m_fluid{};
  std::unique_ptr< NegativeTwoPhaseFlashModel > m_flash{};
  std::unique_ptr< ModelParameters > m_parameters{};
  array1d< integer > m_phaseTypes{};
};

//using NegativeTwoPhaseFlashModel3 = NegativeTwoPhaseFlashModelTestFixture< 3 >;
//using PengRobinson9 = NegativeTwoPhaseFlashModelTestFixture< 9 >;
using SoreideWhitson4 = NegativeTwoPhaseFlashModelTestFixture< 4, EquationOfStateType::SoreideWhitson >;

//TEST_P( NegativeTwoPhaseFlashModel3, testFlash )
//{
//  testFlash( GetParam() );
//}
//TEST_P( NegativeTwoPhaseFlashModel3, testFlashDerivatives )
//{
//  testFlashDerivatives( GetParam() );
//}

//TEST_P( NegativeTwoPhaseFlashModel9, testFlash )
//{
//  testFlash( GetParam() );
//}
//TEST_P( PengRobinson9, testFlashDerivatives )
//{
//  testFlashDerivatives( GetParam() );
//}
TEST_P( SoreideWhitson4, testFlashDerivatives )
{
  testFlashDerivatives( GetParam() );
}

//-------------------------------------------------------------------------------
// Data
//-------------------------------------------------------------------------------
/* UNCRUSTIFY-OFF */
/**
INSTANTIATE_TEST_SUITE_P(
  NegativeTwoPhaseFlashModel, PengRobinson9,
  ::testing::ValuesIn<FlashData<9>>({
    //{1.0000e+07, 293.15, {0.10900000, 0.00300000, 0.23470000, 0.11460000, 0.08790000, 0.04560000, 0.02090000, 0.01510000, 0.36920000}, 0.185888, 0.414112, {0.000000, 0.111648, 0.551139}, {0.000000, 0.731628, 0.000000}}
    //{1.5000e+07, 293.15, {0.13701652, 0.00320982, 0.10906376, 0.08659405, 0.09558042, 0.05469515, 0.02598050, 0.01902140, 0.46883838}, 0.185888, 0.414112, {0.000000, 0.111648, 0.551139}, {0.000000, 0.731628, 0.000000}},
    //{1.0000e+07, 293.15, {0.13701652, 0.00320982, 0.10906370, 0.08659405, 0.09558042, 0.05469515, 0.02598050, 0.01902140, 0.46883844}, 0.185888, 0.414112, {0.000000, 0.111648, 0.551139}, {0.000000, 0.731628, 0.000000}}
    {1.0000e+07, 293.15, {0.13701652, 0.00320982, 0.10906370, 0.00000000, 0.09558042, 0.05469515, 0.02598050, 0.01902140, 0.55543249}, 0.185888, 0.414112, {0.000000, 0.111648, 0.551139}, {0.000000, 0.731628, 0.000000}}
    //{1.0001e+07, 293.15, {0.13701652, 0.00320982, 0.10906376, 0.08659405, 0.09558042, 0.05469515, 0.02598050, 0.01902140, 0.46883838}, 0.185888, 0.414112, {0.000000, 0.111648, 0.551139}, {0.000000, 0.731628, 0.000000}}
    //FlashData<9>{1.0000e+07, 293.15, {0.00396181, 0.00219389, 0.58376764, 0.30883583, 0.06016813, 0.01165533, 0.00271459, 0.02656289, 0.00013988}, 0.185888, 0.414112, {0.000000, 0.111648, 0.551139}, {0.000000, 0.731628, 0.000000}}
  })
);
*/

INSTANTIATE_TEST_SUITE_P(
  NegativeTwoPhaseFlashModel, SoreideWhitson4,
  ::testing::ValuesIn<FlashData<4>>({
    //{1.0e+07, 293.15, {0.000000, 0.000000, 0.000000, 1.000000}, 0.000000, 1.000000, {0.000000, 0.125000, 0.125000}, {0.000000, 0.125000, 0.125000}}
    //{1.0e+07, 293.15, {0.000005, 0.000005, 0.000005, 0.999985}, 0.000000, 1.000000, {0.000000, 0.125000, 0.125000}, {0.000000, 0.125000, 0.125000}}
    //{1.0e+07, 293.15, {0.000000, 1.000000, 0.000000, 0.000000}, 0.000000, 1.000000, {0.000000, 0.125000, 0.125000}, {0.000000, 0.125000, 0.125000}}
    {1.0e+07, 293.15, {0.000005, 0.999985, 0.000005, 0.000005}, 0.000000, 1.000000, {0.000000, 0.125000, 0.125000}, {0.000000, 0.125000, 0.125000}}
  })
);

/* UNCRUSTIFY-ON */

} // testing

} // geos
