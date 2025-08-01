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
struct FluidData< 2 >
{
  static constexpr integer testComponents[numTestComps] = {0, 1, 1};
  static std::unique_ptr< TestFluid< 2 > > createFluid()
  {
    auto fluid = TestFluid< 2 >::create( {Fluid::CH4, Fluid::C5H12} );
    const std::array< real64 const, 2 > bics = { 0.5 };
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
    auto fluid = TestFluid< 4 >::create( {Fluid::CH4, Fluid::CO2, Fluid::H2S, Fluid::H2O} );
    const std::array< real64 const, 6 > bics = { 0.0, 0.0, 0.0, 0.4850, 0.1896, 0.1353 };
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
    if constexpr (EOS == EquationOfStateType::SoreideWhitson)
    {
      eosName = EnumStrings< EquationOfStateType >::toString( EquationOfStateType::PengRobinson );
    }
    equationOfState->m_equationsOfStateNames.emplace_back( eosName );

    auto * criticalVolume = const_cast< CriticalVolume * >(m_parameters->get< CriticalVolume >());
    criticalVolume->m_componentCriticalVolume.resize( numComps );
    TestFluid< NC >::populateArray( criticalVolume->m_componentCriticalVolume, this->m_fluid->criticalVolume );

    auto * flashParameters = const_cast< FlashParameters * >(m_parameters->get< FlashParameters >());
    flashParameters->m_continuousParameters[FlashParameters::STABILITY_TOLERANCE] = 1.0e-12;
    flashParameters->m_continuousParameters[FlashParameters::STABILITY_THRESHOLD] = -1.0e-3;
    flashParameters->m_discreteParameters[FlashParameters::FLASH_MAX_ITERATIONS] = 20;

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
                                       auto const & phaseComponentFractionDerivs,
                                       real64 const scale = 1.0 ){
      integer j = 0;
      for( integer ip = 0; ip < numPhases; ++ip )
      {
        derivs[j++] = scale * phaseFractionDerivs( ip, kc );
        for( integer ic = 0; ic < numComps; ++ic )
        {
          derivs[j++] = scale * phaseComponentFractionDerivs( ip, ic, kc );
        }
      }
    };

    auto const evaluateFlash = [&]( real64 const p, real64 const t, auto const & zmf, auto & values, real64 const scale = 1.0 ){
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
        values[j++] = scale * displacedPhaseFraction[ip];
        for( integer ic = 0; ic < numComps; ++ic )
        {
          values[j++] = scale * displacedPhaseComponentFraction( ip, ic );
        }
      }
    };

    // Test against numerically calculated values
    // --- Pressure derivatives ---
    real64 constexpr pressureScale = 1.0e6;
    concatDerivatives( Deriv::dP, derivatives, dPhaseFraction, dPhaseComponentFraction, pressureScale );
    real64 const dp = 1.0e-4 * pressure;
    geos::testing::internal::testNumericalDerivative< numValues >(
      pressure, dp, derivatives,
      [&]( real64 const p, auto & values ) {
      evaluateFlash( p, temperature, composition.toSliceConst(), values, pressureScale );
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
      if( sumZ < absTol || 1.0 - sumZ < absTol )
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

  void setStabilityThreshold( real64 const threshold )
  {
    auto * flashParameters = const_cast< FlashParameters * >(m_parameters->get< FlashParameters >());
    flashParameters->m_continuousParameters[FlashParameters::STABILITY_THRESHOLD] = threshold;
  }

protected:
  std::unique_ptr< TestFluid< NC > > m_fluid{};
  std::unique_ptr< NegativeTwoPhaseFlashModel > m_flash{};
  std::unique_ptr< ModelParameters > m_parameters{};
  array1d< integer > m_phaseTypes{};
};

using PengRobinson9 = NegativeTwoPhaseFlashModelTestFixture< 9 >;
using SoaveRedlichKwong2 = NegativeTwoPhaseFlashModelTestFixture< 2, EquationOfStateType::SoaveRedlichKwong >;
using SoreideWhitson4 = NegativeTwoPhaseFlashModelTestFixture< 4, EquationOfStateType::SoreideWhitson >;

TEST_P( PengRobinson9, testFlash )
{
  testFlash( GetParam() );
}
TEST_P( SoaveRedlichKwong2, testFlash )
{
  testFlash( GetParam() );
}
TEST_P( SoreideWhitson4, testFlash )
{
  testFlash( GetParam() );
}

TEST_P( PengRobinson9, testFlashDerivatives )
{
  auto const param = GetParam();
  setStabilityThreshold( -1.0e-8 );
  testFlashDerivatives( param );
  setStabilityThreshold( 1.0e-8 );
  testFlashDerivatives( param );
}
TEST_P( SoaveRedlichKwong2, testFlashDerivatives )
{
  auto const param = GetParam();
  setStabilityThreshold( -1.0e-8 );
  testFlashDerivatives( param );
  setStabilityThreshold( 1.0e-8 );
  testFlashDerivatives( param );
}
TEST_P( SoreideWhitson4, testFlashDerivatives )
{
  auto const param = GetParam();
  setStabilityThreshold( -1.0e-8 );
  testFlashDerivatives( param );
  setStabilityThreshold( 1.0e-8 );
  testFlashDerivatives( param );
}

//-------------------------------------------------------------------------------
// Data
//-------------------------------------------------------------------------------

/* UNCRUSTIFY-OFF */

INSTANTIATE_TEST_SUITE_P(
  NegativeTwoPhaseFlashModel, PengRobinson9,
  ::testing::ValuesIn<FlashData<9>>({
    {1.00000e+07, 293.15, { 0.10900000, 0.00300000, 0.23470000, 0.11460000, 0.08790000, 0.04560000, 0.02090000, 0.01510000, 0.36920000 }, 0.78697508, 0.21302492, { 0.13701652, 0.08659405, 0.46883838 }, { 0.00549895, 0.21806200, 0.00110727 }},
    {1.00000e+07, 333.15, { 0.00549900, 0.00222500, 0.69884300, 0.21806200, 0.05952600, 0.01200000, 0.00213100, 0.00061300, 0.00110000 }, 0.00000000, 1.00000000, { 0.00549901, 0.21806223, 0.00110000 }, { 0.00549900, 0.21806200, 0.00110000 }},
    {1.00000e+07, 298.15, { 0.00549900, 0.00222500, 0.69883600, 0.21806200, 0.05952600, 0.01200000, 0.00213100, 0.00061300, 0.00110700 }, 0.00000000, 1.00000000, { 0.00549900, 0.21806200, 0.00110700 }, { 0.00549900, 0.21806200, 0.00110700 }},
    {1.00000e+07, 288.00, { 0.13701652, 0.00320982, 0.10906380, 0.08659405, 0.09558042, 0.05469515, 0.02598050, 0.01902140, 0.46883834 }, 1.00000000, 0.00000000, { 0.13701652, 0.08659405, 0.46883834 }, { 0.13701652, 0.08659405, 0.46883834 }},
    {1.00000e+07, 293.00, { 0.13701652, 0.00320982, 0.10906380, 0.08659405, 0.09558042, 0.05469515, 0.02598050, 0.01902140, 0.46883834 }, 1.00000000, 0.00000000, { 0.13701652, 0.08659405, 0.46883834 }, { 0.13701652, 0.08659405, 0.46883834 }},
    {1.00000e+07, 293.15, { 0.13701652, 0.00320982, 0.10906380, 0.08659405, 0.09558042, 0.05469515, 0.02598050, 0.01902140, 0.46883834 }, 0.99999993, 0.00000007, { 0.13701653, 0.08659404, 0.46883837 }, { 0.00549895, 0.21806197, 0.00110727 }}
  })
);

INSTANTIATE_TEST_SUITE_P(
  NegativeTwoPhaseFlashModel, SoreideWhitson4,
  ::testing::ValuesIn<FlashData<4>>({
    {2.01675e+07, 349.15, { 0.00000000, 0.00180029, 0.00073062, 0.99746900 }, 1.00000000, 0.00000000, { 0.00000000, 0.00180029, 0.99746900 }, { 0.00000000, 0.00180029, 0.99746900 }},
    {1.07250e+07, 349.15, { 0.42400000, 0.05624000, 0.05624000, 0.46351000 }, 0.46391991, 0.53608009, { 0.00109926, 0.00153071, 0.99250187 }, { 0.78998321, 0.10358608, 0.00573284 }},
    {1.06611e+07, 349.15, { 0.42341000, 0.05616000, 0.05616000, 0.46426000 }, 0.46465757, 0.53534243, { 0.00109380, 0.00152509, 0.99252663 }, { 0.78997286, 0.10358214, 0.00575264 }},
    {1.92943e+07, 349.15, { 0.62210000, 0.08283000, 0.08283000, 0.21224000 }, 0.21083142, 0.78916858, { 0.00172629, 0.00206379, 0.99021290 }, { 0.78783679, 0.10440721, 0.00439958 }},
    {2.00000e+07, 349.15, { 0.75000000, 0.10000000, 0.10000000, 0.05000000 }, 0.04630725, 0.95369275, { 0.00176884, 0.00209963, 0.99002996 }, { 0.78633091, 0.10475362, 0.00435616 }},
    {1.99953e+07, 349.15, { 0.96334500, 0.01267200, 0.01266500, 0.01131700 }, 0.00782092, 0.99217908, { 0.00199275, 0.00025091, 0.99699415 }, { 0.97092389, 0.01276992, 0.00354735 }},
    {1.99971e+07, 349.15, { 0.96329100, 0.01267200, 0.01266500, 0.01137300 }, 0.00787742, 0.99212258, { 0.00199287, 0.00025093, 0.99699395 }, { 0.97092271, 0.01277061, 0.00354719 }},
    {1.00000e+07, 293.15, { 0.75000000, 0.05000000, 0.05000000, 0.15000000 }, 0.15072066, 0.84927934, { 0.00150781, 0.00166788, 0.99309883 }, { 0.88283407, 0.05857745, 0.00037619 }},
    {1.00000e+07, 288.15, { 0.00144000, 0.00165800, 0.00007300, 0.99682900 }, 1.00000000, 0.00000000, { 0.00144000, 0.00165800, 0.99682900 }, { 0.00144000, 0.00165800, 0.99682900 }},
    {1.00000e+07, 293.15, { 0.00144000, 0.00165800, 0.00007300, 0.99682900 }, 1.00000000, 0.00000000, { 0.00144000, 0.00165800, 0.99682900 }, { 0.00144000, 0.00165800, 0.99682900 }},
    {1.00000e+07, 353.15, { 0.00144000, 0.00165800, 0.00007300, 0.99682900 }, 0.99958936, 0.00041064, { 0.00108254, 0.00160931, 0.99723579 }, { 0.87159649, 0.12018752, 0.00661262 }},
    {1.00000e+07, 278.15, { 0.88228400, 0.05854300, 0.05882300, 0.00035000 }, 0.00020630, 0.99979370, { 0.00181934, 0.00244237, 0.99161209 }, { 0.88246568, 0.05855458, 0.00014546 }},
    {1.00000e+07, 293.15, { 0.88228400, 0.05854300, 0.05882300, 0.00035000 }, 0.00000000, 1.00000000, { 0.88228400, 0.05854300, 0.00035000 }, { 0.88228400, 0.05854300, 0.00035000 }},
    {1.00000e+07, 303.15, { 0.88228400, 0.05854300, 0.05882300, 0.00035000 }, 0.00000000, 1.00000000, { 0.88228400, 0.05854300, 0.00035000 }, { 0.88228400, 0.05854300, 0.00035000 }},
    {1.00000e+07, 293.15, { 0.00000000, 0.00000500, 0.00000000, 0.99999500 }, 1.00000000, 0.00000000, { 0.00000000, 0.00000500, 0.99999500 }, { 0.00000000, 0.00071931, 0.99928069 }},
    {1.00000e+07, 293.15, { 0.00000000, 0.02572500, 0.00000000, 0.97427500 }, 1.00000000, 0.00000000, { 0.00000000, 0.02572500, 0.97427500 }, { 0.00000000, 0.02572500, 0.97427500 }},
    {1.00000e+07, 293.15, { 0.00000000, 0.02572600, 0.00000000, 0.97427400 }, 0.99999969, 0.00000031, { 0.00000000, 0.02572570, 0.97427430 }, { 0.00000000, 0.99724814, 0.00275186 }},
    {1.00000e+07, 293.15, { 0.99999500, 0.00000000, 0.00000000, 0.00000500 }, 0.00000000, 1.00000000, { 0.99999500, 0.00000000, 0.00000500 }, { 0.99999500, 0.00000000, 0.00000500 }},
    {1.00000e+07, 288.15, { 0.99966400, 0.00000000, 0.00000000, 0.00033600 }, 0.00008980, 0.99991020, { 0.00168896, 0.00000000, 0.99831104 }, { 0.99975363, 0.00000000, 0.00024637 }},
    {1.00000e+07, 293.15, { 0.99966400, 0.00000000, 0.00000000, 0.00033600 }, 0.00000000, 1.00000000, { 0.99966400, 0.00000000, 0.00033600 }, { 0.99966400, 0.00000000, 0.00033600 }},
    {1.00000e+07, 303.15, { 0.99966400, 0.00000000, 0.00000000, 0.00033600 }, 0.00000000, 1.00000000, { 0.99966400, 0.00000000, 0.00033600 }, { 0.99966400, 0.00000000, 0.00033600 }}
  })
);

INSTANTIATE_TEST_SUITE_P(
  NegativeTwoPhaseFlashModel, SoaveRedlichKwong2,
  ::testing::ValuesIn<FlashData<2>>({
    {1.00000e+07, 293.15, { 0.75000000, 0.25000000 }, 0.26178317, 0.73821683, { 0.08177948, 0.91822052, 0.91822052 }, { 0.98696139, 0.01303861, 0.01303861 }},
    {1.00000e+07, 293.15, { 0.08178000, 0.91822000 }, 0.99999943, 0.00000057, { 0.08177948, 0.91822052, 0.91822052 }, { 0.98696139, 0.01303861, 0.01303861 }},
    {1.00000e+07, 293.15, { 0.08177900, 0.91822100 }, 1.00000000, 0.00000000, { 0.08177900, 0.91822100, 0.91822100 }, { 0.08177900, 0.91822100, 0.91822100 }},
    {1.00000e+07, 293.15, { 0.08177800, 0.91822200 }, 1.00000000, 0.00000000, { 0.08177800, 0.91822200, 0.91822200 }, { 0.08177800, 0.91822200, 0.91822200 }},
    {1.00000e+07, 293.15, { 0.98696200, 0.01303800 }, 0.00000000, 1.00000000, { 0.98696200, 0.01303800, 0.01303800 }, { 0.98696200, 0.01303800, 0.01303800 }},
    {1.00000e+07, 293.15, { 0.98696100, 0.01303900 }, 0.00000044, 0.99999956, { 0.08177948, 0.91822052, 0.91822052 }, { 0.98696139, 0.01303861, 0.01303861 }},
    {1.00000e+07, 293.15, { 0.98696000, 0.01304000 }, 0.00000154, 0.99999846, { 0.08177948, 0.91822052, 0.91822052 }, { 0.98696139, 0.01303861, 0.01303861 }}
  })
);

/* UNCRUSTIFY-ON */

} // testing

} // geos
