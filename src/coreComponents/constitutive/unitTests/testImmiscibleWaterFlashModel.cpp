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
#include "constitutive/fluid/multifluid/compositional/parameters/PhaseType.hpp"
#include "constitutive/fluid/multifluid/compositional/models/ImmiscibleWaterFlashModel.hpp"
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
  real64 const,                 // expected aqueous fraction
  Feed< numTestComps > const,   // expected liquid mole fractions
  Feed< numTestComps > const,   // expected vapour mole fractions
  Feed< numTestComps > const    // expected aqueous mole fractions
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
    const stdArray< real64 const, 3 > bics = { 0.25, 0.0, 0.0 };
    fluid->setBinaryCoefficients( bics );
    return fluid;
  }
};

template<>
struct FluidData< 9 >
{
  static constexpr integer testComponents[numTestComps] = {0, 2, 8};
  static std::unique_ptr< TestFluid< 9 > > createFluid()
  {
    auto fluid = TestFluid< 9 >::create( {Fluid::H2O, Fluid::CO2, Fluid::N2, Fluid::C5H12, Fluid::C2H6, Fluid::C3H8, Fluid::C4H10, Fluid::C5H12, Fluid::C10H22} );
    const stdArray< real64, 36 > bics = {
      0.01, 0, 0.003732, 0, 0.01, 0, 0, 0.01, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0.01, 0, 0.028, 0.01, 0.01, 0, 0, 0.01, 0, 0.04532, 0.01, 0.01, 0, 0, 0
    };
    fluid->setBinaryCoefficients( bics );
    return fluid;
  }
};

template< integer NC >
class ImmiscibleWaterFlashModelTestFixture :  public ::testing::TestWithParam< FlashData< NC > >
{
  static constexpr real64 relTol = 1.0e-5;
  static constexpr real64 absTol = 1.0e-7;
  static constexpr integer numPhases = 3;
  static constexpr integer numComps = NC;
  static constexpr integer numDofs = NC + 2;
  using Deriv = geos::constitutive::multifluid::DerivativeOffset;
  using PhasePropSlice = ImmiscibleWaterFlashModelUpdate::PhaseProp::SliceType;
  using PhaseCompSlice = ImmiscibleWaterFlashModelUpdate::PhaseComp::SliceType;

public:
  ImmiscibleWaterFlashModelTestFixture()
    : m_fluid( FluidData< NC >::createFluid() )
  {
    ComponentProperties const & componentProperties = this->m_fluid->getComponentProperties();

    m_parameters = ImmiscibleWaterFlashModel::createParameters( std::move( m_parameters ) );

    auto * equationOfState = const_cast< EquationOfState * >(m_parameters->get< EquationOfState >());
    string const eosName = EnumStrings< EquationOfStateType >::toString( EquationOfStateType::PengRobinson );
    equationOfState->m_equationsOfStateNames.emplace_back( eosName );
    equationOfState->m_equationsOfStateNames.emplace_back( eosName );
    equationOfState->m_equationsOfStateNames.emplace_back( eosName );

    auto * criticalVolume = const_cast< CriticalVolume * >(m_parameters->get< CriticalVolume >());
    TestFluid< NC >::createArray( criticalVolume->m_componentCriticalVolume, this->m_fluid->criticalVolume );

    m_phaseTypes.emplace_back( static_cast< integer >(PhaseType::LIQUID));
    m_phaseTypes.emplace_back( static_cast< integer >(PhaseType::VAPOUR));
    m_phaseTypes.emplace_back( static_cast< integer >(PhaseType::AQUEOUS));

    m_flash = std::make_unique< ImmiscibleWaterFlashModel >( "FlashModel", componentProperties, *m_parameters, m_phaseTypes );
  }

  ~ImmiscibleWaterFlashModelTestFixture() = default;

  void testFlash( FlashData< NC > const & data )
  {
    real64 const pressure = std::get< 0 >( data );
    real64 const temperature = std::get< 1 >( data );
    stackArray1d< real64, numComps > composition;
    TestFluid< NC >::createArray( composition, std::get< 2 >( data ));

    real64 const expectedPhaseFraction[numPhases] = {std::get< 3 >( data ), std::get< 4 >( data ), std::get< 5 >( data )};
    Feed< numTestComps > const expectedPhaseComponentFraction[numPhases] = {std::get< 6 >( data ), std::get< 7 >( data ), std::get< 8 >( data )    };

    stackArray2d< real64, (numPhases-1)*numComps > kValues( numPhases-1, numComps );
    LvArray::forValuesInSlice( kValues.toSlice(), []( real64 & v ){ v = 0.0; } );

    StackArray< real64, 3, numPhases, multifluid::LAYOUT_PHASE > phaseFractionData( 1, 1, numPhases );
    StackArray< real64, 4, numPhases *numDofs, multifluid::LAYOUT_PHASE_DC > dPhaseFractionData( 1, 1, numPhases, numDofs );
    StackArray< real64, 4, numPhases *numComps, multifluid::LAYOUT_PHASE_COMP > phaseComponentFractionData( 1, 1, numPhases, numComps );
    StackArray< real64, 5, numPhases *numComps *numDofs, multifluid::LAYOUT_PHASE_COMP_DC > dPhaseComponentFractionData( 1, 1, numPhases, numComps, numDofs );

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

    StackArray< real64, 3, numPhases, multifluid::LAYOUT_PHASE > phaseFractionData( 1, 1, numPhases );
    StackArray< real64, 4, numPhases *numDofs, multifluid::LAYOUT_PHASE_DC > dPhaseFractionData( 1, 1, numPhases, numDofs );
    StackArray< real64, 4, numPhases *numComps, multifluid::LAYOUT_PHASE_COMP > phaseComponentFractionData( 1, 1, numPhases, numComps );
    StackArray< real64, 5, numPhases *numComps *numDofs, multifluid::LAYOUT_PHASE_COMP_DC > dPhaseComponentFractionData( 1, 1, numPhases, numComps, numDofs );

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
    auto const concatDerivatives = []( integer const kc, auto & derivs, auto const & phaseFractionDerivs, auto const & phaseComponentFractionDerivs ){
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
      StackArray< real64, 3, numPhases, multifluid::LAYOUT_PHASE > displacedPhaseFractionData( 1, 1, numPhases );
      StackArray< real64, 4, numPhases *numDofs, multifluid::LAYOUT_PHASE_DC > displacedPhaseFractionDerivsData( 1, 1, numPhases, numDofs );
      StackArray< real64, 4, numPhases *numComps, multifluid::LAYOUT_PHASE_COMP > displacedPhaseComponentFractionData( 1, 1, numPhases, numComps );
      StackArray< real64, 5, numPhases *numComps *numDofs, multifluid::LAYOUT_PHASE_COMP_DC > displacedPhaseComponentFractionDerivsData( 1, 1, numPhases, numComps, numDofs );

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
    real64 const dT = 1.0e-6 * temperature;
    geos::testing::internal::testNumericalDerivative< numValues >(
      temperature, dT, derivatives,
      [&]( real64 const t, auto & values ) {
      evaluateFlash( pressure, t, composition.toSliceConst(), values );
    } );

    // -- Composition derivatives derivative
    real64 const dz = 1.0e-7;
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
  std::unique_ptr< ImmiscibleWaterFlashModel > m_flash{};
  std::unique_ptr< ModelParameters > m_parameters{};
  array1d< integer > m_phaseTypes{};
};

using ImmiscibleWaterFlashModel3 = ImmiscibleWaterFlashModelTestFixture< 3 >;
using ImmiscibleWaterFlashModel9 = ImmiscibleWaterFlashModelTestFixture< 9 >;

TEST_P( ImmiscibleWaterFlashModel3, testFlash )
{
  testFlash( GetParam() );
}
TEST_P( ImmiscibleWaterFlashModel3, testFlashDerivatives )
{
  testFlashDerivatives( GetParam() );
}

TEST_P( ImmiscibleWaterFlashModel9, testFlash )
{
  testFlash( GetParam() );
}
TEST_P( ImmiscibleWaterFlashModel9, testFlashDerivatives )
{
  testFlashDerivatives( GetParam() );
}

//-------------------------------------------------------------------------------
// Data
//-------------------------------------------------------------------------------
/* UNCRUSTIFY-OFF */

INSTANTIATE_TEST_SUITE_P(
  ImmiscibleWaterFlashModel, ImmiscibleWaterFlashModel3,
  ::testing::ValuesIn<FlashData<3>>({
    {1.0e+05, 293.15, {0.000000, 0.000000, 1.000000}, 0.000000, 0.000000, 1.000000, {0.500000, 0.500000, 0.000000}, {0.500000, 0.500000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    {1.0e+06, 293.15, {0.000000, 0.000000, 1.000000}, 0.000000, 0.000000, 1.000000, {0.500000, 0.500000, 0.000000}, {0.500000, 0.500000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    {1.0e+07, 293.15, {0.000000, 0.000000, 1.000000}, 0.000000, 0.000000, 1.000000, {0.500000, 0.500000, 0.000000}, {0.500000, 0.500000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    {1.0e+05, 313.15, {0.000000, 0.000000, 1.000000}, 0.000000, 0.000000, 1.000000, {0.500000, 0.500000, 0.000000}, {0.500000, 0.500000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    {1.0e+06, 313.15, {0.000000, 0.000000, 1.000000}, 0.000000, 0.000000, 1.000000, {0.500000, 0.500000, 0.000000}, {0.500000, 0.500000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    {1.0e+07, 313.15, {0.000000, 0.000000, 1.000000}, 0.000000, 0.000000, 1.000000, {0.500000, 0.500000, 0.000000}, {0.500000, 0.500000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    {1.0e+05, 353.15, {0.000000, 0.000000, 1.000000}, 0.000000, 0.000000, 1.000000, {0.500000, 0.500000, 0.000000}, {0.500000, 0.500000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    {1.0e+06, 353.15, {0.000000, 0.000000, 1.000000}, 0.000000, 0.000000, 1.000000, {0.500000, 0.500000, 0.000000}, {0.500000, 0.500000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    {1.0e+07, 353.15, {0.000000, 0.000000, 1.000000}, 0.000000, 0.000000, 1.000000, {0.500000, 0.500000, 0.000000}, {0.500000, 0.500000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    {1.0e+05, 293.15, {0.300000, 0.300000, 0.400000}, 0.300287, 0.299713, 0.400000, {0.002474, 0.997526, 0.000000}, {0.998479, 0.001521, 0.000000}, {0.000000, 0.000000, 1.000000}},
    {1.0e+06, 293.15, {0.300000, 0.300000, 0.400000}, 0.307342, 0.292658, 0.400000, {0.024081, 0.975919, 0.000000}, {0.999798, 0.000202, 0.000000}, {0.000000, 0.000000, 1.000000}},
    {1.0e+07, 293.15, {0.300000, 0.300000, 0.400000}, 0.366467, 0.233533, 0.400000, {0.181547, 0.818453, 0.000000}, {0.999725, 0.000275, 0.000000}, {0.000000, 0.000000, 1.000000}},
    {1.0e+05, 313.15, {0.300000, 0.300000, 0.400000}, 0.299054, 0.300946, 0.400000, {0.002386, 0.997614, 0.000000}, {0.994485, 0.005515, 0.000000}, {0.000000, 0.000000, 1.000000}},
    {1.0e+06, 313.15, {0.300000, 0.300000, 0.400000}, 0.306975, 0.293025, 0.400000, {0.023388, 0.976612, 0.000000}, {0.999301, 0.000699, 0.000000}, {0.000000, 0.000000, 1.000000}},
    {1.0e+07, 313.15, {0.300000, 0.300000, 0.400000}, 0.367247, 0.232753, 0.400000, {0.183491, 0.816509, 0.000000}, {0.999398, 0.000602, 0.000000}, {0.000000, 0.000000, 1.000000}},
    {1.0e+05, 353.15, {0.300000, 0.300000, 0.400000}, 0.287017, 0.312983, 0.400000, {0.002244, 0.997756, 0.000000}, {0.956460, 0.043540, 0.000000}, {0.000000, 0.000000, 1.000000}},
    {1.0e+06, 353.15, {0.300000, 0.300000, 0.400000}, 0.305469, 0.294531, 0.400000, {0.022889, 0.977111, 0.000000}, {0.994831, 0.005169, 0.000000}, {0.000000, 0.000000, 1.000000}},
    {1.0e+07, 353.15, {0.300000, 0.300000, 0.400000}, 0.369634, 0.230366, 0.400000, {0.189883, 0.810117, 0.000000}, {0.997597, 0.002403, 0.000000}, {0.000000, 0.000000, 1.000000}},
    {1.0e+05, 293.15, {0.200000, 0.800000, 0.000000}, 0.801682, 0.198318, 0.000000, {0.002474, 0.997526, 0.000000}, {0.998479, 0.001521, 0.000000}, {0.000000, 0.000000, 1.000000}},
    {1.0e+06, 293.15, {0.200000, 0.800000, 0.000000}, 0.819703, 0.180297, 0.000000, {0.024081, 0.975919, 0.000000}, {0.999798, 0.000202, 0.000000}, {0.000000, 0.000000, 1.000000}},
    {1.0e+07, 293.15, {0.200000, 0.800000, 0.000000}, 0.977446, 0.022554, 0.000000, {0.181547, 0.818453, 0.000000}, {0.999725, 0.000275, 0.000000}, {0.000000, 0.000000, 1.000000}},
    {1.0e+05, 313.15, {0.200000, 0.800000, 0.000000}, 0.800812, 0.199188, 0.000000, {0.002386, 0.997614, 0.000000}, {0.994485, 0.005515, 0.000000}, {0.000000, 0.000000, 1.000000}},
    {1.0e+06, 313.15, {0.200000, 0.800000, 0.000000}, 0.819029, 0.180971, 0.000000, {0.023388, 0.976612, 0.000000}, {0.999301, 0.000699, 0.000000}, {0.000000, 0.000000, 1.000000}},
    {1.0e+07, 313.15, {0.200000, 0.800000, 0.000000}, 0.979766, 0.020234, 0.000000, {0.183491, 0.816509, 0.000000}, {0.999398, 0.000602, 0.000000}, {0.000000, 0.000000, 1.000000}},
    {1.0e+05, 353.15, {0.200000, 0.800000, 0.000000}, 0.792756, 0.207244, 0.000000, {0.002244, 0.997756, 0.000000}, {0.956460, 0.043540, 0.000000}, {0.000000, 0.000000, 1.000000}},
    {1.0e+06, 353.15, {0.200000, 0.800000, 0.000000}, 0.817776, 0.182224, 0.000000, {0.022889, 0.977111, 0.000000}, {0.994831, 0.005169, 0.000000}, {0.000000, 0.000000, 1.000000}},
    {1.0e+07, 353.15, {0.200000, 0.800000, 0.000000}, 0.987475, 0.012525, 0.000000, {0.189883, 0.810117, 0.000000}, {0.997597, 0.002403, 0.000000}, {0.000000, 0.000000, 1.000000}}
  })
);

INSTANTIATE_TEST_SUITE_P(
  ImmiscibleWaterFlashModel, ImmiscibleWaterFlashModel9,
  ::testing::ValuesIn<FlashData<9>>({
    {1.0e+06, 293.15, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 0.332132, 0.658868, 0.009000, {0.000000, 0.015946, 0.509203}, {0.000000, 0.803506, 0.000117}, {1.000000, 0.000000, 0.000000}},
    {1.0e+07, 293.15, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 0.488745, 0.502255, 0.009000, {0.000000, 0.163221, 0.346062}, {0.000000, 0.905768, 0.000128}, {1.000000, 0.000000, 0.000000}},
    {1.0e+08, 293.15, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 0.991000, 0.000000, 0.009000, {0.000000, 0.539556, 0.170737}, {0.000000, 0.539556, 0.170737}, {1.000000, 0.000000, 0.000000}},
    {1.0e+06, 313.15, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 0.302693, 0.688307, 0.009000, {0.000000, 0.014899, 0.557960}, {0.000000, 0.770281, 0.000450}, {1.000000, 0.000000, 0.000000}},
    {1.0e+07, 313.15, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 0.470546, 0.520454, 0.009000, {0.000000, 0.158760, 0.359214}, {0.000000, 0.883836, 0.000333}, {1.000000, 0.000000, 0.000000}},
    {1.0e+08, 313.15, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 0.991000, 0.000000, 0.009000, {0.000000, 0.539556, 0.170737}, {0.000000, 0.539556, 0.170737}, {1.000000, 0.000000, 0.000000}},
    {1.0e+06, 353.15, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 0.240571, 0.750429, 0.009000, {0.000000, 0.013515, 0.690435}, {0.000000, 0.708194, 0.004132}, {1.000000, 0.000000, 0.000000}},
    {1.0e+07, 353.15, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 0.431744, 0.559256, 0.009000, {0.000000, 0.154806, 0.389666}, {0.000000, 0.836581, 0.001724}, {1.000000, 0.000000, 0.000000}},
    {1.0e+08, 353.15, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 0.991000, 0.000000, 0.009000, {0.000000, 0.539556, 0.170737}, {0.000000, 0.539556, 0.170737}, {1.000000, 0.000000, 0.000000}},
    {1.0e+06, 293.15, {0.400000, 0.001820, 0.323730, 0.069380, 0.053220, 0.027610, 0.012650, 0.009140, 0.102450}, 0.201090, 0.398910, 0.400000, {0.000000, 0.015946, 0.509241}, {0.000000, 0.803498, 0.000117}, {1.000000, 0.000000, 0.000000}},
    {1.0e+07, 293.15, {0.400000, 0.001820, 0.323730, 0.069380, 0.053220, 0.027610, 0.012650, 0.009140, 0.102450}, 0.295912, 0.304088, 0.400000, {0.000000, 0.163219, 0.346087}, {0.000000, 0.905762, 0.000128}, {1.000000, 0.000000, 0.000000}},
    {1.0e+08, 293.15, {0.400000, 0.001820, 0.323730, 0.069380, 0.053220, 0.027610, 0.012650, 0.009140, 0.102450}, 0.600000, 0.000000, 0.400000, {0.000000, 0.539550, 0.170750}, {0.000000, 0.539550, 0.170750}, {1.000000, 0.000000, 0.000000}},
    {1.0e+06, 313.15, {0.400000, 0.001820, 0.323730, 0.069380, 0.053220, 0.027610, 0.012650, 0.009140, 0.102450}, 0.183267, 0.416733, 0.400000, {0.000000, 0.014899, 0.557997}, {0.000000, 0.770276, 0.000450}, {1.000000, 0.000000, 0.000000}},
    {1.0e+07, 313.15, {0.400000, 0.001820, 0.323730, 0.069380, 0.053220, 0.027610, 0.012650, 0.009140, 0.102450}, 0.284893, 0.315107, 0.400000, {0.000000, 0.158758, 0.359240}, {0.000000, 0.883830, 0.000333}, {1.000000, 0.000000, 0.000000}},
    {1.0e+08, 313.15, {0.400000, 0.001820, 0.323730, 0.069380, 0.053220, 0.027610, 0.012650, 0.009140, 0.102450}, 0.600000, 0.000000, 0.400000, {0.000000, 0.539550, 0.170750}, {0.000000, 0.539550, 0.170750}, {1.000000, 0.000000, 0.000000}},
    {1.0e+06, 353.15, {0.400000, 0.001820, 0.323730, 0.069380, 0.053220, 0.027610, 0.012650, 0.009140, 0.102450}, 0.145660, 0.454340, 0.400000, {0.000000, 0.013515, 0.690461}, {0.000000, 0.708195, 0.004133}, {1.000000, 0.000000, 0.000000}},
    {1.0e+07, 353.15, {0.400000, 0.001820, 0.323730, 0.069380, 0.053220, 0.027610, 0.012650, 0.009140, 0.102450}, 0.261401, 0.338599, 0.400000, {0.000000, 0.154804, 0.389693}, {0.000000, 0.836577, 0.001724}, {1.000000, 0.000000, 0.000000}},
    {1.0e+08, 353.15, {0.400000, 0.001820, 0.323730, 0.069380, 0.053220, 0.027610, 0.012650, 0.009140, 0.102450}, 0.600000, 0.000000, 0.400000, {0.000000, 0.539550, 0.170750}, {0.000000, 0.539550, 0.170750}, {1.000000, 0.000000, 0.000000}},
    {1.0e+06, 313.15, {1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000}, 0.000000, 0.000000, 1.000000, {0.000000, 0.125000, 0.125000}, {0.000000, 0.125000, 0.125000}, {1.000000, 0.000000, 0.000000}},
    {1.0e+08, 353.15, {1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000}, 0.000000, 0.000000, 1.000000, {0.000000, 0.125000, 0.125000}, {0.000000, 0.125000, 0.125000}, {1.000000, 0.000000, 0.000000}}
  })
);

/* UNCRUSTIFY-ON */

} // testing

} // geos
