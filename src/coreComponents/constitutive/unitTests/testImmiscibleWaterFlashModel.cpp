/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

// Source includes
#include "codingUtilities/UnitTestUtilities.hpp"
#include "constitutive/fluid/multifluid/MultiFluidUtils.hpp"
#include "constitutive/fluid/multifluid/compositional/models/EquationOfState.hpp"
#include "constitutive/fluid/multifluid/compositional/models/CriticalVolume.hpp"
#include "constitutive/fluid/multifluid/compositional/models/ImmiscibleWaterFlashModel.hpp"
#include "TestFluid.hpp"
#include "TestFluidUtilities.hpp"

using namespace geos::constitutive;
using namespace geos::constitutive::compositional;

namespace geos
{
namespace testing
{

constexpr int numTestComps = 3;

template< int NC >
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

template< int NC >
struct FluidData {};

template<>
struct FluidData< 3 >
{
  static constexpr integer testComponents[numTestComps] = {0, 1, 2};
  static std::unique_ptr< TestFluid< 3 > > createFluid()
  {
    auto fluid = TestFluid< 3 >::create( {Fluid::C1, Fluid::C10, Fluid::H2O} );
    const std::array< real64 const, 3 > bics = { 0.25, 0.0, 0.0 };
    fluid->setBinaryCoefficients( bics );
    return fluid;
  }
};

template<>
struct FluidData< 9 >
{
  static constexpr integer testComponents[numTestComps] = {0, 1, 8};
  static std::unique_ptr< TestFluid< 9 > > createFluid()
  {
    auto fluid = TestFluid< 9 >::create( {Fluid::H2O, Fluid::CO2, Fluid::N2, Fluid::C5, Fluid::C2, Fluid::C3, Fluid::C4, Fluid::C5, Fluid::C10} );
    const std::array< real64 const, 36 > bics = {
      0.01, 0, 0.003732, 0, 0.01, 0, 0, 0.01, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0.01, 0, 0.028, 0.01, 0.01, 0, 0, 0.01, 0, 0.04532, 0.01, 0.01, 0, 0, 0
    };
    fluid->setBinaryCoefficients( bics );
    return fluid;
  }
};

template< int NC >
class ImmiscibleWaterFlashModelTestFixture :  public ::testing::TestWithParam< FlashData< NC > >
{
  static constexpr real64 relTol = 1.0e-5;
  static constexpr real64 absTol = 1.0e-7;
  static constexpr int numPhases = 3;
  static constexpr int numComps = NC;
  static constexpr int numDofs = NC + 2;
  using Deriv = geos::constitutive::multifluid::DerivativeOffset;
  using PhasePropSlice = ImmiscibleWaterFlashModelUpdate::PhaseProp::SliceType;
  using PhaseCompSlice = ImmiscibleWaterFlashModelUpdate::PhaseComp::SliceType;

public:
  ImmiscibleWaterFlashModelTestFixture()
    : m_fluid( FluidData< NC >::createFluid() )
  {
    ComponentProperties const & componentProperties = this->m_fluid->getComponentProperties();

    m_parameters = ImmiscibleWaterFlashModel::createParameters( std::move( m_parameters ) );

    auto * parameters = const_cast< CriticalVolume * >(m_parameters->get< CriticalVolume >());
    parameters->m_componentCriticalVolume.resize( NC );
    TestFluid< 9 >::populateArray( parameters->m_componentCriticalVolume, this->m_fluid->criticalVolume );

    auto * equationOfState = const_cast< EquationOfState * >(m_parameters->get< EquationOfState >());
    string const eosName = EnumStrings< EquationOfStateType >::toString( EquationOfStateType::PengRobinson );
    equationOfState->m_equationsOfStateNames.emplace_back( eosName );
    equationOfState->m_equationsOfStateNames.emplace_back( eosName );
    string const waterEosName = EnumStrings< EquationOfStateType >::toString( EquationOfStateType::ImmiscibleWater );
    equationOfState->m_equationsOfStateNames.emplace_back( waterEosName );

    m_flash = std::make_unique< ImmiscibleWaterFlashModel >( "FlashModel", componentProperties, *m_parameters );
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

    stackArray1d< real64, numPhases > phaseFraction( numPhases );
    stackArray2d< real64, numPhases *numDofs > dPhaseFraction( numPhases, numDofs );
    stackArray2d< real64, numPhases *numComps > phaseComponentFraction( numPhases, numComps );
    stackArray3d< real64, numPhases *numComps *numDofs > dPhaseComponentFraction( numPhases, numComps, numDofs );

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
      for( integer const ic : FluidData< numComps >::testComponents )
      {
        checkRelativeError( phaseComponentFraction[ip][ic], expectedPhaseComponentFraction[ip][ic], relTol, absTol );
      }
    }

/**
((void)expectedPhaseFraction);
((void)expectedPhaseComponentFraction);

    std::cout << "FlashData<" << numComps << ">{";
    std::cout << std::scientific << std::setprecision( 1 ) << pressure
              << ", "
              << std::fixed << std::setprecision( 2 ) << temperature
              << ", {";
    std::cout << std::fixed << std::setprecision( 6 );
    std::cout << composition[0];
    for( integer ic = 1; ic < numComps; ic++ )
    {
      std::cout << ", " << composition[ic];
    }
    std::cout << "}, "
              << phaseFraction[0] << ", "
              << phaseFraction[1] << ", "
              << phaseFraction[2];
    for( integer phaseIndex = 0; phaseIndex < 3; phaseIndex++ )
    {
      std::cout << ", {" << phaseComponentFraction[phaseIndex][FluidData< numComps >::testComponents[0]];
      for( integer ic = 1; ic < numTestComps; ic++ )
      {
        std::cout << ", " << phaseComponentFraction[phaseIndex][FluidData< numComps >::testComponents[ic]];
      }
      std::cout << "}";
    }
    std::cout << "},\n";
*/
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

    stackArray1d< real64, numPhases > phaseFraction( numPhases );
    stackArray2d< real64, numPhases *numDofs > dPhaseFraction( numPhases, numDofs );
    stackArray2d< real64, numPhases *numComps > phaseComponentFraction( numPhases, numComps );
    stackArray3d< real64, numPhases *numComps *numDofs > dPhaseComponentFraction( numPhases, numComps, numDofs );
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
      stackArray1d< real64, numPhases > displacedPhaseFraction( numPhases );
      stackArray2d< real64, numPhases *numDofs > displacedPhaseFractionDerivs( numPhases, numDofs );
      stackArray2d< real64, numPhases *numComps > displacedPhaseComponentFraction( numPhases, numComps );
      stackArray3d< real64, numPhases *numComps *numDofs > displacedPhaseComponentFractionDerivs( numPhases, numComps, numDofs );

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

std::vector< FlashData< 9 > > generate()
{
  std::vector< FlashData< 9 > > data;
  std::array<real64, 3> zero = {0,0,0};
  std::vector<std::array<real64, 9>> comps{
    {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200},
    {0.40000,	0.00182,	0.32373,	0.06938,	0.05322,	0.02761,	0.01265,	0.00914,	0.10245},
    {1,0,0,0,0,0,0,0,0}
  };
  for (std::array<real64, 9> const & comp : comps)
  {
    for (real64 const t : {20.0, 40.0, 80.0})
    {
  for (real64 const p : {10.0, 100.0, 1000.0})
  {
data.emplace_back(FlashData< 9 >{1.0e5*p, t+273.15, comp, 0.0, 0.0, 0.0, zero, zero, zero});
    }
  }
  }
  return data;
}

INSTANTIATE_TEST_SUITE_P(
  ImmiscibleWaterFlashModel, ImmiscibleWaterFlashModel3,
  ::testing::Values(
    FlashData<3>{1.0e+05, 293.15, {0.000000, 0.000000, 1.000000}, 0.000000, 0.000000, 1.000000, {0.500000, 0.500000, 0.000000}, {0.500000, 0.500000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+06, 293.15, {0.000000, 0.000000, 1.000000}, 0.000000, 0.000000, 1.000000, {0.500000, 0.500000, 0.000000}, {0.500000, 0.500000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+07, 293.15, {0.000000, 0.000000, 1.000000}, 0.000000, 0.000000, 1.000000, {0.500000, 0.500000, 0.000000}, {0.500000, 0.500000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+05, 313.15, {0.000000, 0.000000, 1.000000}, 0.000000, 0.000000, 1.000000, {0.500000, 0.500000, 0.000000}, {0.500000, 0.500000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+06, 313.15, {0.000000, 0.000000, 1.000000}, 0.000000, 0.000000, 1.000000, {0.500000, 0.500000, 0.000000}, {0.500000, 0.500000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+07, 313.15, {0.000000, 0.000000, 1.000000}, 0.000000, 0.000000, 1.000000, {0.500000, 0.500000, 0.000000}, {0.500000, 0.500000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+05, 353.15, {0.000000, 0.000000, 1.000000}, 0.000000, 0.000000, 1.000000, {0.500000, 0.500000, 0.000000}, {0.500000, 0.500000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+06, 353.15, {0.000000, 0.000000, 1.000000}, 0.000000, 0.000000, 1.000000, {0.500000, 0.500000, 0.000000}, {0.500000, 0.500000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+07, 353.15, {0.000000, 0.000000, 1.000000}, 0.000000, 0.000000, 1.000000, {0.500000, 0.500000, 0.000000}, {0.500000, 0.500000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+05, 293.15, {0.300000, 0.300000, 0.400000}, 0.300262, 0.299738, 0.400000, {0.000874, 0.999126, 0.000000}, {1.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+06, 293.15, {0.300000, 0.300000, 0.400000}, 0.302619, 0.297381, 0.400000, {0.008653, 0.991347, 0.000000}, {1.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+07, 293.15, {0.300000, 0.300000, 0.400000}, 0.325574, 0.274426, 0.400000, {0.078551, 0.921449, 0.000000}, {1.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+05, 313.15, {0.300000, 0.300000, 0.400000}, 0.300281, 0.299719, 0.400000, {0.000935, 0.999065, 0.000000}, {1.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+06, 313.15, {0.300000, 0.300000, 0.400000}, 0.302802, 0.297198, 0.400000, {0.009255, 0.990745, 0.000000}, {1.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+07, 313.15, {0.300000, 0.300000, 0.400000}, 0.327540, 0.272460, 0.400000, {0.084080, 0.915920, 0.000000}, {1.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+05, 353.15, {0.300000, 0.300000, 0.400000}, 0.300314, 0.299686, 0.400000, {0.001062, 0.998938, 0.000000}, {0.999982, 0.000018, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+06, 353.15, {0.300000, 0.300000, 0.400000}, 0.303187, 0.296813, 0.400000, {0.010514, 0.989486, 0.000000}, {0.999998, 0.000002, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+07, 353.15, {0.300000, 0.300000, 0.400000}, 0.331645, 0.268355, 0.400000, {0.095418, 0.904582, 0.000000}, {0.999999, 0.000001, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+05, 293.15, {0.200000, 0.800000, 0.000000}, 0.800700, 0.199300, 0.000000, {0.000874, 0.999126, 0.000000}, {1.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+06, 293.15, {0.200000, 0.800000, 0.000000}, 0.806983, 0.193017, 0.000000, {0.008653, 0.991347, 0.000000}, {1.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+07, 293.15, {0.200000, 0.800000, 0.000000}, 0.868198, 0.131802, 0.000000, {0.078551, 0.921449, 0.000000}, {1.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+05, 313.15, {0.200000, 0.800000, 0.000000}, 0.800749, 0.199251, 0.000000, {0.000935, 0.999065, 0.000000}, {1.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+06, 313.15, {0.200000, 0.800000, 0.000000}, 0.807473, 0.192527, 0.000000, {0.009255, 0.990745, 0.000000}, {1.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+07, 313.15, {0.200000, 0.800000, 0.000000}, 0.873439, 0.126561, 0.000000, {0.084080, 0.915920, 0.000000}, {1.000000, 0.000000, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+05, 353.15, {0.200000, 0.800000, 0.000000}, 0.800847, 0.199153, 0.000000, {0.001062, 0.998938, 0.000000}, {0.999982, 0.000018, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+06, 353.15, {0.200000, 0.800000, 0.000000}, 0.808500, 0.191500, 0.000000, {0.010514, 0.989486, 0.000000}, {0.999998, 0.000002, 0.000000}, {0.000000, 0.000000, 1.000000}},
    FlashData<3>{1.0e+07, 353.15, {0.200000, 0.800000, 0.000000}, 0.884386, 0.115614, 0.000000, {0.095418, 0.904582, 0.000000}, {0.999999, 0.000001, 0.000000}, {0.000000, 0.000000, 1.000000}}    
  )
);

INSTANTIATE_TEST_SUITE_P(
  ImmiscibleWaterFlashModel, ImmiscibleWaterFlashModel9,
  ::testing::Values(
    FlashData<9>{1.0e+06, 293.15, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 0.198828, 0.792172, 0.009000, {0.000000, 0.000862, 0.850986}, {0.000000, 0.003571, 0.000000}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+07, 293.15, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 0.307007, 0.683993, 0.009000, {0.000000, 0.003780, 0.551128}, {0.000000, 0.002689, 0.000000}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+08, 293.15, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 0.991000, 0.000000, 0.009000, {0.000000, 0.003027, 0.170737}, {0.000000, 0.003027, 0.170737}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+06, 313.15, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 0.190493, 0.800507, 0.009000, {0.000000, 0.000635, 0.888223}, {0.000000, 0.003596, 0.000000}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+07, 313.15, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 0.288620, 0.702380, 0.009000, {0.000000, 0.003255, 0.586237}, {0.000000, 0.002934, 0.000001}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+08, 313.15, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 0.991000, 0.000000, 0.009000, {0.000000, 0.003027, 0.170737}, {0.000000, 0.003027, 0.170737}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+06, 353.15, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 0.181787, 0.809213, 0.009000, {0.000000, 0.000395, 0.930750}, {0.000000, 0.003619, 0.000003}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+07, 353.15, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 0.261660, 0.729340, 0.009000, {0.000000, 0.002467, 0.646622}, {0.000000, 0.003228, 0.000007}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+08, 353.15, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 0.947424, 0.043576, 0.009000, {0.000000, 0.003062, 0.178339}, {0.000000, 0.002269, 0.005438}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+06, 293.15, {0.400000, 0.001820, 0.323730, 0.069380, 0.053220, 0.027610, 0.012650, 0.009140, 0.102450}, 0.120388, 0.479612, 0.400000, {0.000000, 0.000864, 0.850996}, {0.000000, 0.003578, 0.000000}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+07, 293.15, {0.400000, 0.001820, 0.323730, 0.069380, 0.053220, 0.027610, 0.012650, 0.009140, 0.102450}, 0.185888, 0.414112, 0.400000, {0.000000, 0.003787, 0.551139}, {0.000000, 0.002695, 0.000000}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+08, 293.15, {0.400000, 0.001820, 0.323730, 0.069380, 0.053220, 0.027610, 0.012650, 0.009140, 0.102450}, 0.600000, 0.000000, 0.400000, {0.000000, 0.003033, 0.170750}, {0.000000, 0.003033, 0.170750}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+06, 313.15, {0.400000, 0.001820, 0.323730, 0.069380, 0.053220, 0.027610, 0.012650, 0.009140, 0.102450}, 0.115342, 0.484658, 0.400000, {0.000000, 0.000637, 0.888229}, {0.000000, 0.003604, 0.000000}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+07, 313.15, {0.400000, 0.001820, 0.323730, 0.069380, 0.053220, 0.027610, 0.012650, 0.009140, 0.102450}, 0.174755, 0.425245, 0.400000, {0.000000, 0.003261, 0.586247}, {0.000000, 0.002940, 0.000001}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+08, 313.15, {0.400000, 0.001820, 0.323730, 0.069380, 0.053220, 0.027610, 0.012650, 0.009140, 0.102450}, 0.600000, 0.000000, 0.400000, {0.000000, 0.003033, 0.170750}, {0.000000, 0.003033, 0.170750}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+06, 353.15, {0.400000, 0.001820, 0.323730, 0.069380, 0.053220, 0.027610, 0.012650, 0.009140, 0.102450}, 0.110071, 0.489929, 0.400000, {0.000000, 0.000396, 0.930753}, {0.000000, 0.003626, 0.000003}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+07, 353.15, {0.400000, 0.001820, 0.323730, 0.069380, 0.053220, 0.027610, 0.012650, 0.009140, 0.102450}, 0.158432, 0.441568, 0.400000, {0.000000, 0.002472, 0.646630}, {0.000000, 0.003235, 0.000007}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+08, 353.15, {0.400000, 0.001820, 0.323730, 0.069380, 0.053220, 0.027610, 0.012650, 0.009140, 0.102450}, 0.573644, 0.026356, 0.400000, {0.000000, 0.003068, 0.178345}, {0.000000, 0.002273, 0.005437}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+06, 313.15, {1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000}, 0.000000, 0.000000, 1.000000, {0.000000, 0.125000, 0.125000}, {0.000000, 0.125000, 0.125000}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+08, 353.15, {1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000}, 0.000000, 0.000000, 1.000000, {0.000000, 0.125000, 0.125000}, {0.000000, 0.125000, 0.125000}, {1.000000, 0.000000, 0.000000}}   
  )
);

/* UNCRUSTIFY-ON */

} // testing

} // geos
