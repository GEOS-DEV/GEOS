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

template< integer NC >
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

public:
  NegativeTwoPhaseFlashModelTestFixture()
    : m_fluid( FluidData< NC >::createFluid() )
  {
    ComponentProperties const & componentProperties = this->m_fluid->getComponentProperties();

    m_parameters = NegativeTwoPhaseFlashModel::createParameters( std::move( m_parameters ) );

    auto * equationOfState = const_cast< EquationOfState * >(m_parameters->get< EquationOfState >());
    string const eosName = EnumStrings< EquationOfStateType >::toString( EquationOfStateType::PengRobinson );
    equationOfState->m_equationsOfStateNames.emplace_back( eosName );
    equationOfState->m_equationsOfStateNames.emplace_back( eosName );

    auto * criticalVolume = const_cast< CriticalVolume * >(m_parameters->get< CriticalVolume >());
    criticalVolume->m_componentCriticalVolume.resize( numComps );
    TestFluid< NC >::createArray( criticalVolume->m_componentCriticalVolume, this->m_fluid->criticalVolume );

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
    std::cout << std::fixed << std::setprecision( 8 );
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


/**

   Z: { 0.10900000, 0.00300000, 0.23470000, 0.11460000, 0.08790000, 0.04560000, 0.02090000, 0.01510000, 0.36920000 }
   X: { 0.15807941, 0.00337666, 0.07159711, 0.02384273, 0.10085780, 0.06146075, 0.02939718, 0.00974393, 0.54164443 }
   Y: { 0.00396181, 0.00219389, 0.58376764, 0.30883583, 0.06016813, 0.01165533, 0.00271459, 0.02656289, 0.00013988 }

 */
    std::cout << " Z: " << composition.toSliceConst() << std::endl;
    std::cout << " X: " << phaseComponentFraction[0] << std::endl;
    std::cout << " Y: " << phaseComponentFraction[1] << std::endl;
    std::cout << " V: " << phaseFraction.toSliceConst() << std::endl;
    std::cout << std::scientific;
    std::cout << "DV: " << dPhaseFraction[0] << std::endl;
    std::cout << "DV: " << dPhaseFraction[1] << std::endl;
/**
    // Combine derivatives into a single output
    auto const concatDerivatives = []( integer const kc, auto & derivs, auto const & phaseFractionDerivs, auto const &
       phaseComponentFractionDerivs ){
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
      StackArray< real64, 4, numPhases *numComps, multifluid::LAYOUT_PHASE_COMP > displacedPhaseComponentFractionData( 1, 1, numPhases,
         numComps );
      StackArray< real64, 5, numPhases *numComps *numDofs, multifluid::LAYOUT_PHASE_COMP_DC > displacedPhaseComponentFractionDerivsData( 1,
         1, numPhases, numComps, numDofs );

      auto displacedPhaseFraction = displacedPhaseFractionData[0][0];
      auto displacedPhaseFractionDerivs = displacedPhaseFractionDerivsData[0][0];
      auto displacedPhaseComponentFraction = displacedPhaseComponentFractionData[0][0];
      auto displacedPhaseComponentFractionDerivs = displacedPhaseComponentFractionDerivsData[0][0];

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
 */
  }

protected:
  std::unique_ptr< TestFluid< NC > > m_fluid{};
  std::unique_ptr< NegativeTwoPhaseFlashModel > m_flash{};
  std::unique_ptr< ModelParameters > m_parameters{};
  array1d< integer > m_phaseTypes{};
};

//using NegativeTwoPhaseFlashModel3 = NegativeTwoPhaseFlashModelTestFixture< 3 >;
using NegativeTwoPhaseFlashModel9 = NegativeTwoPhaseFlashModelTestFixture< 9 >;

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
TEST_P( NegativeTwoPhaseFlashModel9, testFlashDerivatives )
{
  testFlashDerivatives( GetParam() );
}

//-------------------------------------------------------------------------------
// Data
//-------------------------------------------------------------------------------
/* UNCRUSTIFY-OFF */

INSTANTIATE_TEST_SUITE_P(
  NegativeTwoPhaseFlashModel, NegativeTwoPhaseFlashModel9,
  ::testing::ValuesIn<FlashData<9>>({
    {1.0000e+07, 293.15, {0.10900000, 0.00300000, 0.23470000, 0.11460000, 0.08790000, 0.04560000, 0.02090000, 0.01510000, 0.36920000}, 0.185888, 0.414112, {0.000000, 0.111648, 0.551139}, {0.000000, 0.731628, 0.000000}},
    {1.5000e+07, 293.15, {0.13701652, 0.00320982, 0.10906376, 0.08659405, 0.09558042, 0.05469515, 0.02598050, 0.01902140, 0.46883838}, 0.185888, 0.414112, {0.000000, 0.111648, 0.551139}, {0.000000, 0.731628, 0.000000}},
    {1.0000e+07, 293.15, {0.13701652, 0.00320982, 0.10906376, 0.08659405, 0.09558042, 0.05469515, 0.02598050, 0.01902140, 0.46883838}, 0.185888, 0.414112, {0.000000, 0.111648, 0.551139}, {0.000000, 0.731628, 0.000000}},
    {0.9999e+07, 293.15, {0.13701652, 0.00320982, 0.10906376, 0.08659405, 0.09558042, 0.05469515, 0.02598050, 0.01902140, 0.46883838}, 0.185888, 0.414112, {0.000000, 0.111648, 0.551139}, {0.000000, 0.731628, 0.000000}}
    //FlashData<9>{1.0000e+07, 293.15, {0.00396181, 0.00219389, 0.58376764, 0.30883583, 0.06016813, 0.01165533, 0.00271459, 0.02656289, 0.00013988}, 0.185888, 0.414112, {0.000000, 0.111648, 0.551139}, {0.000000, 0.731628, 0.000000}}

        /**,
    FlashData<9>{1.0e+06, 293.15, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 0.198828, 0.792172, 0.009000, {0.000000, 0.011551, 0.850986}, {0.000000, 0.672081, 0.000000}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+07, 293.15, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 0.307007, 0.683993, 0.009000, {0.000000, 0.111647, 0.551128}, {0.000000, 0.731621, 0.000000}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+08, 293.15, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 0.991000, 0.000000, 0.009000, {0.000000, 0.539556, 0.170737}, {0.000000, 0.539556, 0.170737}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+06, 313.15, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 0.190493, 0.800507, 0.009000, {0.000000, 0.010985, 0.888223}, {0.000000, 0.665337, 0.000000}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+07, 313.15, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 0.288620, 0.702380, 0.009000, {0.000000, 0.106713, 0.586237}, {0.000000, 0.717418, 0.000001}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+08, 313.15, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 0.991000, 0.000000, 0.009000, {0.000000, 0.539556, 0.170737}, {0.000000, 0.539556, 0.170737}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+06, 353.15, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 0.181787, 0.809213, 0.009000, {0.000000, 0.010453, 0.930750}, {0.000000, 0.658417, 0.000003}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+07, 353.15, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 0.261660, 0.729340, 0.009000, {0.000000, 0.101277, 0.646622}, {0.000000, 0.696794, 0.000007}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+08, 353.15, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 0.947424, 0.043576, 0.009000, {0.000000, 0.529993, 0.178339}, {0.000000, 0.747468, 0.005438}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+06, 293.15, {0.400000, 0.001820, 0.323730, 0.069380, 0.053220, 0.027610, 0.012650, 0.009140, 0.102450}, 0.120388, 0.479612, 0.400000, {0.000000, 0.011551, 0.850996}, {0.000000, 0.672084, 0.000000}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+08, 293.15, {0.400000, 0.001820, 0.323730, 0.069380, 0.053220, 0.027610, 0.012650, 0.009140, 0.102450}, 0.600000, 0.000000, 0.400000, {0.000000, 0.539550, 0.170750}, {0.000000, 0.539550, 0.170750}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+06, 313.15, {0.400000, 0.001820, 0.323730, 0.069380, 0.053220, 0.027610, 0.012650, 0.009140, 0.102450}, 0.115342, 0.484658, 0.400000, {0.000000, 0.010985, 0.888229}, {0.000000, 0.665341, 0.000000}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+07, 313.15, {0.400000, 0.001820, 0.323730, 0.069380, 0.053220, 0.027610, 0.012650, 0.009140, 0.102450}, 0.174755, 0.425245, 0.400000, {0.000000, 0.106714, 0.586247}, {0.000000, 0.717425, 0.000001}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+08, 313.15, {0.400000, 0.001820, 0.323730, 0.069380, 0.053220, 0.027610, 0.012650, 0.009140, 0.102450}, 0.600000, 0.000000, 0.400000, {0.000000, 0.539550, 0.170750}, {0.000000, 0.539550, 0.170750}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+06, 353.15, {0.400000, 0.001820, 0.323730, 0.069380, 0.053220, 0.027610, 0.012650, 0.009140, 0.102450}, 0.110071, 0.489929, 0.400000, {0.000000, 0.010453, 0.930753}, {0.000000, 0.658421, 0.000003}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+07, 353.15, {0.400000, 0.001820, 0.323730, 0.069380, 0.053220, 0.027610, 0.012650, 0.009140, 0.102450}, 0.158432, 0.441568, 0.400000, {0.000000, 0.101278, 0.646630}, {0.000000, 0.696800, 0.000007}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+08, 353.15, {0.400000, 0.001820, 0.323730, 0.069380, 0.053220, 0.027610, 0.012650, 0.009140, 0.102450}, 0.573644, 0.026356, 0.400000, {0.000000, 0.529997, 0.178345}, {0.000000, 0.747479, 0.005437}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+06, 313.15, {1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000}, 0.000000, 0.000000, 1.000000, {0.000000, 0.125000, 0.125000}, {0.000000, 0.125000, 0.125000}, {1.000000, 0.000000, 0.000000}},
    FlashData<9>{1.0e+08, 353.15, {1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000}, 0.000000, 0.000000, 1.000000, {0.000000, 0.125000, 0.125000}, {0.000000, 0.125000, 0.125000}, {1.000000, 0.000000, 0.000000}}*/
  })
);

/* UNCRUSTIFY-ON */

} // testing

} // geos
