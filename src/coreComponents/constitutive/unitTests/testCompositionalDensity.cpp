/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

// Source includes
#include "codingUtilities/UnitTestUtilities.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/CubicEOSPhaseModel.hpp"
#include "constitutive/fluid/multifluid/compositional/models/CompositionalDensity.hpp"
#include "TestFluid.hpp"
#include "TestFluidUtilities.hpp"

using namespace geos::constitutive::compositional;

namespace geos
{
namespace testing
{

template< int NC >
using DensityData = std::tuple<
  real64 const,       // pressure
  real64 const,       // temperature
  Feed< NC > const,   // phase composition
  real64 const,       // expected molar density
  real64 const        // expected mass density
  >;

template< int NC >
struct FluidData {};

template<>
struct FluidData< 9 >
{
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
  static std::vector< DensityData< 9 > > generateTestData()
  {
    std::array< real64 const, 2 > pressures( {1.83959e+06, 1.83959e+08} );
    std::array< real64 const, 2 > temperatures( {2.97150e+02, 3.63000e+02} );
    std::array< Feed< 9 > const, 3 > feeds( {
          Feed< 9 >{0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920},
          Feed< 9 >{0.00826, 0.00544, 0.77032, 0.10456, 0.06177, 0.02459, 0.00884, 0.00472, 0.01149},
          Feed< 9 >{0.00899, 0.00299, 0.53281, 0.11447, 0.08791, 0.04566, 0.02095, 0.01516, 0.17107}
        } );
    std::vector< DensityData< 9 > > testData;
    for( const real64 pressure : pressures )
    {
      for( const real64 temperature : temperatures )
      {
        for( const auto & composition : feeds )
        {
          testData.emplace_back( pressure, temperature, composition, 0.0, 0.0 );
        }
      }
    }
    return testData;
  }
};

template< int NC, typename EOS_TYPE >
class CompositionalDensityTestFixture :  public ::testing::TestWithParam< DensityData< NC > >
{
  static constexpr real64 relTol = 1.0e-5;
  static constexpr real64 absTol = 1.0e-7;
  static constexpr int numComps = NC;
  static constexpr int numDofs = NC + 2;
  using Deriv = geos::constitutive::multifluid::DerivativeOffset;
public:
  CompositionalDensityTestFixture()
    : m_fluid( FluidData< NC >::createFluid() )
  {
    ComponentProperties const & componentProperties = this->m_fluid->getComponentProperties();
    m_density = std::make_unique< CompositionalDensity< EOS_TYPE > >( "PhaseDensity", componentProperties );
  }

  ~CompositionalDensityTestFixture() = default;

  void testDensityValues( DensityData< NC > const & data )
  {
    real64 const pressure = std::get< 0 >( data );
    real64 const temperature = std::get< 1 >( data );
    stackArray1d< real64, numComps > composition;
    TestFluid< NC >::createArray( composition, std::get< 2 >( data ));
    real64 const expectedMolarDensity = std::get< 3 >( data );
    real64 const expectedMassDensity = std::get< 4 >( data );

    auto componentProperties = m_fluid->createKernelWrapper();
    auto kernelWrapper = this->m_density->createKernelWrapper();

    real64 molarDensity = 0.0;
    real64 massDensity = 0.0;
    stackArray2d< real64, numComps *numDofs > dComposition( numComps, numDofs );
    stackArray1d< real64, numDofs > tempDerivs( numDofs );

    kernelWrapper.compute( componentProperties,
                           pressure,
                           temperature,
                           composition,
                           dComposition,
                           molarDensity,
                           tempDerivs,
                           massDensity,
                           tempDerivs,
                           false );

    checkRelativeError( molarDensity, expectedMolarDensity, relTol, absTol );
    checkRelativeError( massDensity, expectedMassDensity, relTol, absTol );
  }

  void testDensityDerivatives( DensityData< NC > const & data )
  {
    real64 const pressure = std::get< 0 >( data );
    real64 const temperature = std::get< 1 >( data );
    stackArray1d< real64, numComps > composition;
    TestFluid< NC >::createArray( composition, std::get< 2 >( data ));

    auto componentProperties = m_fluid->createKernelWrapper();
    auto kernelWrapper = this->m_density->createKernelWrapper();

    real64 molarDensity = 0.0;
    real64 massDensity = 0.0;
    stackArray1d< real64, numDofs > molarDensityDerivs( numDofs );
    stackArray1d< real64, numDofs > massDensityDerivs( numDofs );
    stackArray1d< real64, numComps > liquidComposition( numComps );
    stackArray2d< real64, numComps *numDofs > dLiquidComposition( numComps, numDofs );
    calculatePhaseComposition( pressure,
                               temperature,
                               composition,
                               liquidComposition,
                               dLiquidComposition );
    kernelWrapper.compute( componentProperties,
                           pressure,
                           temperature,
                           liquidComposition,
                           dLiquidComposition,
                           molarDensity,
                           molarDensityDerivs,
                           massDensity,
                           massDensityDerivs,
                           false );

    auto calculateDensity = [&]( real64 const p, real64 const t, auto const & zmf ) -> std::pair< real64, real64 > {
      stackArray1d< real64, numComps > xmf( numComps );
      stackArray2d< real64, numComps *numDofs > dxmf( numComps, numDofs );

      this->calculatePhaseComposition( p, t, zmf, xmf, dxmf );

      real64 densityMolar = 0.0;
      real64 densityMass = 0.0;
      stackArray1d< real64, numDofs > tempDerivs( numDofs );
      kernelWrapper.compute( componentProperties, p, t,
                             xmf, dxmf,
                             densityMolar, tempDerivs, densityMass, tempDerivs, false );
      return {densityMolar, densityMass};
    };

    // Compare against numerical derivatives
    // -- Pressure derivative
    real64 const dp = 1.0e-4 * pressure;
    internal::testNumericalDerivative(
      pressure, dp, molarDensityDerivs[Deriv::dP],
      [&]( real64 const p ) -> real64 {
      return calculateDensity( p, temperature, composition ).first;
    } );
    internal::testNumericalDerivative(
      pressure, dp, massDensityDerivs[Deriv::dP],
      [&]( real64 const p ) -> real64 {
      return calculateDensity( p, temperature, composition ).second;
    } );

    // -- Temperature derivative
    real64 const dT = 1.0e-6 * temperature;
    internal::testNumericalDerivative(
      temperature, dT, molarDensityDerivs[Deriv::dT],
      [&]( real64 const t ) -> real64 {
      return calculateDensity( pressure, t, composition ).first;
    } );
    internal::testNumericalDerivative(
      temperature, dT, massDensityDerivs[Deriv::dT],
      [&]( real64 const t ) -> real64 {
      return calculateDensity( pressure, t, composition ).second;
    } );

    // -- Composition derivatives derivative
    real64 const dz = 1.0e-7;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      internal::testNumericalDerivative(
        0.0, dz, molarDensityDerivs[Deriv::dC + ic],
        [&]( real64 const z ) -> real64 {
        stackArray1d< real64, numComps > zmf( numComps );
        for( integer jc = 0; jc < numComps; ++jc )
        {
          zmf[jc] = composition[jc];
        }
        zmf[ic] += z;
        return calculateDensity( pressure, temperature, zmf ).first;
      } );
      internal::testNumericalDerivative(
        0.0, dz, massDensityDerivs[Deriv::dC + ic],
        [&]( real64 const z ) -> real64 {
        stackArray1d< real64, numComps > zmf( numComps );
        for( integer jc = 0; jc < numComps; ++jc )
        {
          zmf[jc] = composition[jc];
        }
        zmf[ic] += z;
        return calculateDensity( pressure, temperature, zmf ).second;
      } );
    }
  }

private:
  void calculatePhaseComposition( real64 const pressure,
                                  real64 const temperature,
                                  arraySlice1d< real64 const > const & totalComposition,
                                  arraySlice1d< real64 > const & phaseComposition,
                                  arraySlice2d< real64 > const & dPhaseComposition ) const;

protected:
  std::unique_ptr< CompositionalDensity< EOS_TYPE > > m_density{};
  std::unique_ptr< TestFluid< NC > > m_fluid{};
};

template< int NC, typename EOS_TYPE >
void
CompositionalDensityTestFixture< NC, EOS_TYPE >::
calculatePhaseComposition( real64 const pressure,
                           real64 const temperature,
                           arraySlice1d< real64 const > const & totalComposition,
                           arraySlice1d< real64 > const & phaseComposition,
                           arraySlice2d< real64 > const & dPhaseComposition ) const
{
  real64 constexpr V = 0.25;
  real64 constexpr PRef = 1.5e7;
  real64 constexpr TRef = 323.15;
  real64 constexpr cP = 5.0e-10;
  real64 constexpr cT = 6.0e-3;
  real64 sum = 0.0;
  stackArray1d< real64, numDofs > sumDerivs( numDofs );
  sumDerivs.zero();
  for( integer ic = 0; ic < numComps; ic++ )
  {
    real64 const d = 1.0 - 2.0*(ic + 0.5)/numComps;
    real64 const k = exp( cT*d*(temperature + TRef))*exp( cP*d*(pressure - PRef));
    real64 const dk_dP = cP*d*k;
    real64 const dk_dT = cT*d*k;
    real64 const m = 1.0 / (1.0 - V + k*V);
    real64 const xi = totalComposition[ic] * m;
    phaseComposition[ic] = xi;
    dPhaseComposition( ic, Deriv::dP ) = -totalComposition[ic]*V*dk_dP*m*m;
    dPhaseComposition( ic, Deriv::dT ) = -totalComposition[ic]*V*dk_dT*m*m;
    dPhaseComposition( ic, Deriv::dC+ic ) = m;

    sum += xi;
    for( integer kc = 0; kc < numDofs; kc++ )
    {
      sumDerivs[kc] += dPhaseComposition( ic, kc );
    }
  }
  real64 const oneOverSum = 1.0 / sum;
  for( integer ic = 0; ic < numComps; ic++ )
  {
    real64 const xi = phaseComposition[ic];
    phaseComposition[ic] *= oneOverSum;
    for( integer kc = 0; kc < numDofs; kc++ )
    {
      dPhaseComposition( ic, kc ) = (dPhaseComposition( ic, kc )*sum - sumDerivs[kc]*xi)*oneOverSum*oneOverSum;
    }
  }
}


using CompositionalDensity9CompPR = CompositionalDensityTestFixture< 9, CubicEOSPhaseModel< PengRobinsonEOS > >;
using CompositionalDensity9CompSRK = CompositionalDensityTestFixture< 9, CubicEOSPhaseModel< SoaveRedlichKwongEOS > >;

TEST_P( CompositionalDensity9CompPR, testDensityDerivatives )
{
  testDensityDerivatives( GetParam() );
}

TEST_P( CompositionalDensity9CompSRK, testDensityDerivatives )
{
  testDensityDerivatives( GetParam() );
}

INSTANTIATE_TEST_SUITE_P(
  CompositionalDensityTest,
  CompositionalDensity9CompPR,
  ::testing::ValuesIn( {
      DensityData< 9 >{1.839590e+06, 2.971500e+02, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 8.355571e+03, 4.559906e+02 },
      DensityData< 9 >{1.839590e+06, 2.971500e+02, {0.008260, 0.005440, 0.770320, 0.104560, 0.061770, 0.024590, 0.008840, 0.004720, 0.011490}, 7.703898e+02, 2.691914e+01 },
      DensityData< 9 >{1.839590e+06, 2.971500e+02, {0.008990, 0.002990, 0.532810, 0.114470, 0.087910, 0.045660, 0.020950, 0.015160, 0.171070}, 8.337694e+03, 4.567935e+02 },
      DensityData< 9 >{1.839590e+06, 3.630000e+02, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 9.073321e+02, 4.951606e+01 },
      DensityData< 9 >{1.839590e+06, 3.630000e+02, {0.008260, 0.005440, 0.770320, 0.104560, 0.061770, 0.024590, 0.008840, 0.004720, 0.011490}, 6.178234e+02, 2.158813e+01 },
      DensityData< 9 >{1.839590e+06, 3.630000e+02, {0.008990, 0.002990, 0.532810, 0.114470, 0.087910, 0.045660, 0.020950, 0.015160, 0.171070}, 9.197865e+02, 5.039192e+01 },
      DensityData< 9 >{1.839590e+08, 2.971500e+02, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 1.095078e+04, 5.976195e+02 },
      DensityData< 9 >{1.839590e+08, 2.971500e+02, {0.008260, 0.005440, 0.770320, 0.104560, 0.061770, 0.024590, 0.008840, 0.004720, 0.011490}, 2.480270e+04, 8.666618e+02 },
      DensityData< 9 >{1.839590e+08, 2.971500e+02, {0.008990, 0.002990, 0.532810, 0.114470, 0.087910, 0.045660, 0.020950, 0.015160, 0.171070}, 1.087917e+04, 5.960323e+02 },
      DensityData< 9 >{1.839590e+08, 3.630000e+02, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 1.065128e+04, 5.812747e+02 },
      DensityData< 9 >{1.839590e+08, 3.630000e+02, {0.008260, 0.005440, 0.770320, 0.104560, 0.061770, 0.024590, 0.008840, 0.004720, 0.011490}, 2.305823e+04, 8.057060e+02 },
      DensityData< 9 >{1.839590e+08, 3.630000e+02, {0.008990, 0.002990, 0.532810, 0.114470, 0.087910, 0.045660, 0.020950, 0.015160, 0.171070}, 1.058381e+04, 5.798506e+02 }
    } )
  );
INSTANTIATE_TEST_SUITE_P(
  CompositionalDensityTest,
  CompositionalDensity9CompSRK,
  ::testing::ValuesIn( {
      DensityData< 9 >{1.839590e+06, 2.971500e+02, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 7.433979e+03, 4.056963e+02 },
      DensityData< 9 >{1.839590e+06, 2.971500e+02, {0.008260, 0.005440, 0.770320, 0.104560, 0.061770, 0.024590, 0.008840, 0.004720, 0.011490}, 7.629968e+02, 2.666082e+01 },
      DensityData< 9 >{1.839590e+06, 2.971500e+02, {0.008990, 0.002990, 0.532810, 0.114470, 0.087910, 0.045660, 0.020950, 0.015160, 0.171070}, 7.416959e+03, 4.063495e+02 },
      DensityData< 9 >{1.839590e+06, 3.630000e+02, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 8.919848e+02, 4.867851e+01 },
      DensityData< 9 >{1.839590e+06, 3.630000e+02, {0.008260, 0.005440, 0.770320, 0.104560, 0.061770, 0.024590, 0.008840, 0.004720, 0.011490}, 6.133569e+02, 2.143206e+01 },
      DensityData< 9 >{1.839590e+06, 3.630000e+02, {0.008990, 0.002990, 0.532810, 0.114470, 0.087910, 0.045660, 0.020950, 0.015160, 0.171070}, 9.045641e+02, 4.955794e+01 },
      DensityData< 9 >{1.839590e+08, 2.971500e+02, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 9.868675e+03, 5.385656e+02 },
      DensityData< 9 >{1.839590e+08, 2.971500e+02, {0.008260, 0.005440, 0.770320, 0.104560, 0.061770, 0.024590, 0.008840, 0.004720, 0.011490}, 2.257420e+04, 7.887929e+02 },
      DensityData< 9 >{1.839590e+08, 2.971500e+02, {0.008990, 0.002990, 0.532810, 0.114470, 0.087910, 0.045660, 0.020950, 0.015160, 0.171070}, 9.803814e+03, 5.371171e+02 },
      DensityData< 9 >{1.839590e+08, 3.630000e+02, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}, 9.615674e+03, 5.247585e+02 },
      DensityData< 9 >{1.839590e+08, 3.630000e+02, {0.008260, 0.005440, 0.770320, 0.104560, 0.061770, 0.024590, 0.008840, 0.004720, 0.011490}, 2.109248e+04, 7.370183e+02 },
      DensityData< 9 >{1.839590e+08, 3.630000e+02, {0.008990, 0.002990, 0.532810, 0.114470, 0.087910, 0.045660, 0.020950, 0.015160, 0.171070}, 9.554300e+03, 5.234471e+02 }
    } )
  );

} // testing

} // geos
