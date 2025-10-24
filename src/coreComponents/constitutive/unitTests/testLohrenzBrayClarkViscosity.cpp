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
#include "constitutive/fluid/multifluid/compositional/models/CompositionalDensity.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/CriticalVolume.hpp"
#include "constitutive/fluid/multifluid/compositional/models/LohrenzBrayClarkViscosity.hpp"
#include "TestFluid.hpp"
#include "TestFluidUtilities.hpp"

using namespace geos::constitutive::compositional;

namespace geos
{
namespace testing
{

template< int NC >
using ViscosityData = std::tuple<
  integer const,            // Mixing type
  real64 const,             // pressure
  real64 const,             // temperature
  Feed< NC > const,         // phase composition
  real64 const              // expected viscosity
  >;

template< int NC >
struct FluidData {};

template<>
struct FluidData< 9 >
{
  static std::unique_ptr< TestFluid< 9 > > createFluid()
  {
    auto fluid = TestFluid< 9 >::create( {Fluid::H2O, Fluid::CO2, Fluid::N2, Fluid::CH4, Fluid::C2H6, Fluid::C3H8, Fluid::C4H10, Fluid::C5H12, Fluid::C10H22} );
    const stdArray< real64 const, 36 > bics = {
      0.01, 0, 0.003732, 0, 0.01, 0, 0, 0.01, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0.01, 0, 0.028, 0.01, 0.01, 0, 0, 0.01, 0, 0.04532, 0.01, 0.01, 0, 0, 0
    };
    fluid->setBinaryCoefficients( bics );
    return fluid;
  }
};

template< int NC >
class LohrenzBrayClarkViscosityTestFixture :  public ::testing::TestWithParam< ViscosityData< NC > >
{
  static constexpr real64 relTol = 1.0e-5;
  static constexpr real64 absTol = 1.0e-7;
  static constexpr int numComps = NC;
  static constexpr int numDofs = NC + 2;
  using Deriv = geos::constitutive::multifluid::DerivativeOffset;
public:
  LohrenzBrayClarkViscosityTestFixture()
    : m_fluid( FluidData< NC >::createFluid() )
  {
    ComponentProperties const & componentProperties = this->m_fluid->getComponentProperties();

    m_parameters = CompositionalDensity::createParameters( std::make_unique< ModelParameters >() );
    m_parameters = LohrenzBrayClarkViscosity::createParameters( std::move( m_parameters ) );

    auto * parameters = const_cast< CriticalVolume * >(m_parameters->get< CriticalVolume >());
    parameters->m_componentCriticalVolume.resize( NC );
    TestFluid< 9 >::populateArray( parameters->m_componentCriticalVolume, this->m_fluid->criticalVolume );

    auto * equationOfState = const_cast< EquationOfState * >(m_parameters->get< EquationOfState >());
    string const eosName = EnumStrings< EquationOfStateType >::toString( EquationOfStateType::PengRobinson );
    equationOfState->m_equationsOfStateNames.emplace_back( eosName );

    m_density = std::make_unique< CompositionalDensity >( "PhaseDensity", componentProperties, 0, *m_parameters );
    m_viscosity = std::make_unique< LohrenzBrayClarkViscosity >( "PhaseViscosity", componentProperties, 0, *m_parameters );
  }

  ~LohrenzBrayClarkViscosityTestFixture() = default;

  void testViscosity( ViscosityData< NC > const & data )
  {
    auto const mixing_type = static_cast< LohrenzBrayClarkViscosityUpdate::MixingType >(std::get< 0 >( data ));
    real64 const pressure = std::get< 1 >( data );
    real64 const temperature = std::get< 2 >( data );
    stackArray1d< real64, numComps > phaseComposition;
    TestFluid< NC >::createArray( phaseComposition, std::get< 3 >( data ));
    real64 const expectedViscosity = std::get< 4 >( data );

    real64 molarDensity = 0.0;
    real64 massDensity = 0.0;
    real64 viscosity = 0.0;
    stackArray1d< real64, numDofs > tempDerivs( numDofs );

    auto componentProperties = m_fluid->createKernelWrapper();
    auto densityKernelWrapper = m_density->createKernelWrapper();
    auto viscosityKernelWrapper = m_viscosity->createKernelWrapper();

    viscosityKernelWrapper.setMixingType( mixing_type );

    densityKernelWrapper.compute( componentProperties,
                                  pressure,
                                  temperature,
                                  phaseComposition.toSliceConst(),
                                  molarDensity,
                                  tempDerivs.toSlice(),
                                  massDensity,
                                  tempDerivs.toSlice(),
                                  false );

    viscosityKernelWrapper.compute( componentProperties,
                                    pressure,
                                    temperature,
                                    phaseComposition.toSliceConst(),
                                    massDensity,
                                    tempDerivs.toSliceConst(),
                                    viscosity,
                                    tempDerivs.toSlice(),
                                    false );

    checkRelativeError( viscosity, expectedViscosity, relTol, absTol );
  }

  void testViscosityDerivatives( ViscosityData< NC > const & data )
  {
    auto const mixing_type = static_cast< LohrenzBrayClarkViscosityUpdate::MixingType >(std::get< 0 >( data ));
    real64 const pressure = std::get< 1 >( data );
    real64 const temperature = std::get< 2 >( data );
    stackArray1d< real64, numComps > phaseComposition;
    TestFluid< NC >::createArray( phaseComposition, std::get< 3 >( data ));

    auto componentProperties = m_fluid->createKernelWrapper();
    auto densityKernelWrapper = m_density->createKernelWrapper();
    auto viscosityKernelWrapper = m_viscosity->createKernelWrapper();

    viscosityKernelWrapper.setMixingType( mixing_type );

    real64 molarDensity = 0.0;
    real64 massDensity = 0.0;
    real64 viscosity = 0.0;
    stackArray1d< real64, numDofs > molarDensityDerivs( numDofs );
    stackArray1d< real64, numDofs > massDensityDerivs( numDofs );
    stackArray1d< real64, numDofs > viscosityDerivs( numDofs );

    densityKernelWrapper.compute( componentProperties,
                                  pressure,
                                  temperature,
                                  phaseComposition.toSliceConst(),
                                  molarDensity,
                                  molarDensityDerivs.toSlice(),
                                  massDensity,
                                  massDensityDerivs.toSlice(),
                                  false );

    viscosityKernelWrapper.compute( componentProperties,
                                    pressure,
                                    temperature,
                                    phaseComposition.toSliceConst(),
                                    massDensity,
                                    massDensityDerivs.toSliceConst(),
                                    viscosity,
                                    viscosityDerivs.toSlice(),
                                    false );

    auto calculateViscosity = [&]( real64 const p, real64 const t, auto const & zmf ) -> real64 {
      real64 densityMolar = 0.0;
      real64 densityMass = 0.0;
      real64 phaseViscosity = 0.0;
      stackArray1d< real64, numDofs > tempDerivs( numDofs );
      densityKernelWrapper.compute( componentProperties, p, t, zmf.toSliceConst(),
                                    densityMolar, tempDerivs.toSlice(), densityMass, tempDerivs.toSlice(), false );
      viscosityKernelWrapper.compute( componentProperties, p, t, zmf.toSliceConst(),
                                      densityMass, tempDerivs.toSliceConst(), phaseViscosity, tempDerivs.toSlice(), false );
      return phaseViscosity;
    };

    // Viscosity values are very small so we will inflate the values to avoid false positives due
    // to the absolute value check
    real64 constexpr scale = 1.0e6;

    // Compare against numerical derivatives
    // -- Pressure derivative
    real64 const dp = 1.0e-4 * pressure;
    internal::testNumericalDerivative(
      pressure, dp, scale*viscosityDerivs[Deriv::dP],
      [&]( real64 const p ) -> real64 {
      return scale*calculateViscosity( p, temperature, phaseComposition );
    } );

    // -- Temperature derivative
    real64 const dT = 1.0e-6 * temperature;
    internal::testNumericalDerivative(
      temperature, dT, scale*viscosityDerivs[Deriv::dT],
      [&]( real64 const t ) -> real64 {
      return scale*calculateViscosity( pressure, t, phaseComposition );
    } );

    // -- Composition derivatives derivative
    real64 const dz = 1.0e-7;
    for( integer ic = 0; ic < 1; ++ic )
    {
      internal::testNumericalDerivative(
        0.0, dz, scale*viscosityDerivs[Deriv::dC + ic],
        [&]( real64 const z ) -> real64 {
        stackArray1d< real64, numComps > zmf( numComps );
        for( integer jc = 0; jc < numComps; ++jc )
        {
          zmf[jc] = phaseComposition[jc];
        }
        zmf[ic] += z;
        return scale*calculateViscosity( pressure, temperature, zmf );
      } );
    }
  }

protected:
  std::unique_ptr< TestFluid< NC > > m_fluid{};
  std::unique_ptr< CompositionalDensity > m_density{};
  std::unique_ptr< LohrenzBrayClarkViscosity > m_viscosity{};
  std::unique_ptr< ModelParameters > m_parameters{};
};

using LohrenzBrayClarkViscosity9 = LohrenzBrayClarkViscosityTestFixture< 9 >;

TEST_P( LohrenzBrayClarkViscosity9, testViscosity )
{
  testViscosity( GetParam() );
}

TEST_P( LohrenzBrayClarkViscosity9, testViscosityDerivatives )
{
  testViscosityDerivatives( GetParam() );
}

//-------------------------------------------------------------------------------
// Data
//-------------------------------------------------------------------------------

/* UNCRUSTIFY-OFF */

INSTANTIATE_TEST_SUITE_P(LohrenzBrayClarkViscosity, LohrenzBrayClarkViscosity9,
  ::testing::ValuesIn<ViscosityData< 9 >>( {
    {0, 1.839590e+06, 297.15, {0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920}, 1.277311e-05},
    {1, 1.839590e+06, 297.15, {0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920}, 1.110154e-05},
    {2, 1.839590e+06, 297.15, {0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920}, 1.570131e-05},
    {0, 1.839590e+06, 297.15, {0.00826, 0.00544, 0.77032, 0.10456, 0.06177, 0.02459, 0.00884, 0.00472, 0.01149}, 1.598556e-05},
    {1, 1.839590e+06, 297.15, {0.00826, 0.00544, 0.77032, 0.10456, 0.06177, 0.02459, 0.00884, 0.00472, 0.01149}, 1.525278e-05},
    {2, 1.839590e+06, 297.15, {0.00826, 0.00544, 0.77032, 0.10456, 0.06177, 0.02459, 0.00884, 0.00472, 0.01149}, 1.573849e-05},
    {0, 1.839590e+06, 297.15, {0.00899, 0.00299, 0.53281, 0.11447, 0.08791, 0.04566, 0.02095, 0.01516, 0.17107}, 1.275691e-05},
    {1, 1.839590e+06, 297.15, {0.00899, 0.00299, 0.53281, 0.11447, 0.08791, 0.04566, 0.02095, 0.01516, 0.17107}, 1.108572e-05},
    {2, 1.839590e+06, 297.15, {0.00899, 0.00299, 0.53281, 0.11447, 0.08791, 0.04566, 0.02095, 0.01516, 0.17107}, 1.570641e-05},
    {0, 1.839590e+06, 363.00, {0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920}, 1.441390e-05},
    {1, 1.839590e+06, 363.00, {0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920}, 1.261317e-05},
    {2, 1.839590e+06, 363.00, {0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920}, 1.790554e-05},
    {0, 1.839590e+06, 363.00, {0.00826, 0.00544, 0.77032, 0.10456, 0.06177, 0.02459, 0.00884, 0.00472, 0.01149}, 1.843194e-05},
    {1, 1.839590e+06, 363.00, {0.00826, 0.00544, 0.77032, 0.10456, 0.06177, 0.02459, 0.00884, 0.00472, 0.01149}, 1.766697e-05},
    {2, 1.839590e+06, 363.00, {0.00826, 0.00544, 0.77032, 0.10456, 0.06177, 0.02459, 0.00884, 0.00472, 0.01149}, 1.818931e-05},
    {0, 1.839590e+06, 363.00, {0.00899, 0.00299, 0.53281, 0.11447, 0.08791, 0.04566, 0.02095, 0.01516, 0.17107}, 1.438427e-05},
    {1, 1.839590e+06, 363.00, {0.00899, 0.00299, 0.53281, 0.11447, 0.08791, 0.04566, 0.02095, 0.01516, 0.17107}, 1.258368e-05},
    {2, 1.839590e+06, 363.00, {0.00899, 0.00299, 0.53281, 0.11447, 0.08791, 0.04566, 0.02095, 0.01516, 0.17107}, 1.790103e-05},
    {0, 1.839590e+08, 297.15, {0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920}, 4.180223e-04},
    {1, 1.839590e+08, 297.15, {0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920}, 4.163507e-04},
    {2, 1.839590e+08, 297.15, {0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920}, 4.209505e-04},
    {0, 1.839590e+08, 297.15, {0.00826, 0.00544, 0.77032, 0.10456, 0.06177, 0.02459, 0.00884, 0.00472, 0.01149}, 1.823560e-04},
    {1, 1.839590e+08, 297.15, {0.00826, 0.00544, 0.77032, 0.10456, 0.06177, 0.02459, 0.00884, 0.00472, 0.01149}, 1.816232e-04},
    {2, 1.839590e+08, 297.15, {0.00826, 0.00544, 0.77032, 0.10456, 0.06177, 0.02459, 0.00884, 0.00472, 0.01149}, 1.821089e-04},
    {0, 1.839590e+08, 297.15, {0.00899, 0.00299, 0.53281, 0.11447, 0.08791, 0.04566, 0.02095, 0.01516, 0.17107}, 4.200580e-04},
    {1, 1.839590e+08, 297.15, {0.00899, 0.00299, 0.53281, 0.11447, 0.08791, 0.04566, 0.02095, 0.01516, 0.17107}, 4.183868e-04},
    {2, 1.839590e+08, 297.15, {0.00899, 0.00299, 0.53281, 0.11447, 0.08791, 0.04566, 0.02095, 0.01516, 0.17107}, 4.230075e-04},
    {0, 1.839590e+08, 363.00, {0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920}, 3.016154e-04},
    {1, 1.839590e+08, 363.00, {0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920}, 2.998146e-04},
    {2, 1.839590e+08, 363.00, {0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920}, 3.051070e-04},
    {0, 1.839590e+08, 363.00, {0.00826, 0.00544, 0.77032, 0.10456, 0.06177, 0.02459, 0.00884, 0.00472, 0.01149}, 1.225458e-04},
    {1, 1.839590e+08, 363.00, {0.00826, 0.00544, 0.77032, 0.10456, 0.06177, 0.02459, 0.00884, 0.00472, 0.01149}, 1.217808e-04},
    {2, 1.839590e+08, 363.00, {0.00826, 0.00544, 0.77032, 0.10456, 0.06177, 0.02459, 0.00884, 0.00472, 0.01149}, 1.223032e-04},
    {0, 1.839590e+08, 363.00, {0.00899, 0.00299, 0.53281, 0.11447, 0.08791, 0.04566, 0.02095, 0.01516, 0.17107}, 3.034493e-04},
    {1, 1.839590e+08, 363.00, {0.00899, 0.00299, 0.53281, 0.11447, 0.08791, 0.04566, 0.02095, 0.01516, 0.17107}, 3.016487e-04},
    {2, 1.839590e+08, 363.00, {0.00899, 0.00299, 0.53281, 0.11447, 0.08791, 0.04566, 0.02095, 0.01516, 0.17107}, 3.069661e-04}
  } )
);

/* UNCRUSTIFY-ON */

} // testing
} // geos
