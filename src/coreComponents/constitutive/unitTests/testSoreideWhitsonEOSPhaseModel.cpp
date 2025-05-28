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
#include "constitutive/fluid/multifluid/compositional/functions/SoreideWhitsonEOSPhaseModel.hpp"
#include "TestFluid.hpp"
#include "TestFluidUtilities.hpp"

namespace geos
{

using namespace constitutive;
using namespace constitutive::compositional;

namespace testing
{

template< int NC >
struct FluidData {};

template<>
struct FluidData< 3 >
{
  static std::unique_ptr< TestFluid< 3 > > create()
  {
    return TestFluid< 3 >::create( {Fluid::H2O, Fluid::C1, Fluid::H2} );
  }

  static std::array< Feed< 3 >, 3 > constexpr feeds = {
    Feed< 3 >{0.995, 0.000, 0.005},
    Feed< 3 >{0.990, 0.005, 0.005},
    Feed< 3 >{0.970, 0.025, 0.005}
  };
};

template<>
struct FluidData< 4 >
{
  static std::unique_ptr< TestFluid< 4 > > create()
  {
    return TestFluid< 4 >::create( {Fluid::N2, Fluid::C1, Fluid::CO2, Fluid::H2O} );
  }

  static std::array< Feed< 4 >, 3 > constexpr feeds = {
    Feed< 4 >{0.030933, 0.319683, 0.637861, 0.011523},
    Feed< 4 >{0.000000, 0.349686, 0.637891, 0.012423},
    Feed< 4 >{0.000000, 0.349686, 0.650314, 0.000000}
  };
};

template< int NC >
using TestData = std::tuple<
  real64 const,         // Pressure
  real64 const,         // Temperature
  real64 const,         // Salinity
  Feed< NC > const      // Input composition
  >;

template< integer NC, typename EOS_TYPE >
class SoreideWhitsonEOSPhaseModelTestFixture : public ::testing::TestWithParam< TestData< NC > >
{
public:
  static constexpr integer numComps = NC;
  static constexpr integer numDof = NC + 2;
  static constexpr real64 absTol = 1.0e-4;
  static constexpr real64 relTol = 1.0e-5;
  using ParamType = TestData< NC >;
  using EOS = SoreideWhitsonEOSPhaseModel< EOS_TYPE >;
  using Deriv = typename EOS::Deriv;
public:
  SoreideWhitsonEOSPhaseModelTestFixture();
  ~SoreideWhitsonEOSPhaseModelTestFixture() = default;

  void testPureCoefficients( ParamType const & testData );
  void testBinaryInteractionCoefficients( ParamType const & testData );
  void testMixtureCoefficients( ParamType const & testData );
  void testCompressibilityFactor( ParamType const & testData );
  void testLogFugacityCoefficients( ParamType const & testData );

protected:
  std::unique_ptr< TestFluid< NC > > m_fluid{};
};

template< integer NC, typename EOS_TYPE >
SoreideWhitsonEOSPhaseModelTestFixture< NC, EOS_TYPE >::SoreideWhitsonEOSPhaseModelTestFixture():
  m_fluid( FluidData< NC >::create() )
{}

template< int NC >
std::vector< TestData< NC > > generateTestData()
{
  auto const pressures = {1.0e+05, 1.83959e+06, 1.83959e+08};
  auto const temperatures = {297.15, 363.0};
  auto const salinities = {0.0, 1.7};
  std::vector< TestData< NC > > testData;
  for( const auto & composition : FluidData< NC >::feeds )
  {
    for( const real64 pressure : pressures )
    {
      for( const real64 temperature : temperatures )
      {
        for( const real64 salinity : salinities )
        {
          testData.emplace_back( pressure, temperature, salinity, composition );
        }
      }
    }
  }
  return testData;
}

template< integer NC, typename EOS_TYPE >
void
SoreideWhitsonEOSPhaseModelTestFixture< NC, EOS_TYPE >::testPureCoefficients( ParamType const & testData )
{
  auto componentProperties = this->m_fluid->createKernelWrapper();
  real64 const pressure = std::get< 0 >( testData );
  real64 const temperature = std::get< 1 >( testData );
  real64 const salinity = std::get< 2 >( testData );

  real64 aCoefficient = 0.0;
  real64 bCoefficient = 0.0;
  real64 daCoefficient_dp = 0.0;
  real64 dbCoefficient_dp = 0.0;
  real64 daCoefficient_dt = 0.0;
  real64 dbCoefficient_dt = 0.0;

  for( integer ic = 0; ic < numComps; ++ic )
  {
    EOS::computePureCoefficientsAndDerivs( ic,
                                           pressure,
                                           temperature,
                                           componentProperties,
                                           salinity,
                                           aCoefficient,
                                           bCoefficient,
                                           daCoefficient_dp,
                                           dbCoefficient_dp,
                                           daCoefficient_dt,
                                           dbCoefficient_dt );

    real64 const dp = 1.0e-4 * pressure;
    internal::testNumericalDerivative( pressure, dp, daCoefficient_dp,
                                       [&]( real64 p ) -> real64 {
      real64 a = 0.0, b = 0.0;
      EOS::computePureCoefficients( ic, p, temperature, componentProperties, salinity, a, b );
      return a;
    }, absTol, relTol );
    internal::testNumericalDerivative( pressure, dp, dbCoefficient_dp,
                                       [&]( real64 p ) -> real64 {
      real64 a = 0.0, b = 0.0;
      EOS::computePureCoefficients( ic, p, temperature, componentProperties, salinity, a, b );
      return b;
    }, absTol, relTol );

    real64 const dT = 1.0e-6 * temperature;
    internal::testNumericalDerivative( temperature, dT, daCoefficient_dt,
                                       [&]( real64 t ) -> real64 {
      real64 a = 0.0, b = 0.0;
      EOS::computePureCoefficients( ic, pressure, t, componentProperties, salinity, a, b );
      return a;
    }, absTol, relTol );
    internal::testNumericalDerivative( temperature, dT, dbCoefficient_dt,
                                       [&]( real64 t ) -> real64 {
      real64 a = 0.0, b = 0.0;
      EOS::computePureCoefficients( ic, pressure, t, componentProperties, salinity, a, b );
      return b;
    }, absTol, relTol );
  }
}

template< integer NC, typename EOS_TYPE >
void
SoreideWhitsonEOSPhaseModelTestFixture< NC, EOS_TYPE >::testBinaryInteractionCoefficients( ParamType const & testData )
{
  auto componentProperties = this->m_fluid->createKernelWrapper();
  real64 const pressure = std::get< 0 >( testData );
  real64 const temperature = std::get< 1 >( testData );
  real64 const salinity = std::get< 2 >( testData );

  // Test symmetry
  for( integer ic = 0; ic < numComps; ++ic )
  {
    for( integer jc = ic; jc < numComps; ++jc )
    {
      real64 kij = 0.0;
      real64 kji = 0.0;
      real64 dk_dT = 0.0;
      EOS::getBinaryInteractionCoefficient( pressure,
                                            temperature,
                                            componentProperties,
                                            salinity,
                                            ic,
                                            jc,
                                            kij,
                                            dk_dT );
      EOS::getBinaryInteractionCoefficient( pressure,
                                            temperature,
                                            componentProperties,
                                            salinity,
                                            ic,
                                            jc,
                                            kji,
                                            dk_dT );
      EXPECT_NEAR( kij, kji, absTol );
      if( ic == jc )
      {
        EXPECT_NEAR( kij, 0.0, absTol );
      }
    }
  }

  // Test numerical derivatives
  for( integer ic = 0; ic < numComps; ++ic )
  {
    for( integer jc = 0; jc < numComps; ++jc )
    {
      real64 kij = 0.0;
      real64 dkij_dT = 0.0;
      EOS::getBinaryInteractionCoefficient( pressure,
                                            temperature,
                                            componentProperties,
                                            salinity,
                                            ic,
                                            jc,
                                            kij,
                                            dkij_dT );

      real64 const dT = 1.0e-6 * temperature;
      internal::testNumericalDerivative( temperature, dT, dkij_dT,
                                         [&]( real64 t ) -> real64 {
        real64 l_kij, l_dkij_dT;
        EOS::getBinaryInteractionCoefficient( pressure,
                                              t,
                                              componentProperties,
                                              salinity,
                                              ic,
                                              jc,
                                              l_kij,
                                              l_dkij_dT );
        return l_kij;
      }, absTol, relTol );
    }
  }
}

template< integer NC, typename EOS_TYPE >
void
SoreideWhitsonEOSPhaseModelTestFixture< NC, EOS_TYPE >::testMixtureCoefficients( ParamType const & testData )
{
  auto componentProperties = this->m_fluid->createKernelWrapper();
  real64 const pressure = std::get< 0 >( testData );
  real64 const temperature = std::get< 1 >( testData );
  real64 const salinity = std::get< 2 >( testData );
  stackArray1d< real64, numComps > composition;
  TestFluid< NC >::createArray( composition, std::get< 3 >( testData ));

  integer constexpr numValues = 2;

  stackArray1d< real64, numValues > mixtureCoefficient( numValues );
  stackArray2d< real64, numValues *numComps > pureCoefficients( numValues, numComps );
  stackArray2d< real64, numValues *numDof > mixtureCoefficientDerivs( numValues, numDof );

  auto const & aCoefficients = pureCoefficients[0];
  auto const & bCoefficients = pureCoefficients[1];
  auto const & aCoefficientDerivs = mixtureCoefficientDerivs[0];
  auto const & bCoefficientDerivs = mixtureCoefficientDerivs[1];

  EOS::computeMixtureCoefficients( numComps,
                                   pressure,
                                   temperature,
                                   composition.toSliceConst(),
                                   componentProperties,
                                   salinity,
                                   aCoefficients,
                                   bCoefficients,
                                   mixtureCoefficient[0],
                                   mixtureCoefficient[1] );

  EOS::computeMixtureCoefficientDerivs( numComps,
                                        pressure,
                                        temperature,
                                        composition.toSliceConst(),
                                        componentProperties,
                                        salinity,
                                        aCoefficients.toSliceConst(),
                                        bCoefficients.toSliceConst(),
                                        mixtureCoefficient[0],
                                        mixtureCoefficient[1],
                                        aCoefficientDerivs,
                                        bCoefficientDerivs );

  stackArray1d< real64, numValues > derivatives( numValues );

  // Pressure derivatives
  real64 const dp = 1.0e-4 * pressure;
  derivatives[0] = aCoefficientDerivs[Deriv::dP];
  derivatives[1] = bCoefficientDerivs[Deriv::dP];
  internal::testNumericalDerivative< numValues >( pressure, dp, derivatives.toSliceConst(),
                                                  [&]( real64 const p, auto & values )
  {
    SoreideWhitsonEOSPhaseModel< EOS_TYPE >::
    computeMixtureCoefficients( numComps,
                                p,
                                temperature,
                                composition.toSliceConst(),
                                componentProperties,
                                salinity,
                                aCoefficients,
                                bCoefficients,
                                values[0],
                                values[1] );
  }, absTol, relTol );

  // Temperature derivatives
  real64 const dT = 1.0e-6 * temperature;
  derivatives[0] = aCoefficientDerivs[Deriv::dT];
  derivatives[1] = bCoefficientDerivs[Deriv::dT];
  internal::testNumericalDerivative< numValues >( temperature, dT, derivatives.toSliceConst(),
                                                  [&]( real64 const t, auto & values )
  {
    SoreideWhitsonEOSPhaseModel< EOS_TYPE >::
    computeMixtureCoefficients( numComps,
                                pressure,
                                t,
                                composition.toSliceConst(),
                                componentProperties,
                                salinity,
                                aCoefficients,
                                bCoefficients,
                                values[0],
                                values[1] );
  }, absTol, relTol );

  // Composition derivatives
  real64 const dz = 1.0e-7;
  for( integer kc = 0; kc < numComps; kc++ )
  {
    derivatives[0] = aCoefficientDerivs[Deriv::dC+kc];
    derivatives[1] = bCoefficientDerivs[Deriv::dC+kc];
    internal::testNumericalDerivative< numValues >( 0, dz, derivatives.toSliceConst(),
                                                    [&]( real64 const z, auto & values )
    {
      real64 const z_old = composition[kc];
      composition[kc] += z;
      SoreideWhitsonEOSPhaseModel< EOS_TYPE >::
      computeMixtureCoefficients( numComps,
                                  pressure,
                                  temperature,
                                  composition.toSliceConst(),
                                  componentProperties,
                                  salinity,
                                  aCoefficients,
                                  bCoefficients,
                                  values[0],
                                  values[1] );
      composition[kc] = z_old;
    }, absTol, relTol );
  }
}

template< integer NC, typename EOS_TYPE >
void
SoreideWhitsonEOSPhaseModelTestFixture< NC, EOS_TYPE >::testCompressibilityFactor( ParamType const & testData )
{
  auto componentProperties = this->m_fluid->createKernelWrapper();
  real64 const pressure = std::get< 0 >( testData );
  real64 const temperature = std::get< 1 >( testData );
  real64 const salinity = std::get< 2 >( testData );
  stackArray1d< real64, numComps > composition;
  TestFluid< NC >::createArray( composition, std::get< 3 >( testData ));

  real64 compressibilityFactor = 0.0;
  stackArray1d< real64, numDof > compressibilityFactorDerivs( numDof );

  EOS::computeCompressibilityFactorAndDerivs( numComps,
                                              pressure,
                                              temperature,
                                              composition.toSliceConst(),
                                              componentProperties,
                                              salinity,
                                              compressibilityFactor,
                                              compressibilityFactorDerivs.toSlice() );

  auto computeCompressibility = [&]( real64 p, real64 t, auto z ) -> real64 {
    real64 zfactor = 0.0;
    EOS::computeCompressibilityFactor( numComps, p, t, z, componentProperties, salinity, zfactor );
    return zfactor;
  };

  real64 const dp = 1.0e-4 * pressure;
  internal::testNumericalDerivative( pressure, dp, compressibilityFactorDerivs[Deriv::dP],
                                     [&]( real64 p )
  {
    return computeCompressibility( p, temperature, composition.toSliceConst() );
  }, absTol, relTol );

  real64 const dT = 1.0e-6 * temperature;
  internal::testNumericalDerivative( temperature, dT, compressibilityFactorDerivs[Deriv::dT],
                                     [&]( real64 t )
  {
    return computeCompressibility( pressure, t, composition.toSliceConst() );
  }, absTol, relTol );

  real64 const dz = 1.0e-7;
  for( integer kc = 0; kc < numComps; kc++ )
  {
    internal::testNumericalDerivative( 0.0, dz, compressibilityFactorDerivs[Deriv::dC+kc],
                                       [&]( real64 z )
    {
      real64 const z_old = composition[kc];
      composition[kc] += z;
      real64 zfactor = computeCompressibility( pressure, temperature, composition.toSliceConst() );
      composition[kc] = z_old;
      return zfactor;
    }, absTol, relTol );
  }
}


template< integer NC, typename EOS_TYPE >
void
SoreideWhitsonEOSPhaseModelTestFixture< NC, EOS_TYPE >::testLogFugacityCoefficients( ParamType const & testData )
{
  auto componentProperties = this->m_fluid->createKernelWrapper();
  real64 const pressure = std::get< 0 >( testData );
  real64 const temperature = std::get< 1 >( testData );
  real64 const salinity = std::get< 2 >( testData );
  stackArray1d< real64, numComps > composition;
  TestFluid< NC >::createArray( composition, std::get< 3 >( testData ));

  stackArray1d< real64, numComps > logFugacityCoefficients( numComps );
  stackArray2d< real64, numComps *numDof > logFugacityCoefficientDerivs( numComps, numDof );

  // Inflate the values of the log fugacity coefficients to catch errors
  real64 const scale = 1.0e3;

  EOS::computeLogFugacityCoefficients( numComps,
                                       pressure,
                                       temperature,
                                       composition.toSliceConst(),
                                       componentProperties,
                                       salinity,
                                       logFugacityCoefficients.toSlice() );

  EOS::computeLogFugacityCoefficientDerivs( numComps,
                                            pressure,
                                            temperature,
                                            composition.toSliceConst(),
                                            componentProperties,
                                            salinity,
                                            logFugacityCoefficients.toSliceConst(),
                                            logFugacityCoefficientDerivs.toSlice() );

  stackArray1d< real64, numComps > derivatives( numComps );

  auto const concatDerivatives = [&]( integer idof )
  {
    for( integer ic = 0; ic < numComps; ++ic )
    {
      derivatives[ic] = scale * logFugacityCoefficientDerivs( ic, idof );
    }
  };

  // Pressure derivatives
  concatDerivatives( Deriv::dP );
  real64 const dp = 1.0e-4 * pressure;
  internal::testNumericalDerivative< numComps >( pressure, dp, derivatives.toSliceConst(),
                                                 [&]( real64 const p, auto & values )
  {
    EOS::computeLogFugacityCoefficients( numComps,
                                         p,
                                         temperature,
                                         composition.toSliceConst(),
                                         componentProperties,
                                         salinity,
                                         values );
    LvArray::forValuesInSlice( values.toSlice(), [scale]( real64 & a ){ a *= scale; } );
  }, absTol, relTol );

  // Temperature derivatives
  concatDerivatives( Deriv::dT );
  real64 const dT = 1.0e-6 * temperature;
  internal::testNumericalDerivative< numComps >( temperature, dT, derivatives.toSliceConst(),
                                                 [&]( real64 const t, auto & values )
  {
    EOS::computeLogFugacityCoefficients( numComps,
                                         pressure,
                                         t,
                                         composition.toSliceConst(),
                                         componentProperties,
                                         salinity,
                                         values );
    LvArray::forValuesInSlice( values.toSlice(), [scale]( real64 & a ){ a *= scale; } );
  }, absTol, relTol );

  // Composition derivatives
  real64 constexpr dz = 1.0e-6;
  for( integer jc = 0; jc < numComps; ++jc )
  {
    concatDerivatives( Deriv::dC + jc );
    internal::testNumericalDerivative< numComps >( 0.0, dz, derivatives.toSliceConst(),
                                                   [&]( real64 const z, auto & values )
    {
      real64 const zj_old = composition[jc];
      composition[jc] += z;
      EOS::computeLogFugacityCoefficients( numComps,
                                           pressure,
                                           temperature,
                                           composition.toSliceConst(),
                                           componentProperties,
                                           salinity,
                                           values );
      composition[jc] = zj_old;
      LvArray::forValuesInSlice( values.toSlice(), [scale]( real64 & a ){ a *= scale; } );
    }, absTol, relTol );
  }
}

using PengRobinson4 = SoreideWhitsonEOSPhaseModelTestFixture< 4, PengRobinsonEOS >;
using SoaveRedlichKwong3 = SoreideWhitsonEOSPhaseModelTestFixture< 3, SoaveRedlichKwongEOS >;

TEST_P( PengRobinson4, testSWModel )
{
  auto const testParam = GetParam();
  testPureCoefficients( testParam );
  testBinaryInteractionCoefficients( testParam );
  testMixtureCoefficients( testParam );
  testCompressibilityFactor( testParam );
  testLogFugacityCoefficients( testParam );
}

TEST_P( SoaveRedlichKwong3, testSWModel )
{
  auto const testParam = GetParam();
  testPureCoefficients( testParam );
  testBinaryInteractionCoefficients( testParam );
  testMixtureCoefficients( testParam );
  testCompressibilityFactor( testParam );
  testLogFugacityCoefficients( testParam );
}

INSTANTIATE_TEST_SUITE_P( SoreideWhitsonEOSPhaseModelTest, PengRobinson4, ::testing::ValuesIn( generateTestData< 4 >()) );
INSTANTIATE_TEST_SUITE_P( SoreideWhitsonEOSPhaseModelTest, SoaveRedlichKwong3, ::testing::ValuesIn( generateTestData< 3 >()) );

} // namespace testing

} // namespace geos
