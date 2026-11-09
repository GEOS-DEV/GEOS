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
#include "TestFluidFixture.hpp"
#include "TestFluidUtilities.hpp"

using namespace geos::constitutive::compositional;

namespace geos
{
namespace testing
{

template< integer NC >
struct DensityData : public TestDataBase< NC >
{
  real64 expectedMolarDensity; ///< Expected molar density
  real64 expectedMassDensity;  ///< Expected mass density
};

template< TestFluidName FLUID, EquationOfStateType EOS_TYPE >
class CompositionalDensityTestFixture : public TestFluidFixture< FLUID, DensityData >
{
  using Base = TestFluidFixture< FLUID, DensityData >;
  using Base::numComps;
  using Base::numDofs;
  using Base::relTol;
  using Base::absTol;
  using typename Base::Deriv;
public:
  CompositionalDensityTestFixture()
  {
    ComponentProperties const & componentProperties = this->m_fluid->getComponentProperties();
    m_parameters = CompositionalDensity::createParameters( std::make_unique< ModelParameters >() );

    auto equationOfState = const_cast< EquationOfState * >(m_parameters->get< EquationOfState >());
    string const eosName = EnumStrings< EquationOfStateType >::toString( EOS_TYPE );
    equationOfState->m_equationsOfStateNames.emplace_back( eosName );

    m_density = std::make_unique< CompositionalDensity >( "PhaseDensity", componentProperties, 0, *m_parameters );
  }

  ~CompositionalDensityTestFixture() = default;

  void testDensityValues( DensityData< numComps > const & data )
  {
    real64 const pressure = data.pressure;
    real64 const temperature = data.temperature;
    stackArray1d< real64, numComps > phaseComposition;
    TestFluid< numComps >::createArray( phaseComposition, data.composition );
    real64 const expectedMolarDensity = data.expectedMolarDensity;
    real64 const expectedMassDensity = data.expectedMassDensity;

    auto componentProperties = this->m_fluid->createKernelWrapper();
    auto kernelWrapper = this->m_density->createKernelWrapper();

    real64 molarDensity = 0.0;
    real64 massDensity = 0.0;
    stackArray1d< real64, numDofs > tempDerivs( numDofs );

    kernelWrapper.compute( componentProperties,
                           pressure,
                           temperature,
                           phaseComposition.toSliceConst(),
                           molarDensity,
                           tempDerivs.toSlice(),
                           massDensity,
                           tempDerivs.toSlice(),
                           false );

    checkRelativeError( molarDensity, expectedMolarDensity, relTol, absTol );
    checkRelativeError( massDensity, expectedMassDensity, relTol, absTol );
  }

  void testDensityDerivatives( DensityData< numComps > const & data )
  {
    stackArray1d< real64, numComps > phaseComposition;
    TestFluid< numComps >::createArray( phaseComposition, data.composition );

    auto componentProperties = this->m_fluid->createKernelWrapper();
    auto kernelWrapper = m_density->createKernelWrapper();

    // Values: molar density and mass density
    real64 molarDensity = 0.0;
    real64 massDensity = 0.0;
    stackArray2d< real64, 2*numDofs > derivatives( 2, numDofs );
    kernelWrapper.compute( componentProperties,
                           data.pressure,
                           data.temperature,
                           phaseComposition.toSliceConst(),
                           molarDensity,
                           derivatives[0],
                           massDensity,
                           derivatives[1],
                           false );

    stackArray1d< real64, numDofs > tempDerivs( numDofs );
    Base::template checkNumericalDerivatives< 2 >( data, derivatives.toSliceConst(),
                                                   [&]( real64 const p, real64 const t, auto const & z, auto & values )
    {
      kernelWrapper.compute( componentProperties, p, t, z,
                             values[0], tempDerivs.toSlice(), values[1], tempDerivs.toSlice(), false );
    } );
  }

protected:
  std::unique_ptr< CompositionalDensity > m_density{};
  std::unique_ptr< ModelParameters > m_parameters{};
};

using PengRobinson = CompositionalDensityTestFixture< TestFluidName::VolatileOil, EquationOfStateType::PengRobinson >;
using SoaveRedlichKwong = CompositionalDensityTestFixture< TestFluidName::VolatileOil, EquationOfStateType::SoaveRedlichKwong >;

TEST_P( PengRobinson, testDensityDerivatives )
{
  testDensityDerivatives( GetParam() );
}

TEST_P( SoaveRedlichKwong, testDensityDerivatives )
{
  testDensityDerivatives( GetParam() );
}

/* UNCRUSTIFY-OFF */

INSTANTIATE_TEST_SUITE_P(
  CompositionalDensityTest, PengRobinson,
  ::testing::ValuesIn( {
      DensityData< 9 >{{1.839590e+06, 2.971500e+02, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}}, 8.355571e+03, 4.559906e+02},
      DensityData< 9 >{{1.839590e+06, 2.971500e+02, {0.008260, 0.005440, 0.770320, 0.104560, 0.061770, 0.024590, 0.008840, 0.004720, 0.011490}}, 7.703898e+02, 2.691914e+01},
      DensityData< 9 >{{1.839590e+06, 2.971500e+02, {0.008990, 0.002990, 0.532810, 0.114470, 0.087910, 0.045660, 0.020950, 0.015160, 0.171070}}, 8.337694e+03, 4.567935e+02},
      DensityData< 9 >{{1.839590e+06, 3.630000e+02, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}}, 9.073321e+02, 4.951606e+01},
      DensityData< 9 >{{1.839590e+06, 3.630000e+02, {0.008260, 0.005440, 0.770320, 0.104560, 0.061770, 0.024590, 0.008840, 0.004720, 0.011490}}, 6.178234e+02, 2.158813e+01},
      DensityData< 9 >{{1.839590e+06, 3.630000e+02, {0.008990, 0.002990, 0.532810, 0.114470, 0.087910, 0.045660, 0.020950, 0.015160, 0.171070}}, 9.197865e+02, 5.039192e+01},
      DensityData< 9 >{{1.839590e+08, 2.971500e+02, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}}, 1.095078e+04, 5.976195e+02},
      DensityData< 9 >{{1.839590e+08, 2.971500e+02, {0.008260, 0.005440, 0.770320, 0.104560, 0.061770, 0.024590, 0.008840, 0.004720, 0.011490}}, 2.480270e+04, 8.666618e+02},
      DensityData< 9 >{{1.839590e+08, 2.971500e+02, {0.008990, 0.002990, 0.532810, 0.114470, 0.087910, 0.045660, 0.020950, 0.015160, 0.171070}}, 1.087917e+04, 5.960323e+02},
      DensityData< 9 >{{1.839590e+08, 3.630000e+02, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}}, 1.065128e+04, 5.812747e+02},
      DensityData< 9 >{{1.839590e+08, 3.630000e+02, {0.008260, 0.005440, 0.770320, 0.104560, 0.061770, 0.024590, 0.008840, 0.004720, 0.011490}}, 2.305823e+04, 8.057060e+02},
      DensityData< 9 >{{1.839590e+08, 3.630000e+02, {0.008990, 0.002990, 0.532810, 0.114470, 0.087910, 0.045660, 0.020950, 0.015160, 0.171070}}, 1.058381e+04, 5.798506e+02}
    } )
  );
INSTANTIATE_TEST_SUITE_P(
  CompositionalDensityTest, SoaveRedlichKwong,
  ::testing::ValuesIn( {
      DensityData< 9 >{{1.839590e+06, 2.971500e+02, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}}, 7.433979e+03, 4.056963e+02},
      DensityData< 9 >{{1.839590e+06, 2.971500e+02, {0.008260, 0.005440, 0.770320, 0.104560, 0.061770, 0.024590, 0.008840, 0.004720, 0.011490}}, 7.629968e+02, 2.666082e+01},
      DensityData< 9 >{{1.839590e+06, 2.971500e+02, {0.008990, 0.002990, 0.532810, 0.114470, 0.087910, 0.045660, 0.020950, 0.015160, 0.171070}}, 7.416959e+03, 4.063495e+02},
      DensityData< 9 >{{1.839590e+06, 3.630000e+02, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}}, 8.919848e+02, 4.867851e+01},
      DensityData< 9 >{{1.839590e+06, 3.630000e+02, {0.008260, 0.005440, 0.770320, 0.104560, 0.061770, 0.024590, 0.008840, 0.004720, 0.011490}}, 6.133569e+02, 2.143206e+01},
      DensityData< 9 >{{1.839590e+06, 3.630000e+02, {0.008990, 0.002990, 0.532810, 0.114470, 0.087910, 0.045660, 0.020950, 0.015160, 0.171070}}, 9.045641e+02, 4.955794e+01},
      DensityData< 9 >{{1.839590e+08, 2.971500e+02, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}}, 9.868675e+03, 5.385656e+02},
      DensityData< 9 >{{1.839590e+08, 2.971500e+02, {0.008260, 0.005440, 0.770320, 0.104560, 0.061770, 0.024590, 0.008840, 0.004720, 0.011490}}, 2.257420e+04, 7.887929e+02},
      DensityData< 9 >{{1.839590e+08, 2.971500e+02, {0.008990, 0.002990, 0.532810, 0.114470, 0.087910, 0.045660, 0.020950, 0.015160, 0.171070}}, 9.803814e+03, 5.371171e+02},
      DensityData< 9 >{{1.839590e+08, 3.630000e+02, {0.009000, 0.003000, 0.534700, 0.114600, 0.087900, 0.045600, 0.020900, 0.015100, 0.169200}}, 9.615674e+03, 5.247585e+02},
      DensityData< 9 >{{1.839590e+08, 3.630000e+02, {0.008260, 0.005440, 0.770320, 0.104560, 0.061770, 0.024590, 0.008840, 0.004720, 0.011490}}, 2.109248e+04, 7.370183e+02},
      DensityData< 9 >{{1.839590e+08, 3.630000e+02, {0.008990, 0.002990, 0.532810, 0.114470, 0.087910, 0.045660, 0.020950, 0.015160, 0.171070}}, 9.554300e+03, 5.234471e+02}
    } )
  );

/* UNCRUSTIFY-ON */

} // testing

} // geos
