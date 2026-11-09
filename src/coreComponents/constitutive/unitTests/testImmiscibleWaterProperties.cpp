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
#include "constitutive/fluid/multifluid/compositional/parameters/ImmiscibleWaterParameters.hpp"
#include "constitutive/fluid/multifluid/compositional/models/ImmiscibleWaterDensity.hpp"
#include "constitutive/fluid/multifluid/compositional/models/ImmiscibleWaterViscosity.hpp"
#include "TestFluidFixture.hpp"
#include "TestFluidUtilities.hpp"

using namespace geos::constitutive::compositional;

namespace geos
{
namespace testing
{

template< integer NC >
struct TestData : public TestDataBase< NC >
{
  real64 expectedMolarDensity; ///< Expected molar density
  real64 expectedMassDensity;  ///< Expected mass density
  real64 expectedViscosity;    ///< Expected viscosity
};

template< TestFluidName FLUID >
class ImmiscibleWaterPropertiesTestFixture : public TestFluidFixture< FLUID, TestData >
{
  using Base = TestFluidFixture< FLUID, TestData >;
  using Base::numComps;
  using Base::numDofs;
  using Base::relTol;
  using Base::absTol;
  using typename Base::Deriv;
public:
  ImmiscibleWaterPropertiesTestFixture()
  {
    ComponentProperties const & componentProperties = this->m_fluid->getComponentProperties();

    m_parameters = ImmiscibleWaterDensity::createParameters( std::make_unique< ModelParameters >() );
    m_parameters = ImmiscibleWaterViscosity::createParameters( std::move( m_parameters ) );

    auto * waterParameters = const_cast< ImmiscibleWaterParameters * >(m_parameters->get< ImmiscibleWaterParameters >());
    waterParameters->m_waterReferencePressure = 215.0e5;
    waterParameters->m_waterDensity = 1020.0;
    waterParameters->m_waterCompressibility = 4.1483E-10;
    waterParameters->m_waterViscosityCompressibility = 2.0E-11;
    waterParameters->m_waterViscosity = 0.32929e-3;

    m_density = std::make_unique< ImmiscibleWaterDensity >( "PhaseDensity", componentProperties, 0, *m_parameters );
    m_viscosity = std::make_unique< ImmiscibleWaterViscosity >( "PhaseViscosity", componentProperties, 0, *m_parameters );
  }

  ~ImmiscibleWaterPropertiesTestFixture() = default;

  void testProperties( TestData< numComps > const & data )
  {
    real64 const pressure = data.pressure;
    real64 const temperature = data.temperature;
    stackArray1d< real64, numComps > phaseComposition;
    TestFluid< numComps >::createArray( phaseComposition, data.composition );
    real64 const expectedMolarDensity = data.expectedMolarDensity;
    real64 const expectedMassDensity = data.expectedMassDensity;
    real64 const expectedViscosity = data.expectedViscosity;

    real64 molarDensity = 0.0;
    real64 massDensity = 0.0;
    real64 viscosity = 0.0;
    stackArray1d< real64, numDofs > tempDerivs( numDofs );

    auto componentProperties = this->m_fluid->createKernelWrapper();
    auto densityKernelWrapper = m_density->createKernelWrapper();
    auto viscosityKernelWrapper = m_viscosity->createKernelWrapper();

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

    checkRelativeError( molarDensity, expectedMolarDensity, relTol, absTol );
    checkRelativeError( massDensity, expectedMassDensity, relTol, absTol );
    checkRelativeError( viscosity, expectedViscosity, relTol, absTol );
  }

  void testPropertyDerivatives( TestData< numComps > const & data )
  {
    stackArray1d< real64, numComps > phaseComposition;
    TestFluid< numComps >::createArray( phaseComposition, data.composition );

    auto componentProperties = this->m_fluid->createKernelWrapper();
    auto densityKernelWrapper = m_density->createKernelWrapper();
    auto viscosityKernelWrapper = m_viscosity->createKernelWrapper();

    // Values: molar density, mass density and viscosity
    auto computeProperties = [&]( real64 const p, real64 const t, arraySlice1d< real64 const > const & z,
                                  auto & values, auto & derivatives )
    {
      stackArray1d< real64, numDofs > massDensityDerivs( numDofs );
      densityKernelWrapper.compute( componentProperties, p, t, z,
                                    values[0], derivatives[0], values[1], massDensityDerivs.toSlice(), false );
      viscosityKernelWrapper.compute( componentProperties, p, t, z,
                                      values[1], massDensityDerivs.toSliceConst(), values[2], derivatives[2], false );
      for( integer idof = 0; idof < numDofs; ++idof )
      {
        derivatives( 1, idof ) = massDensityDerivs[idof];
      }
    };

    stackArray1d< real64, 3 > values( 3 );
    stackArray2d< real64, 3*numDofs > derivatives( 3, numDofs );
    computeProperties( data.pressure, data.temperature, phaseComposition.toSliceConst(), values, derivatives );

    DerivativeCheckSettings settings;
    settings.absTolerance = absTol;
    settings.relTolerance = relTol;

    // Viscosity values are very small so they are inflated to avoid false positives due
    // to the absolute value check
    stackArray2d< real64, 3*numDofs > tempDerivs( 3, numDofs );
    Base::template checkNumericalDerivatives< 3 >( data, derivatives.toSliceConst(),
                                                   [&]( real64 const p, real64 const t, auto const & z, auto & displacedValues )
    {
      computeProperties( p, t, z, displacedValues, tempDerivs );
    }, settings, { 1.0, 1.0, 1.0e6 } );
  }

protected:
  std::unique_ptr< ImmiscibleWaterDensity > m_density{};
  std::unique_ptr< ImmiscibleWaterViscosity > m_viscosity{};
  std::unique_ptr< ModelParameters > m_parameters{};
};

using ImmiscibleWaterProperties3 = ImmiscibleWaterPropertiesTestFixture< TestFluidName::GasOilWater >;

TEST_P( ImmiscibleWaterProperties3, testProperties )
{
  testProperties( GetParam() );
}

TEST_P( ImmiscibleWaterProperties3, testPropertyDerivatives )
{
  testPropertyDerivatives( GetParam() );
}

//-------------------------------------------------------------------------------
// Data
//-------------------------------------------------------------------------------

/* UNCRUSTIFY-OFF */

INSTANTIATE_TEST_SUITE_P(
  ImmiscibleWaterProperties, ImmiscibleWaterProperties3,
  ::testing::ValuesIn<TestData< 3 >>({
    {{2.0e+06, 293.15, {0.30, 0.30, 0.40}}, 5.616239e+04, 1.011782e+03, 3.291616e-04},
    {{1.5e+07, 293.15, {0.30, 0.30, 0.40}}, 5.646608e+04, 1.017253e+03, 3.292472e-04},
    {{6.0e+07, 293.15, {0.30, 0.30, 0.40}}, 5.753005e+04, 1.036421e+03, 3.295437e-04},
    {{2.0e+06, 353.15, {0.30, 0.30, 0.40}}, 5.616239e+04, 1.011782e+03, 3.291616e-04},
    {{1.5e+07, 353.15, {0.30, 0.30, 0.40}}, 5.646608e+04, 1.017253e+03, 3.292472e-04},
    {{6.0e+07, 353.15, {0.30, 0.30, 0.40}}, 5.753005e+04, 1.036421e+03, 3.295437e-04},
    {{2.0e+06, 573.15, {0.30, 0.30, 0.40}}, 5.616239e+04, 1.011782e+03, 3.291616e-04},
    {{1.5e+07, 573.15, {0.30, 0.30, 0.40}}, 5.646608e+04, 1.017253e+03, 3.292472e-04},
    {{6.0e+07, 573.15, {0.30, 0.30, 0.40}}, 5.753005e+04, 1.036421e+03, 3.295437e-04},
    {{2.0e+06, 293.15, {0.00, 0.00, 1.00}}, 5.616239e+04, 1.011782e+03, 3.291616e-04},
    {{1.5e+07, 293.15, {0.00, 0.00, 1.00}}, 5.646608e+04, 1.017253e+03, 3.292472e-04},
    {{6.0e+07, 293.15, {0.00, 0.00, 1.00}}, 5.753005e+04, 1.036421e+03, 3.295437e-04},
    {{2.0e+06, 353.15, {0.00, 0.00, 1.00}}, 5.616239e+04, 1.011782e+03, 3.291616e-04},
    {{1.5e+07, 353.15, {0.00, 0.00, 1.00}}, 5.646608e+04, 1.017253e+03, 3.292472e-04},
    {{6.0e+07, 353.15, {0.00, 0.00, 1.00}}, 5.753005e+04, 1.036421e+03, 3.295437e-04},
    {{2.0e+06, 573.15, {0.00, 0.00, 1.00}}, 5.616239e+04, 1.011782e+03, 3.291616e-04},
    {{1.5e+07, 573.15, {0.00, 0.00, 1.00}}, 5.646608e+04, 1.017253e+03, 3.292472e-04},
    {{6.0e+07, 573.15, {0.00, 0.00, 1.00}}, 5.753005e+04, 1.036421e+03, 3.295437e-04},
    {{2.0e+06, 293.15, {0.20, 0.80, 0.00}}, 5.616239e+04, 1.011782e+03, 3.291616e-04},
    {{1.5e+07, 293.15, {0.20, 0.80, 0.00}}, 5.646608e+04, 1.017253e+03, 3.292472e-04},
    {{6.0e+07, 293.15, {0.20, 0.80, 0.00}}, 5.753005e+04, 1.036421e+03, 3.295437e-04},
    {{2.0e+06, 353.15, {0.20, 0.80, 0.00}}, 5.616239e+04, 1.011782e+03, 3.291616e-04},
    {{1.5e+07, 353.15, {0.20, 0.80, 0.00}}, 5.646608e+04, 1.017253e+03, 3.292472e-04},
    {{6.0e+07, 353.15, {0.20, 0.80, 0.00}}, 5.753005e+04, 1.036421e+03, 3.295437e-04},
    {{2.0e+06, 573.15, {0.20, 0.80, 0.00}}, 5.616239e+04, 1.011782e+03, 3.291616e-04},
    {{1.5e+07, 573.15, {0.20, 0.80, 0.00}}, 5.646608e+04, 1.017253e+03, 3.292472e-04},
    {{6.0e+07, 573.15, {0.20, 0.80, 0.00}}, 5.753005e+04, 1.036421e+03, 3.295437e-04}
  })
);

/* UNCRUSTIFY-ON */

} // testing

} // geos
