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
#include "constitutive/fluid/multifluid/compositional/models/PhillipsBrineViscosity.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/BrineSalinity.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/PressureTemperatureCoordinates.hpp"
#include "TestFluidFixture.hpp"
#include "TestFluidUtilities.hpp"
#include "functions/FunctionManager.hpp"

#include <conduit.hpp>

using namespace geos::constitutive::compositional;

namespace geos
{
namespace testing
{

template< integer NC >
struct ViscosityData : public TestDataBase< NC >
{
  real64 expectedViscosity; ///< Expected viscosity
};

template< integer NC >
std::ostream & operator<<( std::ostream & os, ViscosityData< NC > const & data )
{
  return TestDataPrinter( os, data ).output( data.expectedViscosity );
}

template< TestFluidName FLUID, int SALINITY = 0 >
class PhillipsBrineViscosityTestFixture : public TestFluidFixture< FLUID, ViscosityData >
{
  using Base = TestFluidFixture< FLUID, ViscosityData >;
  using Base::numComps;
  using Base::numDofs;
  using Base::relTol;
  using Base::absTol;
  using typename Base::Deriv;

public:
  PhillipsBrineViscosityTestFixture()
    : m_parent( "parent", m_node )
  {
    m_functionManager = std::make_unique< FunctionManager >( "FunctionManager", &m_parent );

    ComponentProperties const & componentProperties = this->m_fluid->getComponentProperties();

    m_parameters = PhillipsBrineViscosity::createParameters( std::move( m_parameters ) );

    auto coordinates = const_cast< PressureTemperatureCoordinates * >(m_parameters->get< PressureTemperatureCoordinates >());
    createTemperaturePoints( coordinates->m_temperatureCoordinates );

    auto brineSalinity = const_cast< BrineSalinity * >(m_parameters->get< BrineSalinity >());
    real64 const massFraction = 1.0e-6*SALINITY;
    brineSalinity->m_salinity = massFraction / 58.44e-3;

    m_viscosity = std::make_unique< PhillipsBrineViscosity >( "PhaseViscosity", componentProperties, 0, *m_parameters );
  }

  ~PhillipsBrineViscosityTestFixture() = default;

  void testViscosity( ViscosityData< numComps > const & data )
  {
    real64 const pressure = data.pressure;
    real64 const temperature = data.temperature;
    stackArray1d< real64, numComps > phaseComposition;
    TestFluid< numComps >::createArray( phaseComposition, data.composition );
    real64 const expectedViscosity = data.expectedViscosity;

    real64 const massDensity = 1000.0;
    real64 viscosity = 0.0;
    stackArray1d< real64, numDofs > tempDerivs( numDofs );

    auto componentProperties = this->m_fluid->createKernelWrapper();
    auto viscosityKernelWrapper = m_viscosity->createKernelWrapper();

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

  void testViscosityDerivatives( ViscosityData< numComps > const & data )
  {
    stackArray1d< real64, numComps > phaseComposition;
    TestFluid< numComps >::createArray( phaseComposition, data.composition );

    auto componentProperties = this->m_fluid->createKernelWrapper();
    auto viscosityKernelWrapper = m_viscosity->createKernelWrapper();

    // The mass density is fixed, so its derivatives are zero
    real64 const massDensity = 1000.0;
    stackArray1d< real64, numDofs > massDensityDerivs( numDofs );
    LvArray::forValuesInSlice( massDensityDerivs.toSlice(), []( real64 & a ){ a = 0.0; } );

    real64 viscosity = 0.0;
    stackArray2d< real64, numDofs > derivatives( 1, numDofs );
    viscosityKernelWrapper.compute( componentProperties,
                                    data.pressure,
                                    data.temperature,
                                    phaseComposition.toSliceConst(),
                                    massDensity,
                                    massDensityDerivs.toSliceConst(),
                                    viscosity,
                                    derivatives[0],
                                    false );

    // Viscosity values are very small so they are inflated to avoid false positives due
    // to the absolute value check
    stackArray1d< real64, numDofs > tempDerivs( numDofs );
    Base::template checkNumericalDerivatives< 1 >( data, derivatives.toSliceConst(),
                                                   [&]( real64 const p, real64 const t, auto const & z, auto & values )
    {
      viscosityKernelWrapper.compute( componentProperties, p, t, z,
                                      massDensity, massDensityDerivs.toSliceConst(), values[0], tempDerivs.toSlice(), false );
    }, {}, { 1.0e6 } );
  }

private:
  static void createTemperaturePoints( array1d< real64 > & temperatureCoordinates )
  {
    int constexpr n = 20;
    real64 constexpr minTemperature = 12.0 + constants::zeroDegreesCelsiusInKelvin;
    real64 constexpr maxTemperature = 120.0 + constants::zeroDegreesCelsiusInKelvin;
    real64 constexpr dT = (maxTemperature - minTemperature) / n;
    temperatureCoordinates.resize( n+1 );
    for( integer i = 0; i <= n; i++ )
    {
      temperatureCoordinates[i] = minTemperature + i*dT + 0.01;
    }
  }

protected:
  conduit::Node m_node;
  dataRepository::Group m_parent;
  std::unique_ptr< FunctionManager > m_functionManager{};
  std::unique_ptr< PhillipsBrineViscosity > m_viscosity{};
  std::unique_ptr< ModelParameters > m_parameters{};
};

using PhillipsBrineViscosity4 = PhillipsBrineViscosityTestFixture< TestFluidName::EnhancedOilRecovery >;
using PhillipsBrineViscosity2_100K = PhillipsBrineViscosityTestFixture< TestFluidName::CarbonStorage, 100000 >;

TEST_P( PhillipsBrineViscosity4, testViscosity )
{
  testViscosity( GetParam() );
}

TEST_P( PhillipsBrineViscosity2_100K, testViscosity )
{
  testViscosity( GetParam() );
}

TEST_P( PhillipsBrineViscosity4, testViscosityDerivatives )
{
  testViscosityDerivatives( GetParam() );
}

TEST_P( PhillipsBrineViscosity2_100K, testViscosityDerivatives )
{
  testViscosityDerivatives( GetParam() );
}

//-------------------------------------------------------------------------------
// Data
//-------------------------------------------------------------------------------

/* UNCRUSTIFY-OFF */

INSTANTIATE_TEST_SUITE_P(
  PhillipsBrineViscosity,
  PhillipsBrineViscosity4,
  ::testing::ValuesIn({
    ViscosityData<4>{{1.0e+07, 288.15, {0.0000, 1.0000, 0.0000, 0.0000}}, 1.153800e-03},
    ViscosityData<4>{{1.0e+07, 293.15, {0.0000, 1.0000, 0.0000, 0.0000}}, 1.012544e-03},
    ViscosityData<4>{{1.0e+07, 353.15, {0.0000, 1.0000, 0.0000, 0.0000}}, 3.552664e-04},
    ViscosityData<4>{{1.0e+07, 373.15, {0.0000, 1.0000, 0.0000, 0.0000}}, 2.822282e-04},
    ViscosityData<4>{{1.0e+07, 288.15, {0.0000, 0.9990, 0.0010, 0.0000}}, 1.153800e-03},
    ViscosityData<4>{{1.0e+07, 293.15, {0.0000, 0.9990, 0.0010, 0.0000}}, 1.012544e-03},
    ViscosityData<4>{{1.0e+07, 353.15, {0.0000, 0.9990, 0.0010, 0.0000}}, 3.552664e-04},
    ViscosityData<4>{{1.0e+07, 373.15, {0.0000, 0.9990, 0.0010, 0.0000}}, 2.822282e-04},
    ViscosityData<4>{{1.0e+07, 288.15, {0.0000, 0.9900, 0.0100, 0.0000}}, 1.153800e-03},
    ViscosityData<4>{{1.0e+07, 293.15, {0.0000, 0.9900, 0.0100, 0.0000}}, 1.012544e-03},
    ViscosityData<4>{{1.0e+07, 353.15, {0.0000, 0.9900, 0.0100, 0.0000}}, 3.552664e-04},
    ViscosityData<4>{{1.0e+07, 373.15, {0.0000, 0.9900, 0.0100, 0.0000}}, 2.822282e-04},
    ViscosityData<4>{{1.0e+07, 288.15, {0.0010, 0.9990, 0.0000, 0.0000}}, 1.153800e-03},
    ViscosityData<4>{{1.0e+07, 293.15, {0.0010, 0.9990, 0.0000, 0.0000}}, 1.012544e-03},
    ViscosityData<4>{{1.0e+07, 353.15, {0.0010, 0.9990, 0.0000, 0.0000}}, 3.552664e-04},
    ViscosityData<4>{{1.0e+07, 373.15, {0.0010, 0.9990, 0.0000, 0.0000}}, 2.822282e-04},
    ViscosityData<4>{{1.0e+07, 288.15, {0.0010, 0.9980, 0.0010, 0.0000}}, 1.153800e-03},
    ViscosityData<4>{{1.0e+07, 293.15, {0.0010, 0.9980, 0.0010, 0.0000}}, 1.012544e-03},
    ViscosityData<4>{{1.0e+07, 353.15, {0.0010, 0.9980, 0.0010, 0.0000}}, 3.552664e-04},
    ViscosityData<4>{{1.0e+07, 373.15, {0.0010, 0.9980, 0.0010, 0.0000}}, 2.822282e-04},
    ViscosityData<4>{{1.0e+07, 288.15, {0.0010, 0.9890, 0.0100, 0.0000}}, 1.153800e-03},
    ViscosityData<4>{{1.0e+07, 293.15, {0.0010, 0.9890, 0.0100, 0.0000}}, 1.012544e-03},
    ViscosityData<4>{{1.0e+07, 353.15, {0.0010, 0.9890, 0.0100, 0.0000}}, 3.552664e-04},
    ViscosityData<4>{{1.0e+07, 373.15, {0.0010, 0.9890, 0.0100, 0.0000}}, 2.822282e-04},
    ViscosityData<4>{{1.0e+07, 288.15, {0.0100, 0.9900, 0.0000, 0.0000}}, 1.153800e-03},
    ViscosityData<4>{{1.0e+07, 293.15, {0.0100, 0.9900, 0.0000, 0.0000}}, 1.012544e-03},
    ViscosityData<4>{{1.0e+07, 353.15, {0.0100, 0.9900, 0.0000, 0.0000}}, 3.552664e-04},
    ViscosityData<4>{{1.0e+07, 373.15, {0.0100, 0.9900, 0.0000, 0.0000}}, 2.822282e-04},
    ViscosityData<4>{{1.0e+07, 288.15, {0.0100, 0.9890, 0.0010, 0.0000}}, 1.153800e-03},
    ViscosityData<4>{{1.0e+07, 293.15, {0.0100, 0.9890, 0.0010, 0.0000}}, 1.012544e-03},
    ViscosityData<4>{{1.0e+07, 353.15, {0.0100, 0.9890, 0.0010, 0.0000}}, 3.552664e-04},
    ViscosityData<4>{{1.0e+07, 373.15, {0.0100, 0.9890, 0.0010, 0.0000}}, 2.822282e-04},
    ViscosityData<4>{{1.0e+07, 288.15, {0.0100, 0.9800, 0.0100, 0.0000}}, 1.153800e-03},
    ViscosityData<4>{{1.0e+07, 293.15, {0.0100, 0.9800, 0.0100, 0.0000}}, 1.012544e-03},
    ViscosityData<4>{{1.0e+07, 353.15, {0.0100, 0.9800, 0.0100, 0.0000}}, 3.552664e-04},
    ViscosityData<4>{{1.0e+07, 373.15, {0.0100, 0.9800, 0.0100, 0.0000}}, 2.822282e-04}
  })
);

INSTANTIATE_TEST_SUITE_P(
  PhillipsBrineViscosity,
  PhillipsBrineViscosity2_100K,
  ::testing::ValuesIn({
    ViscosityData<2>{{1.0e+07, 288.15, {1.0000, 0.0000}}, 1.364366e-03},
    ViscosityData<2>{{1.0e+07, 293.15, {1.0000, 0.0000}}, 1.199555e-03},
    ViscosityData<2>{{1.0e+07, 353.15, {1.0000, 0.0000}}, 4.302578e-04},
    ViscosityData<2>{{1.0e+07, 373.15, {1.0000, 0.0000}}, 3.442847e-04},
    ViscosityData<2>{{1.0e+07, 288.15, {0.9990, 0.0010}}, 1.364366e-03},
    ViscosityData<2>{{1.0e+07, 293.15, {0.9990, 0.0010}}, 1.199555e-03},
    ViscosityData<2>{{1.0e+07, 353.15, {0.9990, 0.0010}}, 4.302578e-04},
    ViscosityData<2>{{1.0e+07, 373.15, {0.9990, 0.0010}}, 3.442847e-04},
    ViscosityData<2>{{1.0e+07, 288.15, {0.9900, 0.0100}}, 1.364366e-03},
    ViscosityData<2>{{1.0e+07, 293.15, {0.9900, 0.0100}}, 1.199555e-03},
    ViscosityData<2>{{1.0e+07, 353.15, {0.9900, 0.0100}}, 4.302578e-04},
    ViscosityData<2>{{1.0e+07, 373.15, {0.9900, 0.0100}}, 3.442847e-04}
  })
);

/* UNCRUSTIFY-ON */

} // testing

} // geos
