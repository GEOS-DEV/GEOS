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
#include "TestFluidFixture.hpp"
#include "TestFluidUtilities.hpp"

using namespace geos::constitutive::compositional;

namespace geos
{
namespace testing
{

template< integer NC >
struct ViscosityData : public TestDataBase< NC >
{
  integer mixingType;       ///< Mixing type
  real64 expectedViscosity; ///< Expected viscosity
};

template< integer NC >
std::ostream & operator<<( std::ostream & os, ViscosityData< NC > const & data )
{
  return TestDataPrinter( os, data ).value( data.mixingType )
                                    .output( data.expectedViscosity );
}

template< TestFluidName FLUID >
class LohrenzBrayClarkViscosityTestFixture : public TestFluidFixture< FLUID, ViscosityData >
{
  using Base = TestFluidFixture< FLUID, ViscosityData >;
  using Base::numComps;
  using Base::numDofs;
  using Base::relTol;
  using Base::absTol;
  using typename Base::Deriv;

public:
  LohrenzBrayClarkViscosityTestFixture()
  {
    ComponentProperties const & componentProperties = this->m_fluid->getComponentProperties();

    m_parameters = CompositionalDensity::createParameters( std::make_unique< ModelParameters >() );
    m_parameters = LohrenzBrayClarkViscosity::createParameters( std::move( m_parameters ) );

    auto * parameters = const_cast< CriticalVolume * >(m_parameters->get< CriticalVolume >());
    parameters->m_componentCriticalVolume.resize( numComps );
    TestFluid< numComps >::populateArray( parameters->m_componentCriticalVolume, this->m_fluid->criticalVolume );

    auto * equationOfState = const_cast< EquationOfState * >(m_parameters->get< EquationOfState >());
    string const eosName = EnumStrings< EquationOfStateType >::toString( EquationOfStateType::PengRobinson );
    equationOfState->m_equationsOfStateNames.emplace_back( eosName );

    m_density = std::make_unique< CompositionalDensity >( "PhaseDensity", componentProperties, 0, *m_parameters );
    m_viscosity = std::make_unique< LohrenzBrayClarkViscosity >( "PhaseViscosity", componentProperties, 0, *m_parameters );
  }

  ~LohrenzBrayClarkViscosityTestFixture() = default;

  void testViscosity( ViscosityData< numComps > const & data )
  {
    auto const mixing_type = static_cast< LohrenzBrayClarkViscosityUpdate::MixingType >(data.mixingType );
    real64 const pressure = data.pressure;
    real64 const temperature = data.temperature;
    stackArray1d< real64, numComps > phaseComposition;
    TestFluid< numComps >::createArray( phaseComposition, data.composition );
    real64 const expectedViscosity = data.expectedViscosity;

    real64 molarDensity = 0.0;
    real64 massDensity = 0.0;
    real64 viscosity = 0.0;
    stackArray1d< real64, numDofs > tempDerivs( numDofs );

    auto componentProperties = this->m_fluid->createKernelWrapper();
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

  void testViscosityDerivatives( ViscosityData< numComps > const & data )
  {
    stackArray1d< real64, numComps > phaseComposition;
    TestFluid< numComps >::createArray( phaseComposition, data.composition );

    auto componentProperties = this->m_fluid->createKernelWrapper();
    auto densityKernelWrapper = m_density->createKernelWrapper();
    auto viscosityKernelWrapper = m_viscosity->createKernelWrapper();
    viscosityKernelWrapper.setMixingType( static_cast< LohrenzBrayClarkViscosityUpdate::MixingType >( data.mixingType ) );

    auto computeViscosity = [&]( real64 const p, real64 const t, arraySlice1d< real64 const > const & z,
                                 real64 & viscosity, arraySlice1d< real64 > const & viscosityDerivs )
    {
      real64 molarDensity = 0.0;
      real64 massDensity = 0.0;
      stackArray1d< real64, numDofs > molarDensityDerivs( numDofs );
      stackArray1d< real64, numDofs > massDensityDerivs( numDofs );
      densityKernelWrapper.compute( componentProperties, p, t, z,
                                    molarDensity, molarDensityDerivs.toSlice(), massDensity, massDensityDerivs.toSlice(), false );
      viscosityKernelWrapper.compute( componentProperties, p, t, z,
                                      massDensity, massDensityDerivs.toSliceConst(), viscosity, viscosityDerivs, false );
    };

    real64 viscosity = 0.0;
    stackArray2d< real64, numDofs > derivatives( 1, numDofs );
    computeViscosity( data.pressure, data.temperature, phaseComposition.toSliceConst(), viscosity, derivatives[0] );

    // Viscosity values are very small so they are inflated to avoid false positives due
    // to the absolute value check
    stackArray1d< real64, numDofs > tempDerivs( numDofs );
    Base::template checkNumericalDerivatives< 1 >( data, derivatives.toSliceConst(),
                                                   [&]( real64 const p, real64 const t, auto const & z, auto & values )
    {
      computeViscosity( p, t, z, values[0], tempDerivs.toSlice() );
    }, {}, { 1.0e6 } );
  }

protected:
  std::unique_ptr< CompositionalDensity > m_density{};
  std::unique_ptr< LohrenzBrayClarkViscosity > m_viscosity{};
  std::unique_ptr< ModelParameters > m_parameters{};
};

using LohrenzBrayClarkViscosity9 = LohrenzBrayClarkViscosityTestFixture< TestFluidName::VolatileOil >;

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
    {{1.839590e+06, 297.15, {0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920}}, 0, 1.277311e-05},
    {{1.839590e+06, 297.15, {0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920}}, 1, 1.110154e-05},
    {{1.839590e+06, 297.15, {0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920}}, 2, 1.570131e-05},
    {{1.839590e+06, 297.15, {0.00826, 0.00544, 0.77032, 0.10456, 0.06177, 0.02459, 0.00884, 0.00472, 0.01149}}, 0, 1.598556e-05},
    {{1.839590e+06, 297.15, {0.00826, 0.00544, 0.77032, 0.10456, 0.06177, 0.02459, 0.00884, 0.00472, 0.01149}}, 1, 1.525278e-05},
    {{1.839590e+06, 297.15, {0.00826, 0.00544, 0.77032, 0.10456, 0.06177, 0.02459, 0.00884, 0.00472, 0.01149}}, 2, 1.573849e-05},
    {{1.839590e+06, 297.15, {0.00899, 0.00299, 0.53281, 0.11447, 0.08791, 0.04566, 0.02095, 0.01516, 0.17107}}, 0, 1.275691e-05},
    {{1.839590e+06, 297.15, {0.00899, 0.00299, 0.53281, 0.11447, 0.08791, 0.04566, 0.02095, 0.01516, 0.17107}}, 1, 1.108572e-05},
    {{1.839590e+06, 297.15, {0.00899, 0.00299, 0.53281, 0.11447, 0.08791, 0.04566, 0.02095, 0.01516, 0.17107}}, 2, 1.570641e-05},
    {{1.839590e+06, 363.00, {0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920}}, 0, 1.441390e-05},
    {{1.839590e+06, 363.00, {0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920}}, 1, 1.261317e-05},
    {{1.839590e+06, 363.00, {0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920}}, 2, 1.790554e-05},
    {{1.839590e+06, 363.00, {0.00826, 0.00544, 0.77032, 0.10456, 0.06177, 0.02459, 0.00884, 0.00472, 0.01149}}, 0, 1.843194e-05},
    {{1.839590e+06, 363.00, {0.00826, 0.00544, 0.77032, 0.10456, 0.06177, 0.02459, 0.00884, 0.00472, 0.01149}}, 1, 1.766697e-05},
    {{1.839590e+06, 363.00, {0.00826, 0.00544, 0.77032, 0.10456, 0.06177, 0.02459, 0.00884, 0.00472, 0.01149}}, 2, 1.818931e-05},
    {{1.839590e+06, 363.00, {0.00899, 0.00299, 0.53281, 0.11447, 0.08791, 0.04566, 0.02095, 0.01516, 0.17107}}, 0, 1.438427e-05},
    {{1.839590e+06, 363.00, {0.00899, 0.00299, 0.53281, 0.11447, 0.08791, 0.04566, 0.02095, 0.01516, 0.17107}}, 1, 1.258368e-05},
    {{1.839590e+06, 363.00, {0.00899, 0.00299, 0.53281, 0.11447, 0.08791, 0.04566, 0.02095, 0.01516, 0.17107}}, 2, 1.790103e-05},
    {{1.839590e+08, 297.15, {0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920}}, 0, 4.180223e-04},
    {{1.839590e+08, 297.15, {0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920}}, 1, 4.163507e-04},
    {{1.839590e+08, 297.15, {0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920}}, 2, 4.209505e-04},
    {{1.839590e+08, 297.15, {0.00826, 0.00544, 0.77032, 0.10456, 0.06177, 0.02459, 0.00884, 0.00472, 0.01149}}, 0, 1.823560e-04},
    {{1.839590e+08, 297.15, {0.00826, 0.00544, 0.77032, 0.10456, 0.06177, 0.02459, 0.00884, 0.00472, 0.01149}}, 1, 1.816232e-04},
    {{1.839590e+08, 297.15, {0.00826, 0.00544, 0.77032, 0.10456, 0.06177, 0.02459, 0.00884, 0.00472, 0.01149}}, 2, 1.821089e-04},
    {{1.839590e+08, 297.15, {0.00899, 0.00299, 0.53281, 0.11447, 0.08791, 0.04566, 0.02095, 0.01516, 0.17107}}, 0, 4.200580e-04},
    {{1.839590e+08, 297.15, {0.00899, 0.00299, 0.53281, 0.11447, 0.08791, 0.04566, 0.02095, 0.01516, 0.17107}}, 1, 4.183868e-04},
    {{1.839590e+08, 297.15, {0.00899, 0.00299, 0.53281, 0.11447, 0.08791, 0.04566, 0.02095, 0.01516, 0.17107}}, 2, 4.230075e-04},
    {{1.839590e+08, 363.00, {0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920}}, 0, 3.016154e-04},
    {{1.839590e+08, 363.00, {0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920}}, 1, 2.998146e-04},
    {{1.839590e+08, 363.00, {0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920}}, 2, 3.051070e-04},
    {{1.839590e+08, 363.00, {0.00826, 0.00544, 0.77032, 0.10456, 0.06177, 0.02459, 0.00884, 0.00472, 0.01149}}, 0, 1.225458e-04},
    {{1.839590e+08, 363.00, {0.00826, 0.00544, 0.77032, 0.10456, 0.06177, 0.02459, 0.00884, 0.00472, 0.01149}}, 1, 1.217808e-04},
    {{1.839590e+08, 363.00, {0.00826, 0.00544, 0.77032, 0.10456, 0.06177, 0.02459, 0.00884, 0.00472, 0.01149}}, 2, 1.223032e-04},
    {{1.839590e+08, 363.00, {0.00899, 0.00299, 0.53281, 0.11447, 0.08791, 0.04566, 0.02095, 0.01516, 0.17107}}, 0, 3.034493e-04},
    {{1.839590e+08, 363.00, {0.00899, 0.00299, 0.53281, 0.11447, 0.08791, 0.04566, 0.02095, 0.01516, 0.17107}}, 1, 3.016487e-04},
    {{1.839590e+08, 363.00, {0.00899, 0.00299, 0.53281, 0.11447, 0.08791, 0.04566, 0.02095, 0.01516, 0.17107}}, 2, 3.069661e-04}
  } )
);

/* UNCRUSTIFY-ON */

} // testing
} // geos
