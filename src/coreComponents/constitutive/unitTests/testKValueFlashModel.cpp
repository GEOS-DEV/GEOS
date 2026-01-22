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

#include "constitutive/fluid/multifluid/compositional/models/KValueFlashModel.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/KValueFlashParameters.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/PhaseType.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/EquationOfState.hpp"
#include "constitutive/fluid/multifluid/compositional/CompositionalMultiphaseFluid.hpp"
#include "constitutive/unitTests/TestFluid.hpp"
#include "constitutive/unitTests/TestFluidUtilities.hpp"

#include "functions/FunctionManager.hpp"
#include "functions/TableFunction.hpp"

#include <conduit.hpp>

using namespace geos::constitutive;
using namespace geos::constitutive::compositional;
using namespace geos::testing;

namespace geos
{

template< int numComps >
struct TestData;

template<>
struct TestData< 9 >
{
  static constexpr integer numTestComps = 3;
  static constexpr integer testComponents[numTestComps] = {0, 2, 8};
  static std::unique_ptr< TestFluid< 9 > > createFluid()
  {
    auto fluid = TestFluid< 9 >::create( {Fluid::H2O, Fluid::CO2, Fluid::N2, Fluid::CH4, Fluid::C2H6, Fluid::C3H8, Fluid::C4H10, Fluid::C5H12, Fluid::C10H22} );
    const std::array< real64 const, 36 > bics = {
      0.01, 0, 0.003732, 0, 0.01, 0, 0, 0.01, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 0, 0.01, 0, 0.028, 0.01, 0.01, 0, 0, 0.01, 0, 0.04532, 0.01, 0.01, 0, 0, 0
    };
    fluid->setBinaryCoefficients( bics );
    return fluid;
  }
};

template< integer numPhases, integer numComps >
using FlashData = std::tuple<
  real64 const,                                             // pressure
  real64 const,                                             // temperature
  Feed< numComps > const,                                   // total composition
  Feed< numPhases > const,                                  // expected phase fractions
  Feed< TestData< numComps >::numTestComps *numPhases >     // expected phase composition (2 selected components)
  >;

template< integer numPhases, integer numComps >
class KValueFlashTestFixture : public ::testing::TestWithParam< FlashData< numPhases, numComps > >
{
  static constexpr real64 relTol = 1.0e-5;
  static constexpr real64 absTol = 1.0e-7;
  static constexpr int numDofs = numComps + 2;
  using Deriv = geos::constitutive::multifluid::DerivativeOffset;
  using FlashModelType = KValueFlashModel< numPhases >;
  using FlashModelParamType = KValueFlashParameters< numPhases >;
  using PhasePropSlice = typename FlashModelType::KernelWrapper::PhaseProp::SliceType;
  using PhaseCompSlice = typename FlashModelType::KernelWrapper::PhaseComp::SliceType;
public:
  KValueFlashTestFixture();
  ~KValueFlashTestFixture() override;

  void testFlash( typename KValueFlashTestFixture::ParamType const & data );
  void testNumericalDerivative( typename KValueFlashTestFixture::ParamType const & data );

  class MockFluid;

protected:
  conduit::Node m_node;
  dataRepository::Group m_parent;
  std::unique_ptr< FunctionManager > m_functionManager{};
  std::unique_ptr< TestFluid< numComps > > m_fluid{};
  std::unique_ptr< ModelParameters > m_parameters{};
  std::unique_ptr< FlashModelType > m_flash{};
  array1d< integer > m_phaseTypes{};

private:
  void generateTables( string_array & names, string const fluidName );

  static void writeToFile( string const & fileName, string const & content );

  static void removeFile( string const & fileName );
private:
  string_array m_fileNames;
};

template< integer numPhases, integer numComps >
class KValueFlashTestFixture< numPhases, numComps >::MockFluid : public MultiFluidBase
{
public:
  MockFluid( Group * const parent ): MultiFluidBase( "fluid", parent ) {}
  string getCatalogName() const override { return ""; }
  void checkTablesParameters( real64, real64 ) const override {}
  integer getWaterPhaseIndex() const override { return 0; };

  void setProperties( ComponentProperties const & componentProperties )
  {
    string_array & phaseNames = getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );
    TestFluid< numPhases >::createArray( phaseNames, std::array< string, 2 >{"oil", "gas"} );
    string_array & componentNames = getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
    TestFluid< numPhases >::createArray( componentNames, componentProperties.getComponentName());
  }
};

template< integer numPhases, integer numComps >
KValueFlashTestFixture< numPhases, numComps >::KValueFlashTestFixture()
  : m_parent( "parent", m_node ),
  m_fluid( TestData< numComps >::createFluid() )
{
  m_functionManager = std::make_unique< FunctionManager >( "FunctionManager", &m_parent );

  ComponentProperties const & componentProperties = this->m_fluid->getComponentProperties();

  m_phaseTypes.emplace_back( static_cast< integer >(PhaseType::LIQUID));
  m_phaseTypes.emplace_back( static_cast< integer >(PhaseType::VAPOUR));
  if constexpr (numPhases == 3)
  {
    m_phaseTypes.emplace_back( static_cast< integer >(PhaseType::AQUEOUS));
  }

  string const fluidName = GEOS_FMT( "fluid_{}_{}", numPhases, numComps );
  m_parameters = FlashModelType::createParameters( std::move( m_parameters ));
  FlashModelParamType * parameters = const_cast< FlashModelParamType * >(m_parameters->get< FlashModelParamType >());
  parameters->m_kValueTables.resize( (numPhases-1)*numComps );
  generateTables( parameters->m_kValueTables, fluidName );

  MockFluid mockFluid( &m_parent );
  mockFluid.setProperties( componentProperties );
  m_parameters->postInputInitialization( &mockFluid, componentProperties );

  string const flashName = GEOS_FMT( "{}_flash", fluidName );
  m_flash = std::make_unique< FlashModelType >( flashName, componentProperties, *m_parameters, m_phaseTypes );
}

template< integer numPhases, integer numComps >
KValueFlashTestFixture< numPhases, numComps >::~KValueFlashTestFixture()
{
  for( string const & fileName : m_fileNames )
  {
    removeFile( fileName );
  }
}

// Crookston correlations pressure (bar), temperature (K)
real64 getKValue( integer const phaseIndex, integer const compIndex, real64 const pressure, real64 const temperature )
{
  GEOS_UNUSED_VAR( phaseIndex );

  static std::array< real64, 5 > constexpr crookstonCoefficients[9] = {
    {3.0620e+00, 8.9414e+02, 1.1912e-02, 5.3659e+02, 1.1951e+02},
    {-1.8141e+00, 6.2655e+02, 6.7489e-03, 5.0732e+00, 2.5249e+02},
    {-3.5742e-01, 4.8660e+02, 5.7887e-03, 9.1910e+01, 1.8027e+02},
    {3.6577e+00, 7.1889e+02, 1.4154e-02, 5.7323e+02, 1.1088e+02},
    {8.2845e+00, 1.0144e+03, 5.2968e-02, 9.5800e+02, 9.3548e+01},
    {1.5732e+01, 1.4968e+03, 1.7407e-01, 1.3707e+03, 8.1261e+01},
    {2.5068e+01, 2.2625e+03, 5.3888e-01, 1.7981e+03, 7.1948e+01},
    {3.6243e+01, 4.1985e+03, 1.9055e+00, 2.3575e+03, 5.8025e+01},
    {-2.3346e-01, 7.9356e-01, 6.8406e-03, 7.1217e+02, 1.3827e+02}
  };
  auto const [A, B, C, D, E] = std::apply( []( auto &&... args ) { return std::make_tuple( args ... ); }, crookstonCoefficients[compIndex] );
  real64 const kValue = (A + B/pressure + C*pressure)*LvArray::math::exp( -D/(temperature-E));
  return LvArray::math::max( kValue, 1.0e-8 );
}

template< integer numPhases, integer numComps >
void KValueFlashTestFixture< numPhases, numComps >::generateTables( string_array & names, string const fluidName )
{
  FunctionManager & functionManager = FunctionManager::getInstance();
  std::ostringstream content;
  for( integer ip = 0; ip < numPhases-1; ++ip )
  {
    for( integer ic = 0; ic < numComps; ++ic )
    {
      string const tableName = GEOS_FMT( "{}_KVALUE_{}_{}", fluidName, ip+1, ic+1 );

      // Generate tables whose end-points are not exactly the same
      real64 const th = (ip + 0.5)*numComps + ic + 0.5;
      real64 const c = LvArray::math::cos( th );
      real64 const s = LvArray::math::sin( th );
      // Pressure in bar
      real64 const minPressure = 1.0 - 0.1*c*c;
      real64 const maxPressure = 600.0 + 50.0*c*c;
      integer const NP = static_cast< integer >(6.0 + 8.0*s*s);
      real64 const dp = pow( maxPressure/minPressure, 1.0/NP );

      content.str( "" );
      real64 constexpr BAR_2_PA = 1.0e5;
      real64 pressure = minPressure;
      content << BAR_2_PA * pressure << "\n";
      for( integer i = 1; i <= NP; ++i )
      {
        pressure *= dp;
        content << BAR_2_PA * pressure << "\n";
      }
      string const pressureFileName = GEOS_FMT( "{}_PRESSURE.txt", tableName );
      writeToFile( pressureFileName, content.str());
      m_fileNames.emplace_back( pressureFileName );

      // Temp in degC
      real64 const minTemp = 10.0 + 5.0*s*s;
      real64 const maxTemp = 200.0 + 50.0*s*s;
      integer const NT = static_cast< integer >(6.0 + 8.0*c*c);
      real64 const dt = (maxTemp - minTemp)/NP;
      content.str( "" );
      for( integer j = 0; j <= NT; ++j )
      {
        real64 const temperature = minTemp + j*dt + 273.15;
        content << temperature << "\n";
      }
      string const temperatureFileName = GEOS_FMT( "{}_TEMP.txt", tableName );
      writeToFile( temperatureFileName, content.str());
      m_fileNames.emplace_back( temperatureFileName );

      content.str( "" );
      pressure = minPressure;
      for( integer i = 0; i <= NP; ++i )
      {
        for( integer j = 0; j <= NT; ++j )
        {
          real64 const temperature = minTemp + j*dt + 273.15;
          real64 const kValue = getKValue( ip, ic, pressure, temperature );
          content << kValue << "\n";
        }
        pressure *= dp;
      }
      string const kValueFileName = GEOS_FMT( "{}_KVALUE.txt", tableName );
      writeToFile( kValueFileName, content.str());
      m_fileNames.emplace_back( kValueFileName );

      TableFunction * tableFunction = dynamicCast< TableFunction * >( functionManager.createChild( "TableFunction", tableName ) );

      string_array & inputVarNames = tableFunction->getWrapper< string_array >( dataRepository::keys::inputVarNames ).reference();
      inputVarNames.emplace_back( "pressure" );
      inputVarNames.emplace_back( "temperature" );

      path_array & coordinateFiles = tableFunction->getWrapper< path_array >( TableFunction::viewKeyStruct::coordinateFilesString())
                                       .reference();
      coordinateFiles.emplace_back( Path() );
      coordinateFiles.emplace_back( Path() );
      dynamicCast< string & >( coordinateFiles[0] ) = pressureFileName;
      dynamicCast< string & >( coordinateFiles[1] ) = temperatureFileName;
      Path & voxelFile = tableFunction->getWrapper< Path >( TableFunction::viewKeyStruct::voxelFileString())
                           .reference();
      dynamicCast< string & >( voxelFile ) = kValueFileName;
      tableFunction->setInterpolationMethod( TableFunction::InterpolationType::Linear );
      tableFunction->initializeFunction();

      integer const tableIndex = numComps*ip + ic;

      names[tableIndex] = tableName;
    }
  }
}

template< integer numPhases, integer numComps >
void KValueFlashTestFixture< numPhases, numComps >::writeToFile( string const & fileName, string const & content )
{
  std::ofstream os( fileName );
  ASSERT_TRUE( os.is_open() );
  os << content;
  os.close();
}

template< integer numPhases, integer numComps >
void KValueFlashTestFixture< numPhases, numComps >::removeFile( string const & fileName )
{
  int const ret = std::remove( fileName.c_str() );
  ASSERT_TRUE( ret == 0 );
}

/* --- Start tests --- */
template< integer numPhases, integer numComps >
void KValueFlashTestFixture< numPhases, numComps >::testFlash( typename KValueFlashTestFixture::ParamType const & data )
{
  auto const & componentProperties = this->m_fluid->getComponentProperties().createKernelWrapper();

  auto flashKernelWrapper = this->m_flash->createKernelWrapper();

  real64 const pressure = std::get< 0 >( data );
  real64 const temperature = std::get< 1 >( data );
  stackArray1d< real64, numComps > composition;
  TestFluid< numComps >::createArray( composition, std::get< 2 >( data ));
  stackArray1d< real64, numPhases > expactedPhaseFractions;
  TestFluid< numComps >::createArray( expactedPhaseFractions, std::get< 3 >( data ));
  stackArray1d< real64, TestData< numComps >::numTestComps *numPhases > expactedPhaseCompositions;
  TestFluid< numComps >::createArray( expactedPhaseCompositions, std::get< 4 >( data ));

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

  flashKernelWrapper.compute( componentProperties,
                              pressure,
                              temperature,
                              composition.toSliceConst(),
                              kValues.toSlice(),
                              PhasePropSlice( phaseFraction, dPhaseFraction ),
                              PhaseCompSlice( phaseComponentFraction, dPhaseComponentFraction ) );

  integer constexpr nc = TestData< numComps >::numTestComps;
  for( integer ip = 0; ip < numPhases; ++ip )
  {
    checkRelativeError( phaseFraction[ip], expactedPhaseFractions[ip], relTol, absTol );
    for( integer i = 0; i < nc; ++i )
    {
      integer const ic = TestData< numComps >::testComponents[i];
      checkRelativeError( phaseComponentFraction[ip][ic], expactedPhaseCompositions[ip*nc+i], relTol, absTol );
    }
  }
}

template< integer numPhases, integer numComps >
void KValueFlashTestFixture< numPhases, numComps >::testNumericalDerivative( typename KValueFlashTestFixture::ParamType const & data )
{
  // Number of output values from each flash calculation
  constexpr integer numValues = numPhases * (1 + numComps);

  auto const & componentProperties = this->m_fluid->getComponentProperties().createKernelWrapper();

  auto flashKernelWrapper = this->m_flash->createKernelWrapper();

  real64 const pressure = std::get< 0 >( data );
  real64 const temperature = std::get< 1 >( data );
  stackArray1d< real64, numComps > composition;
  TestFluid< numComps >::createArray( composition, std::get< 2 >( data ));

  stackArray2d< real64, (numPhases-1)*numComps > kValues( numPhases-1, numComps );
  LvArray::forValuesInSlice( kValues.toSlice(), []( real64 & v ){ v = 0.0; } );

  stackArray1d< real64, numValues > derivatives( numValues );

  StackArray< real64, 3, numPhases, multifluid::LAYOUT_PHASE > phaseFractionData( 1, 1, numPhases );
  StackArray< real64, 4, numPhases *numDofs, multifluid::LAYOUT_PHASE_DC > dPhaseFractionData( 1, 1, numPhases, numDofs );
  StackArray< real64, 4, numPhases *numComps, multifluid::LAYOUT_PHASE_COMP > phaseComponentFractionData( 1, 1, numPhases, numComps );
  StackArray< real64, 5, numPhases *numComps *numDofs, multifluid::LAYOUT_PHASE_COMP_DC > dPhaseComponentFractionData( 1, 1, numPhases, numComps, numDofs );

  auto phaseFraction = phaseFractionData[0][0];
  auto dPhaseFraction = dPhaseFractionData[0][0];
  auto phaseComponentFraction = phaseComponentFractionData[0][0];
  auto dPhaseComponentFraction = dPhaseComponentFractionData[0][0];

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
  real64 const dT = 1.0e-4 * temperature;
  geos::testing::internal::testNumericalDerivative< numValues >(
    temperature, dT, derivatives,
    [&]( real64 const t, auto & values ) {
    evaluateFlash( pressure, t, composition.toSliceConst(), values );
  } );

  // -- Composition derivatives derivative
  real64 const dz = 1.0e-5;
  for( integer const ic : TestData< numComps >::testComponents )
  {
    if( composition[ic] < absTol || 1.0-composition[ic] < absTol )
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

using KValueFlashTest_2_9 = KValueFlashTestFixture< 2, 9 >;
using FlashData_2_9 = FlashData< 2, 9 >;
TEST_P( KValueFlashTest_2_9, testFlash )
{
  testFlash( GetParam() );
}

TEST_P( KValueFlashTest_2_9, testFlashNumericalDerivative )
{
  testNumericalDerivative( GetParam() );
}

/* UNCRUSTIFY-OFF */

// Test data
INSTANTIATE_TEST_SUITE_P(
  KValueFlash, KValueFlashTest_2_9,
  ::testing::ValuesIn< FlashData_2_9 >({
    {1.0e+05, 298.15, {0.000363, 0.000007, 0.003471, 0.006007, 0.018423, 0.034034, 0.042565, 0.056120, 0.839010}, {0.912604, 0.087396}, {0.000121, 0.000126, 0.914986, 0.002885, 0.038401, 0.045656}},
    {1.0e+05, 323.15, {0.000363, 0.000007, 0.003471, 0.006007, 0.018423, 0.034034, 0.042565, 0.056120, 0.839010}, {0.913798, 0.086202}, {0.000263, 0.000150, 0.912936, 0.001423, 0.038675, 0.055346}},
    {1.0e+05, 373.15, {0.000363, 0.000007, 0.003471, 0.006007, 0.018423, 0.034034, 0.042565, 0.056120, 0.839010}, {0.993677, 0.006323}, {0.000364, 0.002621, 0.844231, 0.000199, 0.137069, 0.018459}},
    {1.0e+06, 298.15, {0.000363, 0.000007, 0.003471, 0.006007, 0.018423, 0.034034, 0.042565, 0.056120, 0.839010}, {0.936231, 0.063769}, {0.000069, 0.000237, 0.891605, 0.004676, 0.050950, 0.066831}},
    {1.0e+06, 323.15, {0.000363, 0.000007, 0.003471, 0.006007, 0.018423, 0.034034, 0.042565, 0.056120, 0.839010}, {0.923713, 0.076287}, {0.000173, 0.000446, 0.902155, 0.002662, 0.040096, 0.074422}},
    {1.0e+06, 373.15, {0.000363, 0.000007, 0.003471, 0.006007, 0.018423, 0.034034, 0.042565, 0.056120, 0.839010}, {1.000000, 0.000000}, {0.000363, 0.003471, 0.839010, 0.000565, 0.144455, 0.000000}},
    {1.0e+07, 298.15, {0.000363, 0.000007, 0.003471, 0.006007, 0.018423, 0.034034, 0.042565, 0.056120, 0.839010}, {0.956084, 0.043916}, {0.000169, 0.000288, 0.876451, 0.004582, 0.072778, 0.023888}},
    {1.0e+07, 323.15, {0.000363, 0.000007, 0.003471, 0.006007, 0.018423, 0.034034, 0.042565, 0.056120, 0.839010}, {0.986650, 0.013350}, {0.000339, 0.001404, 0.850276, 0.002119, 0.156249, 0.006348}},
    {1.0e+07, 373.15, {0.000363, 0.000007, 0.003471, 0.006007, 0.018423, 0.034034, 0.042565, 0.056120, 0.839010}, {1.000000, 0.000000}, {0.000363, 0.003471, 0.839010, 0.000299, 0.077400, 0.000000}},
    {5.0e+07, 298.15, {0.000363, 0.000007, 0.003471, 0.006007, 0.018423, 0.034034, 0.042565, 0.056120, 0.839010}, {0.912044, 0.087956}, {0.000072, 0.000140, 0.916453, 0.003384, 0.038006, 0.035974}},
    {5.0e+07, 323.15, {0.000363, 0.000007, 0.003471, 0.006007, 0.018423, 0.034034, 0.042565, 0.056120, 0.839010}, {0.993464, 0.006536}, {0.000341, 0.001947, 0.844472, 0.003713, 0.235032, 0.008766}},
    {5.0e+07, 373.15, {0.000363, 0.000007, 0.003471, 0.006007, 0.018423, 0.034034, 0.042565, 0.056120, 0.839010}, {1.000000, 0.000000}, {0.000363, 0.003471, 0.839010, 0.000523, 0.087190, 0.000000}},
    {1.0e+05, 298.15, {0.000000, 0.000001, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.999999}, {1.000000, 0.000000}, {0.000000, 0.000000, 0.999999, 0.000000, 0.000000, 0.049822}},
    {5.0e+07, 373.15, {0.000000, 0.000001, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.999999}, {1.000000, 0.000000}, {0.000000, 0.000000, 0.999999, 0.000000, 0.000000, 0.000000}},
    {1.0e+05, 298.15, {0.100000, 0.000000, 0.899999, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000001}, {0.000000, 1.000000}, {0.004008, 0.002805, 0.993188, 0.100000, 0.899999, 0.000001}},
    {1.0e+05, 298.15, {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000}, {1.000000, 0.000000}, {0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 1.000000}},
    {1.0e+05, 323.15, {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000}, {1.000000, 0.000000}, {0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 1.000000}},
    {5.0e+07, 373.15, {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000}, {1.000000, 0.000000}, {0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 1.000000}},
    {1.0e+05, 298.15, {0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000}, {0.000000, 1.000000}, {0.000000, 1.000000, 0.000000, 0.000000, 1.000000, 0.000000}},
    {1.0e+05, 323.15, {0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000}, {0.000000, 1.000000}, {0.000000, 1.000000, 0.000000, 0.000000, 1.000000, 0.000000}},
    {5.0e+07, 323.15, {0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000}, {0.000000, 1.000000}, {0.000000, 1.000000, 0.000000, 0.000000, 1.000000, 0.000000}},
    {5.0e+07, 373.15, {0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000}, {0.000000, 1.000000}, {0.000000, 1.000000, 0.000000, 0.000000, 1.000000, 0.000000}}
  })
);

/* UNCRUSTIFY-ON */

} // namespace geos
