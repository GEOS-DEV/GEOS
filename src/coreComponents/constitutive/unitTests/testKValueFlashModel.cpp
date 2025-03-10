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
#include "constitutive/fluid/multifluid/compositional/models/KValueFlashParameters.hpp"
#include "constitutive/fluid/multifluid/compositional/models/EquationOfState.hpp"
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
    auto fluid = TestFluid< 9 >::create( {Fluid::H2O, Fluid::CO2, Fluid::N2, Fluid::C1, Fluid::C2, Fluid::C3, Fluid::C4, Fluid::C5, Fluid::C10} );
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

  string const fluidName = GEOS_FMT( "fluid_{}_{}", numPhases, numComps );

  m_parameters = FlashModelType::createParameters( std::move( m_parameters ));
  FlashModelParamType * parameters = const_cast< FlashModelParamType * >(m_parameters->get< FlashModelParamType >());
  parameters->m_kValueTables.resize( (numPhases-1)*numComps );
  generateTables( parameters->m_kValueTables, fluidName );

  MockFluid mockFluid( &m_parent );
  mockFluid.setProperties( componentProperties );
  m_parameters->postInputInitialization( &mockFluid, componentProperties );

  string const flashName = GEOS_FMT( "{}_flash", fluidName );
  m_flash = std::make_unique< FlashModelType >( flashName, componentProperties, *m_parameters );
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
      real64 const minPressure = 1.0 + c*c;
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
  real64 const dT = 1.0e-6 * temperature;
  geos::testing::internal::testNumericalDerivative< numValues >(
    temperature, dT, derivatives,
    [&]( real64 const t, auto & values ) {
    evaluateFlash( pressure, t, composition.toSliceConst(), values );
  } );

  // -- Composition derivatives derivative
  real64 const dz = 1.0e-7;
  for( integer const ic : TestData< numComps >::testComponents )
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

using KValueFlashTest_2_9 = KValueFlashTestFixture< 2, 9 >;
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
  ::testing::Values(
    FlashData<2, 9>(1.0e+05, 298.15, {0.000363, 0.000007, 0.003471, 0.006007, 0.018423, 0.034034, 0.042565, 0.056120, 0.839010}, {0.944005, 0.055995}, {0.000168, 0.000302, 0.887904, 0.003657, 0.056891, 0.014717}),
    FlashData<2, 9>(1.0e+05, 323.15, {0.000363, 0.000007, 0.003471, 0.006007, 0.018423, 0.034034, 0.042565, 0.056120, 0.839010}, {0.943992, 0.056008}, {0.000295, 0.000345, 0.887764, 0.001502, 0.056153, 0.017287}),
    FlashData<2, 9>(1.0e+05, 373.15, {0.000363, 0.000007, 0.003471, 0.006007, 0.018423, 0.034034, 0.042565, 0.056120, 0.839010}, {1.000000, 0.000000}, {0.000363, 0.003471, 0.839010, 0.000363, 0.003471, 0.839010}),
    FlashData<2, 9>(1.0e+06, 298.15, {0.000363, 0.000007, 0.003471, 0.006007, 0.018423, 0.034034, 0.042565, 0.056120, 0.839010}, {0.915848, 0.084152}, {0.000060, 0.000300, 0.914123, 0.003662, 0.037987, 0.021530}),
    FlashData<2, 9>(1.0e+06, 323.15, {0.000363, 0.000007, 0.003471, 0.006007, 0.018423, 0.034034, 0.042565, 0.056120, 0.839010}, {0.925278, 0.074722}, {0.000182, 0.000685, 0.904929, 0.002600, 0.037969, 0.022737}),
    FlashData<2, 9>(1.0e+06, 373.15, {0.000363, 0.000007, 0.003471, 0.006007, 0.018423, 0.034034, 0.042565, 0.056120, 0.839010}, {1.000000, 0.000000}, {0.000363, 0.003471, 0.839010, 0.000363, 0.003471, 0.839010}),
    FlashData<2, 9>(1.0e+07, 298.15, {0.000363, 0.000007, 0.003471, 0.006007, 0.018423, 0.034034, 0.042565, 0.056120, 0.839010}, {0.980687, 0.019313}, {0.000248, 0.000876, 0.855270, 0.006185, 0.135236, 0.013360}),
    FlashData<2, 9>(1.0e+07, 323.15, {0.000363, 0.000007, 0.003471, 0.006007, 0.018423, 0.034034, 0.042565, 0.056120, 0.839010}, {0.985183, 0.014817}, {0.000339, 0.001692, 0.851514, 0.001990, 0.121759, 0.007644}),
    FlashData<2, 9>(1.0e+07, 373.15, {0.000363, 0.000007, 0.003471, 0.006007, 0.018423, 0.034034, 0.042565, 0.056120, 0.839010}, {1.000000, 0.000000}, {0.000363, 0.003471, 0.839010, 0.000363, 0.003471, 0.839010}),
    FlashData<2, 9>(5.0e+07, 298.15, {0.000363, 0.000007, 0.003471, 0.006007, 0.018423, 0.034034, 0.042565, 0.056120, 0.839010}, {0.945674, 0.054326}, {0.000109, 0.000348, 0.886566, 0.004782, 0.057827, 0.011188}),
    FlashData<2, 9>(5.0e+07, 323.15, {0.000363, 0.000007, 0.003471, 0.006007, 0.018423, 0.034034, 0.042565, 0.056120, 0.839010}, {1.000000, 0.000000}, {0.000363, 0.003471, 0.839010, 0.000363, 0.003471, 0.839010}),
    FlashData<2, 9>(5.0e+07, 373.15, {0.000363, 0.000007, 0.003471, 0.006007, 0.018423, 0.034034, 0.042565, 0.056120, 0.839010}, {1.000000, 0.000000}, {0.000363, 0.003471, 0.839010, 0.000363, 0.003471, 0.839010}),
    FlashData<2, 9>(1.0e+05, 298.15, {0.000000, 0.000001, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.999999}, {1.000000, 0.000000}, {0.000000, 0.000000, 0.999999, 0.000000, 0.000000, 0.999999}),
    FlashData<2, 9>(5.0e+07, 373.15, {0.000000, 0.000001, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.999999}, {1.000000, 0.000000}, {0.000000, 0.000000, 0.999999, 0.000000, 0.000000, 0.999999}),
    FlashData<2, 9>(1.0e+05, 298.15, {0.100000, 0.000000, 0.899999, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000001}, {0.000000, 1.000000}, {0.100000, 0.899999, 0.000001, 0.100000, 0.899999, 0.000001}),
    FlashData<2, 9>(5.0e+07, 373.15, {0.100000, 0.000000, 0.899999, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000001}, {0.000001, 0.999999}, {0.069963, 0.052686, 0.877351, 0.100000, 0.900000, 0.000000}),
    FlashData<2, 9>(1.0e+05, 298.15, {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000}, {1.000000, 0.000000}, {0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 1.000000}),
    FlashData<2, 9>(1.0e+05, 323.15, {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000}, {1.000000, 0.000000}, {0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 1.000000}),
    FlashData<2, 9>(5.0e+07, 373.15, {0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000}, {1.000000, 0.000000}, {0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 1.000000}),
    FlashData<2, 9>(1.0e+05, 298.15, {0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000}, {0.000000, 1.000000}, {0.000000, 1.000000, 0.000000, 0.000000, 1.000000, 0.000000}),
    FlashData<2, 9>(1.0e+05, 323.15, {0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000}, {0.000000, 1.000000}, {0.000000, 1.000000, 0.000000, 0.000000, 1.000000, 0.000000}),
    FlashData<2, 9>(5.0e+07, 323.15, {0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000}, {0.000000, 1.000000}, {0.000000, 1.000000, 0.000000, 0.000000, 1.000000, 0.000000}),
    FlashData<2, 9>(5.0e+07, 373.15, {0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000}, {0.000000, 1.000000}, {0.000000, 1.000000, 0.000000, 0.000000, 1.000000, 0.000000})
  )
);

/* UNCRUSTIFY-ON */

} // namespace geos
