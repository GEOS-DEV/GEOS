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

#include "constitutiveTestHelpers.hpp"

#include "constitutive/fluid/multifluid/compositional/models/KValueFlashModel.hpp"
#include "constitutive/fluid/multifluid/compositional/models/KValueFlashParameters.hpp"
#include "constitutive/fluid/multifluid/compositional/models/EquationOfState.hpp"
#include "constitutive/fluid/multifluid/compositional/CompositionalMultiphaseFluid.hpp"
#include "constitutive/unitTests/TestFluid.hpp"
#include "constitutive/unitTests/TestFluidUtilities.hpp"

#include "functions/FunctionManager.hpp"
#include "functions/TableFunction.hpp"

#include "mainInterface/GeosxState.hpp"
#include "mainInterface/initialization.hpp"

using namespace geos::dataRepository;
using namespace geos::testing;
using namespace geos::constitutive;
using namespace geos::constitutive::compositional;

namespace geos
{

template< integer numPhases, integer numComps >
using FlashData = std::tuple<
  real64 const,                 // pressure
  real64 const,                 // temperature
  Feed< numComps > const,       // total composition
  Feed< numPhases > const,      // expected phase fractions
  Feed< 2*numPhases >           // expected phase composition (2 selected components)
  >;

template< int numComps >
struct TestData;

template<>
struct TestData< 9 >
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
};

template< integer numPhases, integer numComps >
class KValueFlashTestFixture : public ConstitutiveTestBase< MultiFluidBase >, public ::testing::WithParamInterface< FlashData< numPhases, numComps > >
{
  static constexpr real64 relTol = 1.0e-5;
  static constexpr real64 absTol = 1.0e-7;
  static constexpr int numDofs = numComps + 2;
  // Selected components for test
  static constexpr int comp0 = 0;
  static constexpr int comp1 = numComps-1;
  using Deriv = geos::constitutive::multifluid::DerivativeOffset;
public:
  KValueFlashTestFixture();
  ~KValueFlashTestFixture() override;

  void testFlash( typename KValueFlashTestFixture::ParamType const & data );

protected:
  std::unique_ptr< TestFluid< numComps > > m_fluid{};
  std::unique_ptr< KValueFlashModel< numPhases > > m_flash{};
  std::unique_ptr< ModelParameters > m_parameters{};

private:
  void generateTables( arraySlice1d< string > const & names, string const fluidName );

  static void writeToFile( string const & fileName, string const & content );

  static void removeFile( string const & fileName );

  static CompositionalTwoPhaseConstantViscosity * makeFluid( string const & name, Group * parent, TestFluid< numComps > const * testFluid );

private:
  string_array m_fileNames;
};

template< integer numPhases, integer numComps >
KValueFlashTestFixture< numPhases, numComps >::KValueFlashTestFixture():
  m_fluid( TestData< numComps >::createFluid() )
{
  using ModelParamType = KValueFlashParameters< numPhases >;

  string const fluidName = GEOS_FMT( "fluid_{}_{}", numPhases, numComps );

  m_parameters = KValueFlashModel< numPhases >::createParameters( std::move( m_parameters ));
  ModelParamType * parameters = const_cast< ModelParamType * >(m_parameters->get< ModelParamType >());
  parameters->m_kValueTables.resize( (numPhases-1)*numComps );
  generateTables( parameters->m_kValueTables.toSlice(), fluidName );

  ComponentProperties const & componentProperties = this->m_fluid->getComponentProperties();

  string const flashName = fluidName + "_flash";
  m_flash = std::make_unique< KValueFlashModel< numPhases > >( flashName, componentProperties, *m_parameters );

  auto & parent = this->m_parent;
  parent.resize( 1 );


  m_model = makeFluid( fluidName, &parent, m_fluid.get() );

  parent.initialize();
  parent.initializePostInitialConditions();

  m_parameters->postInputInitialization( m_model, componentProperties );
}

template< integer numPhases, integer numComps >
KValueFlashTestFixture< numPhases, numComps >::~KValueFlashTestFixture()
{
  for( string const & fileName : m_fileNames )
  {
    removeFile( fileName );
  }
}

template< integer numComps >
struct MakeFluid;

template<>
struct MakeFluid< 9 >
{
  static void populate( CompositionalTwoPhaseConstantViscosity & fluid, TestFluid< 9 > const * testFluid )
  {
    using FluidModel = CompositionalTwoPhaseConstantViscosity;

    string_array & componentNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::componentNamesString() );
    TestFluid< 9 >::createArray( componentNames, testFluid->componentNames );

    array1d< real64 > & molarWeight = fluid.getReference< array1d< real64 > >( MultiFluidBase::viewKeyStruct::componentMolarWeightString() );
    TestFluid< 9 >::createArray( molarWeight, testFluid->molecularWeight );

    array1d< real64 > & criticalPressure = fluid.getReference< array1d< real64 > >( FluidModel::viewKeyStruct::componentCriticalPressureString() );
    TestFluid< 9 >::createArray( criticalPressure, testFluid->criticalPressure );

    array1d< real64 > & criticalTemperature = fluid.getReference< array1d< real64 > >( FluidModel::viewKeyStruct::componentCriticalTemperatureString() );
    TestFluid< 9 >::createArray( criticalTemperature, testFluid->criticalTemperature );

    array1d< real64 > & acentricFactor = fluid.getReference< array1d< real64 > >( FluidModel::viewKeyStruct::componentAcentricFactorString() );
    TestFluid< 9 >::createArray( acentricFactor, testFluid->acentricFactor );
  }
};

template< integer numPhases, integer numComps >
CompositionalTwoPhaseConstantViscosity *
KValueFlashTestFixture< numPhases, numComps >::makeFluid( string const & name, Group * parent, TestFluid< numComps > const * testFluid )
{
  CompositionalTwoPhaseConstantViscosity & compositionalFluid = parent->registerGroup< CompositionalTwoPhaseConstantViscosity >( name );

  Group & fluid = compositionalFluid;

  auto & phaseNames = fluid.getReference< string_array >( MultiFluidBase::viewKeyStruct::phaseNamesString() );
  phaseNames.emplace_back( "gas" );
  phaseNames.emplace_back( "liquid" );

  string const eosName = EnumStrings< EquationOfStateType >::toString( EquationOfStateType::PengRobinson );
  string_array & equationsOfState = fluid.template getReference< string_array >( EquationOfState::viewKeyStruct::equationsOfStateString() );
  equationsOfState.emplace_back( eosName );
  equationsOfState.emplace_back( eosName );

  MakeFluid< numComps >::populate( compositionalFluid, testFluid );

  compositionalFluid.postInputInitializationRecursive();

  return &compositionalFluid;
}

// Crookston correlations pressure (bar), temperature (K)
real64 getKValue( integer const phaseIndex, integer const compIndex, real64 const pressure, real64 const temperature )
{
  ((void)phaseIndex);
  ((void)compIndex);
  real64 const A = 212.0;
  real64 const B = 10714.452833583635;
  real64 const C = 0.0;
  real64 const D = 2222.222222222222;
  real64 const E = 266.6666666666667;
  return (A + B / pressure + C * pressure ) * LvArray::math::exp( -D / (temperature - E ) );
}

template< integer numPhases, integer numComps >
void KValueFlashTestFixture< numPhases, numComps >::generateTables( arraySlice1d< string > const & names, string const fluidName )
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
      real64 const dp = (maxPressure - minPressure)/NP;
      content.str( "" );
      for( integer i = 0; i <= NP; ++i )
      {
        real64 const pressure = 1.0e5 * (minPressure + i*dp);
        content << pressure << "\n";
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
      for( integer i = 0; i <= NP; ++i )
      {
        real64 const pressure = minPressure + i*dp;
        for( integer j = 0; j <= NT; ++j )
        {
          real64 const temperature = minTemp + j*dt + 273.15;
          real64 const kValue = getKValue( ip, ic, pressure, temperature );
          content << kValue << "\n";
        }
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
  GEOS_UNUSED_VAR( data );
  auto flashKernelWrapper = this->m_flash->createKernelWrapper();
  GEOS_UNUSED_VAR( flashKernelWrapper );
}

using KValueFlashTest_2_9 = KValueFlashTestFixture< 2, 9 >;
TEST_P( KValueFlashTest_2_9, testFlash )
{
  testFlash( GetParam() );
}

/* UNCRUSTIFY-OFF */

// Test data
INSTANTIATE_TEST_SUITE_P(
  KValueFlash, KValueFlashTest_2_9,
  ::testing::Values( 
    FlashData<2, 9>( 1.0e+05, 278.15, {0.000363, 0.000007, 0.003471, 0.006007, 0.018423, 0.034034, 0.042565, 0.056120, 0.839010}, {0.000197, 0.999803}, {0.000361, 0.839175, 0.010852, 0.000016} )
  )
);

/* UNCRUSTIFY-ON */

} // namespace geos

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geos::GeosxState state( geos::basicSetup( argc, argv ) );

  int const result = RUN_ALL_TESTS();

  geos::basicCleanup();

  return result;
}
