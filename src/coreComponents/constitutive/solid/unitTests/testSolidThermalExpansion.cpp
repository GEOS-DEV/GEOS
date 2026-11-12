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

/**
 * @file testSolidThermalExpansion.cpp
 */

#include "gtest/gtest.h"

#include "common/GEOS_RAJA_Interface.hpp"
#include "common/initializeEnvironment.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "constitutive/solid/SolidThermalExpansion.hpp"
#include "dataRepository/xmlWrapper.hpp"
#include "functions/FunctionManager.hpp"
#include "functions/TableFunction.hpp"

using namespace geos;
using namespace geos::constitutive;

namespace
{

using Strategy = SolidThermalExpansionUpdates::Strategy;

/// Table used by the tabulated tests: knots and values
stdVector< real64 > const tableKnots = { 0.0, 20.0, 50.0, 90.0 };
stdVector< real64 > const tableValues = { 2.0e-6, 5.0e-6, 6.0e-6, 1.6e-5 };

/// Reference alpha(T) of the table: linear interpolation, constant outside
real64 referenceTableAlpha( real64 const T )
{
  if( T <= tableKnots.front() )
  {
    return tableValues.front();
  }
  if( T >= tableKnots.back() )
  {
    return tableValues.back();
  }
  for( localIndex i = 0; i + 1 < LvArray::integerConversion< localIndex >( tableKnots.size() ); ++i )
  {
    if( T < tableKnots[i + 1] )
    {
      real64 const w = ( T - tableKnots[i] ) / ( tableKnots[i + 1] - tableKnots[i] );
      return ( 1.0 - w ) * tableValues[i] + w * tableValues[i + 1];
    }
  }
  return tableValues.back();
}

/// Independent reference mean value of alpha over [a, b] (composite midpoint rule on a fine grid)
template< typename ALPHA >
real64 referenceMean( ALPHA && alphaOf, real64 const a, real64 const b )
{
  if( !( a < b ) && !( b < a ) )
  {
    return alphaOf( a );
  }
  integer constexpr numSteps = 200000;
  real64 const h = ( b - a ) / numSteps;
  real64 sum = 0.0;
  for( integer i = 0; i < numSteps; ++i )
  {
    sum += alphaOf( a + ( i + 0.5 ) * h );
  }
  return sum / numSteps;
}

/// Test fixture: a ConstitutiveManager and a FunctionManager under a root group
class SolidThermalExpansionTest : public ::testing::Test
{
public:
  SolidThermalExpansionTest():
    m_node(),
    m_root( "root", m_node ),
    m_functionManager( std::make_unique< FunctionManager >( "FunctionManager", &m_root ) )
  {}

  /// Register a TableFunction of temperature
  void addTable( string const & name,
                 stdVector< real64 > const & knots,
                 stdVector< real64 > const & values,
                 TableFunction::InterpolationType const interpolation = TableFunction::InterpolationType::Linear )
  {
    TableFunction & table = m_functionManager->registerGroup< TableFunction >( name );
    array1d< real64_array > coordinates( 1 );
    for( real64 const t : knots )
    {
      coordinates[0].emplace_back( t );
    }
    real64_array tableValuesArray;
    for( real64 const v : values )
    {
      tableValuesArray.emplace_back( v );
    }
    table.setTableCoordinates( coordinates );
    table.setTableValues( tableValuesArray );
    table.setInterpolationMethod( interpolation );
    table.initializeFunction();
  }

  /// Parse a <Constitutive> block in a fresh ConstitutiveManager and run the input checks of a real run
  ConstitutiveManager & parseManager( string const & constitutiveBlock )
  {
    m_constitutiveManagers.emplace_back(
      std::make_unique< ConstitutiveManager >( "constitutive" + std::to_string( m_constitutiveManagers.size() ), &m_root ) );
    ConstitutiveManager & constitutiveManager = *m_constitutiveManagers.back();

    string const inputStream = "<Constitutive>" + constitutiveBlock + "</Constitutive>";
    xmlWrapper::xmlDocument xmlDocument;
    xmlWrapper::xmlResult const xmlResult = xmlDocument.loadString( inputStream );
    EXPECT_TRUE( xmlResult );

    xmlWrapper::xmlNode xmlConstitutiveNode = xmlDocument.getChild( "Constitutive" );
    constitutiveManager.processInputFileRecursive( xmlDocument, xmlConstitutiveNode );
    constitutiveManager.postInputInitializationRecursive();
    return constitutiveManager;
  }

  /// Parse a <Constitutive> block and run the initialization steps of a real run on the model named @p name
  SolidThermalExpansion & parse( string const & constitutiveBlock, string const & name )
  {
    SolidThermalExpansion & model = parseManager( constitutiveBlock ).getConstitutiveRelation< SolidThermalExpansion >( name );
    model.initialize();
    return model;
  }

  /// Run @p func and check that it throws an InputError whose message contains @p expectedMessage
  template< typename FUNC >
  static void expectInputError( FUNC && func, string const & expectedMessage )
  {
    try
    {
      func();
      ADD_FAILURE() << "No InputError thrown, expected: " << expectedMessage;
    }
    catch( InputError const & e )
    {
      EXPECT_NE( string( e.what() ).find( expectedMessage ), string::npos ) << "Unexpected message: " << e.what();
    }
  }

protected:
  conduit::Node m_node;
  dataRepository::Group m_root;
  std::unique_ptr< FunctionManager > m_functionManager;
  stdVector< std::unique_ptr< ConstitutiveManager > > m_constitutiveManagers;
};

/// Evaluate the wrapper on a list of intervals with the given policy
template< typename POLICY >
void evaluate( SolidThermalExpansionUpdates const & wrapper,
               stdVector< std::pair< real64, real64 > > const & intervals,
               array1d< real64 > & meanAlpha,
               array1d< real64 > & alpha )
{
  localIndex const n = LvArray::integerConversion< localIndex >( intervals.size() );
  array1d< real64 > lower( n );
  array1d< real64 > upper( n );
  for( localIndex i = 0; i < n; ++i )
  {
    lower[i] = intervals[i].first;
    upper[i] = intervals[i].second;
  }
  meanAlpha.resize( n );
  alpha.resize( n );

  arrayView1d< real64 const > const lowerView = lower.toViewConst();
  arrayView1d< real64 const > const upperView = upper.toViewConst();
  arrayView1d< real64 > const meanAlphaView = meanAlpha.toView();
  arrayView1d< real64 > const alphaView = alpha.toView();

  forAll< POLICY >( n, [=] GEOS_HOST_DEVICE ( localIndex const i )
  {
    wrapper.compute( 0, lowerView[i], upperView[i], meanAlphaView[i], alphaView[i] );
  } );

  meanAlpha.move( hostMemorySpace, false );
  alpha.move( hostMemorySpace, false );
}

/// Intervals exercising one segment, several knots, both sides of the table, reversed and empty intervals
stdVector< std::pair< real64, real64 > > const testIntervals =
{
  { 25.0, 40.0 },   // inside one segment
  { 40.0, 25.0 },   // reversed
  { 30.0, 30.0 },   // empty interval
  { 20.0, 50.0 },   // exactly between two knots
  { 10.0, 70.0 },   // across two knots
  { 5.0, 85.0 },    // across three knots
  { -30.0, -10.0 }, // below the table
  { 100.0, 130.0 }, // above the table
  { -20.0, 120.0 }, // across the whole table and beyond both ends
  { 90.0, 60.0 },   // from the last knot, reversed
  { 19.999999, 20.000001 } // tiny interval across a knot
};

real64 relativeError( real64 const value, real64 const reference )
{
  return LvArray::math::abs( value - reference ) / LvArray::math::abs( reference );
}

template< typename POLICY >
void checkTabulated( SolidThermalExpansionUpdates const & wrapper )
{
  EXPECT_EQ( wrapper.getStrategy(), Strategy::Tabulated );

  array1d< real64 > meanAlpha, alpha;
  evaluate< POLICY >( wrapper, testIntervals, meanAlpha, alpha );

  for( localIndex i = 0; i < LvArray::integerConversion< localIndex >( testIntervals.size() ); ++i )
  {
    real64 const a = testIntervals[i].first;
    real64 const b = testIntervals[i].second;
    SCOPED_TRACE( "interval [" + std::to_string( a ) + ", " + std::to_string( b ) + "]" );
    EXPECT_LT( relativeError( alpha[i], referenceTableAlpha( b ) ), 1.0e-12 );
    EXPECT_LT( relativeError( meanAlpha[i], referenceMean( referenceTableAlpha, a, b ) ), 1.0e-8 );
  }
}

template< typename POLICY >
void checkLinear( SolidThermalExpansionUpdates const & wrapper,
                  real64 const alpha0, real64 const slope, real64 const T0 )
{
  EXPECT_EQ( wrapper.getStrategy(), Strategy::Linear );

  auto const linearAlpha = [=]( real64 const T ) { return alpha0 + slope * ( T - T0 ); };

  array1d< real64 > meanAlpha, alpha;
  evaluate< POLICY >( wrapper, testIntervals, meanAlpha, alpha );

  for( localIndex i = 0; i < LvArray::integerConversion< localIndex >( testIntervals.size() ); ++i )
  {
    real64 const a = testIntervals[i].first;
    real64 const b = testIntervals[i].second;
    SCOPED_TRACE( "interval [" + std::to_string( a ) + ", " + std::to_string( b ) + "]" );
    EXPECT_LT( relativeError( alpha[i], linearAlpha( b ) ), 1.0e-12 );
    // closed form: the mean of a linear function is its value at the midpoint
    EXPECT_LT( relativeError( meanAlpha[i], linearAlpha( 0.5 * ( a + b ) ) ), 1.0e-12 );
  }
}

/// d/dT [ meanAlpha( Tn, T ) * ( T - Tn ) ] must equal alpha( T )
void checkJacobian( SolidThermalExpansionUpdates const & wrapper, real64 const Tn, real64 const T )
{
  auto const strain = [&]( real64 const temperature )
  {
    real64 meanAlpha = 0.0, alpha = 0.0;
    wrapper.compute( 0, Tn, temperature, meanAlpha, alpha );
    return meanAlpha * ( temperature - Tn );
  };
  real64 const h = 1.0e-4;
  real64 const finiteDifference = ( strain( T + h ) - strain( T - h ) ) / ( 2.0 * h );
  real64 meanAlpha = 0.0, alpha = 0.0;
  wrapper.compute( 0, Tn, T, meanAlpha, alpha );
  SCOPED_TRACE( "Tn = " + std::to_string( Tn ) + ", T = " + std::to_string( T ) );
  EXPECT_LT( relativeError( finiteDifference, alpha ), 1.0e-7 );
}

} // namespace

TEST_F( SolidThermalExpansionTest, linear )
{
  SolidThermalExpansion & model =
    parse( "<SolidThermalExpansion name=\"tec\" referenceDrainedLinearTEC=\"1.0e-5\" dDrainedLinearTEC_dTemperature=\"2.0e-8\" referenceTemperature=\"20.0\"/>", "tec" );
  SolidThermalExpansionUpdates const wrapper = model.createKernelUpdates();

  checkLinear< serialPolicy >( wrapper, 1.0e-5, 2.0e-8, 20.0 );
  checkLinear< parallelDevicePolicy<> >( wrapper, 1.0e-5, 2.0e-8, 20.0 );

  checkJacobian( wrapper, 20.0, 60.0 );
  checkJacobian( wrapper, 80.0, 10.0 );
}

TEST_F( SolidThermalExpansionTest, constant )
{
  SolidThermalExpansion & model = parse( "<SolidThermalExpansion name=\"tec\" referenceDrainedLinearTEC=\"1.5e-5\"/>", "tec" );
  SolidThermalExpansionUpdates const wrapper = model.createKernelUpdates();

  array1d< real64 > meanAlpha, alpha;
  evaluate< parallelDevicePolicy<> >( wrapper, testIntervals, meanAlpha, alpha );
  for( localIndex i = 0; i < meanAlpha.size(); ++i )
  {
    // a constant coefficient must be returned exactly, as in the legacy code path
    EXPECT_DOUBLE_EQ( meanAlpha[i], 1.5e-5 );
    EXPECT_DOUBLE_EQ( alpha[i], 1.5e-5 );
  }
}

TEST_F( SolidThermalExpansionTest, fromSolidBase )
{
  array1d< real64 > perElementAlpha( 3 );
  perElementAlpha[0] = 1.0e-5;
  perElementAlpha[1] = 2.0e-5;
  perElementAlpha[2] = 3.0e-5;

  SolidThermalExpansionUpdates const wrapper( perElementAlpha.toViewConst() );
  EXPECT_EQ( wrapper.getStrategy(), Strategy::FromSolidBase );

  array1d< real64 > meanAlpha( 3 ), alpha( 3 );
  arrayView1d< real64 > const meanAlphaView = meanAlpha.toView();
  arrayView1d< real64 > const alphaView = alpha.toView();
  forAll< parallelDevicePolicy<> >( 3, [=] GEOS_HOST_DEVICE ( localIndex const k )
  {
    wrapper.compute( k, 10.0, 80.0, meanAlphaView[k], alphaView[k] );
  } );
  meanAlpha.move( hostMemorySpace, false );
  alpha.move( hostMemorySpace, false );

  for( localIndex k = 0; k < 3; ++k )
  {
    EXPECT_DOUBLE_EQ( meanAlpha[k], perElementAlpha[k] );
    EXPECT_DOUBLE_EQ( alpha[k], perElementAlpha[k] );
  }
}

TEST_F( SolidThermalExpansionTest, tabulated )
{
  addTable( "tecTable", tableKnots, tableValues );
  SolidThermalExpansion & model = parse( "<SolidThermalExpansion name=\"tec\" drainedLinearTECTableName=\"tecTable\"/>", "tec" );
  SolidThermalExpansionUpdates const wrapper = model.createKernelUpdates();

  checkTabulated< serialPolicy >( wrapper );
  checkTabulated< parallelDevicePolicy<> >( wrapper );

  // Jacobian inside a segment, across knots, outside the table
  checkJacobian( wrapper, 10.0, 35.0 );
  checkJacobian( wrapper, 5.0, 70.0 );
  checkJacobian( wrapper, 85.0, 10.0 );
  checkJacobian( wrapper, 50.0, 110.0 );
}

TEST_F( SolidThermalExpansionTest, tabulatedSingleKnot )
{
  addTable( "tecTable", { 20.0 }, { 3.0e-6 } );
  SolidThermalExpansion & model = parse( "<SolidThermalExpansion name=\"tec\" drainedLinearTECTableName=\"tecTable\"/>", "tec" );
  SolidThermalExpansionUpdates const wrapper = model.createKernelUpdates();

  real64 meanAlpha = 0.0, alpha = 0.0;
  wrapper.compute( 0, 0.0, 50.0, meanAlpha, alpha );
  EXPECT_LT( relativeError( meanAlpha, 3.0e-6 ), 1.0e-14 );
  EXPECT_LT( relativeError( alpha, 3.0e-6 ), 1.0e-14 );
}

TEST_F( SolidThermalExpansionTest, inputErrors )
{
  addTable( "tecTable", tableKnots, tableValues );
  addTable( "nearestTable", tableKnots, tableValues, TableFunction::InterpolationType::Nearest );

  // neither a coefficient nor a table
  expectInputError( [&]() { parse( "<SolidThermalExpansion name=\"e1\"/>", "e1" ); },
                    "Either referenceDrainedLinearTEC or drainedLinearTECTableName must be given" );
  // a table combined with the linear law
  expectInputError( [&]() { parse( "<SolidThermalExpansion name=\"e2\" drainedLinearTECTableName=\"tecTable\" referenceDrainedLinearTEC=\"1e-5\"/>", "e2" ); },
                    "drainedLinearTECTableName cannot be combined with" );
  // a slope without its reference temperature
  expectInputError( [&]() { parse( "<SolidThermalExpansion name=\"e3\" referenceDrainedLinearTEC=\"1e-5\" dDrainedLinearTEC_dTemperature=\"1e-8\"/>", "e3" ); },
                    "referenceTemperature must be given when dDrainedLinearTEC_dTemperature is non-zero" );
  // a table that does not exist
  expectInputError( [&]() { parse( "<SolidThermalExpansion name=\"e4\" drainedLinearTECTableName=\"missingTable\"/>", "e4" ); },
                    "The TableFunction 'missingTable' could not be found" );
  // a table without linear interpolation
  expectInputError( [&]() { parse( "<SolidThermalExpansion name=\"e5\" drainedLinearTECTableName=\"nearestTable\"/>", "e5" ); },
                    "must use linear interpolation" );
}

/// Input checks done by the coupled solid that refers to the model
TEST_F( SolidThermalExpansionTest, coupledSolid )
{
  string const subModels =
    "<BiotPorosity name=\"rockPorosity\" defaultReferencePorosity=\"0.1\" defaultGrainBulkModulus=\"1.0e11\"/>"
    "<ConstantPermeability name=\"rockPerm\" permeabilityComponents=\"{ 1.0e-15, 1.0e-15, 1.0e-15 }\"/>"
    "<SolidThermalExpansion name=\"rockTEC\" referenceDrainedLinearTEC=\"1.0e-5\"/>";
  string const solid =
    "<ElasticIsotropic name=\"rockSolid\" defaultDensity=\"2700\" defaultBulkModulus=\"1.0e10\" defaultShearModulus=\"1.0e10\"";
  string const porousSolid =
    "<PorousElasticIsotropic name=\"rock\" solidModelName=\"rockSolid\" porosityModelName=\"rockPorosity\" "
    "permeabilityModelName=\"rockPerm\"";

  // the model attached to a porous solid is the one used by its kernels
  {
    CoupledSolidBase const & rock =
      parseManager( subModels + solid + "/>" + porousSolid + " solidThermalExpansionModelName=\"rockTEC\"/>" ).getConstitutiveRelation< CoupledSolidBase >( "rock" );
    EXPECT_TRUE( rock.hasSolidThermalExpansionModel() );
    EXPECT_FALSE( rock.hasTemperatureDependentThermalExpansion() );
    EXPECT_TRUE( rock.createThermalExpansionUpdates().getStrategy() == Strategy::Linear );
  }

  // without the attribute, the coefficient of the solid model is used, as before
  {
    CoupledSolidBase const & rock =
      parseManager( subModels + solid + "/>" + porousSolid + "/>" ).getConstitutiveRelation< CoupledSolidBase >( "rock" );
    EXPECT_FALSE( rock.hasSolidThermalExpansionModel() );
    EXPECT_TRUE( rock.createThermalExpansionUpdates().getStrategy() == Strategy::FromSolidBase );
  }

  // a model that does not exist
  expectInputError( [&]() { parseManager( subModels + solid + "/>" + porousSolid + " solidThermalExpansionModelName=\"missingTEC\"/>" ); },
                    "The SolidThermalExpansion model 'missingTEC' could not be found" );

  // the coefficient given both on the solid model and by the attached model
  expectInputError( [&]() { parseManager( subModels + solid + " defaultDrainedLinearTEC=\"1.0e-5\"/>" + porousSolid + " solidThermalExpansionModelName=\"rockTEC\"/>" ); },
                    "is given both by defaultDrainedLinearTEC on the solid model" );

  // the porosity coefficient given together with the model (the porosity then uses the model)
  string const subModelsWithPorosityCoefficient =
    "<BiotPorosity name=\"rockPorosity\" defaultReferencePorosity=\"0.1\" defaultGrainBulkModulus=\"1.0e11\" defaultPorosityTEC=\"1.0e-5\"/>"
    "<ConstantPermeability name=\"rockPerm\" permeabilityComponents=\"{ 1.0e-15, 1.0e-15, 1.0e-15 }\"/>"
    "<SolidThermalExpansion name=\"rockTEC\" referenceDrainedLinearTEC=\"1.0e-5\"/>";
  expectInputError( [&]() {
    parseManager( subModelsWithPorosityCoefficient + solid + "/>" + porousSolid +
                  " solidThermalExpansionModelName=\"rockTEC\"/>" );
  },
                    "defaultPorosityTEC cannot be given on the porosity model 'rockPorosity'" );

  // a coupled solid that does not support the model
  expectInputError( [&]() {
    parseManager( subModels +
                  "<NullModel name=\"nullSolid\"/>"
                  "<PressurePorosity name=\"pressurePorosity\" defaultReferencePorosity=\"0.1\" referencePressure=\"0.0\" compressibility=\"1.0e-9\"/>"
                  "<CompressibleSolidConstantPermeability name=\"rock\" solidModelName=\"nullSolid\" porosityModelName=\"pressurePorosity\" "
                  "permeabilityModelName=\"rockPerm\" solidThermalExpansionModelName=\"rockTEC\"/>" );
  },
                    "solidThermalExpansionModelName is not supported by CompressibleSolidConstantPermeability" );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  geos::setupEnvironment( argc, argv );
  int const result = RUN_ALL_TESTS();
  geos::cleanupEnvironment();
  return result;
}
