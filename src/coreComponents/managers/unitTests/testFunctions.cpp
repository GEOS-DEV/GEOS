/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "gtest/gtest.h"
#include "managers/initialization.hpp"
#include "managers/Functions/FunctionManager.hpp"
#include "managers/Functions/FunctionBase.hpp"
#include "managers/Functions/TableFunction.hpp"
#ifdef GEOSX_USE_MATHPRESSO
#include "managers/Functions/SymbolicFunction.hpp"
#endif

#include <random>

using namespace geosx;



void evaluate1DFunction( FunctionBase * function,
                         arrayView1d< real64 const > const & inputs,
                         arrayView1d< real64 const > const & outputs )
{
  for( localIndex ii=0; ii<inputs.size(); ++ii )
  {
    real64 input = inputs[ii];
    real64 predicted = function->evaluate( &input );
    real64 expected = outputs[ii];

    ASSERT_NEAR( predicted, expected, 1e-10 );
  }
}



TEST( FunctionTests, 1DTable )
{
  FunctionManager * functionManager = &FunctionManager::FunctionManager::instance();

  // 1D table, various interpolation methods
  localIndex Naxis = 4;
  localIndex Ntest = 6;

  // Setup table
  array1d< real64_array > coordinates;
  coordinates.resize( 1 );
  coordinates[0].resize( Naxis );
  coordinates[0][0] = -1.0;
  coordinates[0][1] = 0.0;
  coordinates[0][2] = 2.0;
  coordinates[0][3] = 5.0;

  real64_array values( Naxis );
  values[0] = 1.0;
  values[1] = 3.0;
  values[2] = -5.0;
  values[3] = 7.0;

  TableFunction * table_a = functionManager->createChild( "TableFunction", "table_a" )->groupCast< TableFunction * >();
  table_a->setTableCoordinates( coordinates );
  table_a->setTableValues( values );
  table_a->reInitializeFunction();

  // Setup testing coordinates, expected values
  real64_array testCoordinates( Ntest );
  testCoordinates[0] = -1.1;
  testCoordinates[1] = -0.5;
  testCoordinates[2] = 0.2;
  testCoordinates[3] = 2.0;
  testCoordinates[4] = 4.0;
  testCoordinates[5] = 10.0;

  // Linear Interpolation
  real64_array testExpected( Ntest );
  testExpected[0] = 1.0;
  testExpected[1] = 2.0;
  testExpected[2] = 2.2;
  testExpected[3] = -5.0;
  testExpected[4] = 3.0;
  testExpected[5] = 7.0;
  table_a->setInterpolationMethod( TableFunction::InterpolationType::Linear );
  evaluate1DFunction( table_a, testCoordinates, testExpected );

  // Upper
  testExpected[0] = 1.0;
  testExpected[1] = 3.0;
  testExpected[2] = -5.0;
  testExpected[3] = -5.0;
  testExpected[4] = 7.0;
  testExpected[5] = 7.0;
  table_a->setInterpolationMethod( TableFunction::InterpolationType::Upper );
  evaluate1DFunction( table_a, testCoordinates, testExpected );

  // Lower
  testExpected[0] = 1.0;
  testExpected[1] = 1.0;
  testExpected[2] = 3.0;
  testExpected[3] = 3.0;
  testExpected[4] = -5.0;
  testExpected[5] = 7.0;
  table_a->setInterpolationMethod( TableFunction::InterpolationType::Lower );
  evaluate1DFunction( table_a, testCoordinates, testExpected );

  // Nearest
  testExpected[0] = 1.0;
  testExpected[1] = 1.0;
  testExpected[2] = 3.0;
  testExpected[3] = -5.0;
  testExpected[4] = 7.0;
  testExpected[5] = 7.0;
  table_a->setInterpolationMethod( TableFunction::InterpolationType::Nearest );
  evaluate1DFunction( table_a, testCoordinates, testExpected );

}



TEST( FunctionTests, 2DTable )
{
  FunctionManager * functionManager = &FunctionManager::FunctionManager::instance();

  // 2D table with linear interpolation
  // f(x, y) = 2*x - 3*y + 5
  localIndex Ndim = 2;
  localIndex Nx = 3;
  localIndex Ny = 4;
  localIndex Ntest = 20;
  string inputName = "coordinates";

  // Setup table
  array1d< real64_array > coordinates;
  coordinates.resize( Ndim );
  coordinates[0].resize( Nx );
  coordinates[0][0] = -1.0;
  coordinates[0][1] = 0.0;
  coordinates[0][2] = 2.0;
  coordinates[1].resize( Ny );
  coordinates[1][0] = -1.0;
  coordinates[1][1] = 0.0;
  coordinates[1][2] = 1.0;
  coordinates[1][3] = 2.0;

  real64_array values( Nx * Ny );
  localIndex tablePosition = 0;
  for( localIndex jj=0; jj<Ny; ++jj )
  {
    for( localIndex ii=0; ii<Nx; ++ii )
    {
      real64 x = coordinates[0][ii];
      real64 y = coordinates[1][jj];
      values[tablePosition] = (2.0*x) - (3.0*y) + 5.0;
      ++tablePosition;
    }
  }

  // Set input names
  string_array inputVarNames( 1 );
  inputVarNames[0] = inputName;

  // Initialize the table
  TableFunction * table_b = functionManager->createChild( "TableFunction", "table_b" )->groupCast< TableFunction * >();
  table_b->setTableCoordinates( coordinates );
  table_b->setTableValues( values );
  table_b->setInterpolationMethod( TableFunction::InterpolationType::Linear );
  table_b->setInputVarNames( inputVarNames );
  table_b->reInitializeFunction();

  // Setup a group for testing the batch mode function evaluation
  string groupName = "testGroup";
  dataRepository::Group testGroup( groupName, nullptr );

  real64_array2d testCoordinates;
  testGroup.registerWrapper( inputName, &testCoordinates )->
    setSizedFromParent( 1 )->
    reference().resizeDimension< 1 >( Ndim );
  testGroup.resize( Ntest );

  // Build testing inputs/outputs
  real64_array expected( Ntest );
  real64_array output( Ntest );
  SortedArray< localIndex > set;

  // Build the set
  for( localIndex ii=0; ii<Ntest; ++ii )
  {
    set.insert( ii );
  }

  // Setup the random number generator
  std::default_random_engine generator;
  std::uniform_real_distribution< double > distribution( -0.99, 1.99 );

  // Test the function
  for( localIndex ii=0; ii<Ntest; ++ii )
  {
    for( localIndex jj=0; jj<Ndim; ++jj )
    {
      testCoordinates[ii][jj] = distribution( generator );
    }

    real64 x = testCoordinates[ii][0];
    real64 y = testCoordinates[ii][1];
    expected[ii] = (2.0*x) - (3.0*y) + 5.0;
    set.insert( ii );
  }

  // Evaluate the function in batch mode
  table_b->evaluate( &(testGroup), 0.0, set.toView(), output );

  // Compare results
  for( localIndex ii=0; ii<Ntest; ++ii )
  {
    ASSERT_NEAR( expected[ii], output[ii], 1e-10 );
  }
}


TEST( FunctionTests, 4DTable_multipleInputs )
{
  FunctionManager * functionManager = &FunctionManager::FunctionManager::instance();

  // 3D table with linear interpolation
  // f(x, y, z, t) = 2.0 + 3*x - 5*y + 7*z + 11*t
  localIndex Ndim = 4;
  localIndex Nx = 3;
  localIndex Ny = 4;
  localIndex Nz = 5;
  localIndex Nt = 2;
  localIndex Ntest = 20;
  localIndex Ntimes = 5;
  string coordinatesName = "coordinates";
  string timeName = "time";

  // Setup table
  array1d< real64_array > coordinates;
  coordinates.resize( Ndim );
  coordinates[0].resize( Nx );
  coordinates[0][0] = -1.0;
  coordinates[0][1] = 0.0;
  coordinates[0][2] = 1.0;
  coordinates[1].resize( Ny );
  coordinates[1][0] = -1.0;
  coordinates[1][1] = 0.0;
  coordinates[1][2] = 0.5;
  coordinates[1][3] = 1.0;
  coordinates[2].resize( Nz );
  coordinates[2][0] = -1.0;
  coordinates[2][1] = -0.4;
  coordinates[2][2] = 0.3;
  coordinates[2][3] = 0.5;
  coordinates[2][4] = 1.0;
  coordinates[3].resize( Nt );
  coordinates[3][0] = -1.0;
  coordinates[3][1] = 1.0;

  real64_array values( Nx * Ny * Nz * Nt );
  localIndex tablePosition = 0;
  for( localIndex mm=0; mm<Nt; ++mm )
  {
    for( localIndex kk=0; kk<Nz; ++kk )
    {
      for( localIndex jj=0; jj<Ny; ++jj )
      {
        for( localIndex ii=0; ii<Nx; ++ii )
        {
          real64 x = coordinates[0][ii];
          real64 y = coordinates[1][jj];
          real64 z = coordinates[2][kk];
          real64 t = coordinates[3][mm];
          values[tablePosition] = 2.0 + (3*x) - (5*y) + (7*z) + (11*t);
          ++tablePosition;
        }
      }
    }
  }

  // Set variable names
  string_array inputVarNames( 2 );
  inputVarNames[0] = coordinatesName;
  inputVarNames[1] = timeName;

  // Initialize the table
  TableFunction * table_c = functionManager->createChild( "TableFunction", "table_c" )->groupCast< TableFunction * >();
  table_c->setTableCoordinates( coordinates );
  table_c->setTableValues( values );
  table_c->setInterpolationMethod( TableFunction::InterpolationType::Linear );
  table_c->setInputVarNames( inputVarNames );
  table_c->reInitializeFunction();

  // Setup a group for testing the batch mode function evaluation
  string groupName = "testGroup";
  dataRepository::Group testGroup( groupName, nullptr );

  real64_array2d testCoordinates;
  testGroup.registerWrapper( coordinatesName, &testCoordinates )->
    setSizedFromParent( 1 )->
    reference().resizeDimension< 1 >( Ndim - 1 );
  testGroup.resize( Ntest );

  // Build testing inputs/outputs
  real64_array expected( Ntest );
  real64_array output( Ntest );
  SortedArray< localIndex > set;

  // Fill out the set
  for( localIndex ii=0; ii<Ntest; ++ii )
  {
    set.insert( ii );
  }

  // Setup a random number generator
  std::default_random_engine generator;
  std::uniform_real_distribution< double > distribution( -0.99, 0.99 );

  // Build the inputs
  for( localIndex ii=0; ii<Ntimes; ++ii )
  {
    real64 t = distribution( generator );

    for( localIndex jj=0; jj<Ntest; ++jj )
    {
      for( localIndex kk=0; kk<Ndim-1; ++kk )
      {
        testCoordinates[jj][kk] = distribution( generator );
      }

      real64 x = testCoordinates[jj][0];
      real64 y = testCoordinates[jj][1];
      real64 z = testCoordinates[jj][2];
      expected[jj] = 2.0 + (3*x) - (5*y) + (7*z) + (11*t);
    }

    // Evaluate the function in batch mode
    table_c->evaluate( &(testGroup), t, set.toView(), output );

    // Compare results
    for( localIndex jj=0; jj<Ntest; ++jj )
    {
      ASSERT_NEAR( expected[jj], output[jj], 1e-10 );
    }
  }
}



#ifdef GEOSX_USE_MATHPRESSO

TEST( FunctionTests, 4DTable_symbolic )
{
  FunctionManager * functionManager = &FunctionManager::FunctionManager::instance();

  // Symbolic function with four inputs
  string expression = "1.0+(2.0*a)-(3.0*b*b)+(5.0*c*c*c)-(7.0*d*d*d*d)";
  localIndex Ntest = 20;

  // Set variable names
  string_array inputVarNames( 4 );
  string nameA = "a";
  string nameB = "b";
  string nameC = "c";
  string nameD = "d";
  inputVarNames[0] = nameA;
  inputVarNames[1] = nameB;
  inputVarNames[2] = nameC;
  inputVarNames[3] = nameD;

  // Initialize the table
  SymbolicFunction * table_d = functionManager->createChild( "SymbolicFunction", "table_d" )->groupCast< SymbolicFunction * >();
  table_d->setSymbolicExpression( expression );
  table_d->setInputVarNames( inputVarNames );
  table_d->setSymbolicVariableNames( inputVarNames );
  table_d->initializeFunction();

  // Setup a group for testing the batch mode function evaluation
  string groupName = "testGroup";
  dataRepository::Group testGroup( groupName, nullptr );
  real64_array inputA;
  real64_array inputB;
  real64_array inputC;
  real64_array inputD;
  testGroup.registerWrapper( nameA, &inputA )->setSizedFromParent( 1 );
  testGroup.registerWrapper( nameB, &inputB )->setSizedFromParent( 1 );
  testGroup.registerWrapper( nameC, &inputC )->setSizedFromParent( 1 );
  testGroup.registerWrapper( nameD, &inputD )->setSizedFromParent( 1 );
  testGroup.resize( Ntest );

  // Build testing inputs/outputs
  real64_array expected( Ntest );
  real64_array output( Ntest );
  SortedArray< localIndex > set;

  // Fill out the set
  for( localIndex ii=0; ii<Ntest; ++ii )
  {
    set.insert( ii );
  }

  // Setup a random number generator
  std::default_random_engine generator;
  std::uniform_real_distribution< double > distribution( -1.0, 1.0 );

  // Build the inputs
  for( localIndex ii=0; ii<Ntest; ++ii )
  {
    real64 a = distribution( generator );
    real64 b = distribution( generator );
    real64 c = distribution( generator );
    real64 d = distribution( generator );
    inputA[ii] = a;
    inputB[ii] = b;
    inputC[ii] = c;
    inputD[ii] = d;

    expected[ii] = 1.0+(2.0*a)-(3.0*b*b)+(5.0*c*c*c)-(7.0*d*d*d*d);
  }

  // Evaluate the function in batch mode
  table_d->evaluate( &(testGroup), 0.0, set.toView(), output );

  // Compare results
  for( localIndex jj=0; jj<Ntest; ++jj )
  {
    ASSERT_NEAR( expected[jj], output[jj], 1e-10 );
  }
}

#endif


int main( int argc, char * * argv )
{
  basicSetup( argc, argv );

  ::testing::InitGoogleTest( &argc, argv );

  int const result = RUN_ALL_TESTS();

  basicCleanup();

  return result;
}
