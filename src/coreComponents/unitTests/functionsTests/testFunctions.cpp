/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "codingUtilities/UnitTestUtilities.hpp"
#include "gtest/gtest.h"
#include "mainInterface/initialization.hpp"
#include "functions/FunctionManager.hpp"
#include "functions/FunctionBase.hpp"
#include "functions/TableFunction.hpp"
#include "functions/MultivariableTableFunction.hpp"
#include "functions/MultivariableTableFunctionKernels.hpp"
#include "mainInterface/GeosxState.hpp"

#ifdef GEOSX_USE_MATHPRESSO
  #include "functions/SymbolicFunction.hpp"
#endif

#include <random>

using namespace geosx;


void evaluate1DFunction( FunctionBase & function,
                         arrayView1d< real64 const > const & inputs,
                         arrayView1d< real64 const > const & outputs )
{
  for( localIndex ii=0; ii<inputs.size(); ++ii )
  {
    real64 input = inputs[ii];
    real64 predicted = function.evaluate( &input );
    real64 expected = outputs[ii];

    ASSERT_NEAR( predicted, expected, 1e-10 );
  }
}

void checkDirectionalDerivative( real64 const (&input)[4],
                                 real64 (& perturbedInput)[4],
                                 real64 const & val,
                                 real64 & perturbedVal,
                                 real64 const (&derivatives)[4],
                                 real64 (& perturbedDerivatives)[4],
                                 real64 const & perturb,
                                 real64 const & relTol,
                                 localIndex const direction,
                                 TableFunction::KernelWrapper kernelWrapper )
{
  LvArray::tensorOps::copy< 4 >( perturbedInput, input );
  real64 const dInput = perturb * ( input[direction] + perturb );
  perturbedInput[direction] += dInput;
  perturbedVal = kernelWrapper.compute( perturbedInput, perturbedDerivatives );

  geosx::testing::checkRelativeError( derivatives[direction], (perturbedVal-val)/dInput, relTol, geosx::testing::DEFAULT_ABS_TOL );
}

TEST( FunctionTests, 1DTable )
{
  FunctionManager * functionManager = &FunctionManager::getInstance();

  // 1D table, various interpolation methods
  localIndex const Naxis = 4;
  localIndex const Ntest = 6;

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

  TableFunction & table_a = dynamicCast< TableFunction & >( *functionManager->createChild( "TableFunction", "table_a" ) );
  table_a.setTableCoordinates( coordinates );
  table_a.setTableValues( values );
  table_a.reInitializeFunction();

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
  table_a.setInterpolationMethod( TableFunction::InterpolationType::Linear );
  table_a.reInitializeFunction();
  evaluate1DFunction( table_a, testCoordinates, testExpected );

  // Upper
  testExpected[0] = 1.0;
  testExpected[1] = 3.0;
  testExpected[2] = -5.0;
  testExpected[3] = -5.0;
  testExpected[4] = 7.0;
  testExpected[5] = 7.0;
  table_a.setInterpolationMethod( TableFunction::InterpolationType::Upper );
  table_a.reInitializeFunction();
  evaluate1DFunction( table_a, testCoordinates, testExpected );

  // Lower
  testExpected[0] = 1.0;
  testExpected[1] = 1.0;
  testExpected[2] = 3.0;
  testExpected[3] = 3.0;
  testExpected[4] = -5.0;
  testExpected[5] = 7.0;
  table_a.setInterpolationMethod( TableFunction::InterpolationType::Lower );
  table_a.reInitializeFunction();
  evaluate1DFunction( table_a, testCoordinates, testExpected );

  // Nearest
  testExpected[0] = 1.0;
  testExpected[1] = 1.0;
  testExpected[2] = 3.0;
  testExpected[3] = -5.0;
  testExpected[4] = 7.0;
  testExpected[5] = 7.0;
  table_a.setInterpolationMethod( TableFunction::InterpolationType::Nearest );
  table_a.reInitializeFunction();
  evaluate1DFunction( table_a, testCoordinates, testExpected );

}



TEST( FunctionTests, 2DTable )
{
  FunctionManager * functionManager = &FunctionManager::getInstance();

  // 2D table with linear interpolation
  // f(x, y) = 2*x - 3*y + 5
  localIndex const Ndim = 2;
  localIndex const Nx = 3;
  localIndex const Ny = 4;
  localIndex const Ntest = 20;
  string const inputName = "coordinates";

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
      real64 const x = coordinates[0][ii];
      real64 const y = coordinates[1][jj];
      values[tablePosition] = (2.0*x) - (3.0*y) + 5.0;
      ++tablePosition;
    }
  }

  // Set input names
  string_array inputVarNames( 1 );
  inputVarNames[0] = inputName;

  // Initialize the table
  TableFunction & table_b = dynamicCast< TableFunction & >( *functionManager->createChild( "TableFunction", "table_b" ) );
  table_b.setTableCoordinates( coordinates );
  table_b.setTableValues( values );
  table_b.setInterpolationMethod( TableFunction::InterpolationType::Linear );
  table_b.setInputVarNames( inputVarNames );
  table_b.reInitializeFunction();

  // Setup a group for testing the batch mode function evaluation
  conduit::Node node;
  dataRepository::Group testGroup( "testGroup", node );

  real64_array2d testCoordinates;
  testGroup.registerWrapper( inputName, &testCoordinates ).
    setSizedFromParent( 1 ).
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

    real64 const x = testCoordinates[ii][0];
    real64 const y = testCoordinates[ii][1];
    expected[ii] = (2.0*x) - (3.0*y) + 5.0;
    set.insert( ii );
  }

  // Evaluate the function in batch mode
  table_b.evaluate( testGroup, 0.0, set.toView(), output );

  // Compare results
  for( localIndex ii=0; ii<Ntest; ++ii )
  {
    ASSERT_NEAR( expected[ii], output[ii], 1e-10 );
  }
}


TEST( FunctionTests, 4DTable_multipleInputs )
{
  FunctionManager * functionManager = &FunctionManager::getInstance();

  // 4D table with linear interpolation
  // f(x, y, z, t) = 2.0 + 3*x - 5*y + 7*z + 11*t
  localIndex const Ndim = 4;
  localIndex const Nx = 3;
  localIndex const Ny = 4;
  localIndex const Nz = 5;
  localIndex const Nt = 2;
  localIndex const Ntest = 20;
  localIndex const Ntimes = 5;
  string const coordinatesName = "coordinates";
  string const timeName = "time";

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
          real64 const x = coordinates[0][ii];
          real64 const y = coordinates[1][jj];
          real64 const z = coordinates[2][kk];
          real64 const t = coordinates[3][mm];
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
  TableFunction & table_c = dynamicCast< TableFunction & >( *functionManager->createChild( "TableFunction", "table_c" ) );
  table_c.setTableCoordinates( coordinates );
  table_c.setTableValues( values );
  table_c.setInterpolationMethod( TableFunction::InterpolationType::Linear );
  table_c.setInputVarNames( inputVarNames );
  table_c.reInitializeFunction();

  // Setup a group for testing the batch mode function evaluation
  conduit::Node node;
  dataRepository::Group testGroup( "testGroup", node );

  real64_array2d testCoordinates;
  testGroup.registerWrapper( coordinatesName, &testCoordinates ).
    setSizedFromParent( 1 ).
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
    real64 const t = distribution( generator );

    for( localIndex jj=0; jj<Ntest; ++jj )
    {
      for( localIndex kk=0; kk<Ndim-1; ++kk )
      {
        testCoordinates[jj][kk] = distribution( generator );
      }

      real64 const x = testCoordinates[jj][0];
      real64 const y = testCoordinates[jj][1];
      real64 const z = testCoordinates[jj][2];
      expected[jj] = 2.0 + (3*x) - (5*y) + (7*z) + (11*t);
    }

    // Evaluate the function in batch mode
    table_c.evaluate( testGroup, t, set.toView(), output );

    // Compare results
    for( localIndex jj=0; jj<Ntest; ++jj )
    {
      ASSERT_NEAR( expected[jj], output[jj], 1e-10 );
    }
  }
}

TEST( FunctionTests, 4DTable_derivatives )
{
  FunctionManager * functionManager = &FunctionManager::getInstance();

  // 4D table with linear interpolation
  // f(x, y, z, t) = 2.0 + 3*x - 5*y*y + 7*z*z*z + 11*t*t*t*t
  localIndex const Ndim = 4;
  localIndex const Nx = 3;
  localIndex const Ny = 4;
  localIndex const Nz = 5;
  localIndex const Nt = 4;
  string const coordinatesName = "coordinates";
  string const timeName = "time";

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
  coordinates[3][1] = 0.34;
  coordinates[3][2] = 0.5;
  coordinates[3][3] = 1.0;

  real64_array values( Nx * Ny * Nz * Nt );
  for( localIndex mm=0, tablePosition=0; mm<Nt; ++mm )
  {
    for( localIndex kk=0; kk<Nz; ++kk )
    {
      for( localIndex jj=0; jj<Ny; ++jj )
      {
        for( localIndex ii=0; ii<Nx; ++ii, ++tablePosition )
        {
          real64 const x = coordinates[0][ii];
          real64 const y = coordinates[1][jj];
          real64 const z = coordinates[2][kk];
          real64 const t = coordinates[3][mm];

          values[tablePosition] = 2.0 + (3*x) - (5*y*y) + (7*z*z*z) + (11*t*t*t*t);
        }
      }
    }
  }

  // Set variable names
  string_array inputVarNames( 2 );
  inputVarNames[0] = coordinatesName;
  inputVarNames[1] = timeName;

  // Initialize the table
  TableFunction & table_d = dynamicCast< TableFunction & >( *functionManager->createChild( "TableFunction", "table_d" ) );
  table_d.setTableCoordinates( coordinates );
  table_d.setTableValues( values );
  table_d.setInterpolationMethod( TableFunction::InterpolationType::Linear );
  table_d.setInputVarNames( inputVarNames );
  table_d.reInitializeFunction();

  real64 const relTol = 1e-4;
  real64 const eps = std::numeric_limits< real64 >::epsilon();
  real64 const perturb = std::sqrt( eps );

  real64 const start = coordinates[0][0]-0.1; // start outside the table, to check derivatives there too
  real64 const end = coordinates[0][Nx-1]+0.1; // end outside the table
  localIndex const nSamples = 7; // try not to fall on table coordinate, otherwise the finite-difference approximation won't work
  real64 const delta = (end-start)/nSamples;

  real64 val = 0.0;
  real64 perturbedVal = 0.0;
  real64 input[4] = { start, start, start, start };
  real64 perturbedInput[4]{};
  real64 derivatives[4]{};
  real64 perturbedDerivatives[4]{};

  TableFunction::KernelWrapper kernelWrapper = table_d.createKernelWrapper();
  for( localIndex mm=0; mm<nSamples; ++mm, input[3] += delta )
  {
    input[2] = start;
    for( localIndex kk=0; kk<nSamples; ++kk, input[2] += delta )
    {
      input[1] = start;
      for( localIndex jj=0; jj<nSamples; ++jj, input[1] += delta )
      {
        input[0] = start;
        for( localIndex ii=0; ii<nSamples; ++ii, input[0] += delta )
        {
          // evaluate once to get the analytical derivatives
          val = kernelWrapper.compute( input, derivatives );

          // check derivatives

          for( localIndex direction = 0; direction < Ndim; ++direction )
          {
            checkDirectionalDerivative( input, perturbedInput,
                                        val, perturbedVal,
                                        derivatives, perturbedDerivatives,
                                        perturb, relTol, direction, kernelWrapper );
          }
        }
      }
    }
  }
}

#ifdef GEOSX_USE_MATHPRESSO

TEST( FunctionTests, 4DTable_symbolic )
{
  FunctionManager * functionManager = &FunctionManager::getInstance();

  // Symbolic function with four inputs
  string const expression = "1.0+(2.0*a)-(3.0*b*b)+(5.0*c*c*c)-(7.0*d*d*d*d)";
  localIndex const Ntest = 20;

  // Set variable names
  string_array inputVarNames( 4 );
  string const nameA = "a";
  string const nameB = "b";
  string const nameC = "c";
  string const nameD = "d";
  inputVarNames[0] = nameA;
  inputVarNames[1] = nameB;
  inputVarNames[2] = nameC;
  inputVarNames[3] = nameD;

  // Initialize the table
  SymbolicFunction & table_e = dynamicCast< SymbolicFunction & >( *functionManager->createChild( "SymbolicFunction", "table_e" ) );
  table_e.setSymbolicExpression( expression );
  table_e.setInputVarNames( inputVarNames );
  table_e.setSymbolicVariableNames( inputVarNames );
  table_e.initializeFunction();

  // Setup a group for testing the batch mode function evaluation
  conduit::Node node;
  dataRepository::Group testGroup( "testGroup", node );

  real64_array inputA;
  real64_array inputB;
  real64_array inputC;
  real64_array inputD;
  testGroup.registerWrapper( nameA, &inputA ).setSizedFromParent( 1 );
  testGroup.registerWrapper( nameB, &inputB ).setSizedFromParent( 1 );
  testGroup.registerWrapper( nameC, &inputC ).setSizedFromParent( 1 );
  testGroup.registerWrapper( nameD, &inputD ).setSizedFromParent( 1 );
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
    real64 const a = distribution( generator );
    real64 const b = distribution( generator );
    real64 const c = distribution( generator );
    real64 const d = distribution( generator );
    inputA[ii] = a;
    inputB[ii] = b;
    inputC[ii] = c;
    inputD[ii] = d;

    expected[ii] = 1.0+(2.0*a)-(3.0*b*b)+(5.0*c*c*c)-(7.0*d*d*d*d);
  }

  // Evaluate the function in batch mode
  table_e.evaluate( testGroup, 0.0, set.toView(), output );

  // Compare results
  for( localIndex jj=0; jj<Ntest; ++jj )
  {
    ASSERT_NEAR( expected[jj], output[jj], 1e-10 );
  }
}

#endif

template< integer NUM_DIMS, integer NUM_OPS >
void testMutivariableFunction( MultivariableTableFunction & function,
                               arrayView1d< real64 const > const & inputs,
                               arrayView1d< real64 const > const & expectedValues,
                               arrayView1d< real64 const > const & expectedDerivatives )
{
  localIndex const numElems = inputs.size() / NUM_DIMS;

  ASSERT_EQ( numElems * NUM_DIMS, inputs.size());
  ASSERT_EQ( numElems * NUM_OPS, expectedValues.size());

  real64_array evaluatedValues( numElems * NUM_OPS );
  real64_array evaluatedDerivatives( numElems * NUM_OPS * NUM_OPS );
  arrayView1d< real64 > evaluatedValuesView = evaluatedValues.toView(),
                        evaluatedDerivativesView = evaluatedDerivatives.toView();


  MultivariableTableFunctionStaticKernel< NUM_DIMS, NUM_OPS > kernel( function.getAxisMinimums(),
                                                                      function.getAxisMaximums(),
                                                                      function.getAxisPoints(),
                                                                      function.getAxisSteps(),
                                                                      function.getAxisStepInvs(),
                                                                      function.getAxisHypercubeMults(),
                                                                      function.getHypercubeData(),
                                                                      inputs,
                                                                      evaluatedValuesView,
                                                                      evaluatedDerivativesView );
  // Loop over cells on the device.
  forAll< parallelDevicePolicy< > >( numElems, [=] GEOSX_HOST_DEVICE
                                       ( localIndex const elemIndex )
  {
    kernel.compute( elemIndex );
  } );

  // Perform checks.
  forAll< serialPolicy >( numElems, [=] ( localIndex const elemIndex )
  {
    ASSERT_NEAR( expectedValues[elemIndex], evaluatedValuesView[elemIndex], 1e-3 );
  } );

  // Perform checks.
  forAll< serialPolicy >( numElems * NUM_DIMS, [=] ( localIndex const elemDerivIndex )
  {
    ASSERT_NEAR( expectedDerivatives[elemDerivIndex], evaluatedDerivativesView[elemDerivIndex], 1e-2 );
  } );
}

TEST( FunctionTests, 1DMultivariableTable )
{
  FunctionManager * functionManager = &FunctionManager::getInstance();

  // 1D table
  localIndex constexpr nDims = 1;
  localIndex constexpr nOps = 1;
  localIndex const nTest = 7;

  // Setup table
  real64_array axisMins;
  real64_array axisMaxs;
  integer_array axisPoints;
  axisMins.resize( nDims );
  axisMaxs.resize( nDims );
  axisPoints.resize( nDims );
  axisMins[0] = 1;
  axisMaxs[0] = 100;
  axisPoints[0] = 100;

  real64_array values( axisPoints[0] );
  // set this up so that values are equal to coordinate: f(x) = x
  for( auto i = 0; i < axisPoints[0]; i++ )
    values[i] = i + 1;

  MultivariableTableFunction & table_f = dynamicCast< MultivariableTableFunction & >( *functionManager->createChild( "MultivariableTableFunction", "table_f" ) );
  table_f.setTableCoordinates( nDims, nOps, axisMins, axisMaxs, axisPoints );
  table_f.setTableValues( values );
  table_f.initializeFunction();


  // Setup testing coordinates, expected values
  real64_array testCoordinates( nTest );
  testCoordinates[0] = -1.1;
  testCoordinates[1] = 0.0;
  testCoordinates[2] = 1.0;
  testCoordinates[3] = 2.8;
  testCoordinates[4] = 99.3;
  testCoordinates[5] = 100.0;
  testCoordinates[6] = 11234.534;

  real64_array testExpectedValues( testCoordinates );
  real64_array testExpectedDerivatives( nTest * nDims * nOps );
  for( auto i = 0; i < nTest * nDims * nOps; i++ )
    testExpectedDerivatives[i] = 1;


  testMutivariableFunction< nDims, nOps >( table_f, testCoordinates, testExpectedValues, testExpectedDerivatives );
}

real64 operator1 ( real64 const x, real64 const y ) { return 2 * x + 3 * y * y; }
real64 dOperator1_dx ( real64 const x, real64 const y ) { return 2; }
real64 dOperator1_dy ( real64 const x, real64 const y ) { return 6 * y; }
real64 operator2 ( real64 const x, real64 const y ) { return 2 * x * x + 3 / (y + 1); }
real64 dOperator2_dx ( real64 const x, real64 const y ) { return 4 * x; }
real64 dOperator2_dy ( real64 const x, real64 const y ) { return -3 / ((y + 1) * (y + 1)); }
real64 operator3 ( real64 const x, real64 const y ) { return 2 * x + 3; }
real64 dOperator3_dx ( real64 const x, real64 const y ) { return 2; }
real64 dOperator3_dy ( real64 const x, real64 const y ) { return 0; }

TEST( FunctionTests, 2DMultivariableTable )
{
  FunctionManager * functionManager = &FunctionManager::getInstance();


  // 1D table
  localIndex constexpr nDims = 2;
  localIndex constexpr nOps = 3;
  localIndex const nTest = 3;

  // Setup table
  real64_array axisMins( nDims );
  real64_array axisMaxs( nDims );
  integer_array axisPoints( nDims );


  axisMins[0] = 1;
  axisMins[1] = 0;
  axisMaxs[0] = 2;
  axisMaxs[1] = 1;
  axisPoints[0] = 1000;
  axisPoints[1] = 1100;

  real64_array axisSteps( nDims );
  for( auto i = 0; i < nDims; i++ )
    axisSteps[i] = (axisMaxs[i] - axisMins[i]) /  axisPoints[i];

  real64_array values( axisPoints[0] * axisPoints[1] * nOps );
  for( auto i = 0; i < axisPoints[0]; i++ )
    for( auto j = 0; j < axisPoints[1]; j++ )
    {
      auto x = axisMins[0] + i * axisSteps[0];
      auto y = axisMins[1] + j * axisSteps[1];

      values[( i * axisPoints[1] + j ) * nOps] = operator1( x, y );
      values[( i * axisPoints[1] + j ) * nOps + 1] = operator2( x, y );
      values[( i * axisPoints[1] + j ) * nOps + 2] = operator3( x, y );
    }

  MultivariableTableFunction & table_g = dynamicCast< MultivariableTableFunction & >( *functionManager->createChild( "MultivariableTableFunction", "table_f" ) );
  table_g.setTableCoordinates( nDims, nOps, axisMins, axisMaxs, axisPoints );
  table_g.setTableValues( values );
  table_g.initializeFunction();


  // Setup testing coordinates, expected values
  real64_array testCoordinates( nTest * nDims );
  testCoordinates[0] = 1.2334;
  testCoordinates[1] = 0.1232;
  testCoordinates[2] = 1.7342;
  testCoordinates[3] = 0.2454;
  testCoordinates[4] = 2.0;
  testCoordinates[5] = 0.7745;

  real64_array testExpectedValues( nTest * nOps );
  real64_array testExpectedDerivatives( nTest * nOps * nDims );
  for( auto i = 0; i < nTest; i++ )
  {
    testExpectedValues[i * nOps] = operator1( testCoordinates[i * nDims], testCoordinates[i * nDims + 1] );
    testExpectedValues[i * nOps + 1] = operator2( testCoordinates[i * nDims], testCoordinates[i * nDims + 1] );
    testExpectedValues[i * nOps + 2] = operator3( testCoordinates[i * nDims], testCoordinates[i * nDims + 1] );
    testExpectedDerivatives[i * nOps * nDims] = dOperator1_dx( testCoordinates[i * nDims], testCoordinates[i * nDims + 1] );
    testExpectedDerivatives[i * nOps * nDims + 1] = dOperator1_dy( testCoordinates[i * nDims], testCoordinates[i * nDims + 1] );
    testExpectedDerivatives[i * nOps * nDims + 2] = dOperator2_dx( testCoordinates[i * nDims], testCoordinates[i * nDims + 1] );
    testExpectedDerivatives[i * nOps * nDims + 3] = dOperator2_dy( testCoordinates[i * nDims], testCoordinates[i * nDims + 1] );
    testExpectedDerivatives[i * nOps * nDims + 4] = dOperator3_dx( testCoordinates[i * nDims], testCoordinates[i * nDims + 1] );
    testExpectedDerivatives[i * nOps * nDims + 5] = dOperator3_dy( testCoordinates[i * nDims], testCoordinates[i * nDims + 1] );
  }


  testMutivariableFunction< nDims, nOps >( table_g, testCoordinates, testExpectedValues, testExpectedDerivatives );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geosx::GeosxState state( geosx::basicSetup( argc, argv ) );

  int const result = RUN_ALL_TESTS();

  geosx::basicCleanup();

  return result;
}
