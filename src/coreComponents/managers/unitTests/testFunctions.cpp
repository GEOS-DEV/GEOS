/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "gtest/gtest.h"
#include "managers/initialization.hpp"
#include "managers/Functions/FunctionManager.hpp"
#include "managers/Functions/FunctionBase.hpp"
#include "managers/Functions/TableFunction.hpp"
#include "managers/Functions/SymbolicFunction.hpp"

using namespace geosx;



void evaluate1DFunction(FunctionBase * function,
                        real64_array inputs,
                        real64_array outputs)
{
  for (localIndex ii=0; ii<inputs.size(); ++ii)
  {
    real64 input = inputs[ii];
    real64 predicted = function->Evaluate(&input);
    real64 expected = outputs[ii];
    
    EXPECT_DOUBLE_EQ( predicted, expected );
  }
}



TEST( FunctionTests, 1DTable )
{
  FunctionManager * functionManager = &FunctionManager::FunctionManager::Instance();

  // 1D table, various interpolation methods
  localIndex Naxis = 4;
  localIndex Ntest = 6;

  // Setup table
  array1d< real64_array > coordinates;
  coordinates.resize(1);
  coordinates[0].resize(Naxis);
  coordinates[0][0] = -1.0;
  coordinates[0][1] = 0.0;
  coordinates[0][2] = 2.0;
  coordinates[0][3] = 5.0;

  real64_array values(Naxis);
  values[0] = 1.0;
  values[1] = 3.0;
  values[2] = -5.0;
  values[3] = 7.0;

  TableFunction * table_a = functionManager->CreateChild( "TableFunction", "table_a" )->group_cast< TableFunction * >();
  table_a->setTableCoordinates(coordinates);
  table_a->setTableValues(values);
  table_a->reInitializeFunction();

  // Setup testing coordinates, expected values
  real64_array testCoordinates(Ntest);
  testCoordinates[0] = -1.1;
  testCoordinates[1] = -0.5;
  testCoordinates[2] = 0.2;
  testCoordinates[3] = 2.0;
  testCoordinates[4] = 4.0;
  testCoordinates[5] = 10.0;

  // Linear Interpolation
  real64_array testExpected(Ntest);
  testExpected[0] = 1.0;
  testExpected[1] = 2.0;
  testExpected[2] = 2.2;
  testExpected[3] = -5.0;
  testExpected[4] = 3.0;
  testExpected[5] = 7.0;
  table_a->setInterpolationMethod("linear");
  evaluate1DFunction( table_a, testCoordinates, testExpected);

  // Upper
  testExpected[0] = 1.0;
  testExpected[1] = 3.0;
  testExpected[2] = -5.0;
  testExpected[3] = -5.0;
  testExpected[4] = 7.0;
  testExpected[5] = 7.0;
  table_a->setInterpolationMethod("upper");
  evaluate1DFunction( table_a, testCoordinates, testExpected);

  // Lower
  testExpected[0] = 1.0;
  testExpected[1] = 1.0;
  testExpected[2] = 3.0;
  testExpected[3] = 3.0;
  testExpected[4] = -5.0;
  testExpected[5] = 7.0;
  table_a->setInterpolationMethod("lower");
  evaluate1DFunction( table_a, testCoordinates, testExpected);

  // Nearest
  testExpected[0] = 1.0;
  testExpected[1] = 1.0;
  testExpected[2] = 3.0;
  testExpected[3] = -5.0;
  testExpected[4] = 7.0;
  testExpected[5] = 7.0;
  table_a->setInterpolationMethod("nearest");
  evaluate1DFunction( table_a, testCoordinates, testExpected);
   
}


// #ifdef GEOSX_USE_MATHPRESSO
// #endif


int main( int argc, char * * argv )
{
  basicSetup( argc, argv );

  ::testing::InitGoogleTest( &argc, argv );
  
  int const result = RUN_ALL_TESTS();

  basicCleanup();

  return result;
}
