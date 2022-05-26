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

#include <gtest/gtest.h>

#include "dataRepository/xmlWrapper.hpp"

const char IGNORE_OUTPUT[] = ".*";

using namespace geosx;

TEST( testXmlWrapper, array3d_errors )
{
  string input;

  // This should work
  {
    input = " { { {0,1,2},{3,4,5} }, { {6,7,8},{9,10,11} }, { {12,13,14},{15,16,17} } , { {18,19,20},{21,22,23} } }";
    array3d< int > array;
    xmlWrapper::stringToInputVariable( array, input );
  }
  // This should fail the num('{')==num('}') test
  {
    input = " { { {0,1,2},{3,4,5} }, { {6,7,8},{9,10,11} }, { {12,13,14},{15,16,17} } , { {18,19,20},{21,22,23} } ";
    array3d< int > array;
    EXPECT_DEATH_IF_SUPPORTED( xmlWrapper::stringToInputVariable( array, input ), IGNORE_OUTPUT );
  }
  {
    input = " { { {0,1,2},{3,4,5} }, { {6,7,8},{9,10,11} }, { {12,13,14},{15,16,17}  , { {18,19,20},{21,22,23} } }";
    array3d< int > array;
    EXPECT_DEATH_IF_SUPPORTED( xmlWrapper::stringToInputVariable( array, input ), IGNORE_OUTPUT );
  }
  {
    input = " { { {0,1,2},{3,4,5} }, { {6,7,8},{9,10,11} }, { {12,13,14,{15,16,17} } , { {18,19,20},{21,22,23} } }";
    array3d< int > array;
    EXPECT_DEATH_IF_SUPPORTED( xmlWrapper::stringToInputVariable( array, input ), IGNORE_OUTPUT );
  }
  {
    input = " { { {0,1,2},{3,4,5} }, { {6,7,8,{9,10,11} }, { {12,13,14},{15,16,17} } , { {18,19,20},{21,22,23} } }";
    array3d< int > array;
    EXPECT_DEATH_IF_SUPPORTED( xmlWrapper::stringToInputVariable( array, input ), IGNORE_OUTPUT );
  }
  {
    input = " { { 0,1,2},{3,4,5} }, { {6,7,8},{9,10,11} }, { {12,13,14},{15,16,17} } , { {18,19,20},{21,22,23} } }";
    array3d< int > array;
    EXPECT_DEATH_IF_SUPPORTED( xmlWrapper::stringToInputVariable( array, input ), IGNORE_OUTPUT );
  }
  {
    input = "  { {0,1,2},{3,4,5} }, { {6,7,8},{9,10,11} }, { {12,13,14},{15,16,17} } , { {18,19,20},{21,22,23} } ";
    array3d< int > array;
    EXPECT_DEATH_IF_SUPPORTED( xmlWrapper::stringToInputVariable( array, input ), IGNORE_OUTPUT );
  }

  {
    input = " { { {,1,2},{3,4,5} }, { {6,7,8},{9,10,11} }, { {12,13,14},{15,16,17} } , { {18,19,20},{21,22,23} } }";
    array3d< int > array;
    EXPECT_DEATH_IF_SUPPORTED( xmlWrapper::stringToInputVariable( array, input ), IGNORE_OUTPUT );
  }

  {
    input = " { { {},{3,4,5} }, { {6,7,8},{9,10,11} }, { {12,13,14},{15,16,17} } , { {18,19,20},{21,22,23} } }";
    array3d< int > array;
    EXPECT_DEATH_IF_SUPPORTED( xmlWrapper::stringToInputVariable( array, input ), IGNORE_OUTPUT );
  }
  {
    input = " { { {0,1,2}}{ } }";
    array3d< int > array;
    EXPECT_DEATH_IF_SUPPORTED( xmlWrapper::stringToInputVariable( array, input ), IGNORE_OUTPUT );
  }


}

TEST( testXmlWrapper, array3d )
{
//  string input = " { { {0,1,2},{3,4,5} }, { {6,7,8},{9,10,11} }, { {12,13,14},{15,16,17} } , { {18,19,20},{21,22,23} }
// }";
  string input;
  input += "{ ";
  localIndex numI = 4;
  localIndex numJ = 5;
  localIndex numK = 3;
  for( localIndex i=0; i<4; ++i )
  {
    input += "{ ";
    for( localIndex j=0; j<5; ++j )
    {
      input += "{ ";
      for( localIndex k=0; k<3; ++k )
      {
        input += std::to_string( i*2+j*3+k*4 );
        if( k<(numK-1) )
        {
          input += " , ";
        }
      }
      input += " }";
      if( j<(numJ-1) )
      {
        input += " , ";
      }
    }
    input += " }";
    if( i<(numI-1) )
    {
      input += " , ";
    }
  }
  input += " }";

  array3d< int > array;
  xmlWrapper::stringToInputVariable( array, input );

  ASSERT_EQ( array.size( 0 ), numI );
  ASSERT_EQ( array.size( 1 ), numJ );
  ASSERT_EQ( array.size( 2 ), numK );

  for( localIndex i=0; i<array.size( 0 ); ++i )
  {
    for( localIndex j=0; j<array.size( 1 ); ++j )
    {
      for( localIndex k=0; k<array.size( 2 ); ++k )
      {
        ASSERT_EQ( array[i][j][k], i*2+j*3+k*4 );
      }
    }
  }
}

template< typename T, int N >
void testTensorWellFormed( string const & input, Tensor< T, N > const expected )
{
  Tensor< T, N > output;
  xmlWrapper::stringToInputVariable( output, input );
  EXPECT_EQ( output, expected );
}

template< typename T, int N >
void testTensorIllFormed( string const & input )
{
  Tensor< T, N > output;
  EXPECT_THROW( xmlWrapper::stringToInputVariable( output, input ), InputError );
}

TEST( testXmlWrapper, TensorWellFormed )
{
  testTensorWellFormed< real64, 3 >( "{1.0,2.0,3.0}", { 1.0, 2.0, 3.0 } );
  testTensorWellFormed< real64, 3 >( "{ 1.0, 2.0, 3.0 }", { 1.0, 2.0, 3.0 } );
  testTensorWellFormed< real64, 3 >( "  {  1.0  , 2.0,3.0  }    ", { 1.0, 2.0, 3.0 } );
  testTensorWellFormed< real64, 6 >( "{1.0,2.0,3.0,4.0,5.0,6.0}", { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 } );
  testTensorWellFormed< real64, 6 >( "{ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 }", { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 } );
  testTensorWellFormed< real64, 6 >( "  { 1.0 ,  2.0 ,3.0 , 4.0,   5.0,6.0}  ", { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 } );
}

TEST( testXmlWrapper, TensorIllFormed )
{
  testTensorIllFormed< real64, 3 >( "1.0, 2.0, 3.0 }" ); // missing opening brance
  testTensorIllFormed< real64, 3 >( "{ 1.0, 2.0, 3.0" ); // missing closing brace
  testTensorIllFormed< real64, 3 >( "{1.0 2.0, 3.0 }" ); // missing a comma
  testTensorIllFormed< real64, 3 >( "{ 1.0, 2.0 }" ); // too few values
  testTensorIllFormed< real64, 3 >( "{ 1.0, 2.0, 3.0, 4.0 }" ); // too many values
  testTensorIllFormed< real64, 3 >( "{ 1.0, 2.0, 3.0 }  ,4.0" ); // extra characters
  testTensorIllFormed< real64, 3 >( "{ 1.0, O.1, 3.0 }" ); // invalid floating point value
  testTensorIllFormed< real64, 3 >( "{ 1.0, 2.O, 3.0 }" ); // invalid floating point value
}


// Attribute format tests
class real64AttributeTestFixture :public ::testing::TestWithParam<std::tuple<string, real64, bool>>
{
protected:
    string attributeString;
    real64 expectedValue;
    real64 parsedValue;
    bool failureFlag;

    void parseString()
    {
      GEOSX_LOG_RANK_0( "Parsing string: " << attributeString );
      xmlWrapper::stringToInputVariable< real64 >(parsedValue, attributeString);
      ASSERT_NEAR( expectedValue, parsedValue, 1e-10 );
    }
};

TEST_P(real64AttributeTestFixture, testParsing)
{
    auto testParams = GetParam();
    attributeString = std::get<0>(testParams);
    expectedValue = std::get<1>(testParams);
    failureFlag = std::get<2>(testParams);

    if (failureFlag)
    {
      EXPECT_THROW( parseString(), InputError );
    }
    else
    {
      parseString();
    }
}

INSTANTIATE_TEST_CASE_P(
        real64AttributeTests,
        real64AttributeTestFixture,
        ::testing::Values(std::make_tuple("1", 1,  false),
                          std::make_tuple("-23", -23,  false),
                          std::make_tuple("4.5", 4.5, false),
                          std::make_tuple("4.", 4.0, false),
                          std::make_tuple("6e1", 6e1, false),
                          std::make_tuple("7e-2", 7e-2, false),
                          std::make_tuple("8.765e0", 8.765, false),
                          std::make_tuple("-1.2e-3", -1.2e-3, false),
                          std::make_tuple("3.5E+4", 3.5e+4, false),
                          std::make_tuple("6.1E05", 6.1e05, false),
                          std::make_tuple("-7.1e06", -7.1e06, false),
                          std::make_tuple("8.1e+02", 8.1e+02, false),
                          std::make_tuple("9.1e-01", 9.1e-01, false),
                          std::make_tuple("alpha", 0, true),
                          std::make_tuple("1beta234", 0, true),
                          std::make_tuple("1.234gamma", 0, true),
                          std::make_tuple("1.2.3", 0, true),
                          std::make_tuple("1e2.3 ", 0, true),
                          std::make_tuple("1 ", 0, true),
                          std::make_tuple("1e", 0, true),
                          std::make_tuple("1e-", 0, true),
                          std::make_tuple("1e+", 0, true)));


class real32AttributeTestFixture :public ::testing::TestWithParam<std::tuple<string, real32, bool>>
{
protected:
    string attributeString;
    real32 expectedValue;
    real32 parsedValue;
    bool failureFlag;

    void parseString()
    {
      GEOSX_LOG_RANK_0( "Parsing string: " << attributeString );
      xmlWrapper::stringToInputVariable< real32 >(parsedValue, attributeString);
      ASSERT_NEAR( expectedValue, parsedValue, 1e-10 );
    }
};

TEST_P(real32AttributeTestFixture, testParsing)
{
    auto testParams = GetParam();
    attributeString = std::get<0>(testParams);
    expectedValue = std::get<1>(testParams);
    failureFlag = std::get<2>(testParams);

    if (failureFlag)
    {
      EXPECT_THROW( parseString(), InputError );
    }
    else
    {
      parseString();
    }
}

INSTANTIATE_TEST_CASE_P(
        real32AttributeTests,
        real32AttributeTestFixture,
        ::testing::Values(std::make_tuple("1", 1,  false),
                          std::make_tuple("-23", -23,  false),
                          std::make_tuple("4.5", 4.5, false),
                          std::make_tuple("4.", 4.0, false),
                          std::make_tuple("6e1", 6e1, false),
                          std::make_tuple("7e-2", 7e-2, false),
                          std::make_tuple("8.765e0", 8.765, false),
                          std::make_tuple("-1.2e-3", -1.2e-3, false),
                          std::make_tuple("3.5E+4", 3.5e+4, false),
                          std::make_tuple("6.1E05", 6.1e05, false),
                          std::make_tuple("-7.1e06", -7.1e06, false),
                          std::make_tuple("8.1e+02", 8.1e+02, false),
                          std::make_tuple("9.1e-01", 9.1e-01, false),
                          std::make_tuple("alpha", 0, true),
                          std::make_tuple("1beta234", 0, true),
                          std::make_tuple("1.234gamma", 0, true),
                          std::make_tuple("1.2.3", 0, true),
                          std::make_tuple("1e2.3 ", 0, true),
                          std::make_tuple("1 ", 0, true),
                          std::make_tuple("1e", 0, true),
                          std::make_tuple("1e-", 0, true),
                          std::make_tuple("1e+", 0, true)));

class integerAttributeTestFixture :public ::testing::TestWithParam<std::tuple<string, integer, bool>>
{
protected:
    string attributeString;
    integer expectedValue;
    integer parsedValue;
    bool failureFlag;

    void parseString()
    {
      GEOSX_LOG_RANK_0( "Parsing string: " << attributeString );
      xmlWrapper::stringToInputVariable< integer >(parsedValue, attributeString);
      ASSERT_TRUE( expectedValue == parsedValue );
    }
};

TEST_P(integerAttributeTestFixture, testParsing)
{
    auto testParams = GetParam();
    attributeString = std::get<0>(testParams);
    expectedValue = std::get<1>(testParams);
    failureFlag = std::get<2>(testParams);

    if (failureFlag)
    {
      EXPECT_THROW( parseString(), InputError );
    }
    else
    {
      parseString();
    }
}

INSTANTIATE_TEST_CASE_P(
        integerAttributeTests,
        integerAttributeTestFixture,
        ::testing::Values(std::make_tuple("1", 1,  false),
                          std::make_tuple("-23", -23,  false),
                          std::make_tuple("4.5", 0, true),
                          std::make_tuple("4.", 0, true),
                          std::make_tuple("alpha", 0, true),
                          std::make_tuple("1beta234", 0, true),
                          std::make_tuple("1234gamma", 0, true),
                          std::make_tuple("1 ", 0, true),
                          std::make_tuple("1 2", 0, true)));



int main( int argc, char * argv[] )
{
  logger::InitializeLogger();

  int result = 0;
  testing::InitGoogleTest( &argc, argv );
  result = RUN_ALL_TESTS();

  logger::FinalizeLogger();

  return result;
}
