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

#include <gtest/gtest.h>

#include "dataRepository/xmlWrapper.hpp"
#include "common/format/EnumStrings.hpp"

#include <filesystem>
#include <fstream>

using namespace geos;

TEST( testXmlWrapper, array3d_errors )
{
  Regex const & array3DIntRegex = rtTypes::getTypeRegex< array3d< int > >();
  array3d< int > array;

  {
    stdVector< string > workingInputs = {
      // testing without spaces with various array sizes
      "{{{0}}}",
      "{{{0,1,2}}}",
      "{{{0,1,2},{3,4,5},{3,4,5}}}",
      "{{{0,1,2},{3,4,5},{3,4,5}},{{6,7,8},{9,10,11},{12,13,14}},{{15,16,17},{18,19,20},{21,22,23}}}",
      //testing with spaces
      " { { { 0 } } } ",
      " { { { 0, 1, 2 } } } ",
      " { { { 0 , 1 , 2 } } } ",
      " { { { 0, 1, 2 }, { 3, 4, 5 }, { 3, 4, 5 } } } ",
      " { { { 0 , 1 , 2 } , { 3 , 4 , 5 } , { 3 , 4 , 5 } } } ",
      " { { { 0, 1, 2 }, { 3, 4, 5 }, { 3, 4, 5 } }, { { 6, 7, 8 }, { 9, 10, 11 }, { 12, 13, 14 } }, { { 15, 16, 17 }, { 18, 19, 20 }, { 21, 22, 23 } } } ",
      " { { { 0 , 1 , 2 } , { 3 , 4 , 5 } , { 3 , 4 , 5 } } , { { 6 , 7 , 8 } , { 9 , 10 , 11 } , { 12 , 13 , 14 } } , { { 15 , 16 , 17 } , { 18 , 19 , 20 } , { 21 , 22 , 23 } } } ",
      "   {   {   {   0   }   }   }   ",
      "   {   {   {   0,   1,   2   }   }   }   ",
      "   {   {   {   0   ,   1   ,   2   }   }   }   ",
      "   {   {   {   0,   1,   2   },   {   3,   4,   5   },   {   3,   4,   5   }   }   }   ",
      "   {   {   {   0   ,   1   ,   2   }   ,   {   3   ,   4   ,   5   }   ,   {   3   ,   4   ,   5   }   }   }   ",
      "   {   {   {   0,   1,   2   },   {   3,   4,   5   },   {   3,   4,   5   }   },   {   {   6,   7,   8   },   {   9,   10,   11   },   {   12,   13,   14   }   },   {   {   15,   16,   17   },   {   18,   19,   20   },   {   21,   22,   23   }   }   }   ",
      "   {   {   {   0   ,   1   ,   2   }   ,   {   3   ,   4   ,   5   }   ,   {   3   ,   4   ,   5   }   }   ,   {   {   6   ,   7   ,   8   }   ,   {   9   ,   10   ,   11   }   ,   {   12   ,   13   ,   14   }   }   ,   {   {   15   ,   16   ,   17   }   ,   {   18   ,   19   ,   20   }   ,   {   21   ,   22   ,   23   }   }   }   ",
    };
    for( string const & input : workingInputs )
    {
      EXPECT_NO_THROW( xmlWrapper::stringToInputVariable( array, input, array3DIntRegex ) );
    }
  }
  {
    stdVector< string > erroneousInputs = {
      // fordbiden characters
      "{{{0,1,2},{3, hello,5},{3,4,5}},{{6,7,8},{9,10,11},{12,13,14}},{{15,16,17},{18,19,20},{21,22,23}}}",
      "{{{0,1,2 + 2},{3, 2 * 2,5},{3,4,5}},{{6,7,8},{9,10,11},{12,13,14}},{{15,16,17},{18,19,20},{21,22,23}}}",
      // illegal empty entries
      "{{{ }}}",
      "{{{0, ,2}}}",
      "{{{0,1,2}, ,{3,4,5}}}",
      "{{{0,1,2},{3,4,5},{3,4,5}}, ,{{15,16,17},{18,19,20},{21,22,23}}}",
      // unclosed brackets (these should fail the num('{')==num('}') test)
      " { { {0,1,2},{3,4,5} }, { {6,7,8},{9,10,11} }, { {12,13,14},{15,16,17} } , { {18,19,20},{21,22,23} } ",
      " { { {0,1,2},{3,4,5} }, { {6,7,8},{9,10,11} }, { {12,13,14},{15,16,17}  , { {18,19,20},{21,22,23} } }",
      " { { {0,1,2},{3,4,5} }, { {6,7,8},{9,10,11} }, { {12,13,14,{15,16,17} } , { {18,19,20},{21,22,23} } }",
      " { { {0,1,2},{3,4,5} }, { {6,7,8,{9,10,11} }, { {12,13,14},{15,16,17} } , { {18,19,20},{21,22,23} } }",
      " { { 0,1,2},{3,4,5} }, { {6,7,8},{9,10,11} }, { {12,13,14},{15,16,17} } , { {18,19,20},{21,22,23} } }",
      "  { {0,1,2},{3,4,5} }, { {6,7,8},{9,10,11} }, { {12,13,14},{15,16,17} } , { {18,19,20},{21,22,23} } ",
      " { { {,1,2},{3,4,5} }, { {6,7,8},{9,10,11} }, { {12,13,14},{15,16,17} } , { {18,19,20},{21,22,23} } }",
      " { { {},{3,4,5} }, { {6,7,8},{9,10,11} }, { {12,13,14},{15,16,17} } , { {18,19,20},{21,22,23} } }",
      " { { {0,1,2}}{ } }",
      // comma mistakes (arbitrary replaced with spaces)
      "{{{0,1,2},{3,4   5},{3,4,5}},{{6,7,8},{9,10,11},{12,13,14}},{{15,16,17},{18,19,20},{21,22,23}}}",
      "{{{0,1,2},{3,4,5}   {3,4,5}},{{6,7,8},{9,10,11},{12,13,14}},{{15,16,17},{18,19,20},{21,22,23}}}",
      "{{{0,1,2},{3,4,5},{3,4,5}}   {{6,7,8},{9,10,11},{12,13,14}},{{15,16,17},{18,19,20},{21,22,23}}}",
    };
    for( string const & input : erroneousInputs )
    {
      EXPECT_ANY_THROW( xmlWrapper::stringToInputVariable( array, input, array3DIntRegex ) );
    }
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
  EXPECT_NO_THROW( xmlWrapper::stringToInputVariable( array, input,
                                                      rtTypes::getTypeRegex< array3d< int > >() ) );

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
  EXPECT_NO_THROW( xmlWrapper::stringToInputVariable( output, input,
                                                      rtTypes::getTypeRegex< Tensor< T, N > >() ) );
  EXPECT_EQ( output, expected );
}

template< typename T, int N >
void testTensorIllFormed( string const & input )
{
  Tensor< T, N > output;
  EXPECT_THROW( xmlWrapper::stringToInputVariable( output, input,
                                                   rtTypes::getTypeRegex< Tensor< T, N > >() ), InputError );
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


// Templated attribute test
template< typename T >
class AttributeReadTestFixture : public ::testing::TestWithParam< std::tuple< string, T, bool > >
{
public:

  std::tuple< string, T, bool > testParams;
  string attributeString;
  T expectedValue;
  T parsedValue;
  bool failureFlag;

  void compareValues()
  {
    EXPECT_EQ( expectedValue, parsedValue );
  }

  void parseString()
  {
    SCOPED_TRACE( "Parsing string: " + attributeString );
    xmlWrapper::stringToInputVariable< T >( parsedValue, attributeString,
                                            rtTypes::getTypeRegex< T >() );
    compareValues();
  }

  void test()
  {
    attributeString = std::get< 0 >( testParams );
    expectedValue = std::get< 1 >( testParams );
    failureFlag = std::get< 2 >( testParams );

    if( failureFlag )
    {
      EXPECT_THROW( parseString(), InputError );
    }
    else
    {
      EXPECT_NO_THROW( parseString() );
    }
  }
};


class real64AttributeTestFixture : public AttributeReadTestFixture< real64 >
{
  void compareValues()
  {
    EXPECT_DOUBLE_EQ( expectedValue, parsedValue );
  }
};

TEST_P( real64AttributeTestFixture, testParsing )
{
  testParams = GetParam();
  this->test();
}

INSTANTIATE_TEST_SUITE_P(
  real64AttributeTests,
  real64AttributeTestFixture,
  ::testing::Values( std::make_tuple( "1", 1, false ),
                     std::make_tuple( "1 ", 1, false ),
                     std::make_tuple( " 1", 1, false ),
                     std::make_tuple( " 1 ", 1, false ),
                     std::make_tuple( "-23", -23, false ),
                     std::make_tuple( "4.5", 4.5, false ),
                     std::make_tuple( "4.", 4.0, false ),
                     std::make_tuple( "6e1", 6e1, false ),
                     std::make_tuple( "7e-2", 7e-2, false ),
                     std::make_tuple( "8.765e0", 8.765, false ),
                     std::make_tuple( "-1.2e-3", -1.2e-3, false ),
                     std::make_tuple( "3.5E+4", 3.5e+4, false ),
                     std::make_tuple( "6.1E05", 6.1e05, false ),
                     std::make_tuple( "-7.1e06", -7.1e06, false ),
                     std::make_tuple( "8.1e+02", 8.1e+02, false ),
                     std::make_tuple( "9.1e-01", 9.1e-01, false ),
                     std::make_tuple( "alpha", 0, true ),
                     std::make_tuple( "1beta234", 0, true ),
                     std::make_tuple( "1.234gamma", 0, true ),
                     std::make_tuple( "1.2.3", 0, true ),
                     std::make_tuple( "1e2.3 ", 0, true ),
                     std::make_tuple( "1e", 0, true ),
                     std::make_tuple( "1e-", 0, true ),
                     std::make_tuple( "1e+", 0, true )));


class real32AttributeTestFixture : public AttributeReadTestFixture< real32 >
{
  void compareValues()
  {
    EXPECT_FLOAT_EQ( expectedValue, parsedValue );
  }
};

TEST_P( real32AttributeTestFixture, testParsing )
{
  testParams = GetParam();
  this->test();
}

INSTANTIATE_TEST_SUITE_P(
  real32AttributeTests,
  real32AttributeTestFixture,
  ::testing::Values( std::make_tuple( "1", 1, false ),
                     std::make_tuple( "1 ", 1, false ),
                     std::make_tuple( " 1", 1, false ),
                     std::make_tuple( " 1 ", 1, false ),
                     std::make_tuple( "-23", -23, false ),
                     std::make_tuple( "4.5", 4.5, false ),
                     std::make_tuple( "4.", 4.0, false ),
                     std::make_tuple( "6e1", 6e1, false ),
                     std::make_tuple( "7e-2", 7e-2, false ),
                     std::make_tuple( "8.765e0", 8.765, false ),
                     std::make_tuple( "-1.2e-3", -1.2e-3, false ),
                     std::make_tuple( "3.5E+4", 3.5e+4, false ),
                     std::make_tuple( "6.1E05", 6.1e05, false ),
                     std::make_tuple( "-7.1e06", -7.1e06, false ),
                     std::make_tuple( "8.1e+02", 8.1e+02, false ),
                     std::make_tuple( "9.1e-01", 9.1e-01, false ),
                     std::make_tuple( "alpha", 0, true ),
                     std::make_tuple( "1beta234", 0, true ),
                     std::make_tuple( "1.234gamma", 0, true ),
                     std::make_tuple( "1.2.3", 0, true ),
                     std::make_tuple( "1e2.3 ", 0, true ),
                     std::make_tuple( "1e", 0, true ),
                     std::make_tuple( "1e-", 0, true ),
                     std::make_tuple( "1e+", 0, true )));


class integerAttributeTestFixture : public AttributeReadTestFixture< integer > {};

TEST_P( integerAttributeTestFixture, testParsing )
{
  testParams = GetParam();
  this->test();
}

INSTANTIATE_TEST_SUITE_P(
  integerAttributeTests,
  integerAttributeTestFixture,
  ::testing::Values( std::make_tuple( "1", 1, false ),
                     std::make_tuple( "1 ", 1, false ),
                     std::make_tuple( " 1", 1, false ),
                     std::make_tuple( " 1 ", 1, false ),
                     std::make_tuple( "-23", -23, false ),
                     std::make_tuple( "4.5", 0, true ),
                     std::make_tuple( "4.", 0, true ),
                     std::make_tuple( "alpha", 0, true ),
                     std::make_tuple( "1beta234", 0, true ),
                     std::make_tuple( "1234gamma", 0, true ),
                     std::make_tuple( "1 2", 0, true )));


enum class TestEnum { None, Default, Value, Value2 };
ENUM_STRINGS( TestEnum, "None", "Default", "Value", "Value2" );

class enumAttributeTestFixture : public AttributeReadTestFixture< TestEnum > {};

TEST_P( enumAttributeTestFixture, testParsing )
{
  testParams = GetParam();
  this->test();
}

INSTANTIATE_TEST_SUITE_P(
  enumAttributeTests,
  enumAttributeTestFixture,
  ::testing::Values( std::make_tuple( "None", TestEnum::None, false ),
                     std::make_tuple( "Default", TestEnum::Default, false ),
                     std::make_tuple( "Value", TestEnum::Value, false ),
                     std::make_tuple( "Value2", TestEnum::Value2, false ),
                     std::make_tuple( "0", TestEnum( 0 ), true ),
                     std::make_tuple( "4.", TestEnum( 0 ), true ),
                     std::make_tuple( "alpha", TestEnum( 0 ), true ),
                     std::make_tuple( "Val", TestEnum( 0 ), true ),
                     std::make_tuple( "Def ault", TestEnum( 0 ), true ),
                     std::make_tuple( "None123", TestEnum( 0 ), true ) ) );


TEST( testXmlWrapper, testGroupNamesFormats )
{
  struct GroupNameTest
  {
    Regex const & m_regex;
    string m_valueToTest;
    GroupNameTest( Regex const & regex, string_view valueToTest ):
      m_regex( regex ), m_valueToTest( valueToTest ) {}
  };

  Regex const & groupNameRegex = rtTypes::getTypeRegex< string >( rtTypes::CustomTypes::groupName );
  string groupName;

  {
    stdVector< GroupNameTest > workingInputs = {
      GroupNameTest( groupNameRegex, "testname" ),
      GroupNameTest( groupNameRegex, "testname01" ),
      GroupNameTest( groupNameRegex, "test_name" ),
      GroupNameTest( groupNameRegex, "test-name" ),
      GroupNameTest( groupNameRegex, " testname " ),
      GroupNameTest( groupNameRegex, " \t testname \n " ),
    };
    for( GroupNameTest const & input : workingInputs )
    {
      EXPECT_NO_THROW( xmlWrapper::stringToInputVariable( groupName, input.m_valueToTest, input.m_regex ) )
        << "Parsing input '"<< input.m_valueToTest
        << "' with regex '" << input.m_regex.m_regexStr << "' didn't throw an InputError as expected.";
      string_view stringTrimed = stringutilities::trimSpaces( input.m_valueToTest );
      EXPECT_STREQ( groupName.c_str(), std::string{stringTrimed}.c_str()  );
    }
  }
  {
    stdVector< GroupNameTest > erroneousInputs = {
      //empty entries
      GroupNameTest( groupNameRegex, "" ),
      GroupNameTest( groupNameRegex, " " ),
      GroupNameTest( groupNameRegex, "\t" ),
      //white spaces
      GroupNameTest( groupNameRegex, "test name" ),
      GroupNameTest( groupNameRegex, "test\tname" ),
      //fordbiden characters
      GroupNameTest( groupNameRegex, "test/name" ),
      GroupNameTest( groupNameRegex, "test:name" ),
      GroupNameTest( groupNameRegex, "test;name" ),
      GroupNameTest( groupNameRegex, "test\\name" ),
      GroupNameTest( groupNameRegex, "{[arrayElement]}" ),
      GroupNameTest( groupNameRegex, "{name1,name2,name3}" ),
    };
    for( GroupNameTest const & input : erroneousInputs )
    {
      EXPECT_THROW( xmlWrapper::stringToInputVariable( groupName, input.m_valueToTest, input.m_regex ),
                    InputError ) << "Parsing input '"<< input.m_valueToTest
                                 << "' with regex '" << input.m_regex.m_regexStr << "' didn't throw an InputError as expected.";
    }
  }
}
TEST( testXmlWrapper, testGroupNamesArrayFormats )
{
  struct GroupNameTest
  {
    Regex const & m_regex;
    string m_valueToTest;
    GroupNameTest( Regex const & regex, string_view valueToTest ):
      m_regex( regex ), m_valueToTest( valueToTest ) {}
  };

  Regex const & groupNameRefArrayRegex = rtTypes::getTypeRegex< string >( rtTypes::CustomTypes::groupNameRefArray );
  string groupName;

  {
    stdVector< GroupNameTest > workingInputs = {
      GroupNameTest( groupNameRefArrayRegex, "{}" ),
      GroupNameTest( groupNameRefArrayRegex, " {} " ),
      GroupNameTest( groupNameRefArrayRegex, " \t {} \n " ),
      GroupNameTest( groupNameRefArrayRegex, "{groupName}" ),
      GroupNameTest( groupNameRefArrayRegex, "{123name}" ),
      GroupNameTest( groupNameRefArrayRegex, "{name123}" ),
      GroupNameTest( groupNameRefArrayRegex, "{a-Z0-9./*[]-_,}" ),
      GroupNameTest( groupNameRefArrayRegex, "{name.with-special_chars}" ),
      GroupNameTest( groupNameRefArrayRegex, "{path/to/resource*}" ),
      GroupNameTest( groupNameRefArrayRegex, "{[arrayElement]}" ),
      GroupNameTest( groupNameRefArrayRegex, "{name1,name2,name3}" ),
    };
    for( GroupNameTest const & input : workingInputs )
    {
      EXPECT_NO_THROW( xmlWrapper::stringToInputVariable( groupName, input.m_valueToTest, input.m_regex ) );
    }
  }
  {
    stdVector< GroupNameTest > erroneousInputs = {
      GroupNameTest( groupNameRefArrayRegex, "" ),
      GroupNameTest( groupNameRefArrayRegex, " " ),
      GroupNameTest( groupNameRefArrayRegex, " \t " ),
      GroupNameTest( groupNameRefArrayRegex, "{\t}" ),
      GroupNameTest( groupNameRefArrayRegex, "{test:name}" ),
      GroupNameTest( groupNameRefArrayRegex, "{test;name}" ),
      GroupNameTest( groupNameRefArrayRegex, "{test\\name}" ),
      GroupNameTest( groupNameRefArrayRegex, "{test|name}" ),
      GroupNameTest( groupNameRefArrayRegex, "{test^name}" ),
      GroupNameTest( groupNameRefArrayRegex, "{test$name}" ),
      GroupNameTest( groupNameRefArrayRegex, "{test&name}" ),
      GroupNameTest( groupNameRefArrayRegex, "{test#name}" ),
      GroupNameTest( groupNameRefArrayRegex, "{test!name}" ),
      GroupNameTest( groupNameRefArrayRegex, "{test%name}" ),
      GroupNameTest( groupNameRefArrayRegex, "{test@name}" ),
      GroupNameTest( groupNameRefArrayRegex, "{test(name}" ),
      GroupNameTest( groupNameRefArrayRegex, "{test)name}" ),
      GroupNameTest( groupNameRefArrayRegex, "{test=name}" ),
      GroupNameTest( groupNameRefArrayRegex, "{test+name}" ),
      GroupNameTest( groupNameRefArrayRegex, "{test<name}" ),
      GroupNameTest( groupNameRefArrayRegex, "{test>name}" ),
      GroupNameTest( groupNameRefArrayRegex, "{test\tname}" ),
      GroupNameTest( groupNameRefArrayRegex, "{test\nname}" ),
      GroupNameTest( groupNameRefArrayRegex, "{test\rname}" ),
      GroupNameTest( groupNameRefArrayRegex, "groupName" ),
      GroupNameTest( groupNameRefArrayRegex, "{groupName" ),
      GroupNameTest( groupNameRefArrayRegex, "groupName}" ),
      GroupNameTest( groupNameRefArrayRegex, "{groupName}} " ),
      GroupNameTest( groupNameRefArrayRegex, "{}groupName" ),
      GroupNameTest( groupNameRefArrayRegex, "groupName{}" ),
      GroupNameTest( groupNameRefArrayRegex, "{hello, \t\n\r ,world}" ),
      GroupNameTest( groupNameRefArrayRegex, "test{groupName}" ),
      GroupNameTest( groupNameRefArrayRegex, "{groupName}test" ),
      GroupNameTest( groupNameRefArrayRegex, "{groupName} test" ),
      GroupNameTest( groupNameRefArrayRegex, "{element with space, another}" ),
      GroupNameTest( groupNameRefArrayRegex, "{element, element with space, 123.456, a-b}" ),
      GroupNameTest( groupNameRefArrayRegex, "{ space at ends } " ),
      GroupNameTest( groupNameRefArrayRegex, "{valuewith,,commas }" ),
      GroupNameTest( groupNameRefArrayRegex, "{ value with , commas }" ),
      GroupNameTest( groupNameRefArrayRegex, "{{groupname}}" ),
    };
    for( GroupNameTest const & input : erroneousInputs )
    {
      EXPECT_THROW( xmlWrapper::stringToInputVariable( groupName, input.m_valueToTest, input.m_regex ),
                    InputError ) << "Parsing input '"<< input.m_valueToTest
                                 << "' with regex '" << input.m_regex.m_regexStr << "' didn't throw an InputError as expected.";
    }
  }
}

class CollectIncludedTest : public ::testing::Test
{
protected:
  std::filesystem::path m_tempDir;

  void SetUp() override
  {
    m_tempDir = std::filesystem::temp_directory_path() / "geos_collectIncluded_test";
    std::filesystem::create_directories( m_tempDir );
  }

  void TearDown() override
  {
    std::filesystem::remove_all( m_tempDir );
  }

  string filePath( string const & filename )
  {
    return ( m_tempDir / filename ).string();
  }

  void writeXML( string const & filename, string const & content )
  {
    std::ofstream f( filePath( filename ) );
    f << content;
  }
};


TEST_F( CollectIncludedTest, collectIncluded_noIncludes )
{
  writeXML( "base.xml", "<Problem>"
                        "</Problem>" );

  auto result = xmlWrapper::collectIncluded( filePath( "base.xml" ) );

  EXPECT_TRUE( result.empty() );
}

TEST_F( CollectIncludedTest, collectIncluded_singleInclude )
{
  writeXML( "child.xml", "<Problem>"
                         "</Problem>" );
  writeXML( "base.xml", "<Problem>"
                        " <Included>"
                        "   <File name=\"child.xml\"/>"
                        " </Included>"
                        "</Problem>" );

  auto result = xmlWrapper::collectIncluded( filePath( "base.xml" ) );

  EXPECT_NE( result.find( filePath( "child.xml" ) ), result.end() );
}

TEST_F( CollectIncludedTest, collectIncluded_multipleIncludes )
{
  writeXML( "child1.xml", "<Problem>"
                          "</Problem>" );
  writeXML( "child2.xml", "<Problem>"
                          "</Problem>" );
  writeXML( "base.xml", "<Problem>"
                        " <Included>"
                        "   <File name=\"child1.xml\"/>"
                        "   <File name=\"child2.xml\"/>"
                        " </Included>"
                        "</Problem>" );

  auto result = xmlWrapper::collectIncluded( filePath( "base.xml" ) );

  EXPECT_NE( result.find( filePath( "child1.xml" ) ), result.end() );
  EXPECT_NE( result.find( filePath( "child2.xml" ) ), result.end() );
}

TEST_F( CollectIncludedTest, collectIncluded_emptyNameAttribute )
{
  writeXML( "base.xml", "<Problem>"
                        " <Included>"
                        "   <File name=\"\"/>"
                        " </Included>"
                        "</Problem>" );

  std::set< string > result;

  EXPECT_ANY_THROW( xmlWrapper::collectIncluded( filePath( "base.xml" ), result ) );

  EXPECT_TRUE( result.empty() );
}

TEST_F( CollectIncludedTest, collectIncluded_existingEntriesKept )
{
  std::set< string > existingCollection;
  existingCollection.insert( "/somewhere/thereisanalreadyexistingxmlfile.xml" );
  writeXML( "base.xml", "<Problem>"
                        "</Problem>" );

  xmlWrapper::collectIncluded( filePath( "base.xml" ), existingCollection );

  EXPECT_NE( existingCollection.find( "/somewhere/thereisanalreadyexistingxmlfile.xml" ),
             existingCollection.end() );
}


TEST_F( CollectIncludedTest, collectIncludedRecursive_noIncludes )
{
  writeXML( "base.xml", "<Problem>"
                        "</Problem>" );

  auto result = xmlWrapper::collectIncludedRecursive( filePath( "base.xml" ) );

  EXPECT_EQ( result.size(), 1 );  // size 1 because collectIncludedRecursive collects the base file
}

TEST_F( CollectIncludedTest, collectIncludedRecursive_singleInclude )
{
  writeXML( "child.xml", "<Problem>"
                         "</Problem>" );
  writeXML( "base.xml", "<Problem>"
                        " <Included>"
                        "   <File name=\"child.xml\"/>"
                        " </Included>"
                        "</Problem>" );

  auto result = xmlWrapper::collectIncludedRecursive( filePath( "base.xml" ) );

  EXPECT_NE( result.find( filePath( "child.xml" ) ), result.end() );
}

TEST_F( CollectIncludedTest, collectIncludedRecursive_multipleIncludes )
{
  writeXML( "child1.xml", "<Problem>"
                          "</Problem>" );
  writeXML( "child2.xml", "<Problem>"
                          "</Problem>" );
  writeXML( "base.xml", "<Problem>"
                        " <Included>"
                        "   <File name=\"child1.xml\"/>"
                        "   <File name=\"child2.xml\"/>"
                        " </Included>"
                        "</Problem>" );

  auto result = xmlWrapper::collectIncludedRecursive( filePath( "base.xml" ) );

  EXPECT_NE( result.find( filePath( "child1.xml" ) ), result.end() );
  EXPECT_NE( result.find( filePath( "child2.xml" ) ), result.end() );
}

TEST_F( CollectIncludedTest, collectIncludedRecursive_simpleRecursive )
{
  writeXML( "child.xml", "<Problem>"
                         "</Problem>" );

  writeXML( "middle.xml", "<Problem>"
                          " <Included>"
                          "   <File name=\"child.xml\"/>"
                          " </Included>"
                          "</Problem>" );

  writeXML( "base.xml", "<Problem>"
                        " <Included>"
                        "   <File name=\"middle.xml\"/>"
                        " </Included>"
                        "</Problem>" );

  auto result = xmlWrapper::collectIncludedRecursive( filePath( "base.xml" ) );

  EXPECT_NE( result.find( filePath( "middle.xml" ) ), result.end() );
  EXPECT_NE( result.find( filePath( "child.xml" ) ), result.end() );
}

TEST_F( CollectIncludedTest, collectIncludedRecursive_cyclePrevention )
{
  writeXML( "cycle.xml", "<Problem>"
                         " <Included>"
                         "   <File name=\"cycle.xml\"/>"
                         " </Included>"
                         "</Problem>" );

  auto result = xmlWrapper::collectIncludedRecursive( filePath( "cycle.xml" ) );

  EXPECT_NE( result.find( filePath( "cycle.xml" ) ), result.end() );
}

TEST_F( CollectIncludedTest, collectIncludedRecursive_noDuplicates )
{
  writeXML( "base.xml", "<Problem>"
                        " <Included>"
                        "   <File name=\"base.xml\"/>"
                        " </Included>"
                        "</Problem>" );

  auto result = xmlWrapper::collectIncludedRecursive( filePath( "base.xml" ) );

  EXPECT_EQ( result.size(), 1 );  // collectIncludedRecursive collects the base file
}

TEST_F( CollectIncludedTest, collectIncludedRecursive_existingEntriesKept )
{
  std::set< string > existingCollection;
  existingCollection.insert( "/somewhere/thereisanalreadyexistingxmlfile.xml" );
  writeXML( "base.xml", "<Problem>"
                        "</Problem>" );

  xmlWrapper::collectIncludedRecursive( filePath( "base.xml" ), existingCollection );

  EXPECT_NE( existingCollection.find( "/somewhere/thereisanalreadyexistingxmlfile.xml" ),
             existingCollection.end() );
}


int main( int argc, char * argv[] )
{
  logger::InitializeLogger();

  int result = 0;
  testing::InitGoogleTest( &argc, argv );
  result = RUN_ALL_TESTS();

  logger::FinalizeLogger();

  return result;
}
