/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
    LvArray::Array< int, 3, localIndex > array;
    xmlWrapper::StringToInputVariable( array, input );
  }
  // This should fail the num('{')==num('}') test
  {
    input = " { { {0,1,2},{3,4,5} }, { {6,7,8},{9,10,11} }, { {12,13,14},{15,16,17} } , { {18,19,20},{21,22,23} } ";
    LvArray::Array< int, 3, localIndex > array;
    EXPECT_DEATH_IF_SUPPORTED( xmlWrapper::StringToInputVariable( array, input ), IGNORE_OUTPUT );
  }
  {
    input = " { { {0,1,2},{3,4,5} }, { {6,7,8},{9,10,11} }, { {12,13,14},{15,16,17}  , { {18,19,20},{21,22,23} } }";
    LvArray::Array< int, 3, localIndex > array;
    EXPECT_DEATH_IF_SUPPORTED( xmlWrapper::StringToInputVariable( array, input ), IGNORE_OUTPUT );
  }
  {
    input = " { { {0,1,2},{3,4,5} }, { {6,7,8},{9,10,11} }, { {12,13,14,{15,16,17} } , { {18,19,20},{21,22,23} } }";
    LvArray::Array< int, 3, localIndex > array;
    EXPECT_DEATH_IF_SUPPORTED( xmlWrapper::StringToInputVariable( array, input ), IGNORE_OUTPUT );
  }
  {
    input = " { { {0,1,2},{3,4,5} }, { {6,7,8,{9,10,11} }, { {12,13,14},{15,16,17} } , { {18,19,20},{21,22,23} } }";
    LvArray::Array< int, 3, localIndex > array;
    EXPECT_DEATH_IF_SUPPORTED( xmlWrapper::StringToInputVariable( array, input ), IGNORE_OUTPUT );
  }
  {
    input = " { { 0,1,2},{3,4,5} }, { {6,7,8},{9,10,11} }, { {12,13,14},{15,16,17} } , { {18,19,20},{21,22,23} } }";
    LvArray::Array< int, 3, localIndex > array;
    EXPECT_DEATH_IF_SUPPORTED( xmlWrapper::StringToInputVariable( array, input ), IGNORE_OUTPUT );
  }
  {
    input = "  { {0,1,2},{3,4,5} }, { {6,7,8},{9,10,11} }, { {12,13,14},{15,16,17} } , { {18,19,20},{21,22,23} } ";
    LvArray::Array< int, 3, localIndex > array;
    EXPECT_DEATH_IF_SUPPORTED( xmlWrapper::StringToInputVariable( array, input ), IGNORE_OUTPUT );
  }

  {
    input = " { { {,1,2},{3,4,5} }, { {6,7,8},{9,10,11} }, { {12,13,14},{15,16,17} } , { {18,19,20},{21,22,23} } }";
    LvArray::Array< int, 3, localIndex > array;
    EXPECT_DEATH_IF_SUPPORTED( xmlWrapper::StringToInputVariable( array, input ), IGNORE_OUTPUT );
  }

  {
    input = " { { {},{3,4,5} }, { {6,7,8},{9,10,11} }, { {12,13,14},{15,16,17} } , { {18,19,20},{21,22,23} } }";
    LvArray::Array< int, 3, localIndex > array;
    EXPECT_DEATH_IF_SUPPORTED( xmlWrapper::StringToInputVariable( array, input ), IGNORE_OUTPUT );
  }
  {
    input = " { { {0,1,2}}{ } }";
    LvArray::Array< int, 3, localIndex > array;
    EXPECT_DEATH_IF_SUPPORTED( xmlWrapper::StringToInputVariable( array, input ), IGNORE_OUTPUT );
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
  for( localIndex i=0 ; i<4 ; ++i )
  {
    input += "{ ";
    for( localIndex j=0 ; j<5 ; ++j )
    {
      input += "{ ";
      for( localIndex k=0 ; k<3 ; ++k )
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

  LvArray::Array< int, 3, localIndex > array;
  xmlWrapper::StringToInputVariable( array, input );

  ASSERT_EQ( array.size( 0 ), numI );
  ASSERT_EQ( array.size( 1 ), numJ );
  ASSERT_EQ( array.size( 2 ), numK );

  for( localIndex i=0 ; i<array.size( 0 ) ; ++i )
  {
    for( localIndex j=0 ; j<array.size( 1 ) ; ++j )
    {
      for( localIndex k=0 ; k<array.size( 2 ) ; ++k )
      {
        ASSERT_EQ( array[i][j][k], i*2+j*3+k*4 );
      }
    }
  }
}

int main( int argc, char * argv[] )
{
  logger::InitializeLogger();

  int result = 0;
  testing::InitGoogleTest( &argc, argv );
  result = RUN_ALL_TESTS();

  logger::FinalizeLogger();

#ifdef USE_CHAI
  chai::ArrayManager::finalize();
#endif

  return result;
}
