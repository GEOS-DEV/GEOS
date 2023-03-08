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

// Source includes
#include "../StringUtilities.hpp"

// TPL includes
#include <gtest/gtest.h>

using namespace geosx;

struct TokenizeTest
{
  string const strToTest;
  string const delims;
  bool const treatConsecutiveDelimAsOne;
  bool const preTrimStr;
  std::vector< string > const expected;

  TokenizeTest( string const & strToTest_,
                string const & delims_,
                bool treatConsecutiveDelimAsOne_,
                bool preTrimStr_,
                std::vector< string > const & expected_ ):
    strToTest( strToTest_ ),
    delims( delims_ ),
    treatConsecutiveDelimAsOne( treatConsecutiveDelimAsOne_ ),
    preTrimStr( preTrimStr_ ),
    expected( expected_ )
  {}

  friend std::ostream & operator<<( std::ostream & os, TokenizeTest const & test )
  {
    os << "TokenizeTests( \"" << test.strToTest
       << "\", \"" << test.delims << '\"'
       << ", delimsAsOne=" << test.treatConsecutiveDelimAsOne
       << ", preTrim=" << test.preTrimStr << " )";
    return os;
  }
};

struct TokenizeBySpacesTest
{
  string const strToTest;
  std::vector< string > const expected;

  TokenizeBySpacesTest( string const & strToTest_,
                        std::vector< string > const & expected_ ):
    strToTest( strToTest_ ),
    expected( expected_ )
  {}

  friend std::ostream & operator<<( std::ostream & os, TokenizeBySpacesTest const & test )
  {
    os << "TokenizeBySpacesTest( \"" << test.strToTest << "\" )";
    return os;
  }
};

TEST( testStringUtilities, tokenize )
{

  // Path tokenizing test
  {
    map< string, std::pair< std::vector< string >, std::vector< string > > > const
    entries =
    {
      { "//entry0//entry1//entry2", { { "", "entry0", "entry1", "entry2" },
          { "", "", "entry0", "", "entry1", "", "entry2" } } },
      { "entry0//entry1/entry2", { { "entry0", "entry1", "entry2" },
          { "entry0", "", "entry1", "entry2" } } }
    };

    for( auto const & entry : entries )
    {
      string const & key = entry.first;
      std::vector< string > const & values0 = entry.second.first;
      std::vector< string > const & values1 = entry.second.second;

      std::vector< string > tokens0 = stringutilities::tokenize( key, "/", true );
      std::vector< string > tokens1 = stringutilities::tokenize( key, "/", false );


      EXPECT_TRUE( tokens0==values0 );
      EXPECT_TRUE( tokens1==values1 );
    }
  }

  // Various strings tokenizing test
  {
    std::vector< TokenizeTest > const tokenizeTests = {
      TokenizeTest( string( "a|b" ), "|", false, false, { "a", "b" } ),
      TokenizeTest( string( "a|b" ), "|", true, false, { "a", "b" } ),
      TokenizeTest( string( "a|b" ), "|", false, true, { "a", "b" } ),
      TokenizeTest( string( "a|b" ), "|", true, true, { "a", "b" } ),

      TokenizeTest( string( "|a|b|" ), "|", false, false, { "", "a", "b", "" } ),
      TokenizeTest( string( "|a|b|" ), "|", true, false, { "", "a", "b", "" } ),
      TokenizeTest( string( "|a|b|" ), "|", false, true, { "a", "b" } ),
      TokenizeTest( string( "|a|b|" ), "|", true, true, { "a", "b" } ),

      TokenizeTest( string( "|a||b|" ), "|", false, false, { "", "a", "", "b", "" } ),
      TokenizeTest( string( "|a||b|" ), "|", true, false, { "", "a", "b", "" } ),
      TokenizeTest( string( "|a||b|" ), "|", false, true, { "a", "", "b" } ),
      TokenizeTest( string( "|a||b|" ), "|", true, true, { "a", "b" } ),

      TokenizeTest( string( "a|||b" ), "|", false, false, { "a", "", "", "b" } ),
      TokenizeTest( string( "a|||b" ), "|", true, false, { "a", "b" } ),
      TokenizeTest( string( "a|||b" ), "|", false, true, { "a", "", "", "b" } ),
      TokenizeTest( string( "a|||b" ), "|", true, true, { "a", "b" } ),

      TokenizeTest( string( "||a|||b||" ), "|", false, false, { "", "", "a", "", "", "b", "", "" } ),
      TokenizeTest( string( "||a|||b||" ), "|", true, false, { "", "a", "b", "" } ),
      TokenizeTest( string( "||a|||b||" ), "|", false, true, { "a", "", "", "b" } ),
      TokenizeTest( string( "||a|||b||" ), "|", true, true, { "a", "b" } ),

      TokenizeTest( string( "|" ), "|", false, false, { "", "" } ),
      TokenizeTest( string( "|" ), "|", true, false, { "", "" } ),
      TokenizeTest( string( "|" ), "|", false, true, { } ),
      TokenizeTest( string( "|" ), "|", true, true, { } ),

      TokenizeTest( string( "|||" ), "|", false, false, { "", "", "", "" } ),
      TokenizeTest( string( "|||" ), "|", true, false, { "", "" } ),
      TokenizeTest( string( "|||" ), "|", false, true, { } ),
      TokenizeTest( string( "|||" ), "|", true, true, { } ),

      TokenizeTest( string( "ab" ), "|", false, false, { "ab" } ),
      TokenizeTest( string( "ab" ), "|", true, false, { "ab" } ),
      TokenizeTest( string( "ab" ), "|", false, true, { "ab" } ),
      TokenizeTest( string( "ab" ), "|", true, true, { "ab" } ),

      TokenizeTest( string( "" ), "|", false, false, { } ),
      TokenizeTest( string( "" ), "|", true, false, { } ),
      TokenizeTest( string( "" ), "|", false, true, { } ),
      TokenizeTest( string( "" ), "|", true, true, { } ),

      TokenizeTest( string( "" ), "", false, false, { } ),
      TokenizeTest( string( "" ), "", true, false, { } ),
      TokenizeTest( string( "" ), "", false, true, { } ),
      TokenizeTest( string( "" ), "", true, true, { } ),

      TokenizeTest( string( "ab" ), "", false, false, { "ab" } ),
      TokenizeTest( string( "ab" ), "", true, false, { "ab" } ),
      TokenizeTest( string( "ab" ), "", false, true, { "ab" } ),
      TokenizeTest( string( "ab" ), "", true, true, { "ab" } ),

      TokenizeTest( string( "| a|  b |" ), "| ", false, false, { "", "", "a", "", "", "b", "", "" } ),
      TokenizeTest( string( "| a|  b |" ), "| ", true, false, { "", "a", "b", "" } ),
      TokenizeTest( string( "| a|  b |" ), "| ", false, true, { "a", "", "", "b" } ),
      TokenizeTest( string( "| a|  b |" ), "| ", true, true, { "a", "b" } ),

      TokenizeTest( string( " |a ||b| " ), "| ", false, false, { "", "", "a", "", "", "b", "", "" } ),
      TokenizeTest( string( " |a ||b| " ), "| ", true, false, { "", "a", "b", "" } ),
      TokenizeTest( string( " |a ||b| " ), "| ", false, true, { "a", "", "", "b" } ),
      TokenizeTest( string( " |a ||b| " ), "| ", true, true, { "a", "b" } ),
    };

    for( TokenizeTest test : tokenizeTests )
    {
      std::vector< string > const r = stringutilities::tokenize( test.strToTest,
                                                                 test.delims,
                                                                 test.treatConsecutiveDelimAsOne,
                                                                 test.preTrimStr );

      if( r != test.expected )
      {
        string const result = "{ '" + stringutilities::join( r, "' , '" ) + "' }";
        string const expected = "{ '" + stringutilities::join( test.expected, "' , '" ) + "' }";
        FAIL() << test << " failed: " << result <<" results instead of " << expected << ".";
      }
    }

  }

  // Spaces tokenizing test
  {
    std::vector< TokenizeBySpacesTest > const tokenizeBSTests = {
      TokenizeBySpacesTest( string( "a b" ), { "a", "b" } ),
      TokenizeBySpacesTest( string( " a b " ), { "a", "b" } ),
      TokenizeBySpacesTest( string( "  a  b  " ), { "a", "b" } ),
      TokenizeBySpacesTest( string( "\ta \tb\n\r" ), { "a", "b" } ),
      TokenizeBySpacesTest( string( " " ), { } ),
      TokenizeBySpacesTest( string( "   " ), { } ),
      TokenizeBySpacesTest( string( "ab" ), { "ab" } ),
      TokenizeBySpacesTest( string( "\t ab \n\r" ), { "ab" } ),
    };

    for( TokenizeBySpacesTest test : tokenizeBSTests )
    {
      std::vector< string > const r = stringutilities::tokenizeBySpaces( test.strToTest );

      if( r != test.expected )
      {
        string const result = "{ '" + stringutilities::concat( r, "' , '" ) + "' }";
        string const expected = "{ '" + stringutilities::concat( test.expected, "' , '" ) + "' }";
        FAIL() << test << " failed: " << result << " results instead of " << expected << ".";
      }
    }
  }

}
