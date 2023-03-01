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
  string m_strToTest;
  string m_delims;
  bool m_treatConsecutiveDelimAsOne;
  bool m_preTrimStr;
  std::vector< string > m_results;

  TokenizeTest( string const & strToTest,
                string const & delims,
                bool treatConsecutiveDelimAsOne,
                bool preTrimStr,
                std::vector< string > const & results ):
    m_strToTest( strToTest ),
    m_delims( delims ),
    m_treatConsecutiveDelimAsOne( treatConsecutiveDelimAsOne ),
    m_preTrimStr( preTrimStr ),
    m_results( results )
  {}

  friend std::ostream & operator<<( std::ostream & os, TokenizeTest const & test )
  {
    os << "TokenizeTests( \"" << test.m_strToTest
       << "\", \"" << test.m_delims << '\"'
       << ", delimsAsOne=" << test.m_treatConsecutiveDelimAsOne
       << ", preTrim=" << test.m_preTrimStr << " )";
    return os;
  }
};

struct TokenizeBySpacesTest
{
  string m_strToTest;
  std::vector< string > m_results;

  TokenizeBySpacesTest( string const & strToTest,
                        std::vector< string > const & results ):
    m_strToTest( strToTest ),
    m_results( results )
  {}

  friend std::ostream & operator<<( std::ostream & os, TokenizeBySpacesTest const & test )
  {
    os << "TokenizeBySpacesTest( \"" << test.m_strToTest << "\" )";
    return os;
  }
};

TEST( testStringUtilities, tokenize )
{

  // Path tokenizing test
  {
    map< string, std::pair< std::vector< string >, std::vector< string > > >
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
    std::vector< TokenizeTest > tokenizeTests = {
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
      std::vector< string > r = stringutilities::tokenize( test.m_strToTest,
                                                           test.m_delims,
                                                           test.m_treatConsecutiveDelimAsOne,
                                                           test.m_preTrimStr );

      if( r.size() != test.m_results.size() )
      {
        FAIL() << test << " failed: "
               << r.size() << " results instead of " << test.m_results.size();
      }
      else
      {
        for( size_t i = 0; i < r.size(); ++i )
        {
          if( r[i] != test.m_results[i] )
          {
            FAIL() << test << " failed: "
                   << "result no." << i << " was '" << r[i] << "' instead of '" << test.m_results[i] <<"'.";
          }
        }
      }
    }

  }

  // Spaces tokenizing test
  {
    std::vector< TokenizeBySpacesTest > tokenizeBSTests = {
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
      std::vector< string > r = stringutilities::tokenizeBySpaces( test.m_strToTest );

      if( r.size() != test.m_results.size() )
      {
        FAIL() << test << " failed: "
               << r.size() << " results instead of " << test.m_results.size();
      }
      else
      {
        for( size_t i = 0; i < r.size(); ++i )
        {
          if( r[i] != test.m_results[i] )
          {
            FAIL() << test << " failed: "
                   << "result no." << i << " was '" << r[i] << "' instead of '" << test.m_results[i] <<"'.";
          }
        }
      }
    }
  }

}
