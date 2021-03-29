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

#include "StringUtilities.hpp"

#include "common/Logger.hpp"
#include <algorithm>
#include <fstream>

// TODO: this should not be here! Move what is not CPG-specific to codingUtilities/StringUtilities
// TODO: std::string -> string

namespace geosx
{

namespace CPMeshStringUtilities
{

void eclipseDataBufferToVector( std::string & inputBuffer, std::vector< double > & v )
{
  std::string const star( "*" );
  int positionOfStar( 0 );
  int nRepetitions( 0 );
  double repeatedValue( 0. );

  // split the inputBuffer in chunks of information (spaces are used to split)
  std::istringstream splitBuffer( inputBuffer );
  std::string chunk;

  // Traverse through all chunks
  do
  {
    // Read a chunk from the buffer
    splitBuffer >> chunk;
    if( chunk.size() > 0 )
    {
      // does this word contain a star "*"
      positionOfStar = int(chunk.find( star ));
      if( positionOfStar>0 )
      {
        std::string const nTimesValueIsRepeated( chunk.substr( 0, positionOfStar ));
        std::string const valueAsString( chunk.substr( positionOfStar+1, chunk.size()));
        nRepetitions = std::stoi( nTimesValueIsRepeated );
        repeatedValue = std::stod( valueAsString );
        std::fill_n( back_inserter( v ), nRepetitions, repeatedValue ); // append nRepetitions of the repeatedValue as double
      }
      else
      {
        // The value is only present once (no repetition)
        repeatedValue = std::stod( chunk ); // conversion to double
        v.push_back( repeatedValue ); // append the repeatedValue once
      }
    }
    // While there is more to read
  }
  while (splitBuffer);
}

void eclipseDataBufferToVector( std::string & inputBuffer, std::vector< int > & v )
{
  std::string const star( "*" );
  int positionOfStar( 0 );
  int nRepetitions( 0 );
  int repeatedValue( 0 );

  // split the inputBuffer in chunks of information (spaces are used to split)
  std::istringstream splitBuffer( inputBuffer );

  // Traverse through all chunks
  do
  {
    // Read a word
    std::string word;
    splitBuffer >> word;
    if( word.size()>0 )
    {
      // does this word contain a "*"
      positionOfStar = int(word.find( star ));
      if( positionOfStar>0 )
      {
        std::string const multiplierAsString( word.substr( 0, positionOfStar ));
        std::string const valueAsString( word.substr( positionOfStar+1, word.size()));
        nRepetitions = std::stoi( multiplierAsString );
        repeatedValue = std::stoi( valueAsString );
        std::fill_n( back_inserter( v ), nRepetitions, repeatedValue ); // append nRepetitions of the repeatedValue as int
      }
      else
      {
        // The value is only present once (no repetition)
        repeatedValue = std::stoi( word ); // conversion to integer
        v.push_back( repeatedValue ); // append the repeatedValue once
      }
    }
    // While there is more to read
  }
  while (splitBuffer);
}

std::string fileToString( const std::string filePath )
{
  std::ifstream meshFile;
  std::string fileContent( "A" );

  //Open file
  meshFile.open( filePath );
  GEOSX_THROW_IF( !meshFile.is_open(),
                  "File could not be open",
                  InputError );

  //Transfer file content into string for easing broadcast
  std::stringstream buffer;
  buffer << meshFile.rdbuf();
  fileContent = buffer.str();

  //Close file
  meshFile.close();
  return fileContent;
}

void trim( std::string & str )
{
  if( str.size() > 0 )
  {
    int first = int(str.find_first_not_of( ' ' ));
    int last = int(str.find_last_not_of( ' ' ));
    if( first < 0 )
    {
      str.resize( 0 );
    }
    else
    {
      str = str.substr( first, (last - first + 1));
    }

  }
}

bool removeStringAndFollowingContentFromLine( std::string toBeRemoved, std::string & line )
{
  // does the line have a "toBeRemoved" character(s)
  std::size_t pos = line.find( toBeRemoved );
  std::size_t initial_line_length = line.size();
  bool res = false;
  if( pos != std::string::npos )
  {
    res = true;

    std::size_t end_line_position = line.find( '\n' );
    if( end_line_position != std::string::npos )
    {
      // remove the character and everything afterwards
      line = line.substr( 0, pos )+line.substr( end_line_position+1, initial_line_length );
    }
    else
    {
      line = line.substr( 0, pos );
    }
  }
  return res;
}

void removeTab( std::string & v )
{
  std::replace( v.begin(), v.end(), '\t', ' ' );
}

void removeEndOfLine( std::string & v )
{
  std::replace( v.begin(), v.end(), '\r', ' ' );
  std::replace( v.begin(), v.end(), '\n', ' ' );
}

void removeExtraSpaces( std::string & v )
{
  v.erase( std::unique( v.begin(), v.end(),
                        []( char a, char b ) { return a == ' ' && b == ' '; } ), v.end());
}

} // end namespace CPMeshStringUtilities

} // end namespace geosx
