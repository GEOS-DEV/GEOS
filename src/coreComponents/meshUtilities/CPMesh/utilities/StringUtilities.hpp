/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOSX_MESHUTILITIES_CPMESH_STRINGUTILITIES_HPP_
#define GEOSX_MESHUTILITIES_CPMESH_STRINGUTILITIES_HPP_

#include <sstream>
#include <vector>

// TODO: this should not be here! Move what is not CPG-specific to codingUtilities/StringUtilities

namespace geosx
{

namespace CPMeshStringUtilities
{

void eclipseDataBufferToVector( std::string & inputBuffer, std::vector< double > & outputVector );
void eclipseDataBufferToVector( std::string & inputBuffer, std::vector< int > & outputVector );
std::string fileToString( const std::string filePath );
void trim( std::string & str );
bool removeStringAndFollowingContentFromLine( std::string toBeRemoved, std::string & line );
void removeTab( std::string & v );
void removeEndOfLine( std::string & v );
void removeExtraSpaces( std::string & v );

template< typename T >
void fromStringTo( std::string & data, std::vector< T > & v )
{
  v.reserve( data.size());
  std::istringstream iss( data );
  T sub;
  while( iss >> sub )
  {
    v.push_back( sub );
  }
}

template< typename T >
void fromStringTo( std::string & data, T & v )
{
  std::istringstream iss( data );
  iss >> v;
}

} // end namespace CPMeshStringUtilities

} // end namespace geosx

#endif
