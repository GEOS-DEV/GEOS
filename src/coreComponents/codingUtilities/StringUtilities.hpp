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

/**
 * @file StringUtilities.hpp
 */

#ifndef GEOSX_CODINGUTILITIES_STRINGUTILITIES_HPP_
#define GEOSX_CODINGUTILITIES_STRINGUTILITIES_HPP_

#include <cxxabi.h>
#include <cstring>
#include <memory>
#include <sstream>
#include <algorithm>
#include <map>


#include "common/DataTypes.hpp"
#include "cxx-utilities/src/IntegerConversion.hpp"

namespace geosx
{
namespace stringutilities
{

/// Overloaded function to check equality between strings and char arrays
/// Mainly used to avoid char*==char* mistakes
inline bool streq( std::string const & strA, std::string const & strB )
{ return strA == strB; }

inline bool streq( std::string const & strA, char const * const strB )
{ return strA == strB; }

inline bool streq( char const * const strA, std::string const & strB )
{ return strA == strB; }

inline bool streq( char const * const strA, char const * const strB )
{ return !strcmp( strA, strB ); }

/// Subdivide string by delimiters
string_array Tokenize( std::string const & str, std::string const & delimiters );

string EncodeBase64( unsigned char const * const bytes,
                     integer dataSize );

}
}

#endif /* GEOSX_CODINGUTILITIES_STRINGUTILITIES_HPP_ */
