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

/**
 * @file DataTypes.cpp
 */


#include "DataTypes.hpp"
#include "Logger.hpp"
#include "LvArray/src/system.hpp"
#include "codingUtilities/StringUtilities.hpp"

namespace geos
{
#ifdef GEOSX_USE_MPI
MPI_Comm MPI_COMM_GEOSX;
#else
int MPI_COMM_GEOSX = 0;
#endif

void printTypeSummary()
{
  GEOS_LOG_RANK_0( "real64 is alias of " <<LvArray::system::demangle( typeid(real64).name() ) );
  GEOS_LOG_RANK_0( "localIndex is alias of " <<LvArray::system::demangle( typeid(localIndex).name() ) );
  GEOS_LOG_RANK_0( "globalIndex is alias of "<<LvArray::system::demangle( typeid(globalIndex).name()) );
}

string rtTypes::getTypeName( std::type_index const key )
{
  static const std::unordered_map< std::type_index, string > type_names =
  {
    {std::type_index( typeid(integer)), "integer"},
    {std::type_index( typeid(real32)), "real32"},
    {std::type_index( typeid(real64)), "real64"},
    {std::type_index( typeid(localIndex)), "localIndex"},
    {std::type_index( typeid(globalIndex)), "globalIndex"},
    {std::type_index( typeid(R1Tensor)), "R1Tensor"},
    {std::type_index( typeid(R1Tensor32)), "R1Tensor32"},
    {std::type_index( typeid(R2SymTensor)), "R2SymTensor"},
    {std::type_index( typeid(integer_array)), "integer_array"},
    {std::type_index( typeid(real32_array)), "real32_array"},
    {std::type_index( typeid(real64_array)), "real64_array"},
    {std::type_index( typeid(localIndex_array)), "localIndex_array"},
    {std::type_index( typeid(globalIndex_array)), "globalIndex_array"},
    {std::type_index( typeid(integer_array2d)), "integer_array2d"},
    {std::type_index( typeid(real32_array2d)), "real32_array2d"},
    {std::type_index( typeid(real64_array2d)), "real64_array2d"},
    {std::type_index( typeid(localIndex_array2d)), "localIndex_array2d"},
    {std::type_index( typeid(globalIndex_array2d)), "globalIndex_array2d"},
    {std::type_index( typeid(integer_array3d)), "integer_array3d"},
    {std::type_index( typeid(real32_array3d)), "real32_array3d"},
    {std::type_index( typeid(real64_array3d)), "real64_array3d"},
    {std::type_index( typeid(localIndex_array3d)), "localIndex_array3d"},
    {std::type_index( typeid(globalIndex_array3d)), "globalIndex_array3d"},
    {std::type_index( typeid(real64_array4d)), "real64_array4d"},
    {std::type_index( typeid(string)), "string"},
    {std::type_index( typeid(Path)), "path"},
    {std::type_index( typeid(string_array)), "string_array"},
    {std::type_index( typeid(path_array)), "path_array"},
  };

  // If the data type is not defined here, return type_info.name()
  auto const iter = type_names.find( key );
  if( iter != type_names.end() )
  {
    return iter->second;
  }
  else
  {
    return LvArray::system::demangle( key.name());
  }
}

/**
 * @brief Build Array regexes.
 * @param subPattern
 * @param dimension
 * @return
 *
 * @note The sub pattern is the base object you are targeting.  It can either
 *       be a simple type or a lower-dimensional array. Sub-elements and
 *       axes are given as a comma-separated list enclosed in a curly brace.
 *       For example, a 2D string array would look like: {{"a", "b"}, {"c", "d"}}
 */
string constructArrayRegex( string_view subPattern, integer dimension )
{
  string subPatternStr = dimension > 1 ?
                         constructArrayRegex( subPattern, dimension-1 ) :
                         string( subPattern );
  // Add trailing space if is not already done
  if( !stringutilities::endsWith( subPatternStr, "\\s*" ) )
    subPatternStr+="\\s*";
  // Allow the bottom-level to be empty
  return dimension > 1 ?
    "\\{\\s*(" + subPatternStr + ",\\s*)*" + subPatternStr + "\\}" :
    "\\s*\\{\\s*((" + subPatternStr + "\\s*,\\s*)*" + subPatternStr + ")?\\}\\s*";
}

rtTypes::RegexMapType rtTypes::createBasicTypesRegexMap()
{
  // Define the component regexes:
  // Regex to match an unsigned int (123, etc.)
  // string_view const ru = "[\\d]+";// unused

  // Regex to match an signed int (-123, 455, +789, etc.)
  string_view const intRegex = "[+-]?[\\d]+";

  // Regex to match a float (1, +2.3, -.4, 5.6e7, 8E-9, etc.)
  // Explanation of parts:
  // [+-]?[\\d]*  matches an optional +/- at the beginning, any numbers preceding the decimal
  // ([\\d]\\.?|\\.[\\d]) matches the decimal region of the number (0, 1., 2.3, .4)
  // [\\d]*  matches any number of numbers following the decimal
  // ([eE][-+]?[\\d]+|\\s*)  matches an optional scientific notation number
  // Note: the xsd regex implementation does not allow an empty branch, so use allow whitespace at the end
  string_view const realRegex = "[+-]?[\\d]*([\\d]\\.?|\\.[\\d])[\\d]*([eE][-+]?[\\d]+|\\s*)";

  // Regex to match a R1Tensor
  string const R1Regex = "\\s*\\{\\s*(" + string( realRegex ) + ",\\s*){2}" + string( realRegex ) + "\\s*\\}";
  // Regex to match a R2SymTensor
  string const R2Regex = "\\s*\\{\\s*(" + string( realRegex ) + ",\\s*){5}" + string( realRegex ) + "\\s*\\}";

  // Regex to match a string that can't be empty and does not contain any whitespaces nor the characters ,{}
  string_view const strRegex = "[^,\\{\\}\\s]+\\s*";
  // Regex to match a string that does not contain any whitespaces nor the characters ,{}
  string_view const strRegexE = "[^,\\{\\}\\s]*\\s*";

  // Regex to match a path: a string that can't be empty and does not contain any space nor the characters *?<>|:",
  string_view const pathRegex = "[^*?<>\\|:\";,\\s]+\\s*";
  // Regex to match a path: a string that does not contain any space nor the characters *?<>|:",
  string_view const pathRegexE = "[^*?<>\\|:\";,\\s]*\\s*";

  // Regex to match a group name: it can't be empty and contains only upper/lower letters, digits, and the .-_ characters.
  string_view const groupNameRegex = "[a-zA-Z0-9.\\-_]+";
  // Regex to match an optionnal group name reference: it can be empty, contains only upper/lower letters, digits, the .-_ characters, and
  // the / character for paths.
  string_view const groupNameRefRegexE = "[a-zA-Z0-9.\\-_\\/]*";


  // Build master list of regexes
  RegexMapType regexMap =
  {
    {"integer", string( intRegex )},
    {"localIndex", string( intRegex )},
    {"globalIndex", string( intRegex )},
    {"real32", string( realRegex )},
    {"real64", string( realRegex )},
    {"R1Tensor", string( R1Regex )},
    {"R1Tensor32", string( R1Regex )},
    {"R2SymTensor", string( R2Regex )},
    {"integer_array", constructArrayRegex( intRegex, 1 )},
    {"localIndex_array", constructArrayRegex( intRegex, 1 )},
    {"globalIndex_array", constructArrayRegex( intRegex, 1 )},
    {"real32_array", constructArrayRegex( realRegex, 1 )},
    {"real64_array", constructArrayRegex( realRegex, 1 )},
    {"integer_array2d", constructArrayRegex( intRegex, 2 )},
    {"localIndex_array2d", constructArrayRegex( intRegex, 2 )},
    {"globalIndex_array2d", constructArrayRegex( intRegex, 2 )},
    {"real32_array2d", constructArrayRegex( realRegex, 2 )},
    {"real64_array2d", constructArrayRegex( realRegex, 2 )},
    {"integer_array3d", constructArrayRegex( intRegex, 3 )},
    {"localIndex_array3d", constructArrayRegex( intRegex, 3 )},
    {"globalIndex_array3d", constructArrayRegex( intRegex, 3 )},
    {"real32_array3d", constructArrayRegex( realRegex, 3 )},
    {"real64_array3d", constructArrayRegex( realRegex, 3 )},
    {"real64_array4d", constructArrayRegex( realRegex, 4 )},
    {"string", string( strRegexE )},
    {"path", string( pathRegexE )},
    {"string_array", constructArrayRegex( strRegex, 1 )},
    {"path_array", constructArrayRegex( pathRegex, 1 )},

    {string( CustomTypes::mapPair ), string( strRegexE )},
    {string( CustomTypes::plotLevel ), string( intRegex )},
    {string( CustomTypes::groupName ), string( groupNameRegex )},
    {string( CustomTypes::groupNameRef ), string( groupNameRefRegexE )},
    {string( CustomTypes::groupNameRefArray ), constructArrayRegex( groupNameRefRegexE, 1 )},
  };
  return regexMap;
}



}
