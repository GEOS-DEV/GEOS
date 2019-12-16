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
 * @file xmlWrapper.hpp
 */

#ifndef GEOSX_DATAREPOSITORY_XMLWRAPPER_HPP_
#define GEOSX_DATAREPOSITORY_XMLWRAPPER_HPP_

// Source includes
#include "common/DataTypes.hpp"
#include "dataRepository/DefaultValue.hpp"
#include "cxx-utilities/src/ArrayUtilities.hpp"

// TPL includes
#include <pugixml.hpp>

// System includes
#include <algorithm>
#include <sstream>

namespace geosx
{

namespace dataRepository
{
class Group;
}


/**
 * @class xmlWrapper
 *
 * This class wraps/provides facilities to process entries from an xml file into the appropriate
 * data types. xmlWrapper provides some aliases that will wrap the underlying xml package being
 * used to extract data from the xml file, and a set of static member functions that facilitate
 * the parsing of string data from the xml into the variables that will hold those values in the
 * code.
 */
class xmlWrapper
{
public:

  /// Alias for the type of xml document.
  using xmlDocument = pugi::xml_document;

  /// Alias for the type of the result from an xml parse attempt.
  using xmlResult = pugi::xml_parse_result;

  /// Alias for the type of an xml node.
  using xmlNode = pugi::xml_node;

  /// Alias for the type of an xml attribute.
  using xmlAttribute = pugi::xml_attribute;

  /// Alias for the type variant of an xml node.
  using xmlTypes = pugi::xml_node_type;

  /**
   * @brief constexpr variable to hold name for inserting the file path into the xml file.
   *
   * This is used because we would like the option to hold the file path in the xml structure.
   */
  static constexpr auto filePathString = "filePath";

  /**
   * @name Constructors/destructor/assignment.
   */
  ///@{
  xmlWrapper() = delete;
  ~xmlWrapper() = delete;
  xmlWrapper( xmlWrapper const & ) = delete;
  xmlWrapper( xmlWrapper && ) = delete;
  xmlWrapper & operator=( xmlWrapper const & ) = delete;
  xmlWrapper & operator=( xmlWrapper && ) = delete;
  ///@}

  /**
   * @brief Function to add xml nodes from included files.
   * @param targetNode the node for which to look for included children specifications
   *
   * This function looks for a "Included" node under the targetNode, loops over all subnodes under the "Included"
   * node, and then parses the file specified in those subnodes taking all the nodes in the file and adding them to
   * the targetNode.
   */
  static void addIncludedXML( xmlNode & targetNode );

  /**
   * @name String to variable parsing.
   *
   * These functions take in @p value and parse that string based on the type of
   * @p target. The function implementation should provide sufficient error checking
   * in the case that @p value is formatted incorrectly for the type specified in @p target.
   */
  ///@{

  /**
   * @brief Parse a string and fill a variable with the value(s) in the string.
   * @tparam T the type of variable fill with string value
   * @param[out] target the object to read values into
   * @param[in]  value  the string that contains the data to be parsed into target
   */
  template< typename T >
  static void StringToInputVariable( T & target, string value );

  /**
   * @copybrief StringToInputVariable(T &, string)
   * @param[out] target the object to read values into
   * @param[in]  value  the string that contains the data to be parsed into target
   */
  static void StringToInputVariable( R1Tensor & target, string value );

  /**
   * @copybrief StringToInputVariable(T &, string)
   * @tparam T    data type of the array
   * @tparam NDIM number of dimensions of the array
   * @param[out] array the array to read values into
   * @param[in]  value the string that contains the data to be parsed into target
   */
  template< typename T, int NDIM >
  static void StringToInputVariable( Array< T, NDIM > & array, string value );

  ///@}


  /**
   * @name Attribute extraction from XML nodes.
   */
  ///@{

  /**
   * @brief Extract attribute in an xml tree, and translate its value into a typed variable.
   * @tparam T             the type of variable fill with xml attribute.
   * @param[out] rval      the variable to fill with value
   * @param[in] name       the name of the xml attribute to process
   * @param[in] targetNode the xml node that should contain the attribute
   * @param[in] required   whether or not the value is required
   * @return boolean value indicating whether the value was successfully read from XML.
   */
  template< typename T >
  static bool ReadAttributeAsType( T & rval,
                                   string const & name,
                                   xmlNode const & targetNode,
                                   bool const required );

  /**
   * @copybrief ReadAttributeAsType(T &, string const &, xmlNode const &, bool const);
   * @tparam T             the type of variable fill with xml attribute.
   * @tparam T_DEF         type of the default value for @p rval
   * @param[out] rval      the variable to fill with value
   * @param[in] name       the name of the xml attribute to process
   * @param[in] targetNode the xml node that should contain the attribute
   * @param[in] defVal     default value of @p rval (or entries of @p rval, if it is an array)
   * @return boolean value indicating whether the value was successfully read from XML.
   */
  template< typename T, typename T_DEF = T >
  static bool ReadAttributeAsType( T & rval,
                                   string const & name,
                                   xmlNode const & targetNode,
                                   T_DEF const & defVal );

  /**
   * @copybrief ReadAttributeAsType(T &, string const &, xmlNode const &, bool const);
   * @tparam T             the type of variable fill with xml attribute.
   * @param[out] rval      the variable to fill with value
   * @param[in] name       the name of the xml attribute to process
   * @param[in] targetNode the xml node that should contain the attribute
   * @param[in] defVal     default value of @p rval (or entries of @p rval, if it is an array)
   * @return boolean value indicating whether the value was successfully read from XML.
   */
  template< typename T >
  static std::enable_if_t< !dataRepository::DefaultValue< T >::has_default_value, bool >
  ReadAttributeAsType( T & rval,
                       string const & name,
                       xmlNode const & targetNode,
                       dataRepository::DefaultValue< T > const & defVal )
  {
    GEOSX_UNUSED_VAR( defVal );
    return ReadAttributeAsType( rval, name, targetNode, false );
  }

  /**
   * @copybrief ReadAttributeAsType(T &, string const &, xmlNode const &, bool const);
   * @tparam T             the type of variable fill with xml attribute.
   * @param[out] rval      the variable to fill with value
   * @param[in] name       the name of the xml attribute to process
   * @param[in] targetNode the xml node that should contain the attribute
   * @param[in] defVal     default value of @p rval (or entries of @p rval, if it is an array)
   * @return boolean value indicating whether the value was successfully read from XML.
   */
  template< typename T >
  static typename std::enable_if_t< dataRepository::DefaultValue< T >::has_default_value, bool >
  ReadAttributeAsType( T & rval,
                       string const & name,
                       xmlNode const & targetNode,
                       dataRepository::DefaultValue< T > const & defVal )
  {
    return ReadAttributeAsType( rval, name, targetNode, defVal.value );
  }

  ///@}

};


template< typename T >
void xmlWrapper::StringToInputVariable( T & target, string inputValue )
{
  std::istringstream ss( inputValue );
  ss>>target;
}

template< typename T, int NDIM >
void xmlWrapper::StringToInputVariable( Array< T, NDIM > & array, string valueString )
{
  cxx_utilities::stringToArray( array, valueString );
}

template< typename T >
bool xmlWrapper::ReadAttributeAsType( T & rval,
                                      string const & name,
                                      xmlNode const & targetNode,
                                      bool const required )
{
  pugi::xml_attribute xmlatt = targetNode.attribute( name.c_str() );

  bool const success = !(xmlatt.empty() && required);
  //GEOSX_ERROR_IF( xmlatt.empty() && required, "Input variable " + name + " is required in " + targetNode.path() );

  if( success )
  {
    // parse the string/attribute into a value
    StringToInputVariable( rval, xmlatt.value() );
  }
  return success;
}


template< typename T, typename T_DEF >
bool xmlWrapper::ReadAttributeAsType( T & rval,
                                      string const & name,
                                      xmlNode const & targetNode,
                                      T_DEF const & defVal )
{
  pugi::xml_attribute xmlatt = targetNode.attribute( name.c_str() );
  if( !xmlatt.empty() )
  {
    // parse the string/attribute into a value
    StringToInputVariable( rval, xmlatt.value() );
  }
  else
  {
    // set the value to the default value
    rval = defVal;
  }
  return true;
}



} /* namespace geosx */

#endif /*GEOSX_DATAREPOSITORY_XMLWRAPPER_HPP_ */
