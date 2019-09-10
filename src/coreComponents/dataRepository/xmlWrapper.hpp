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

/**
 * @file xmlWrapper.hpp
 */

#ifndef _FILEIO_XMLWRAPPER_HPP_
#define _FILEIO_XMLWRAPPER_HPP_

#include <algorithm>
#include <sstream>

#include "pugixml.hpp"

#include "common/DataTypes.hpp"
#include "dataRepository/DefaultValue.hpp"
#include "ArrayUtilities.hpp"

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
  /// @typedef alias for the type of xml document
  using xmlDocument = pugi::xml_document;

  /// @typedef alias for the type of the result from an xml parse attempt
  using xmlResult = pugi::xml_parse_result;

  /// @typedef alias for the type of an xml node
  using xmlNode = pugi::xml_node;

  /// @typedef alias for the type of an xml attribute
  using xmlAttribute = pugi::xml_attribute;

  /// @typedef alias for the type variant of an xml node
  using xmlTypes = pugi::xml_node_type;

  /// constexpr variable to hold name for inserting the file path into the xml file. This is used
  /// because we would like the option to hold the file path in the xml structure.
  static constexpr auto filePathString = "filePath";

  /**
   * @name FUNCTION GROUP for rule of five functions...which are all deleted in this case.
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
   * Function to add xml nodes from included files.
   * @param targetNode the node for which to look for included children specifications
   *
   * This function looks for a "Included" node under the targetNode, loops over all subnodes under the "Included"
   * node, and then parses the file specified in those subnodes taking all the nodes in the file and adding them to
   * the targetNode.
   */
  static void addIncludedXML( xmlNode & targetNode );

  /**
   * @name FUNCTION GROUP for StringToInputVariable()
   * @brief Functions to parse a string, and fill a variable with the value/s in the string.
   * @tparam T the type of variable fill with string value
   * @param[out] target the object to read values into
   * @param[in] value the string that contains the data to be parsed into target
   *
   * This function takes in @p value and parses that string based on the type of
   * @p target. The function implementation should provide sufficient error
   * checking in the case that @p value is formatted incorrectly for the type specified in
   * @p target.
   */
  ///@{
  template< typename T >
  static void StringToInputVariable( T & target, string value );

  static void StringToInputVariable( R1Tensor & target, string value );

  template< typename T, int NDIM >
  static void StringToInputVariable( LvArray::Array< T, NDIM, localIndex > & array, string value );
  ///@}


  /**
   * @name FUNCTION GROUP for ReadAttributeAsType()
   * @brief Functions to extract attributes in an xml tree, and translate those values into a
   *        typed variable.
   * @tparam T the type of variable fill with xml attribute.
   * @tparam T_DEF the default value of @p T, or in the case where @p T is an array,
   *                the entries of T.
   * @param rval the variable to fill.
   * @param name the name of the xml attribute to process
   * @param targetNode The xml node that should contain the attribute
   * @param required whether or not the value is required
   */
  ///@{

  template< typename T >
  static void ReadAttributeAsType( T & rval,
                                   string const & name,
                                   xmlNode const & targetNode,
                                   bool const required );

  template< typename T, typename T_DEF = T >
  static void ReadAttributeAsType( T & rval,
                                   string const & name,
                                   xmlNode const & targetNode,
                                   T_DEF const & defVal );

  template< typename T >
  static typename std::enable_if_t< !dataRepository::DefaultValue< T >::has_default_value >
  ReadAttributeAsType( T & rval,
                       string const & name,
                       xmlNode const & targetNode,
                       dataRepository::DefaultValue< T > const & defVal )
  {
    ReadAttributeAsType( rval, name, targetNode, false );
  }

  template< typename T >
  static typename std::enable_if_t< dataRepository::DefaultValue< T >::has_default_value >
  ReadAttributeAsType( T & rval,
                       string const & name,
                       xmlNode const & targetNode,
                       dataRepository::DefaultValue< T > const & defVal )
  {
    ReadAttributeAsType( rval, name, targetNode, defVal.value );
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
void xmlWrapper::StringToInputVariable( LvArray::Array< T, NDIM, localIndex > & array, string valueString )
{
  cxx_utilities::stringToArray( array, valueString );
}

template< typename T >
void xmlWrapper::ReadAttributeAsType( T & rval,
                                      string const & name,
                                      xmlNode const & targetNode,
                                      bool const required )
{
  pugi::xml_attribute xmlatt = targetNode.attribute( name.c_str() );

  GEOS_ERROR_IF( xmlatt.empty() && required, "Input variable " + name + " is required in " + targetNode.path() );

  // parse the string/attribute into a value
  StringToInputVariable( rval, xmlatt.value() );
}


template< typename T, typename T_DEF >
void xmlWrapper::ReadAttributeAsType( T & rval,
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
}



} /* namespace geosx */

#endif /*_FILEIO_XMLWRAPPER_HPP_ */
