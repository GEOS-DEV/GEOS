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
 * @file InputParsing.hpp
 */

#ifndef GEOS_DATAREPOSITORY_INPUTPARSING_HPP_
#define GEOS_DATAREPOSITORY_INPUTPARSING_HPP_

#include "xmlWrapper.hpp"
#include "KeyNames.hpp"
// #include "InputExtensionGroup.hpp"

namespace geos
{

namespace inputParsing
{

using input_document_type = xmlWrapper::xmlDocument;

/**
 * @brief constexpr variable to hold name for inserting the file path into the xml file.
 *
 * This is used because we would like the option to hold the file path in the xml structure.
 * The name is uglified with underscores to avoid collisions with real attribute names.
 */
constexpr char const filePathString[] = "__filePath__";

/**
 * @brief constexpr variable to hold node character offset from the start of the xml file.
 *
 * This is used because we would like the option to hold the offset in the xml structure.
 * The name is uglified with underscores to avoid collisions with real attribute names.
 */
constexpr char const charOffsetString[] = "__charOffset__";

inline string mergeInputDocuments( string_array const & inputFileList, string const & outputDir = {} )
{
  // TODO: handle different document types
  return xmlWrapper::mergeInputDocuments( inputFileList, outputDir );
}

template < typename Document >
void processIncludes( Document & document, int level = 0 )
{
  if constexpr( std::is_same_v< Document, xmlWrapper::xmlDocument > )
  {
    document.processIncludes( document.getFirstChild(), level );
  }
}

template < typename Document >
bool isDocMetadataAttribute( string const & attributeName )
{
  if constexpr( std::is_same_v< Document, xmlWrapper::xmlDocument > )
  {
    return xmlWrapper::isFileMetadataAttribute( attributeName );
  }
  return false;
}

template < typename Document, typename T, typename U >
bool readAttributeAsType( T & rval, string const & name, Regex const & regex, typename Document::node_type const & docNode, U const & lastArg )
{
  if constexpr( std::is_same_v< Document, xmlWrapper::xmlDocument > )
  {
    return xmlWrapper::readAttributeAsType( rval, name, regex, docNode, lastArg );
  }
  return false;
}

template < typename Document >
void processInputException( std::exception const & ex, string const & wrapperName, typename Document::node_type const & docNode, typename Document::node_pos_type const & docNodePos )
{
  if constexpr( std::is_same_v< Document, xmlWrapper::xmlDocument > )
  {
    return xmlWrapper::processInputException( ex, wrapperName, docNode, docNodePos );
  }
}

} // namespace inputParsing

} // namespace geos

#endif