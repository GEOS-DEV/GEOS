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
 * @file xmlWrapper.cpp
 */

#include "xmlWrapper.hpp"

#include "codingUtilities/StringUtilities.hpp"
#include "common/MpiWrapper.hpp"
#include "dataRepository/KeyNames.hpp"

namespace geosx
{
using namespace dataRepository;

namespace xmlWrapper
{

template< typename T, int SIZE >
void stringToInputVariable( Tensor< T, SIZE > & target, string const & inputValue )
{
  std::istringstream ss( inputValue );
  auto const errorMsg = [&]( auto const & msg )
  {
    return GEOSX_FMT( "{} for Tensor<{}> at position {} in input: {}", msg, SIZE, ss.tellg(), inputValue );
  };

  // Read the head
  ss >> std::ws;
  GEOSX_THROW_IF_NE_MSG( ss.peek(), '{', errorMsg( "Missing opening brace" ), InputError );
  ss.ignore(); // skip the opening brace

  // Read the values
  T value;
  int count = 0;
  while( ss >> value )
  {
    GEOSX_THROW_IF_GE_MSG( count, SIZE, errorMsg( "Too many values" ), InputError );
    target[count++] = value;
    ss >> std::ws;
    if( count < SIZE )
    {
      GEOSX_THROW_IF_NE_MSG( ss.peek(), ',', errorMsg( ss.peek() == '}' ? "Not enough values" : "Missing comma separator" ), InputError );
      ss.ignore(); // skip the comma
    }
  }
  GEOSX_THROW_IF_LT_MSG( count, SIZE, errorMsg( "Not enough values" ), InputError );
  ss.clear();

  // Read the tail
  GEOSX_THROW_IF_NE_MSG( ss.peek(), '}', errorMsg( "Missing closing brace" ), InputError );
  ss.ignore(); // skip the closing brace
  ss >> std::ws;
  GEOSX_THROW_IF( ss.peek() != std::char_traits< char >::eof(), errorMsg( "Unparsed characters" ), InputError );
}

template void stringToInputVariable( Tensor< real64, 3 > & target, string const & inputValue );
template void stringToInputVariable( Tensor< real64, 6 > & target, string const & inputValue );

void addIncludedXML( xmlNode & targetNode, int const level )
{
  GEOSX_THROW_IF( level > 100, "XML include level limit reached, please check input for include loops", InputError );

  xmlNode const rootNode = targetNode.root();
  string const currentFilePath = rootNode.child( filePathString ).attribute( filePathString ).value();

  // Schema currently allows a single unique <Included>, but a non-validating file may include multiple
  for( xmlNode includedNode : targetNode.children( includedListTag ) )
  {
    for( xmlNode fileNode : includedNode.children() )
    {
      // Extract the file name and construct full includedDirPath
      string const includedFilePath = [&]()
      {
        GEOSX_THROW_IF_NE_MSG( string( fileNode.name() ), includedFileTag,
                               GEOSX_FMT( "<{}> must only contain <{}> tags", includedListTag, includedFileTag ),
                               InputError );
        xmlAttribute const nameAttr = fileNode.attribute( "name" );
        string const fileName = nameAttr.value();
        GEOSX_THROW_IF( !nameAttr || fileName.empty(),
                        GEOSX_FMT( "<{}> tag must have a non-empty 'name' attribute", includedFileTag ),
                        InputError );
        return isAbsolutePath( fileName ) ? fileName : joinPath( splitPath( currentFilePath ).first, fileName );
      }();

      xmlDocument includedXmlDocument;
      xmlResult const result = includedXmlDocument.load_file( includedFilePath.c_str() );
      GEOSX_THROW_IF( !result, GEOSX_FMT( "Errors found while parsing included XML file {}\nDescription: {}\nOffset: {}",
                                          includedFilePath, result.description(), result.offset ), InputError );

      // All included files must contain a root node that must match the target node.
      // Currently, schema only allows <Included> tags at the top level (inside <Problem>).
      // We then proceed to merge each nested node from included file with the one in main.

      xmlNode includedRootNode = includedXmlDocument.first_child();
      GEOSX_THROW_IF_NE_MSG( string( includedRootNode.name() ), string( targetNode.name() ),
                             "Included document root does not match the including XML node", InputError );

      // Process potential includes in the included file to allow nesting
      includedXmlDocument.append_child( filePathString ).append_attribute( filePathString ).set_value( includedFilePath.c_str() );
      addIncludedXML( includedRootNode, level + 1 );

      // Add each top level tag of imported document to current
      // This may result in repeated XML blocks, which will be implicitly merged when processed
      for( xmlNode importedNode : includedRootNode.children() )
      {
        // Save path to current file on the node, in order to handle relative paths later
        if( !importedNode.attribute( filePathString ) )
        {
          importedNode.append_attribute( filePathString ).set_value( includedFilePath.c_str() );
        }
        targetNode.append_copy( importedNode );
      }
    }
  }

  // Just in case, remove <Included> tags that have been processed
  while( targetNode.remove_child( includedListTag ) )
  {}
}

string buildMultipleInputXML( string_array const & inputFileList,
                              string const & outputDir )
{
  if( inputFileList.empty() )
  {
    return {};
  }
  if( inputFileList.size() == 1 )
  {
    return inputFileList[0];
  }

  string inputFileName = joinPath( outputDir, "composite_input.xml" );

  // Write the composite xml file on one rank
  if( MpiWrapper::commRank() == 0 )
  {
    xmlWrapper::xmlDocument compositeTree;
    xmlWrapper::xmlNode compositeRoot = compositeTree.append_child( dataRepository::keys::ProblemManager );
    xmlWrapper::xmlNode includedRoot = compositeRoot.append_child( includedListTag );

    for( auto & fileName: inputFileList )
    {
      xmlWrapper::xmlNode fileNode = includedRoot.append_child( includedFileTag );
      fileNode.append_attribute( "name" ) = fileName.c_str();
    }

    compositeTree.save_file( inputFileName.c_str() );
  }

  // Everybody else has to wait before attempting to read
  MpiWrapper::barrier();

  return inputFileName;
}

} /* namespace xmlWrapper */

} /* namespace geosx */
