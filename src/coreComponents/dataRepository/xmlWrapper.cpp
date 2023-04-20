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

#include <regex>

#include "xmlWrapper.hpp"

#include "codingUtilities/StringUtilities.hpp"
#include "common/MpiWrapper.hpp"
#include "dataRepository/KeyNames.hpp"

namespace geos
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
    return GEOS_FMT( "{} for Tensor<{}> at position {} in input: {}", msg, SIZE, ss.tellg(), inputValue );
  };

  // Read the head
  ss >> std::ws;
  GEOS_THROW_IF_NE_MSG( ss.peek(), '{', errorMsg( "Missing opening brace" ), InputError );
  ss.ignore(); // skip the opening brace

  // Read the values
  T value;
  int count = 0;
  while( ss >> value )
  {
    GEOS_THROW_IF_GE_MSG( count, SIZE, errorMsg( "Too many values" ), InputError );
    target[count++] = value;
    ss >> std::ws;
    if( count < SIZE )
    {
      GEOS_THROW_IF_NE_MSG( ss.peek(), ',', errorMsg( ss.peek() == '}' ? "Not enough values" : "Missing comma separator" ), InputError );
      ss.ignore(); // skip the comma
    }
  }
  GEOS_THROW_IF_LT_MSG( count, SIZE, errorMsg( "Not enough values" ), InputError );
  ss.clear();

  // Read the tail
  GEOS_THROW_IF_NE_MSG( ss.peek(), '}', errorMsg( "Missing closing brace" ), InputError );
  ss.ignore(); // skip the closing brace
  ss >> std::ws;
  GEOS_THROW_IF( ss.peek() != std::char_traits< char >::eof(), errorMsg( "Unparsed characters" ), InputError );
}

template void stringToInputVariable( Tensor< real32, 3 > & target, string const & inputValue );
template void stringToInputVariable( Tensor< real64, 3 > & target, string const & inputValue );
template void stringToInputVariable( Tensor< real64, 6 > & target, string const & inputValue );

/**
 * @brief Adds the filePath and character offset infos on the node in filePathString
 * and charOffsetString attributes. This function allow to keep track of the source
 * filename & offset of each node.
 * @param node the target node to add the informations on.
 * @param filePath the absolute path of the xml file containing the node.
 */
void addNodeFileInfos( xmlNode node, string const & filePath )
{
  // we keep the file path and the character offset on each node so we keep track of these
  // informations, even if the nodes are manipulated within the xml hierarchy.
  node.append_attribute( filePathString ).set_value( filePath.c_str() );
  node.append_attribute( charOffsetString ).set_value( node.offset_debug() );

  for( xmlNode subNode : node.children() )
  {
    addNodeFileInfos( subNode, filePath );
  }
}
/**
 * @brief Returns true if the addNodeFileInfos() command has been called of the specified node.
 */
bool xmlDocument::hasNodeFileInfos() const
{ return !first_child().attribute( filePathString ).empty(); }

void xmlDocument::addIncludedXML( xmlNode & targetNode, int const level )
{
  GEOS_THROW_IF( level > 100, "XML include level limit reached, please check input for include loops", InputError );

  string const currentFilePath = targetNode.attribute( filePathString ).value();

  // Schema currently allows a single unique <Included>, but a non-validating file may include multiple
  for( xmlNode includedNode : targetNode.children( includedListTag ) )
  {
    for( xmlNode fileNode : includedNode.children() )
    {
      // Extract the file name and construct full includedDirPath
      string const includedFilePath = [&]()
      {
        GEOS_THROW_IF_NE_MSG( string( fileNode.name() ), includedFileTag,
                              GEOS_FMT( "<{}> must only contain <{}> tags", includedListTag, includedFileTag ),
                              InputError );
        xmlAttribute const nameAttr = fileNode.attribute( "name" );
        string const fileName = nameAttr.value();
        GEOS_THROW_IF( !nameAttr || fileName.empty(),
                       GEOS_FMT( "<{}> tag must have a non-empty 'name' attribute", includedFileTag ),
                       InputError );
        return isAbsolutePath( fileName ) ? fileName : joinPath( splitPath( currentFilePath ).first, fileName );
      }();

      xmlDocument includedXmlDocument;
      xmlResult const result = includedXmlDocument.load_file( includedFilePath.c_str(),
                                                              hasNodeFileInfos() );
      GEOS_THROW_IF( !result, GEOS_FMT( "Errors found while parsing included XML file {}\n"
                                        "Description: {}\nOffset: {}",
                                        includedFilePath, result.description(), result.offset ),
                     InputError );


      // All included files must contain a root node that must match the target node.
      // Currently, schema only allows <Included> tags at the top level (inside <Problem>).
      // We then proceed to merge each nested node from included file with the one in main.

      xmlNode includedRootNode = includedXmlDocument.first_child();
      GEOS_THROW_IF_NE_MSG( string( includedRootNode.name() ), string( targetNode.name() ),
                            "Included document root does not match the including XML node", InputError );

      // Process potential includes in the included file to allow nesting
      addIncludedXML( includedRootNode, level + 1 );

      // Add each top level tag of imported document to current
      // This may result in repeated XML blocks, which will be implicitly merged when processed
      for( xmlNode importedNode : includedRootNode.children() )
      {
        targetNode.append_copy( importedNode );
      }

      m_originalBuffers[includedXmlDocument.getFilePath()] = includedXmlDocument.getOriginalBuffer();
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

bool isFileMetadataAttribute( string const & name )
{
  static const std::set< string > fileMetadataAttributes {
    "name", "xmlns:xsi", "xsi:noNamespaceSchemaLocation", xmlWrapper::filePathString, xmlWrapper::charOffsetString
  };
  return fileMetadataAttributes.find( name ) != fileMetadataAttributes.end();
}

const size_t xmlDocument::npos = string::npos;
size_t documentId=0;

xmlDocument::xmlDocument():
  pugi::xml_document(),
  m_rootFilePath( "CodeIncludedXML" + std::to_string( documentId++ ) )
{}

xmlResult xmlDocument::load_string( const pugi::char_t * contents, bool loadNodeFileInfos,
                                    unsigned int options )
{
  xmlResult result = pugi::xml_document::load_string( contents, options );

  // keeping a copy of original buffer to allow line retrieval
  if( loadNodeFileInfos )
  {
    new (&m_originalBuffers) map< string, string >();
    m_originalBuffers[m_rootFilePath] = string( contents );

    addNodeFileInfos( first_child(), m_rootFilePath );
  }

  return result;
}
xmlResult xmlDocument::load_file( const char * path, bool loadNodeFileInfos,
                                  unsigned int options, pugi::xml_encoding encoding )
{
  xmlResult result = pugi::xml_document::load_file( path, options, encoding );
  m_rootFilePath = getAbsolutePath( path );

  // keeping a copy of original buffer to allow line retrieval
  if( loadNodeFileInfos )
  {
    std::ifstream t( path );
    std::stringstream buffer;
    buffer << t.rdbuf();

    new (&m_originalBuffers) map< string, string >();
    m_originalBuffers[m_rootFilePath] = string( buffer.str() );

    addNodeFileInfos( first_child(), getAbsolutePath( m_rootFilePath ) );
  }

  return result;
}
xmlResult xmlDocument::load_buffer( const void * contents, size_t size, bool loadNodeFileInfos,
                                    unsigned int options, pugi::xml_encoding encoding )
{
  xmlResult result = pugi::xml_document::load_buffer( contents, size, options, encoding );

  //keeping a copy of original buffer
  if( loadNodeFileInfos )
  {
    new (&m_originalBuffers) map< string, string >();
    m_originalBuffers[m_rootFilePath] = string( ( char const * )contents, size );

    addNodeFileInfos( first_child(), m_rootFilePath );
  }

  return result;
}

string const & xmlDocument::getFilePath() const
{ return m_rootFilePath; }

string const & xmlDocument::getOriginalBuffer() const
{ return m_originalBuffers.find( m_rootFilePath )->second; }

string const * xmlDocument::getOriginalBuffer( string const & filePath ) const
{
  map< string, string >::const_iterator it = m_originalBuffers.find( filePath );
  return it != m_originalBuffers.cend() ? &it->second : nullptr;
}

map< string, string > const & xmlDocument::getOriginalBuffers() const
{ return m_originalBuffers; }

xmlNodePos xmlDocument::getNodePosition( xmlNode const & node ) const
{
  size_t line = npos;
  size_t offsetInLine = npos;
  size_t offset = npos;
  xmlAttribute filePathAtt = node.attribute( filePathString );
  xmlAttribute charOffsetAtt = node.attribute( charOffsetString );
  string filePath;

  if( filePathAtt && charOffsetAtt )
  {
    filePath = string( filePathAtt.value() );
    offset = std::atoi( charOffsetAtt.value() );

    auto sourceBuffer = m_originalBuffers.find( filePath );
    string nodeName = node.name();

    if( sourceBuffer!=m_originalBuffers.end() &&
        offset > 0 && offset + nodeName.size() < sourceBuffer->second.size() &&
        sourceBuffer->second.substr( offset, nodeName.size() ) == nodeName )
    {
      line = std::count( sourceBuffer->second.begin(), sourceBuffer->second.begin() + offset, '\n' ) + 1;
      offsetInLine = offset - sourceBuffer->second.rfind( '\n', offset );
    }
    else
    {
      std::ostringstream oss;
      oss << "Node from file=" << filePath << ", path=" << node.path() << ", offset=" << offset;
      oss << " could not be found: ";
      oss << ( sourceBuffer == m_originalBuffers.end() ?
               ( filePathAtt.empty() ? "Source file not found." : "Source file infos not loaded") :
               ( offset > 0 && offset + nodeName.size() < sourceBuffer->second.size() ?
                 "The offset doesn't lead to this node." :
                 "Offset out of bounds." )
               );
      oss << "\nSource files list:";
      for( auto const & buffer : m_originalBuffers )
        oss << "\n  " << buffer.first << " buffer size = " << std::to_string( buffer.second.size() );

      throw InputError( oss.str() );
    }
  }

  return xmlNodePos( *this, filePath, line, offsetInLine, offset );
}

xmlNodePos::xmlNodePos( xmlDocument const & document_, string const & filePath_, size_t line_,
                        size_t offsetInLine_, size_t offset_ ):
  document( document_ ),
  filePath( filePath_ ),
  line( line_ ),
  offsetInLine( offsetInLine_ ),
  offset( offset_ )
{}

bool xmlNodePos::isFound() const
{ return line != xmlDocument::npos; }

xmlAttributePos xmlNodePos::getAttributeLine( string const & attName ) const
{
  string const * buffer = document.getOriginalBuffer( filePath );
  size_t tagEnd = xmlDocument::npos;
  size_t attOffset = xmlDocument::npos;
  size_t attLine = xmlDocument::npos;
  size_t attOffsetInLine = xmlDocument::npos;

  if( isFound() && buffer != nullptr && offset < buffer->size() )
  {
    tagEnd = buffer->find( '>', offset );
    if( tagEnd != string::npos )
    {
      std::smatch m;
      // we search for a string which is the attribute name followed by an '=', eventually separated by spaces
      if( std::regex_search( buffer->cbegin() + offset, buffer->cbegin() + tagEnd,
                             m, std::regex( attName + "\\s*=" ) ) )
      {
        attOffset = m.position() + offset;
        attLine = line + std::count( buffer->cbegin() + offset, buffer->cbegin() + attOffset, '\n' );
        attOffsetInLine = attOffset - buffer->rfind( '\n', attOffset );
      }
    }
  }

  return xmlAttributePos( filePath, attLine, attOffsetInLine, attOffset );
}

xmlAttributePos::xmlAttributePos( string const & filePath_, size_t line_, size_t offsetInLine_,
                                  size_t offset_ ):
  filePath( filePath_ ),
  line( line_ ),
  offsetInLine( offsetInLine_ ),
  offset( offset_ )
{}

bool xmlAttributePos::isFound() const
{ return offset != xmlDocument::npos; }

} /* namespace xmlWrapper */

} /* namespace geos */
