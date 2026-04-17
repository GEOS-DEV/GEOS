/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ArchiveInputDeck.cpp
 */

#include "ArchiveInputDeck.hpp"

#include "common/MpiWrapper.hpp"
#include "common/Path.hpp"
#include "common/format/Format.hpp"
#include "common/initializeEnvironment.hpp"
#include "common/logger/Logger.hpp"
#include "dataRepository/Group.hpp"
#include "dataRepository/xmlWrapper.hpp"
#include "mainInterface/ProblemManager.hpp"

#include <conduit.hpp>

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <system_error>

namespace geos
{

using namespace dataRepository;

namespace archiveInputDeck
{

namespace
{

string makeTimestamp()
{
  auto const now = std::chrono::system_clock::now();
  auto const time_t_now = std::chrono::system_clock::to_time_t( now );
  std::ostringstream timestampStream;
  timestampStream << std::put_time( std::localtime( &time_t_now ), "%Y%m%d_%H%M%S" );
  return timestampStream.str();
}

void stripMetadataAttributes( xmlWrapper::xmlNode node )
{
  node.remove_attribute( xmlWrapper::filePathString );
  node.remove_attribute( xmlWrapper::charOffsetString );

  for( xmlWrapper::xmlNode child : node.children() )
  {
    stripMetadataAttributes( child );
  }
}

void reorderTags( xmlWrapper::xmlNode rootNode, string_array const & tagOrder )
{
  xmlWrapper::xmlNode lastInserted;
  for( string const & tagName : tagOrder )
  {
    xmlWrapper::xmlNode tag = rootNode.child( tagName.c_str() );
    if( !tag )
    {
      continue;
    }

    lastInserted ? rootNode.insert_move_after( tag, lastInserted )
                  : rootNode.append_move( tag );

    lastInserted = tag;
  }

  // ProblemManager's order list doesn't provide every XML tags available in GEOS
  // so we put the missing ones below the ones it provides.
  // And sort them alphabetically
  stdVector< string > missingTags;

  for( xmlWrapper::xmlNode const & tag : rootNode.children() )
  {
    string const & tagName = tag.name();

    if( std::find( tagOrder.begin(), tagOrder.end(), tag.name() ) == tagOrder.end() )
    {
      missingTags.push_back( tagName );
    }
  }

  std::sort( missingTags.begin(), missingTags.end() );

  for( string const & tagName : missingTags )
  {
    xmlWrapper::xmlNode tag = rootNode.child( tagName.c_str() );

    if( tag )
    {
      rootNode.append_move( tag );
    }
  }
}

void sortAttributes( xmlWrapper::xmlNode node )
{
  stdVector< std::pair< string, string > > attributes;
  for( xmlWrapper::xmlAttribute attr = node.first_attribute();
       attr;
       attr = attr.next_attribute() )
  {
    attributes.emplace_back( attr.name(), attr.value() );
  }

  std::sort( attributes.begin(),
             attributes.end(),
             []( std::pair< string, string > const & a,
                 std::pair< string, string > const & b )
  {
    // name attribute should be the first attribute, and not sorted alphabetically
    bool const aIsName = ( a.first == "name" );
    bool const bIsName = ( b.first == "name" );
    if( aIsName != bIsName )
    {
      return aIsName;
    }

    // other attributes are sorted alphabetically
    return a.first < b.first;
  } );

  // pugi doesn't have any move_attribute method yet, so we have to
  // copy and remove attributes
  while( node.remove_attribute( node.first_attribute() ) )
  {}
  for( auto const & attr : attributes )
  {
    node.append_attribute( attr.first.c_str() ).set_value( attr.second.c_str() );
  }

  for( xmlWrapper::xmlNode child : node.children() )
  {
    sortAttributes( child );
  }
}

xmlWrapper::xmlDocument flattenXMLs( string_array const & fileNames )
{
  xmlWrapper::xmlDocument flatDoc;
  xmlWrapper::xmlNode root = flatDoc.appendChild( "Problem" );

  for( string const & fileName : fileNames )
  {
    xmlWrapper::xmlDocument doc;
    xmlWrapper::xmlResult const result = doc.loadFile( fileName, true );
    GEOS_THROW_IF( !result,
                   GEOS_FMT( "Could not load XML file '{}': {}", fileName, result.description() ),
                   InputError );
    xmlWrapper::xmlNode docRoot = doc.getFirstChild();

    doc.addIncludedXML( docRoot );

    for( xmlWrapper::xmlNode & node : docRoot.children() )
    {
      root.append_copy( node );
    }
  }

  return flatDoc;
}


}


void archiveInputDeck( CommandLineOptions const & opts )
{
  if( opts.inputFileNames.empty() || opts.outputDirectory.empty() )
  {
    return;
  }

  if( MpiWrapper::commRank() != 0 )
  {
    return;
  }

  conduit::Node tempRoot;
  ProblemManager tempPM( tempRoot );

  string_array xmlTagOrder;
  tempPM.initializationOrder( xmlTagOrder );

  string const timestamp = makeTimestamp();
  string const archiveDir = joinPath( opts.outputDirectory, "archive_inputFiles", timestamp );
  makeDirsForPath( archiveDir + "/" );

  xmlWrapper::xmlDocument flatDoc = flattenXMLs( opts.inputFileNames );
  xmlWrapper::xmlNode root = flatDoc.getFirstChild();

  stripMetadataAttributes( root );
  reorderTags( root, xmlTagOrder );
  sortAttributes( root );

  flatDoc.saveFile( joinPath( archiveDir, "input.xml" ) );
}


} /* namespace archiveInputDeck */

} /* namespace geos */
