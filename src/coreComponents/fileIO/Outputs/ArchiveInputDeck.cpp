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

#include "common/Path.hpp"
#include "common/format/Format.hpp"
#include "common/logger/Logger.hpp"
#include "dataRepository/xmlWrapper.hpp"

#include <algorithm>
#include <chrono>
#include <filesystem>

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

    }

    return prefix + relPath;
}

}

void archiveInputDeck( string_array const & inputFileNames,
                       string const & outputDirectory )
{
  if ( inputFileNames.empty() )
  {
    return;
  }


string archiveInputDeck( string_array const & inputFileNames,
                         string const & outputDirectory,
                         string_array const & xmlTagOrder )
{
  if( inputFileNames.empty() || outputDirectory.empty() )
  {
    return {};
  }

  string const timestamp = makeTimestamp();
  string const archiveDir = joinPath( outputDirectory, "archive_inputFiles", timestamp );
  makeDirsForPath( archiveDir + "/" );

  xmlWrapper::xmlDocument flatDoc;
  xmlWrapper::xmlNode root = flatDoc.appendChild( "Problem" );

  for( string const & fileName : inputFileNames )
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

  stripMetadataAttributes( root );
  reorderTags( root, xmlTagOrder );

  flatDoc.saveFile( joinPath( archiveDir, "input.xml" ) );

  return archiveDir;
}

} /* namespace archiveInputDeck */

} /* namespace geos */
