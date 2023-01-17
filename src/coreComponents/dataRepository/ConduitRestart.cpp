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
 * @file ConduitRestart.cpp
 */

// Source includes
#include "ConduitRestart.hpp"
#include "common/MpiWrapper.hpp"
#include "common/TimingMacros.hpp"
#include "common/Path.hpp"

// TPL includes
#include <conduit_relay.hpp>

namespace geosx
{
namespace dataRepository
{

string writeRootFile( conduit::Node & root, string const & rootPath )
{
  string const completeRootPath = rootPath;
  string const rootFileName = splitPath( completeRootPath ).second;

  if( MpiWrapper::commRank() == 0 )
  {
    makeDirsForPath( completeRootPath );

    root[ "protocol/name" ] = "hdf5";
    root[ "protocol/version" ] = CONDUIT_VERSION;

    root[ "number_of_files" ] = MpiWrapper::commSize();
    root[ "file_pattern" ] = rootFileName + "/rank_%07d.hdf5";

    root[ "number_of_trees" ] = 1;
    root[ "tree_pattern" ] = "/";

    conduit::relay::io::save( root, completeRootPath + ".root", "hdf5" );
  }

  MpiWrapper::barrier( MPI_COMM_GEOSX );
  return GEOSX_FMT( "{}/rank_{:07}.hdf5", completeRootPath.data(), MpiWrapper::commRank() );
}


string readRootNode( string const & rootPath )
{
  string rankFilePattern;
  if( MpiWrapper::commRank() == 0 )
  {
    conduit::Node node;
    conduit::relay::io::load( rootPath + ".root", "hdf5", node );

    int const nFiles = node.fetch_existing( "number_of_files" ).value();
    GEOSX_THROW_IF_NE( nFiles, MpiWrapper::commSize(), InputError );

    string const filePattern = node.fetch_existing( "file_pattern" ).as_string();
    string const rootDirName = splitPath( rootPath ).first;

    rankFilePattern = rootDirName + "/" + filePattern;
    GEOSX_LOG_RANK_VAR( rankFilePattern );
  }

  MpiWrapper::broadcast( rankFilePattern, 0 );

  char buffer[ 1024 ];
  GEOSX_ERROR_IF_GE( std::snprintf( buffer, 1024, rankFilePattern.data(), MpiWrapper::commRank() ), 1024 );
  return buffer;
}

void writeTree( string const & path, conduit::Node & root )
{
  GEOSX_MARK_FUNCTION;

  conduit::Node rootFileNode;
  string const filePathForRank = writeRootFile( rootFileNode, path );
  GEOSX_LOG_RANK( "Writing out restart file at " << filePathForRank );
  conduit::relay::io::save( root, filePathForRank, "hdf5" );
}

void loadTree( string const & path, conduit::Node & root )
{
  GEOSX_MARK_FUNCTION;
  string const filePathForRank = readRootNode( path );
  GEOSX_LOG_RANK( "Reading in restart file at " << filePathForRank );
  conduit::relay::io::load( filePathForRank, "hdf5", root );
}

} /* end namespace dataRepository */
} /* end namespace geosx */
