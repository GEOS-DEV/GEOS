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
 * @file ConduitRestart.cpp
 */

// Source includes
#include "ConduitRestart.hpp"
#include "mpiCommunications/MpiWrapper.hpp"
#include "common/TimingMacros.hpp"
#include "fileIO/utils/utils.hpp"

// TPL includes
#include <conduit_relay.hpp>
#include <axom/core/utilities/FileUtilities.hpp>
#include <fmt/fmt.hpp>

namespace geosx
{
namespace dataRepository
{

conduit::Node rootConduitNode;


std::string writeRootNode( std::string const & rootPath )
{
  std::string rootDirName, rootFileName;
  splitPath( rootPath, rootDirName, rootFileName );

  if( MpiWrapper::Comm_rank() == 0 )
  {
    conduit::Node root;
    root[ "number_of_files" ] = MpiWrapper::Comm_size();
    root[ "file_pattern" ] = rootFileName + "/rank_%07d.hdf5";

    conduit::relay::io::save( root, rootPath + ".root", "hdf5" );

    axom::utilities::filesystem::makeDirsForPath( rootPath );
  }

  MpiWrapper::Barrier( MPI_COMM_GEOSX );

  return fmt::sprintf( rootPath + "/rank_%07d.hdf5", MpiWrapper::Comm_rank() );
}


std::string readRootNode( std::string const & rootPath )
{
  std::string rankFilePattern;
  if( MpiWrapper::Comm_rank() == 0 )
  {
    conduit::Node node;
    conduit::relay::io::load( rootPath + ".root", "hdf5", node );

    int const nFiles = node.fetch_child( "number_of_files" ).value();
    GEOS_ERROR_IF_NE( nFiles, MpiWrapper::Comm_size() );

    std::string const filePattern = node.fetch_child( "file_pattern" ).as_string();

    std::string rootDirName, rootFileName;
    splitPath( rootPath, rootDirName, rootFileName );

    rankFilePattern = rootDirName + "/" + filePattern;
    GEOS_LOG_RANK_VAR( rankFilePattern );
  }

  MpiWrapper::Broadcast( rankFilePattern, 0 );
  return fmt::sprintf( rankFilePattern, MpiWrapper::Comm_rank() );
}

/* Write out a restart file. */
void writeTree( std::string const & path )
{
  GEOSX_MARK_FUNCTION;

  std::string const filePathForRank = writeRootNode( path );
  GEOS_LOG_RANK( "Writing out restart file at " << filePathForRank );
  conduit::relay::io::save( rootConduitNode, filePathForRank, "hdf5" );
}


void loadTree( std::string const & path )
{
  GEOSX_MARK_FUNCTION;
  std::string const filePathForRank = readRootNode( path );
  GEOS_LOG_RANK( "Reading in restart file at " << filePathForRank );
  conduit::relay::io::load( filePathForRank, "hdf5", rootConduitNode );
}

} /* end namespace dataRepository */
} /* end namespace geosx */
