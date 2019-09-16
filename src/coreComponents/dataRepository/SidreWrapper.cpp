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
 * @file SidreWrapper.cpp
 */

#include "SidreWrapper.hpp"
#include "common/TimingMacros.hpp"

#include <string>
#include <cstdio>


namespace geosx
{
namespace dataRepository
{

#ifdef GEOSX_USE_ATK
using namespace axom::sidre;
#endif

SidreWrapper::SidreWrapper()
{}

SidreWrapper::~SidreWrapper()
{}

#ifdef GEOSX_USE_ATK
DataStore & SidreWrapper::dataStore()
{
  static DataStore datastore;
  return datastore;
}
#endif


/* Write out a restart file. */
void SidreWrapper::writeTree( int MPI_PARAM( num_files ),
                              const std::string & path,
                              const std::string & protocol,
                              MPI_Comm MPI_PARAM( comm ) )
{
#ifdef GEOSX_USE_ATK
  GEOSX_MARK_FUNCTION;
#ifdef GEOSX_USE_MPI
  axom::sidre::IOManager ioManager( comm );
  ioManager.write( SidreWrapper::dataStore().getRoot(), num_files, path, protocol );
#else
  SidreWrapper::dataStore().getRoot()->save( path, protocol );
#endif
#endif
}


void SidreWrapper::reconstructTree( const std::string & root_path,
                                    const std::string & protocol,
                                    MPI_Comm MPI_PARAM( comm ) )
{
#ifdef GEOSX_USE_ATK
  GEOSX_MARK_FUNCTION;
  if( !SidreWrapper::dataStore().hasAttribute( "__sizedFromParent__" ))
  {
    SidreWrapper::dataStore().createAttributeScalar( "__sizedFromParent__", -1 );
  }

#ifdef GEOSX_USE_MPI
  axom::sidre::IOManager ioManager( comm );
  ioManager.read( SidreWrapper::dataStore().getRoot(), root_path, protocol );
#else
  SidreWrapper::dataStore().getRoot()->load( root_path, protocol );
#endif
#endif
}


/* Load sidre external data. */
void SidreWrapper::loadExternalData( const std::string & root_path,
                                     MPI_Comm MPI_PARAM( comm ) )
{
#ifdef GEOSX_USE_ATK
  GEOSX_MARK_FUNCTION;
#ifdef GEOSX_USE_MPI

  axom::sidre::IOManager ioManager( comm );
  ioManager.loadExternalData( SidreWrapper::dataStore().getRoot(), root_path );
#else
  SidreWrapper::dataStore().getRoot()->loadExternalData( root_path );
#endif
#endif
}

} /* end namespace dataRepository */
} /* end namespace geosx */
