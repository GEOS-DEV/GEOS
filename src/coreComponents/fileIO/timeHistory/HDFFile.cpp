/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2020-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "HDFFile.hpp"

#include "common/MpiWrapper.hpp"

#include <hdf5.h>

namespace geos
{


HDFFile::HDFFile( string const & fnm, bool deleteExisting, bool parallelAccess, MPI_Comm comm ):
  m_filename( ),
  m_fileId( 0 ),
  m_faplId( 0 ),
  m_mpioFapl( parallelAccess ),
  m_comm( comm )
{
  int rnk = MpiWrapper::commRank( comm );
#ifdef GEOSX_USE_MPI
  if( m_mpioFapl )
  {
    m_faplId = H5Pcreate( H5P_FILE_ACCESS );
    H5Pset_fapl_mpio( m_faplId, m_comm, MPI_INFO_NULL );
    m_filename = fnm + ".hdf5";
  }
  else
  {
    m_faplId = H5P_DEFAULT;
    m_filename = fnm + "." + std::to_string( rnk ) + ".hdf5";
  }
#else
  m_faplId = H5P_DEFAULT;
  m_filename = fnm + ".hdf5";
#endif
  // check if file already exists
  htri_t exists = 0;
  H5E_BEGIN_TRY
  {
    exists = H5Fis_hdf5( m_filename.c_str() );
  }
  H5E_END_TRY
  // if there is an non-hdf file with the same name,
  // and we're either not using parallel access or we're rank 0
  if( exists == 0 && ( !m_mpioFapl || rnk == 0 ) )
  {
    remove( m_filename.c_str() );
  }
  if( exists > 0 && !deleteExisting )
  {
    m_fileId = H5Fopen( m_filename.c_str(), H5F_ACC_RDWR, m_faplId );
  }
  else if( exists >= 0 && deleteExisting )
  {
    m_fileId = H5Fcreate( m_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, m_faplId );
  }
  else if( exists < 0 )
  {
    m_fileId = H5Fcreate( m_filename.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, m_faplId );
  }
}

HDFFile::~HDFFile()
{
  if( m_mpioFapl )
  {
    H5Pclose( m_faplId );
  }
  H5Fflush( m_fileId, H5F_SCOPE_GLOBAL );
  H5Fclose( m_fileId );
}

bool HDFFile::hasDataset( const string & name ) const
{
  int exists = 0;
  H5E_BEGIN_TRY
  {
    exists = H5Gget_objinfo( *this, name.c_str(), 0, NULL );
  }
  H5E_END_TRY
  return ( exists == 0 );
}

}
