/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CPMeshParser.cpp
 */

#include "CPMeshParser.hpp"

#include "meshUtilities/CPMesh/utilities/StringUtilities.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include <iostream>

namespace geosx
{

namespace CPMesh
{

CPMeshParser::CPMeshParser( string const & name )
  : m_meshName( name )
{}


void CPMeshParser::readNumberOfCells( Path const & filePath,
                                      localIndex & nX,
                                      localIndex & nY,
                                      localIndex & nZ ) const
{
  // Note: This is definitely not what we want to do (reading the full file into a string)
  //       For now, to get a quick start, I ported this from PAMELA. I will improve that in a second step
  string fileContent = CPMeshStringUtilities::fileToString( filePath );
  std::istringstream meshFile;
  meshFile.str( fileContent );

  nX = 0;
  nY = 0;
  nZ = 0;

  // first read the dimensions of the mesh
  std::string line, buffer;
  while( getline( meshFile, line ) )
  {
    CPMeshStringUtilities::removeStringAndFollowingContentFromLine( "--", line );
    CPMeshStringUtilities::removeExtraSpaces( line );
    CPMeshStringUtilities::removeEndOfLine( line );
    CPMeshStringUtilities::removeTab( line );
    CPMeshStringUtilities::trim( line );
    if( line == "SPECGRID" || line == "DIMENS" )
    {
      std::vector< localIndex > bufInt;
      buffer = extractDataBelowKeyword( meshFile );
      CPMeshStringUtilities::fromStringTo( buffer, bufInt );

      nX = bufInt[0];
      nY = bufInt[1];
      nZ = bufInt[2];
    }
  }
  GEOSX_THROW_IF( nX <= 0 || nY <= 0 || nZ <= 0,
                  "Neither SPECGRID nor DIMENS was found",
                  InputError );

}

void CPMeshParser::readMesh( Path const & filePath,
                             CPMeshDimensions const & cPMeshDims )
{
  // Note: This is definitely not what we want to do (reading the full file into a string)
  //       For now, to get a quick start, I ported this from PAMELA. I will improve that in a second step23
  string fileContent = CPMeshStringUtilities::fileToString( filePath );
  std::istringstream meshFile;
  meshFile.str( fileContent );

  bool foundCOORD = false;
  bool foundZCORN = false;
  bool foundACTNUM = false;
  bool foundPERMX = false;
  bool foundPERMY = false;
  bool foundPERMZ = false;

  std::string line;
  while( getline( meshFile, line ) )
  {
    CPMeshStringUtilities::removeStringAndFollowingContentFromLine( "--", line );
    CPMeshStringUtilities::removeExtraSpaces( line );
    CPMeshStringUtilities::removeEndOfLine( line );
    CPMeshStringUtilities::removeTab( line );
    CPMeshStringUtilities::trim( line );

    // TODO: shorten, there must be a way to factorize all these functions
    //       at least, we shoud be able to merge readLocalACTNUM and readLocalPROP

    // at this point, SPECGRID has already been read, we skip it
    if( line == "COORD" )
    {
      foundCOORD = true;
      readLocalCOORD( meshFile, cPMeshDims );
    }
    else if( line == "ZCORN" )
    {
      foundZCORN = true;
      readLocalZCORN( meshFile, cPMeshDims );
    }
    else if( line == "ACTNUM" )
    {
      foundACTNUM = true;
      readLocalACTNUM( meshFile, cPMeshDims );
    }
    else if( line == "PERMX" )
    {
      foundPERMX = true;
      readLocalPROP( meshFile, cPMeshDims, m_permx );
    }
    else if( line == "PERMY" )
    {
      foundPERMY = true;
      readLocalPROP( meshFile, cPMeshDims, m_permy );
    }
    else if( line == "PERMZ" )
    {
      foundPERMZ = true;
      readLocalPROP( meshFile, cPMeshDims, m_permz );
    }
    else if( line == "PORO" )
    {
      readLocalPROP( meshFile, cPMeshDims, m_poro );
    }
  }

  if( !foundACTNUM )
  {
    fillLocalACTNUM( cPMeshDims );
  }

  // TODO: error out if one of them is missing
  if( foundPERMX && foundPERMY && foundPERMZ )
  {
    m_perm.resizeDimension< 0, 1 >( m_permx.size(), 3 );
    for( localIndex i = 0; i < m_permx.size(); ++i )
    {
      m_perm( i, 0 ) = m_permx( i );
      m_perm( i, 1 ) = m_permy( i );
      m_perm( i, 2 ) = m_permz( i );
    }
  }

  GEOSX_THROW_IF( !foundCOORD || !foundZCORN,
                  "At least one of the following keywords was not found: COORD, ZCORN",
                  InputError );


}

void CPMeshParser::readLocalCOORD( std::istringstream & meshFile,
                                   CPMeshDimensions const & cPMeshDims )
{
  localIndex const nX = cPMeshDims.nX();
  localIndex const nY = cPMeshDims.nY();

  localIndex const maxSize = 6*(nX+1)*(nY+1); // TODO: change
  std::vector< real64 > work;
  work.reserve( maxSize );

  string buffer = extractDataBelowKeyword( meshFile ); // TODO: change to save only what is necessary
  CPMeshStringUtilities::eclipseDataBufferToVector( buffer, work ); // TODO: change to save only what is necessary

  localIndex const nXLocal = cPMeshDims.nXLocal();
  localIndex const nYLocal = cPMeshDims.nYLocal();

  localIndex const iMinLocal = cPMeshDims.iMinLocal();
  localIndex const jMinLocal = cPMeshDims.jMinLocal();

  m_coord.resize( 6*(nXLocal+1)*(nYLocal+1) );

  // Note: Below I mimic what I will ultimately do: not storing the full file in a string,
  //       but instead going through the file and saving only what the MPI rank needs
  for( localIndex j = 0; j < nYLocal+1; ++j )
  {
    for( localIndex i = 0; i < nXLocal+1; ++i )
    {
      localIndex const indPillarLocal = j*(nXLocal+1)+i;
      localIndex const indPillar = (j+jMinLocal)*(nX+1)+(i+iMinLocal);
      for( localIndex pos = 0; pos < 6; ++pos )
      {
        m_coord( 6*indPillarLocal+pos ) = work[ 6*indPillar+pos ]; // TODO: change to avoid doing that
      }
    }
  }
}

void CPMeshParser::readLocalZCORN( std::istringstream & meshFile,
                                   CPMeshDimensions const & cPMeshDims )
{
  localIndex const nX = cPMeshDims.nX();
  localIndex const nY = cPMeshDims.nY();
  localIndex const nZ = cPMeshDims.nZ();

  localIndex const maxSize = 8*nX*nY*nZ; // TODO: change
  std::vector< real64 > work;
  work.reserve( maxSize );

  string buffer = extractDataBelowKeyword( meshFile ); // TODO: change to save only what is necessary
  CPMeshStringUtilities::eclipseDataBufferToVector( buffer, work ); // TODO: change to save only what is necessary

  localIndex const nXLocal = cPMeshDims.nXLocal();
  localIndex const nYLocal = cPMeshDims.nYLocal();
  localIndex const nZLocal = cPMeshDims.nZLocal();

  localIndex const iMinLocal = cPMeshDims.iMinLocal();
  localIndex const jMinLocal = cPMeshDims.jMinLocal();
  localIndex const kMinLocal = 0;

  m_zcorn.resize( 8*nXLocal*nYLocal*nZLocal );

  // Note: Below I mimic what I will ultimately do: not storing the full file in a string,
  //       but instead going through the file and saving only what the MPI rank needs
  for( localIndex k = 0; k < nZLocal; ++k )
  {
    for( localIndex j = 0; j < nYLocal; ++j )
    {
      for( localIndex i = 0; i < nXLocal; ++i )
      {
        localIndex const i_xmLocal = k*8*nXLocal*nYLocal + j*4*nXLocal + 2*i;
        localIndex const i_xpLocal = i_xmLocal + 2*nXLocal;
        localIndex const i_xm = (k+kMinLocal)*8*nX*nY + (j+jMinLocal)*4*nX + 2*(i+iMinLocal);
        localIndex const i_xp = i_xm + 2*nX;

        m_zcorn( i_xmLocal ) = work[ i_xm ];
        m_zcorn( i_xmLocal + 1 ) = work[ i_xm + 1 ];
        m_zcorn( i_xpLocal ) = work[ i_xp ];
        m_zcorn( i_xpLocal + 1 ) = work[ i_xp + 1 ];
        m_zcorn( i_xmLocal + 4*nXLocal*nYLocal ) = work[ i_xm + 4*nX*nY ];
        m_zcorn( i_xmLocal + 4*nXLocal*nYLocal + 1 ) = work[ i_xm + 4*nX*nY + 1 ];
        m_zcorn( i_xpLocal + 4*nXLocal*nYLocal ) = work[ i_xp + 4*nX*nY ];
        m_zcorn( i_xpLocal + 4*nXLocal*nYLocal + 1 ) = work[ i_xp + 4*nX*nY + 1 ];
      }
    }
  }
}

void CPMeshParser::readLocalACTNUM( std::istringstream & meshFile,
                                    CPMeshDimensions const & cPMeshDims )
{
  localIndex const nX = cPMeshDims.nX();
  localIndex const nY = cPMeshDims.nY();
  localIndex const nZ = cPMeshDims.nZ();

  localIndex const maxSize = nX*nY*nZ; // TODO: change
  std::vector< real64 > work;
  work.reserve( maxSize );

  string buffer = extractDataBelowKeyword( meshFile ); // TODO: change to save only what is necessary
  CPMeshStringUtilities::eclipseDataBufferToVector( buffer, work ); // TODO: change to save only what is necessary
  std::replace( work.begin(), work.end(), 2, 0 );
  std::replace( work.begin(), work.end(), 3, 0 );

  localIndex const nXLocal = cPMeshDims.nXLocal();
  localIndex const nYLocal = cPMeshDims.nYLocal();
  localIndex const nZLocal = cPMeshDims.nZLocal();

  localIndex const iMinLocal = cPMeshDims.iMinLocal();
  localIndex const jMinLocal = cPMeshDims.jMinLocal();
  localIndex const kMinLocal = 0;

  m_actnum.resize( nXLocal*nYLocal*nZLocal );

  // Note: Below I mimic what I will ultimately do: not storing the full file in a string,
  //       but instead going through the file and saving only what the MPI rank needs
  for( localIndex k = 0; k < nZLocal; ++k )
  {
    for( localIndex j = 0; j < nYLocal; ++j )
    {
      for( localIndex i = 0; i < nXLocal; ++i )
      {
        localIndex const eiLocal = k*nXLocal*nYLocal + j*nXLocal + i;
        localIndex const ei = (k+kMinLocal)*nX*nY + (j+jMinLocal)*nX + i+iMinLocal;
        m_actnum( eiLocal ) = work[ ei ];
      }
    }
  }
}

void CPMeshParser::readLocalPROP( std::istringstream & meshFile,
                                  CPMeshDimensions const & cPMeshDims,
                                  array1d< real64 > & prop )
{
  localIndex const nX = cPMeshDims.nX();
  localIndex const nY = cPMeshDims.nY();
  localIndex const nZ = cPMeshDims.nZ();

  localIndex const maxSize = nX*nY*nZ; // TODO: change
  std::vector< real64 > work;
  work.reserve( maxSize );

  string buffer = extractDataBelowKeyword( meshFile ); // TODO: change to save only what is necessary
  CPMeshStringUtilities::eclipseDataBufferToVector( buffer, work ); // TODO: change to save only what is necessary

  localIndex const nXLocal = cPMeshDims.nXLocal();
  localIndex const nYLocal = cPMeshDims.nYLocal();
  localIndex const nZLocal = cPMeshDims.nZLocal();

  localIndex const iMinLocal = cPMeshDims.iMinLocal();
  localIndex const jMinLocal = cPMeshDims.jMinLocal();
  localIndex const kMinLocal = 0;

  prop.resize( nXLocal*nYLocal*nZLocal );

  // Note: Below I mimic what I will ultimately do: not storing the full file in a string,
  //       but instead going through the file and saving only what the MPI rank needs

  for( localIndex k = 0; k < nZLocal; ++k )
  {
    for( localIndex j = 0; j < nYLocal; ++j )
    {
      for( localIndex i = 0; i < nXLocal; ++i )
      {
        localIndex const eiLocal = k*nXLocal*nYLocal + j*nXLocal + i;
        localIndex const ei = (k+kMinLocal)*nX*nY + (j+jMinLocal)*nX + i+iMinLocal;
        prop( eiLocal ) = work[ ei ];
      }
    }
  }
}

void CPMeshParser::fillLocalACTNUM( CPMeshDimensions const & cPMeshDims )
{
  localIndex const nXLocal = cPMeshDims.nXLocal();
  localIndex const nYLocal = cPMeshDims.nYLocal();
  localIndex const nZLocal = cPMeshDims.nZLocal();

  m_actnum.resize( nXLocal*nYLocal*nZLocal );
  m_actnum.setValues< serialPolicy >( 1 );
}

std::string CPMeshParser::extractDataBelowKeyword( std::istringstream & stringBlock ) const
{
  char keywordEnd = '/';
  std::string chunk;
  std::streampos pos;
  std::vector< std::string > res;
  getline( stringBlock, chunk, keywordEnd );
  stringBlock.clear();

  CPMeshStringUtilities::removeStringAndFollowingContentFromLine( "--", chunk );
  CPMeshStringUtilities::removeTab( chunk );
  CPMeshStringUtilities::removeEndOfLine( chunk );
  CPMeshStringUtilities::trim( chunk );
  res.push_back( chunk );
  stringBlock.ignore( 10, '\n' );
  getline( stringBlock, chunk );
  CPMeshStringUtilities::removeTab( chunk );
  CPMeshStringUtilities::removeEndOfLine( chunk );
  return res[0];
}

REGISTER_CATALOG_ENTRY( CPMeshParser, CPMeshParser, string const & )

} // namespace CPMesh

} // end namespace geosx
