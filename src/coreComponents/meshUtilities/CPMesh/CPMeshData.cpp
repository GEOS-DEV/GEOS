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
 * @file CPMeshData.cpp
 */

#include "CPMeshData.hpp"

#include "mpiCommunications/MpiWrapper.hpp"

namespace geosx
{

namespace CPMesh
{

CPMeshData::CPMeshData( string const & name )
  :
  m_nX( 0 ),
  m_nY( 0 ),
  m_nZ( 0 ),
  m_iMin( 0 ),
  m_jMin( 0 ),
  m_iMax( 0 ),
  m_jMax( 0 ),
  m_iMinOverlap( 0 ),
  m_jMinOverlap( 0 ),
  m_iMaxOverlap( 0 ),
  m_jMaxOverlap( 0 ),
  m_nOwnedVertices( 0 ),
  m_nLocalVertices( 0 ),
  m_nOwnedActiveCells( 0 ),
  m_nLocalActiveCells( 0 ),
  m_meshName( name )
{}

void CPMeshData::defineDomainBoundaries( localIndex const nX, localIndex const nY, localIndex const nZ )
{
  m_nX = nX;
  m_nY = nY;
  m_nZ = nZ;
}

void CPMeshData::definePartitionBoundaries( localIndex const iMin, localIndex const jMin,
                                            localIndex const iMax, localIndex const jMax )
{
  m_iMin = iMin;
  m_jMin = jMin;
  m_iMax = iMax;
  m_jMax = jMax;
}

void CPMeshData::definePartitionOverlaps( localIndex const iMinOverlap, localIndex const jMinOverlap,
                                          localIndex const iMaxOverlap, localIndex const jMaxOverlap )
{
  m_iMinOverlap = iMinOverlap;
  m_jMinOverlap = jMinOverlap;
  m_iMaxOverlap = iMaxOverlap;
  m_jMaxOverlap = jMaxOverlap;
}


REGISTER_CATALOG_ENTRY( CPMeshData, CPMeshData, string const & )

} // namespace CPMesh

} // end namespace geosx
