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
 * @file FaceElementStencil.cpp
 */

#include "FaceElementStencil.hpp"

namespace geosx
{

FaceElementStencil::FaceElementStencil():
  StencilBase< FaceElementStencil_Traits, FaceElementStencil >()
{}

void FaceElementStencil::move( LvArray::MemorySpace const space )
{
  StencilBase< FaceElementStencil_Traits, FaceElementStencil >::move( space );
  m_cellCenterToEdgeCenters.move( space, true );
}

void FaceElementStencil::add( localIndex const numPts,
                              localIndex const * const elementRegionIndices,
                              localIndex const * const elementSubRegionIndices,
                              localIndex const * const elementIndices,
                              real64 const * const weights,
                              localIndex const connectorIndex )
{
  GEOSX_ERROR_IF( numPts >= MAX_STENCIL_SIZE, "Maximum stencil size exceeded" );

  typename decltype( m_connectorIndices )::iterator iter = m_connectorIndices.find( connectorIndex );
  if( iter==m_connectorIndices.end() )
  {
    m_elementRegionIndices.appendArray( elementRegionIndices, elementRegionIndices + numPts );
    m_elementSubRegionIndices.appendArray( elementSubRegionIndices, elementSubRegionIndices + numPts );
    m_elementIndices.appendArray( elementIndices, elementIndices + numPts );
    m_weights.appendArray( weights, weights + numPts );

    m_connectorIndices[connectorIndex] = m_weights.size() - 1;
  }
  else
  {
    localIndex const stencilIndex = iter->second;
    m_elementRegionIndices.clearArray( stencilIndex );
    m_elementSubRegionIndices.clearArray( stencilIndex );
    m_elementIndices.clearArray( stencilIndex );
    m_weights.clearArray( stencilIndex );

    m_elementRegionIndices.appendToArray( stencilIndex, elementRegionIndices, elementRegionIndices + numPts );
    m_elementSubRegionIndices.appendToArray( stencilIndex, elementSubRegionIndices, elementSubRegionIndices + numPts );
    m_elementIndices.appendToArray( stencilIndex, elementIndices, elementIndices + numPts );
    m_weights.appendToArray( stencilIndex, weights, weights + numPts );
  }
}

void FaceElementStencil::add( localIndex const numPts,
                              R1Tensor const * const cellCenterToEdgeCenter,
                              localIndex const connectorIndex )
{
  GEOSX_ERROR_IF( numPts >= MAX_STENCIL_SIZE, "Maximum stencil size exceeded" );

  typename decltype( m_connectorIndices )::iterator iter = m_connectorIndices.find( connectorIndex );
  if( iter==m_connectorIndices.end() )
  {
    GEOSX_ERROR( "Wrong connectorIndex" );
  }
  else
  {
    localIndex const stencilIndex = iter->second;
    if( stencilIndex < m_cellCenterToEdgeCenters.size())
    {
      m_cellCenterToEdgeCenters.clearArray( stencilIndex );
      m_cellCenterToEdgeCenters.appendToArray( stencilIndex, cellCenterToEdgeCenter, cellCenterToEdgeCenter + numPts );
    }
    else
    {
      m_cellCenterToEdgeCenters.appendArray( cellCenterToEdgeCenter, cellCenterToEdgeCenter + numPts );
    }
  }
}

void FaceElementStencilWrapper::computeTransmissibility( localIndex iconn,
                                                         real64 transmissibility )
{
  localIndex const er0  =  m_elementRegionIndices[iconn][0];
  localIndex const esr0 =  m_elementSubRegionIndices[iconn][0];
  localIndex const ei0  =  m_elementIndices[iconn][0];

  localIndex const er1  =  m_elementRegionIndices[iconn][1];
  localIndex const esr1 =  m_elementSubRegionIndices[iconn][1];
  localIndex const ei1  =  m_elementIndices[iconn][1];

  real64 const t0 = m_weights[iconn][0] * m_permeability[er0][esr0][ei0];
  real64 const t1 = m_weights[iconn][1] * m_permeability[er1][esr1][ei1];

  real64 const harmonicWeight   = t0*t1 / (t0+t1);
  real64 const arithmeticWeight = (t0+t1)/2;

  transmissibility = meanPermCoeff * harmonicWeight + (1 - meanPermCoeff) * arithmeticWeight;
}

void FaceElementStencilWrapper::dTrans_dPressure( localIndex iconn,
                                                  real64 (&dTrans_dPressure )[2] )
{
  localIndex const er0  =  m_elementRegionIndices[iconn][0];
  localIndex const esr0 =  m_elementSubRegionIndices[iconn][0];
  localIndex const ei0  =  m_elementIndices[iconn][0];

  localIndex const er1  =  m_elementRegionIndices[iconn][1];
  localIndex const esr1 =  m_elementSubRegionIndices[iconn][1];
  localIndex const ei1  =  m_elementIndices[iconn][1];

  real64 const dt0 = m_weights[iconn][0] * m_dPerm_dPressure[er0][esr0][ei0];
  real64 const dt1 = m_weights[iconn][1] * m_dPerm_dPressure[er1][esr1][ei1];

  dTrans_dPressure[0];
  dTrans_dPressure[1];
}




} /* namespace geosx */
