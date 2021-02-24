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
 * @file CellElementStencilTPFA.cpp
 */


#include "CellElementStencilTPFA.hpp"
#include "codingUtilities/Utilities.hpp"

namespace geosx
{

CellElementStencilTPFA::CellElementStencilTPFA():
  StencilBase< CellElementStencilTPFA_Traits, CellElementStencilTPFA >()
{}


void CellElementStencilTPFA::add( localIndex const numPts,
                                  localIndex const * const elementRegionIndices,
                                  localIndex const * const elementSubRegionIndices,
                                  localIndex const * const elementIndices,
                                  real64 const * const weights,
                                  localIndex const connectorIndex )
{
  GEOSX_ERROR_IF_NE_MSG( numPts, 2, "Number of cells in TPFA stencil should be 2" );

  localIndex const oldSize = m_elementRegionIndices.size( 0 );
  localIndex const newSize = oldSize + 1;
  m_elementRegionIndices.resize( newSize, numPts );
  m_elementSubRegionIndices.resize( newSize, numPts );
  m_elementIndices.resize( newSize, numPts );
  m_weights.resize( newSize, numPts );

  for( localIndex a=0; a<numPts; ++a )
  {
    m_elementRegionIndices( oldSize, a ) = elementRegionIndices[a];
    m_elementSubRegionIndices( oldSize, a ) = elementSubRegionIndices[a];
    m_elementIndices( oldSize, a ) = elementIndices[a];
    m_weights( oldSize, a ) = weights[a];
  }
  m_connectorIndices[connectorIndex] = oldSize;
}

  WeightContainerViewConstType & CellElementStencilTPFA::
  computeEffectiveWeights( CoefficientAccessor< arrayView2d< real64 const > > const & coefficient ) const
  {
    WeightContainerType effectiveWeights;

    for( localIndex iconn = 0; iconn < size(); iconn++ )
    {
      localIndex const er0  =  m_elementRegionIndices[iconn][0];
      localIndex const esr0 =  m_elementSubRegionIndices[iconn][0];
      localIndex const ei0  =  m_elementIndices[iconn][0];

      localIndex const er1  =  m_elementRegionIndices[iconn][1];
      localIndex const esr1 =  m_elementSubRegionIndices[iconn][1];
      localIndex const ei1  =  m_elementIndices[iconn][1];

      // need the normal direction here.
      real64 const t0 = m_weights[iconn][0] * coefficient[er0][esr0][ei0] * aperture[er0][esr0][ei0];
      real64 const t1 = m_weights[iconn][1] * coefficient[er1][esr1][ei1] * aperture[er1][esr1][ei1];

      real64 const harmonicWeight   = t0*t1 / (t0+t1);
      real64 const arithmeticWeight = (t0+t1)/2;

      real64 const c = 1.00;

      real64 const weight = c * harmonicWeight + (1 - c) * arithmeticWeight;

      effectiveWeights[iconn][0] = weight;
      effectiveWeights[iconn][1] = weight;
    }

    return effectiveWeights;
  }

} /* namespace geosx */
