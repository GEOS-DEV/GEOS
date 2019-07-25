/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file CellElementStencilTPFA.cpp
 */


#include "CellElementStencilTPFA.hpp"
#include "codingUtilities/Utilities.hpp"

namespace geosx
{

CellElementStencilTPFA::CellElementStencilTPFA():
  StencilBase<CellElementStencilTPFA_Traits,CellElementStencilTPFA>()
{
}


void CellElementStencilTPFA::add( localIndex const numPts,
                                  localIndex const * const elementRegionIndices,
                                  localIndex const * const elementSubRegionIndices,
                                  localIndex const * const elementIndices,
                                  real64 const * const weights,
                                  localIndex const connectorIndex )
{
  GEOS_ERROR_IF( numPts!=2, "number of cells in TPFA stencil should be 2");

  localIndex const oldSize = m_elementRegionIndices.size(0);
  localIndex const newSize = oldSize + 1;
  m_elementRegionIndices.resize( newSize, numPts );
  m_elementSubRegionIndices.resize( newSize, numPts );
  m_elementIndices.resize( newSize, numPts );
  m_weights.resize( newSize, numPts );

  for( localIndex a=0 ; a<numPts ; ++a )
  {
    m_elementRegionIndices(oldSize,a) = elementRegionIndices[a];
    m_elementSubRegionIndices(oldSize,a) = elementSubRegionIndices[a];
    m_elementIndices(oldSize,a) = elementIndices[a];
    m_weights(oldSize,a) = weights[a];
  }
  m_connectorIndices[connectorIndex] = oldSize;
}

} /* namespace geosx */
