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
 * @file CellElementStencilMPFA.cpp
 */


#include "CellElementStencilMPFA.hpp"

namespace geosx
{

CellElementStencilMPFA::CellElementStencilMPFA():
  StencilBase<CellElementStencilMPFA_Traits,CellElementStencilMPFA>()
{
}


void CellElementStencilMPFA::reserve( localIndex const size )
{
  m_elementRegionIndices.reserve( size * 9 );
  m_elementSubRegionIndices.reserve( size * 9 );
  m_elementIndices.reserve( size * 9 );
  m_weights.reserve( size * 9 );
}

void CellElementStencilMPFA::add( localIndex const numPts,
                                  localIndex const * const elementRegionIndices,
                                  localIndex const * const elementSubRegionIndices,
                                  localIndex const * const elementIndices,
                                  real64 const * const weights,
                                  localIndex const connectorIndex )
{
  GEOS_ERROR_IF( numPts >= MAX_STENCIL_SIZE, "Maximum stencil size exceeded" );

  m_elementRegionIndices.appendArray( elementRegionIndices, numPts );
  m_elementSubRegionIndices.appendArray( elementSubRegionIndices, numPts );
  m_elementIndices.appendArray( elementIndices, numPts );
  m_weights.appendArray( weights, numPts );

  m_connectorIndices[connectorIndex] = m_elementRegionIndices.size()-1;
}

} /* namespace geosx */
