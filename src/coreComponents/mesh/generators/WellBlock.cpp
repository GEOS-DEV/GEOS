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
 * @file DomainPartition.cpp
 */

#include "mesh/Perforation.hpp"
#include "mesh/generators/WellBlock.hpp"
#include "mesh/generators/InternalWellGenerator.hpp"




namespace geosx
{

WellBlock::WellBlock( string const & name, Group * const parent ):
  WellBlockABC( name, parent )
{}

void WellBlock::setPrevElemIndices( arrayView1d< arrayView1d< globalIndex const > const > prevElemIndices) {
  int size = prevElemIndices.size();
  m_prevElemId.resize( size );
  for ( int i = 0; i < size; i++ )
  {
    m_prevElemId[i] = prevElemIndices[i];
  }
}
}
