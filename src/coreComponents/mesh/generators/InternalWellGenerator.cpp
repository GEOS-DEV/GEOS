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

/*
 * @file InternalWellGenerator.cpp
 *
 */

#include "InternalWellGenerator.hpp"

namespace geos
{
using namespace dataRepository;

InternalWellGenerator::InternalWellGenerator( string const & name, Group * const parent ):
  WellGeneratorBase( name, parent )
{
  enableLogLevelInput();

  registerWrapper( viewKeyStruct::polylineNodeCoordsString(), &m_polyNodeCoords ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Physical coordinates of the well polyline nodes" );

  registerWrapper( viewKeyStruct::polylineSegmentConnString(), &m_segmentToPolyNodeMap ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Connectivity of the polyline segments" );
}

void InternalWellGenerator::postProcessInput()
{
  GEOS_THROW_IF( m_polyNodeCoords.size( 1 ) != m_nDims,
                 "Invalid number of physical coordinates in " << viewKeyStruct::polylineNodeCoordsString() << " for well " << getName(),
                 InputError );

  GEOS_THROW_IF( m_segmentToPolyNodeMap.size( 1 ) != 2,
                 "Invalid size in " << viewKeyStruct::polylineSegmentConnString() << " for well " << getName(),
                 InputError );

  GEOS_THROW_IF( m_polyNodeCoords.size( 0 )-1 != m_segmentToPolyNodeMap.size( 0 ),
                 "Incompatible sizes of " << viewKeyStruct::polylineNodeCoordsString() << " and " << viewKeyStruct::polylineSegmentConnString() << " in well " << getName(),
                 InputError );

  // TODO: add more checks here
  // TODO: check that the connectivity of the well is valid
  // TODO: check that with no branching we can go from top to bottom and touch all the elements
}


REGISTER_CATALOG_ENTRY( WellGeneratorBase, InternalWellGenerator, string const &, Group * const )
}
