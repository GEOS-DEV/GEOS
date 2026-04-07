/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
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
  registerWrapper( viewKeyStruct::polylineNodeCoordsString(), &m_polyNodeCoords ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Physical coordinates of the well polyline nodes" );

  registerWrapper( viewKeyStruct::polylineSegmentConnString(), &m_segmentToPolyNodeMap ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Connectivity of the polyline segments" );
}

void InternalWellGenerator::postInputInitialization()
{
  GEOS_THROW_IF( m_polyNodeCoords.size( 1 ) != m_nDims,
                 GEOS_FMT( "InternalWell {}: Invalid number of physical coordinates.",
                           getWrapperDataContext( viewKeyStruct::polylineNodeCoordsString() ) ),
                 InputError, getWrapperDataContext( viewKeyStruct::polylineNodeCoordsString() ) );

  GEOS_THROW_IF( m_segmentToPolyNodeMap.size( 1 ) != 2,
                 GEOS_FMT( "InternalWell {}: Invalid size.",
                           getWrapperDataContext( viewKeyStruct::polylineSegmentConnString() ) ),
                 InputError, getWrapperDataContext( viewKeyStruct::polylineSegmentConnString() ) );

  GEOS_THROW_IF( m_polyNodeCoords.size( 0 )-1 != m_segmentToPolyNodeMap.size( 0 ),
                 GEOS_FMT( "Incompatible sizes of {} and {}",
                           getWrapperDataContext( viewKeyStruct::polylineNodeCoordsString() ),
                           getWrapperDataContext( viewKeyStruct::polylineSegmentConnString() ) ),
                 InputError, getWrapperDataContext( viewKeyStruct::polylineSegmentConnString() ) );

  // TODO: add more checks here
  // TODO: check that the connectivity of the well is valid
  // TODO: check that with no branching we can go from top to bottom and touch all the elements
}


REGISTER_CATALOG_ENTRY( MeshComponentBase, InternalWellGenerator, string const &, Group * const )
}
