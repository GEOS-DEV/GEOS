/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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

/*
 * @file SegmentManager.cpp
 *
 */

#include "SegmentManager.hpp"
#include "Segment.hpp"

#include "managers/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"

namespace geosx
{

using namespace dataRepository;

SegmentManager::SegmentManager(string const & name, ManagedGroup * const parent)
  : ObjectManagerBase(name, parent)
{
  RegisterViewWrapper( viewKeyStruct::segElementRegionString, &m_segElementRegion, false );
  RegisterViewWrapper( viewKeyStruct::segElementSubregionString, &m_segElementSubregion, false );
  RegisterViewWrapper( viewKeyStruct::segElementIndexString, &m_segElementIndex, false );
  RegisterViewWrapper( viewKeyStruct::segIndexString, &m_segIndex, false );

  RegisterViewWrapper( viewKeyStruct::gravityDepthString, &m_gravityDepth, false );
}

SegmentManager::~SegmentManager()
{

}

ManagedGroup * SegmentManager::CreateChild(string const & childKey, string const & childName)
{
  if ( childKey == groupKeyStruct::segmentString )
  {
    m_segmentList.push_back( childName );
    return RegisterGroup<Segment>( childName );
  }
  else
  {
    GEOS_ERROR( "Unrecognized node: " << childKey );
  }
  return nullptr;
}

const string SegmentManager::getCatalogName() const
{
  return keys::segments;
}

void SegmentManager::InitializePreSubGroups( ManagedGroup * const problemManager )
{

}

void SegmentManager::InitializePostInitialConditions_PreSubGroups( ManagedGroup * const problemManager )
{
  DomainPartition const * domain = problemManager->GetGroup<DomainPartition>( keys::domain );
  MeshLevel const * mesh = domain->getMeshBody(0)->getMeshLevel(0);
  PrecomputeData( mesh );
}

void SegmentManager::PrecomputeData( MeshLevel const * mesh )
{

}

Segment const * SegmentManager::getSegment( localIndex iseg ) const
{
  return this->GetGroup<Segment>( m_segmentList[iseg] );
}

Segment * SegmentManager::getSegment( localIndex iseg )
{
  return this->GetGroup<Segment>( m_segmentList[iseg] );
}


} //namespace geosx
