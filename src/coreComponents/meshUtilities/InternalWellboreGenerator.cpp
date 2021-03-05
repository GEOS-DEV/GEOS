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
 * @file InternalWellboreGenerator.cpp
 */

#include "InternalWellboreGenerator.hpp"
#include "managers/DomainPartition.hpp"
#include "mpiCommunications/PartitionBase.hpp"

namespace geosx
{
using namespace dataRepository;

InternalWellboreGenerator::InternalWellboreGenerator( string const & name, Group * const parent ):
  InternalMeshGenerator( name, parent )
{}

void InternalWellboreGenerator::generateMesh( DomainPartition & domain )
{
  InternalMeshGenerator::generateMesh( domain );

  Group & meshBodies = domain.getGroup( string( "MeshBodies" ));
  MeshBody & meshBody = meshBodies.registerGroup< MeshBody >( this->getName() );
  MeshLevel & meshLevel0 = meshBody.registerGroup< MeshLevel >( string( "Level0" ));
  NodeManager & nodeManager = meshLevel0.getNodeManager();
  Group & nodeSets = nodeManager.sets();

  // Wellbore nodesets
  SortedArray< localIndex > & rnegNodes = nodeSets.registerWrapper< SortedArray< localIndex > >( string( "rneg" ) ).reference();
  SortedArray< localIndex > & rposNodes = nodeSets.registerWrapper< SortedArray< localIndex > >( string( "rpos" ) ).reference();
  SortedArray< localIndex > & tnegNodes = nodeSets.registerWrapper< SortedArray< localIndex > >( string( "tneg" ) ).reference();
  SortedArray< localIndex > & tposNodes = nodeSets.registerWrapper< SortedArray< localIndex > >( string( "tpos" ) ).reference();

  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();
  for( int localNodeIndex=0; localNodeIndex<nodeManager.size(); ++localNodeIndex )
  {
    real64 const xCoord = X( localNodeIndex, 0 );
    real64 const yCoord = X( localNodeIndex, 1 );
    real64 rCoord = sqrt( xCoord * xCoord + yCoord * yCoord );

    if( isEqual( rCoord, m_min[0], m_positionTolerance ) )
    {
      rnegNodes.insert( localNodeIndex );
    }

    if( isEqual( rCoord, m_max[0], m_positionTolerance ) )
    {
      rposNodes.insert( localNodeIndex );
    }

    real64 tCoord;

    if( yCoord>=0 )
    {
      tCoord = acos( xCoord/rCoord );
    }
    else
    {
      tCoord = 2*M_PI - acos( xCoord/rCoord );
    }

    tCoord *= 180/M_PI;

    if( isEqual( tCoord, m_min[1], m_positionTolerance ) )
    {
      tnegNodes.insert( localNodeIndex );
    }
    if( isEqual( tCoord, m_max[1], m_positionTolerance ) )
    {
      tposNodes.insert( localNodeIndex );
    }
  }
}

REGISTER_CATALOG_ENTRY( MeshGeneratorBase, InternalWellboreGenerator, string const &, Group * const )
} /* namespace geosx */
