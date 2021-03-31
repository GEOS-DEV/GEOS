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
 * @file CornerPointMeshGenerator.cpp
 */

#include "CornerPointMeshGenerator.hpp"

#include "managers/DomainPartition.hpp"
#include "mesh/MeshBody.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace CPMesh;

CornerPointMeshGenerator::CornerPointMeshGenerator( string const & name, Group * const parent ):
  MeshGeneratorBase( name, parent )
{
  registerWrapper( viewKeyStruct::filePathString(), &m_filePath ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Path to the mesh file" );

  string const builderName = "cpMeshBuilder";
  m_cPMeshBuilder = std::make_unique< CPMesh::CPMeshBuilder >( builderName );
}

CornerPointMeshGenerator::~CornerPointMeshGenerator()
{}

void CornerPointMeshGenerator::postProcessInput()
{
  m_cPMeshBuilder->buildMesh( m_filePath );
}

Group * CornerPointMeshGenerator::createChild( string const & GEOSX_UNUSED_PARAM( childKey ),
                                               string const & GEOSX_UNUSED_PARAM( childName ) )
{
  return nullptr;
}

void CornerPointMeshGenerator::generateMesh( DomainPartition & domain )
{
  Group & meshBodies = domain.getGroup( string( "MeshBodies" ));
  MeshBody & meshBody = meshBodies.registerGroup< MeshBody >( this->getName() );

  MeshLevel & meshLevel0 = meshBody.registerGroup< MeshLevel >( string( "Level0" ));
  NodeManager & nodeManager = meshLevel0.getNodeManager();
  CellBlockManager & cellBlockManager = domain.getGroup< CellBlockManager >( keys::cellManager );

  // Step 0: transfer the neighbor list

  // TODO: change the name, we don't use metis
  domain.getMetisNeighborList() = m_cPMeshBuilder->neighborsList();

  // we can start constructing the mesh

  // Step 1: fill vertex information

  arrayView2d< real64 const > vertices = m_cPMeshBuilder->vertices();
  arrayView1d< globalIndex const > vertexToGlobalVertex = m_cPMeshBuilder->vertexToGlobalVertex();
  localIndex const nVertices = vertices.size( 0 );
  nodeManager.resize( nVertices );
  arrayView1d< globalIndex > const & vertexLocalToGlobal = nodeManager.localToGlobalMap();
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();

  Group & vertexSets = nodeManager.sets();
  SortedArray< localIndex > & allVertices =
    vertexSets.registerWrapper< SortedArray< localIndex > >( string( "all" ) ).reference();

  real64 xMin[3] = { std::numeric_limits< real64 >::max() };
  real64 xMax[3] = { std::numeric_limits< real64 >::min() };

  for( localIndex iVertex = 0; iVertex < nVertices; ++iVertex )
  {
    X( iVertex, 0 ) = vertices( iVertex, 0 );
    X( iVertex, 1 ) = vertices( iVertex, 1 );
    X( iVertex, 2 ) = vertices( iVertex, 2 );
    allVertices.insert( iVertex );
    vertexLocalToGlobal( iVertex ) = vertexToGlobalVertex( iVertex );

    for( int dim = 0; dim < 3; dim++ )
    {
      if( X( iVertex, dim ) > xMax[dim] )
      {
        xMax[dim] = X( iVertex, dim );
      }
      if( X( iVertex, dim ) < xMin[dim] )
      {
        xMin[dim] = X( iVertex, dim );
      }
    }
  }
  LvArray::tensorOps::subtract< 3 >( xMax, xMin );
  meshBody.setGlobalLengthScale( LvArray::tensorOps::l2Norm< 3 >( xMax ) );

  // Step 2: fill cell information

  CellBlock * cellBlock = &cellBlockManager.getGroup( keys::cellBlocks ).registerGroup< CellBlock >( "DEFAULT_HEX" );

  // temporary accessors while we only support the conforming case
  // this is not what we ultimately want, which is instead:
  //  -> map from elem to faces
  //  -> map from face to nodes
  // what is below is a temporary mess (that maybe, should be hidden in CPMeshData)
  arrayView1d< localIndex const > activeCellInsidePartitionToActiveCell = m_cPMeshBuilder->activeCellInsidePartitionToActiveCell();
  arrayView1d< globalIndex const > activeCellToGlobalCell = m_cPMeshBuilder->activeCellToGlobalCell();
  arrayView1d< localIndex const > activeCellToCell = m_cPMeshBuilder->activeCellToCell();
  arrayView1d< localIndex const > cellToCPVertices = m_cPMeshBuilder->cellToCPVertices();
  arrayView1d< localIndex const > cPVertexToVertex = m_cPMeshBuilder->cPVertexToVertex();

  localIndex const nActiveCellsInsidePartition = activeCellInsidePartitionToActiveCell.size();
  cellBlock->setElementType( "C3D8" );
  cellBlock->resize( nActiveCellsInsidePartition );

  arrayView1d< globalIndex > cellLocalToGlobal = cellBlock->localToGlobalMap();
  auto & cellToVertex = cellBlock->nodeList(); // TODO: remove auto
  cellToVertex.resize( nActiveCellsInsidePartition, 8 );

  for( localIndex iActiveCellInside = 0; iActiveCellInside < nActiveCellsInsidePartition; ++iActiveCellInside )
  {
    // this filtering is needed because the CPMeshBuilder uses a layer of cells around each MPI partition
    // to facilitate the treatment of non-matching faces at the boundary between MPI partitions.
    // But, it should be hidden in CPMeshBuilder/CPMeshData. I will take care of this when we handle the non-conforming case.
    localIndex const iActiveCell = activeCellInsidePartitionToActiveCell( iActiveCellInside );
    localIndex const iFirstCPVertex = cellToCPVertices( activeCellToCell( iActiveCell ) );

    // temporary code while we don't support the non-conforming case
    cellToVertex( iActiveCellInside, 0 ) = cPVertexToVertex( iFirstCPVertex );
    cellToVertex( iActiveCellInside, 1 ) = cPVertexToVertex( iFirstCPVertex + 1 );
    cellToVertex( iActiveCellInside, 2 ) = cPVertexToVertex( iFirstCPVertex + 3 );
    cellToVertex( iActiveCellInside, 3 ) = cPVertexToVertex( iFirstCPVertex + 2 );
    cellToVertex( iActiveCellInside, 4 ) = cPVertexToVertex( iFirstCPVertex + 4 );
    cellToVertex( iActiveCellInside, 5 ) = cPVertexToVertex( iFirstCPVertex + 5 );
    cellToVertex( iActiveCellInside, 6 ) = cPVertexToVertex( iFirstCPVertex + 7 );
    cellToVertex( iActiveCellInside, 7 ) = cPVertexToVertex( iFirstCPVertex + 6 );

    cellLocalToGlobal( iActiveCellInside ) = activeCellToGlobalCell( iActiveCell );
  }

  // Step 3: fill property information

  if( cellBlock != nullptr && cellBlock->size() > 0 )
  {

    // Step 3.a: fill porosity in active cells
    arrayView1d< real64 const > porosityField = m_cPMeshBuilder->porosityField();
    if( !porosityField.empty() )
    {
      arrayView1d< real64 > referencePorosity = cellBlock->addProperty< array1d< real64 > >( "referencePorosity" ).toView();
      for( localIndex iActiveCellInside = 0; iActiveCellInside < nActiveCellsInsidePartition; ++iActiveCellInside )
      {
        localIndex const iActiveCell = activeCellInsidePartitionToActiveCell( iActiveCellInside );
        referencePorosity( iActiveCellInside ) = porosityField( activeCellToCell( iActiveCell ) );
      }
    }

    // Step 3.b: fill permeability in active cells
    arrayView2d< real64 const > permeabilityField = m_cPMeshBuilder->permeabilityField();
    if( !permeabilityField.empty() )
    {
      array2d< real64 > & permeability = cellBlock->addProperty< array2d< real64 > >( "permeability" );
      permeability.resizeDimension< 1 >( 3 );
      for( localIndex iActiveCellInside = 0; iActiveCellInside < nActiveCellsInsidePartition; ++iActiveCellInside )
      {
        // this filtering is needed because the CPMeshBuilder uses a layer of cells around each MPI partition
        localIndex const iActiveCell = activeCellInsidePartitionToActiveCell( iActiveCellInside );
        for( localIndex dim = 0; dim < 3; dim++ )
        {
          permeability( iActiveCellInside, dim ) = permeabilityField( activeCellToCell( iActiveCell ), dim );
        }
      }
    }
  }
}

REGISTER_CATALOG_ENTRY( MeshGeneratorBase, CornerPointMeshGenerator, string const &, Group * const )
}
