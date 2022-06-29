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
 * @file CornerPointMeshGenerator.cpp
 */

#include "CornerPointMeshGenerator.hpp"
#include "mesh/generators/CellBlockManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshBody.hpp"
#include "mesh/ElementType.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include <numeric>

namespace geosx
{
using namespace dataRepository;
using namespace cornerPointMesh;

CornerPointMeshGenerator::CornerPointMeshGenerator( string const & name, Group * const parent ):
  ExternalMeshGeneratorBase( name, parent ),
  m_permeabilityUnitInInputFile( PermeabilityUnit::Millidarcy ),
  m_coordinatesUnitInInputFile( CoordinatesUnit::Meter ),
  m_toSquareMeter( 1.0 ),
  m_toMeter( 1.0 )
{

  registerWrapper( viewKeyStruct::permeabilityUnitInInputFileString(), &m_permeabilityUnitInInputFile ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( PermeabilityUnit::Millidarcy ).
    setDescription( "Flag to specify the unit of permeability in the input file. \n"
                    "Two options are available: Millidarcy (default) or SquareMeter. \n"
                    "If Millidarcy is chosen, a conversion factor is applied as GEOSX internally works with square meters." );

  registerWrapper( viewKeyStruct::coordinatesUnitInInputFileString(), &m_coordinatesUnitInInputFile ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( CoordinatesUnit::Meter ).
    setDescription( "Flag to specify the unit of the coordinates in the input file. \n"
                    "Two options are available: Meter (default) or Foot. \n"
                    "If Foot is chosen, a conversion factor is applied as GEOSX internally works with meters." );

  string const builderName = "cpMeshBuilder";
  m_cpMeshBuilder = std::make_unique< cornerPointMesh::CornerPointMeshBuilder >( builderName );
}

CornerPointMeshGenerator::~CornerPointMeshGenerator()
{}

void CornerPointMeshGenerator::postProcessInput()
{
  if( m_permeabilityUnitInInputFile == PermeabilityUnit::SquareMeter )
  {
    m_toSquareMeter = 1.0;
  }
  else if( m_permeabilityUnitInInputFile == PermeabilityUnit::Millidarcy )
  {
    m_toSquareMeter = 9.869232667160128e-16;
  }
  else
  {
    GEOSX_THROW( "This unit for permeability is not supported", InputError );
  }

  if( m_coordinatesUnitInInputFile == CoordinatesUnit::Meter )
  {
    m_toMeter = 1.0;
  }
  else if( m_coordinatesUnitInInputFile == CoordinatesUnit::Foot )
  {
    m_toMeter = 0.3048;
  }
  else
  {
    GEOSX_THROW( "This unit for coordinates is not supported", InputError );
  }

  // Isaac: call buildMesh method to get all maps needed
  m_cpMeshBuilder->buildMesh( m_filePath );
}

Group * CornerPointMeshGenerator::createChild( string const & GEOSX_UNUSED_PARAM( childKey ),
                                               string const & GEOSX_UNUSED_PARAM( childName ) )
{
  return nullptr;
}

void CornerPointMeshGenerator::generateMesh( DomainPartition & domain )
{
  MeshBody & meshBody = domain.getMeshBodies().registerGroup< MeshBody >( this->getName() );
  meshBody.createMeshLevel( 0 );
  CellBlockManager & cellBlockManager = meshBody.registerGroup< CellBlockManager >( keys::cellManager );

  // Step 0: transfer the neighbor list

  // TODO: change the name, we don't use metis
  domain.getMetisNeighborList() = m_cpMeshBuilder->neighborsList();

  // we can start constructing the mesh

  // Step 1: fill vertex information

  arrayView2d< real64 const > vertexPositions = m_cpMeshBuilder->vertexPositions();
  arrayView1d< globalIndex const > vertexToGlobalVertex = m_cpMeshBuilder->vertexToGlobalVertex();
  localIndex const nVertices = vertexPositions.size( 0 );
  cellBlockManager.setNumNodes( nVertices );
  arrayView1d< globalIndex > const & vertexLocalToGlobal = cellBlockManager.getNodeLocalToGlobal();
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & X = cellBlockManager.getNodePositions();

  // Generate the "all" set
  array1d< localIndex > allNodes( nVertices );
  std::iota( allNodes.begin(), allNodes.end(), 0 );
  SortedArray< localIndex > & allNodeSet = cellBlockManager.getNodeSets()[ "all" ];
  allNodeSet.insert( allNodes.begin(), allNodes.end() );

  real64 xMin[3] = { std::numeric_limits< real64 >::max() };
  real64 xMax[3] = { std::numeric_limits< real64 >::min() };

  for( localIndex iVertex = 0; iVertex < nVertices; ++iVertex )
  {
    for( integer i = 0; i < 3; ++i )
    {
      X( iVertex, i ) = ( m_toMeter * vertexPositions( iVertex, i ) + m_translate[i] ) * m_scale[i];
      xMax[i] = std::max( xMax[i], X( iVertex, i ) );
      xMin[i] = std::min( xMin[i], X( iVertex, i ) );
    }
    vertexLocalToGlobal( iVertex ) = vertexToGlobalVertex( iVertex );
  }
  LvArray::tensorOps::subtract< 3 >( xMax, xMin );
  meshBody.setGlobalLengthScale( LvArray::tensorOps::l2Norm< 3 >( xMax ) );

  // Step 2: fill cell information for each region

  ArrayOfArraysView< localIndex const > regionId = m_cpMeshBuilder->regionId();

  for( localIndex er = 0; er < regionId.size(); ++er )
  {
    string const cellBlockName = "DEFAULT_HEX_"+std::to_string( er );
    CellBlock & cellBlock = cellBlockManager.registerCellBlock( cellBlockName );

    if( regionId.sizeOfArray( er ) == 0 )
    {
      continue;
    }

    // temporary accessors while we only support the conforming case
    // this is not what we ultimately want, which is instead:
    //  -> map from elem to faces
    //  -> map from face to nodes
    // what is below is a temporary mess (that maybe, should be hidden in CPMeshData)
    arraySlice1d< localIndex const > ownedActiveCellsInRegion = regionId[er];
    arrayView1d< localIndex const > ownedActiveCellToActiveCell = m_cpMeshBuilder->ownedActiveCellToActiveCell();
    arrayView1d< globalIndex const > ownedActiveCellToGlobalCell = m_cpMeshBuilder->ownedActiveCellToGlobalCell();

    arrayView1d< localIndex const > activeCellToCell = m_cpMeshBuilder->activeCellToCell();
    arrayView1d< localIndex const > cellToCPVertices = m_cpMeshBuilder->cellToCPVertices();
    arrayView1d< localIndex const > cpVertexToVertex = m_cpMeshBuilder->cpVertexToVertex();

    localIndex const nOwnedActiveCellsInRegion = ownedActiveCellsInRegion.size();
    cellBlock.setElementType( ElementType::Hexahedron );
    cellBlock.resize( nOwnedActiveCellsInRegion );

    arrayView1d< globalIndex > cellLocalToGlobal = cellBlock.localToGlobalMap();
    arrayView2d< localIndex, cells::NODE_MAP_USD > const cellToVertex = cellBlock.getElemToNode(); // CPG has an arrayofarray
    //cellToVertex.resize( nOwnedActiveCellsInRegion, 8 );

    for( localIndex iOwnedActiveCellInRegion = 0; iOwnedActiveCellInRegion < nOwnedActiveCellsInRegion; ++iOwnedActiveCellInRegion )
    {
      localIndex const iOwnedActiveCell = ownedActiveCellsInRegion( iOwnedActiveCellInRegion );
      // this filtering is needed because the CPMeshBuilder uses a layer of cells around each MPI partition
      // to facilitate the treatment of non-matching faces at the boundary between MPI partitions.
      // But, it should be hidden in CPMeshBuilder/CPMeshData. I will take care of this when we handle the non-conforming case.
      localIndex const iActiveCell = ownedActiveCellToActiveCell( iOwnedActiveCell );
      localIndex const iFirstCPVertex = cellToCPVertices( activeCellToCell( iActiveCell ) );

      // temporary code while we don't support the non-conforming case
      for( localIndex i = 0; i < 8; ++i )
      {
        cellToVertex( iOwnedActiveCellInRegion, i ) = cpVertexToVertex( iFirstCPVertex + i );
      }
      cellLocalToGlobal( iOwnedActiveCellInRegion ) = ownedActiveCellToGlobalCell( iOwnedActiveCellInRegion );
    }
  }
  cellBlockManager.buildMaps();
}

void CornerPointMeshGenerator::freeResources()
{
  // do something smart here
}

namespace
{

std::vector< localIndex > findFieldNamesForImport( arrayView1d< string const > const & srcFieldNames )
{
  std::vector< localIndex > arrays; // save index
  std::vector< string > validFieldNames = { "PORO", "PERM" }; // Now we only import PORO and PERM field from GRDECL files
  for( string const & sourceName : srcFieldNames )
  {
    auto fieldNameIndex = std::find( validFieldNames.begin(), validFieldNames.end(), sourceName );

    bool const isValid = fieldNameIndex != validFieldNames.end() ? true : false;

    GEOSX_THROW_IF( isValid == false,
                    GEOSX_FMT( "Source field '{}' not found in dataset", sourceName ),
                    InputError );

    localIndex fieldPos = fieldNameIndex -  validFieldNames.begin();
    arrays.push_back( fieldPos );
  }

  return arrays;
}

} // namespace

void CornerPointMeshGenerator::importFields( DomainPartition & domain ) const
{
  GEOSX_LOG_RANK_0( "Importing field data from mesh dataset" );
  GEOSX_ASSERT_MSG( m_cpMeshBuilder, "Must call generateMesh() before importFields()" );

  ElementRegionManager & elemManager = domain.getMeshBody( this->getName() ).getMeshLevel( 0 ).getElemManager();
  ArrayOfArraysView< localIndex const > regionId = m_cpMeshBuilder->regionId();

  std::vector< localIndex > const indexArrays = findFieldNamesForImport( m_fieldsToImport );

  elemManager.forElementSubRegionsComplete< CellElementSubRegion >( [&]( localIndex const er,
                                                                         localIndex const,
                                                                         ElementRegionBase &,
                                                                         CellElementSubRegion & subRegion )
  {
    if( regionId.sizeOfArray( er ) == 0 || indexArrays.size() == 0 )
    {
      return;
    }

    arraySlice1d< localIndex const > ownedActiveCellsInRegion = regionId[er];
    arrayView1d< localIndex const > ownedActiveCellToActiveCell = m_cpMeshBuilder->ownedActiveCellToActiveCell();
    arrayView1d< localIndex const > activeCellToCell = m_cpMeshBuilder->activeCellToCell();

    localIndex const nOwnedActiveCellsInRegion = ownedActiveCellsInRegion.size();

    // Step 3: fill property information

    // TODO: here, just copy over what is done in PAMELAMeshGenerator, if it works
    // Writing properties
    for( std::size_t i = 0; i < indexArrays.size(); ++i )
    {
      // Get source
      // Find destination
      string const wrapperName = m_fieldNamesInGEOSX[indexArrays[i]];
      if( subRegion.hasWrapper( wrapperName ) )
      {
        // Step 3.a: fill porosity in active cells
        if( m_fieldsToImport[i] == "PORO" )
        {
          arrayView1d< real64 const > porosityField = m_cpMeshBuilder->porosityField();
          if( !porosityField.empty() )
          {
            arrayView1d< real64 > & referencePorosity = subRegion.getReference< array1d< real64 > >( wrapperName );
            for( localIndex iOwnedActiveCellInRegion = 0; iOwnedActiveCellInRegion < nOwnedActiveCellsInRegion; ++iOwnedActiveCellInRegion )
            {
              localIndex const iOwnedActiveCell = ownedActiveCellsInRegion( iOwnedActiveCellInRegion );
              localIndex const iActiveCell = ownedActiveCellToActiveCell( iOwnedActiveCell );
              localIndex const iCell = activeCellToCell( iActiveCell );
              referencePorosity( iOwnedActiveCellInRegion ) = porosityField( iCell );
            }
          }
        }

        if( m_fieldsToImport[i] == "PERM" )
        {
          // Step 3.b: fill permeability in active cells
          arrayView2d< real64 const > permeabilityField = m_cpMeshBuilder->permeabilityField();
          if( !permeabilityField.empty())
          {
            arrayView3d< real64 > & permeability = subRegion.getReference< array3d< real64 > >( wrapperName );
            for( localIndex iOwnedActiveCellInRegion = 0; iOwnedActiveCellInRegion < nOwnedActiveCellsInRegion; ++iOwnedActiveCellInRegion )
            {
              localIndex const iOwnedActiveCell = ownedActiveCellsInRegion( iOwnedActiveCellInRegion );
              localIndex const iActiveCell = ownedActiveCellToActiveCell( iOwnedActiveCell );
              localIndex const iCell = activeCellToCell( iActiveCell );
              for( int q = 0; q < permeability.size( 1 ); ++q )
              {
                for( localIndex dim = 0; dim < 3; dim++ )
                {
                  permeability( iOwnedActiveCellInRegion, q, dim ) =
                    LvArray::math::max( 1e-19, m_toSquareMeter * permeabilityField( iCell, dim ) );
                }
              }
            }
          }
        }
      }
    }
  } );

  CommunicationTools::getInstance().synchronizeFields( { { "elems", m_fieldNamesInGEOSX } },
                                                       domain.getMeshBody( this->getName() ).getMeshLevel( 0 ),
                                                       domain.getNeighbors(),
                                                       false );
}

REGISTER_CATALOG_ENTRY( MeshGeneratorBase, CornerPointMeshGenerator, string const &, Group * const )
}
