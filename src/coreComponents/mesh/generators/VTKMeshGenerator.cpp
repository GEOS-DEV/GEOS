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
 * @file VTKMeshGenerator.cpp
 */

#include "VTKMeshGenerator.hpp"

#include "mesh/generators/VTKFaceBlockUtilities.hpp"
#include "mesh/generators/VTKMeshGeneratorTools.hpp"
#include "mesh/generators/CellBlockManager.hpp"
#include "mesh/MeshBody.hpp"
#include "common/DataTypes.hpp"
#include "common/DataLayouts.hpp"
#include "common/MpiWrapper.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"

namespace geosx
{
using namespace dataRepository;


VTKMeshGenerator::VTKMeshGenerator( string const & name,
                                    Group * const parent )
  : ExternalMeshGeneratorBase( name, parent )
{
  registerWrapper( viewKeyStruct::regionAttributeString(), &m_attributeName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( "attribute" ).
    setDescription( "Name of the VTK cell attribute to use as region marker" );

  registerWrapper( viewKeyStruct::nodesetNamesString(), &m_nodesetNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Names of the VTK nodesets to import" );

  registerWrapper( viewKeyStruct::faceBlockNamesString(), &m_faceBlockNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Names of the VTK faces to import" );

  registerWrapper( viewKeyStruct::partitionRefinementString(), &m_partitionRefinement ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1 ).
    setDescription( "Number of partitioning refinement iterations (defaults to 1, recommended value)."
                    "A value of 0 disables graph partitioning and keeps simple kd-tree partitions (not recommended). "
                    "Values higher than 1 may lead to slightly improved partitioning, but yield diminishing returns." );

  registerWrapper( viewKeyStruct::partitionMethodString(), &m_partitionMethod ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Method (library) used to partition the mesh" );

  registerWrapper( viewKeyStruct::useGlobalIdsString(), &m_useGlobalIds ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Controls the use of global IDs in the input file for cells and points."
                    " If set to 0 (default value), the GlobalId arrays in the input mesh are used if available, and generated otherwise."
                    " If set to a negative value, the GlobalId arrays in the input mesh are not used, and generated global Ids are automatically generated."
                    " If set to a positive value, the GlobalId arrays in the input mesh are used and required, and the simulation aborts if they are not available" );
}

void VTKMeshGenerator::generateMesh( DomainPartition & domain )
{
  // TODO refactor void MeshGeneratorBase::generateMesh( DomainPartition & domain )
  GEOSX_MARK_FUNCTION;

  MPI_Comm const comm = MPI_COMM_GEOSX;
  vtkSmartPointer< vtkMultiProcessController > controller = vtk::getController();
  vtkMultiProcessController::SetGlobalController( controller );

  GEOSX_LOG_RANK_0( GEOSX_FMT( "{} '{}': reading mesh from {}", catalogName(), getName(), m_filePath ) );
  {
    GEOSX_LOG_LEVEL_RANK_0( 2, "  reading the dataset..." );
    vtkSmartPointer< vtkDataSet > loadedMesh = vtk::loadMesh( m_filePath );
    GEOSX_LOG_LEVEL_RANK_0( 2, "  redistributing mesh..." );
    m_vtkMesh = vtk::redistributeMesh( *loadedMesh, comm, m_partitionMethod, m_partitionRefinement, m_useGlobalIds );
    GEOSX_LOG_LEVEL_RANK_0( 2, "  finding neighbor ranks..." );
    std::vector< vtkBoundingBox > boxes = vtk::exchangeBoundingBoxes( *m_vtkMesh, comm );
    std::vector< int > const neighbors = vtk::findNeighborRanks( std::move( boxes ) );
    domain.getMetisNeighborList().insert( neighbors.begin(), neighbors.end() );
    GEOSX_LOG_LEVEL_RANK_0( 2, "  done!" );
  }

  GEOSX_LOG_RANK_0( GEOSX_FMT( "{} '{}': generating GEOSX mesh data structure", catalogName(), getName() ) );

  MeshBody & meshBody = domain.getMeshBodies().registerGroup< MeshBody >( this->getName() );
  meshBody.createMeshLevel( 0 );

  CellBlockManager & cellBlockManager = meshBody.registerGroup< CellBlockManager >( keys::cellManager );

  GEOSX_LOG_LEVEL_RANK_0( 2, "  preprocessing..." );
  m_cellMap = vtk::buildCellMap( *m_vtkMesh, m_attributeName );

  GEOSX_LOG_LEVEL_RANK_0( 2, "  writing nodes..." );
  real64 const globalLength = writeNodes( getLogLevel(), *m_vtkMesh, m_nodesetNames, cellBlockManager, this->m_translate, this->m_scale );
  meshBody.setGlobalLengthScale( globalLength );

  GEOSX_LOG_LEVEL_RANK_0( 2, "  writing cells..." );
  writeCells( getLogLevel(), *m_vtkMesh, m_cellMap, cellBlockManager );

  GEOSX_LOG_LEVEL_RANK_0( 2, "  writing surfaces..." );
  writeSurfaces( getLogLevel(), *m_vtkMesh, m_cellMap, cellBlockManager );

  GEOSX_LOG_LEVEL_RANK_0( 2, "  building connectivity maps..." );
  cellBlockManager.buildMaps();

  for( string const & faceBlockName: m_faceBlockNames )
  {
    importFractureNetwork( faceBlockName, m_vtkMesh, cellBlockManager );
  }

  GEOSX_LOG_LEVEL_RANK_0( 2, "  done!" );
  vtk::printMeshStatistics( *m_vtkMesh, m_cellMap, comm );
}

string VTKMeshGenerator::getSourceName( localIndex index ) const {
  std::vector< vtkDataArray * > const srcArrays = vtk::findArraysForImport( *m_vtkMesh, m_fieldsToImport );
  vtkDataArray * vtkArray = srcArrays[index];
  return vtkArray->GetName();
}

void VTKMeshGenerator::importFieldsOnArray( string const regionName, string const wrapperName, int fieldNameIdx, WrapperBase & wrapper, bool importMaterial ) const
{
  GEOSX_LOG_RANK_0( GEOSX_FMT( "{} '{}': importing field data from mesh dataset", catalogName(), getName() ) );
  GEOSX_ASSERT_MSG( m_vtkMesh, "Must call generateMesh() before importFields()" );

  std::vector< vtkDataArray * > const srcArrays = vtk::findArraysForImport( *m_vtkMesh, m_fieldsToImport );
  vtkDataArray * vtkArray = srcArrays[fieldNameIdx];

  for( auto const & typeRegions : m_cellMap )
  {
    // Restrict data import to 3D cells
    if( getElementDim( typeRegions.first ) == 3 )
    {
      for( auto const & regionCells: typeRegions.second )
      {
        string const cellBlockName = vtk::buildCellBlockName( typeRegions.first, regionCells.first );
        // We don't know how the user mapped cell blocks to regions, so we must check all of them
        if( regionName != cellBlockName ) continue;

        GEOSX_LOG_LEVEL_RANK_0( 1, GEOSX_FMT( "Importing field {} -> {}",
                                              vtkArray->GetName(), wrapperName ) );
        if( importMaterial )
        {
          vtk::importMaterialField( regionCells.second, vtkArray, wrapper );
        }
        else
        {
          vtk::importRegularField( regionCells.second, vtkArray, wrapper );
        }
      }
    }
  }
}

void VTKMeshGenerator::freeResources()
{
  m_vtkMesh = nullptr;
  m_cellMap.clear();
}

REGISTER_CATALOG_ENTRY( MeshGeneratorBase, VTKMeshGenerator, string const &, Group * const )

} // namespace geosx
