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
#include "common/DataTypes.hpp"
#include "common/DataLayouts.hpp"
#include "common/MpiWrapper.hpp"

namespace geos
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

  registerWrapper( viewKeyStruct::mainBlockNameString(), &m_mainBlockName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( "main" ).
    setDescription( "For multi-block files, name of the 3d mesh block." );

  registerWrapper( viewKeyStruct::faceBlockNamesString(), &m_faceBlockNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "For multi-block files, names of the face mesh block." );

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

void VTKMeshGenerator::fillCellBlockManager( CellBlockManager & cellBlockManager, array1d< int > const & )
{
  // TODO refactor void MeshGeneratorBase::generateMesh( DomainPartition & domain )
  GEOS_MARK_FUNCTION;

  MPI_Comm const comm = MPI_COMM_GEOSX;
  vtkSmartPointer< vtkMultiProcessController > controller = vtk::getController();
  vtkMultiProcessController::SetGlobalController( controller );

  logger.rank0Log( GEOS_FMT( "{} '{}': reading mesh from {}", catalogName(), getName(), m_filePath ) );
  {
    GEOS_LOG_LEVEL_RANK_0( 2, "  reading the dataset..." );
    vtkSmartPointer< vtkDataSet > loadedMesh = vtk::loadMesh( m_filePath, m_mainBlockName );
    GEOS_LOG_LEVEL_RANK_0( 2, "  redistributing mesh..." );
    m_vtkMesh = vtk::redistributeMesh( loadedMesh, comm, m_partitionMethod, m_partitionRefinement, m_useGlobalIds );
    GEOS_LOG_LEVEL_RANK_0( 2, "  finding neighbor ranks..." );
    std::vector< vtkBoundingBox > boxes = vtk::exchangeBoundingBoxes( *m_vtkMesh, comm );
    std::vector< int > const neighbors = vtk::findNeighborRanks( std::move( boxes ) );
    m_spatialPartition.setMetisNeighborList( std::move( neighbors ) );
    GEOS_LOG_LEVEL_RANK_0( 2, "  done!" );
  }
  logger.rank0Log( GEOS_FMT( "{} '{}': generating GEOSX mesh data structure", catalogName(), getName() ) );


  GEOS_LOG_LEVEL_RANK_0( 2, "  preprocessing..." );
  m_cellMap = vtk::buildCellMap( *m_vtkMesh, m_attributeName );

  GEOS_LOG_LEVEL_RANK_0( 2, "  writing nodes..." );
  cellBlockManager.setGlobalLength( writeNodes( getLogLevel(), *m_vtkMesh, m_nodesetNames, cellBlockManager, this->m_translate, this->m_scale ) );

  GEOS_LOG_LEVEL_RANK_0( 2, "  writing cells..." );
  writeCells( getLogLevel(), *m_vtkMesh, m_cellMap, cellBlockManager );

  GEOS_LOG_LEVEL_RANK_0( 2, "  writing surfaces..." );
  writeSurfaces( getLogLevel(), *m_vtkMesh, m_cellMap, cellBlockManager );

  GEOS_LOG_LEVEL_RANK_0( 2, "  building connectivity maps..." );
  cellBlockManager.buildMaps();

  for( string const & faceBlockName: m_faceBlockNames )
  {
    m_faceBlockMeshes[faceBlockName] = importFractureNetwork( m_filePath, faceBlockName, m_vtkMesh, cellBlockManager );
  }

  GEOS_LOG_LEVEL_RANK_0( 2, "  done!" );
  vtk::printMeshStatistics( *m_vtkMesh, m_cellMap, comm );
}

void VTKMeshGenerator::importVolumicFieldOnArray( string const & cellBlockName,
                                                  string const & meshFieldName,
                                                  bool isMaterialField,
                                                  dataRepository::WrapperBase & wrapper ) const
{
  for( auto const & typeRegions : m_cellMap )
  {
    // Restrict data import to 3D cells
    if( getElementDim( typeRegions.first ) == 3 )
    {
      for( auto const & regionCells: typeRegions.second )
      {
        string const currentCellBlockName = vtk::buildCellBlockName( typeRegions.first, regionCells.first );
        // We don't know how the user mapped cell blocks to regions, so we must check all of them
        if( cellBlockName != currentCellBlockName )
          continue;

        vtkDataArray * vtkArray = vtk::findArrayForImport( *m_vtkMesh, meshFieldName );
        if( isMaterialField )
        {
          return vtk::importMaterialField( regionCells.second, vtkArray, wrapper );
        }
        else
        {
          return vtk::importRegularField( regionCells.second, vtkArray, wrapper );
        }
      }
    }
  }

  GEOS_ERROR( "Could not import field \"" << meshFieldName << "\" from cell block \"" << cellBlockName << "\"." );
}


void VTKMeshGenerator::importSurfacicFieldOnArray( string const & faceBlockName,
                                                   string const & meshFieldName,
                                                   dataRepository::WrapperBase & wrapper ) const
{
  // If the field was not imported from cell blocks, we now look for it in the face blocks.
  // This surely can be improved by clearly stating which field should be on which block.
  // Note that there is no additional work w.r.t. the cells on which we want to import the fields,
  // because the face blocks are heterogeneous.
  // We always take the whole data, we do not select cell type by cell type.
  for( auto const & p: m_faceBlockMeshes )
  {
    vtkSmartPointer< vtkDataSet > faceMesh = p.second;
    if( vtk::hasArray( *faceMesh, meshFieldName ) )
    {
      vtkDataArray * vtkArray = vtk::findArrayForImport( *faceMesh, meshFieldName );
      return vtk::importRegularField( vtkArray, wrapper );
    }
  }

  GEOS_ERROR( "Could not import field \"" << meshFieldName << "\" from face block \"" << faceBlockName << "\"." );
}


void VTKMeshGenerator::importFieldOnArray( Block block,
                                           string const & blockName,
                                           string const & meshFieldName,
                                           bool isMaterialField,
                                           dataRepository::WrapperBase & wrapper ) const
{
  GEOS_ASSERT_MSG( m_vtkMesh, "Must call generateMesh() before importFields()" );

  switch( block )
  {
    case MeshGeneratorBase::Block::VOLUMIC:
      return importVolumicFieldOnArray( blockName, meshFieldName, isMaterialField, wrapper );
    case MeshGeneratorBase::Block::SURFACIC:
    case MeshGeneratorBase::Block::LINEIC:
      return importSurfacicFieldOnArray( blockName, meshFieldName, wrapper );
  }
}

void VTKMeshGenerator::freeResources()
{
  m_vtkMesh = nullptr;
  m_cellMap.clear();
  m_faceBlockMeshes.clear();
}

REGISTER_CATALOG_ENTRY( MeshGeneratorBase, VTKMeshGenerator, string const &, Group * const )

} // namespace geos
