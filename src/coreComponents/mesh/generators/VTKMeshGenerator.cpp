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

/**
 * @file VTKMeshGenerator.cpp
 */

#include "VTKMeshGenerator.hpp"

#include "common/DataTypes.hpp"
#include "mesh/ExternalDataSourceManager.hpp"
#include "mesh/LogLevelsInfo.hpp"
#include "mesh/generators/VTKFaceBlockUtilities.hpp"
#include "mesh/generators/VTKMeshDebug.hpp"
#include "mesh/generators/VTKMeshGeneratorTools.hpp"
#include "mesh/generators/CellBlockManager.hpp"
#include "mesh/mpiCommunications/SpatialPartition.hpp"
#include "mesh/generators/Region.hpp"

#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkAppendFilter.h>
#include <vtkDataSet.h>
#include <vtkCellData.h>

#include <cmath>

namespace geos
{
using namespace dataRepository;

VTKMeshGenerator::VTKMeshGenerator( string const & name,
                                    Group * const parent )
  : ExternalMeshGeneratorBase( name, parent ),
  m_dataSource( nullptr )
{
  getWrapperBase( viewKeyStruct::filePathString()).
    setInputFlag( InputFlags::OPTIONAL );

  registerWrapper( viewKeyStruct::regionAttributeString(), &m_regionAttributeName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( "attribute" ).
    setDescription( "Name of the VTK cell attribute to use as region marker" );

  registerWrapper( viewKeyStruct::structuredIndexAttributeString(), &m_structuredIndexAttributeName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of the VTK cell attribute containing structured cell index (e.g. Cartesian IJK)" );

  registerWrapper( viewKeyStruct::nodesetNamesString(), &m_nodesetNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Names of the VTK nodesets to import" );

  registerWrapper( viewKeyStruct::mainBlockNameString(), &m_mainBlockName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( "main" ).
    setDescription( "For multi-block files, name of the 3d mesh block." );

  registerWrapper( viewKeyStruct::faceBlockNamesString(), &m_faceBlockNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "For multi-block files, names of the face mesh block." );

  registerWrapper( viewKeyStruct::partitionRefinementString(), &m_partitionRefinement ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1 ).
    setDescription( "Number of partitioning refinement iterations (defaults to 1, recommended value)."
                    "A value of 0 disables graph partitioning and keeps the initial scatter partition. "
                    "Values higher than 1 may lead to slightly improved partitioning, but yield diminishing returns." );

  registerWrapper( viewKeyStruct::partitionMethodString(), &m_partitionMethod ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Method (library) used to refine mesh partitioning" );

  registerWrapper( viewKeyStruct::partitionModelString(), &m_partitionModel ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( vtk::PartitionModel::legacy ).
    setDescription( "Partition topology/objective: 'legacy' preserves the existing scatter/refinement path; "
                    "'hybrid' uses exact shared faces and points on a serial root mesh and scatters 3D cells once" );

  registerWrapper( viewKeyStruct::partitionFVMCommunicationWeightString(),
                   &m_hybridPartitionOptions.fvmCommunicationWeight ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1.0 ).
    setDescription( "Relative hybrid objective weight for exact cut-face communication" );

  registerWrapper( viewKeyStruct::partitionFEMCommunicationWeightString(),
                   &m_hybridPartitionOptions.femCommunicationWeight ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1.0 ).
    setDescription( "Relative hybrid objective weight for exact shared-node replication" );

  registerWrapper( viewKeyStruct::partitionNeighborPenaltyString(),
                   &m_hybridPartitionOptions.neighborPenalty ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.0 ).
    setDescription( "Optional weak hybrid objective penalty for each distinct neighboring-rank pair" );

  registerWrapper( viewKeyStruct::partitionImbalanceString(),
                   &m_hybridPartitionOptions.imbalance ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.05 ).
    setDescription( "Relative imbalance tolerances for FVM work, FEM work, and resident memory" );

  registerWrapper( viewKeyStruct::partitionSeedString(), &m_hybridPartitionOptions.seed ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 2022 ).
    setDescription( "Deterministic METIS seed for hybrid root-local partitioning" );

  registerWrapper( viewKeyStruct::partitionFVMWeightFieldString(),
                   &m_hybridPartitionOptions.fvmWeightField ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Optional scalar VTK cell array overriding topology-derived FVM work" );

  registerWrapper( viewKeyStruct::partitionFEMWeightFieldString(),
                   &m_hybridPartitionOptions.femWeightField ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Optional scalar VTK cell array overriding topology-derived FEM work" );

  registerWrapper( viewKeyStruct::partitionMemoryWeightFieldString(),
                   &m_hybridPartitionOptions.memoryWeightField ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Optional scalar VTK cell array overriding topology-derived resident-memory work" );

  registerWrapper( viewKeyStruct::partitionRootGraphMemoryLimitMBString(),
                   &m_hybridPartitionOptions.rootGraphMemoryLimitMB ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Root hybrid-graph peak-memory limit in MiB; zero is unlimited and falls back explicitly when exceeded" );

  registerWrapper( viewKeyStruct::partitionDiagnosticsString(),
                   &m_hybridPartitionOptions.diagnostics ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Enable detailed per-rank hybrid partition load and topology diagnostics" );

  registerWrapper( viewKeyStruct::scatterMethodString(), &m_scatterMethod ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( vtk::ScatterMethod::rcb ).
    setDescription( "Method for initial mesh scatter from rank 0 to all ranks: "
                    "contiguous (cell ID ranges, no geometry), "
                    "cartesian (regular grid using -x/-y/-z partitions), "
                    "rcb (recursive coordinate bisection, default), "
                    "kdtree (VTK built-in kd-tree; automatically falls back to rcb when fractures are present)" );

  registerWrapper( viewKeyStruct::partitionFractureWeightString(), &m_partitionFractureWeight ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Additional weight to fracture-connected super-cells during partitioning" );

  registerWrapper( viewKeyStruct::useGlobalIdsString(), &m_useGlobalIds ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Controls the use of global IDs in the input file for cells and points."
                    " If set to 0 (default value), the GlobalId arrays in the input mesh are used if available, and generated otherwise."
                    " If set to a negative value, the GlobalId arrays in the input mesh are not used, and generated global Ids are automatically generated."
                    " If set to a positive value, the GlobalId arrays in the input mesh are used and required, and the simulation aborts if they are not available" );

  registerWrapper( viewKeyStruct::dataSourceString(), &m_dataSourceName ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of the VTK data source" );

  addLogLevel< logInfo::VTKSteps >();
}

void VTKMeshGenerator::postInputInitialization()
{
  ExternalMeshGeneratorBase::postInputInitialization();

  GEOS_ERROR_IF( m_filePath.empty() && m_dataSourceName.empty(),
                 GEOS_FMT( "Either {} or {} must be specified.",
                           viewKeyStruct::filePathString(), viewKeyStruct::dataSourceString() ),
                 getDataContext() );

  GEOS_ERROR_IF( !m_filePath.empty() && !m_dataSourceName.empty(),
                 GEOS_FMT( "Access to the mesh via file and data source are mutually exclusive. "
                           "You can't set both {} and {} at the same time.",
                           viewKeyStruct::filePathString(), viewKeyStruct::dataSourceString() ),
                 getDataContext() );

  if( !m_dataSourceName.empty())
  {
    ExternalDataSourceManager & externalDataManager = getGroupByPath< ExternalDataSourceManager >( "/Problem/ExternalDataSource" );

    m_dataSource = externalDataManager.getGroupPointer< VTKHierarchicalDataSource >( m_dataSourceName );

    GEOS_THROW_IF( m_dataSource == nullptr,
                   GEOS_FMT( "VTK Data Object Source not found: {}", m_dataSourceName ),
                   InputError, getDataContext() );

    m_dataSource->open();
  }

  GEOS_ERROR_IF( m_partitionRefinement < 0,
                 "partitionRefinement must be nonnegative", getDataContext() );
  if( m_partitionModel == vtk::PartitionModel::hybrid )
  {
#ifndef GEOS_USE_METIS
    GEOS_ERROR( "partitionModel='hybrid' requires GEOS to be built with ENABLE_METIS=ON",
                getDataContext() );
#endif
    GEOS_ERROR_IF( m_partitionMethod != vtk::PartitionMethod::parmetis,
                   "partitionModel='hybrid' requires partitionMethod='parmetis'; PT-Scotch remains available "
                   "through partitionModel='legacy'", getDataContext() );
    GEOS_ERROR_IF( !std::isfinite( m_hybridPartitionOptions.fvmCommunicationWeight ) ||
                   !std::isfinite( m_hybridPartitionOptions.femCommunicationWeight ) ||
                   !std::isfinite( m_hybridPartitionOptions.neighborPenalty ) ||
                   m_hybridPartitionOptions.fvmCommunicationWeight < 0.0 ||
                   m_hybridPartitionOptions.femCommunicationWeight < 0.0 ||
                   m_hybridPartitionOptions.neighborPenalty < 0.0,
                   "Hybrid partition communication weights and neighbor penalty must be finite and nonnegative",
                   getDataContext() );
    GEOS_ERROR_IF( m_hybridPartitionOptions.fvmCommunicationWeight == 0.0 &&
                   m_hybridPartitionOptions.femCommunicationWeight == 0.0,
                   "At least one hybrid partition communication weight must be positive", getDataContext() );
    GEOS_ERROR_IF_NE_MSG( m_hybridPartitionOptions.imbalance.size(), 3,
                          "partitionImbalance must contain exactly three values" );
    for( real64 const tolerance : m_hybridPartitionOptions.imbalance )
    {
      GEOS_ERROR_IF( !std::isfinite( tolerance ) || tolerance < 0.0,
                     "partitionImbalance values must be finite and nonnegative", getDataContext() );
    }
    GEOS_ERROR_IF( m_hybridPartitionOptions.rootGraphMemoryLimitMB < 0,
                   "partitionRootGraphMemoryLimitMB must be nonnegative", getDataContext() );
    GEOS_ERROR_IF( m_hybridPartitionOptions.diagnostics < 0,
                   "partitionDiagnostics must be nonnegative", getDataContext() );
  }

}

void VTKMeshGenerator::fillCellBlockManager( CellBlockManager & cellBlockManager, SpatialPartition & partition )
{
  // TODO refactor void MeshGeneratorBase::generateMesh( DomainPartition & domain )
  GEOS_MARK_FUNCTION;

  MPI_Comm const comm = MPI_COMM_GEOS;
  vtk::meshDebug::log( comm, "VTKMeshGenerator::fillCellBlockManager begin" );
  vtkSmartPointer< vtkMultiProcessController > controller = vtk::getController();
  vtkMultiProcessController::SetGlobalController( controller );

  array1d< int > const & partitions = partition.getPartitions();

  if( m_scatterMethod == vtk::ScatterMethod::cartesian )
  {
    int const product = partitions[0] * partitions[1] * partitions[2];
    GEOS_ERROR_IF( product != MpiWrapper::commSize( comm ),
                   GEOS_FMT( "scatterMethod=\"cartesian\" requires -x * -y * -z = MPI size. "
                             "Got {}x{}x{} = {} but MPI size is {}.",
                             partitions[0], partitions[1], partitions[2],
                             product, MpiWrapper::commSize( comm ) ) );
  }

  GEOS_LOG_LEVEL_RANK_0( logInfo::VTKSteps, "  redistributing mesh..." );
  {
    vtk::AllMeshes allMeshes;

    GEOS_LOG_LEVEL_RANK_0( logInfo::VTKSteps, GEOS_FMT( "{} '{}': reading the dataset...", catalogName(), getName() ) );

    if( !m_filePath.empty())
    {
      GEOS_LOG_RANK_0( GEOS_FMT( "{} '{}': reading mesh from {}", catalogName(), getName(), m_filePath ) );
      vtk::meshDebug::logf( comm,
                            "VTKMeshGenerator loadAllMeshes begin file=%s",
                            m_filePath.c_str() );
      allMeshes = vtk::loadAllMeshes( m_filePath, m_mainBlockName, m_faceBlockNames );
      vtk::meshDebug::logDataSet( comm,
                                  "VTKMeshGenerator loadAllMeshes main result",
                                  allMeshes.getMainMesh() );
      vtk::meshDebug::logf( comm,
                            "VTKMeshGenerator loadAllMeshes faceBlockCount=%lld",
                            static_cast< long long >( allMeshes.getFaceBlocks().size() ) );
    }
    else if( !m_dataSourceName.empty())
    {
      if( MpiWrapper::commRank() == 0 )
      {
        stdVector< vtkSmartPointer< vtkPartitionedDataSet > > vtkPartitions;
        vtkNew< vtkAppendFilter > appender;
        appender->MergePointsOn();
        for( auto & [key, value] : getSubGroups())
        {
          Region const & region = getGroup< Region >( key );

          string path = region.getWrapper< string >( Region::viewKeyStruct::pathInRepositoryString()).reference();
          integer region_id = region.getWrapper< integer >( Region::viewKeyStruct::idString()).reference();

          GEOS_LOG_RANK_0( GEOS_FMT( "{} '{}': reading partition from {}", catalogName(), getName(), path ) );
          vtkPartitionedDataSet * p = m_dataSource->search( path );

          //load the grid
          vtkDataObject * block = p->GetPartition( 0 );
          if( block->IsA( "vtkDataSet" ) )
          {
            vtkSmartPointer< vtkDataSet > dataset = vtkDataSet::SafeDownCast( block );

            vtkIntArray * arr = vtkIntArray::New();
            arr->SetName( m_regionAttributeName.c_str());
            arr->SetNumberOfComponents( 1 );
            arr->SetNumberOfTuples( dataset->GetNumberOfCells());

            arr->FillValue( region_id );

            dataset->GetCellData()->AddArray( arr );
            appender->AddInputDataObject( dataset );
          }
        }
        appender->Update();
        vtkUnstructuredGrid * result = vtkUnstructuredGrid::SafeDownCast( appender->GetOutputDataObject( 0 ) );
        allMeshes.setMainMesh( result );

        //DEBUG code
        vtkNew< vtkXMLUnstructuredGridWriter > writer;
        writer->SetFileName( "tmp_output.vtu" );
        writer->SetInputData( result );
        writer->Write();
      }
      else
      {
        vtkUnstructuredGrid * result = vtkUnstructuredGrid::New();
        allMeshes.setMainMesh( result );
      }
    }

    GEOS_LOG_LEVEL_RANK_0( logInfo::VTKSteps,
                           GEOS_FMT( "{} '{}': redistributing mesh...", catalogName(), getName() ) );
    vtk::meshDebug::logDataSet( comm,
                                "VTKMeshGenerator redistributeMeshes input",
                                allMeshes.getMainMesh() );
    vtk::meshDebug::log( comm, "VTKMeshGenerator redistributeMeshes begin" );
    vtk::HybridPartitionOptions hybridOptions = m_hybridPartitionOptions;
    hybridOptions.refinementPasses = m_partitionRefinement;
    hybridOptions.fractureWeight = m_partitionFractureWeight;
    vtk::AllMeshes redistributedMeshes = vtk::redistributeMeshes( getLogLevel(),
                                                                   allMeshes.getMainMesh(),
                                                                   allMeshes.getFaceBlocks(),
                                                                  comm,
                                                                  m_scatterMethod,
                                                                  partitions.toViewConst(),
                                                                  m_partitionMethod,
                                                                  m_partitionRefinement,
                                                                  m_partitionFractureWeight,
                                                                  m_useGlobalIds,
                                                                  m_structuredIndexAttributeName,
                                                                  m_partitionModel,
                                                                  hybridOptions );
    m_vtkMesh = redistributedMeshes.getMainMesh();
    m_faceBlockMeshes = redistributedMeshes.getFaceBlocks();
    vtk::meshDebug::logDataSet( comm,
                                "VTKMeshGenerator redistributeMeshes main result",
                                m_vtkMesh );
    vtk::meshDebug::logf( comm,
                          "VTKMeshGenerator redistributeMeshes faceBlockResultCount=%lld",
                          static_cast< long long >( m_faceBlockMeshes.size() ) );
    GEOS_LOG_LEVEL_RANK_0( logInfo::VTKSteps, GEOS_FMT( "{} '{}': finding neighbor ranks...", catalogName(), getName() ) );
    stdVector< int > neighbors;
    if( redistributedMeshes.hasExactNeighborRanks() )
    {
      neighbors = redistributedMeshes.getExactNeighborRanks();
      vtk::meshDebug::logf( comm,
                            "VTKMeshGenerator using exact hybrid neighbors count=%lld",
                            static_cast< long long >( neighbors.size() ) );
    }
    else
    {
      vtk::meshDebug::log( comm, "VTKMeshGenerator exchangeBoundingBoxes begin" );
      stdVector< vtkBoundingBox > boxes = vtk::exchangeBoundingBoxes( *m_vtkMesh, MPI_COMM_GEOS );
      vtk::meshDebug::logf( comm,
                            "VTKMeshGenerator exchangeBoundingBoxes end boxes=%lld",
                            static_cast< long long >( boxes.size() ) );
      neighbors = vtk::findNeighborRanks( std::move( boxes ) );
    }
    partition.setMetisNeighborList( std::move( neighbors ) );
    GEOS_LOG_LEVEL_RANK_0( logInfo::VTKSteps, GEOS_FMT( "{} '{}': done!", catalogName(), getName() ) );
  }
  GEOS_LOG_RANK_0( GEOS_FMT( "{} '{}': generating GEOS mesh data structure", catalogName(), getName() ) );


  GEOS_LOG_LEVEL_RANK_0( logInfo::VTKSteps, GEOS_FMT( "{} '{}': preprocessing...", catalogName(), getName() ) );
  vtk::meshDebug::log( comm, "VTKMeshGenerator buildCellMap begin" );
  m_cellMap = vtk::buildCellMap( *m_vtkMesh, m_regionAttributeName );
  vtk::meshDebug::logf( comm,
                        "VTKMeshGenerator buildCellMap end elementTypeCount=%lld",
                        static_cast< long long >( m_cellMap.size() ) );

  GEOS_LOG_LEVEL_RANK_0( logInfo::VTKSteps, GEOS_FMT( "{} '{}': writing nodes...", catalogName(), getName() ) );
  vtk::meshDebug::log( comm, "VTKMeshGenerator writeNodes begin" );
  writeNodes( getLogLevel(), *m_vtkMesh, m_nodesetNames, cellBlockManager, m_translate, m_scale );
  vtk::meshDebug::log( comm, "VTKMeshGenerator writeNodes end" );

  GEOS_LOG_LEVEL_RANK_0( logInfo::VTKSteps, GEOS_FMT( "{} '{}': writing cells...", catalogName(), getName() ) );
  vtk::meshDebug::log( comm, "VTKMeshGenerator writeCells begin" );
  writeCells( getLogLevel(), *m_vtkMesh, m_cellMap, m_structuredIndexAttributeName, cellBlockManager );
  vtk::meshDebug::log( comm, "VTKMeshGenerator writeCells end" );

  GEOS_LOG_LEVEL_RANK_0( logInfo::VTKSteps, GEOS_FMT( "{} '{}': writing surfaces...", catalogName(), getName() ) );
  vtk::meshDebug::log( comm, "VTKMeshGenerator writeSurfaces begin" );
  writeSurfaces( getLogLevel(), *m_vtkMesh, m_cellMap, cellBlockManager );
  vtk::meshDebug::log( comm, "VTKMeshGenerator writeSurfaces end" );

  GEOS_LOG_LEVEL_RANK_0( logInfo::VTKSteps, GEOS_FMT( "{} '{}': building connectivity maps...", catalogName(), getName() ) );
  vtk::meshDebug::log( comm, "VTKMeshGenerator buildMaps begin" );
  cellBlockManager.buildMaps();
  vtk::meshDebug::log( comm, "VTKMeshGenerator buildMaps end" );

  vtk::meshDebug::log( comm, "VTKMeshGenerator getGlobalLengthAndOffset begin" );
  auto lengthAndOffset = getGlobalLengthAndOffset( *m_vtkMesh );
  cellBlockManager.setGlobalLength( lengthAndOffset.first );
  cellBlockManager.setGlobalOffset( lengthAndOffset.second );
  vtk::meshDebug::logf( comm,
                        "VTKMeshGenerator getGlobalLengthAndOffset end globalLength=%lld globalOffset=%lld",
                        static_cast< long long >( lengthAndOffset.first ),
                        static_cast< long long >( lengthAndOffset.second ) );

  for( auto const & [name, mesh]: m_faceBlockMeshes )
  {
    vtk::importFractureNetwork( name, mesh, m_vtkMesh, cellBlockManager );
  }

  GEOS_LOG_LEVEL_RANK_0( logInfo::VTKSteps, GEOS_FMT( "{} '{}': done!", catalogName(), getName() ) );
  vtk::meshDebug::log( comm, "VTKMeshGenerator::fillCellBlockManager end" );
  vtk::printMeshStatistics( *m_vtkMesh, m_cellMap, MPI_COMM_GEOS );
}

void VTKMeshGenerator::importVolumicFieldOnArray( string const & cellBlockName,
                                                  string const & meshFieldName,
                                                  bool isMaterialField,
                                                  dataRepository::WrapperBase & wrapper ) const
{
  for( auto const & typeRegions: m_cellMap )
  {
    // Restrict data import to 3D cells
    if( getElementDim( typeRegions.first ) == 3 )
    {
      for( auto const & regionCells: typeRegions.second )
      {
        string const currentCellBlockName = vtk::buildCellBlockName( typeRegions.first, regionCells.first );
        // We don't know how the user mapped cell blocks to regions, so we must check all of them
        if( cellBlockName != currentCellBlockName )
        {
          continue;
        }

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

  GEOS_ERROR( GEOS_FMT( "Could not import field \"{}\" from cell block \"{}\".", meshFieldName, cellBlockName ),
              getDataContext()  );
}


void VTKMeshGenerator::importSurfacicFieldOnArray( string const & faceBlockName,
                                                   string const & meshFieldName,
                                                   dataRepository::WrapperBase & wrapper ) const
{
  // Note that there is no additional work w.r.t. the cells on which we want to import the fields,
  // because the face blocks are heterogeneous.
  // We always take the whole data, we do not select cell type by cell type.
  vtkSmartPointer< vtkDataSet > faceMesh = m_faceBlockMeshes.at( faceBlockName );

  // I've noticed that there may be some issues when reading empty arrays (empty, not nulls).
  // It looks like we may be reading above the limits of the array; ghosting is surely at stake here.
  if( faceMesh->GetNumberOfCells() == 0 )
  {
    return;
  }

  if( vtk::hasArray( *faceMesh, meshFieldName ) )
  {
    vtkDataArray * vtkArray = vtk::findArrayForImport( *faceMesh, meshFieldName );
    return vtk::importRegularField( vtkArray, wrapper );
  }

  GEOS_ERROR( GEOS_FMT( "Could not import field \"{}\" from face block \"{}\".", meshFieldName, faceBlockName ),
              getDataContext()  );
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
