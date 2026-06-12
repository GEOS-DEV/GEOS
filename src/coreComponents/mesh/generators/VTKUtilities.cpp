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

#include "common/format/table/TableData.hpp"
#include "common/format/table/TableFormatter.hpp"
#include "common/format/table/TableLayout.hpp"
#include "common/TypeDispatch.hpp"

#include "mesh/generators/CollocatedNodes.hpp"
#include "mesh/generators/VTKMeshGeneratorTools.hpp"
#include "mesh/generators/VTKUtilities.hpp"
#include "mesh/MeshFields.hpp"
#include "mesh/generators/VTKSuperCellPartitioning.hpp"

#ifdef GEOS_USE_PARMETIS
#include "mesh/generators/ParMETISInterface.hpp"
#endif
#ifdef GEOS_USE_SCOTCH
#include "mesh/generators/PTScotchInterface.hpp"
#endif

#include <vtkArrayDispatch.h>
#include <vtkBoundingBox.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataSetReader.h>
#include <vtkExtractCells.h>
#include <vtkGenerateGlobalIds.h>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkInformationStringKey.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkNew.h>
#include <vtkPartitionedDataSet.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkRectilinearGrid.h>
#include <vtkRectilinearGridReader.h>
#include <vtkRedistributeDataSetFilter.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridReader.h>
#include <vtkStructuredPoints.h>
#include <vtkStructuredPointsReader.h>
#include <vtkThreshold.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLImageDataReader.h>
#include <vtkXMLMultiBlockDataReader.h>
#include <vtkXMLPImageDataReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPPolyDataReader.h>
#include <vtkXMLPRectilinearGridReader.h>
#include <vtkXMLPStructuredGridReader.h>
#include <vtkXMLPUnstructuredGridReader.h>
#include <vtkXMLRectilinearGridReader.h>
#include <vtkXMLStructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkAppendFilter.h>

#ifdef GEOS_USE_MPI
#include <vtkMPIController.h>
#include <vtkMPI.h>
#else
#include <vtkDummyController.h>
#endif

#include <numeric>
#include <type_traits>


namespace geos
{
using namespace dataRepository;

namespace vtk
{

/**
 * @brief Supported VTK Mesh file extensions
 */
enum class VTKMeshExtension : integer
{
  vtm,  ///< XML multi-block container
  vtk,  ///< Legacy serial format
  vtu,  ///< XML serial vtkUnstructuredGrid (unstructured)
  vtr,  ///< XML serial vtkRectilinearGrid (structured)
  vts,  ///< XML serial vtkStructuredGrid (structured)
  vti,  ///< XML serial vtkImageData (structured)
  vtp,  ///< XML serial vtkPolyData
  pvtu, ///< XML parallel vtkUnstructuredGrid (unstructured)
  pvtr, ///< XML parallel vtkRectilinearGrid (structured)
  pvts, ///< XML parallel vtkStructuredGrid (structured)
  pvti, ///< XML parallel vtkImageData (structured)
  pvtp, ///< XML parallel vtkPolyData
};

/// Strings for VTKMeshGenerator::VTKMeshExtension enumeration
ENUM_STRINGS( VTKMeshExtension,
              "vtm",
              "vtk",
              "vtu",
              "vtr",
              "vts",
              "vti",
              "vtp",
              "pvtu",
              "pvtr",
              "pvts",
              "pvti",
              "pvtp" );

/**
 * @brief Supported VTK legacy dataset types
 */
enum class VTKLegacyDatasetType : integer
{
  structuredPoints, ///< Structured points (structured)
  structuredGrid,   ///< Structured grid (structured)
  unstructuredGrid, ///< Unstructured grid (unstructured)
  rectilinearGrid,  ///< Rectilinear grid (structured)
  polyData,         ///< PolyData
};

/// Strings for VTKMeshGenerator::VTKLegacyDatasetType enumeration
ENUM_STRINGS( VTKLegacyDatasetType,
              "structuredPoints",
              "structuredGrid",
              "unstructuredGrid",
              "rectilinearGrid",
              "polyData" );

/**
 * @brief Gathers all the data from all ranks, merge them, sort them, and remove duplicates.
 * @tparam T Type of the exchanged data.
 * @param data The data to be exchanged.
 * @return The merged data.
 * @note This function makes MPI calls.
 */
template< typename T >
stdVector< T > collectUniqueValues( stdVector< T > const & data )
{
  // Exchange the sizes of the data across all ranks.
  array1d< int > dataSizes( MpiWrapper::commSize() );
  MpiWrapper::allGather( LvArray::integerConversion< int >( data.size() ), dataSizes, MPI_COMM_GEOS );
  // `totalDataSize` contains the total data size across all the MPI ranks.
  int const totalDataSize = std::accumulate( dataSizes.begin(), dataSizes.end(), 0 );

  // Once the MPI exchange is done, `allData` will contain all the data of all the MPI ranks.
  // We want all ranks to get all the data. But each rank may have a different size of information.
  // Therefore, we use `allgatherv` that does not impose the same size across ranks like `allgather` does.
  stdVector< T > allData( totalDataSize );
  // `displacements` is the offset (relative to the receive buffer) to store the data for each rank.
  stdVector< int > displacements( MpiWrapper::commSize(), 0 );
  std::partial_sum( dataSizes.begin(), dataSizes.end() - 1, displacements.begin() + 1 );
  MpiWrapper::allgatherv( data.data(), data.size(), allData.data(), dataSizes.data(), displacements.data(), MPI_COMM_GEOS );

  // Finalizing by sorting, removing duplicates and trimming the result vector at the proper size.
  std::sort( allData.begin(), allData.end() );
  auto newEnd = std::unique( allData.begin(), allData.end() );
  allData.erase( newEnd, allData.end() );

  return allData;
}

/**
 * @brief Check it the vtk grid is a (supported) structured mesh
 * @param[in] mesh a vtk grid
 * @return @p true if @p mesh is structured; @p false otherwise
 */
bool isMeshStructured( vtkSmartPointer< vtkDataSet > mesh )
{
  return mesh->IsA( "vtkStructuredPoints" )
         || mesh->IsA( "vtkStructuredGrid" )
         || mesh->IsA( "vtkRectilinearGrid" )
         || mesh->IsA( "vtkImageData" );
}


/**
 * @brief Generate global point and cell ids
 *
 * @param[in] mesh a vtk grid
 * @return the vtk grid with global ids attributes
 */
vtkSmartPointer< vtkDataSet >
generateGlobalIDs( vtkSmartPointer< vtkDataSet > mesh )
{
  GEOS_MARK_FUNCTION;

  vtkNew< vtkGenerateGlobalIds > generator;
  generator->SetInputDataObject( mesh );
  generator->Update();
  return vtkDataSet::SafeDownCast( generator->GetOutputDataObject( 0 ) );
}


/**
 * @brief Build the element to nodes mappings for the 3D mesh only
 * @tparam INDEX_TYPE The indexing type
 * @param mesh The 3D mesh
 * @return The mapping for 3D cells only (fractures excluded)
 * @note Fractures are partitioned separately based on 3D neighbors, not included in ParMETIS graph
 */
template< typename INDEX_TYPE >
ArrayOfArrays< INDEX_TYPE, INDEX_TYPE >
buildElemToNodes( vtkSmartPointer< vtkDataSet > mesh )
{
  GEOS_MARK_FUNCTION;

  localIndex const numCells = LvArray::integerConversion< localIndex >( mesh->GetNumberOfCells() );

  array1d< INDEX_TYPE > nodeCounts( numCells );

  forAll< parallelHostPolicy >( numCells, [&mesh, nodeCounts = nodeCounts.toView()] ( localIndex const cellIdx )
  {
    nodeCounts[cellIdx] = LvArray::integerConversion< INDEX_TYPE >( mesh->GetCellSize( cellIdx ) );
  } );

  ArrayOfArrays< INDEX_TYPE, INDEX_TYPE > elemToNodes;
  elemToNodes.template resizeFromCapacities< parallelHostPolicy >( numCells, nodeCounts.data() );

  vtkIdTypeArray const & globalPointId = *vtkIdTypeArray::FastDownCast( mesh->GetPointData()->GetGlobalIds() );

  forAll< parallelHostPolicy >( numCells, [&mesh, &globalPointId, elemToNodes = elemToNodes.toView()] ( localIndex const cellIdx )
  {
    vtkSmartPointer< vtkIdList > const ptIds = vtkSmartPointer< vtkIdList >::New();
    mesh->GetCellPoints( cellIdx, ptIds );
    vtkIdType const numPts = ptIds->GetNumberOfIds();
    for( vtkIdType a = 0; a < numPts; ++a )
    {
      vtkIdType const pointIdx = globalPointId.GetValue( ptIds->GetId( a ) );
      elemToNodes.emplaceBack( cellIdx, LvArray::integerConversion< INDEX_TYPE >( pointIdx ) );
    }
  } );

  return elemToNodes;
}


/**
 * @brief Split a mesh by partitioning it
 *
 * @tparam PART_INDEX the type of the partition indexes
 * @param mesh an input mesh (can be empty)
 * @param numParts The number of partitions to create
 * @param part an array of target partitions for each element in local mesh (can be empty)
 * @return A collection of vtkUnstructuredGrid, one for each target partition/rank.
 * If no data is to be split for a given partition, then the returned vtkUnstructuredGrid for this rank will be empty.
 * @details Splits @p mesh according to the @p part information. Each target partition gets separated into
 * its own vtkUnstructuredGrid returned as part of the vtkPartitionedDataSet.
 *
 * @note This function handles empty inputs gracefully:
 * - If @p part is empty, returns @p numParts empty partitions
 * - If @p mesh has no cells, returns @p numParts empty partitions
 */
template< typename PART_INDEX >
vtkSmartPointer< vtkPartitionedDataSet >
splitMeshByPartition( vtkSmartPointer< vtkDataSet > mesh,
                      int const numParts,
                      arrayView1d< PART_INDEX const > const & part )
{
  vtkNew< vtkPartitionedDataSet > result;
  result->SetNumberOfPartitions( LvArray::integerConversion< unsigned int >( numParts ) );

  // Handle empty input: create empty partitions for all ranks
  if( part.empty() || mesh->GetNumberOfCells() == 0 )
  {
    for( int p = 0; p < numParts; ++p )
    {
      result->SetPartition( LvArray::integerConversion< unsigned int >( p ),
                            vtkNew< vtkUnstructuredGrid >() );
    }
    return result;
  }

  // Count cells per partition
  array1d< localIndex > cellCounts( numParts );
  forAll< parallelHostPolicy >( part.size(), [part, cellCounts = cellCounts.toView()] ( localIndex const cellIdx )
  {
    RAJA::atomicInc< parallelHostAtomic >( &cellCounts[LvArray::integerConversion< localIndex >( part[cellIdx] )] );
  } );

  // Build cell lists per partition
  ArrayOfArrays< vtkIdType > cellsLists;
  cellsLists.resizeFromCapacities< serialPolicy >( numParts, cellCounts.data() );

  forAll< parallelHostPolicy >( part.size(), [part, cellsLists = cellsLists.toView()] ( localIndex const cellIdx )
  {
    cellsLists.emplaceBackAtomic< parallelHostAtomic >( LvArray::integerConversion< localIndex >( part[cellIdx] ),
                                                        LvArray::integerConversion< vtkIdType >( cellIdx ) );
  } );

  // Extract cells for each partition
  vtkNew< vtkExtractCells > extractor;
  extractor->SetInputDataObject( mesh );

  for( localIndex p = 0; p < numParts; ++p )
  {
    arraySlice1d< vtkIdType const > const cells = cellsLists[p];

    if( cells.size() > 0 )
    {
      extractor->SetCellIds( cells.dataIfContiguous(), LvArray::integerConversion< vtkIdType >( cells.size() ) );
      extractor->Update();

      vtkNew< vtkUnstructuredGrid > ug;
      ug->ShallowCopy( extractor->GetOutputDataObject( 0 ) );
      result->SetPartition( LvArray::integerConversion< unsigned int >( p ), ug );
    }
    else
    {
      // No cells for this partition - set empty grid
      result->SetPartition( LvArray::integerConversion< unsigned int >( p ),
                            vtkNew< vtkUnstructuredGrid >() );
    }
  }

  return result;
}

vtkSmartPointer< vtkMultiProcessController > getController()
{
#ifdef GEOS_USE_MPI
  vtkNew< vtkMPIController > controller;
  vtkMPICommunicatorOpaqueComm vtkGeosxComm( &MPI_COMM_GEOS );
  vtkNew< vtkMPICommunicator > communicator;
  communicator->InitializeExternal( &vtkGeosxComm );
  controller->SetCommunicator( communicator );
#else
  vtkNew< vtkDummyController > controller;
#endif
  return controller;
}

VTKLegacyDatasetType getVTKLegacyDatasetType( vtkSmartPointer< vtkDataSetReader > const vtkGridReader,
                                              Path const & filePath )
{
  vtkGridReader->SetFileName( filePath.c_str() );

  if( vtkGridReader->IsFileStructuredPoints() )
  {
    return VTKLegacyDatasetType::structuredPoints;
  }
  else if( vtkGridReader->IsFileStructuredGrid() )
  {
    return VTKLegacyDatasetType::structuredGrid;
  }
  else if( vtkGridReader->IsFileUnstructuredGrid() )
  {
    return VTKLegacyDatasetType::unstructuredGrid;
  }
  else if( vtkGridReader->IsFileRectilinearGrid() )
  {
    return VTKLegacyDatasetType::rectilinearGrid;
  }
  else if( vtkGridReader->IsFilePolyData())
  {
    return VTKLegacyDatasetType::polyData;
  }
  else
  {
    GEOS_ERROR( GEOS_FMT( "Unsupported legacy VTK dataset format.\nLegacy supported formats are: {}.",
                          EnumStrings< VTKLegacyDatasetType >::concat( ", " ) ) );
  }
  return {};
}

vtkSmartPointer< vtkDataSet >
loadMesh( Path const & filePath,
          string const & blockName,
          int readerRank = 0 )
{
  string const extension = filePath.extension();

  auto const parallelRead = [&]( auto const vtkGridReader )
  {
    using GridType = TYPEOFPTR( vtkGridReader->GetOutput() );
    vtkGridReader->SetFileName( filePath.c_str() );
    vtkGridReader->UpdateInformation();
    vtkGridReader->UpdatePiece( MpiWrapper::commRank(), MpiWrapper::commSize(), 0 );
    return vtkSmartPointer< GridType >( vtkGridReader->GetOutput() );
  };

  auto const serialRead = [&]( auto const vtkGridReader )
  {
    using GridType = TYPEOFPTR( vtkGridReader->GetOutput() );
    if( MpiWrapper::commRank() == readerRank )
    {
      vtkGridReader->SetFileName( filePath.c_str() );
      vtkGridReader->Update();
      return vtkSmartPointer< GridType >( vtkGridReader->GetOutput() );
    }
    else
    {
      return vtkSmartPointer< GridType >::New();
    }
  };

  switch( EnumStrings< VTKMeshExtension >::fromString( extension ) )
  {
    case VTKMeshExtension::vtm:
    {
      if( MpiWrapper::commRank() == readerRank )
      {
        // The multi-block format is a container of multiple datasets (or even of other containers).
        // We must navigate this multi-block to extract the relevant information.
        auto reader = vtkSmartPointer< vtkXMLMultiBlockDataReader >::New();
        reader->SetFileName( filePath.c_str() );
        reader->Update();
        vtkCompositeDataSet * compositeDataSet = reader->GetOutput();
        if( !compositeDataSet->IsA( "vtkMultiBlockDataSet" ) )
        {
          GEOS_ERROR( GEOS_FMT( "Unsupported vtk multi-block format in file \"{}\".\n{}",
                                filePath,
                                generalMeshErrorAdvice ) );
        }
        vtkMultiBlockDataSet * multiBlockDataSet = vtkMultiBlockDataSet::SafeDownCast( compositeDataSet );

        // Looking for the _first_ block that matches the requested name.
        // No check is performed to validate that there is not name duplication.
        for( unsigned int i = 0; i < multiBlockDataSet->GetNumberOfBlocks(); ++i )
        {
          string const dataSetName = multiBlockDataSet->GetMetaData( i )->Get( multiBlockDataSet->NAME() );
          if( dataSetName == blockName )
          {
            vtkDataObject * block = multiBlockDataSet->GetBlock( i );
            if( block->IsA( "vtkDataSet" ) )
            {
              vtkSmartPointer< vtkDataSet > mesh = vtkDataSet::SafeDownCast( block );
              return mesh;
            }
          }
        }
        GEOS_ERROR( GEOS_FMT( "Could not find mesh \"{}\" in multi-block vtk file \"{}\".\n{}",
                              blockName,
                              filePath,
                              generalMeshErrorAdvice ) );
        return {};
      }
      else
      {
        return vtkSmartPointer< vtkUnstructuredGrid >::New();
      }
    }
    case VTKMeshExtension::vtk:
    {
      VTKLegacyDatasetType const datasetType = getVTKLegacyDatasetType( vtkSmartPointer< vtkDataSetReader >::New(),
                                                                        filePath );
      switch( datasetType )
      {
        case VTKLegacyDatasetType::structuredPoints: return serialRead( vtkSmartPointer< vtkStructuredPointsReader >::New() );
        case VTKLegacyDatasetType::structuredGrid:   return serialRead( vtkSmartPointer< vtkStructuredGridReader >::New() );
        case VTKLegacyDatasetType::unstructuredGrid: return serialRead( vtkSmartPointer< vtkUnstructuredGridReader >::New() );
        case VTKLegacyDatasetType::rectilinearGrid:  return serialRead( vtkSmartPointer< vtkRectilinearGridReader >::New() );
        case VTKLegacyDatasetType::polyData:         return serialRead( vtkSmartPointer< vtkPolyDataReader >::New() );
      }
      break;
    }
    case VTKMeshExtension::vtu: return serialRead( vtkSmartPointer< vtkXMLUnstructuredGridReader >::New() );
    case VTKMeshExtension::vtr: return serialRead( vtkSmartPointer< vtkXMLRectilinearGridReader >::New() );
    case VTKMeshExtension::vts: return serialRead( vtkSmartPointer< vtkXMLStructuredGridReader >::New() );
    case VTKMeshExtension::vti: return serialRead( vtkSmartPointer< vtkXMLImageDataReader >::New() );
    case VTKMeshExtension::vtp: return serialRead( vtkSmartPointer< vtkXMLPolyDataReader >::New() );
    case VTKMeshExtension::pvtu:
    {
      return parallelRead( vtkSmartPointer< vtkXMLPUnstructuredGridReader >::New() );
      // TODO: Apply vtkStaticCleanUnstructuredGrid, once it is included in the next VTK release.
      //       https://gitlab.kitware.com/vtk/vtk/-/blob/master/Filters/Core/vtkStaticCleanUnstructuredGrid.h
      //       This removes duplicate points, either present in the dataset, or resulting from merging pieces.
    }
    case VTKMeshExtension::pvts: return parallelRead( vtkSmartPointer< vtkXMLPStructuredGridReader >::New() );
    case VTKMeshExtension::pvtr: return parallelRead( vtkSmartPointer< vtkXMLPRectilinearGridReader >::New() );
    case VTKMeshExtension::pvti: return parallelRead( vtkSmartPointer< vtkXMLPImageDataReader >::New() );
    case VTKMeshExtension::pvtp: return parallelRead( vtkSmartPointer< vtkXMLPPolyDataReader >::New() );
    default:
    {
      GEOS_ERROR( GEOS_FMT( "{} is not a recognized extension for VTKMesh. Please use .{}",
                            extension,
                            EnumStrings< VTKMeshExtension >::concat( ", ." ) ) );
      break;
    }
  }

  return {};
}


AllMeshes loadAllMeshes( Path const & filePath,
                         string const & mainBlockName,
                         string_array const & faceBlockNames )
{
  vtkSmartPointer< vtkDataSet > main = loadMesh( filePath, mainBlockName );
  stdMap< string, vtkSmartPointer< vtkDataSet > > faces;

  // Load fractures on rank 0 (same as main mesh for 2D cells)
  // This allows building fracture-to-3D connectivity before redistribution
  for( string const & faceBlockName: faceBlockNames )
  {
    faces.insert( { faceBlockName, loadMesh( filePath, faceBlockName, 0 ) } );
  }

  return AllMeshes( main, faces );
}


/**
 * @brief Partition the mesh using cell graph methods (ParMETIS or PTScotch)
 *
 * @param[in] input a collection of vtk main (3D) and fracture mesh
 * @param[in] method the partitioning method
 * @param[in] comm the MPI communicator
 * @param[in] numParts the number of partitions
 * @param[in] minCommonNodes the minimum number of shared nodes for adding a graph edge
 * @param[in] numRefinements the number of refinements for PTScotch
 * @return the cell partitioning array
 */
array1d< int64_t >
partitionByCellGraph( vtkSmartPointer< vtkDataSet > mesh3D,
                      PartitionMethod const method,
                      MPI_Comm const comm,
                      int const numParts,
                      int const minCommonNodes,
                      int const numRefinements )
{
  GEOS_MARK_FUNCTION;

  pmet_idx_t const numElems = mesh3D->GetNumberOfCells();
  pmet_idx_t const numRanks = MpiWrapper::commSize( comm );

  // Note: elemDist contains 3D cells only. Fractures are assigned later based on 3D neighbors.
  array1d< pmet_idx_t > const elemDist( numRanks + 1 );
  {
    array1d< pmet_idx_t > elemCounts;
    MpiWrapper::allGather( numElems, elemCounts, comm );
    std::partial_sum( elemCounts.begin(), elemCounts.end(), elemDist.begin() + 1 );
  }

  // Build element-to-node connectivity for the 3D mesh only
  ArrayOfArrays< pmet_idx_t, pmet_idx_t > const elemToNodes = buildElemToNodes< pmet_idx_t >( mesh3D );
  ArrayOfArrays< pmet_idx_t, pmet_idx_t > graph;
#ifdef GEOS_USE_PARMETIS
  graph = parmetis::meshToDual( elemToNodes.toViewConst(), elemDist, comm, minCommonNodes );
#else
  GEOS_THROW( "GEOS must be built with ParMETIS support (ENABLE_PARMETIS=ON)"
              "to use any graph partitioning method for parallel mesh distribution", InputError );
#endif

  switch( method )
  {
    case PartitionMethod::parmetis:
    {
#ifdef GEOS_USE_PARMETIS
      return parmetis::partition( graph.toViewConst(), elemDist, numParts, comm, numRefinements );
#else
      GEOS_THROW( "GEOS must be built with ParMETIS support (ENABLE_PARMETIS=ON) to use 'parmetis' partitioning method", InputError );
      return {};
#endif
    }
    case PartitionMethod::ptscotch:
    {
#ifdef GEOS_USE_SCOTCH
      GEOS_WARNING_IF( numRefinements > 0, "Partition refinement is not supported by 'ptscotch' partitioning method" );
      return ptscotch::partition( graph.toViewConst(), numParts, comm );
#else
      GEOS_THROW( "GEOS must be built with Scotch support (ENABLE_SCOTCH=ON) to use 'ptscotch' partitioning method", InputError );
      return {};
#endif
    }
    default:
    {
      GEOS_THROW( "Unknown partition method", InputError );
      return{};
    }
  }
  return {};
}

/**
 * @brief Redistribute the 3D mesh using cell graph methods (ParMETIS or PTScotch)
 * @param[in] mesh3D The 3D main mesh to partition
 * @param[in] method The partitioning method
 * @param[in] comm The MPI communicator
 * @param[in] numRefinements The number of refinement iterations
 * @return Redistributed 3D mesh
 * @note Fractures are NOT partitioned here - they will be assigned later
 */
vtkSmartPointer< vtkDataSet >
redistributeByCellGraph( vtkSmartPointer< vtkDataSet > mesh3D,
                         PartitionMethod const method,
                         MPI_Comm const comm,
                         int const numRefinements )
{
  GEOS_MARK_FUNCTION;

  int const numRanks = MpiWrapper::commSize( comm );

  // Partition the 3D main mesh only
  array1d< int64_t > newPartitions = partitionByCellGraph( mesh3D, method, comm, numRanks, 3, numRefinements );

  // Split and redistribute the 3D mesh
  vtkSmartPointer< vtkPartitionedDataSet > const splitMesh =
    splitMeshByPartition( mesh3D, numRanks, newPartitions.toViewConst() );

  return vtk::redistribute( *splitMesh, comm );
}


/**
 * @brief Partition mesh while keeping super-cells (cell groups) intact
 *
 * Builds a weighted graph where vertices represent super-cells, then partitions
 * using ParMETIS/PTScotch. Ensures cells with the same SuperCellId stay together.
 *
 * @param mesh Input mesh with "SuperCellId" cell data array
 * @param method Partitioning algorithm (parmetis or ptscotch)
 * @param comm MPI communicator
 * @param numRefinementIterations Number of ParMETIS refinement passes
 * @return Redistributed mesh with super-cells kept intact
 */
vtkSmartPointer< vtkDataSet >
redistributeBySuperCellGraph(
  vtkSmartPointer< vtkDataSet > mesh,
  PartitionMethod const method,
  MPI_Comm comm,
  int const numRefinementIterations,
  int const fractureWeight )
{
  GEOS_MARK_FUNCTION;

  int const rank = MpiWrapper::commRank( comm );
  int const numRanks = MpiWrapper::commSize( comm );

  // -----------------------------------------------------------------------
  // Step 1: Build base cell graph (standard adjacency)
  // -----------------------------------------------------------------------
  vtkSmartPointer< vtkUnstructuredGrid > ugrid =
    vtkUnstructuredGrid::SafeDownCast( mesh );

  GEOS_ERROR_IF( !ugrid, "Mesh is not vtkUnstructuredGrid" );

  ArrayOfArrays< pmet_idx_t, pmet_idx_t > baseCellGraph;
  array1d< pmet_idx_t > baseElemDist( numRanks + 1 );
  {
    // Build element distribution based on local cell counts
    pmet_idx_t const numLocalCells = ugrid->GetNumberOfCells();

    array1d< pmet_idx_t > cellCounts;
    MpiWrapper::allGather( numLocalCells, cellCounts, comm );

    baseElemDist[0] = 0;
    std::partial_sum( cellCounts.begin(), cellCounts.end(), baseElemDist.begin() + 1 );

    // Build element-to-node connectivity local to each rank
    ArrayOfArrays< pmet_idx_t, pmet_idx_t > const elemToNodes =
      buildElemToNodes< pmet_idx_t >( ugrid );

    // Build dual graph (cell-to-cell via shared nodes)
#ifdef GEOS_USE_PARMETIS
    int const minCommonNodes = 3; // Minimum shared nodes for edge
    baseCellGraph = parmetis::meshToDual( elemToNodes.toViewConst(), baseElemDist, comm, minCommonNodes );
#else
    GEOS_THROW( "GEOS must be built with ParMETIS support (ENABLE_PARMETIS=ON) "
                "to build cell graphs for partitioning", InputError );
#endif
  }

  // -----------------------------------------------------------------------
  // Step 2: Reconstruct super-cell info on all ranks (from SuperCellId array)
  // -----------------------------------------------------------------------
  SuperCellInfo localSuperCellInfo = reconstructSuperCellInfo( ugrid, fractureWeight );

  // -----------------------------------------------------------------------
  // Step 3: Build super-cell graph
  // -----------------------------------------------------------------------
  auto [superCellGraph, superVertexWeights] = buildSuperCellGraph(
    ugrid,
    baseCellGraph,
    baseElemDist,
    localSuperCellInfo,
    comm
    );

  // -----------------------------------------------------------------------
  // Step 4: Compute super-cell element distribution
  // -----------------------------------------------------------------------
  array1d< pmet_idx_t > superElemDist( numRanks + 1 );
  {
    pmet_idx_t const localSuperCellCount = LvArray::integerConversion< pmet_idx_t >( superCellGraph.size() );

    array1d< pmet_idx_t > superCellCounts;
    MpiWrapper::allGather( localSuperCellCount, superCellCounts, comm );

    superElemDist[0] = 0;
    std::partial_sum( superCellCounts.begin(), superCellCounts.end(),
                      superElemDist.begin() + 1 );
  }

  // -----------------------------------------------------------------------
  // Step 5: Validate graph before partitioning
  // -----------------------------------------------------------------------
  validateSuperCellGraph(
    superCellGraph,
    superElemDist.toViewConst(),
    superVertexWeights.toViewConst(),
    comm
    );

  // -----------------------------------------------------------------------
  // Step 6: Partition super-cell graph using ParMETIS/PTScotch
  // -----------------------------------------------------------------------
  array1d< int64_t > superCellPartitioning;

  if( method == PartitionMethod::parmetis )
  {
    GEOS_LOG_RANK_0( "Partitioning super-cell graph with ParMETIS..." );

    // Basic validation
    pmet_idx_t const minGraphSize = MpiWrapper::min( static_cast< pmet_idx_t >( superCellGraph.size() ), comm );
    GEOS_ERROR_IF( minGraphSize == 0, "At least one rank has an empty super-cell graph!" );

    superCellPartitioning = parmetis::partitionWeighted(
      superCellGraph.toViewConst(),
      superVertexWeights.toViewConst(),
      superElemDist.toViewConst(),
      numRanks,
      comm,
      numRefinementIterations
      );
  }
  else
  {
    GEOS_ERROR( GEOS_FMT( "Unsupported partition method: {}",
                          EnumStrings< PartitionMethod >::toString( method ) ) );
  }

  // -----------------------------------------------------------------------
  // Step 7: Map SuperCellId to local super-cell indices for partition unpacking
  // -----------------------------------------------------------------------

  vtkIdTypeArray * superCellIdArray =
    vtkIdTypeArray::SafeDownCast( ugrid->GetCellData()->GetArray( "SuperCellId" ) );

  GEOS_ERROR_IF( !superCellIdArray,
                 GEOS_FMT( "SuperCellId array not found on rank {}", rank ) );

  stdMap< vtkIdType, localIndex > superCellIdToLocalIdx;

// Build ordered list of local super-cell IDs
  stdVector< vtkIdType > orderedSuperCellIds;
  orderedSuperCellIds.reserve( localSuperCellInfo.superCellToOriginalCells.size() );

  for( auto const & [scId, cells] : localSuperCellInfo.superCellToOriginalCells )
  {
    orderedSuperCellIds.push_back( scId );
  }

  // Sort to ensure consistent ordering
  std::sort( orderedSuperCellIds.begin(), orderedSuperCellIds.end() );

  localIndex const numLocalSuperCells = LvArray::integerConversion< localIndex >( orderedSuperCellIds.size() );
  for( localIndex i = 0; i < numLocalSuperCells; ++i )
  {
    superCellIdToLocalIdx.insert( { orderedSuperCellIds[i], i } );
  }

  // -----------------------------------------------------------------------
  // Step 8: Unpack super-cell partitioning to individual cells
  // -----------------------------------------------------------------------
  array1d< int64_t > cellPartitioning = expandSuperCellPartitioningToCells(
    ugrid,
    superCellPartitioning,
    superCellIdToLocalIdx,
    comm
    );

  // -----------------------------------------------------------------------
  // Step 9: Verify no super-cells were split across ranks
  // -----------------------------------------------------------------------
  stdMap< vtkIdType, int64_t > superCellToRank;

  vtkIdType const numCells = ugrid->GetNumberOfCells();
  for( vtkIdType i = 0; i < numCells; ++i )
  {
    vtkIdType const scId = superCellIdArray->GetValue( i );
    int64_t const targetRank = cellPartitioning[i];

    auto [it, inserted] = superCellToRank.insert( {scId, targetRank} );

    GEOS_ERROR_IF( !inserted && it->second != targetRank,
                   GEOS_FMT( "Super-cell {} is split: assigned to both rank {} and rank {}",
                             scId, it->second, targetRank ) );
  }

  // -----------------------------------------------------------------------
  // Step 10: Redistribute mesh according to cell partitioning
  // -----------------------------------------------------------------------
  // Split mesh according to partitioning
  vtkSmartPointer< vtkPartitionedDataSet > splitMesh =
    splitMeshByPartition( ugrid, numRanks, cellPartitioning.toViewConst() );

  // Redistribute using VTK
  vtkSmartPointer< vtkDataSet > redistributed = vtk::redistribute( *splitMesh, comm );

  // -----------------------------------------------------------------------
  // Step 11: Report final distribution statistics
  // -----------------------------------------------------------------------
  array1d< vtkIdType > cellsPerRank;
  vtkIdType const localCells = redistributed->GetNumberOfCells();
  MpiWrapper::allGather( localCells, cellsPerRank, comm );

  if( rank == 0 )
  {
    vtkIdType totalCells = 0;
    vtkIdType minCells = std::numeric_limits< vtkIdType >::max();
    vtkIdType maxCells = 0;

    for( int r = 0; r < numRanks; ++r )
    {
      totalCells += cellsPerRank[r];
      minCells = std::min( minCells, cellsPerRank[r] );
      maxCells = std::max( maxCells, cellsPerRank[r] );
    }

    real64 const avgCells = static_cast< real64 >( totalCells ) / numRanks;
    real64 const imbalance = (numRanks > 1) ? (maxCells - avgCells) / avgCells * 100.0 : 0.0;

    GEOS_LOG_RANK_0( GEOS_FMT( "Final 3D cells redistribution: {} total, {} avg/rank, [{}-{}] range, {:.1f}% imbalance",
                               totalCells, static_cast< vtkIdType >(avgCells),
                               minCells, maxCells, imbalance ) );
  }

  return redistributed;
}


/**
 * @brief Scatter the mesh by blocks  (no geometric information involved, assumes rank 0 has the full mesh)
 *
 * @param[in] mesh a vtk grid
 * @return the vtk grid redistributed
 */
vtkSmartPointer< vtkDataSet >
scatterByBlock( vtkDataSet & mesh )
{
  GEOS_MARK_FUNCTION;

  int const rank = MpiWrapper::commRank();
  int const size = MpiWrapper::commSize();

  // Count total cells across all ranks
  vtkIdType localCells = mesh.GetNumberOfCells();
  vtkIdType totalCells = MpiWrapper::allReduce( localCells, MpiWrapper::Reduction::Sum, MPI_COMM_GEOS );

  // Handle edge cases
  if( totalCells == 0 )
  {
    vtkNew< vtkUnstructuredGrid > emptyMesh;
    return emptyMesh;
  }

  if( size == 1 )
  {
    vtkNew< vtkUnstructuredGrid > copy;
    copy->DeepCopy( &mesh );
    return copy;
  }

  // Verify rank 0 has the complete mesh for redistribution
  if( rank == 0 && localCells != totalCells )
  {
    GEOS_ERROR( GEOS_FMT( "Rank 0 must have the complete mesh. Rank 0 has {} cells but total is {}",
                          localCells,
                          totalCells ) );
  }

  // Scatter cells by contiguous blocks
  vtkIdType cellsPerRank = totalCells / size;
  vtkIdType remainder = totalCells % size;

  // Create partitioned dataset
  vtkNew< vtkPartitionedDataSet > localParts;
  if( rank == 0 )
  {
    // Rank 0 has the full mesh, extract cells for each rank
    for( int r = 0; r < size; ++r )
    {
      vtkIdType rankStart = r * cellsPerRank + std::min( (vtkIdType)r, remainder );
      vtkIdType rankEnd = rankStart + cellsPerRank + (r < remainder ? 1 : 0);

      // Validate cell range
      GEOS_ERROR_IF( rankStart< 0 || rankEnd > totalCells,
                     GEOS_FMT( "Invalid cell range for rank {}: [{}, {}) with total cells {}",
                               r,
                               rankStart,
                               rankEnd,
                               totalCells ) );

      if( rankEnd > rankStart )
      {
        // Add cells for this rank
        vtkNew< vtkExtractCells > extractor;
        extractor->SetInputDataObject( &mesh );
        extractor->AddCellRange( rankStart, rankEnd - 1 );
        extractor->Update();
        vtkUnstructuredGrid * extracted = extractor->GetOutput();
        localParts->SetPartition( r, extracted );
      }
      else
      {
        // Create empty partition for ranks with no cells
        vtkNew< vtkUnstructuredGrid > emptyPartition;
        localParts->SetPartition( r, emptyPartition );
      }
    }
  }
  else
  {
    // Other ranks have an empty mesh, but we still need to create the
    // partitioned data set structure.
    localParts->SetNumberOfPartitions( size );
    for( int r = 0; r < size; ++r )
    {
      vtkNew< vtkUnstructuredGrid > emptyPartition;
      localParts->SetPartition( r, emptyPartition );
    }
  }

  //Send cells to appropriate ranks
  vtkSmartPointer< vtkUnstructuredGrid > result = vtk::redistribute( *localParts, MPI_COMM_GEOS );

  // Final validation
  vtkIdType finalLocalCells = result->GetNumberOfCells();
  vtkIdType finalTotalCells = MpiWrapper::allReduce( finalLocalCells, MpiWrapper::Reduction::Sum, MPI_COMM_GEOS );

  GEOS_ERROR_IF( finalTotalCells != totalCells,
                 GEOS_FMT( "Block redistribution lost cells: started with {}, ended with {}",
                           totalCells,
                           finalTotalCells ) );

  return result;
}


/**
 * @brief Classify cells by dimension
 *
 * Builds index arrays for 2D and 3D cells while ignoring lower-dimensional elements.
 *
 * @param[in] mesh Input mesh containing cells of mixed dimensions
 * @param[out] cells3DIndices Indices of 3D cells in original mesh
 * @param[out] cells2DIndices Indices of 2D cells in original mesh
 */
static void classifyCellsByDimension( vtkDataSet & mesh,
                                      array1d< vtkIdType > & cells3DIndices,
                                      array1d< vtkIdType > & cells2DIndices )
{
  GEOS_MARK_FUNCTION;

  vtkIdType const numCells = mesh.GetNumberOfCells();

  cells3DIndices.clear();
  cells2DIndices.clear();

  // Single pass: classify and populate
  for( vtkIdType i = 0; i < numCells; ++i )
  {
    int const dim = vtkCellTypes::GetDimension( mesh.GetCellType( i ) );
    if( dim == 3 )
    {
      cells3DIndices.emplace_back( i );
    }
    else if( dim == 2 )
    {
      cells2DIndices.emplace_back( i );
    }
  }

  GEOS_LOG_RANK_0( GEOS_FMT( "Classified mesh: {} 3D cells, {} 2D cells (from {} total)",
                             cells3DIndices.size(), cells2DIndices.size(), numCells ) );
}

/**
 * @brief Build mapping from 2D cells to their neighboring 3D cells using indices
 *
 * @param[in] mesh Original mesh
 * @param[in] cells2DIndices Indices of 2D cells in original mesh
 * @param[in] cells3DIndices Indices of 3D cells in original mesh
 * @return Mapping from 2D cell index (in cells2DIndices) to global IDs of neighboring 3D cells
 */
static ArrayOfArrays< localIndex, int64_t >
build2DTo3DNeighbors( vtkDataSet & mesh,
                      arrayView1d< vtkIdType const > cells2DIndices,
                      arrayView1d< vtkIdType const > cells3DIndices )
{
  GEOS_MARK_FUNCTION;

  // Retrieve global cell ID array
  vtkDataArray * globalCellIds = mesh.GetCellData()->GetGlobalIds();
  GEOS_ERROR_IF( globalCellIds == nullptr,
                 "Global cell IDs must be present in mesh for 2D-3D neighbor mapping" );

  // Build reverse lookup: original mesh index to global cell ID (3D cells only)
  stdUnorderedMap< vtkIdType, int64_t > meshIdxToGlobalId3D;
  meshIdxToGlobalId3D.reserve( cells3DIndices.size() );

  for( vtkIdType meshIdx : cells3DIndices )
  {
    int64_t const globalId = static_cast< int64_t >( globalCellIds->GetTuple1( meshIdx ) );
    meshIdxToGlobalId3D.emplace( meshIdx, globalId );
  }

  ArrayOfArrays< localIndex, int64_t > neighbors2Dto3D;
  neighbors2Dto3D.reserve( cells2DIndices.size() );

  // Topology statistics
  localIndex numStandalone = 0;
  localIndex numBoundary = 0;       // 1 neighbor
  localIndex numInternal = 0;       // 2 neighbors
  localIndex numJunction = 0;       // >2 neighbors

  vtkNew< vtkIdList > neighborCells;
  vtkNew< vtkIdList > pointIds2D;

  // Build neighbor list for each 2D cell
  for( vtkIdType meshIdx2D : cells2DIndices )
  {
    mesh.GetCellPoints( meshIdx2D, pointIds2D );

    // Find all cells sharing ALL nodes with this 2D cell (exact face match)
    neighborCells->Reset();
    mesh.GetCellNeighbors( meshIdx2D, pointIds2D, neighborCells );

    // Filter for 3D neighbors and retrieve their global IDs
    array1d< int64_t > neighbor3DGlobalIds;
    neighbor3DGlobalIds.reserve( neighborCells->GetNumberOfIds() );

    for( vtkIdType n = 0; n < neighborCells->GetNumberOfIds(); ++n )
    {
      vtkIdType const neighborIdx = neighborCells->GetId( n );

      // Check if neighbor is a 3D cell and get its global ID
      auto it = meshIdxToGlobalId3D.find( neighborIdx );
      if( it != meshIdxToGlobalId3D.end() )
      {
        neighbor3DGlobalIds.emplace_back( it->second );
      }
      // Non-3D neighbors (2D/1D/0D) are silently skipped
    }

    // Update topology statistics
    localIndex const numNeighbors = neighbor3DGlobalIds.size();
    switch( numNeighbors )
    {
      case 0:  numStandalone++; break;
      case 1:  numBoundary++;   break;
      case 2:  numInternal++;   break;
      default: numJunction++;   break;
    }

    neighbors2Dto3D.appendArray( neighbor3DGlobalIds.begin(), neighbor3DGlobalIds.end() );
  }

  // Print diagnostic summary
  if( numStandalone > 0 || numJunction > 0 )
  {
    // Show detailed breakdown when anomalies exist
    GEOS_LOG_RANK_0( "\n2D-to-3D Neighbor Topology" );
    GEOS_LOG_RANK_0( GEOS_FMT( " Total 2D cells:           {}", cells2DIndices.size() ) );
    GEOS_LOG_RANK_0( GEOS_FMT( " Standalone (0 neighbors): {}", numStandalone ) );
    GEOS_LOG_RANK_0( GEOS_FMT( " Boundary   (1 neighbor):  {}", numBoundary ) );
    GEOS_LOG_RANK_0( GEOS_FMT( " Internal   (2 neighbors): {}", numInternal ) );
    GEOS_LOG_RANK_0( GEOS_FMT( " Junction   (>2 neighbors):{}\n", numJunction ) );
  }
  else
  {
    // Condensed output for normal cases
    GEOS_LOG_RANK_0( GEOS_FMT( "2D cells: {} (1-neighbor/boundary: {}, 2-neighbor/internal: {})",
                               cells2DIndices.size(), numBoundary, numInternal ) );
  }

  // Standalone 2D cells indicate mesh topology errors
  GEOS_ERROR_IF( numStandalone > 0,
                 GEOS_FMT( "{} orphaned 2D cells detected with no 3D neighbors. "
                           "These may be artifacts or detached surfaces. "
                           "Please clean the mesh or verify the geometry.", numStandalone ) );

  return neighbors2Dto3D;
}

/**
 * @brief Assign cells to partitions based on their 3D neighbor locations
 *
 * Generic function for assigning 2D cells or fracture elements to partitions.
 * Uses deterministic tie-breaking (minimum neighbor global ID) to ensure
 * each cell is co-located with at least one of its 3D neighbors.
 *
 * @param[in] neighborsTo3D Neighbor mapping (cell index -> 3D global IDs)
 * @param[in] local3DGlobalIds Global IDs of 3D cells on this rank
 * @param[in] local3DPartitions Partition assignments for local 3D cells
 * @param[in] comm MPI communicator
 * @return Partition assignments for each input cell
 */
static array1d< int >
assignCellsBasedOn3DNeighbors( ArrayOfArrays< localIndex, int64_t > const & neighbors2Dto3D,
                               arrayView1d< int64_t const > local3DGlobalIds,
                               arrayView1d< int const > local3DPartitions,
                               MPI_Comm const comm )
{
  GEOS_MARK_FUNCTION;

  int const numRanks = MpiWrapper::commSize( comm );

  // Build local partition lookup
  stdUnorderedMap< int64_t, int > localPartitionMap;
  localPartitionMap.reserve( local3DGlobalIds.size() );

  for( localIndex i = 0; i < local3DGlobalIds.size(); ++i )
  {
    localPartitionMap.emplace( local3DGlobalIds[i], local3DPartitions[i] );
  }

  // Identify which 3D global IDs we need from other ranks
  SortedArray< int64_t > missingGlobalIds;

  for( localIndex i = 0; i < neighbors2Dto3D.size(); ++i )
  {
    for( int64_t globalId : neighbors2Dto3D[i] )
    {
      if( localPartitionMap.count( globalId ) == 0 )
      {
        missingGlobalIds.insert( globalId );
      }
    }
  }

  // Gather all requested IDs across ranks
  stdVector< int64_t > missingGlobalIdsVec( missingGlobalIds.begin(), missingGlobalIds.end() );
  stdVector< int64_t > allRequestedIdsVec = collectUniqueValues( missingGlobalIdsVec );

  array1d< int64_t > allRequestedIds( allRequestedIdsVec.size() );
  std::copy( allRequestedIdsVec.begin(), allRequestedIdsVec.end(), allRequestedIds.begin() );

  // Each rank contributes partition info for IDs it owns that others need
  array1d< int64_t > contributedGlobalIds;
  array1d< int > contributedPartitions;

  contributedGlobalIds.reserve( allRequestedIds.size() / numRanks );
  contributedPartitions.reserve( allRequestedIds.size() / numRanks );

  for( int64_t requestedId : allRequestedIds )
  {
    auto it = localPartitionMap.find( requestedId );
    if( it != localPartitionMap.end() )
    {
      contributedGlobalIds.emplace_back( requestedId );
      contributedPartitions.emplace_back( it->second );
    }
  }

  // All-gather using the convenience wrapper (handles displacements internally)
  array1d< int64_t > allGlobalIds;
  array1d< int > allPartitions;

  MpiWrapper::allGatherv( contributedGlobalIds.toViewConst(), allGlobalIds, comm );
  MpiWrapper::allGatherv( contributedPartitions.toViewConst(), allPartitions, comm );

  // Build complete partition map from gathered data
  stdUnorderedMap< int64_t, int > completePartitionMap( localPartitionMap );
  completePartitionMap.reserve( localPartitionMap.size() + allGlobalIds.size() );

  for( localIndex i = 0; i < allGlobalIds.size(); ++i )
  {
    completePartitionMap.emplace( allGlobalIds[i], allPartitions[i] );
  }

  // Assign 2D cell partitions using complete map
  array1d< int > partitions2D( neighbors2Dto3D.size() );

  for( localIndex i = 0; i < neighbors2Dto3D.size(); ++i )
  {
    auto neighbors = neighbors2Dto3D[i];

    // Deterministic tie-breaking: minimum global ID
    int64_t const minGlobalId = *std::min_element( neighbors.begin(), neighbors.end() );

    // Look up partition
    auto it = completePartitionMap.find( minGlobalId );
    GEOS_ERROR_IF( it == completePartitionMap.end(),
                   GEOS_FMT( "Partition for 3D neighbor with global ID {} not found", minGlobalId ) );

    partitions2D[i] = it->second;
  }

  return partitions2D;
}

/**
 * @brief Extract cells by indices
 *
 * @param[in] mesh Source mesh
 * @param[in] indices Cell indices to extract
 * @return Extracted cells as unstructured grid (shallow copy)
 */
static vtkSmartPointer< vtkUnstructuredGrid >
extractCellsByIndices( vtkDataSet & mesh,
                       arrayView1d< vtkIdType const > indices )
{
  GEOS_MARK_FUNCTION;

  if( indices.empty() )
  {
    return vtkSmartPointer< vtkUnstructuredGrid >::New();
  }

  vtkNew< vtkIdList > idList;
  idList->SetNumberOfIds( indices.size() );
  for( localIndex i = 0; i < indices.size(); ++i )
  {
    idList->SetId( i, indices[i] );
  }

  vtkNew< vtkExtractCells > extractor;
  extractor->SetInputDataObject( &mesh );
  extractor->SetCellList( idList );
  extractor->Update();

  vtkSmartPointer< vtkUnstructuredGrid > result = vtkSmartPointer< vtkUnstructuredGrid >::New();
  result->ShallowCopy( extractor->GetOutput() );

  return result;
}

/**
 * @brief Extract global IDs for all cells in a mesh
 *
 * @param[in] mesh Mesh containing global IDs
 * @return Array of global IDs for all cells
 */
static array1d< int64_t >
extractGlobalIds( vtkDataSet & mesh )
{
  vtkDataArray * globalIds = mesh.GetCellData()->GetGlobalIds();
  GEOS_ERROR_IF( globalIds == nullptr, "Global IDs not found in mesh" );

  vtkIdType const numCells = mesh.GetNumberOfCells();
  array1d< int64_t > result( numCells );

  for( vtkIdType i = 0; i < numCells; ++i )
  {
    result[i] = static_cast< int64_t >( globalIds->GetTuple1( i ) );
  }

  return result;
}


/**
 * @brief Redistribute 2D cells and fractures, then merge with already-redistributed 3D cells
 *
 * This function completes the 2D/3D/fracture redistribution workflow:
 * 1. Assigns 2D cells to partitions based on their 3D neighbor locations
 * 2. Redistributes 2D cells to appropriate ranks
 * 3. Assigns fracture cells to partitions based on their 3D neighbor locations
 * 4. Redistributes fracture cells to appropriate ranks
 * 5. Merges local 2D and 3D cells on each rank into a unified mesh
 *
 * Both 2D cells and fractures are assigned to ranks using the same strategy:
 * - Find all 3D neighbors of each 2D/fracture element
 * - Assign to the rank owning the 3D neighbor with minimum global ID (deterministic)
 * - Uses efficient MPI communication (only gathers needed 3D partition info)
 *
 * @param[in] redistributed3D Already partitioned 3D cells with global IDs
 * @param[in] originalMesh Original mesh containing 2D cells (on rank 0)
 * @param[in] cells2DIndices Indices of 2D cells in original mesh
 * @param[in] neighbors2Dto3D Pre-computed 2D-to-3D neighbor mapping (fracture element -> 3D cell global IDs)
 * @param[in] unpartitionedFractures Unpartitioned fracture meshes (on rank 0, empty on others)
 * @param[in] fractureNeighbors Pre-computed fracture-to-3D neighbor mappings (on rank 0, empty on others)
 * @param[in] fractureNames Ordered list of fracture names (consistent across all ranks)
 * @param[in] comm MPI communicator
 * @return Complete AllMeshes object with merged main mesh (3D+2D) and redistributed fractures
 */
static AllMeshes
redistribute2DAndMergeWith3D( vtkSmartPointer< vtkDataSet > redistributed3D,
                              vtkSmartPointer< vtkDataSet > originalMesh,
                              arrayView1d< vtkIdType const > cells2DIndices,
                              ArrayOfArrays< localIndex, int64_t > const & neighbors2Dto3D,
                              stdMap< string, vtkSmartPointer< vtkDataSet > > const & unpartitionedFractures,
                              stdMap< string, ArrayOfArrays< localIndex, int64_t > > const & fractureNeighbors,
                              stdVector< string > const & fractureNames,
                              MPI_Comm const comm )
{
  GEOS_MARK_FUNCTION;

  int const rank = MpiWrapper::commRank( comm );
  int const numRanks = MpiWrapper::commSize( comm );

  // -----------------------------------------------------------------------
  // Step 1: Assign and redistribute 2D cells
  // -----------------------------------------------------------------------
  array1d< int64_t > local3DGlobalIds = extractGlobalIds( *redistributed3D );
  array1d< int > partitions3D( redistributed3D->GetNumberOfCells() );
  partitions3D.setValues< parallelHostPolicy >( rank );

  bool const hasLocal2DCells = (rank == 0 && !cells2DIndices.empty());

  array1d< int > partitions2D = assignCellsBasedOn3DNeighbors(
    hasLocal2DCells ? neighbors2Dto3D : ArrayOfArrays< localIndex, int64_t >{},
    local3DGlobalIds.toViewConst(),
    partitions3D.toViewConst(),
    comm );

  vtkSmartPointer< vtkUnstructuredGrid > cells2D = hasLocal2DCells
      ? extractCellsByIndices( *originalMesh, cells2DIndices )
      : vtkSmartPointer< vtkUnstructuredGrid >::New();

  vtkSmartPointer< vtkPartitionedDataSet > split2D =
    splitMeshByPartition( cells2D, numRanks, partitions2D.toViewConst() );

  vtkSmartPointer< vtkUnstructuredGrid > redistributed2D =
    vtk::redistribute( *split2D, comm );

  // Conservation check
  vtkIdType const total2DCells = MpiWrapper::sum( redistributed2D->GetNumberOfCells(), comm );
  vtkIdType expected2DCells = cells2DIndices.size();
  MpiWrapper::broadcast( expected2DCells, 0, comm );

  GEOS_ERROR_IF( total2DCells != expected2DCells,
                 GEOS_FMT( "2D cell redistribution failed: expected {} cells, got {} cells",
                           expected2DCells, total2DCells ) );

  // -----------------------------------------------------------------------
  // Step 2: Redistribute fractures
  // -----------------------------------------------------------------------
  stdMap< string, vtkSmartPointer< vtkDataSet > > redistributedFractures;

  for( string const & fractureName : fractureNames )
  {
    vtkIdType expectedFractureCells = 0;

    ArrayOfArrays< localIndex, int64_t > localFractureNeighbors;
    vtkSmartPointer< vtkDataSet > unpartitionedFracture;

    if( rank == 0 )
    {
      auto fracIt = unpartitionedFractures.find( fractureName );
      unpartitionedFracture = (fracIt != unpartitionedFractures.end()) ? fracIt->second : nullptr;

      if( unpartitionedFracture && unpartitionedFracture->GetNumberOfCells() > 0 )
      {
        expectedFractureCells = unpartitionedFracture->GetNumberOfCells();

        auto fracNeighborIt = fractureNeighbors.find( fractureName );
        GEOS_ERROR_IF( fracNeighborIt == fractureNeighbors.end(),
                       GEOS_FMT( "Fracture '{}' not found in fractureNeighbors map", fractureName ) );

        localFractureNeighbors = fracNeighborIt->second;
      }
    }

    array1d< int > partitionsFracture = assignCellsBasedOn3DNeighbors(
      localFractureNeighbors,  // Empty on non-root ranks
      local3DGlobalIds.toViewConst(),
      partitions3D.toViewConst(),
      comm );

    // Split and redistribute
    vtkSmartPointer< vtkPartitionedDataSet > splitFracture;

    if( rank == 0 && unpartitionedFracture && unpartitionedFracture->GetNumberOfCells() > 0 )
    {
      splitFracture = splitMeshByPartition( unpartitionedFracture, numRanks,
                                            partitionsFracture.toViewConst() );
    }
    else
    {
      splitFracture = vtkSmartPointer< vtkPartitionedDataSet >::New();
      splitFracture->SetNumberOfPartitions( numRanks );
      for( int r = 0; r < numRanks; ++r )
      {
        splitFracture->SetPartition( r, vtkSmartPointer< vtkUnstructuredGrid >::New() );
      }
    }

    vtkSmartPointer< vtkUnstructuredGrid > localFracture =
      vtk::redistribute( *splitFracture, comm );

    vtkIdType const totalFractureCells = MpiWrapper::sum( localFracture->GetNumberOfCells(), comm );
    MpiWrapper::broadcast( expectedFractureCells, 0, comm );

    GEOS_ERROR_IF( totalFractureCells != expectedFractureCells,
                   GEOS_FMT( "Fracture '{}' redistribution lost cells: expected {}, got {}",
                             fractureName, expectedFractureCells, totalFractureCells ) );

    redistributedFractures.insert( { fractureName, localFracture } );
  }

  // -----------------------------------------------------------------------
  // Step 3: Merge local 2D and 3D cells
  // -----------------------------------------------------------------------
  vtkSmartPointer< vtkUnstructuredGrid > mergedMesh = vtkSmartPointer< vtkUnstructuredGrid >::New();

  if( redistributed2D->GetNumberOfCells() > 0 )
  {
    vtkNew< vtkAppendFilter > appendFilter;
    appendFilter->AddInputData( redistributed3D );
    appendFilter->AddInputData( redistributed2D );
    appendFilter->MergePointsOn();
    appendFilter->Update();

    mergedMesh->ShallowCopy( appendFilter->GetOutput() );
  }
  else
  {
    mergedMesh->ShallowCopy( redistributed3D );
  }

  return AllMeshes( mergedMesh, redistributedFractures );
}

/**
 * @brief Find 3D cells whose faces exactly match a fracture element
 *
 * A 3D cell matches if it has a face that shares all nodes with the fracture element,
 * accounting for collocated nodes at split interfaces.
 *
 * @param fractureNodeIds Local node IDs of the fracture element
 * @param collocatedNodes Mapping from local node ID to all collocated global IDs
 * @param nodesToCells Reverse map from global node ID to cells containing that node
 * @return Global IDs of matching 3D cells (typically 0-2 neighbors for fractures)
 */
stdVector< vtkIdType > findMatchingCellsForFractureElement(
  vtkIdList * fractureNodeIds,
  CollocatedNodes const & collocatedNodes,
  stdMap< vtkIdType, std::set< vtkIdType > > const & nodesToCells )
{
  vtkIdType const numFractureNodes = fractureNodeIds->GetNumberOfIds();

  // Build set of ALL collocated nodes for this fracture element
  std::unordered_set< vtkIdType > fractureCollocatedNodes;
  fractureCollocatedNodes.reserve( numFractureNodes * 2 );  // Typical case: 2 versions per node

  for( vtkIdType j = 0; j < numFractureNodes; ++j )
  {
    vtkIdType const localNodeIdx = fractureNodeIds->GetId( j );
    stdVector< vtkIdType > const & ns = collocatedNodes[ localNodeIdx ];
    fractureCollocatedNodes.insert( ns.begin(), ns.end() );
  }

  // Build map: candidate cellId -> set of its nodes that match fracture's collocated nodes
  stdUnorderedMap< vtkIdType, std::unordered_set< vtkIdType > > cellToMatchedNodes;
  cellToMatchedNodes.reserve( fractureCollocatedNodes.size() );

  for( vtkIdType const & collocatedNode : fractureCollocatedNodes )
  {
    auto it = nodesToCells.find( collocatedNode );
    if( it != nodesToCells.cend() )
    {
      for( vtkIdType const & cellId : it->second )
      {
        cellToMatchedNodes.get_inserted( cellId ).insert( collocatedNode );
      }
    }
  }

  // Filter to cells that form a valid matching face
  stdVector< vtkIdType > matchingCells;
  matchingCells.reserve( 2 );  // Most fractures have 0-2 neighbors

  for( auto const & cellToNodes : cellToMatchedNodes )
  {
    vtkIdType const cellId = cellToNodes.first;
    auto const & matchedNodes = cellToNodes.second;

    // Must match exactly the number of fracture nodes
    if( matchedNodes.size() != static_cast< std::size_t >( numFractureNodes ) )
    {
      continue;
    }

    // Verify each fracture node has at least one collocated version in the matched set
    // (ensures we matched a true face, not just any numFractureNodes nodes)
    bool allFractureNodesRepresented = true;

    for( vtkIdType j = 0; j < numFractureNodes; ++j )
    {
      vtkIdType const localNodeIdx = fractureNodeIds->GetId( j );
      stdVector< vtkIdType > const & nodeCollocated = collocatedNodes[ localNodeIdx ];

      // Check if ANY collocated version of this node is in the matched set
      bool nodeRepresented = std::any_of(
        nodeCollocated.begin(),
        nodeCollocated.end(),
        [&matchedNodes]( vtkIdType collocNode ) { return matchedNodes.count( collocNode ) > 0; }
        );

      if( !nodeRepresented )
      {
        allFractureNodesRepresented = false;
        break;
      }
    }

    if( allFractureNodesRepresented )
    {
      matchingCells.push_back( cellId );

      // Early exit - most fractures have ≤2 neighbors (boundary or internal)
      if( matchingCells.size() >= 2 )
      {
        return matchingCells;
      }
    }
  }

  return matchingCells;
}

/**
 * @brief Build fracture-to-3D neighbor connectivity
 *
 * @param originalMesh The original mesh containing both 3D cells and fractures
 * @param fractureMesh The fracture mesh (separate mesh on rank 0)
 * @param cells3DIndices Indices of 3D cells in the original mesh
 * @return Mapping: fracture element index -> global IDs of neighboring 3D cells
 */
static ArrayOfArrays< localIndex, int64_t >
buildFractureTo3DNeighbors( vtkDataSet & originalMesh,
                            vtkSmartPointer< vtkDataSet > fractureMesh,
                            arrayView1d< vtkIdType const > cells3DIndices )
{
  GEOS_MARK_FUNCTION;

  vtkIdType const numFractureElems = fractureMesh->GetNumberOfCells();

  vtkDataArray * meshGlobalCellIds = originalMesh.GetCellData()->GetGlobalIds();
  vtkDataArray * meshGlobalNodeIds = originalMesh.GetPointData()->GetGlobalIds();

  GEOS_ERROR_IF( meshGlobalCellIds == nullptr, "Original mesh must have cell GlobalIds" );
  GEOS_ERROR_IF( meshGlobalNodeIds == nullptr, "Original mesh must have node GlobalIds" );

  // -----------------------------------------------------------------------
  // Step 1: Build reverse lookup: original mesh index -> global cell ID (3D cells only)
  // -----------------------------------------------------------------------
  stdUnorderedMap< vtkIdType, int64_t > mesh3DIdxToGlobalId;
  mesh3DIdxToGlobalId.reserve( cells3DIndices.size() );

  for( localIndex i = 0; i < cells3DIndices.size(); ++i )
  {
    vtkIdType const origIdx = cells3DIndices[i];
    int64_t const globalId = static_cast< int64_t >( meshGlobalCellIds->GetTuple1( origIdx ) );
    mesh3DIdxToGlobalId.emplace( origIdx, globalId );
  }

  // -----------------------------------------------------------------------
  // Step 2: Build node-to-3D-cell connectivity for the original mesh
  // -----------------------------------------------------------------------
  stdMap< vtkIdType, std::set< vtkIdType > > nodeGlobalIdToCells3D;

  vtkSmartPointer< vtkIdList > const ptIds = vtkSmartPointer< vtkIdList >::New();
  for( localIndex i = 0; i < cells3DIndices.size(); ++i )
  {
    vtkIdType const origCellIdx = cells3DIndices[i];
    originalMesh.GetCellPoints( origCellIdx, ptIds );

    for( vtkIdType p = 0; p < ptIds->GetNumberOfIds(); ++p )
    {
      vtkIdType const nodeLocalId = ptIds->GetId( p );
      vtkIdType const nodeGlobalId = static_cast< vtkIdType >( meshGlobalNodeIds->GetTuple1( nodeLocalId ) );

      nodeGlobalIdToCells3D.get_inserted( nodeGlobalId ).insert( origCellIdx );
    }
  }

  // -----------------------------------------------------------------------
  // Step 3: Build collocated nodes for fracture mesh
  // -----------------------------------------------------------------------
  CollocatedNodes collocatedNodes( "fracture", fractureMesh, false );

  // -----------------------------------------------------------------------
  // Step 4: Build fracture-to-3D neighbor mapping
  // -----------------------------------------------------------------------
  ArrayOfArrays< localIndex, int64_t > result;
  result.reserve( numFractureElems );

  vtkSmartPointer< vtkIdList > const fracPtIds = vtkSmartPointer< vtkIdList >::New();
  for( vtkIdType fracElemIdx = 0; fracElemIdx < numFractureElems; ++fracElemIdx )
  {
    fractureMesh->GetCellPoints( fracElemIdx, fracPtIds );

    // Find matching 3D cells using exact node matching
    stdVector< vtkIdType > matchingOriginalCells = findMatchingCellsForFractureElement(
      fracPtIds,
      collocatedNodes,
      nodeGlobalIdToCells3D );

    // Convert original mesh indices to global cell IDs
    array1d< int64_t > neighbor3DGlobalIds;
    neighbor3DGlobalIds.reserve( matchingOriginalCells.size() );

    for( vtkIdType origCellIdx : matchingOriginalCells )
    {
      auto it = mesh3DIdxToGlobalId.find( origCellIdx );
      if( it != mesh3DIdxToGlobalId.end() )
      {
        neighbor3DGlobalIds.emplace_back( it->second );
      }
    }

    // Error on orphaned fractures
    GEOS_ERROR_IF( neighbor3DGlobalIds.empty(),
                   GEOS_FMT( "Fracture element {} has no 3D neighbors", fracElemIdx ) );

    result.appendArray( neighbor3DGlobalIds.begin(), neighbor3DGlobalIds.end() );
  }

  return result;
}


vtkSmartPointer< vtkUnstructuredGrid >
threshold( vtkDataSet & mesh,
           string const & arrayName,
           int const component,
           double const lo,
           double const hi )
{
  vtkNew< vtkThreshold > threshold;
  threshold->SetInputDataObject( &mesh );
  threshold->SetInputArrayToProcess( 0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_CELLS, arrayName.c_str() );
  threshold->SetComponentModeToUseSelected();
  threshold->SetSelectedComponent( component );
  threshold->SetThresholdFunction( vtkThreshold::THRESHOLD_BETWEEN );
  threshold->SetLowerThreshold( lo );
  threshold->SetUpperThreshold( hi );
  threshold->Update();
  return threshold->GetOutput();
}

template< typename VTK_TYPES, typename FUNC >
void dispatchArray( vtkDataSetAttributes & data,
                    string const & arrayName,
                    FUNC && func )
{
  vtkDataArray * const array = data.GetArray( arrayName.c_str() );
  GEOS_THROW_IF( array == nullptr, GEOS_FMT( "VTK array '{}' not found", arrayName ), InputError );
  bool const result = vtkArrayDispatch::DispatchByValueType< VTK_TYPES >::Execute( array, [&]( auto const * const typedArray )
  {
    vtkDataArrayAccessor< TYPEOFPTR( typedArray ) > const accessor( typedArray );
    func( accessor );
  } );
  GEOS_THROW_IF( !result,
                 GEOS_FMT( "VTK dispatch failed, array '{}' is not of expected type", arrayName ),
                 InputError );
}

std::array< std::pair< int, int >, 2 >
findGlobalIndexBounds( vtkDataSet & mesh,
                       MPI_Comm const & comm,
                       string const & indexArrayName )
{
  RAJA::ReduceMin< parallelHostReduce, int > minIdx0{ std::numeric_limits< int >::max() }, minIdx1{ std::numeric_limits< int >::max() };
  RAJA::ReduceMax< parallelHostReduce, int > maxIdx0{ std::numeric_limits< int >::min() }, maxIdx1{ std::numeric_limits< int >::min() };
  dispatchArray< vtkArrayDispatch::Integrals >( *mesh.GetCellData(), indexArrayName, [&]( auto const index )
  {
    forAll< parallelHostPolicy >( mesh.GetNumberOfCells(), [=]( vtkIdType const i )
    {
      auto const idx0 = index.Get( i, 0 );
      minIdx0.min( idx0 );
      maxIdx0.max( idx0 );
      auto const idx1 = index.Get( i, 1 );
      minIdx1.min( idx1 );
      maxIdx1.max( idx1 );
    } );
  } );

  int const minIdxLocal[2] = { minIdx0.get(), minIdx1.get() };
  int const maxIdxLocal[2] = { maxIdx0.get(), maxIdx1.get() };
  int minIdxGlobal[2], maxIdxGlobal[2];
  MpiWrapper::min< int >( minIdxLocal, minIdxGlobal, comm );
  MpiWrapper::max< int >( maxIdxLocal, maxIdxGlobal, comm );
  return { std::make_pair( minIdxGlobal[0], maxIdxGlobal[0] ), std::make_pair( minIdxGlobal[1], maxIdxGlobal[1] ) };
}

AllMeshes
redistributeByAreaGraphAndLayer( AllMeshes & input,
                                 PartitionMethod const method,
                                 string const & indexArrayName,
                                 MPI_Comm const comm,
                                 int const numPartZ,
                                 int const numRefinements )
{
  GEOS_MARK_FUNCTION;

  localIndex const numCells = LvArray::integerConversion< localIndex >( input.getMainMesh()->GetNumberOfCells() );
  int const numProcs = MpiWrapper::commSize( comm );
  int const numPartA = numProcs / numPartZ;
  int const numProcsRemainder = numProcs % numPartZ;
  GEOS_ERROR_IF_NE_MSG( numProcsRemainder, 0, "Number of ranks must evenly divide the number of z-partitions" );

  // Compute conversion from cell z-index to partition z-index
  std::array< std::pair< int, int >, 2 > const idxLimits = findGlobalIndexBounds( *input.getMainMesh(), comm, indexArrayName );
  double const cellPerPartZInv = static_cast< double >( numPartZ ) / ( idxLimits[1].second - idxLimits[1].first + 1 );
  auto const computePartIndexZ = [minZ = idxLimits[1].first, cellPerPartZInv]( auto const zidx )
  {
    return static_cast< int >( ( zidx - minZ + 0.5 ) * cellPerPartZInv );
  };

  // Extract cells in z-layer "0"
  vtkSmartPointer< vtkUnstructuredGrid > layer0 = threshold( *input.getMainMesh(), indexArrayName, 1, idxLimits[1].first, idxLimits[1].first );

  // Ranks that have the layer form a subcomm and compute the area partitioning
  bool const haveLayer0 = layer0->GetNumberOfCells() > 0;
  MPI_Comm subComm = MpiWrapper::commSplit( comm, haveLayer0 ? 0 : MPI_UNDEFINED, MpiWrapper::commRank( comm ) );
  array1d< int64_t > layer0Parts;
  if( haveLayer0 )
  {
    AllMeshes layer0input( layer0, {} ); // fracture mesh not supported yet
    layer0Parts = partitionByCellGraph( layer0input.getMainMesh(), method, subComm, numPartA, 3, numRefinements );
    MpiWrapper::commFree( subComm );
  }

  // pack the area index into unused high 32 bits of the 64-bit partition index
  dispatchArray< vtkArrayDispatch::Integrals >( *layer0->GetCellData(), indexArrayName,
                                                [&]( auto const index )
  {
    forAll< parallelHostPolicy >( layer0Parts.size(), [dst = layer0Parts.toView(), index,
                                                       minA = idxLimits[0].first]( localIndex const i )
    {
      auto const aidx = index.Get( i, 0 );

      GEOS_ASSERT_GT( std::numeric_limits< int32_t >::max(), dst[i] );
      if constexpr ( std::is_signed_v< decltype(aidx) > )
        GEOS_ASSERT_GE( aidx, 0 );

      dst[i] |= static_cast< int64_t >( aidx - minA ) << 32;
    } );
  } );

  // Distribute the partitioning of layer 0 (or columns) to everyone
  array1d< int64_t > columnParts;
  MpiWrapper::allGatherv( layer0Parts.toViewConst(), columnParts, comm );

  // Sort w.r.t. area index. If this becomes too expensive, may need a diff approach.
  // The goal is to enable fast lookup of area partition index by area index.
  std::sort( columnParts.begin(), columnParts.end() );

  auto const computePartIndexA = [minA = idxLimits[0].first, partsA = columnParts.toViewConst()]( auto const aidx )
  {
    return partsA[aidx - minA] & 0xFFFFFFFF;
  };

  // Finally, we can compute the target partitioning by combining area and z-partition data
  array1d< int > const newParts( numCells );
  dispatchArray< vtkArrayDispatch::Integrals >( *input.getMainMesh()->GetCellData(), indexArrayName,
                                                [&]( auto const index )
  {
    forAll< parallelHostPolicy >( numCells, [computePartIndexA, computePartIndexZ,
                                             newParts = newParts.toView(),
                                             index, numPartZ]( localIndex const i )
    {
      int const partIdxA = computePartIndexA( index.Get( i, 0 ) );
      int const partIdxZ = computePartIndexZ( index.Get( i, 1 ) );
      newParts[i] = partIdxA * numPartZ + partIdxZ;
    } );
  } );
  vtkSmartPointer< vtkPartitionedDataSet > const splitMesh = splitMeshByPartition( input.getMainMesh(), numProcs, newParts.toViewConst() );
  return AllMeshes( vtk::redistribute( *splitMesh, MPI_COMM_GEOS ), {} );
}

/**
 * @brief Redistributes the mesh using a Kd-Tree
 *
 * @param[in] mesh a vtk grid
 * @return the vtk grid redistributed
 */
vtkSmartPointer< vtkDataSet >
redistributeByKdTree( vtkDataSet & mesh )
{
  GEOS_MARK_FUNCTION;

  // Count input cells for verification
  vtkIdType localInputCells = mesh.GetNumberOfCells();
  vtkIdType globalInputCells = MpiWrapper::allReduce( localInputCells, MpiWrapper::Reduction::Sum, MPI_COMM_GEOS );

  // Use a VTK filter which employs a kd-tree partition internally
  vtkNew< vtkRedistributeDataSetFilter > rdsf;
  rdsf->SetInputDataObject( &mesh );
  rdsf->SetNumberOfPartitions( MpiWrapper::commSize() );
  rdsf->Update();

  vtkSmartPointer< vtkDataSet > result = vtkDataSet::SafeDownCast( rdsf->GetOutputDataObject( 0 ) );

  // Verify we didn't lose any cells
  vtkIdType localOutputCells = result->GetNumberOfCells();
  vtkIdType globalOutputCells = MpiWrapper::allReduce( localOutputCells, MpiWrapper::Reduction::Sum, MPI_COMM_GEOS );

  if( globalOutputCells != globalInputCells )
  {
    if( MpiWrapper::commRank() == 0 )
    {
      GEOS_WARNING( GEOS_FMT( "VTK KdTree redistribution lost {} elements! Falling back to block redistribution.",
                              globalInputCells - globalOutputCells ) );
    }
    return scatterByBlock( mesh );
  }

  return result;
}

stdVector< int >
findNeighborRanks( stdVector< vtkBoundingBox > boundingBoxes )
{
  int const numParts = LvArray::integerConversion< int >( boundingBoxes.size() );
  int const thisRank = MpiWrapper::commRank();

  // Inflate boxes to detect intersections more reliably.
  //
  // Pure relative scaling (ScaleAboutCenter) fails for small or thin partitions in Release builds:
  // two face-adjacent partitions share an exact planar boundary, so their bounding boxes only
  // *touch* rather than overlap. After scaling by 1% about the center, the expanded extents of
  // one box may still not penetrate the other in that dimension if the box is very thin relative
  // to floating-point rounding under -O2/-O3.
  //
  // Fix: compute a per-box absolute inflation in each dimension as
  //   max( relative_expand, absoluteMinExpand )
  // where absoluteMinExpand is a small fixed distance that guarantees overlap even for
  // degenerate (zero-length) extents along a split plane.
  double constexpr relativeExpandFactor = 0.01; // 1% of the box length in each dimension
  double constexpr absoluteMinExpand    = 1.0e-6; // absolute fallback (mesh units)

  for( vtkBoundingBox & box : boundingBoxes )
  {
    double const * minPt = box.GetMinPoint();
    double const * maxPt = box.GetMaxPoint();
    double newMin[3], newMax[3];
    for( int d = 0; d < 3; ++d )
    {
      double const half = 0.5 * ( maxPt[d] - minPt[d] );
      double const expand = std::max( relativeExpandFactor * half, absoluteMinExpand );
      newMin[d] = minPt[d] - expand;
      newMax[d] = maxPt[d] + expand;
    }
    box.SetMinPoint( newMin );
    box.SetMaxPoint( newMax );
  }

  stdVector< int > neighbors;
  for( int i = 0; i < numParts; ++i )
  {
    if( i != thisRank && boundingBoxes[thisRank].Intersects( boundingBoxes[i] ) )
    {
      neighbors.push_back( i );
    }
  }

  return neighbors;
}


vtkSmartPointer< vtkDataSet > manageGlobalIds( vtkSmartPointer< vtkDataSet > mesh, int useGlobalIds, bool isFractured )
{
  auto hasGlobalIds = []( vtkSmartPointer< vtkDataSet > m ) -> bool
  {
    return m->GetPointData()->GetGlobalIds() != nullptr && m->GetCellData()->GetGlobalIds() != nullptr;
  };

  {
    // Add global ids on the fly if needed
    int const me = hasGlobalIds( mesh );
    int const everyone = MpiWrapper::allReduce( me, MpiWrapper::Reduction::Max, MPI_COMM_GEOS );

    if( everyone and not me )
    {
      mesh->GetPointData()->SetGlobalIds( vtkIdTypeArray::New() );
      mesh->GetCellData()->SetGlobalIds( vtkIdTypeArray::New() );
    }
  }

  vtkSmartPointer< vtkDataSet > output;
  bool const globalIdsAvailable = hasGlobalIds( mesh );
  if( useGlobalIds > 0 && !globalIdsAvailable )
  {
    GEOS_ERROR( "Global IDs strictly required (useGlobalId > 0) but unavailable. Set useGlobalIds to 0 to build them automatically." );
  }
  else if( useGlobalIds >= 0 && globalIdsAvailable )
  {
    output = mesh;
    vtkIdTypeArray const * const globalCellId = vtkIdTypeArray::FastDownCast( output->GetCellData()->GetGlobalIds() );
    vtkIdTypeArray const * const globalPointId = vtkIdTypeArray::FastDownCast( output->GetPointData()->GetGlobalIds() );
    GEOS_ERROR_IF( globalCellId->GetNumberOfComponents() != 1 && globalCellId->GetNumberOfTuples() != output->GetNumberOfCells(),
                   GEOS_FMT( "Global cell IDs are invalid. Check the array or enable automatic generation (useGlobalId < 0).\n{}",
                             generalMeshErrorAdvice ) );
    GEOS_ERROR_IF( globalPointId->GetNumberOfComponents() != 1 && globalPointId->GetNumberOfTuples() != output->GetNumberOfPoints(),
                   GEOS_FMT( "Global cell IDs are invalid. Check the array or enable automatic generation (useGlobalId < 0).\n{}",
                             generalMeshErrorAdvice ) );

    GEOS_LOG_RANK_0( "Using global Ids defined in VTK mesh" );
  }
  else
  {
    GEOS_ERROR_IF( isFractured,
                   GEOS_FMT( "Automatic generation of global IDs for fractured meshes is disabled. Please split with  mesh_doctor. \n{}",
                             generalMeshErrorAdvice ) );

    GEOS_LOG_RANK_0( "Generating global Ids from VTK mesh" );
    output = generateGlobalIDs( mesh );
  }

  return output;
}

/**
 * @brief This function tries to make sure that no MPI rank is empty
 *
 * @param[in] mesh a vtk grid
 * @param[in] comm the MPI communicator
 * @return the vtk grid redistributed
 */
vtkSmartPointer< vtkDataSet >
ensureNoEmptyRank( vtkSmartPointer< vtkDataSet > mesh,
                   MPI_Comm const comm )
{
  GEOS_MARK_FUNCTION;

  // step 1: figure out who is a donor and who is a recipient
  localIndex const numElems = LvArray::integerConversion< localIndex >( mesh->GetNumberOfCells() );
  integer const numProcs = MpiWrapper::commSize( comm );

  array1d< localIndex > elemCounts( numProcs );
  MpiWrapper::allGather( numElems, elemCounts, comm );

  SortedArray< integer > recipientRanks;
  array1d< integer > donorRanks;
  recipientRanks.reserve( numProcs );
  donorRanks.reserve( numProcs );

  for( integer iRank = 0; iRank < numProcs; ++iRank )
  {
    if( elemCounts[iRank] == 0 )
    {
      recipientRanks.insert( iRank );
    }
    else if( elemCounts[iRank] > 1 ) // need at least two elems to be a donor
    {
      donorRanks.emplace_back( iRank );
    }
  }

  // step 2: at this point, we need to determine the donors and which cells they donate

  // First we sort the donor in order of the number of elems they contain
  std::stable_sort( donorRanks.begin(), donorRanks.end(),
                    [&elemCounts] ( auto i1, auto i2 )
  { return elemCounts[i1] > elemCounts[i2]; } );

  // Then, if my position is "i" in the donorRanks array, I will send half of my elems to the i-th recipient
  integer const myRank = MpiWrapper::commRank();
  auto const pos = std::find( donorRanks.begin(), donorRanks.end(), myRank );
  bool const isDonor = ( pos != donorRanks.end() );

  // step 3: my rank was selected to donate cells, let's proceed
  // we need to make a distinction between two configurations

  array1d< localIndex > newParts( numElems );
  newParts.setValues< parallelHostPolicy >( myRank );

  // step 3.1: donorRanks.size() >= recipientRanks.size()
  // we use a strategy that preserves load balancing
  if( isDonor && donorRanks.size() >= recipientRanks.size() )
  {
    auto const myPosition = std::distance( donorRanks.begin(), pos );
    if( myPosition < recipientRanks.size() )
    {
      integer const recipientRank = recipientRanks[myPosition];
      for( localIndex iElem = numElems/2; iElem < numElems; ++iElem )
      {
        newParts[iElem] = recipientRank; // I donate half of my cells
      }
    }
  }
  // step 3.2: donorRanks.size() < recipientRanks.size()
  // this is the unhappy path: we don't care anymore about load balancing at this stage
  // we just want the simulation to run and count on ParMetis/PTScotch to restore load balancing
  else if( isDonor && donorRanks.size() < recipientRanks.size() )
  {
    auto const myPosition = std::distance( donorRanks.begin(), pos );
    localIndex firstRecipientPosition = 0;
    for( integer iRank = 0; iRank < myPosition; ++iRank )
    {
      firstRecipientPosition += elemCounts[donorRanks[iRank]] - 1;
    }
    if( firstRecipientPosition < recipientRanks.size() )
    {
      bool const isLastDonor = myPosition == donorRanks.size() - 1;
      localIndex const lastRecipientPosition = firstRecipientPosition + numElems - 1;
      GEOS_THROW_IF( isLastDonor && ( lastRecipientPosition < recipientRanks.size() ),
                     "The current implementation is unable to guarantee that all ranks have at least one element",
                     geos::RuntimeError );

      for( localIndex iElem = 1; iElem < numElems; ++iElem ) // I only keep my first element
      {
        // this is the brute force approach
        // each donor donates all its elems except the first one
        localIndex const recipientPosition = firstRecipientPosition + iElem - 1;
        if( recipientPosition < recipientRanks.size() )
        {
          newParts[iElem] = recipientRanks[recipientPosition];
        }
      }
    }
  }

  GEOS_WARNING_IF( donorRanks.size() < recipientRanks.size(),
                   "We strongly encourage the use of partitionRefinement > 5 for this number of MPI ranks" );

  vtkSmartPointer< vtkPartitionedDataSet > const splitMesh = splitMeshByPartition( mesh, numProcs, newParts.toViewConst() );
  return vtk::redistribute( *splitMesh, MPI_COMM_GEOS );
}

AllMeshes
redistributeMeshes( integer const logLevel,
                    vtkSmartPointer< vtkDataSet > loadedMesh,
                    stdMap< string, vtkSmartPointer< vtkDataSet > > & namesToFractures,
                    MPI_Comm const comm,
                    PartitionMethod const method,
                    int const partitionRefinement,
                    int const partitionFractureWeight,
                    int const useGlobalIds,
                    string const & structuredIndexAttributeName,
                    int const numPartZ )
{
  GEOS_MARK_FUNCTION;
  int const numRanks = MpiWrapper::commSize( comm );
  int const rank = MpiWrapper::commRank( comm );

  stdVector< vtkSmartPointer< vtkDataSet > > fractures;
  for( auto & nameToFracture: namesToFractures )
  {
    fractures.push_back( nameToFracture.second );
  }

  // Generate global IDs for vertices and cells, if needed
  vtkSmartPointer< vtkDataSet > mesh = manageGlobalIds( loadedMesh, useGlobalIds, !std::empty( fractures ) );

  // Ensure mesh is always a valid VTK object, even if empty
  if( !mesh || (mesh->GetNumberOfCells() == 0 && mesh->GetNumberOfPoints() == 0) )
  {
    mesh = vtkSmartPointer< vtkUnstructuredGrid >::New();
  }

  // Verify fractures are on rank 0
  if( rank != 0 )
  {
    for( auto const & [name, fracture]: namesToFractures )
    {
      GEOS_UNUSED_VAR( name );
      GEOS_ASSERT_EQ( fracture->GetNumberOfCells(), 0 );
    }
  }

  // -----------------------------------------------------------------------
  // Step 1: Classify cells by dimension
  // -----------------------------------------------------------------------
  array1d< vtkIdType > cells3DIndices, cells2DIndices;
  classifyCellsByDimension( *mesh, cells3DIndices, cells2DIndices );

  // -----------------------------------------------------------------------
  // Step 2: Build 2D-to-3D neighbor mapping
  // -----------------------------------------------------------------------
  ArrayOfArrays< localIndex, int64_t > neighbors2Dto3D;
  if( !cells2DIndices.empty() )
  {
    neighbors2Dto3D = build2DTo3DNeighbors( *mesh,
                                            cells2DIndices.toViewConst(),
                                            cells3DIndices.toViewConst() );
  }

  // -----------------------------------------------------------------------
  // Step 3: Build fracture-to-3D neighbor mappings
  // -----------------------------------------------------------------------
  stdMap< string, ArrayOfArrays< localIndex, int64_t > > fractureNeighbors;
  stdVector< string > fractureNames;

  // Collect fracture names on rank 0
  if( rank == 0 )
  {
    for( auto const & [fractureName, fractureMesh]: namesToFractures )
    {
      fractureNames.push_back( fractureName );
    }
  }

  // Broadcast fracture names to all ranks
  {
    int numFractures = LvArray::integerConversion< int >( fractureNames.size() );
    MpiWrapper::broadcast( numFractures, 0, comm );

    if( rank != 0 )
    {
      fractureNames.resize( numFractures );
    }

    for( string & name : fractureNames )
    {
      MpiWrapper::broadcast( name, 0, comm );
    }
  }

  // Initialize empty neighbor arrays
  for( auto const & fractureName : fractureNames )
  {
    fractureNeighbors.get_inserted( fractureName ).resize( 0, 0 );
  }

  // Build actual neighbors only on rank 0
  if( rank == 0 )
  {
    for( auto const & [fractureName, fractureMesh]: namesToFractures )
    {
      if( fractureMesh && fractureMesh->GetNumberOfCells() > 0 )
      {
        fractureNeighbors[fractureName] = buildFractureTo3DNeighbors(
          *mesh,
          fractureMesh,
          cells3DIndices.toViewConst() );
      }
    }
  }

  // -----------------------------------------------------------------------
  // Step 4: Extract cells and tag super-cells
  // -----------------------------------------------------------------------
  SuperCellInfo superCellInfo;
  bool hasSuperCells = false;

  vtkSmartPointer< vtkUnstructuredGrid > cells3D = extractCellsByIndices( *mesh, cells3DIndices );

  if( rank == 0 && !fractureNames.empty() )
  {
    superCellInfo = tagCellsWithSuperCellIds( cells3D, fractureNeighbors, partitionFractureWeight );

    // Verify SuperCellId array was actually created
    vtkIdTypeArray * scArray =
      vtkIdTypeArray::SafeDownCast( cells3D->GetCellData()->GetArray( "SuperCellId" ) );
    hasSuperCells = (scArray != nullptr);

    if( hasSuperCells )
    {
      GEOS_LOG_RANK_0( GEOS_FMT( "Tagged {} super-cells from fracture connectivity",
                                 superCellInfo.atomicSuperCells.size() ));
    }
  }

  // Broadcast whether we have super-cells
  {
    int hasSuperCellsInt = hasSuperCells ? 1 : 0;
    MpiWrapper::broadcast( hasSuperCellsInt, 0, comm );
    hasSuperCells = (hasSuperCellsInt > 0);
  }

  // -----------------------------------------------------------------------
  // Step 5: Initial redistribution of 3D cells
  // -----------------------------------------------------------------------
  vtkIdType const minCellsOnAnyRank = MpiWrapper::min( cells3D->GetNumberOfCells(), comm );

  vtkSmartPointer< vtkDataSet > redistributed3D;

  if( minCellsOnAnyRank == 0 )
  {
    if( hasSuperCells )
    {
      redistributed3D = redistributeBySuperCellBlocks( cells3D, comm );
    }
    else
    {
      GEOS_LOG_RANK_0( "Initial redistribution using KD-tree..." );
      redistributed3D = redistributeByKdTree( *cells3D );

      if( MpiWrapper::min( redistributed3D->GetNumberOfCells(), comm ) == 0 )
      {
        redistributed3D = ensureNoEmptyRank( redistributed3D, comm );
      }
    }
  }
  else
  {
    // All ranks already have cells - no redistribution needed
    redistributed3D = cells3D;
  }

  // Check all ranks have cells after redistribution
  GEOS_ERROR_IF( redistributed3D->GetNumberOfCells() == 0, "Rank has no cells after initial redistribution." );


  // -----------------------------------------------------------------------
  // Step 6: Refined partitioning
  // -----------------------------------------------------------------------
  if( partitionRefinement > 0 )
  {
    if( hasSuperCells )
    {
      GEOS_LOG_RANK_0( "Refining partition with super-cell constraints..." );

      redistributed3D = redistributeBySuperCellGraph(
        redistributed3D,
        method,
        comm,
        partitionRefinement - 1,
        partitionFractureWeight );
    }
    else if( !structuredIndexAttributeName.empty() )
    {
      // Use structured mesh layered partitioning (no fractures/super-cells)
      GEOS_LOG_RANK_0( "Refining partition with structured mesh layers..." );

      // Wrap in AllMeshes for compatibility with redistributeByAreaGraphAndLayer
      AllMeshes tempWrapper;
      tempWrapper.setMainMesh( redistributed3D );
      tempWrapper.setFaceBlocks( stdMap< string, vtkSmartPointer< vtkDataSet > >() );  // Empty fractures

      tempWrapper = redistributeByAreaGraphAndLayer(
        tempWrapper,
        method,
        structuredIndexAttributeName,
        comm,
        numPartZ,
        partitionRefinement - 1 );

      redistributed3D = tempWrapper.getMainMesh();
    }
    else
    {
      // Use standard cell graph partitioning
      GEOS_LOG_RANK_0( "Refining partition with standard cell graph..." );

      redistributed3D = redistributeByCellGraph(
        redistributed3D,
        method,
        comm,
        partitionRefinement - 1 );
    }
  }

  // -----------------------------------------------------------------------
  // Step 7: Merge 2D cells and fractures back with redistributed 3D cells
  // -----------------------------------------------------------------------
  AllMeshes finalResult = redistribute2DAndMergeWith3D(
    redistributed3D,
    mesh,
    cells2DIndices.toViewConst(),
    neighbors2Dto3D,
    namesToFractures,
    fractureNeighbors,
    fractureNames,
    comm );

  // -----------------------------------------------------------------------
  // Step 8: Final logging
  // -----------------------------------------------------------------------
  if( logLevel >= 5 )
  {
    vtkIdType local2DCells = 0;
    vtkIdType local3DCells = 0;
    vtkIdType localFractureCells = 0;

    vtkSmartPointer< vtkDataSet > finalMesh = finalResult.getMainMesh();
    for( vtkIdType i = 0; i < finalMesh->GetNumberOfCells(); ++i )
    {
      int const dim = finalMesh->GetCell( i )->GetCellDimension();
      if( dim == 2 )
        local2DCells++;
      else if( dim == 3 )
        local3DCells++;
    }

    for( auto const & [faceName, faceMesh] : finalResult.getFaceBlocks() )
    {
      localFractureCells += faceMesh->GetNumberOfCells();
    }

    array1d< vtkIdType > all2D, all3D, allFracture;
    MpiWrapper::allGather( local2DCells, all2D, comm );
    MpiWrapper::allGather( local3DCells, all3D, comm );
    MpiWrapper::allGather( localFractureCells, allFracture, comm );

    if( rank == 0 )
    {
      GEOS_LOG_RANK_0( "\n-------------------------------------------------------------" );
      GEOS_LOG_RANK_0( "| Rank  |  3D Cells |  2D Cells | Fractures |     Total     |" );
      GEOS_LOG_RANK_0( "|-----------------------------------------------------------|" );
      vtkIdType sum2D = 0, sum3D = 0, sumFracture = 0;
      for( int r = 0; r < numRanks; ++r )
      {
        sum2D += all2D[r];
        sum3D += all3D[r];
        sumFracture += allFracture[r];

        GEOS_LOG_RANK_0( GEOS_FMT( "| {:>5} | {:9} | {:9} | {:9} | {:13} |",
                                   r, all3D[r], all2D[r], allFracture[r],
                                   all3D[r] + all2D[r] + allFracture[r] ) );
      }
      GEOS_LOG_RANK_0( "|-----------------------------------------------------------|" );
      GEOS_LOG_RANK_0( GEOS_FMT( "| Total | {:9} | {:9} | {:9} | {:13} |",
                                 sum3D, sum2D, sumFracture,
                                 sum3D + sum2D + sumFracture ) );
      GEOS_LOG_RANK_0( "-------------------------------------------------------------" );
    }
  }

  return finalResult;
}


/**
 * @brief Identify the GEOSX type of the polyhedron
 *
 * @param cell The vtk cell VTK_POLYHEDRON
 * @return The geos element type associated to VTK_POLYHEDRON
 */
geos::ElementType buildGeosxPolyhedronType( vtkCell * const cell )
{
  GEOS_ERROR_IF_NE_MSG( cell->GetCellType(), VTK_POLYHEDRON, "Input for polyhedronType() must be a VTK_POLYHEDRON." );

  localIndex const numFaces = cell->GetNumberOfFaces();
  localIndex numTriangles = 0;
  localIndex numQuads = 0;

  // Compute the number of triangles and quads
  for( localIndex iFace = 0; iFace < numFaces; ++iFace )
  {
    localIndex numFaceNodes = cell->GetFace( iFace )->GetNumberOfPoints();
    if( numFaceNodes == 3 )
    {
      numTriangles++;
    }
    if( numFaceNodes == 4 )
    {
      numQuads++;
    }
  }

  if( numTriangles == 4 && numFaces == 4 )
  {
    return geos::ElementType::Tetrahedron;
  }

  if( numQuads == 6 && numFaces == 6 )
  {
    return geos::ElementType::Hexahedron;
  }

  if( numTriangles == 2 && numQuads == 3 && numFaces == 5 )
  {
    return geos::ElementType::Wedge;
  }

  if( numTriangles == 4 && numQuads == 1 && numFaces == 5 )
  {
    return geos::ElementType::Pyramid;
  }

  if( numFaces - numQuads != 2 )
  {
    return geos::ElementType::Polyhedron;
  }

  // Check if the polyhedron is a prism
  // quadsPoints contains points defining all the quads
  // noQuadsPoints contains points defining all the faces which are not quad
  set< localIndex > quadsPoints;
  set< localIndex > noQuadsPoints;
  for( localIndex iFace = 0; iFace < numFaces; ++iFace )
  {
    vtkCell *cellFace = cell->GetFace( iFace );
    if( cellFace->GetNumberOfPoints() == 4 )
    {
      for( localIndex iPoint = 0; iPoint < 4; ++iPoint )
      {
        quadsPoints.insert( cellFace->GetPointId( iPoint ) );
      }
    }
    else
    {
      for( localIndex iPoint = 0; iPoint < cellFace->GetNumberOfPoints(); ++iPoint )
      {
        noQuadsPoints.insert( cellFace->GetPointId( iPoint ) );
      }
    }
  }

  if( quadsPoints != noQuadsPoints )
  {
    return geos::ElementType::Polyhedron;
  }

  // The polyhedron is a prism
  switch( numQuads )
  {
    case 5:  return geos::ElementType::Prism5;
    case 6:  return geos::ElementType::Prism6;
    case 7:  return geos::ElementType::Prism7;
    case 8:  return geos::ElementType::Prism8;
    case 9:  return geos::ElementType::Prism9;
    case 10: return geos::ElementType::Prism10;
    case 11: return geos::ElementType::Prism11;
    default:
    {
      GEOS_ERROR( GEOS_FMT( "Prism with {} sides is not supported.\n{}", numQuads, generalMeshErrorAdvice ) );
      return{};
    }
  }
}

/**
 * @brief Get the GEOSX element type
 * @param[in] cell The vtk cell type
 * @return The GEOSX element type
 */
ElementType convertVtkToGeosxElementType( vtkCell *cell )
{
  switch( cell->GetCellType() )
  {
    case VTK_VERTEX:           return ElementType::Vertex;
    case VTK_LINE:             return ElementType::Line;
    case VTK_TRIANGLE:         return ElementType::Triangle;
    case VTK_QUAD:             return ElementType::Quadrilateral;
    case VTK_POLYGON:          return ElementType::Polygon;
    case VTK_TETRA:            return ElementType::Tetrahedron;
    case VTK_PYRAMID:          return ElementType::Pyramid;
    case VTK_WEDGE:            return ElementType::Wedge;
    case VTK_VOXEL:            return ElementType::Hexahedron;
    case VTK_HEXAHEDRON:       return ElementType::Hexahedron;
    case VTK_PENTAGONAL_PRISM: return ElementType::Prism5;
    case VTK_HEXAGONAL_PRISM:  return ElementType::Prism6;
    case VTK_POLYHEDRON:       return buildGeosxPolyhedronType( cell );
    default:
    {
      GEOS_ERROR( GEOS_FMT( "{} is not a recognized cell type to be used with the VTKMeshGenerator.\n{}",
                            cell->GetCellType(),
                            generalMeshErrorAdvice ) );
      return {};
    }
  }
}

/**
 * @brief Split and arrange the cells of a grid by type
 *
 * @param[in] mesh a vtk grid
 * @return a map of cells grouped by type
 */
stdMap< ElementType, stdVector< vtkIdType > >
splitCellsByType( vtkDataSet & mesh )
{
  stdMap< ElementType, stdVector< vtkIdType > > typeToCells;
  vtkIdType const numCells = mesh.GetNumberOfCells();

  // Count the number of each cell type
  std::array< size_t, numElementTypes() > cellTypeCounts{};
  for( vtkIdType c = 0; c < numCells; c++ )
  {
    ElementType const elemType = convertVtkToGeosxElementType( mesh.GetCell( c ) );
    ++cellTypeCounts[static_cast< integer >( elemType )];
  }

  // Allocate space to hold cell id lists by type
  std::array< stdVector< vtkIdType >, numElementTypes() > cellListsByType;
  for( integer t = 0; t < numElementTypes(); ++t )
  {
    cellListsByType[t].reserve( cellTypeCounts[t] );
  }

  // Collect cell lists for each type (using array for speed in a hot loop)
  for( vtkIdType c = 0; c < numCells; c++ )
  {
    ElementType const elemType = convertVtkToGeosxElementType( mesh.GetCell( c ) );
    cellListsByType[static_cast< integer >( elemType )].push_back( c );
  }

  // Convert from array to map while also performing some checks
  for( integer t = 0; t < numElementTypes(); ++t )
  {
    // Avoid creating unneeded map entries that will show up in statistics
    if( cellListsByType[t].empty() )
    {
      continue;
    }

    ElementType const type = static_cast< ElementType >( t );
    switch( getElementDim( type ) )
    {
      case 0:
      case 1:
      {
        // Ignore vertex/line elements for now; maybe later we can import well polylines here
        break;
      }
      case 2:
      {
        // Merge all 2D elements together as polygons (we don't track their shapes).
        stdVector< vtkIdType > & surfaceCells = typeToCells.get_inserted( ElementType::Polygon );
        surfaceCells.insert( surfaceCells.end(), cellListsByType[t].begin(), cellListsByType[t].end() );
        break;
      }
      case 3:
      {
        // Collect 3D elements as is
        typeToCells.emplace( type, std::move( cellListsByType[t] ) );
        break;
      }
      default:
      {
        GEOS_ERROR( GEOS_FMT( "Invalid element dimension: {}", getElementDim( type ) ) );
      }
    }
  }

  return typeToCells;
}

/**
 * @brief Split and arrange the cells of a grid according to an attribute
 *
 * @param[in] typeToCells a map of cells grouped by type
 * @param attributeDataArray an attribute
 * @return a map of cell lists grouped by type
 */
CellMapType
splitCellsByTypeAndAttribute( stdMap< ElementType, stdVector< vtkIdType > > & typeToCells,
                              vtkDataArray * const attributeDataArray )
{
  CellMapType typeToAttributeToCells;
  for( auto & t2c : typeToCells )
  {
    ElementType const elemType = t2c.first;
    stdVector< vtkIdType > & cells = t2c.second;
    stdUnorderedMap< int, stdVector< vtkIdType > > & attributeToCells = typeToAttributeToCells.get_inserted( elemType );

    if( attributeDataArray == nullptr )
    {
      attributeToCells.emplace( -1, std::move( cells ) );
    }
    else
    {
      GEOS_ERROR_IF_NE_MSG( attributeDataArray->GetNumberOfComponents(), 1,
                            "Invalid number of components in attribute array" );
      vtkArrayDispatch::Dispatch::Execute( attributeDataArray, [&]( auto const * const attributeArray )
      {
        using ArrayType = TYPEOFPTR( attributeArray );
        vtkDataArrayAccessor< ArrayType > attribute( attributeArray );
        stdUnorderedMap< int, size_t > cellCounts;
        for( vtkIdType c: cells )
        {
          int const region = static_cast< int >( attribute.Get( c, 0 ) );
          ++cellCounts.get_inserted( region );
        }
        for( auto const & count : cellCounts )
        {
          attributeToCells.get_inserted( count.first ).reserve( count.second );
        }
        for( vtkIdType c: cells )
        {
          int const region = static_cast< int >( attribute.Get( c, 0 ) );
          attributeToCells.get_inserted( region ).push_back( c );
        }
      } );
    }
  }
  return typeToAttributeToCells;
}

/**
 * @brief Gather all element types encountered on any rank and enrich the local collection
 *
 * @param[in] cellMap a map of cell lists grouped by type
 */
void extendCellMapWithRemoteKeys( CellMapType & cellMap )
{
  // Gather all element types encountered on any rank and enrich the local collection
  stdVector< ElementType > allElementTypes = collectUniqueValues( mapKeys( cellMap ) );
  stdVector< int > allCellAttributes;
  for( auto const & typeRegions : cellMap )
  {
    if( getElementDim( typeRegions.first ) == 3 )
    {
      stdVector< int > const attrs = mapKeys( typeRegions.second );
      allCellAttributes.insert( allCellAttributes.end(), attrs.begin(), attrs.end() );
    }
  }
  allCellAttributes = collectUniqueValues( allCellAttributes );

  for( ElementType const elemType : allElementTypes )
  {
    if( getElementDim( elemType ) == 3 )
    {
      for( int attrValue: allCellAttributes )
      {
        // This code inserts an empty element list if one was not present
        cellMap.get_inserted( elemType ).get_inserted( attrValue );
      }
    }
  }

  // Treat surfaces separately - and avoid inadvertently creating a map entry for polygons
  stdVector< int > const surfaceAttributes = cellMap.count( ElementType::Polygon ) > 0
                                             ? mapKeys( cellMap.at( ElementType::Polygon ) )
                                             : stdVector< int >();
  stdVector< int > allSurfaceAttributes = collectUniqueValues( surfaceAttributes );
  for( int attrValue: allSurfaceAttributes )
  {
    cellMap.get_inserted( ElementType::Polygon ).get_inserted( attrValue );
  }
}

/**
 * @brief Get the tetrahedron node ordering from a VTK_POLYHEDRON
 *
 * @param cell The vtk cell, type VTK_POLYHEDRON
 * @return The node ordering
 */
stdVector< localIndex > getTetrahedronNodeOrderingFromPolyhedron( vtkCell * const cell )
{
  GEOS_ERROR_IF_NE_MSG( cell->GetCellType(), VTK_POLYHEDRON, "Input must be a VTK_POLYHEDRON." );

  real64 vectorA[3];
  real64 vectorB[3];
  real64 vectorC[3];
  real64 vectorN[3];
  vtkPoints *cellPoint = cell->GetPoints();

  // Check if the orientation is correct
  for( localIndex i = 0; i < 3; ++i )
  {
    vectorA[i] = cellPoint->GetPoint( 1 )[i] - cellPoint->GetPoint( 0 )[i];
    vectorB[i] = cellPoint->GetPoint( 2 )[i] - cellPoint->GetPoint( 0 )[i];
    vectorC[i] = cellPoint->GetPoint( 3 )[i] - cellPoint->GetPoint( 0 )[i];
  }

  LvArray::tensorOps::crossProduct( vectorN, vectorA, vectorB );

  if( LvArray::tensorOps::AiBi< 3 >( vectorN, vectorC ) > 0 )
  {
    // The orientation is correct
    return { 0, 1, 2, 3 };
  }
  else
  {
    // The orientation is incorrect, renumber nodes
    return { 0, 2, 1, 3 };
  }
}

/**
 * @brief Get the hexahedron node ordering from a VTK_POLYHEDRON
 *
 * @param cell The vtk cell, type VTK_POLYHEDRON
 * @return The node ordering
 *
 * It could be possible to use getPrismNodeOrderingFromPolyhedron< 4 > with additional
 * permutations. But at this point computationalGeometry::prismVolume< NUM_SIDES >
 * is not ready.
 */
stdVector< localIndex > getHexahedronNodeOrderingFromPolyhedron( vtkCell * const cell )
{
  GEOS_ERROR_IF_NE_MSG( cell->GetCellType(), VTK_POLYHEDRON, "Input must be a VTK_POLYHEDRON." );

  localIndex iFace;
  stdVector< localIndex > nodeOrder( 8 );

  // Generate global to local map
  stdUnorderedMap< localIndex, localIndex > G2L;
  for( localIndex iPoint = 0; iPoint < 8; ++iPoint )
  {
    G2L.insert( {cell->GetPointId( iPoint ), iPoint} );
  }

  // Assuming the input parameters are correct, take the first quad
  iFace = 0;

  // Get global pointIds for the first quad
  vtkCell *cellFace = cell->GetFace( iFace );
  nodeOrder[0] = cellFace->GetPointId( 0 );
  nodeOrder[1] = cellFace->GetPointId( 1 );
  nodeOrder[2] = cellFace->GetPointId( 3 );
  nodeOrder[3] = cellFace->GetPointId( 2 );

  // Generate the global pointIds for the opposit quad using the edges connecting the two bases
  for( localIndex iEdge = 0; iEdge < 12; ++iEdge )
  {
    auto edgeNode0 = cell->GetEdge( iEdge )->GetPointId( 0 );
    auto edgeNode1 = cell->GetEdge( iEdge )->GetPointId( 1 );
    auto it0 = std::find( &nodeOrder[0], &nodeOrder[4], edgeNode0 );
    auto it1 = std::find( &nodeOrder[0], &nodeOrder[4], edgeNode1 );

    if( it0 != &nodeOrder[4] && it1 == &nodeOrder[4] )
    {
      nodeOrder[ 4 + std::distance( &nodeOrder[0], it0 ) ] = edgeNode1;
    }

    if( it0 == &nodeOrder[4] && it1 != &nodeOrder[4] )
    {
      nodeOrder[ 4 + std::distance( &nodeOrder[0], it1 ) ] = edgeNode0;
    }
  }

  // Convert global numbering to local numbering
  for( localIndex iPoint = 0; iPoint < 8; ++iPoint )
  {
    nodeOrder[iPoint] = G2L.at( nodeOrder[iPoint] );
  }

  real64 Xlocal[8][3];
  for( localIndex iPoint = 0; iPoint < 8; ++iPoint )
  {
    std::copy_n( cell->GetPoints()->GetPoint( nodeOrder[iPoint] ), 3, Xlocal[iPoint] );
  }
  real64 const cellVolume = computationalGeometry::hexahedronVolume( Xlocal );

  // If cell volume is negative swap the quads
  if( cellVolume < 0 )
  {
    std::rotate( nodeOrder.begin(), nodeOrder.begin()+4, nodeOrder.end());
  }

  return nodeOrder;
}

/**
 * @brief Get the wedge node ordering from a VTK_POLYHEDRON
 *
 * @param cell The vtk cell, type VTK_POLYHEDRON
 * @return The node ordering
 *
 * It could be possible to use getPrismNodeOrderingFromPolyhedron< 3 > with additional
 * permutations. But at this point computationalGeometry::prismVolume< NUM_SIDES >
 * is not ready.
 */
stdVector< localIndex > getWedgeNodeOrderingFromPolyhedron( vtkCell * const cell )
{
  GEOS_ERROR_IF_NE_MSG( cell->GetCellType(), VTK_POLYHEDRON, "Input must be a VTK_POLYHEDRON." );

  localIndex iFace;
  stdVector< localIndex > nodeTri0( 3 );
  stdVector< localIndex > nodeTri1( 3 );
  stdVector< localIndex > nodeOrder( 6 );

  // Generate global to local map
  stdUnorderedMap< localIndex, localIndex > G2L;
  for( localIndex iPoint = 0; iPoint < 6; ++iPoint )
  {
    G2L.get_inserted( cell->GetPointId( iPoint )) = iPoint;
  }

  // Assuming the input parameters are correct, identify one of the triangles
  localIndex const numFaces = cell->GetNumberOfFaces();
  for( iFace = 0; iFace < numFaces; ++iFace )
  {
    if( cell->GetFace( iFace )->GetNumberOfPoints() == 3 )
    {
      break;
    }
  }

  GEOS_ERROR_IF( iFace == numFaces, GEOS_FMT( "Invalid wedge.\n{}", generalMeshErrorAdvice ) );

  // Get global pointIds for the first triangle
  for( localIndex i = 0; i < 3; ++i )
  {
    nodeTri0[i] = cell->GetFace( iFace )->GetPointId( i );
  }

  // Generate the global pointIds for the second triangle using the edges connecting the two triangles
  for( localIndex iEdge = 0; iEdge < 9; ++iEdge )
  {
    auto edgeNode0 = cell->GetEdge( iEdge )->GetPointId( 0 );
    auto edgeNode1 = cell->GetEdge( iEdge )->GetPointId( 1 );
    auto it0 = std::find( nodeTri0.begin(), nodeTri0.end(), edgeNode0 );
    auto it1 = std::find( nodeTri0.begin(), nodeTri0.end(), edgeNode1 );

    if( it0 != nodeTri0.end() && it1 == nodeTri0.end() )
    {
      nodeTri1[std::distance( nodeTri0.begin(), it0 )] = edgeNode1;
    }

    if( it0 == nodeTri0.end() && it1 != nodeTri0.end() )
    {
      nodeTri1[std::distance( nodeTri0.begin(), it1 )] = edgeNode0;
    }
  }

  // Convert global numbering to local numbering
  for( int iPoint = 0; iPoint < 3; ++iPoint )
  {
    nodeTri0[iPoint] = G2L.at( nodeTri0[iPoint] );
    nodeTri1[iPoint] = G2L.at( nodeTri1[iPoint] );
  }

  // Set nodes numbering
  nodeOrder[0] = nodeTri0[0];
  nodeOrder[1] = nodeTri1[0];
  nodeOrder[2] = nodeTri0[1];
  nodeOrder[3] = nodeTri1[1];
  nodeOrder[4] = nodeTri0[2];
  nodeOrder[5] = nodeTri1[2];

  // Compute cell volume
  real64 Xlocal[6][3];
  for( localIndex iPoint = 0; iPoint < 6; ++iPoint )
  {
    std::copy_n( cell->GetPoints()->GetPoint( nodeOrder[iPoint] ), 3, Xlocal[iPoint] );
  }
  real64 const cellVolume = computationalGeometry::wedgeVolume( Xlocal );

  // If cell volume is negative reorder nodes
  if( cellVolume < 0 )
  {
    nodeOrder[0] = nodeTri0[0];
    nodeOrder[1] = nodeTri1[0];
    nodeOrder[2] = nodeTri0[2];
    nodeOrder[3] = nodeTri1[2];
    nodeOrder[4] = nodeTri0[1];
    nodeOrder[5] = nodeTri1[1];
  }

  return nodeOrder;
}

/**
 * @brief Get the pyramid node ordering from a VTK_POLYHEDRON
 *
 * @param cell The vtk cell, type VTK_POLYHEDRON
 * @return The node ordering
 */
stdVector< localIndex > getPyramidNodeOrderingFromPolyhedron( vtkCell * const cell )
{
  GEOS_ERROR_IF_NE_MSG( cell->GetCellType(), VTK_POLYHEDRON, "Input must be a VTK_POLYHEDRON." );

  localIndex iPoint;
  localIndex iFace;
  stdVector< localIndex > nodeOrder( 5 );

  // Generate global to local map
  stdUnorderedMap< localIndex, localIndex > G2L;
  for( iPoint = 0; iPoint < 5; ++iPoint )
  {
    G2L.get_inserted( cell->GetPointId( iPoint )) = iPoint;
  }

  // Assuming the input parameters are correct, identify the base
  localIndex const numFaces = cell->GetNumberOfFaces();
  for( iFace = 0; iFace < numFaces; ++iFace )
  {
    if( cell->GetFace( iFace )->GetNumberOfPoints() == 4 )
    {
      break;
    }
  }

  GEOS_ERROR_IF( iFace == numFaces, GEOS_FMT( "Invalid pyramid.\n{}", generalMeshErrorAdvice ) );

  // Get global pointIds for the base
  vtkCell * cellFace = cell->GetFace( iFace );
  nodeOrder[0] = cellFace->GetPointId( 0 );
  nodeOrder[1] = cellFace->GetPointId( 1 );
  nodeOrder[2] = cellFace->GetPointId( 3 );
  nodeOrder[3] = cellFace->GetPointId( 2 );

  // Identify the missing node
  iPoint = 0;
  while( std::find( &nodeOrder[0], &nodeOrder[4], cell->GetPointId( iPoint ) ) != &nodeOrder[4] )
  {
    iPoint++;
  }

  // Add it to nodeOrder
  nodeOrder[4] = cell->GetPointId( iPoint );

  // Convert global numbering to local numbering
  for( iPoint = 0; iPoint < 5; ++iPoint )
  {
    nodeOrder[iPoint] = G2L.at( nodeOrder[iPoint] );
  }

  // Compute cell volume
  real64 Xlocal[5][3];
  for( iPoint = 0; iPoint < 5; ++iPoint )
  {
    std::copy_n( cell->GetPoints()->GetPoint( nodeOrder[iPoint] ), 3, Xlocal[iPoint] );
  }
  real64 const cellVolume = computationalGeometry::pyramidVolume( Xlocal );

  // If cell volume is negative swap two nodes from the base
  if( cellVolume < 0 )
  {
    std::swap( nodeOrder[1], nodeOrder[2] );
  }

  return nodeOrder;
}

/**
 * @brief Get the prism node ordering from a VTK_POLYHEDRON
 *
 * @tparam NUM_SIDES number of sides of the prism
 * @param cell The vtk cell, type VTK_POLYHEDRON
 * @return The node ordering
 */
template< integer NUM_SIDES >
stdVector< localIndex > getPrismNodeOrderingFromPolyhedron( vtkCell * const cell )
{
  GEOS_ERROR_IF_NE_MSG( cell->GetCellType(), VTK_POLYHEDRON, "Input must be a VTK_POLYHEDRON." );

  localIndex iFace;
  stdVector< localIndex > nodeOrder( 2*NUM_SIDES );

  // Generate global to local map
  stdUnorderedMap< localIndex, localIndex > G2L;
  for( localIndex iPoint = 0; iPoint < cell->GetNumberOfPoints(); ++iPoint )
  {
    G2L.get_inserted( cell->GetPointId( iPoint )) = iPoint;
  }

  // Assuming the input parameters are correct, identify one of the bases
  localIndex const numFaces = cell->GetNumberOfFaces();
  for( iFace = 0; iFace < numFaces; ++iFace )
  {
    if( cell->GetFace( iFace )->GetNumberOfPoints() == NUM_SIDES )
    {
      break;
    }
  }

  GEOS_ERROR_IF( iFace == numFaces, GEOS_FMT( "Invalid prism.\n{}", generalMeshErrorAdvice ) );

  // Get global pointIds for the first base
  vtkCell *cellFace = cell->GetFace( iFace );
  for( localIndex iPoint = 0; iPoint < NUM_SIDES; ++iPoint )
  {
    nodeOrder[iPoint] = cellFace->GetPointId( iPoint );
  }

  // Generate the global pointIds for the second base using the edges connecting the two bases
  localIndex numEdges = cell->GetNumberOfEdges();
  for( localIndex iEdge = 0; iEdge < numEdges; ++iEdge )
  {
    auto edgeNode0 = cell->GetEdge( iEdge )->GetPointId( 0 );
    auto edgeNode1 = cell->GetEdge( iEdge )->GetPointId( 1 );
    auto it0 = std::find( &nodeOrder[0], &nodeOrder[NUM_SIDES], edgeNode0 );
    auto it1 = std::find( &nodeOrder[0], &nodeOrder[NUM_SIDES], edgeNode1 );

    if( it0 != &nodeOrder[NUM_SIDES] && it1 == &nodeOrder[NUM_SIDES] )
    {
      nodeOrder[NUM_SIDES + std::distance( &nodeOrder[0], it0 )] = edgeNode1;
    }

    if( it0 == &nodeOrder[NUM_SIDES] && it1 != &nodeOrder[NUM_SIDES] )
    {
      nodeOrder[NUM_SIDES + std::distance( &nodeOrder[0], it1 )] = edgeNode0;
    }
  }

  // Convert global numbering to local numbering
  for( localIndex iPoint = 0; iPoint < 2*NUM_SIDES; ++iPoint )
  {
    nodeOrder[iPoint] = G2L.at( nodeOrder[iPoint] );
  }

  // Compute cell volume
  real64 Xlocal[2*NUM_SIDES][3];
  for( localIndex iPoint = 0; iPoint < 2*NUM_SIDES; ++iPoint )
  {
    std::copy_n( cell->GetPoints()->GetPoint( nodeOrder[iPoint] ), 3, Xlocal[iPoint] );
  }
  real64 const cellVolume = computationalGeometry::prismVolume< NUM_SIDES >( Xlocal );

  // If cell volume is negative swap the bases
  if( cellVolume < 0 )
  {
    std::rotate( nodeOrder.begin(), nodeOrder.begin()+NUM_SIDES, nodeOrder.end());
  }

  return nodeOrder;
}

CellMapType buildCellMap( vtkDataSet & mesh, string const & attributeName )
{

  // First, pass through all VTK cells and split them int sub-lists based on type.
  stdMap< ElementType, stdVector< vtkIdType > > typeToCells = splitCellsByType( mesh );

  // Now, actually split into groups according to region attribute, if present
  vtkDataArray * const attributeDataArray =
    vtkDataArray::FastDownCast( mesh.GetCellData()->GetAbstractArray( attributeName.c_str() ) );

  CellMapType cellMap =
    splitCellsByTypeAndAttribute( typeToCells, attributeDataArray );

  // Gather all element types encountered on any rank and enrich the local collection
  extendCellMapWithRemoteKeys( cellMap );

  return cellMap;
}

bool vtkToGeosxNodeOrderingExists( ElementType const elemType )
{
  switch( elemType )
  {
    case ElementType::Vertex:
    case ElementType::Line:
    case ElementType::Triangle:
    case ElementType::Quadrilateral:
    case ElementType::Tetrahedron:
    case ElementType::Pyramid:
    case ElementType::Wedge:
    case ElementType::Hexahedron:
    case ElementType::Prism5:
    case ElementType::Prism6:
    {
      return true;
    }
    default:
    {
      return false;
    }
  }
}

stdVector< int > getVtkToGeosxNodeOrdering( ElementType const elemType )
{
  switch( elemType )
  {
    case ElementType::Vertex:        return { 0 };
    case ElementType::Line:          return { 0, 1 };
    case ElementType::Triangle:      return { 0, 1, 2 };
    case ElementType::Quadrilateral: return { 0, 1, 2, 3 };  // TODO check (QUAD vs. PIXEL)
    case ElementType::Tetrahedron:   return { 0, 1, 2, 3 };
    case ElementType::Pyramid:       return { 0, 1, 3, 2, 4 };
    case ElementType::Wedge:         return { 0, 3, 2, 5, 1, 4 };
    case ElementType::Hexahedron:    return { 0, 1, 3, 2, 4, 5, 7, 6 };
    case ElementType::Prism5:        return { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    case ElementType::Prism6:        return { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
    default:
    {
      GEOS_ERROR( GEOS_FMT( "Cannot get vtk to geos node ordering based on geos element type {}", elemType ) );
      break;
    }
  }
  return {};
}

stdVector< int > getVtkToGeosxNodeOrdering( VTKCellType const vtkType )
{
  switch( vtkType )
  {
    case VTK_VERTEX:           return getVtkToGeosxNodeOrdering( ElementType::Vertex );
    case VTK_LINE:             return getVtkToGeosxNodeOrdering( ElementType::Line );
    case VTK_TRIANGLE:         return getVtkToGeosxNodeOrdering( ElementType::Triangle );
    case VTK_QUAD:             return getVtkToGeosxNodeOrdering( ElementType::Quadrilateral );
    case VTK_TETRA:            return getVtkToGeosxNodeOrdering( ElementType::Tetrahedron );
    case VTK_PYRAMID:          return getVtkToGeosxNodeOrdering( ElementType::Pyramid );
    case VTK_WEDGE:            return getVtkToGeosxNodeOrdering( ElementType::Wedge );
    case VTK_VOXEL:            return { 0, 1, 2, 3, 4, 5, 6, 7 };
    case VTK_HEXAHEDRON:       return getVtkToGeosxNodeOrdering( ElementType::Hexahedron );
    case VTK_PENTAGONAL_PRISM: return getVtkToGeosxNodeOrdering( ElementType::Prism5 );
    case VTK_HEXAGONAL_PRISM:  return getVtkToGeosxNodeOrdering( ElementType::Prism6 );
    default:
    {
      GEOS_ERROR( GEOS_FMT( "Cannot get vtk to geos node ordering based on vtk cell type {}",
                            static_cast< int >( vtkType ) ) );
      break;
    }
  }
  return {};
}

stdVector< int > getVtkToGeosxPolyhedronNodeOrdering( ElementType const elemType,
                                                      vtkCell *cell )
{
  GEOS_ERROR_IF_NE_MSG( cell->GetCellType(), VTK_POLYHEDRON, "Input must be a VTK_POLYHEDRON." );
  switch( elemType )
  {
    case ElementType::Tetrahedron: return getTetrahedronNodeOrderingFromPolyhedron( cell );
    case ElementType::Pyramid:     return getPyramidNodeOrderingFromPolyhedron( cell );
    case ElementType::Wedge:       return getWedgeNodeOrderingFromPolyhedron( cell );
    case ElementType::Hexahedron:  return getHexahedronNodeOrderingFromPolyhedron( cell );
    case ElementType::Prism5:      return getPrismNodeOrderingFromPolyhedron< 5 >( cell );
    case ElementType::Prism6:      return getPrismNodeOrderingFromPolyhedron< 6 >( cell );
    case ElementType::Prism7:      return getPrismNodeOrderingFromPolyhedron< 7 >( cell );
    case ElementType::Prism8:      return getPrismNodeOrderingFromPolyhedron< 8 >( cell );
    case ElementType::Prism9:      return getPrismNodeOrderingFromPolyhedron< 9 >( cell );
    case ElementType::Prism10:     return getPrismNodeOrderingFromPolyhedron< 10 >( cell );
    case ElementType::Prism11:     return getPrismNodeOrderingFromPolyhedron< 11 >( cell );
    default:
    {
      GEOS_ERROR( "Unsupported VTK polyhedral cell" );
      break;
    }
  }
  return {};
}

/**
 * @brief Fill @p cellBlock with the appropriate nodes and local/global mappings.
 * @param[in] cellIds the cell indexes of cell type \p cellType within this region
 * @param[in] mesh the vtkUnstructuredGrid or vtkStructuredGrid that is loaded
 * @param[in,out] cellBlock The cell block to be written
 */
void fillCellBlock( vtkDataSet & mesh,
                    Span< vtkIdType const > const cellIds,
                    CellBlock & cellBlock )
{
  GEOS_ASSERT_EQ( LvArray::integerConversion< localIndex >( cellIds.size() ), cellBlock.size() );

  localIndex const numNodesPerElement = cellBlock.numNodesPerElement();
  arrayView2d< localIndex, cells::NODE_MAP_USD > const cellToVertex = cellBlock.getElemToNode();
  arrayView1d< globalIndex > const & localToGlobal = cellBlock.localToGlobalMap();
  vtkIdTypeArray const * const globalCellId = vtkIdTypeArray::FastDownCast( mesh.GetCellData()->GetGlobalIds() );
  GEOS_ERROR_IF( !cellIds.empty() && globalCellId == nullptr, "Global cell IDs have not been generated" );

  localIndex cellCount = 0;
  auto const writeCell = [&]( vtkIdType const c, vtkCell * const cell, auto const & nodeOrder )
  {
    for( localIndex v = 0; v < numNodesPerElement; v++ )
    {
      cellToVertex[cellCount][v] = LvArray::integerConversion< localIndex >( cell->GetPointId( nodeOrder[v] ) );
    }
    localToGlobal[cellCount++] = globalCellId->GetValue( c );
  };

  // Writing connectivity and Local to Global
  ElementType const elemType = cellBlock.getElementType();
  stdVector< int > const nodeOrderFixed = vtkToGeosxNodeOrderingExists( elemType )
                                          ? getVtkToGeosxNodeOrdering( elemType )
                                          : stdVector< int >();
  stdVector< int > const nodeOrderVoxel = ( elemType == ElementType::Hexahedron)
                                          ? getVtkToGeosxNodeOrdering( VTK_VOXEL )
                                          : stdVector< int >();

  for( vtkIdType c: cellIds )
  {
    vtkCell * const cell = mesh.GetCell( c );
    VTKCellType const vtkType = static_cast< VTKCellType >( cell->GetCellType() );
    switch( vtkType )
    {
      case VTK_POLYHEDRON:
      {
        writeCell( c, cell, getVtkToGeosxPolyhedronNodeOrdering( elemType, cell ) );
        break;
      }
      case VTK_VOXEL:
      {
        writeCell( c, cell, nodeOrderVoxel );
        break;
      }
      default:
      {
        writeCell( c, cell, nodeOrderFixed );
        break;
      }
    }
  }
}

void importMaterialField( stdVector< vtkIdType > const & cellIds,
                          vtkDataArray * vtkArray,
                          WrapperBase & wrapper )
{
  // Scalar material fields are stored as 2D arrays, vector/tensor are 3D
  using ImportTypes = types::ListofTypeList< types::ArrayTypes< types::RealTypes, types::DimsRange< 2, 3 > > >;
  types::dispatch( ImportTypes{}, [&]( auto tupleOfTypes )
  {
    using ArrayType = camp::first< decltype( tupleOfTypes ) >;
    Wrapper< ArrayType > & wrapperT = Wrapper< ArrayType >::cast( wrapper );
    auto const view = wrapperT.reference().toView();

    localIndex const numComponentsSrc = LvArray::integerConversion< localIndex >( vtkArray->GetNumberOfComponents() );
    localIndex const numComponentsDst = wrapperT.numArrayComp() / view.size( 1 );
    GEOS_ERROR_IF_NE_MSG( numComponentsDst, numComponentsSrc,
                          GEOS_FMT( "Mismatch in number of components for field {}", vtkArray->GetName() ) );

    vtkArrayDispatch::DispatchByValueType< vtkArrayDispatch::Reals >::Execute( vtkArray, [&]( auto const * srcArray )
    {
      vtkDataArrayAccessor< TYPEOFPTR( srcArray ) > data( srcArray );
      localIndex cellCount = 0;
      for( vtkIdType cellIdx : cellIds )
      {
        for( localIndex q = 0; q < view.size( 1 ); ++q )
        {
          // The same value is copied for all quadrature points.
          LvArray::forValuesInSlice( view[cellCount][q], [&, componentIdx = 0]( auto & val ) mutable
          {
            val = data.Get( cellIdx, componentIdx++ );
          } );
        }
        ++cellCount;
      }
    } );
  }, wrapper );
}

void importRegularField( stdVector< vtkIdType > const & cellIds,
                         vtkDataArray * vtkArray,
                         WrapperBase & wrapper )
{
  using ImportTypes = types::ListofTypeList< types::ArrayTypes< types::RealTypes, types::DimsRange< 1, 2 > > >;
  types::dispatch( ImportTypes{}, [&]( auto tupleOfTypes )
  {
    using ArrayType = camp::first< decltype( tupleOfTypes ) >;
    Wrapper< ArrayType > & wrapperT = Wrapper< ArrayType >::cast( wrapper );
    auto const view = wrapperT.reference().toView();

    localIndex const numComponentsSrc = LvArray::integerConversion< localIndex >( vtkArray->GetNumberOfComponents() );
    localIndex const numComponentsDst = wrapperT.numArrayComp();
    GEOS_ERROR_IF_NE_MSG( numComponentsDst, numComponentsSrc,
                          GEOS_FMT( "Mismatch in number of components for field {}", vtkArray->GetName() ) );

    vtkArrayDispatch::DispatchByValueType< vtkArrayDispatch::Reals >::Execute( vtkArray, [&]( auto const * srcArray )
    {
      vtkDataArrayAccessor< TYPEOFPTR( srcArray ) > data( srcArray );
      localIndex cellCount = 0;
      for( vtkIdType cellIdx : cellIds )
      {
        LvArray::forValuesInSlice( view[cellCount], [&, componentIdx = 0]( auto & val ) mutable
        {
          val = data.Get( cellIdx, componentIdx++ );
        } );
        ++cellCount;
      }
    } );
  }, wrapper );
}


void importRegularField( vtkDataArray * vtkArray,
                         WrapperBase & wrapper )
{
  stdVector< vtkIdType > cellIds( wrapper.size() );
  std::iota( cellIds.begin(), cellIds.end(), 0 );
  return importRegularField( cellIds, vtkArray, wrapper );
}


void printMeshStatistics( vtkDataSet &,
                          CellMapType const & cellMap,
                          MPI_Comm const comm )
{
  auto accumulateElemsCount = []( stdMap< ElementType, globalIndex > & elemsTarget ) -> globalIndex
  {
    return std::accumulate(
      std::begin( elemsTarget ), std::end( elemsTarget ), globalIndex{0},
      []( std::size_t const previous, auto const & elems )
    { return previous + elems.second; } );
  };

  int const rank = MpiWrapper::commRank( comm );
  int const size = MpiWrapper::commSize( comm );

  stdMap< ElementType, globalIndex > totalLocalElems;
  stdMap< ElementType, globalIndex > minLocalElemsCounts;
  stdMap< ElementType, globalIndex > avgLocalElemsCounts;
  stdMap< ElementType, globalIndex > maxLocalElemsCounts;

  for( auto const & typeToCells : cellMap )
  {
    localIndex const localElemsOfType =
      std::accumulate( typeToCells.second.begin(), typeToCells.second.end(), localIndex{},
                       []( auto const s, auto const & region ) { return s + region.second.size(); } );

    totalLocalElems.insert( {typeToCells.first, MpiWrapper::sum( globalIndex{ localElemsOfType }, comm )} );
    minLocalElemsCounts.insert( {typeToCells.first, MpiWrapper::min( localElemsOfType )} );
    avgLocalElemsCounts.insert( {typeToCells.first, LvArray::integerConversion< localIndex >( MpiWrapper::sum( localElemsOfType ) / size )} );
    maxLocalElemsCounts.insert( {typeToCells.first, MpiWrapper::max( localElemsOfType )} );
  }

  if( rank == 0 )
  {

    auto sumOfElemsType = accumulateElemsCount( totalLocalElems );
    auto sumOfMinElemsType = accumulateElemsCount( minLocalElemsCounts );
    auto sumOfAvgElemsType = accumulateElemsCount( avgLocalElemsCounts );
    auto sumOfMaxElemsType = accumulateElemsCount( maxLocalElemsCounts );

    TableLayout const elemsLayout( "Load balancing, element / rank",
                                   { TableLayout::Column()
                                       .setName( "Element type" )
                                       .setValuesAlignment( TableLayout::Alignment::left ),
                                     "Total over ranks",
                                     "minimum",
                                     "average",
                                     "maximum" } );
    TableData elemsData;
    for( auto const & typeCount: totalLocalElems )
    {
      elemsData.addRow( typeCount.first, typeCount.second,
                        minLocalElemsCounts[typeCount.first],
                        avgLocalElemsCounts[typeCount.first],
                        maxLocalElemsCounts[typeCount.first] );
    }
    elemsData.addRow( "total elements", sumOfElemsType, sumOfMinElemsType, sumOfAvgElemsType, sumOfMaxElemsType );
    TableTextFormatter elemsText( elemsLayout );
    GEOS_LOG_RANK_0( elemsText.toString( elemsData ));
  }
}

vtkDataArray *
findArrayForImport( vtkDataSet & mesh,
                    string const & sourceName )
{
  vtkCellData & cellData = *mesh.GetCellData();

  vtkAbstractArray * const curArray = cellData.GetAbstractArray( sourceName.c_str() );
  GEOS_THROW_IF( curArray == nullptr,
                 GEOS_FMT( "Source field '{}' not found in dataset", sourceName ),
                 InputError );

  int const dataType = curArray->GetDataType();
  GEOS_ERROR_IF( dataType != VTK_FLOAT && dataType != VTK_DOUBLE,
                 GEOS_FMT( "Source field '{}' has unsupported type: {} (expected floating point type)",
                           sourceName, curArray->GetDataTypeAsString() ) );
  return vtkDataArray::SafeDownCast( curArray );
}

bool hasArray( vtkDataSet & mesh, string const & sourceName )
{
  return mesh.GetCellData()->GetAbstractArray( sourceName.c_str() ) != nullptr;
}

/**
 * @brief Builds the cell block name.
 * @param[in] type The element name.
 * @param[in] regionId The region Id.
 * @return The name.
 * @warning This name will be visible in the input file... Consider refactoring with great care.
 */
string buildCellBlockName( ElementType const type, int const regionId )
{
  GEOS_ERROR_IF_LT_MSG( regionId, -1, "Invalid region id" );
  string const cellTypeName = getElementTypeName( type );
  return regionId != -1 ? std::to_string( regionId ) + "_" + cellTypeName : cellTypeName;
}

} // namespace vtk

namespace
{

/**
 * @brief Extract node indices from a binary mask array using parallel scan
 *
 * @tparam ValueType The underlying data type (vtkTypeUInt8, vtkTypeUInt16, vtkTypeUInt32, vtkTypeUInt64)
 * @param[in] rawData Pointer to the raw mask array data (values should be 0 or 1)
 * @param[in] numPoints Total number of points in the mesh
 * @param[out] targetNodeset The sorted array to populate with node indices
 *
 * @note Input array should contain only 0 (not in set) or 1 (in set). Values other than
 *       1 will be treated as 0.
 */
template< typename ValueType >
void extractNodesetFromMask( ValueType const * rawData,
                             localIndex const numPoints,
                             SortedArray< localIndex > & targetNodeset )
{
  // Early exit for empty mesh
  if( numPoints == 0 )
  {
    return;
  }

  // Allocate temporary array for positions (initially stores flags)
  array1d< localIndex > positions( numPoints );

  // Step 1: Extract binary flags (0 or 1) from mask (parallel)
  forAll< parallelHostPolicy >( numPoints, [rawData, &positions]( localIndex const j )
  {
    positions[j] = ( rawData[j] == 1 );  // Expect exactly 1 for membership
  } );

  // Save last flag value before scan overwrites it
  localIndex const lastFlag = positions[numPoints - 1];

  // Step 2: In-place exclusive scan to compute write positions (parallel)
  RAJA::exclusive_scan_inplace< parallelHostPolicy >(
    RAJA::make_span( positions.data(), numPoints ),
    RAJA::operators::plus< localIndex >{}
    );

  // Step 3: Compute total count from scan result
  localIndex const lastPosition = positions[numPoints - 1];
  localIndex const count = lastPosition + lastFlag;

  if( count == 0 )
  {
    return;  // No nodes in this nodeset
  }

  // Step 4: Allocate exact-size output array
  array1d< localIndex > nodeIndices( count );

  // Step 5: Scatter flagged node indices to compacted array (parallel)
  forAll< parallelHostPolicy >( numPoints, [rawData, &nodeIndices, &positions]( localIndex const j )
  {
    if( rawData[j] == 1 )
    {
      nodeIndices[ positions[j] ] = j;
    }
  } );

  // Step 6: Insert into sorted nodeset container exploiting sortedness and uniqueness of nodeIndices
  targetNodeset.insert( nodeIndices.begin(), nodeIndices.end() );
}

}

/**
 * @brief Build node sets from binary mask arrays
 *
 * @param[in] logLevel The log level
 * @param[in] mesh The vtk grid that is loaded
 * @param[in] nodesetNames An array of the node sets names
 * @param[in] cellBlockManager The instance that stores the node sets.
 *
 * @note Nodeset arrays should be binary masks (0 or 1 values) stored as unsigned integers.
 *       UInt8 is recommended for optimal memory usage.
 */
void importNodesets( integer const logLevel,
                     vtkDataSet & mesh,
                     string_array & nodesetNames,
                     CellBlockManager & cellBlockManager )
{
  auto & nodeSets = cellBlockManager.getNodeSets();
  localIndex const numPoints = LvArray::integerConversion< localIndex >( mesh.GetNumberOfPoints() );

  for( size_t i = 0; i < nodesetNames.size(); ++i )
  {
    string const & nodesetName = nodesetNames[i];

    GEOS_LOG_RANK_0_IF( logLevel >= 2, "    Processing nodeset: " + nodesetName );

    vtkAbstractArray * const curArray = mesh.GetPointData()->GetAbstractArray( nodesetName.c_str() );

    GEOS_THROW_IF( curArray == nullptr,
                   GEOS_FMT( "Nodeset '{}' not found in mesh point data", nodesetName ),
                   InputError );

    // Get the target nodeset container
    SortedArray< localIndex > & targetNodeset = nodeSets.get_inserted( nodesetName );

    // Get array metadata
    int const dataType = curArray->GetDataType();
    string const dataTypeName = curArray->GetDataTypeAsString();

    vtkDataArray * dataArray = vtkDataArray::FastDownCast( curArray );
    void const * rawData = dataArray->GetVoidPointer( 0 );

    // Dispatch based on data type
    switch( dataType )
    {
      case VTK_TYPE_UINT8:
        extractNodesetFromMask( static_cast< vtkTypeUInt8 const * >( rawData ), numPoints, targetNodeset );
        break;
      case VTK_TYPE_UINT16:
        extractNodesetFromMask( static_cast< vtkTypeUInt16 const * >( rawData ), numPoints, targetNodeset );
        break;
      case VTK_TYPE_UINT32:
        extractNodesetFromMask( static_cast< vtkTypeUInt32 const * >( rawData ), numPoints, targetNodeset );
        break;
      case VTK_TYPE_UINT64:
        extractNodesetFromMask( static_cast< vtkTypeUInt64 const * >( rawData ), numPoints, targetNodeset );
        break;
      default:
        GEOS_THROW( GEOS_FMT( "Nodeset '{}': unsupported type '{}' (use UInt8/16/32/64)", nodesetName, dataTypeName ), InputError );
    }
  }
}

void writeNodes( integer const logLevel,
                 vtkDataSet & mesh,
                 string_array & nodesetNames,
                 CellBlockManager & cellBlockManager,
                 const geos::R1Tensor & translate,
                 const geos::R1Tensor & scale )
{
  localIndex const numPts = LvArray::integerConversion< localIndex >( mesh.GetNumberOfPoints() );
  cellBlockManager.setNumNodes( numPts );

  // Writing the points
  arrayView1d< globalIndex > const nodeLocalToGlobal = cellBlockManager.getNodeLocalToGlobal();
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const X = cellBlockManager.getNodePositions();

  std::unordered_set< globalIndex > nodeGlobalIds;
  nodeGlobalIds.reserve( numPts );

  vtkIdTypeArray const & globalPointId = *vtkIdTypeArray::FastDownCast( mesh.GetPointData()->GetGlobalIds() );
  forAll< serialPolicy >( numPts, [&, X, nodeLocalToGlobal]( localIndex const k )
  {
    double point[3];
    mesh.GetPoint( k, point );
    LvArray::tensorOps::add< 3 >( point, translate );
    LvArray::tensorOps::hadamardProduct< 3 >( X[k], point, scale );
    globalIndex const pointGlobalID = globalPointId.GetValue( k );
    nodeLocalToGlobal[k] = pointGlobalID;

    // TODO: remove this check once the input mesh is cleaned of duplicate points via a filter
    //       and make launch policy parallel again
    GEOS_ERROR_IF( nodeGlobalIds.count( pointGlobalID ) > 0,
                   GEOS_FMT( "At least one duplicate point detected (globalID = {}).\n"
                             "Potential fixes :\n- Make sure partitionRefinement is set to 1 or higher.\n"
                             "- {}",
                             pointGlobalID, generalMeshErrorAdvice ) );
    nodeGlobalIds.insert( pointGlobalID );
  } );

  // Generate the "all" set
  array1d< localIndex > allNodes( numPts );
  std::iota( allNodes.begin(), allNodes.end(), 0 );
  SortedArray< localIndex > & allNodeSet = cellBlockManager.getNodeSets().get_inserted( "all" );
  allNodeSet.insert( allNodes.begin(), allNodes.end() );

  // Import remaining nodesets
  importNodesets( logLevel, mesh, nodesetNames, cellBlockManager );
}

void writeStructuredIndex( vtkDataSet & mesh,
                           string const & indexArrayName,
                           Span< vtkIdType const > const cellIds,
                           CellBlock & cellBlock )
{
  GEOS_ASSERT_EQ( LvArray::integerConversion< localIndex >( cellIds.size() ), cellBlock.size() );

  vtkDataArray * const srcDataArray = mesh.GetCellData()->GetArray( indexArrayName.c_str() );
  GEOS_THROW_IF( srcDataArray == nullptr, "Structured index array not found", InputError );
  integer const numComp = srcDataArray->GetNumberOfComponents();

  cellBlock.addProperty< fields::StructuredIndex::type >( fields::StructuredIndex::key() )
    .resizeDimension< 1 >( numComp );
  auto const dstIndex =
    cellBlock.getReference< fields::StructuredIndex::type >( fields::StructuredIndex::key() ).toView();

  vtkArrayDispatch::DispatchByValueType< vtkArrayDispatch::Integrals >::Execute( srcDataArray, [&]( auto const * const srcArray )
  {
    vtkDataArrayAccessor< TYPEOFPTR( srcArray ) > const srcIndex( srcArray );
    forAll< parallelHostPolicy >( cellBlock.size(), [numComp, dstIndex, srcIndex, cellIds]( localIndex const i )
    {
      for( integer c = 0; c < numComp; ++c )
      {
        dstIndex( i, c ) = static_cast< fields::StructuredIndex::dataType >( srcIndex.Get( cellIds[i], c ) );
      }
    } );
  } );
}

static void writePassthroughField( vtkDataSet & mesh,
                                    Span< vtkIdType const > const cellIds,
                                    string const & fieldName,
                                    CellBlock & cellBlock )
{
  vtkDataArray * srcArray = mesh.GetCellData()->GetArray( fieldName.c_str() );
  if( srcArray == nullptr )
  {
    GEOS_WARNING( GEOS_FMT( "Passthrough field '{}' not found in VTK CellData; skipping.", fieldName ) );
    return;
  }

  localIndex const numElems = LvArray::integerConversion< localIndex >( cellIds.size() );
  integer const numComp = srcArray->GetNumberOfComponents();

  // Map VTK float types to real64; everything else (integers, id types) to integer.
  bool const isRealArray = ( srcArray->GetDataType() == VTK_FLOAT || srcArray->GetDataType() == VTK_DOUBLE );

  if( isRealArray )
  {
    if( numComp == 1 )
    {
      array1d< real64 > & dst = cellBlock.addProperty< array1d< real64 > >( fieldName );
      dst.resize( numElems );
      cellBlock.getWrapperBase( fieldName ).setPlotLevel( dataRepository::PlotLevel::LEVEL_1 );
      forAll< parallelHostPolicy >( numElems, [&dst, srcArray, cellIds]( localIndex const i )
      {
        dst[i] = static_cast< real64 >( srcArray->GetTuple1( cellIds[i] ) );
      } );
    }
    else
    {
      array2d< real64 > & dst = cellBlock.addProperty< array2d< real64 > >( fieldName );
      dst.resize( numElems, numComp );
      cellBlock.getWrapperBase( fieldName ).setPlotLevel( dataRepository::PlotLevel::LEVEL_1 );
      forAll< parallelHostPolicy >( numElems, [&dst, srcArray, cellIds, numComp]( localIndex const i )
      {
        for( integer c = 0; c < numComp; ++c )
          dst( i, c ) = static_cast< real64 >( srcArray->GetComponent( cellIds[i], c ) );
      } );
    }
  }
  else
  {
    if( numComp == 1 )
    {
      array1d< integer > & dst = cellBlock.addProperty< array1d< integer > >( fieldName );
      dst.resize( numElems );
      cellBlock.getWrapperBase( fieldName ).setPlotLevel( dataRepository::PlotLevel::LEVEL_1 );
      forAll< parallelHostPolicy >( numElems, [&dst, srcArray, cellIds]( localIndex const i )
      {
        dst[i] = static_cast< integer >( srcArray->GetTuple1( cellIds[i] ) );
      } );
    }
    else
    {
      array2d< integer > & dst = cellBlock.addProperty< array2d< integer > >( fieldName );
      dst.resize( numElems, numComp );
      cellBlock.getWrapperBase( fieldName ).setPlotLevel( dataRepository::PlotLevel::LEVEL_1 );
      forAll< parallelHostPolicy >( numElems, [&dst, srcArray, cellIds, numComp]( localIndex const i )
      {
        for( integer c = 0; c < numComp; ++c )
          dst( i, c ) = static_cast< integer >( srcArray->GetComponent( cellIds[i], c ) );
      } );
    }
  }
}

void writeCells( integer const logLevel,
                 vtkDataSet & mesh,
                 vtk::CellMapType const & cellMap,
                 string const & structuredIndexAttributeName,
                 CellBlockManager & cellBlockManager,
                 Span< string const > passthroughFieldNames )
{
  if( !passthroughFieldNames.empty() )
  {
    string fieldList;
    for( string const & n : passthroughFieldNames )
      fieldList += "'" + n + "' ";
    GEOS_LOG_RANK_0( GEOS_FMT( "VTKMesh: propagating passthrough fields: {}", fieldList ) );
  }

  // Creates a new cell block for each region and for each type of cell.
  for( auto const & typeRegions : cellMap )
  {
    ElementType const elemType = typeRegions.first;
    if( getElementDim( elemType ) != 3 )
    {
      continue;
    }
    stdUnorderedMap< int, stdVector< vtkIdType > > const & regionIdToCellIds = typeRegions.second;
    for( auto const & regionCells : regionIdToCellIds )
    {
      int const regionId = regionCells.first;
      stdVector< vtkIdType > const & cellIds = regionCells.second;

      string const cellBlockName = vtk::buildCellBlockName( elemType, regionId );
      GEOS_LOG_RANK_0_IF( logLevel >= 1, "Importing cell block " << cellBlockName );

      // Create and resize the cell block.
      CellBlock & cellBlock = cellBlockManager.registerCellBlock( cellBlockName, regionId );
      cellBlock.setElementType( elemType );
      cellBlock.resize( LvArray::integerConversion< localIndex >( cellIds.size() ) );

      vtk::fillCellBlock( mesh, cellIds, cellBlock );

      if( !structuredIndexAttributeName.empty() )
      {
        writeStructuredIndex( mesh, structuredIndexAttributeName, cellIds, cellBlock );
      }

      for( string const & fieldName : passthroughFieldNames )
      {
        writePassthroughField( mesh, cellIds, fieldName, cellBlock );
      }
    }
  }
}

void writeSurfaces( integer const logLevel,
                    vtkDataSet & mesh,
                    const geos::vtk::CellMapType & cellMap,
                    CellBlockManager & cellBlockManager )
{
  if( cellMap.count( ElementType::Polygon ) == 0 )
  {
    return;
  }
  stdMap< string, SortedArray< localIndex > > & nodeSets = cellBlockManager.getNodeSets();

  for( auto const & surfaceCells: cellMap.at( ElementType::Polygon ) )
  {
    int const surfaceId = surfaceCells.first;
    stdVector< vtkIdType > const & cellIds = surfaceCells.second;
    string const surfaceName = std::to_string( surfaceId );
    GEOS_LOG_RANK_0_IF( logLevel >= 1, "Importing surface " << surfaceName );

    // Get or create all surfaces (even those which are empty in this rank)
    SortedArray< localIndex > & curNodeSet = nodeSets.get_inserted( surfaceName );

    for( vtkIdType const c : cellIds )
    {
      vtkCell * const currentCell = mesh.GetCell( c );
      for( int v = 0; v < currentCell->GetNumberOfPoints(); ++v )
      {
        curNodeSet.insert( currentCell->GetPointId( v ) );
      }
    }
  }
}

std::pair< real64, real64 > getGlobalLengthAndOffset( vtkDataSet & mesh )
{
  constexpr real64 minReal = LvArray::NumericLimits< real64 >::min;
  constexpr real64 maxReal = LvArray::NumericLimits< real64 >::max;
  real64 xMin[3] = { maxReal, maxReal, maxReal };
  real64 xMax[3] = { minReal, minReal, minReal };

  vtkBoundingBox bb( mesh.GetBounds() );
  if( bb.IsValid() )
  {
    bb.GetMinPoint( xMin );
    bb.GetMaxPoint( xMax );
  }

  MpiWrapper::min< real64 >( xMin, xMin, MPI_COMM_GEOS );
  MpiWrapper::max< real64 >( xMax, xMax, MPI_COMM_GEOS );
  LvArray::tensorOps::subtract< 3 >( xMax, xMin );

  // global length
  real64 size[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( xMax );
  LvArray::tensorOps::subtract< 3 >( size, xMin );
  // global offset
  real64 offset[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( xMin );
  LvArray::tensorOps::add< 3 >( offset, xMax );
  LvArray::tensorOps::scale< 3 >( offset, 0.5 );

  return { LvArray::tensorOps::l2Norm< 3 >( size ), LvArray::tensorOps::l2Norm< 3 >( offset ) };
}

} // namespace geos
