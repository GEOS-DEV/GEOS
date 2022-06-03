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

#include "common/DataTypes.hpp"
#include "common/DataLayouts.hpp"
#include "common/MpiWrapper.hpp"
#include "common/TypeDispatch.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshBody.hpp"
#include "mesh/generators/CellBlockManager.hpp"
#include "mesh/generators/VTKMeshGeneratorTools.hpp"
#include "mesh/generators/ParMETISInterface.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"

#include <vtkArrayDispatch.h>
#include <vtkBoundingBox.h>
#include <vtkCellData.h>
#include <vtkExtractCells.h>
#include <vtkGenerateGlobalIds.h>
#include <vtkPartitionedDataSet.h>
#include <vtkPointData.h>
#include <vtkRedistributeDataSetFilter.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLPUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>

#ifdef GEOSX_USE_MPI
#include <vtkMPIController.h>
#include <vtkMPI.h>
#else
#include <vtkDummyController.h>
#endif

#include <numeric>
#include <unordered_set>
#include <algorithm>
#include "mesh/CellElementSubRegion.hpp"
#include "common/TypeDispatch.hpp"
#include "mesh/MeshLevel.hpp"
#include "mesh/generators/CellBlockUtilities.hpp"

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

  registerWrapper( viewKeyStruct::partitionRefinementString(), &m_partitionRefinement ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1 ).
    setDescription( "Number of partitioning refinement iterations (defaults to 1, recommended value)."
                    "A value of 0 disables graph partitioning and keeps simple kd-tree partitions (not recommended). "
                    "Values higher than 1 may lead to slightly improved partitioning, but yield diminishing returns." );
}

namespace vtk
{

/**
 * @brief Return a VTK controller for multiprocessing.
 */
vtkSmartPointer< vtkMultiProcessController > getController()
{
#ifdef GEOSX_USE_MPI
  vtkNew< vtkMPIController > controller;
  vtkMPICommunicatorOpaqueComm vtkGeosxComm( &MPI_COMM_GEOSX );
  vtkNew< vtkMPICommunicator > communicator;
  communicator->InitializeExternal( &vtkGeosxComm );
  controller->SetCommunicator( communicator );
#else
  vtkNew< vtkDummyController > controller;
#endif
  return controller;
}

/**
 * @brief Load the VTK file into the VTK data structure
 * @param[in] filePath the Path of the file to load
 */
vtkSmartPointer< vtkUnstructuredGrid >
loadMesh( Path const & filePath )
{
  string const extension = filePath.extension();
  vtkSmartPointer< vtkUnstructuredGrid > loadedMesh;

  if( extension == "pvtu" )
  {
    auto const vtkUgReader = vtkSmartPointer< vtkXMLPUnstructuredGridReader >::New();
    vtkUgReader->SetFileName( filePath.c_str() );
    vtkUgReader->UpdateInformation();
    vtkUgReader->UpdatePiece( MpiWrapper::commRank(), MpiWrapper::commSize(), 0 );
    loadedMesh = vtkUgReader->GetOutput();

    // TODO: Apply vtkStaticCleanUnstructuredGrid, once it is included in the next VTK release.
    //       https://gitlab.kitware.com/vtk/vtk/-/blob/master/Filters/Core/vtkStaticCleanUnstructuredGrid.h
    //       This removes duplicate points, either present in the dataset, or resulting from merging pieces.
  }
  else
  {
    if( MpiWrapper::commRank() == 0 )
    {
      auto const read = [&]( auto const vtkUgReader )
      {
        vtkUgReader->SetFileName( filePath.c_str() );
        vtkUgReader->Update();
        return vtkUgReader->GetOutput();
      };

      if( extension == "vtk" )
      {
        loadedMesh = read( vtkSmartPointer< vtkUnstructuredGridReader >::New() );
      }
      else if( extension == "vtu" )
      {
        loadedMesh = read( vtkSmartPointer< vtkXMLUnstructuredGridReader >::New() );
      }
      else
      {
        GEOSX_ERROR( extension << " is not a recognized extension for VTKMesh. Please use .vtk, .vtu or .pvtu." );
      }
    }
    else
    {
      loadedMesh = vtkSmartPointer< vtkUnstructuredGrid >::New();
    }
  }

  return loadedMesh;
}

template< typename INDEX_TYPE, typename POLICY >
ArrayOfArrays< INDEX_TYPE, INDEX_TYPE >
buildElemToNodesImpl( vtkUnstructuredGrid & mesh )
{
  localIndex const numCells = LvArray::integerConversion< localIndex >( mesh.GetNumberOfCells() );
  array1d< INDEX_TYPE > nodeCounts( numCells );
  vtkCellArray & cells = *mesh.GetCells();

  // GetCellSize() is always thread-safe, can run in parallel
  forAll< parallelHostPolicy >( numCells, [nodeCounts = nodeCounts.toView(), &cells] ( localIndex const cellIdx )
  {
    nodeCounts[cellIdx] = LvArray::integerConversion< INDEX_TYPE >( cells.GetCellSize( cellIdx ) );
  } );

  ArrayOfArrays< INDEX_TYPE, INDEX_TYPE > elemToNodes;
  elemToNodes.template resizeFromCapacities< parallelHostPolicy >( numCells, nodeCounts.data() );

  vtkIdTypeArray const & globalPointId = *vtkIdTypeArray::FastDownCast( mesh.GetPointData()->GetGlobalIds() );

  // GetCellAtId() is conditionally thread-safe, use POLICY argument
  forAll< POLICY >( numCells, [&cells, &globalPointId, elemToNodes = elemToNodes.toView()] ( localIndex const cellIdx )
  {
    vtkIdType numPts;
    vtkIdType const * points;
    cells.GetCellAtId( cellIdx, numPts, points );
    for( int a = 0; a < numPts; ++a )
    {
      vtkIdType const pointIdx = globalPointId.GetValue( points[a] );
      elemToNodes.emplaceBack( cellIdx, LvArray::integerConversion< INDEX_TYPE >( pointIdx ) );
    }
  } );

  return elemToNodes;
}

template< typename INDEX_TYPE >
ArrayOfArrays< INDEX_TYPE, INDEX_TYPE >
buildElemToNodes( vtkUnstructuredGrid & mesh )
{
  // According to VTK docs, IsStorageShareable() indicates whether pointers extracted via
  // vtkCellArray::GetCellAtId() are pointers into internal storage rather than temp buffer
  // and thus results can be used in a thread-safe way.
  return mesh.GetCells()->IsStorageShareable()
       ? buildElemToNodesImpl< INDEX_TYPE, parallelHostPolicy >( mesh )
       : buildElemToNodesImpl< INDEX_TYPE, serialPolicy >( mesh );
}

template< typename PART_INDEX >
vtkSmartPointer< vtkPartitionedDataSet >
splitMeshByPartition( vtkUnstructuredGrid & mesh,
                      PART_INDEX const numParts,
                      arrayView1d< PART_INDEX const > const & part )
{
  array1d< localIndex > cellCounts( numParts );
  forAll< parallelHostPolicy >( part.size(), [part, cellCounts = cellCounts.toView()] ( localIndex const cellIdx )
  {
    RAJA::atomicInc< parallelHostAtomic >( &cellCounts[part[cellIdx]] );
  } );

  ArrayOfArrays< vtkIdType > cellsLists;
  cellsLists.resizeFromCapacities< serialPolicy >( numParts, cellCounts.data() );

  forAll< parallelHostPolicy >( part.size(), [part, cellsLists = cellsLists.toView()] ( localIndex const cellIdx )
  {
    cellsLists.emplaceBackAtomic< parallelHostAtomic >( LvArray::integerConversion< localIndex >( part[cellIdx] ),
                                                        LvArray::integerConversion< vtkIdType >( cellIdx ) );
  } );

  vtkNew< vtkPartitionedDataSet > result;
  result->SetNumberOfPartitions( LvArray::integerConversion< unsigned int >( numParts ) );

  vtkNew< vtkExtractCells > extractor;
  extractor->SetInputDataObject( &mesh );

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
  }
  return result;
}

vtkSmartPointer< vtkUnstructuredGrid >
redistributeByCellGraph( vtkUnstructuredGrid & mesh,
                         MPI_Comm const comm,
                         int const numRefinements )
{
  GEOSX_MARK_FUNCTION;

  // Use int64_t here to match ParMETIS' idx_t
  ArrayOfArrays< int64_t, int64_t > const elemToNodes = buildElemToNodes< int64_t >( mesh );
  int64_t const numProcs = MpiWrapper::commSize( comm );
  array1d< int64_t > const newParts = parmetis::partMeshKway( elemToNodes.toViewConst(), comm, 3, numRefinements );
  vtkSmartPointer< vtkPartitionedDataSet > const splitMesh = splitMeshByPartition( mesh, numProcs, newParts.toViewConst() );
  return vtk::redistribute( *splitMesh, MPI_COMM_GEOSX );
}

vtkSmartPointer< vtkUnstructuredGrid >
redistributeByKdTree( vtkUnstructuredGrid & mesh )
{
  GEOSX_MARK_FUNCTION;

  // Use a VTK filter which employs a kd-tree partition internally
  vtkNew< vtkRedistributeDataSetFilter > rdsf;
  rdsf->SetInputDataObject( &mesh );
  rdsf->SetNumberOfPartitions( MpiWrapper::commSize() );
  rdsf->Update();
  return vtkUnstructuredGrid::SafeDownCast( rdsf->GetOutputDataObject( 0 ) );
}

/**
 * @brief Compute the rank neighbor candidate list.
 * @param[in] boundingBoxes the bounding boxes used by the VTK partitioner for all ranks
 * @return the list of neighboring MPI ranks, will be updated
 */
std::vector< int >
findNeighborRanks( std::vector< vtkBoundingBox > boundingBoxes )
{
  int const numParts = LvArray::integerConversion< int >( boundingBoxes.size() );
  int const thisRank = MpiWrapper::commRank();

  // Inflate boxes to detect intersections more reliably
  double constexpr inflateFactor = 1.01;
  for( vtkBoundingBox & box : boundingBoxes )
  {
    box.ScaleAboutCenter( inflateFactor );
  }

  std::vector< int > neighbors;
  for( int i = 0; i < numParts; ++i )
  {
    if( i != thisRank && boundingBoxes[thisRank].Intersects( boundingBoxes[i] ) )
    {
      neighbors.push_back( i );
    }
  }

  return neighbors;
}

/**
 * @brief Generate global point/cell IDs and redistribute the mesh among MPI ranks.
 * @param[in] loadedMesh the mesh that was loaded on one or several MPI ranks
 * @param[in] comm the MPI communicator
 * @param[in] partitionRefinement number of graph partitioning refinement cycles
 */
vtkSmartPointer< vtkUnstructuredGrid >
redistributeMesh( vtkUnstructuredGrid & loadedMesh,
                  MPI_Comm const comm,
                  int const partitionRefinement )
{
  GEOSX_MARK_FUNCTION;

  // Generate global IDs for vertices and cells
  vtkNew< vtkGenerateGlobalIds > generator;
  generator->SetInputDataObject( &loadedMesh );
  generator->Update();
  vtkSmartPointer< vtkUnstructuredGrid > mesh =
    vtkUnstructuredGrid::SafeDownCast( generator->GetOutputDataObject( 0 ) );

  // Determine if redistribution is required
  vtkIdType const minCellsOnAnyRank = MpiWrapper::min( loadedMesh.GetNumberOfCells(), comm );
  if( minCellsOnAnyRank == 0 )
  {
    // Redistribute the mesh over all ranks using simple octree partitions
    mesh = redistributeByKdTree( *mesh );
  }

  // Redistribute the mesh again using higher-quality graph partitioner
  if( partitionRefinement > 0 )
  {
    mesh = redistributeByCellGraph( *mesh, comm, partitionRefinement - 1 );
  }

  return mesh;
}

/**
 * @brief Gathers all the data from all ranks, merge them, sort them, and remove duplicates.
 * @tparam T Type of the exchanged data.
 * @param data The data to be exchanged.
 * @return The merged data.
 * @note This function makes MPI calls.
 */
template< typename T >
std::vector< T > collectUniqueValues( std::vector< T > const & data )
{
  // Exchange the sizes of the data across all ranks.
  array1d< int > dataSizes( MpiWrapper::commSize() );
  MpiWrapper::allGather( LvArray::integerConversion< int >( data.size() ), dataSizes, MPI_COMM_GEOSX );
  // `totalDataSize` contains the total data size across all the MPI ranks.
  int const totalDataSize = std::accumulate( dataSizes.begin(), dataSizes.end(), 0 );

  // Once the MPI exchange is done, `allData` will contain all the data of all the MPI ranks.
  // We want all ranks to get all the data. But each rank may have a different size of information.
  // Therefore, we use `allgatherv` that does not impose the same size across ranks like `allgather` does.
  std::vector< T > allData( totalDataSize );
  // `displacements` is the offset (relative to the receive buffer) to store the data for each rank.
  std::vector< int > displacements( MpiWrapper::commSize(), 0 );
  std::partial_sum( dataSizes.begin(), dataSizes.end() - 1, displacements.begin() + 1 );
  MpiWrapper::allgatherv( data.data(), data.size(), allData.data(), dataSizes.data(), displacements.data(), MPI_COMM_GEOSX );

  // Finalizing by sorting, removing duplicates and trimming the result vector at the proper size.
  std::sort( allData.begin(), allData.end() );
  auto newEnd = std::unique( allData.begin(), allData.end() );
  allData.erase( newEnd, allData.end() );

  return allData;
}

/**
 * @brief Get the GEOSX element type
 * @param[in] cellType The vtk cell type
 * @return The GEOSX element type
 */
ElementType convertVtkToGeosxElementType( VTKCellType const cellType )
{
  switch( cellType )
  {
    case VTK_VERTEX:           return ElementType::Vertex;
    case VTK_LINE:             return ElementType::Line;
    case VTK_TRIANGLE:         return ElementType::Triangle;
    case VTK_QUAD:             return ElementType::Quadrilateral;
    case VTK_POLYGON:          return ElementType::Polygon;
    case VTK_TETRA:            return ElementType::Tetrahedron;
    case VTK_PYRAMID:          return ElementType::Pyramid;
    case VTK_WEDGE:            return ElementType::Wedge;
    case VTK_HEXAHEDRON:       return ElementType::Hexahedron;
    case VTK_PENTAGONAL_PRISM: return ElementType::Prism5;
    case VTK_HEXAGONAL_PRISM:  return ElementType::Prism6;
    case VTK_POLYHEDRON:       return ElementType::Polyhedron;
    default:
    {
      GEOSX_ERROR( cellType << " is not a recognized cell type to be used with the VTKMeshGenerator" );
    }
  }
  return ElementType::Polyhedron;
}

std::map< ElementType, std::vector< vtkIdType > >
splitCellsByType( vtkUnstructuredGrid & mesh )
{
  std::map< ElementType, std::vector< vtkIdType > > typeToCells;
  vtkIdType const numCells = mesh.GetNumberOfCells();

  // Count the number of each cell type
  std::array< size_t, numElementTypes() > cellTypeCounts{};
  for( vtkIdType c = 0; c < numCells; c++ )
  {
    ElementType const elemType = convertVtkToGeosxElementType( static_cast< VTKCellType >( mesh.GetCellType( c ) ) );
    ++cellTypeCounts[static_cast< integer >( elemType )];
  }

  // Allocate space to hold cell id lists by type
  std::array< std::vector< vtkIdType >, numElementTypes() > cellListsByType;
  for( integer t = 0; t < numElementTypes(); ++t )
  {
    cellListsByType[t].reserve( cellTypeCounts[t] );
  }

  // Collect cell lists for each type (using array for speed in a hot loop)
  for( vtkIdType c = 0; c < numCells; c++ )
  {
    ElementType const elemType = convertVtkToGeosxElementType( static_cast< VTKCellType >( mesh.GetCellType( c ) ) );
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
        std::vector< vtkIdType > & surfaceCells = typeToCells[ ElementType::Polygon ];
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
        GEOSX_ERROR( "Invalid element dimension: " << getElementDim( type ) );
      }
    }
  }

  return typeToCells;
}

VTKMeshGenerator::CellMapType
splitCellsByTypeAndAttribute( std::map< ElementType, std::vector< vtkIdType > > & typeToCells,
                              vtkDataArray * const attributeDataArray )
{
  VTKMeshGenerator::CellMapType typeToAttributeToCells;
  for( auto & t2c : typeToCells )
  {
    ElementType const elemType = t2c.first;
    std::vector< vtkIdType > & cells = t2c.second;
    std::unordered_map< int, std::vector< vtkIdType > > & attributeToCells = typeToAttributeToCells[elemType];

    if( attributeDataArray == nullptr )
    {
      attributeToCells.emplace( -1, std::move( cells ) );
    }
    else
    {
      GEOSX_ERROR_IF_NE_MSG( attributeDataArray->GetNumberOfComponents(), 1,
                             "Invalid number of components in attribute array" );
      vtkArrayDispatch::Dispatch::Execute( attributeDataArray, [&]( auto const * attributeArray )
      {
        using ArrayType = TYPEOFPTR( attributeArray );
        vtkDataArrayAccessor< ArrayType > attribute( attributeArray );
        std::unordered_map< int, size_t > cellCounts;
        for( vtkIdType c: cells )
        {
          int const region = static_cast< int >( attribute.Get( c, 0 ) );
          ++cellCounts[region];
        }
        for( auto const & count : cellCounts )
        {
          attributeToCells[count.first].reserve( count.second );
        }
        for( vtkIdType c: cells )
        {
          int const region = static_cast< int >( attribute.Get( c, 0 ) );
          attributeToCells[region].push_back( c );
        }
      } );
    }
  }
  return typeToAttributeToCells;
}

void extendCellMapWithRemoteKeys( VTKMeshGenerator::CellMapType & cellMap )
{
  // Gather all element types encountered on any rank and enrich the local collection
  std::vector< ElementType > allElementTypes = collectUniqueValues( mapKeys( cellMap ) );
  std::vector< int > allCellAttributes;
  for( auto const & typeRegions : cellMap )
  {
    if( getElementDim( typeRegions.first ) == 3 )
    {
      std::vector< int > const attrs = mapKeys( typeRegions.second );
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
        cellMap[elemType][attrValue];
      }
    }
  }

  // Treat surfaces separately - and avoid inadvertently creating a map entry for polygons
  std::vector< int > const surfaceAttributes = cellMap.count( ElementType::Polygon ) > 0
                                             ? mapKeys( cellMap.at( ElementType::Polygon ) )
                                             : std::vector< int >();
  std::vector< int > allSurfaceAttributes = collectUniqueValues( surfaceAttributes );
  for( int attrValue: allSurfaceAttributes )
  {
    cellMap[ElementType::Polygon][attrValue];
  }
}

///////////////////////////////////////////////////////////////////////////////
int PolyhedronType( vtkCell* const cell )
{
int numberOfFaces = cell->GetNumberOfFaces();
int numberOfQuads = 0;
int numberOfTriangles = 0;
set < int > ListQuadsPoints;
set < int > ListNoQuadsPoints;

if( cell->GetCellType() != 42 ) return 1;     // Stop with error
//if( cell->GetCellType() != 42 ) GEOSX_ERROR( "Input for PolyhedronType() must be a polyhedron." );

// Compute the number of triangles
for( int iFace = 0 ; iFace < numberOfFaces ; ++iFace )
 if( cell->GetFace(iFace)->GetNumberOfPoints() == 3 ) numberOfTriangles++;

// Compute the number of quads
for( int iFace = 0 ; iFace < numberOfFaces ; ++iFace )
 if( cell->GetFace(iFace)->GetNumberOfPoints() == 4 ) numberOfQuads++;

// VTK_TETRA (=10) 
if(numberOfTriangles == 4 && numberOfFaces == 4) return 10;

// VTK_HEXAHEDRON (=12)
if(numberOfQuads == 6 && numberOfFaces == 6) return 12;

// VTK_WEDGE (=13)
if(numberOfTriangles == 2 && numberOfQuads == 3 && numberOfFaces == 5) return 13;

// VTK_PYRAMID (=14)
if(numberOfTriangles == 4 && numberOfQuads == 1 && numberOfFaces == 5) return 14;

// Check if the polyhedron is a prism
if( numberOfFaces - numberOfQuads == 2)
 {
  for( int iFace = 0 ; iFace < numberOfFaces ; ++iFace )
   if( cell->GetFace(iFace)->GetNumberOfPoints() == 4)
    for( int iPoint = 0 ; iPoint < 4 ; ++iPoint )
     ListQuadsPoints.insert( cell->GetFace(iFace)->GetPointId(iPoint) );
   else
    for( int iPoint = 0 ; iPoint < cell->GetFace(iFace)->GetNumberOfPoints() ; ++iPoint )
     ListNoQuadsPoints.insert( cell->GetFace(iFace)->GetPointId(iPoint) );
 }
else return 0;     // General polyhedron

if(ListQuadsPoints == ListNoQuadsPoints)   // The polyhedron is a prism
 {
  if( numberOfQuads == 5) return 15;   // VTK_PENTAGONAL_PRISM (=15)
  if( numberOfQuads == 6) return 16;   // VTK_HEXAGONAL_PRISM (=16)
  return (-numberOfQuads);             // N_PRISM
 }
else return 0;     // General polyhedron
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
std::vector< int > GetTetraNodeOrderingFromPolyhedron( vtkCell* const cell )
{
std::vector < int > nodeOrder(4);

// Check if the orientation is correct
real64 vectorA[3];
real64 vectorB[3];
real64 vectorC[3];
real64 vectorN[3];

for( int i = 0; i < 3; ++i )        // HOW TO AVOID FOR LOOP?!?!?
 {
   vectorA[i] = cell->GetPoints()->GetPoint(1)[i] - cell->GetPoints()->GetPoint(0)[i];
   vectorB[i] = cell->GetPoints()->GetPoint(2)[i] - cell->GetPoints()->GetPoint(0)[i];
   vectorC[i] = cell->GetPoints()->GetPoint(3)[i] - cell->GetPoints()->GetPoint(0)[i];
 }

LvArray::tensorOps::crossProduct( vectorN, vectorA, vectorB );

if( LvArray::tensorOps::AiBi< 3 >( vectorN, vectorC ) > 0 ) // The orientation is correctly
 {
  nodeOrder[0] = 0;
  nodeOrder[1] = 1;
  nodeOrder[2] = 2;
  nodeOrder[3] = 3;
 }
else // The orientation is incorrectly
 {
  nodeOrder[0] = 0;
  nodeOrder[1] = 2;
  nodeOrder[2] = 1;
  nodeOrder[3] = 3;
 }

return nodeOrder;
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
std::vector< int > GetHexahedronNodeOrderingFromPolyhedron( vtkCell* const cell )
{
int iFace;
std::vector < int > nodeOrder(8);
std::map < int , int > G2L;

// Generate global to local map
for( int iPoint = 0 ; iPoint < 8 ; ++iPoint ) G2L.insert({cell->GetPointId(iPoint),iPoint});

// Assuming the input parameters are correct, take the first quad
iFace = 0;

// Get global pointIds for the first quad
nodeOrder[0] = cell->GetFace(iFace)->GetPointId(0);
nodeOrder[1] = cell->GetFace(iFace)->GetPointId(1);
nodeOrder[2] = cell->GetFace(iFace)->GetPointId(3);
nodeOrder[3] = cell->GetFace(iFace)->GetPointId(2);

// Using the edges, generate the global pointIds for the opposit quad
for( int iEdge = 0 ; iEdge < 12 ; ++iEdge )
 {
  auto edgeNode0 = cell->GetEdge(iEdge)->GetPointId(0);
  auto edgeNode1 = cell->GetEdge(iEdge)->GetPointId(1);
  auto it0 = std::find( &nodeOrder[0] , &nodeOrder[4] , edgeNode0 );
  auto it1 = std::find( &nodeOrder[0] , &nodeOrder[4] , edgeNode1 );

  if(it0 != &nodeOrder[4] && it1 == &nodeOrder[4])
   nodeOrder[ 4 + std::distance(&nodeOrder[0],it0) ] = edgeNode1;

  if(it0 == &nodeOrder[4] && it1 != &nodeOrder[4])
   nodeOrder[ 4 + std::distance(&nodeOrder[0],it1) ] = edgeNode0;
 }

// Convert global numbering to local numbering
for ( int iPoint = 0; iPoint < 8; ++iPoint ) nodeOrder[iPoint] = G2L.find(nodeOrder[iPoint])->second;

// Check if the quads are correctly positioned
real64 vectorA[3];
real64 vectorB[3];
real64 vectorC[3];
real64 vectorN[3];

for( int i = 0; i < 3; ++i )        // HOW TO AVOID FOR LOOP?!?!?
 {
   vectorA[i] = cell->GetPoints()->GetPoint(nodeOrder[1])[i] - cell->GetPoints()->GetPoint(nodeOrder[0])[i];
   vectorB[i] = cell->GetPoints()->GetPoint(nodeOrder[2])[i] - cell->GetPoints()->GetPoint(nodeOrder[0])[i];
   vectorC[i] = cell->GetPoints()->GetPoint(nodeOrder[4])[i] - cell->GetPoints()->GetPoint(nodeOrder[0])[i];
 }

LvArray::tensorOps::crossProduct( vectorN, vectorA, vectorB );

// If incorrect swap the quads
if( LvArray::tensorOps::AiBi< 3 >( vectorN, vectorC ) < 0 )
 std::rotate( nodeOrder.begin(), nodeOrder.begin()+4, nodeOrder.end());

return nodeOrder;
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
std::vector< int > GetWedgeNodeOrderingFromPolyhedron( vtkCell* const cell )
{
int iFace;
std::vector < int > nodeTri0(3);
std::vector < int > nodeTri1(3);
std::vector < int > nodeOrder(6);
std::map < int , int > G2L;

// Generate global to local map
for( int iPoint = 0 ; iPoint < 6 ; ++iPoint ) G2L.insert({cell->GetPointId(iPoint),iPoint});

// Assuming the input parameters are correct, identify one of the triangles
iFace = 0; while( cell->GetFace(iFace)->GetNumberOfPoints() != 3 ) iFace++;

// Get global pointIds for the first triangle
nodeTri0[0] = cell->GetFace(iFace)->GetPointId(0);
nodeTri0[1] = cell->GetFace(iFace)->GetPointId(1);
nodeTri0[2] = cell->GetFace(iFace)->GetPointId(2);

// Using the edges, generate the global pointIds for the second triangle
for( int iEdge = 0 ; iEdge < 9 ; ++iEdge )
{
  auto edgeNode0 = cell->GetEdge(iEdge)->GetPointId(0);
  auto edgeNode1 = cell->GetEdge(iEdge)->GetPointId(1);
  auto it0 = std::find( &nodeTri0[0] , &nodeTri0[3] , edgeNode0 );
  auto it1 = std::find( &nodeTri0[0] , &nodeTri0[3] , edgeNode1 );

  if( it0 != &nodeTri0[3] && it1 == &nodeTri0[3] )
   nodeTri1[std::distance(&nodeTri0[0],it0)] = edgeNode1;

  if( it0 == &nodeTri0[3] && it1 != &nodeTri0[3] )
   nodeTri1[std::distance(&nodeTri0[0],it1)] = edgeNode0;
}

// Convert global numbering to local numbering
for ( int iPoint = 0; iPoint < 3; ++iPoint ) 
 {
  nodeTri0[iPoint] = G2L.find(nodeTri0[iPoint])->second;
  nodeTri1[iPoint] = G2L.find(nodeTri1[iPoint])->second;
 }

// Check if the bases are correctly positioned
real64 vectorA[3];
real64 vectorB[3];
real64 vectorC[3];
real64 vectorN[3];

for( int i = 0; i < 3; ++i )        // HOW TO AVOID FOR LOOP?!?!?
 {
   vectorA[i] = cell->GetPoints()->GetPoint(nodeTri0[1])[i] - cell->GetPoints()->GetPoint(nodeTri0[0])[i];
   vectorB[i] = cell->GetPoints()->GetPoint(nodeTri0[2])[i] - cell->GetPoints()->GetPoint(nodeTri0[0])[i];
   vectorC[i] = cell->GetPoints()->GetPoint(nodeTri1[0])[i] - cell->GetPoints()->GetPoint(nodeTri0[0])[i];
 }

LvArray::tensorOps::crossProduct( vectorN, vectorA, vectorB );

if( LvArray::tensorOps::AiBi< 3 >( vectorN, vectorC ) > 0 ) // Bases are correctly positioned
 {
  nodeOrder[0]=nodeTri0[0];
  nodeOrder[1]=nodeTri1[0];
  nodeOrder[2]=nodeTri0[1];
  nodeOrder[3]=nodeTri1[1];
  nodeOrder[4]=nodeTri0[2];
  nodeOrder[5]=nodeTri1[2];
 }
else // Bases are incorrectly positioned
 {
  nodeOrder[0]=nodeTri0[0];
  nodeOrder[1]=nodeTri1[0];
  nodeOrder[2]=nodeTri0[2];
  nodeOrder[3]=nodeTri1[2];
  nodeOrder[4]=nodeTri0[1];
  nodeOrder[5]=nodeTri1[1];   
 }

return nodeOrder;
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
std::vector< int > GetPyramidNodeOrderingFromPolyhedron( vtkCell* const cell )
{
int iPoint;
int iFace;
std::vector < int > nodeOrder(5);
std::map < int , int > G2L;

// Generate global to local map
for(iPoint = 0 ; iPoint < 5 ; ++iPoint ) G2L.insert({cell->GetPointId(iPoint),iPoint});

// Assuming the input parameters are correct, identify the base
iFace = 0; while( cell->GetFace(iFace)->GetNumberOfPoints() != 4 ) iFace++;

// Get global pointIds for the base
nodeOrder[0] = cell->GetFace(iFace)->GetPointId(0);
nodeOrder[1] = cell->GetFace(iFace)->GetPointId(1);
nodeOrder[2] = cell->GetFace(iFace)->GetPointId(3);
nodeOrder[3] = cell->GetFace(iFace)->GetPointId(2);

// Identify the missing node
iPoint = 0; while(std::find( &nodeOrder[0] , &nodeOrder[4] , cell->GetPointId(iPoint) ) != &nodeOrder[4]) iPoint++;

// Add it to "nodeOrder"
nodeOrder[4] = cell->GetPointId(iPoint);

// Convert global numbering to local numbering
for (iPoint = 0; iPoint < 5; ++iPoint ) nodeOrder[iPoint] = G2L.find(nodeOrder[iPoint])->second;

// Check the orientation of the pyramid
real64 vectorA[3];
real64 vectorB[3];
real64 vectorC[3];
real64 vectorN[3];

for( int i = 0; i < 3; ++i )        // HOW TO AVOID FOR LOOP?!?!?
 {
   vectorA[i] = cell->GetPoints()->GetPoint(nodeOrder[1])[i] - cell->GetPoints()->GetPoint(nodeOrder[0])[i];
   vectorB[i] = cell->GetPoints()->GetPoint(nodeOrder[2])[i] - cell->GetPoints()->GetPoint(nodeOrder[0])[i];
   vectorC[i] = cell->GetPoints()->GetPoint(nodeOrder[4])[i] - cell->GetPoints()->GetPoint(nodeOrder[0])[i];
 }

LvArray::tensorOps::crossProduct( vectorN, vectorA, vectorB );

// If incorrect swap two nodes from the base
if( LvArray::tensorOps::AiBi< 3 >( vectorN, vectorC ) < 0 ) std::swap(nodeOrder[1],nodeOrder[2]);

return nodeOrder;
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
std::vector< int > GetPrismNodeOrderingFromPolyhedron( vtkCell* const cell , int numberOfSides )
{
int iFace;
std::vector < int > nodeOrder(2*numberOfSides);
std::map < int , int > G2L;

// Generate global to local map
for( int iPoint = 0 ; iPoint < cell->GetNumberOfPoints() ; ++iPoint ) G2L.insert({cell->GetPointId(iPoint),iPoint});

// Assuming the input parameters are correct, identify one of the bases
iFace = 0; while( cell->GetFace(iFace)->GetNumberOfPoints() != numberOfSides ) iFace++;

// Get global pointIds for the first base
for( int iPoint = 0 ; iPoint < numberOfSides ; ++iPoint ) nodeOrder[iPoint] = cell->GetFace(iFace)->GetPointId(iPoint);

// Using the edges, generate the global pointIds for the second base
for( int iEdge = 0 ; iEdge < cell->GetNumberOfEdges() ; ++iEdge )
{
  auto edgeNode0 = cell->GetEdge(iEdge)->GetPointId(0);
  auto edgeNode1 = cell->GetEdge(iEdge)->GetPointId(1);
  auto it0 = std::find( &nodeOrder[0] , &nodeOrder[numberOfSides] , edgeNode0 );
  auto it1 = std::find( &nodeOrder[0] , &nodeOrder[numberOfSides] , edgeNode1 );

  if( it0 != &nodeOrder[numberOfSides] && it1 == &nodeOrder[numberOfSides] )
   nodeOrder[numberOfSides + std::distance(&nodeOrder[0],it0)] = edgeNode1;

  if( it0 == &nodeOrder[numberOfSides] && it1 != &nodeOrder[numberOfSides] )
   nodeOrder[numberOfSides + std::distance(&nodeOrder[0],it1)] = edgeNode0;
}

// Convert global numbering to local numbering
for ( int iPoint = 0; iPoint < 2*numberOfSides; ++iPoint ) nodeOrder[iPoint] = G2L.find(nodeOrder[iPoint])->second;

// Check if the bases are correctly positioned
real64 vectorA[3];
real64 vectorB[3];
real64 vectorC[3];
real64 vectorN[3];

for( int i = 0; i < 3; ++i )        // HOW TO AVOID FOR LOOP?!?!?
 {
   vectorA[i] = cell->GetPoints()->GetPoint(nodeOrder[1])[i] - cell->GetPoints()->GetPoint(nodeOrder[0])[i];
   vectorB[i] = cell->GetPoints()->GetPoint(nodeOrder[(numberOfSides+1)/2])[i] - cell->GetPoints()->GetPoint(nodeOrder[0])[i];
   vectorC[i] = cell->GetPoints()->GetPoint(nodeOrder[numberOfSides])[i] - cell->GetPoints()->GetPoint(nodeOrder[0])[i];
 }

LvArray::tensorOps::crossProduct( vectorN, vectorA, vectorB );

// If incorrect swap the bases
if( LvArray::tensorOps::AiBi< 3 >( vectorN, vectorC ) < 0 )
 std::rotate( nodeOrder.begin(), nodeOrder.begin()+numberOfSides, nodeOrder.end());

return nodeOrder;
}
///////////////////////////////////////////////////////////////////////////////

/**
 * @brief Collect lists of VTK cell indices organized by type and attribute value.
 * @param[in] mesh the vtkUnstructuredGrid that is loaded
 * @param[in] attributeName name of the VTK data array containing the attribute, if any
 * @return A map from element type to a map of attribute to the associated cell ids for the current rank.
 *         The map contains entries for all types and attribute values across all MPI ranks,
 *         even if there are no cells on current rank (then the list will be empty).
 */
VTKMeshGenerator::CellMapType
buildCellMap( vtkUnstructuredGrid & mesh,
              string const & attributeName )
{

std::vector < int >  nodeOrder;
int polyhedronType;

for ( int i = 0 ; i < mesh.GetNumberOfCells() ; ++i )
 {
  cout << "Cell(" << i << ") = " << mesh.GetCell(i)->GetCellType() << "\n";

  if (mesh.GetCell(i)->GetCellType() == 42) polyhedronType = PolyhedronType( mesh.GetCell(i) );
  else    polyhedronType = 0 ;

 cout << "polyhedronType = " << polyhedronType << "\n";

  if( polyhedronType == 10 )    // TETRA
   {
    nodeOrder = GetTetraNodeOrderingFromPolyhedron( mesh.GetCell(i) );

    cout << "Bases points using local numbering:" << "\n"; 
    for (const auto item : nodeOrder) cout << item << ", ";
    cout << "\n\n";

  real64 Xlocal[12][3];
  vtkCell *cell = mesh.GetCell(i);
  real64 volume;

  for (int j = 0 ; j < cell->GetNumberOfPoints() ; ++j)
   {
     Xlocal[j][0]=cell->GetPoints()->GetPoint(j)[0];
     Xlocal[j][1]=cell->GetPoints()->GetPoint(j)[1];
     Xlocal[j][2]=cell->GetPoints()->GetPoint(j)[2];
   }

   for (int j = 0 ; j < cell->GetNumberOfPoints() ; ++j)
    {
    int k = nodeOrder[j];
    cout << Xlocal[k][0] << " " << Xlocal[k][1] << " " << Xlocal[k][2] << " " << k << " " << j << "\n";
    }    
   
   } 


  if( polyhedronType == 13 )    // WEDGE
   {
    nodeOrder = GetWedgeNodeOrderingFromPolyhedron( mesh.GetCell(i) );

    cout << "Bases points using local numbering:" << "\n"; 
    for (const auto item : nodeOrder) cout << item << ", ";
    cout << "\n\n";

  real64 Xlocal[12][3];
  vtkCell *cell = mesh.GetCell(i);
  real64 volume;

  for (int j = 0 ; j < cell->GetNumberOfPoints() ; ++j)
   {
     Xlocal[j][0]=cell->GetPoints()->GetPoint(j)[0];
     Xlocal[j][1]=cell->GetPoints()->GetPoint(j)[1];
     Xlocal[j][2]=cell->GetPoints()->GetPoint(j)[2];
   }

   for (int j = 0 ; j < cell->GetNumberOfPoints() ; ++j)
    {
    int k = nodeOrder[j];
    cout << Xlocal[k][0] << " " << Xlocal[k][1] << " " << Xlocal[k][2] << " " << k << " " << j << "\n";
    }    
   
   } 

  if( polyhedronType == 14 )    // PYRAMID
   {
    nodeOrder = GetPyramidNodeOrderingFromPolyhedron( mesh.GetCell(i) );

    cout << "Bases points using local numbering:" << "\n"; 
    for (const auto item : nodeOrder) cout << item << ", ";
    cout << "\n\n";

  real64 Xlocal[12][3];
  vtkCell *cell = mesh.GetCell(i);
  real64 volume;

  for (int j = 0 ; j < cell->GetNumberOfPoints() ; ++j)
   {
     Xlocal[j][0]=cell->GetPoints()->GetPoint(j)[0];
     Xlocal[j][1]=cell->GetPoints()->GetPoint(j)[1];
     Xlocal[j][2]=cell->GetPoints()->GetPoint(j)[2];
   }

   for (int j = 0 ; j < cell->GetNumberOfPoints() ; ++j)
    {
    int k = nodeOrder[j];
    cout << Xlocal[k][0] << " " << Xlocal[k][1] << " " << Xlocal[k][2] << " " << k << " " << j << "\n";
    }    
   
   } 


  if( polyhedronType == 12 )     // HEXA
   {
    nodeOrder = GetHexahedronNodeOrderingFromPolyhedron( mesh.GetCell(i) );

    cout << "Bases points using local numbering:" << "\n"; 
    for (const auto item : nodeOrder) cout << item << ", ";
    cout << "\n\n";

  real64 Xlocal[12][3];
  vtkCell *cell = mesh.GetCell(i);
  real64 volume;

  for (int j = 0 ; j < cell->GetNumberOfPoints() ; ++j)
   {
     Xlocal[j][0]=cell->GetPoints()->GetPoint(j)[0];
     Xlocal[j][1]=cell->GetPoints()->GetPoint(j)[1];
     Xlocal[j][2]=cell->GetPoints()->GetPoint(j)[2];
   }

   for (int j = 0 ; j < cell->GetNumberOfPoints() ; ++j)
    {
    int k = nodeOrder[j];
    cout << Xlocal[k][0] << " " << Xlocal[k][1] << " " << Xlocal[k][2] << " " << k << " " << j << "\n";
    }    


   }
  if( polyhedronType < 0 ) 
  {
   nodeOrder = GetPrismNodeOrderingFromPolyhedron(mesh.GetCell(i), -polyhedronType );

   cout << "Bases points using local numbering:" << "\n"; 
   for (const auto item : nodeOrder) cout << item << ", ";
   cout << "\n\n";

  real64 Xlocal[12][3];
  vtkCell *cell = mesh.GetCell(i);
  real64 volume;

  for (int j = 0 ; j < cell->GetNumberOfPoints() ; ++j)
   {
     Xlocal[j][0]=cell->GetPoints()->GetPoint(j)[0];
     Xlocal[j][1]=cell->GetPoints()->GetPoint(j)[1];
     Xlocal[j][2]=cell->GetPoints()->GetPoint(j)[2];
   }

   for (int j = 0 ; j < cell->GetNumberOfPoints() ; ++j)
    {
    int k = nodeOrder[j];
    cout << Xlocal[k][0] << " " << Xlocal[k][1] << " " << Xlocal[k][2] << " " << k << " " << j << "\n";
    }
   //volume = geosx::computationalGeometry::prismVolume< 5 >( Xlocal );
   //cout << "Volume = " << volume << "\n";

   }

  if( polyhedronType == 15 || polyhedronType == 16 ) 
  //if( polyhedronType == 15 ) 
  {
    cout << "Prism56" << "\n";
   nodeOrder = GetPrismNodeOrderingFromPolyhedron(mesh.GetCell(i), polyhedronType-10 );

   cout << "Bases points using local numbering:" << "\n"; 
   for (const auto item : nodeOrder) cout << item << ", ";
   cout << "\n\n";

  real64 Xlocal[12][3];
  vtkCell *cell = mesh.GetCell(i);
  real64 volume;

  for (int j = 0 ; j < cell->GetNumberOfPoints() ; ++j)
   {
     Xlocal[j][0]=cell->GetPoints()->GetPoint(j)[0];
     Xlocal[j][1]=cell->GetPoints()->GetPoint(j)[1];
     Xlocal[j][2]=cell->GetPoints()->GetPoint(j)[2];
   }

   for (int j = 0 ; j < cell->GetNumberOfPoints() ; ++j)
   {
    int k = nodeOrder[j];
    cout << Xlocal[k][0] << " " << Xlocal[k][1] << " " << Xlocal[k][2] << " " << k << " " << j << "\n";
   }
   
   //volume = geosx::real64 geosx::computationalGeometry::prismVolume< 5 >( Xlocal );
   //cout << "Volume = " << volume << "\n";
   }
 }
exit(0);

  // First, pass through all VTK cells and split them int sub-lists based on type.
  std::map< ElementType, std::vector< vtkIdType > > typeToCells = splitCellsByType( mesh );

  // Now, actually split into groups according to region attribute, if present
  vtkDataArray * const attributeDataArray =
    vtkDataArray::FastDownCast( mesh.GetCellData()->GetAbstractArray( attributeName.c_str() ) );

  VTKMeshGenerator::CellMapType cellMap =
    splitCellsByTypeAndAttribute( typeToCells, attributeDataArray );

  // Gather all element types encountered on any rank and enrich the local collection
  extendCellMapWithRemoteKeys( cellMap );

  return cellMap;
}

std::vector< int > getGeosxToVtkNodeOrdering( ElementType const elemType )
{
  switch( elemType )
  {
    case ElementType::Vertex:        return { 0 };
    case ElementType::Line:          return { 0, 1 };
    case ElementType::Triangle:      return { 0, 1, 2 };
    case ElementType::Quadrilateral: return { 0, 1, 2, 3 }; // TODO check
    case ElementType::Polygon:       return { 0, 1, 2, 3, 4, 5, 6, 7, 8 }; // TODO
    case ElementType::Tetrahedron:   return { 0, 1, 2, 3 };
    case ElementType::Pyramid:       return { 0, 1, 3, 2, 4 };
    case ElementType::Wedge:         return { 0, 3, 2, 5, 1, 4 };
    case ElementType::Hexahedron:    return { 0, 1, 3, 2, 4, 5, 7, 6 };
    case ElementType::Prism5:        return { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    case ElementType::Prism6:        return { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 };
    case ElementType::Polyhedron:    return { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 }; // TODO
  }
  return {};
}

/**
 * @brief Fill @p cellBlock with the appropriate nodes and local/global mappings.
 * @param[in] elemType the vtk cell type for cells of the CellBlock being written
 * @param[in] cellIds the cell indexes of cell type \p cellType within this region
 * @param[in] mesh the vtkUnstructuredGrid that is loaded
 * @param[in,out] cellBlock The cell block to be written
 */
void fillCellBlock( vtkUnstructuredGrid & mesh,
                    ElementType const elemType,
                    std::vector< vtkIdType > const & cellIds,
                    CellBlock & cellBlock )
{
  localIndex const numNodesPerElement = cellBlock.numNodesPerElement();
  arrayView2d< localIndex, cells::NODE_MAP_USD > const cellToVertex = cellBlock.getElemToNode();
  arrayView1d< globalIndex > const & localToGlobal = cellBlock.localToGlobalMap();

  std::vector< int > const nodeOrder = getGeosxToVtkNodeOrdering( elemType );
  vtkIdTypeArray const * const globalCellId = vtkIdTypeArray::FastDownCast( mesh.GetCellData()->GetGlobalIds() );
  GEOSX_ERROR_IF( !cellIds.empty() && globalCellId == nullptr, "Global cell IDs have not been generated" );

  // Writing connectivity and Local to Global
  localIndex cellCount = 0;
  for( vtkIdType c: cellIds )
  {
    vtkCell * const currentCell = mesh.GetCell( c );
    for( localIndex v = 0; v < numNodesPerElement; v++ )
    {
      cellToVertex[cellCount][v] = currentCell->GetPointId( nodeOrder[v] );
    }
    localToGlobal[cellCount++] = globalCellId->GetValue( c );
  }
}

/**
 * @brief Returns a string describing the element.
 * @param[in] type The element type.
 * @return The name.
 * @warning This information will be visible in the input file... Consider refactoring with great care.
 */
string getElementTypeName( ElementType const type )
{
  switch( type )
  {
    case ElementType::Hexahedron:  return "hexahedra";
    case ElementType::Tetrahedron: return "tetrahedra";
    case ElementType::Wedge:       return "wedges";
    case ElementType::Pyramid:     return "pyramids";
    case ElementType::Prism5:      return "pentagonalPrisms";
    case ElementType::Prism6:      return "hexagonalPrisms";
    case ElementType::Polyhedron:  return "polyhedra";
    default:
    {
      GEOSX_ERROR( "Element type '" << type << "' is not supported" );
      return {};
    }
  }
}

/**
 * @brief Builds the cell block name.
 * @param[in] elementName The element name.
 * @param[in] regionId The region Id.
 * @return The name.
 * @warning This name will be visible in the input file... Consider refactoring with great care.
 */
string buildCellBlockName( ElementType const type, int const regionId )
{
  GEOSX_ERROR_IF_LT_MSG( regionId, -1, "Invalid region id" );
  string const cellTypeName = getElementTypeName( type );
  return regionId != -1 ? std::to_string( regionId ) + "_" + cellTypeName : cellTypeName;
}

/**
 * @brief Collect a set of material field names registered in a subregion.
 * @param subRegion the target subregion
 * @return a set of wrapper names
 */
std::unordered_set< string > getMaterialWrapperNames( ElementSubRegionBase const & subRegion )
{
  using namespace constitutive;
  std::unordered_set< string > materialWrapperNames;
  subRegion.getConstitutiveModels().forSubGroups< ConstitutiveBase >( [&]( ConstitutiveBase const & material )
  {
    material.forWrappers( [&]( WrapperBase const & wrapper )
    {
      if( wrapper.sizedFromParent() )
      {
        materialWrapperNames.insert( ConstitutiveBase::makeFieldName( material.getName(), wrapper.getName() ) );
      }
    } );
  } );
  return materialWrapperNames;
}

/**
 * @brief Imports 2d and 3d arrays from @p vtkArray to @p wrapper, only for @p cellIds
 * @param cellIds The cells for which we should copy the data.
 * @param vtkArray The source.
 * @param wrapper The destination.
 */
void importMaterialField( std::vector< vtkIdType > const & cellIds,
                          vtkDataArray * vtkArray,
                          WrapperBase & wrapper )
{
  // Scalar material fields are stored as 2D arrays, vector/tensor are 3D
  using ImportTypes = types::ArrayTypes< types::RealTypes, types::DimsRange< 2, 3 > >;
  types::dispatch( ImportTypes{}, wrapper.getTypeId(), true, [&]( auto array )
  {
    using ArrayType = decltype( array );
    Wrapper< ArrayType > & wrapperT = Wrapper< ArrayType >::cast( wrapper );
    auto const view = wrapperT.reference().toView();

    localIndex const numComponentsSrc = LvArray::integerConversion< localIndex >( vtkArray->GetNumberOfComponents() );
    localIndex const numComponentsDst = wrapperT.numArrayComp() / view.size( 1 );
    GEOSX_ERROR_IF_NE_MSG( numComponentsDst, numComponentsSrc,
                           "Mismatch in number of components for field " << vtkArray->GetName() );

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
  } );
}

/**
 * @brief Imports 1d and 2d arrays from @p vtkArray to @p wrapper, only for @p cellIds
 * @param cellIds The cells for which we should copy the data.
 * @param vtkArray The source.
 * @param wrapper The destination.
 */
void importRegularField( std::vector< vtkIdType > const & cellIds,
                         vtkDataArray * vtkArray,
                         WrapperBase & wrapper )
{
  using ImportTypes = types::ArrayTypes< types::RealTypes, types::DimsRange< 1, 2 > >;
  types::dispatch( ImportTypes{}, wrapper.getTypeId(), true, [&]( auto dstArray )
  {
    using ArrayType = decltype( dstArray );
    Wrapper< ArrayType > & wrapperT = Wrapper< ArrayType >::cast( wrapper );
    auto const view = wrapperT.reference().toView();

    localIndex const numComponentsSrc = LvArray::integerConversion< localIndex >( vtkArray->GetNumberOfComponents() );
    localIndex const numComponentsDst = wrapperT.numArrayComp();
    GEOSX_ERROR_IF_NE_MSG( numComponentsDst, numComponentsSrc,
                           "Mismatch in number of components for field " << vtkArray->GetName() );

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
  } );
}

void printMeshStatistics( vtkUnstructuredGrid & mesh,
                          VTKMeshGenerator::CellMapType const & cellMap,
                          MPI_Comm const comm )
{
  int const rank = MpiWrapper::commRank( comm );
  int const size = MpiWrapper::commSize( comm );

  vtkIdTypeArray const & globalPointId = *vtkIdTypeArray::FastDownCast( mesh.GetPointData()->GetGlobalIds() );
  RAJA::ReduceMax< parallelHostReduce, globalIndex > maxGlobalNode( -1 );
  forAll< parallelHostPolicy >( mesh.GetNumberOfPoints(), [&globalPointId, maxGlobalNode]( vtkIdType const k )
  {
    maxGlobalNode.max( globalPointId.GetValue( k ) );
  } );
  globalIndex const numGlobalNodes = MpiWrapper::max( maxGlobalNode.get(), comm ) + 1;

  localIndex numLocalElems = 0;
  globalIndex numGlobalElems = 0;
  std::map< ElementType, globalIndex > elemCounts;

  for( auto const & typeToCells : cellMap )
  {
    localIndex const localElemsOfType =
      std::accumulate( typeToCells.second.begin(), typeToCells.second.end(), localIndex{},
                       []( auto const s, auto const & region ) { return s + region.second.size(); } );
    numLocalElems += localElemsOfType;

    globalIndex const globalElemsOfType = MpiWrapper::sum( globalIndex{ localElemsOfType }, comm );
    numGlobalElems += globalElemsOfType;
    elemCounts[typeToCells.first] = globalElemsOfType;
  }

  localIndex const minLocalElems = MpiWrapper::min( numLocalElems );
  localIndex const maxLocalElems = MpiWrapper::max( numLocalElems );
  localIndex const avgLocalElems = LvArray::integerConversion< localIndex >( numGlobalElems / size );

  if( rank == 0 )
  {
    int const widthGlobal = static_cast< int >( std::log10( std::max( numGlobalElems, numGlobalNodes ) ) + 1 );
    GEOSX_LOG( GEOSX_FMT( "Number of nodes: {:>{}}", numGlobalNodes, widthGlobal ) );
    GEOSX_LOG( GEOSX_FMT( "Number of elems: {:>{}}", numGlobalElems, widthGlobal ) );
    for( auto const & typeCount: elemCounts )
    {
      GEOSX_LOG( GEOSX_FMT( "{:>15}: {:>{}}", toString( typeCount.first ), typeCount.second, widthGlobal ) );
    }

    int const widthLocal = static_cast< int >( std::log10( maxLocalElems ) + 1 );
    GEOSX_LOG( GEOSX_FMT( "Load balancing: {1:>{0}} {2:>{0}} {3:>{0}}\n"
                          "(element/rank): {4:>{0}} {5:>{0}} {6:>{0}}",
                          widthLocal, "min", "avg", "max",
                          minLocalElems, avgLocalElems, maxLocalElems ) );
  }
}

/**
 * @brief Collect the data to be imported.
 * @return A list of pointers to VTK data arrays.
 */
std::vector< vtkDataArray * >
findArraysForImport( vtkUnstructuredGrid & mesh,
                     arrayView1d< string const > const & srcFieldNames )
{
  std::vector< vtkDataArray * > arrays;
  vtkCellData & cellData = *mesh.GetCellData();

  for( string const & sourceName : srcFieldNames )
  {
    vtkAbstractArray * const curArray = cellData.GetAbstractArray( sourceName.c_str() );
    GEOSX_THROW_IF( curArray == nullptr,
                    GEOSX_FMT( "Source field '{}' not found in dataset", sourceName ),
                    InputError );

    int const dataType = curArray->GetDataType();
    GEOSX_ERROR_IF( dataType != VTK_FLOAT && dataType != VTK_DOUBLE,
                    GEOSX_FMT( "Source field '{}' has unsupported type: {} (expected floating point type)",
                               sourceName, curArray->GetDataTypeAsString() ) );
    arrays.push_back( vtkDataArray::SafeDownCast( curArray ) );
  }

  return arrays;
}

} // namespace vtk

void VTKMeshGenerator::
  importFieldOnCellElementSubRegion( int const regionId,
                                     ElementType const elemType,
                                     std::vector< vtkIdType > const & cellIds,
                                     ElementRegionManager & elemManager,
                                     arrayView1d< string const > const & fieldNames,
                                     std::vector< vtkDataArray * > const & srcArrays,
                                     FieldIdentifiers & fieldsToBeSync ) const
{
  string const cellBlockName = vtk::buildCellBlockName( elemType, regionId );

  elemManager.forElementSubRegionsComplete< CellElementSubRegion >( [&]( localIndex,
                                                                         localIndex,
                                                                         ElementRegionBase const & region,
                                                                         CellElementSubRegion & subRegion )
  {
    // We don't know how the user mapped cell blocks to regions, so we must check all of them
    if( subRegion.getName() != cellBlockName )
    {
      return;
    }
    std::unordered_set< string > const materialWrapperNames = vtk::getMaterialWrapperNames( subRegion );

    // Writing properties
    for( localIndex i = 0; i < fieldNames.size(); ++i )
    {
      // Get source
      vtkDataArray * vtkArray = srcArrays[i];

      // Find destination
      string const wrapperName = fieldNames[i];
      if( !subRegion.hasWrapper( wrapperName ) )
      {
        // Skip - the user may have not enabled a particular physics model/solver on this dstRegion.
        GEOSX_LOG_LEVEL_RANK_0( 1, GEOSX_FMT( "Skipping import of {} -> {} on {}/{} (field not found)",
                                              vtkArray->GetName(), wrapperName, region.getName(), subRegion.getName() ) );

        continue;
      }

      // Now that we know that the subRegion has this wrapper, we can add the wrapperName to the list of fields to synchronize
      fieldsToBeSync.addElementFields( {wrapperName}, {region.getName()} );

      WrapperBase & wrapper = subRegion.getWrapperBase( wrapperName );

      GEOSX_LOG_LEVEL_RANK_0( 1, GEOSX_FMT( "Importing field {} -> {} on {}/{}",
                                            vtkArray->GetName(), wrapperName, region.getName(), subRegion.getName() ) );

      if( materialWrapperNames.count( wrapperName ) > 0 && wrapper.numArrayDims() > 1 )
      {
        vtk::importMaterialField( cellIds, vtkArray, wrapper );
      }
      else
      {
        vtk::importRegularField( cellIds, vtkArray, wrapper );
      }
    }
  } );
}

void VTKMeshGenerator::importFields( DomainPartition & domain ) const
{
  GEOSX_LOG_RANK_0( GEOSX_FMT( "{} '{}': importing field data from mesh dataset", catalogName(), getName() ) );
  GEOSX_ASSERT_MSG( m_vtkMesh, "Must call generateMesh() before importFields()" );

  // TODO Having CellElementSubRegion and ConstitutiveBase... here in a pure geometric module is problematic.
  ElementRegionManager & elemManager = domain.getMeshBody( this->getName() ).getMeshLevel( 0 ).getElemManager();

  std::vector< vtkDataArray * > const srcArrays = vtk::findArraysForImport( *m_vtkMesh, m_fieldsToImport );

  FieldIdentifiers fieldsToBeSync;

  for( auto const & typeRegions : m_cellMap )
  {
    // Restrict data import to 3D cells
    if( getElementDim( typeRegions.first ) == 3 )
    {
      for( auto const & regionCells: typeRegions.second )
      {
        importFieldOnCellElementSubRegion( regionCells.first,
                                           typeRegions.first,
                                           regionCells.second,
                                           elemManager,
                                           m_fieldNamesInGEOSX,
                                           srcArrays,
                                           fieldsToBeSync );
      }
    }
  }

  CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                       domain.getMeshBody( this->getName() ).getMeshLevel( 0 ),
                                                       domain.getNeighbors(),
                                                       false );
}

real64 VTKMeshGenerator::writeNodes( CellBlockManager & cellBlockManager ) const
{
  localIndex const numPts = LvArray::integerConversion< localIndex >( m_vtkMesh->GetNumberOfPoints() );
  cellBlockManager.setNumNodes( numPts );

  // Writing the points
  arrayView1d< globalIndex > const nodeLocalToGlobal = cellBlockManager.getNodeLocalToGlobal();
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const X = cellBlockManager.getNodePositions();

  std::unordered_set< globalIndex > nodeGlobalIds;
  nodeGlobalIds.reserve( numPts );

  vtkIdTypeArray const & globalPointId = *vtkIdTypeArray::FastDownCast( m_vtkMesh->GetPointData()->GetGlobalIds() );
  forAll< serialPolicy >( numPts, [&, X, nodeLocalToGlobal]( localIndex const k )
  {
    double point[3];
    m_vtkMesh->GetPoint( k, point );
    LvArray::tensorOps::add< 3 >( point, m_translate );
    LvArray::tensorOps::hadamardProduct< 3 >( X[k], point, m_scale );
    globalIndex const pointGlobalID = globalPointId.GetValue( k );
    nodeLocalToGlobal[k] = pointGlobalID;

    // TODO: remove this check once the input mesh is cleaned of duplicate points via a filter
    //       and make launch policy parallel again
    GEOSX_ERROR_IF( nodeGlobalIds.count( pointGlobalID ) > 0,
                    GEOSX_FMT( "Duplicate point detected: globalID = {}\n"
                               "Consider cleaning the dataset in Paraview using 'Clean to grid' filter.\n"
                               "Make sure partitionRefinement is set to 1 or higher (this may help).",
                               pointGlobalID ) );
    nodeGlobalIds.insert( pointGlobalID );
  } );

  // Generate the "all" set
  array1d< localIndex > allNodes( numPts );
  std::iota( allNodes.begin(), allNodes.end(), 0 );
  SortedArray< localIndex > & allNodeSet = cellBlockManager.getNodeSets()[ "all" ];
  allNodeSet.insert( allNodes.begin(), allNodes.end() );

  constexpr real64 minReal = LvArray::NumericLimits< real64 >::min;
  constexpr real64 maxReal = LvArray::NumericLimits< real64 >::max;
  real64 xMin[3] = { maxReal, maxReal, maxReal };
  real64 xMax[3] = { minReal, minReal, minReal };

  vtkBoundingBox bb( m_vtkMesh->GetBounds() );
  if( bb.IsValid() )
  {
    bb.GetMinPoint( xMin );
    bb.GetMaxPoint( xMax );
  }

  MpiWrapper::min< real64 >( xMin, xMin, MPI_COMM_GEOSX );
  MpiWrapper::max< real64 >( xMax, xMax, MPI_COMM_GEOSX );
  LvArray::tensorOps::subtract< 3 >( xMax, xMin );
  return LvArray::tensorOps::l2Norm< 3 >( xMax );
}

/**
 * @brief Build all the cell blocks.
 * @param[in] cellBlockManager The instance that stores the cell blocks.
 */
void VTKMeshGenerator::writeCells( CellBlockManager & cellBlockManager ) const
{
  // Creates a new cell block for each region and for each type of cell.
  for( auto const & typeRegions : m_cellMap )
  {
    ElementType const elemType = typeRegions.first;
    if( getElementDim( elemType ) != 3 )
    {
      continue;
    }
    std::unordered_map< int, std::vector< vtkIdType > > const & regionIdToCellIds = typeRegions.second;
    for( auto const & regionCells : regionIdToCellIds )
    {
      int const regionId = regionCells.first;
      std::vector< vtkIdType > const & cellIds = regionCells.second;

      string const cellBlockName = vtk::buildCellBlockName( elemType, regionId );
      GEOSX_LOG_LEVEL_RANK_0( 1, "Importing cell block " << cellBlockName );

      // Create and resize the cell block.
      CellBlock & cellBlock = cellBlockManager.registerCellBlock( cellBlockName );
      cellBlock.setElementType( elemType );
      cellBlock.resize( LvArray::integerConversion< localIndex >( cellIds.size() ) );

      vtk::fillCellBlock( *m_vtkMesh, elemType, cellIds, cellBlock );
    }
  }
}

/**
 * @brief Build the "surface" node sets from the surface information.
 * @param[in] mesh The vtkUnstructuredGrid that is loaded
 * @param[in] surfacesIdsToCellsIds Map from the surfaces index to the list of cells in this surface in this rank.
 * @param[out] cellBlockManager The instance that stores the node sets.
 * @note @p surfacesIdsToCellsIds will contain all the surface ids across all the MPI ranks, but only its cell ids.
 * If the current MPI rank has no cell id for a given surface, then an empty set will be created.
 */
void VTKMeshGenerator::writeSurfaces( CellBlockManager & cellBlockManager ) const
{
  if( m_cellMap.count( ElementType::Polygon ) == 0 )
  {
    return;
  }
  std::map< string, SortedArray< localIndex > > & nodeSets = cellBlockManager.getNodeSets();

  for( auto const & surfaceCells: m_cellMap.at( ElementType::Polygon ) )
  {
    int const surfaceId = surfaceCells.first;
    std::vector< vtkIdType > const & cellIds = surfaceCells.second;
    string const surfaceName = std::to_string( surfaceId );
    GEOSX_LOG_LEVEL_RANK_0( 1, "Importing surface " << surfaceName );

    // Get or create all surfaces (even those which are empty in this rank)
    SortedArray< localIndex > & curNodeSet = nodeSets[ surfaceName ];

    for( vtkIdType const c : cellIds )
    {
      vtkCell * const currentCell = m_vtkMesh->GetCell( c );
      for( int v = 0; v < currentCell->GetNumberOfPoints(); ++v )
      {
        curNodeSet.insert( currentCell->GetPointId( v ) );
      }
    }
  }
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
    vtkSmartPointer< vtkUnstructuredGrid > loadedMesh = vtk::loadMesh( m_filePath );
    GEOSX_LOG_LEVEL_RANK_0( 2, "  redistributing mesh..." );
    m_vtkMesh = vtk::redistributeMesh( *loadedMesh, comm, m_partitionRefinement );
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
  real64 const globalLength = writeNodes( cellBlockManager );
  meshBody.setGlobalLengthScale( globalLength );

  GEOSX_LOG_LEVEL_RANK_0( 2, "  writing cells..." );
  writeCells( cellBlockManager );

  GEOSX_LOG_LEVEL_RANK_0( 2, "  writing surfaces..." );
  writeSurfaces( cellBlockManager );

  GEOSX_LOG_LEVEL_RANK_0( 2, "  building connectivity maps..." );
  cellBlockManager.buildMaps();

  GEOSX_LOG_LEVEL_RANK_0( 2, "  done!" );
  vtk::printMeshStatistics( *m_vtkMesh, m_cellMap, comm );
}

void VTKMeshGenerator::freeResources()
{
  m_vtkMesh = nullptr;
  m_cellMap.clear();
}

REGISTER_CATALOG_ENTRY( MeshGeneratorBase, VTKMeshGenerator, string const &, Group * const )

} // namespace geosx
