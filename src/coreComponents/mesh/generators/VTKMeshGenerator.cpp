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

#include "mesh/DomainPartition.hpp"
#include "mesh/generators/CellBlockManager.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/MeshBody.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"

#include "common/MpiWrapper.hpp"
#include "common/TypeDispatch.hpp"

#include "common/DataTypes.hpp"
#include "common/DataLayouts.hpp"

#include <vtkArrayDispatch.h>
#include <vtkBoundingBox.h>
#include <vtkCellData.h>
#include <vtkDIYExplicitAssigner.h>
#include <vtkGenerateGlobalIds.h>
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
    setDescription( "Translate the coordinates of the vertices by a given vector" );
}

namespace
{

/**
 * @brief Return a VTK controller for multiprocessing.
 */
vtkSmartPointer< vtkMultiProcessController > getVTKController()
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
loadVTKMesh( Path const & filePath )
{
  // TODO maybe code a method in Path to get the file extension?
  string const extension = stringutilities::tokenize( filePath, "." ).back();
  vtkSmartPointer< vtkUnstructuredGrid > loadedMesh;

  if( extension == "pvtu" )
  {
    vtkSmartPointer< vtkXMLPUnstructuredGridReader > vtkUgReader = vtkSmartPointer< vtkXMLPUnstructuredGridReader >::New();
    vtkUgReader->SetFileName( filePath.c_str() );
    vtkUgReader->UpdateInformation();
    int const numberOfPieces = vtkUgReader->GetNumberOfPieces();
    if( MpiWrapper::commSize() == 1 )
    {
      vtkUgReader->Update();
    }
    else if( MpiWrapper::commRank() < numberOfPieces )
    {
      vtkUgReader->UpdatePiece( MpiWrapper::commRank(), 2, 0 ); // TODO What's that numPieces = 2?
    }
    loadedMesh = vtkUgReader->GetOutput();
  }
  else
  {
    if( MpiWrapper::commRank() == 0 )
    {
      auto read = [&]( auto vtkUgReader ) // auto can be multiple types in the same function
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
        GEOSX_ERROR( extension << " is not a recognized extension for using the VTK reader with GEOSX. Please use .vtk or .vtu" );
      }
    }
    else
    {
      loadedMesh = vtkSmartPointer< vtkUnstructuredGrid >::New();
    }
  }

  return loadedMesh;
}

/**
 * @brief Compute assignments of mesh partitions produced by vtkRedistributeDataSetFilter to ranks.
 * @param numParts number of partitions (cuts)
 * @param numRanks number of MPI ranks
 * @return a vector of rank assignments
 *
 * This function merges partitions starting from the back, until the desired number is reached.
 * It is a stripped down version of vtkDIYKdTreeUtilities::ComputeAssignments(). We cannot call
 * that function directly because it's part of a VTK private module (i.e. header not installed).
 * We rely on copying the assignment procedure for lack of a better solution.
 * In particular, there is no API to extract this information out of vtkRedistributeDataSetFilter.
 */
std::vector< int > computePartitionAssignments( int const numParts, int const numRanks )
{
  std::vector< int > assignments;
  assignments.reserve( numParts );
  std::div_t const d = std::div( numParts, numRanks );

  for( int i = 0; i < numRanks; ++i )
  {
    // This puts more parts on higher ranks, fewer parts on lower.
    int const numLocal = d.quot + ( i < numRanks - d.rem ? 0 : 1 );
    for( int j = 0; j < numLocal; ++j )
    {
      assignments.push_back( i );
    }
  }

  return assignments;
}


/**
 * @brief Redistribute the mesh among the available MPI ranks
 * @details this method will also generate global ids for points and cells in the VTK Mesh
 * @param[in] loadedMesh the mesh that was loaded on one or several MPI ranks
 * @param[out] cuts the bounding boxes used by the VTK partitioner needed to compute neighboring ranks
 */
vtkSmartPointer< vtkUnstructuredGrid >
redistributeMesh( vtkUnstructuredGrid & loadedMesh,
                  std::vector< vtkBoundingBox > & cuts )
{
  // Generate global IDs for vertices and cells
  vtkNew< vtkGenerateGlobalIds > generator;
  generator->SetInputDataObject( &loadedMesh );
  generator->Update();

  // Redistribute data all over the available ranks
  vtkNew< vtkRedistributeDataSetFilter > rdsf;
  rdsf->SetInputDataObject( generator->GetOutputDataObject( 0 ) );
  rdsf->SetNumberOfPartitions( MpiWrapper::commSize() );
  rdsf->Update();

  cuts = rdsf->GetCuts();
  return vtkUnstructuredGrid::SafeDownCast( rdsf->GetOutputDataObject( 0 ) );
}

/**
 * @brief Compute the potential rank neighbor list
 * @param[in] cuts the bounding boxes used by the VTK partitioner for all ranks
 * @details Compute the rank neighbors. The assumption is that 2 ranks are neighbors if
 * the corresponding bounding boxes intersect.
 */
std::unordered_set< int > computeMPINeighborRanks( std::vector< vtkBoundingBox > cuts )
{
  int const numParts = LvArray::integerConversion< int >( cuts.size() );
  int const numRanks = MpiWrapper::commSize();
  int const thisRank = MpiWrapper::commRank();

  double constexpr inflateFactor = 1.1;

  std::vector< int > const assignments = computePartitionAssignments( numParts, numRanks );
  std::vector< vtkBoundingBox > myCuts;
  for( int i = 0; i < numParts; ++i )
  {
    cuts[i].ScaleAboutCenter( inflateFactor );
    if( assignments[i] == thisRank )
    {
      myCuts.push_back( cuts[i] );
    }
  }

  std::unordered_set< int > neighbors;
  for( int i = 0; i < numParts; ++i )
  {
    if( assignments[i] != thisRank && std::any_of( myCuts.begin(), myCuts.end(),
                                                   [&]( auto const & b ) { return b.Intersects( cuts[i] ); } ) )
    {
      neighbors.insert( assignments[i] );
    }
  }

  return neighbors;
}

/**
 * @brief Copy the VTK mesh nodes into the nodeManager of GEOSX
 * @param[in] nodeManager the NodeManager of the domain in which the points will be copied.
 * @param[in] mesh the vtkUnstructuredGrid that is loaded
 * @return the global length of the mesh (diagonal of the bounding box)
 */
real64 writeMeshNodes( vtkUnstructuredGrid & mesh,
                       R1Tensor const & scale,
                       R1Tensor const & translate,
                       CellBlockManager & cellBlockManager )
{
  vtkIdType const numPts = mesh.GetNumberOfPoints();
  cellBlockManager.setNumNodes( numPts );

  // Writing the points
  arrayView1d< globalIndex > const & nodeLocalToGlobal = cellBlockManager.getNodeLocalToGlobal();
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & X = cellBlockManager.getNodePositions();

  constexpr real64 minReal = LvArray::NumericLimits< real64 >::min;
  constexpr real64 maxReal = LvArray::NumericLimits< real64 >::max;
  real64 xMin[3] = { maxReal, maxReal, maxReal };
  real64 xMax[3] = { minReal, minReal, minReal };

  vtkIdTypeArray const * const globalPointId = vtkIdTypeArray::FastDownCast( mesh.GetPointData()->GetGlobalIds() );
  GEOSX_ERROR_IF( numPts > 0 && globalPointId == nullptr, "Global point IDs have not been generated" );

  for( vtkIdType v = 0; v < numPts; ++v )
  {
    double const * const point = mesh.GetPoint( v );
    nodeLocalToGlobal[v] = globalPointId->GetValue( v );
    for( integer i = 0; i < 3; ++i )
    {
      X( v, i ) = ( point[i] + translate[i] ) * scale[i];
      xMax[i] = std::max( xMax[i], X( v, i ) );
      xMin[i] = std::min( xMin[i], X( v, i ) );
    }
  }

  // Generate the "all" set
  array1d< localIndex > allNodes( numPts );
  std::iota( allNodes.begin(), allNodes.end(), 0 );
  SortedArray< localIndex > & allNodeSet = cellBlockManager.getNodeSets()[ "all" ];
  allNodeSet.insert( allNodes.begin(), allNodes.end() );

  MpiWrapper::allReduce( xMin, xMin, 3, MPI_MIN, MPI_COMM_GEOSX );
  MpiWrapper::allReduce( xMax, xMax, 3, MPI_MAX, MPI_COMM_GEOSX );
  LvArray::tensorOps::subtract< 3 >( xMax, xMin );
  return LvArray::tensorOps::l2Norm< 3 >( xMax );
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
    case VTK_VERTEX:     return ElementType::Vertex;
    case VTK_LINE:       return ElementType::Line;
    case VTK_TRIANGLE:   return ElementType::Triangle;
    case VTK_QUAD:       return ElementType::Quadrilateral;
    case VTK_POLYGON:    return ElementType::Polygon;
    case VTK_TETRA:      return ElementType::Tetrahedron;
    case VTK_PYRAMID:    return ElementType::Pyramid;
    case VTK_WEDGE:      return ElementType::Prism;
    case VTK_HEXAHEDRON: return ElementType::Hexahedron;
    case VTK_POLYHEDRON: return ElementType::Polyhedron;
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
    case ElementType::Tetrahedron:   return { 1, 0, 2, 3 };
    case ElementType::Pyramid:       return { 0, 1, 2, 3, 4, 5 }; // TODO check
    case ElementType::Prism:         return { 0, 3, 1, 4, 2, 5 };
    case ElementType::Hexahedron:    return { 0, 1, 3, 2, 4, 5, 7, 6 };
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
    case ElementType::Prism:       return "wedges";
    case ElementType::Pyramid:     return "pyramids";
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
      for( vtkIdType c: cellIds )
      {
        for( localIndex q = 0; q < view.size( 1 ); ++q )
        {
          // The same value is copied for all quadrature points.
          int componentIdx = 0;
          LvArray::forValuesInSlice( view[cellCount][q], [&]( auto & val )
          {
            val = data.Get( c, componentIdx++ );
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
      for( vtkIdType c: cellIds )
      {
        int componentIdx = 0;
        LvArray::forValuesInSlice( view[cellCount], [&]( auto & val )
        {
          val = data.Get( c, componentIdx++ );
        } );
        ++cellCount;
      }
    } );
  } );
}

void printMeshStatistics( CellBlockManager & cellBlockManager,
                          VTKMeshGenerator::CellMapType const & cellMap )
{
  // TODO cellBlockManager should be const, but it's API is weird...
  arrayView1d< globalIndex const > const nodeGlobalIndices = cellBlockManager.getNodeLocalToGlobal();
  auto const maxIter = std::max_element( nodeGlobalIndices.begin(), nodeGlobalIndices.end() );
  globalIndex const maxLocalNode = ( maxIter != nodeGlobalIndices.end() ) ? *maxIter : -1;
  globalIndex const numGlobalNodes = MpiWrapper::max( maxLocalNode ) + 1;

  globalIndex numGlobalElems = 0;
  globalIndex maxTypeElems = 0;
  std::map< ElementType, globalIndex > elemCounts;
  for( auto const & typeToCells : cellMap )
  {
    localIndex const localElems =
      std::accumulate( typeToCells.second.begin(), typeToCells.second.end(), localIndex{},
                       []( auto const s, auto const & region ) { return s + region.second.size(); } );
    globalIndex const globalElems = MpiWrapper::sum( LvArray::integerConversion< globalIndex >( localElems ) );
    elemCounts[typeToCells.first] = globalElems;
    numGlobalElems += globalElems;
    maxTypeElems = std::max( maxTypeElems, globalElems );
  }
  int const width = static_cast< int >( std::log10( maxTypeElems ) + 1 );
  std::ostringstream elemStats;
  for( auto const & typeCount : elemCounts )
  {
    elemStats << GEOSX_FMT( "\n{:>15}: {:>{}}", toString( typeCount.first ), typeCount.second, width );
  }
  GEOSX_LOG_RANK_0( GEOSX_FMT( "Number of nodes: {}\nNumber of elems: {}{}",
                               numGlobalNodes, numGlobalElems, elemStats.str() ) );
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

} // namespace

void VTKMeshGenerator::
  importFieldOnCellElementSubRegion( int const regionId,
                                     ElementType const elemType,
                                     std::vector< vtkIdType > const & cellIds,
                                     ElementRegionManager & elemManager,
                                     arrayView1d< string const > const & fieldNames,
                                     std::vector< vtkDataArray * > const & srcArrays ) const
{
  string const cellBlockName = buildCellBlockName( elemType, regionId );

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
    std::unordered_set< string > const materialWrapperNames = getMaterialWrapperNames( subRegion );

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
      WrapperBase & wrapper = subRegion.getWrapperBase( wrapperName );

      GEOSX_LOG_LEVEL_RANK_0( 1, GEOSX_FMT( "Importing field {} -> {} on {}/{}",
                                            vtkArray->GetName(), wrapperName, region.getName(), subRegion.getName() ) );

      if( materialWrapperNames.count( wrapperName ) > 0 && wrapper.numArrayDims() > 1 )
      {
        importMaterialField( cellIds, vtkArray, wrapper );
      }
      else
      {
        importRegularField( cellIds, vtkArray, wrapper );
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

  std::vector< vtkDataArray * > const srcArrays = findArraysForImport( *m_vtkMesh, m_fieldsToImport );

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
                                           srcArrays );
      }
    }
  }

  CommunicationTools::getInstance().synchronizeFields( { { "elems", m_fieldNamesInGEOSX } },
                                                       domain.getMeshBody( this->getName() ).getMeshLevel( 0 ),
                                                       domain.getNeighbors(),
                                                       false );
}

/**
 * @brief Build all the cell blocks.
 * @param[in] mesh the vtkUnstructuredGrid that is loaded
 * @param[in] regionsHex map from region index to the hexahedron indexes in this region
 * @param[in] regionsTetra map from region index to the tetra indexes in this region
 * @param[in] regionsWedges map from region index to the wedges indexes in this region
 * @param[in] regionsPyramids map from region index to the pyramids indexes in this region
 * @param[out] cellBlockManager The instance that stores the cell blocks.
 */
void VTKMeshGenerator::buildCellBlocks( CellBlockManager & cellBlockManager ) const
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
    for( auto const & regionCells: regionIdToCellIds )
    {
      int const regionId = regionCells.first;
      std::vector< vtkIdType > const & cellIds = regionCells.second;

      string const cellBlockName = buildCellBlockName( elemType, regionId );
      GEOSX_LOG_LEVEL_RANK_0( 1, "Importing cell block " << cellBlockName );

      // Create and resize the cell block.
      CellBlock & cellBlock = cellBlockManager.registerCellBlock( cellBlockName );
      cellBlock.setElementType( elemType );
      cellBlock.resize( LvArray::integerConversion< localIndex >( cellIds.size() ) );

      fillCellBlock( *m_vtkMesh, elemType, cellIds, cellBlock );
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
void VTKMeshGenerator::buildSurfaces( CellBlockManager & cellBlockManager ) const
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

  vtkSmartPointer< vtkMultiProcessController > controller = getVTKController();
  vtkMultiProcessController::SetGlobalController( controller );

  GEOSX_LOG_RANK_0( GEOSX_FMT( "{} '{}': reading mesh from {}", catalogName(), getName(), m_filePath ) );
  {
    vtkSmartPointer< vtkUnstructuredGrid > loadedMesh = loadVTKMesh( m_filePath );
    std::vector< vtkBoundingBox > cuts;
    m_vtkMesh = redistributeMesh( *loadedMesh, cuts );

    // TODO Check that the neighbor information set is bulletproof
    std::unordered_set< int > const neighbors = computeMPINeighborRanks( std::move( cuts ) );
    domain.getMetisNeighborList().insert( neighbors.begin(), neighbors.end() );
  }

  GEOSX_LOG_RANK_0( GEOSX_FMT( "{} '{}': generating GEOSX mesh data structure", catalogName(), getName() ) );

  MeshBody & meshBody = domain.getMeshBodies().registerGroup< MeshBody >( this->getName() );
  meshBody.createMeshLevel( 0 );

  CellBlockManager & cellBlockManager = meshBody.registerGroup< CellBlockManager >( keys::cellManager );

  real64 const globalLength = writeMeshNodes( *m_vtkMesh, m_scale, m_translate, cellBlockManager );
  meshBody.setGlobalLengthScale( globalLength );

  m_cellMap = buildCellMap( *m_vtkMesh, m_attributeName );
  printMeshStatistics( cellBlockManager, m_cellMap );

  buildCellBlocks( cellBlockManager );
  buildSurfaces( cellBlockManager );

  // TODO Check the memory usage that seems prohibitive - Do we need to build all connections?
  cellBlockManager.buildMaps();
}

void VTKMeshGenerator::freeResources()
{
  m_vtkMesh = nullptr;
  m_cellMap.clear();
}

REGISTER_CATALOG_ENTRY( MeshGeneratorBase, VTKMeshGenerator, string const &, Group * const )

} // namespace geosx
