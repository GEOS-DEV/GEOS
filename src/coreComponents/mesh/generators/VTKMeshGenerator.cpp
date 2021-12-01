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
#include "mesh/mpiCommunications/PartitionBase.hpp"
#include "mesh/mpiCommunications/SpatialPartition.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/MeshBody.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"

#include "common/MpiWrapper.hpp"
#include "common/TypeDispatch.hpp"

#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLPUnstructuredGridReader.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkRedistributeDataSetFilter.h>
#include <vtkGenerateGlobalIds.h>
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



VTKMeshGenerator::VTKMeshGenerator( string const & name, Group * const parent ):
  MeshGeneratorBase( name, parent )
{

  registerWrapper( viewKeyStruct::filePathString(), &m_filePath ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "path to the mesh file" );
}

void VTKMeshGenerator::postProcessInput()
{}

Group * VTKMeshGenerator::createChild( string const & GEOSX_UNUSED_PARAM( childKey ), string const & GEOSX_UNUSED_PARAM( childName ) )
{
  return nullptr;
}

template< typename VTK_ARRAY >
VTK_ARRAY const * getDataArrayOptional( vtkFieldData & source, char const * name )
{
  // returns nullptr if either lookup or cast fail
  return VTK_ARRAY::FastDownCast( source.GetAbstractArray( name ) );
}

template< typename VTK_ARRAY >
VTK_ARRAY const & getDataArray( vtkFieldData & source, char const * name )
{
  VTK_ARRAY const * const array = getDataArrayOptional< VTK_ARRAY >( source, name );
  GEOSX_ERROR_IF( array == nullptr, "VTK array '" << name << "' not found or has unexpected data type" );
  return *array;
}

/**
 * @brief Check if two boxes are intersected
 * @param[in] box1 the first box [xMin xMax yMin yMax zMin zMax]
 * @param[in] box2 the second box [xMin xMax yMin yMax zMin zMax]
 * @return true if the boxes are intersecting, false otherwise
 */
bool checkIntersection( double box1[6], double box2[6] )
{
  return box1[0] < box2[1] && box1[1]> box2[0] && box1[2] < box2[3] && box1[3]> box2[2] && box1[4] < box2[5] && box1[5]> box2[4];
}

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
vtkSmartPointer< vtkUnstructuredGrid > loadVTKMesh( Path const & filePath )
{
  string_array filePathTokenized = stringutilities::tokenize( filePath, "." ); //TODO maybe code a method in Path to get the file extension?
  string extension = filePathTokenized[filePathTokenized.size() - 1];
  vtkSmartPointer< vtkUnstructuredGrid > loadedMesh = vtkSmartPointer< vtkUnstructuredGrid >::New();

  auto read = [&]( auto vtkUgReader ) { // auto can be multiple types in the same function
    vtkUgReader->SetFileName( filePath.c_str() );
    vtkUgReader->Update();
    loadedMesh = vtkUgReader->GetOutput();
  };

  if( extension == "pvtu" )
  {
    vtkSmartPointer< vtkXMLPUnstructuredGridReader > vtkUgReader = vtkSmartPointer< vtkXMLPUnstructuredGridReader >::New();
    vtkUgReader->SetFileName( filePath.c_str() );
    vtkUgReader->UpdateInformation();
    int numberOfPieces = vtkUgReader->GetNumberOfPieces();
    if( MpiWrapper::commSize() == 1 )
    {
      vtkUgReader->Update();
    }
    else
    {
      if( MpiWrapper::commRank() < numberOfPieces )
      {
        vtkUgReader->UpdatePiece( MpiWrapper::commRank(), 2, 0 );
      }
    }
    loadedMesh = vtkUgReader->GetOutput();
  }
  else
  {
    if( MpiWrapper::commRank() == 0 )
    {
      if( extension == "vtk" )
      {
        read( vtkSmartPointer< vtkUnstructuredGridReader >::New());
      }
      else if( extension == "vtu" )
      {
        read( vtkSmartPointer< vtkXMLUnstructuredGridReader >::New());
      }
      else
      {
        GEOSX_ERROR( extension << " is not a recognized extension for using the VTK reader with GEOSX. Please use .vtk or .vtu" );
      }
    }
  }

  return loadedMesh;
}

/**
 * @brief Redistribute the mesh among the available MPI ranks
 * @details this method will also generate global ids for points and cells in the VTK Mesh
 * @param[in] loadedMesh the mesh that was loaded on one or several MPI ranks
 */
vtkSmartPointer< vtkUnstructuredGrid > redistributeMesh( vtkUnstructuredGrid & loadedMesh )
{
  // Redistribute data all over the available ranks
  vtkNew< vtkRedistributeDataSetFilter > rdsf;
  rdsf->SetInputDataObject( &loadedMesh );
  rdsf->SetNumberOfPartitions( MpiWrapper::commSize() );
  GEOSX_LOG_RANK( "11111" );
  rdsf->GenerateGlobalCellIdsOn();
  GEOSX_LOG_RANK( "22222" );
  rdsf->Update();
  GEOSX_LOG_RANK( "33333" );
  MpiWrapper::barrier();

  // Generate global IDs for vertices and cells
  vtkNew< vtkGenerateGlobalIds > generator;
  generator->SetInputDataObject( vtkUnstructuredGrid::SafeDownCast( rdsf->GetOutputDataObject( 0 )));
  generator->Update();
  vtkSmartPointer< vtkUnstructuredGrid > mesh = vtkUnstructuredGrid::SafeDownCast( generator->GetOutputDataObject( 0 ));
  GEOSX_LOG_RANK( "number of cells after " << mesh->GetNumberOfCells());
  return mesh;
}

/**
 * @brief Copy the VTK mesh nodes into the nodeManager of GEOSX
 * @param[in] nodeManager the NodeManager of the domain in which the poiints will be copied.
 * @param[in] mesh the vtkUnstructuredGrid that is loaded
 * @return the global length of the mesh (diagonal of the bounding box)
 */
double writeMeshNodes( NodeManager & nodeManager, vtkUnstructuredGrid & mesh )
{
  nodeManager.resize( mesh.GetNumberOfPoints() );
  arrayView1d< globalIndex > const & nodeLocalToGlobal = nodeManager.localToGlobalMap();

  // Writing the points
  GEOSX_ERROR_IF( mesh.GetNumberOfPoints() == 0, "Mesh is empty, aborting" );

  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();
  Group & nodeSets = nodeManager.sets();
  SortedArray< localIndex > & allNodes  = nodeSets.registerWrapper< SortedArray< localIndex > >( "all" ).reference();
  real64 xMax[3] = { std::numeric_limits< real64 >::min(), std::numeric_limits< real64 >::min(), std::numeric_limits< real64 >::min() };
  real64 xMin[3] = { std::numeric_limits< real64 >::max(), std::numeric_limits< real64 >::max(), std::numeric_limits< real64 >::max()  };

  vtkIdTypeArray const & globalPointIdDataArray = getDataArray< vtkIdTypeArray >( *mesh.GetPointData(), "GlobalPointIds" );

  for( vtkIdType v = 0; v < mesh.GetNumberOfPoints(); v++ )
  {
    double * point = mesh.GetPoint( v );
    nodeLocalToGlobal[v] = globalPointIdDataArray.GetValue( v );
    for( integer i = 0; i < 3; i++ )
    {
      X( v, i ) = point[i];
      xMax[i] = std::max( xMax[i], point[i] );
      xMin[i] = std::min( xMin[i], point[i] );
    }
    allNodes.insert( v );
  }
  LvArray::tensorOps::subtract< 3 >( xMax, xMin );
  return LvArray::tensorOps::l2Norm< 3 >( xMax );
}

/**
 * @brief Compute the potential rank neighbor list
 * @param[in] rankBounds the boundary of the mesh in this rank
 * @param[in] globalLength the globalLength of the model
 * @details Compute the metisNeighbor. This method computes the bounding box
 * of each domains. If these boundings boxes are crossing, it is possible that the corresponsing
 * domains are neighbors.
 */
std::set< int > computePotentialNeighborLists( double const rankBounds[], double const globalLength )
{
  double haloSize = globalLength/1000.;
  double rankBoundsHalo[6];
  /// compute the halo
  for( unsigned int i = 0; i < 6; i++ )
  {
    rankBoundsHalo[i] = rankBounds[i];
    rankBoundsHalo[i] += i%2?-haloSize:haloSize;
  }

  std::vector< double > allBounds( 6 * MpiWrapper::commSize());
  MpiWrapper::allgather< double, double >( rankBounds, 6, allBounds.data(), 6, MPI_COMM_GEOSX );

  //Checking the intersection
  std::set< int > rankNeighbor;
  for( int i = 0; i < MpiWrapper::commSize(); i++ )
  {
    if( i == MpiWrapper::commRank())
    {
      continue;
    }
    std::vector< double >::const_iterator first = allBounds.begin() + 6*i;
    std::vector< double >::const_iterator last = allBounds.begin() + 6*(i+1);
    std::vector< double > curBound( first, last );
    if( checkIntersection( rankBoundsHalo, curBound.data()))
    {
      rankNeighbor.insert( i );
    }
  }
  return rankNeighbor;
}

template< typename TYPE >
std::vector< TYPE > communicateAttributes( std::vector< TYPE > const & attributeVector )
{
  std::vector< TYPE > allAttributes;
  array1d< int > attributeSizes( MpiWrapper::commSize());
  MpiWrapper::allGather( LvArray::integerConversion< int >( attributeVector.size()), attributeSizes, MPI_COMM_GEOSX );
  int const totalNbAttributeId = std::accumulate( attributeSizes.begin(), attributeSizes.end(), 0 );
  allAttributes.resize( totalNbAttributeId );
  std::vector< int > displacements( MpiWrapper::commSize(), 0 );
  std::partial_sum( attributeSizes.begin(), attributeSizes.end() - 1, displacements.begin() + 1 );
  MpiWrapper::allgatherv( attributeVector.data(), attributeVector.size(), allAttributes.data(), attributeSizes.data(), displacements.data(), MPI_COMM_GEOSX );

  std::sort( allAttributes.begin(), allAttributes.end());
  auto last = std::unique( allAttributes.begin(), allAttributes.end());
  allAttributes.erase( last, allAttributes.end());
  return allAttributes;
}

/**
 * @brief This method is used to preprocess the the VTK mesh and count the number of cells, facets, regions
 * and surfaces over the current MPI rank
 * @param[out] regionsHex map from region index to the hexahedron indexes
 * @param[out] regionsTetra map from region index to the tetra indexes
 * @param[out] regionsWedges map from region index to the wedges indexes
 * @param[out] regionsPyramids map from region index to the pyramids indexes
 * @param[out] regions a set containing all the region indexes
 * @param[out] surfaces a set containing all the surface indexes, from this MPI rank
 * @param[out] allSurfaces a vector containing all the surfaces among all the MPI rank
 * @param[in] mesh the vtkUnstructuredGrid that is loaded
 */
void countCellsAndFaces( std::map< int, std::vector< vtkIdType > > & regionsHex, std::map< int, std::vector< vtkIdType > > & regionsTetra,
                         std::map< int, std::vector< vtkIdType > > & regionsWedges, std::map< int, std::vector< vtkIdType > > & regionsPyramids,
                         std::map< int, std::vector< vtkIdType > > & surfaces, std::vector< int > & allSurfaces,
                         vtkUnstructuredGrid & mesh )
{
  vtkIntArray const * const attributeDataArray = getDataArrayOptional< vtkIntArray >( *mesh.GetCellData(), "attribute" );

  auto countCell = [&]( localIndex c, std::map< int, std::vector< vtkIdType > > & regionsElem ) { // auto can be multiple types in the same
                                                                                                  // function
    if( attributeDataArray != nullptr )
    {
      int region = attributeDataArray->GetValue( c );
      GEOSX_LOG_RANK_0_IF( region < 0, "Attribute value " << region << " not supported, please use value >=0 to describe regions" );
      regionsElem[region].push_back( c );
    }
    else
    {
      regionsElem[-1].push_back( c );
    }
  };


  for( vtkIdType c = 0; c < mesh.GetNumberOfCells(); c++ )
  {
    if( mesh.GetCellType( c ) == VTK_HEXAHEDRON )
    {
      countCell( c, regionsHex );
    }
    else if( mesh.GetCellType( c ) == VTK_TETRA )
    {
      countCell( c, regionsTetra );
    }
    else if( mesh.GetCellType( c ) == VTK_WEDGE )
    {
      countCell( c, regionsWedges );
    }
    else if( mesh.GetCellType( c ) == VTK_PYRAMID )
    {
      countCell( c, regionsPyramids );
    }
    else if( mesh.GetCellType( c ) == VTK_TRIANGLE || mesh.GetCellType( c ) == VTK_QUAD )
    {
      if( attributeDataArray != nullptr )
      {
        surfaces[attributeDataArray->GetValue( c )].push_back( c );
      }
    }
  }

  // Communicate all the cells blocks
  std::vector< int > cellBlockRegionIndex( regionsHex.size() + regionsTetra.size() +
                                           regionsPyramids.size() + regionsWedges.size() );

  std::vector< ElementType > cellBlockElementType( regionsHex.size() + regionsTetra.size() +
                                                   regionsPyramids.size() + regionsWedges.size() );


  unsigned int count =0;
  auto fillCellBlockNames = [&]( std::map< int, std::vector< vtkIdType > > & regionsElem, ElementType type ) {
    for( auto const & region : regionsElem )
    {
      cellBlockRegionIndex[count] = region.first;
      cellBlockElementType[count++] = type;
    }
  };
  fillCellBlockNames( regionsHex, ElementType::Hexahedron );
  fillCellBlockNames( regionsTetra, ElementType::Tetrahedron );
  fillCellBlockNames( regionsWedges, ElementType::Prism );
  fillCellBlockNames( regionsPyramids, ElementType::Pyramid );

  // Comunicate all the region names
  // TODO if in a mesh, there is for instance a region with tetra, and a region without tetra,
  // they will both contain a cell block tetrahedron
  std::vector< int > allCellBlockRegionIndex = communicateAttributes( cellBlockRegionIndex );
  std::vector< ElementType > allCellBlockElementType = communicateAttributes( cellBlockElementType );
  for( unsigned int i = 0; i < allCellBlockRegionIndex.size(); i++ )
  {
    for( unsigned int j = 0; j < allCellBlockElementType.size(); j++ )
    {
      string name = allCellBlockRegionIndex[i] != -1 ? std::to_string( allCellBlockRegionIndex[i] ) + "_" : "";
      switch( allCellBlockElementType[j] )
      {
        case ElementType::Tetrahedron:
          regionsTetra[allCellBlockRegionIndex[i]];
          break;
        case ElementType::Hexahedron:
          regionsHex[allCellBlockRegionIndex[i]];
          break;
        case ElementType::Prism:
          regionsWedges[allCellBlockRegionIndex[i]];
          break;
        case ElementType::Pyramid:
          regionsPyramids[allCellBlockRegionIndex[i]];
          break;
        default:
          GEOSX_ERROR( allCellBlockElementType[j] << " is not a recognized cell type" );
      }
    }
  }

  // Communicate all the surfaces attributes to the ranks
  std::vector< int > surfaceVector;
  surfaceVector.reserve( surfaces.size());
  for( auto & surface : surfaces )
  {
    surfaceVector.push_back( surface.first );
  }
  allSurfaces = communicateAttributes( surfaceVector );
}

/**
 * @brief Find the properties to be imported
 * @details all the float and double vtkArrays will be imported
 * @return a vector containing all the vtkDataArray that can be imported
 */
std::vector< vtkDataArray * > findArrayToBeImported( vtkUnstructuredGrid & mesh )
{
  std::vector< vtkDataArray * > arrayToBeImported;
  vtkCellData * cellData = mesh.GetCellData();
  for( int i = 0; i < cellData->GetNumberOfArrays(); i++ )
  {
    vtkAbstractArray * curArray = cellData->GetAbstractArray( i );
    // We support only float and double type. Good assumption?
    int dataType = curArray->GetDataType();
    if( !curArray->IsNumeric() )
    {
      GEOSX_LOG_RANK_0( curArray->GetName() << " is not a numeric array" );
      continue;
    }
    if( dataType != VTK_FLOAT && dataType != VTK_DOUBLE )
    {
      GEOSX_LOG_RANK_0( "Underlying data Type " << curArray->GetDataTypeAsString() << " for array "
                                                << curArray->GetName() << " is not supported by GEOSX and will be ignored "
                                                << "(Only double and float are supported)" );
      continue;
    }
    GEOSX_LOG_RANK_0( "Importing array " << curArray->GetName() );
    arrayToBeImported.push_back( vtkDataArray::SafeDownCast( curArray ) );
  }
  return arrayToBeImported;
}

/**
 * @brief Get the number of points of a cell knowing its vtk cell type
 * @param[in] cellType the vtk cell type
 * @return the number of points contained in the cell
 */
localIndex getNumberOfPoints( int cellType )
{
  localIndex numberOfPoints = 0;
  switch( cellType )
  {
    case VTK_HEXAHEDRON:
      numberOfPoints = 8;
      break;
    case VTK_TETRA:
      numberOfPoints = 4;
      break;
    case VTK_WEDGE:
      numberOfPoints = 6;
      break;
    case VTK_PYRAMID:
      numberOfPoints = 5;
      break;
    default:
      GEOSX_ERROR( cellType << " is not a recognized cell type to be used with the VTKMeshGenerator" );
  }
  return numberOfPoints;
}

/**
 * @brief Get the GEOSX element type
 * @param[in] cellType the vtk cell type
 * @return the GEOSX element type
 */
ElementType getElementType( int cellType )
{
  ElementType elementType = ElementType::Polyhedron;
  switch( cellType )
  {
    case VTK_HEXAHEDRON:
      elementType = ElementType::Hexahedron;
      break;
    case VTK_TETRA:
      elementType = ElementType::Tetrahedron;
      break;
    case VTK_WEDGE:
      elementType = ElementType::Prism;
      break;
    case VTK_PYRAMID:
      elementType = ElementType::Pyramid;
      break;
    default:
      GEOSX_ERROR( cellType << " is not a recognized cell type to be used with the VTKMeshGenerator" );
  }
  return elementType;
}


/**
 * @brief Write the hexahedron vertices
 * @details The node ordering from VTK differs from the node ordering in GEOSX
 * @param[in,out] cellToVertex list of nodes organized per cells
 * @param[in] cellIds the hex indexes within this region
 * @param[in,out] localToGLobal the local index to global index map
 * @param[in] mesh the vtkUnstructuredGrid that is loaded
 */
void writeHexahedronVertices( CellBlock::NodeMapType & cellToVertex, std::vector< vtkIdType > const & cellIds,
                              arrayView1d< globalIndex > const & localToGlobal, vtkUnstructuredGrid & mesh )
{
  localIndex cellCount = 0;
  vtkIdTypeArray const & globalCellIdsDataArray = getDataArray< vtkIdTypeArray >( *mesh.GetCellData(), "GlobalCellIds" );
  for( vtkIdType c : cellIds )
  {
    vtkCell * currentCell = mesh.GetCell( c );
    cellToVertex[cellCount][0] = currentCell->GetPointId( 0 );
    cellToVertex[cellCount][1] = currentCell->GetPointId( 1 );
    cellToVertex[cellCount][2] = currentCell->GetPointId( 3 );
    cellToVertex[cellCount][3] = currentCell->GetPointId( 2 );
    cellToVertex[cellCount][4] = currentCell->GetPointId( 4 );
    cellToVertex[cellCount][5] = currentCell->GetPointId( 5 );
    cellToVertex[cellCount][6] = currentCell->GetPointId( 7 );
    cellToVertex[cellCount][7] = currentCell->GetPointId( 6 );
    localToGlobal[cellCount] = globalCellIdsDataArray.GetValue( c );
    cellCount++;
  }
}

/**
 * @brief Write the wedge vertices
 * @details The node ordering from VTK differs from the node ordering in GEOSX
 * @param[in,out] cellToVertex list of nodes organized per cells
 * @param[in] cellIds the wedges indexes within this region
 * @param[in,out] localToGLobal the local index to global index map
 * @param[in] mesh the vtkUnstructuredGrid that is loaded
 */
void writeWedgeVertices( CellBlock::NodeMapType & cellToVertex, std::vector< vtkIdType > const & cellIds,
                         arrayView1d< globalIndex > const & localToGlobal, vtkUnstructuredGrid & mesh )
{
  localIndex cellCount = 0;
  vtkIdTypeArray const & globalCellIdsDataArray = getDataArray< vtkIdTypeArray >( *mesh.GetCellData(), "GlobalCellIds" );
  for( vtkIdType c : cellIds )
  {
    vtkCell * currentCell = mesh.GetCell( c );
    cellToVertex[cellCount][0] = currentCell->GetPointId( 0 );
    cellToVertex[cellCount][1] = currentCell->GetPointId( 3 );
    cellToVertex[cellCount][2] = currentCell->GetPointId( 1 );
    cellToVertex[cellCount][3] = currentCell->GetPointId( 4 );
    cellToVertex[cellCount][4] = currentCell->GetPointId( 2 );
    cellToVertex[cellCount][4] = currentCell->GetPointId( 5 );
    localToGlobal[cellCount] = globalCellIdsDataArray.GetValue( c );
    cellCount++;
  }
}

/**
 * @brief Write a CellBlock of a given cell type
 * @param[in] name the name of the cellBlock to be written
 * @param[in] cellIds the cell indexes of cell type \p cellType within this region
 * @param[in] region_id the id of the region
 * @param[in] cellType the vtk cell type for cells of the CellBlock being written
 * @param[in] cellBlockManager the CellBlockManager
 * @param[in] arraysTobeImported the list of arrays to be imported
 * @param[in] mesh the vtkUnstructuredGrid that is loaded
 */
void writeCellBlock( string const & name, std::vector< vtkIdType > const & cellIds, int region_id, int cellType,
                     CellBlockManager & cellBlockManager,
                     std::vector< vtkDataArray * > const & arraysTobeImported,
                     vtkUnstructuredGrid & mesh )
{

  localIndex numCells = cellIds.size();
  assert( region_id >= -1 );
  string const cellBlockName = region_id != -1  ? std::to_string( region_id ) + "_" + name : name;

  CellBlock & cellBlock = cellBlockManager.getGroup( keys::cellBlocks ).getGroup< CellBlock >( cellBlockName );
  int numPointsPerCell = getNumberOfPoints( cellType );
  CellBlock::NodeMapType & cellToVertex = cellBlock.nodeList();
  arrayView1d< globalIndex > const & localToGlobal = cellBlock.localToGlobalMap();

  /// Writing connectivity and Local to Global
  cellToVertex.resize( numCells, numPointsPerCell );
  if( cellType == VTK_HEXAHEDRON ) // Special case for hexahedron because of the ordering
  {
    writeHexahedronVertices( cellToVertex, cellIds, localToGlobal, mesh );
    return;
  }
  if( cellType == VTK_WEDGE ) // Special case for hexahedron because of the ordering
  {
    writeWedgeVertices( cellToVertex, cellIds, localToGlobal, mesh );
    return;
  }
  vtkIdTypeArray const & globalCellIdsDataArray = getDataArray< vtkIdTypeArray >( *mesh.GetCellData(), "GlobalCellIds" );
  localIndex cellCount = 0;
  for( vtkIdType c : cellIds )
  {
    vtkCell * currentCell = mesh.GetCell( c );
    for( localIndex v = 0; v < numPointsPerCell; v++ )
    {
      cellToVertex[cellCount][v] = currentCell->GetPointId( v );
    }
    localToGlobal[cellCount] = globalCellIdsDataArray.GetValue( c );
    cellCount++;
  }
}

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

void importFieldOnCellBlock( string const & name, std::vector< vtkIdType > const & cellIds, int region_id,
                             ElementRegionManager & elemManager,
                             std::vector< vtkDataArray * > const & arraysTobeImported )
{
  assert( region_id >= -1 );
  string const cellBlockName = region_id != -1  ? std::to_string( region_id ) + "_" + name : name;
  elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
  {
    if( subRegion.getName() == cellBlockName )
    {
      /// Writing properties
      std::unordered_set< string > const materialWrapperNames = getMaterialWrapperNames( subRegion );
      for( vtkDataArray * vtkArray : arraysTobeImported )
      {
        WrapperBase & wrapper = subRegion.getWrapperBase( vtkArray->GetName() );
        if( materialWrapperNames.count( vtkArray->GetName() ) > 0 && wrapper.numArrayDims() > 1 )
        {
          // Scalar material fields are stored as 2D arrays, vector/tensor are 3D
          GEOSX_LOG_RANK_0( "going in with : " << vtkArray->GetName());
          using ImportTypes = types::ArrayTypes< types::RealTypes, types::DimsRange< 2, 3 > >;
          types::dispatch( ImportTypes{}, wrapper.getTypeId(), true, [&]( auto array )
          {
            using ArrayType = decltype( array );
            Wrapper< ArrayType > & wrapperT = Wrapper< ArrayType >::cast( wrapper );
            auto const view = wrapperT.reference().toView();

            localIndex cellCount = 0;
            for( vtkIdType c : cellIds )
            {
              int comp = 0;
              for( localIndex q = 0; q < view.size( 1 ); ++q )
              {
                // The same value is copied for all quadrature points.
                LvArray::forValuesInSlice( view[cellCount][q], [&]( auto & val )
                {
                  val = vtkArray->GetComponent( c, comp );
                } );
              }
              cellCount++;
              comp++;
            }
          } );
        }
        else
        {
          using ImportTypes = types::ArrayTypes< types::RealTypes, types::DimsRange< 1, 2 > >;
          types::dispatch( ImportTypes{}, wrapper.getTypeId(), true, [&]( auto array )
          {
            using ArrayType = decltype( array );
            Wrapper< ArrayType > & wrapperT = Wrapper< ArrayType >::cast( wrapper );
            auto const view = wrapperT.reference().toView();
            localIndex cellCount = 0;
            for( vtkIdType c : cellIds )
            {
              LvArray::forValuesInSlice( view[cellCount], [&]( auto & val )
              {
                val = vtkArray->GetComponent( c, 0 );
              } );
              cellCount++;
            }
          } );
        }
      }
    }
  } );
}
void VTKMeshGenerator::importFields( DomainPartition & domain ) const
{

  GEOSX_LOG_RANK_0( "import fields" );
  ElementRegionManager & elemManager = domain.getMeshBody( this->getName() ).getMeshLevel( 0 ).getElemManager();
  auto importFieldsForElementType= [&]( string const & name, std::map< int, std::vector< vtkIdType > > const & cellBlocks ) {
    for( auto & cellBlock : cellBlocks )
    {
      GEOSX_LOG_RANK_0( "importing a cell block" );
      importFieldOnCellBlock( name, cellBlock.second, cellBlock.first, elemManager, m_arraysToBeImported );
    }
  };

  importFieldsForElementType( "hexahedron", m_regionsHex );
  importFieldsForElementType( "tetrahedron", m_regionsTetra );
  importFieldsForElementType( "wedges", m_regionsWedges );
  importFieldsForElementType( "pyramids", m_regionsPyramids );

  std::map< string, string_array > fieldNames;
  for( vtkDataArray * vtkArray : m_arraysToBeImported )
  {
    fieldNames["elems"].emplace_back( vtkArray->GetName() );
    std::cout << "fieldNames" << fieldNames["elems"];
  }

  CommunicationTools::getInstance().synchronizeFields( fieldNames, domain.getMeshBody( this->getName() ).getMeshLevel( 0 ), domain.getNeighbors(), false );
}

void registerCellBlock( string const & name, std::vector< vtkIdType > const & cellIds, int region_id, int cellType,
                        CellBlockManager & cellBlockManager )
{
  localIndex numCells = cellIds.size();
  assert( region_id >= -1 );
  string const cellBlockName = region_id != -1  ? std::to_string( region_id ) + "_" + name : name;

  CellBlock & cellBlock = cellBlockManager.getGroup( keys::cellBlocks ).registerGroup< CellBlock >( cellBlockName );
  cellBlock.setElementType( getElementType( cellType ) );
  cellBlock.resize( numCells );
}

/**
 * @brief Write all the cell blocks
 * @param[in] domain the domain in which the cell blocks will be written
 * @param[in] regionsHex map from region index to the hexahedron indexes in this region
 * @param[in] regionsTetra map from region index to the tetra indexes in this region
 * @param[in] regionsWedges map from region index to the wedges indexes in this region
 * @param[in] regionsPyramids map from region index to the pyramids indexes in this region
 * @param[in] regions a set containing all the region indexes
 * @param[in] arraysToBeImported a vector containing all the vtkDataArray that can be imported
 * @param[in] mesh the vtkUnstructuredGrid that is loaded
 */
void writeCellBlocks( DomainPartition & domain,
                      std::map< int, std::vector< vtkIdType > > & regionsHex, std::map< int, std::vector< vtkIdType > > & regionsTetra,
                      std::map< int, std::vector< vtkIdType > > & regionsWedges, std::map< int, std::vector< vtkIdType > > & regionsPyramids,
                      std::vector< vtkDataArray * > & arraysToBeImported,
                      vtkUnstructuredGrid & mesh )
{
  CellBlockManager & cellBlockManager = domain.getGroup< CellBlockManager >( keys::cellManager );

  auto writeCellBlockForElementType= [&]( string const & name, std::map< int, std::vector< vtkIdType > > & cellBlocks, int vtkType ) {
    for( auto & cellBlock : cellBlocks )
    {
      registerCellBlock( name, cellBlock.second, cellBlock.first, vtkType, cellBlockManager );
      writeCellBlock( name, cellBlock.second, cellBlock.first, vtkType, cellBlockManager, arraysToBeImported, mesh );
    }
  };

  writeCellBlockForElementType( "hexahedron", regionsHex, VTK_HEXAHEDRON );
  writeCellBlockForElementType( "tetrahedron", regionsTetra, VTK_TETRA );
  writeCellBlockForElementType( "wedges", regionsWedges, VTK_WEDGE );
  writeCellBlockForElementType( "pyramids", regionsPyramids, VTK_PYRAMID );

}

/**
 * @param[in] nodeManager the NodeManager of the domain in which the poiints will be copied.
 * @param[in] allSurfaces the surfaces id to be imported
 * @param[in] surfaces map from the surfaces index to the list of cells in this surface in this rank
 * @param[in] mesh the vtkUnstructuredGrid that is loaded
 */
void writeSurfaces( NodeManager & nodeManager, std::vector< int > const & allSurfaces,
                    std::map< int, std::vector< vtkIdType > > surfaces, vtkUnstructuredGrid & mesh )
{
  Group & nodeSets = nodeManager.sets();
  /// Register all surfaces (even those which are empty in this rank)
  for( int surface : allSurfaces )
  {
    nodeSets.registerWrapper< SortedArray< localIndex > >( std::to_string( surface ) ).reference();
  }
  for( auto surface : surfaces )
  {
    SortedArray< localIndex > & curNodeSet  = nodeSets.getWrapper< SortedArray< localIndex > >( std::to_string( surface.first ) ).reference();
    for( vtkIdType c : surface.second )
    {
      vtkCell * currentCell = mesh.GetCell( c );
      for( localIndex v = 0; v < currentCell->GetNumberOfPoints(); v++ )
      {
        curNodeSet.insert( currentCell->GetPointId( v ) );
      }
    }
  }
}
void VTKMeshGenerator::generateMesh( DomainPartition & domain )
{
  int mpiSize = MpiWrapper::commSize();
  GEOSX_ERROR_IF( (mpiSize & (mpiSize - 1)) != 0, "MPI size is not a power of 2. Can't be used with the VTKMeshGenerator" );
  vtkSmartPointer< vtkMultiProcessController > controller = getVTKController();
  vtkMultiProcessController::SetGlobalController( controller );

  vtkSmartPointer< vtkUnstructuredGrid > loadedMesh = loadVTKMesh( m_filePath );
  MpiWrapper::barrier();

  m_vtkMesh =redistributeMesh( *loadedMesh );
  MpiWrapper::barrier();

  GEOSX_LOG_RANK_0( "Successfully redistributed the mesh on the " << controller->GetNumberOfProcesses() << " MPI processes" );

  Group & meshBodies = domain.getMeshBodies();
  MeshBody & meshBody = meshBodies.registerGroup< MeshBody >( this->getName() );

  MeshLevel & meshLevel0 = meshBody.registerGroup< MeshLevel >( string( "Level0" ));
  NodeManager & nodeManager = meshLevel0.getNodeManager();

  double globalLength = writeMeshNodes( nodeManager, *m_vtkMesh );
  meshBody.setGlobalLengthScale( globalLength );

  domain.getMetisNeighborList() = computePotentialNeighborLists( m_vtkMesh->GetBounds(), globalLength );

  std::map< int, std::vector< vtkIdType > > surfaces; //local to this rank
  std::vector< int > allSurfaces; // global to the whole mesh
  std::vector< string > allCellBlockNames;

  countCellsAndFaces( m_regionsHex, m_regionsTetra, m_regionsWedges, m_regionsPyramids, surfaces, allSurfaces, *m_vtkMesh );

  m_arraysToBeImported = findArrayToBeImported( *m_vtkMesh );

  writeCellBlocks( domain, m_regionsHex, m_regionsTetra, m_regionsWedges, m_regionsPyramids, m_arraysToBeImported, *m_vtkMesh );

  importFields( domain );

  writeSurfaces( nodeManager, allSurfaces, surfaces, *m_vtkMesh );

  MpiWrapper::barrier();
  GEOSX_LOG_RANK_0( "Mesh was loaded successfully" );
}
REGISTER_CATALOG_ENTRY( MeshGeneratorBase, VTKMeshGenerator, string const &, Group * const )
}
