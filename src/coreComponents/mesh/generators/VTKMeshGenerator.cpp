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
#include "mesh/MeshBody.hpp"

#include "common/MpiWrapper.hpp"

#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLPUnstructuredGridReader.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkRedistributeDataSetFilter.h>
#include <vtkGenerateGlobalIds.h>

#include <numeric>
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
{
}

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
bool checkIntersection(double box1[6], double box2[6])
{
    //GEOSX_LOG_RANK(box2[0] << " and " box1[0] and)
    return box1[0] < box2[1] && box1[1]> box2[0] && box1[2] < box2[3] && box1[3]> box2[2] && box1[4] < box2[5] && box1[5]> box2[4];
}

/**
 * @brief Return a VTK controller for multiprocessing.
 */
vtkSmartPointer<vtkMultiProcessController> getVTKController()
{
  #ifdef GEOSX_USE_MPI
    vtkNew<vtkMPIController> controller;
    vtkMPICommunicatorOpaqueComm vtkGeosxComm(&MPI_COMM_GEOSX);
    vtkNew<vtkMPICommunicator> communicator;
    communicator->InitializeExternal( &vtkGeosxComm);
    controller->SetCommunicator(communicator);
  #else
    vtkNew<vtkDummyController> controller;
  #endif
  return controller;
}

/**
 * @brief Load the VTK file into the VTK data structure
 * @param[in] filePath the Path of the file to load
 */
vtkSmartPointer<vtkUnstructuredGrid> loadVTKMesh( Path const & filePath )
{
  string_array filePathTokenized = stringutilities::tokenize( filePath, "."); //TODO maybe code a method in Path to get the file extension?
  string extension = filePathTokenized[filePathTokenized.size() - 1];
  vtkSmartPointer< vtkUnstructuredGrid > loadedMesh = vtkSmartPointer< vtkUnstructuredGrid >::New();

  auto read = [&](auto vtkUgReader) { // auto can be multiple types in the same function
    vtkUgReader->SetFileName( filePath.c_str() );
    vtkUgReader->Update();
    loadedMesh = vtkUgReader->GetOutput();
  };

  if( extension == "pvtu" )
  {
    read(vtkSmartPointer<vtkXMLPUnstructuredGridReader>::New());
  }
  else
  {
    if( MpiWrapper::commRank() == 0 )
    {
    if( extension == "vtk")
      {
        read(vtkSmartPointer<vtkUnstructuredGridReader>::New());
      }
      else if( extension == "vtu")
      {
        read(vtkSmartPointer<vtkXMLUnstructuredGridReader>::New());
      }
      else
      {
        GEOSX_ERROR( extension << " is not a recognized extension for using the VTK reader with GEOSX. Please use .vtk or .vtu" );
      }
    } 
  }
  GEOSX_LOG_RANK( loadedMesh->GetNumberOfCells());

  return loadedMesh;
}

/**
 * @brief Redistribute the mesh among the available MPI ranks
 * @details this method will also generate global ids for points and cells in the VTK Mesh
 * @param[in] loadedMesh the mesh that was loaded on one or several MPI ranks
 */
vtkSmartPointer< vtkUnstructuredGrid > redistributeMesh( vtkUnstructuredGrid & loadedMesh)
{
  // Redistribute data all over the available ranks
  vtkNew<vtkRedistributeDataSetFilter> rdsf;
  rdsf->SetInputDataObject( &loadedMesh );
  rdsf->SetNumberOfPartitions( MpiWrapper::commSize() );
  rdsf->GenerateGlobalCellIdsOn();
  rdsf->Update();
  MpiWrapper::barrier();

  // Generate global IDs for vertices and cells
  vtkNew<vtkGenerateGlobalIds> generator;
  generator->SetInputDataObject(vtkUnstructuredGrid::SafeDownCast( rdsf->GetOutputDataObject(0)));
  generator->Update();
  vtkSmartPointer< vtkUnstructuredGrid > mesh = vtkUnstructuredGrid::SafeDownCast( generator->GetOutputDataObject(0));
  return mesh;
}

/**
 * @brief Copy the VTK mesh nodes into the nodeManager of GEOSX
 * @param[in] nodeManager the NodeManager of the domain in which the poiints will be copied.
 * @return the global length of the mesh (diagonal of the bounding box)
 */
double writeMeshNodes(NodeManager & nodeManager, vtkUnstructuredGrid & mesh)
{
  nodeManager.resize( mesh.GetNumberOfPoints() );
  arrayView1d< globalIndex > const & nodeLocalToGlobal = nodeManager.localToGlobalMap();

  // Writing the points
  GEOSX_ERROR_IF( mesh.GetNumberOfPoints() == 0, "Mesh is empty, aborting");

  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();
  Group & nodeSets = nodeManager.sets();
  SortedArray< localIndex > & allNodes  = nodeSets.registerWrapper< SortedArray< localIndex > >( "all" ).reference();
  real64 xMax[3] = { std::numeric_limits< real64 >::min(), std::numeric_limits< real64 >::min(), std::numeric_limits< real64 >::min() };
  real64 xMin[3] = { std::numeric_limits< real64 >::max(), std::numeric_limits< real64 >::max() , std::numeric_limits< real64 >::max()  };

  vtkIdTypeArray const & globalPointIdDataArray = getDataArray< vtkIdTypeArray >( *mesh.GetPointData(), "GlobalPointIds" );

  for( vtkIdType v = 0; v < mesh.GetNumberOfPoints(); v++)
  {
    double * point = mesh.GetPoint( v );
    nodeLocalToGlobal[v] = globalPointIdDataArray.GetValue( v );
    for( integer i = 0; i < 3; i++)
    {
      X(v,i) = point[i];
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
  for(unsigned int i = 0; i < 6;i++)
  {
    rankBoundsHalo[i] = rankBounds[i];
    rankBoundsHalo[i] += i%2?-haloSize:haloSize;
  }

  std::vector<double> allBounds(6 * MpiWrapper::commSize());
  MpiWrapper::allgather<double, double>(rankBounds, 6, allBounds.data(),6,MPI_COMM_GEOSX);

  //Checking the intersection
  std::set< int > rankNeighbor;
  for(int i = 0; i < MpiWrapper::commSize(); i++)
  {
    if( i == MpiWrapper::commRank())
    {
      continue;
    }
    std::vector< double >::const_iterator first = allBounds.begin() + 6*i;
    std::vector< double >::const_iterator last = allBounds.begin() + 6*(i+1);
    std::vector< double > curBound(first,last);
    if(checkIntersection(rankBoundsHalo, curBound.data()))
    {
      rankNeighbor.insert(i);
    }
  }
  return rankNeighbor;
}

/**
 * @brief This method is used to preprocess the the VTK mesh and count the number of cells, facets, regions
 * and surfaces over the current MPI rank
 * @param[out] numHex number of hexahedra
 * @param[out] numTet number of tetra
 * @param[out] numWedge number of wedges
 * @param[out] numPyr number of pyramids
 * @param[out] regionsHex map from region index to the number of hexahedron in this region
 * @param[out] regionsTetra map from region index to the number of tetra in this region
 * @param[out] regionsWedges map from region index to the number of wedges in this region
 * @param[out] regionsPyramids map from region index to the number of pyramids in this region
 * @param[out] regions a set containing all the region indexes
 * @param[out] surfaces a set containing all the surface indexes, from this MPI rank
 * @param[out] allSurfaces a vector containing all the surfaces among all the MPI rank
 */
void countCellsAndFaces( localIndex & numHex, localIndex & numTet, localIndex & numWedge, localIndex & numPyr,
                         std::map<int,localIndex> & regionsHex, std::map<int,localIndex> & regionsTetra,
                         std::map<int,localIndex> & regionsWedges, std::map<int,localIndex> & regionsPyramids,
                         std::set< int > & regions, std::set< int > & surfaces, std::vector< int > & allSurfaces,
                         vtkUnstructuredGrid & mesh)
{
  vtkIntArray const * const attributeDataArray = getDataArrayOptional< vtkIntArray >( *mesh.GetCellData(), "attribute" );

  auto countCell = [&](localIndex c, localIndex & numElem, std::map<int,localIndex> & regionsElem ) { // auto can be multiple types in the same function
    numElem++;
    if (attributeDataArray != nullptr)
    {
      int region = attributeDataArray->GetValue(c);
      GEOSX_LOG_RANK_0_IF( region < 0,"Attribute value " << region << " not supported, please use value >=0 to describe regions");
      //regions.insert(region);
      ++regionsElem[region];
    }
  };


  for( vtkIdType c = 0; c < mesh.GetNumberOfCells(); c++)
  {
    if( mesh.GetCellType(c) == VTK_HEXAHEDRON )
    {
      countCell(c, numHex, regionsHex);
    }
    else if( mesh.GetCellType(c) == VTK_TETRA )
    {
      countCell(c, numTet, regionsTetra);
    }
    else if( mesh.GetCellType(c) == VTK_WEDGE )
    {
      countCell(c, numWedge, regionsWedges);
    }
    else if( mesh.GetCellType(c) == VTK_PYRAMID )
    {
      countCell(c, numPyr, regionsPyramids);
    }
    else if( mesh.GetCellType(c) == VTK_TRIANGLE || mesh.GetCellType(c) == VTK_QUAD )
    {
      if (attributeDataArray != nullptr)
      {
        surfaces.insert(attributeDataArray->GetValue(c));
      }
    }
  }

  // Fill the region set
  auto fillRegionSet = [&](std::map<int,localIndex> & regionsElem ) {
  for( auto const & region : regionsElem )
  {
    regions.insert(region.first);
  }
  };
  fillRegionSet(regionsHex);
  fillRegionSet(regionsTetra);
  fillRegionSet(regionsPyramids);
  fillRegionSet(regionsWedges);


  // Communicate all the surfaces attributes to the ranks
  const std::vector< int > surfaceVector(surfaces.cbegin(), surfaces.cend());
  array1d< int > surfaceSizes( MpiWrapper::commSize());
  MpiWrapper::allGather( LvArray::integerConversion<int>(surfaceVector.size()), surfaceSizes, MPI_COMM_GEOSX );
  int const totalNbSurfaceId = std::accumulate( surfaceSizes.begin(), surfaceSizes.end(), 0 );
  GEOSX_LOG_RANK("total nb surface id " << totalNbSurfaceId);
  allSurfaces.resize(totalNbSurfaceId);
  std::vector< int > displacements( MpiWrapper::commSize(),0);
  std::partial_sum( surfaceSizes.begin(), surfaceSizes.end() - 1, displacements.begin() + 1 );
  MpiWrapper::allgatherv(surfaceVector.data(),  surfaceVector.size(), allSurfaces.data(),surfaceSizes.data(), displacements.data(), MPI_COMM_GEOSX);

  std::sort(allSurfaces.begin(), allSurfaces.end());
  auto last = std::unique( allSurfaces.begin(), allSurfaces.end());
  allSurfaces.erase(last, allSurfaces.end());
}

/**
 * @brief Find the properties to be imported
 * @details all the float and double vtkArrays will be imported
 * @return a vector containing all the vtkDataArray that can be imported
 */
std::vector< vtkDataArray * > findArrayToBeImported(vtkUnstructuredGrid & mesh)
{
  std::vector< vtkDataArray * > arrayToBeImported;
  vtkCellData * cellData = mesh.GetCellData();
  for( int i = 0; i < cellData->GetNumberOfArrays(); i++)
  {
    vtkAbstractArray * curArray = cellData->GetAbstractArray(i);
    // We support only float and double type. Good assumption?
    int dataType = curArray->GetDataType();
    if( !curArray->IsNumeric() )
    {
      GEOSX_LOG_RANK_0( curArray->GetName() << " is not a numeric array");
      continue;
    }
    if( dataType != VTK_FLOAT && dataType != VTK_DOUBLE )
    {
      GEOSX_LOG_RANK_0( "Underlying data Type " << curArray->GetDataTypeAsString() << " for array "
                        << curArray->GetName() << " is not supported by GEOSX and will be ignored "
                        << "(Only double and float are supported)");
      continue;
    }
    GEOSX_LOG_RANK_0("Importing array " << curArray->GetName() );
    arrayToBeImported.push_back( vtkDataArray::SafeDownCast(curArray) );
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
  switch( cellType)
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
      GEOSX_ERROR( cellType << " is not a recognized cell type to be used with the VTKMeshGenerator");
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
  switch( cellType)
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
      GEOSX_ERROR( cellType << " is not a recognized cell type to be used with the VTKMeshGenerator");
  }
  return elementType;
}


/**
 * @brief Write the hexahedron vertices
 * @details The node ordering from VTK differs from the node ordering in GEOSX
 * @param[in,out] cellToVertex list of nodes organized per cells
 */
void writeHexahedronVertices( CellBlock::NodeMapType & cellToVertex, int region_id,
                              arrayView1d< globalIndex > const & localToGlobal, vtkUnstructuredGrid & mesh )
{
  localIndex cellCount = 0;
  vtkIntArray const * const attributeDataArray = getDataArrayOptional< vtkIntArray >( *mesh.GetCellData(), "attribute" );
  vtkIdTypeArray const & globalCellIdsDataArray = getDataArray< vtkIdTypeArray >( *mesh.GetCellData(), "GlobalCellIds" );
  for( vtkIdType c = 0; c < mesh.GetNumberOfCells(); c++)
  {
    vtkCell * currentCell = mesh.GetCell(c);
    if( currentCell->GetCellType() == VTK_HEXAHEDRON )
    {
      if( attributeDataArray != nullptr && attributeDataArray->GetValue(c) != region_id) continue;
		  cellToVertex[cellCount][0] = currentCell->GetPointId(0);
		  cellToVertex[cellCount][1] = currentCell->GetPointId(1);
		  cellToVertex[cellCount][2] = currentCell->GetPointId(3);
		  cellToVertex[cellCount][3] = currentCell->GetPointId(2);
		  cellToVertex[cellCount][4] = currentCell->GetPointId(4);
		  cellToVertex[cellCount][5] = currentCell->GetPointId(5);
		  cellToVertex[cellCount][6] = currentCell->GetPointId(7);
		  cellToVertex[cellCount][7] = currentCell->GetPointId(6);
      localToGlobal[cellCount] = globalCellIdsDataArray.GetValue(c);
      cellCount++;
    }
  }
}

/**
 * @brief Write the wedge vertices
 * @details The node ordering from VTK differs from the node ordering in GEOSX
 * @param[in,out] cellToVertex list of nodes organized per cells
 */
void writeWedgeVertices( CellBlock::NodeMapType & cellToVertex, int region_id,
                         arrayView1d< globalIndex > const & localToGlobal, vtkUnstructuredGrid & mesh )
{
  localIndex cellCount = 0;
  vtkIntArray const * const attributeDataArray = getDataArrayOptional< vtkIntArray >( *mesh.GetCellData(), "attribute" );
  vtkIdTypeArray const & globalCellIdsDataArray = getDataArray< vtkIdTypeArray >( *mesh.GetCellData(), "GlobalCellIds" );
  for( vtkIdType c = 0; c < mesh.GetNumberOfCells(); c++)
  {
    vtkCell * currentCell = mesh.GetCell(c);
    if( currentCell->GetCellType() == VTK_WEDGE )
    {
      if( attributeDataArray != nullptr && attributeDataArray->GetValue(c) != region_id) continue;
		  cellToVertex[cellCount][0] = currentCell->GetPointId(0);
		  cellToVertex[cellCount][1] = currentCell->GetPointId(3);
		  cellToVertex[cellCount][2] = currentCell->GetPointId(1);
		  cellToVertex[cellCount][3] = currentCell->GetPointId(4);
		  cellToVertex[cellCount][4] = currentCell->GetPointId(2);
		  cellToVertex[cellCount][4] = currentCell->GetPointId(5);
      localToGlobal[cellCount] = globalCellIdsDataArray.GetValue(c);
      cellCount++;
    }
  }
}

/**
 * @brief Write a CellBlock of a given cell type
 * @param[in] name the name of the cellBlock to be written
 * @param[in] numcells number of cells the CellBlock will contain
 * @param[in] region_id the id of the region
 * @param[in] cellType the vtk cell type for cells of the CellBlock being written
 * @param[in] cellBlockManager the CellBlockManager
 * @param[in] arraysTobeImported the list of arrays to be imported
 */
void writeCellBlock( string const & name, localIndex numCells, int region_id, int cellType,
                                       CellBlockManager & cellBlockManager,
                                       std::vector< vtkDataArray * > const & arraysTobeImported,
                                       vtkUnstructuredGrid & mesh )
{

  if( numCells > 0)
  {
    assert(region_id >= -1 );
    string const cellBlockName = region_id != -1  ? std::to_string(region_id) + "_" + name : name;

    GEOSX_LOG_RANK_0( "Writing " << numCells << " cells to CellBlock " << cellBlockName);
    CellBlock & cellBlock = cellBlockManager.getGroup( keys::cellBlocks ).registerGroup< CellBlock >( cellBlockName );
    int numPointsPerCell = getNumberOfPoints( cellType);
    cellBlock.setElementType( getElementType(cellType) );
    CellBlock::NodeMapType & cellToVertex = cellBlock.nodeList();
    cellBlock.resize( numCells );
    arrayView1d< globalIndex > const & localToGlobal = cellBlock.localToGlobalMap();

    vtkIntArray const * const attributeDataArray = getDataArrayOptional< vtkIntArray >( *mesh.GetCellData(), "attribute" );


    /// Writing properties
    for( vtkDataArray * array : arraysTobeImported )
    {
      int dimension = array->GetNumberOfComponents();
      if( dimension == 1 )
      {
        real64_array & property = cellBlock.addProperty< real64_array >( array->GetName() );
        localIndex cellCount = 0;
        for( vtkIdType c = 0; c < mesh.GetNumberOfCells(); c++)
        {
          vtkCell * currentCell = mesh.GetCell(c);
          if( currentCell->GetCellType() == cellType )
          {
            if( attributeDataArray != nullptr && attributeDataArray->GetValue(c) != region_id) continue;
            property[cellCount++] = array->GetTuple1(c);
          }
        }
      }
      else
      {
        array2d< real64 > & property = cellBlock.addProperty< array2d< real64 > >( array->GetName() );
        property.resizeDimension< 1 >( dimension );
        localIndex cellCount = 0;
        for( vtkIdType c = 0; c < mesh.GetNumberOfCells(); c++)
        {
          vtkCell * currentCell = mesh.GetCell(c);
          if( currentCell->GetCellType() == cellType )
          {
            if( attributeDataArray != nullptr && attributeDataArray->GetValue(c) != region_id ) continue;
            for(int i = 0; i < dimension; i++ )
            {
              property(cellCount,i) = array->GetComponent(c,i);
            }
            cellCount++;
          }
        }
      }
    }

    /// Writing connectivity and Local to Global
    cellToVertex.resize( numCells, numPointsPerCell );
    if( cellType == VTK_HEXAHEDRON ) // Special case for hexahedron because of the ordering
    {
      writeHexahedronVertices( cellToVertex, region_id, localToGlobal, mesh );
      return;
    }
    if( cellType == VTK_WEDGE ) // Special case for hexahedron because of the ordering
    {
      writeWedgeVertices( cellToVertex, region_id, localToGlobal, mesh );
      return;
    }
    vtkIdTypeArray const & globalCellIdsDataArray = getDataArray< vtkIdTypeArray >( *mesh.GetCellData(), "GlobalCellIds" );
    localIndex cellCount = 0;
    for( vtkIdType c = 0; c < mesh.GetNumberOfCells(); c++)
    {
      vtkCell * currentCell = mesh.GetCell(c);
      if( currentCell->GetCellType() == cellType )
      {
        if( attributeDataArray != nullptr && attributeDataArray->GetValue(c) != region_id) continue;
        for( localIndex v = 0; v < numPointsPerCell; v++)
        {
          cellToVertex[cellCount][v] = currentCell->GetPointId(v);
        }
        localToGlobal[cellCount] = globalCellIdsDataArray.GetValue(c);
        cellCount++;
      }
    }
  }
}


/**
 * @brief Write all the cell blocks
 * @param[in] domain the domain in which the cell blocks will be written
 * @param[in] numHex number of hexahedra
 * @param[in] numTet number of tetra
 * @param[in] numWedge number of wedges
 * @param[in] numPyr number of pyramids
 * @param[in] regionsHex map from region index to the number of hexahedron in this region
 * @param[in] regionsTetra map from region index to the number of tetra in this region
 * @param[in] regionsWedges map from region index to the number of wedges in this region
 * @param[in] regionsPyramids map from region index to the number of pyramids in this region
 * @param[in] regions a set containing all the region indexes
 * @param[in] arraysToBeImported a vector containing all the vtkDataArray that can be imported
 */
void writeCellBlocks( DomainPartition& domain, localIndex numHex, localIndex numTet, localIndex numWedge, localIndex numPyr,
                      std::map<int,localIndex> & regionsHex, std::map<int,localIndex> & regionsTetra,
                      std::map<int,localIndex> & regionsWedges, std::map<int,localIndex> & regionsPyramids,
                      std::set< int > & regions, std::vector< vtkDataArray * > & arraysToBeImported,
                      vtkUnstructuredGrid & mesh)
{
  CellBlockManager & cellBlockManager = domain.getGroup< CellBlockManager >( keys::cellManager );

  if( !regions.empty() )
  {
    for (int region : regions)
    {
      std::map< int, localIndex>::iterator it;
      it = regionsHex.find(region);
      if( it != regionsHex.end())
      {
        writeCellBlock("hexahedron", it->second, it->first, VTK_HEXAHEDRON, cellBlockManager, arraysToBeImported, mesh);
      }
      it = regionsTetra.find(region);
      if( it != regionsTetra.end())
      {
        writeCellBlock("tetrahedron", it->second, it->first, VTK_TETRA, cellBlockManager, arraysToBeImported, mesh);
      }
      it = regionsWedges.find(region);
      if( it != regionsWedges.end())
      {
        writeCellBlock("wedges", it->second, it->first, VTK_WEDGE, cellBlockManager, arraysToBeImported,mesh);
      }
      it = regionsPyramids.find(region);
      if( it != regionsPyramids.end())
        writeCellBlock("pyramids", it->second, it->first, VTK_PYRAMID, cellBlockManager, arraysToBeImported, mesh);
    }
  }
  else
  {
      writeCellBlock("hexahedron", numHex, -1, VTK_HEXAHEDRON, cellBlockManager, arraysToBeImported, mesh);
      writeCellBlock("tetrahedron", numTet, -1, VTK_TETRA, cellBlockManager, arraysToBeImported, mesh);
      writeCellBlock("wedges", numWedge, -1, VTK_WEDGE, cellBlockManager, arraysToBeImported, mesh);
      writeCellBlock("pyramids", numPyr, -1, VTK_PYRAMID, cellBlockManager, arraysToBeImported, mesh);
  }
}


/**
 * @param[in] nodeManager the NodeManager of the domain in which the poiints will be copied.
 * @param[in] allSurfaces the surfaces id to be imported
 */
void writeSurfaces( NodeManager & nodeManager, std::vector<int> const & allSurfaces, vtkUnstructuredGrid & mesh )
{
  Group & nodeSets = nodeManager.sets();
  vtkIntArray const * const attributeDataArray = getDataArrayOptional< vtkIntArray >( *mesh.GetCellData(), "attribute" );
  for( int surface : allSurfaces)
  {
    SortedArray< localIndex > & curNodeSet  = nodeSets.registerWrapper< SortedArray< localIndex > >( std::to_string( surface) ).reference();
    for( vtkIdType c = 0; c < mesh.GetNumberOfCells(); c++)
    {
      vtkCell * currentCell = mesh.GetCell(c);
      if( mesh.GetCellType(c) == VTK_TRIANGLE ||  mesh.GetCellType(c) == VTK_QUAD)
      {
        if( attributeDataArray->GetValue(c) == surface )
        {
          for( localIndex v = 0; v < currentCell->GetNumberOfPoints(); v++)
          {
            curNodeSet.insert(currentCell->GetPointId(v) );
          }
        } 
      }
    }
  }
}
void VTKMeshGenerator::generateMesh( DomainPartition & domain )
{
  vtkSmartPointer<vtkMultiProcessController> controller = getVTKController();
  vtkMultiProcessController::SetGlobalController(controller);

  vtkSmartPointer<vtkUnstructuredGrid> loadedMesh = loadVTKMesh(m_filePath);
  MpiWrapper::barrier();

  m_vtkMesh =redistributeMesh(*loadedMesh);
  MpiWrapper::barrier();

  GEOSX_LOG_RANK_0(" Successfully redistributed the mesh on the " << controller->GetNumberOfProcesses() << " MPI processes");

  Group & meshBodies = domain.getMeshBodies();
  MeshBody & meshBody = meshBodies.registerGroup< MeshBody >( this->getName() );

  MeshLevel & meshLevel0 = meshBody.registerGroup< MeshLevel >( string( "Level0" ));
  NodeManager & nodeManager = meshLevel0.getNodeManager();

  double globalLength = writeMeshNodes( nodeManager, *m_vtkMesh );
  meshBody.setGlobalLengthScale( globalLength);

  domain.getMetisNeighborList() = computePotentialNeighborLists( m_vtkMesh->GetBounds(), globalLength );


  localIndex numHex =0;
  localIndex numTet = 0;
  localIndex numWedge = 0;
  localIndex numPyr = 0;

  // region id link to their number of cells
  std::map<int,localIndex> regionsHex;
  std::map<int,localIndex> regionsTetra;
  std::map<int,localIndex> regionsWedges;
  std::map<int,localIndex> regionsPyramids;

  // surface and greionids
  std::set<int> regions;
  std::set<int> surfaces; //local to this rank
  std::vector<int> allSurfaces; // global to the whole mesh

  countCellsAndFaces( numHex, numTet, numWedge, numPyr, regionsHex, regionsTetra, regionsWedges, regionsPyramids, regions, surfaces, allSurfaces, *m_vtkMesh);

  std::vector< vtkDataArray * > arraysToBeImported = findArrayToBeImported(*m_vtkMesh);

  writeCellBlocks( domain, numHex, numTet, numWedge, numPyr,regionsHex, regionsTetra, regionsWedges, regionsPyramids, regions, arraysToBeImported, *m_vtkMesh );

  writeSurfaces( nodeManager, allSurfaces, *m_vtkMesh);

  GEOSX_LOG_RANK_0("Mesh was loaded successfully");
}
REGISTER_CATALOG_ENTRY( MeshGeneratorBase, VTKMeshGenerator, string const &, Group * const )
}
