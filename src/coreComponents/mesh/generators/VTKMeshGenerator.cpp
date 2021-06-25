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

#include "Elements/Element.hpp"
#include "MeshDataWriters/Variable.hpp"
#include "mesh/DomainPartition.hpp"

#include <math.h>

#include "mesh/mpiCommunications/PartitionBase.hpp"
#include "mesh/mpiCommunications/SpatialPartition.hpp"
#include "common/MpiWrapper.hpp"
#include "Mesh/MeshFactory.hpp"

#include "MeshDataWriters/MeshParts.hpp"

#include "mesh/MeshBody.hpp"

#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLPUnstructuredGridReader.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkRedistributeDataSetFilter.h>
#include <vtkGenerateGlobalIds.h>
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

VTKMeshGenerator::~VTKMeshGenerator()
{}

void VTKMeshGenerator::postProcessInput()
{
  GEOSX_LOG_RANK_0("debut");
}

Group * VTKMeshGenerator::createChild( string const & GEOSX_UNUSED_PARAM( childKey ), string const & GEOSX_UNUSED_PARAM( childName ) )
{
  return nullptr;
}

bool CheckIntersection(double box1[6], double box2[6])
{
    //GEOSX_LOG_RANK(box2[0] << " and " box1[0] and)
    return box1[0] < box2[1] && box1[1]> box2[0] && box1[2] < box2[3] && box1[3]> box2[2] && box1[4] < box2[5] && box1[5]> box2[4];
}

void VTKMeshGenerator::generateMesh( DomainPartition & domain )
{
  vtkSmartPointer<vtkMultiProcessController> controller = GetVTKController();
  vtkMultiProcessController::SetGlobalController(controller);

  vtkSmartPointer<vtkUnstructuredGrid> loadedMesh = LoadVTKMesh();
  MpiWrapper::barrier();

  RedistributeMesh(loadedMesh);
  MpiWrapper::barrier();

  GEOSX_LOG_RANK_0(" Successfully redistributed the mesh on the " << controller->GetNumberOfProcesses() << " MPI processes");

  Group & meshBodies = domain.getGroup( string( "MeshBodies" ));
  MeshBody & meshBody = meshBodies.registerGroup< MeshBody >( this->getName() );

  MeshLevel & meshLevel0 = meshBody.registerGroup< MeshLevel >( string( "Level0" ));
  NodeManager & nodeManager = meshLevel0.getNodeManager();

  double globalLength = WriteMeshNodes( nodeManager );
  meshBody.setGlobalLengthScale( globalLength);

  ComputePotentialNeighborLists( domain );

  localIndex nbHex =0;
  localIndex nbTet = 0;
  localIndex nbWedge = 0;
  localIndex nbPyr = 0;

  // region id link to their number of cells
  std::map<int,localIndex> regions_hex;
  std::map<int,localIndex> regions_tetra;
  std::map<int,localIndex> regions_wedges;
  std::map<int,localIndex> regions_pyramids;

  // surface and greionids
  std::set<int> regions;
  std::set<int> surfaces; //local to this rank
  std::vector<int> allSurfaces; // global to the whole mesh

  CountCellsAndFaces( nbHex, nbTet, nbWedge, nbPyr, regions_hex, regions_tetra, regions_wedges, regions_pyramids, regions, surfaces, allSurfaces);

  std::vector< vtkDataArray * > arraysToBeImported = FindArrayToBeImported();

  WriteCellBlocks( domain, nbHex, nbTet, nbWedge, nbPyr,regions_hex, regions_tetra, regions_wedges, regions_pyramids, regions, arraysToBeImported );

  WriteSurfaces( nodeManager, allSurfaces);

  GEOSX_LOG_RANK_0("Mesh was loaded successfully");
}

vtkSmartPointer<vtkMultiProcessController> VTKMeshGenerator::GetVTKController()
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

vtkSmartPointer<vtkUnstructuredGrid> VTKMeshGenerator::LoadVTKMesh()
{
  string_array filePathTokenized = stringutilities::tokenize( m_filePath, "."); //TODO maybe code a method in Path to get the file extension?
  string extension = filePathTokenized[filePathTokenized.size() - 1];
  vtkSmartPointer< vtkUnstructuredGrid > loadedMesh = vtkSmartPointer< vtkUnstructuredGrid >::New();
  if( extension == "pvtu" )
  {
    vtkSmartPointer<vtkXMLPUnstructuredGridReader> vtkUgReader = vtkSmartPointer<vtkXMLPUnstructuredGridReader>::New();
	  vtkUgReader->SetFileName( m_filePath.c_str() );
    vtkUgReader->Update();
	  loadedMesh = vtkUgReader->GetOutput();
  }
  else
  {
    if( MpiWrapper::commRank() == 0 )
    {
    if( extension == "vtk")
      {
        vtkSmartPointer<vtkUnstructuredGridReader> vtkUgReader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
	      vtkUgReader->SetFileName( m_filePath.c_str() );
        vtkUgReader->Update();
	      loadedMesh = vtkUgReader->GetOutput();
      }
      else if( extension == "vtu")
      {
        vtkSmartPointer<vtkXMLUnstructuredGridReader> vtkUgReader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
	      vtkUgReader->SetFileName( m_filePath.c_str() );
        vtkUgReader->Update();
	      loadedMesh = vtkUgReader->GetOutput();
      }
      else
      {
        GEOSX_ERROR( extension << " is not a recognized extension for using the VTK reader with GEOSX. Please use .vtk or .vtu" );
      }
    } 
    else // Empty mesh for those ranks. It will be filled afther
    {
      loadedMesh = vtkSmartPointer<vtkUnstructuredGrid>::New();
    }
  }
  GEOSX_LOG_RANK( loadedMesh->GetNumberOfCells());

  return loadedMesh;
}

void VTKMeshGenerator::RedistributeMesh( vtkUnstructuredGrid * loadedMesh)
{
  // Redistribute data all over the available ranks
  vtkNew<vtkRedistributeDataSetFilter> rdsf;
  rdsf->SetInputDataObject( loadedMesh );
  rdsf->SetNumberOfPartitions( MpiWrapper::commSize() );
  rdsf->GenerateGlobalCellIdsOn();
  rdsf->Update();
  MpiWrapper::barrier();

  // Generate global IDs for vertices and cells
  vtkNew<vtkGenerateGlobalIds> generator;
  generator->SetInputDataObject(vtkUnstructuredGrid::SafeDownCast( rdsf->GetOutputDataObject(0)));
  generator->Update();
  m_vtkMesh = vtkUnstructuredGrid::SafeDownCast( generator->GetOutputDataObject(0));
  GEOSX_LOG_RANK( m_vtkMesh->GetNumberOfCells());
}

double VTKMeshGenerator::WriteMeshNodes(NodeManager & nodeManager)
{
  nodeManager.resize( m_vtkMesh->GetNumberOfPoints() );
  arrayView1d< globalIndex > const & nodeLocalToGlobal = nodeManager.localToGlobalMap();

  // Writing the points
  GEOSX_ERROR_IF( m_vtkMesh->GetNumberOfPoints() == 0, "Mesh is empty, aborting");
  GEOSX_LOG_RANK_0("Writing " << m_vtkMesh->GetNumberOfPoints() << " vertices");

  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();
  Group & nodeSets = nodeManager.sets();
  SortedArray< localIndex > & allNodes  = nodeSets.registerWrapper< SortedArray< localIndex > >( string( "all" ) ).reference();
  real64 xMax[3] = { std::numeric_limits< real64 >::min() };
  real64 xMin[3] = { std::numeric_limits< real64 >::max() };

  vtkPointData * pointData = m_vtkMesh->GetPointData();
  int globalPointIdIndex = -1;
  vtkAbstractArray * globalPointIdArray = pointData->GetAbstractArray("GlobalPointIds", globalPointIdIndex);
  GEOSX_ERROR_IF(globalPointIdIndex >=0 && globalPointIdArray->GetDataType() != VTK_ID_TYPE, "Array name \"vtkGlobalPointIds\" is not present or do not have " 
                  << "the underlying type has to be an id type");
  vtkIdTypeArray * globalPointIdDataArray =  nullptr;
  if( globalPointIdIndex >= 0)
  {
    globalPointIdDataArray = vtkIdTypeArray::SafeDownCast( globalPointIdArray );
  }
  for( vtkIdType v = 0; v < m_vtkMesh->GetNumberOfPoints(); v++)
  {
    double * point = m_vtkMesh->GetPoint( v );
    nodeLocalToGlobal[v] = globalPointIdDataArray->GetValue( v );
    for( integer i = 0; i < 3; i++)
    {
      X(v,i) = point[i];
      if( X( v, i ) > xMax[i] )
      {
        xMax[i] = X( v, i );
      }
      if( X( v, i ) < xMin[i] )
      {
        xMin[i] = X( v, i );
      }
    }
    allNodes.insert( v );
  }
  LvArray::tensorOps::subtract< 3 >( xMax, xMin );
  return LvArray::tensorOps::l2Norm< 3 >( xMax );
}

void VTKMeshGenerator::ComputePotentialNeighborLists( DomainPartition & domain)
{
  double * rankBounds = m_vtkMesh->GetBounds();
  /// compute the halo 
  for(unsigned int i = 0; i < 6;i++)
  {
    if( i%2)
   {
     rankBounds[i]-=0.01;
   } 
   else
   {
     rankBounds[i]+=0.01;
   }
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
    if(CheckIntersection(rankBounds, curBound.data()))
    {
      domain.getMetisNeighborList().insert(i);
    }
  }
}

vtkIntArray * VTKMeshGenerator::GetAttributeDataArray()
{
  vtkCellData * cellData = m_vtkMesh->GetCellData();
  int attributeIndex = -1;
  vtkAbstractArray * attributeArray = cellData->GetAbstractArray("attribute", attributeIndex);
  GEOSX_ERROR_IF(attributeIndex >=0 && attributeArray->GetDataType() != VTK_INT, "Array name \"attribute\" is reserved to store the regions " 
                  << "used by GEOSX. The underlying type has to be an INT");
  vtkIntArray * attributeDataArray =  nullptr;
  if( attributeIndex >= 0)
  {
    attributeDataArray = vtkIntArray::SafeDownCast( attributeArray );
  }
  return attributeDataArray;
}
  vtkIdTypeArray * VTKMeshGenerator::GetCellGlobalIdDataArray()
  {
    int globalCellIdsIndex = -1;
    vtkCellData * cellData = m_vtkMesh->GetCellData();
    vtkAbstractArray * globalCellIdsArray = cellData->GetAbstractArray("GlobalCellIds", globalCellIdsIndex);
    GEOSX_ERROR_IF(globalCellIdsIndex >=0 && globalCellIdsArray->GetDataType() != VTK_ID_TYPE, "Array name \"GlobalCellIds\" is not present or do not have " 
                  << " the underlying type id type");
    vtkIdTypeArray * globalCellIdsDataArray =  nullptr;
    if( globalCellIdsIndex >= 0)
    {
      globalCellIdsDataArray = vtkIdTypeArray::SafeDownCast( globalCellIdsArray );
    }
    return globalCellIdsDataArray;
  }
void VTKMeshGenerator::CountCellsAndFaces( localIndex & nbHex, localIndex & nbTet, localIndex & nbWedge, localIndex & nbPyr,
                                           std::map<int,localIndex> & regions_hex, std::map<int,localIndex> & regions_tetra,
                                           std::map<int,localIndex> & regions_wedges, std::map<int,localIndex> & regions_pyramids,
                                           std::set< int > & regions, std::set< int > & surfaces, std::vector< int > & allSurfaces)
{
  vtkIntArray * attributeDataArray =  GetAttributeDataArray();

  for( vtkIdType c = 0; c < m_vtkMesh->GetNumberOfCells(); c++)
  {
    if( m_vtkMesh->GetCellType(c) == VTK_HEXAHEDRON )
    {
      nbHex++;
      if (attributeDataArray != nullptr)
      {
        //TODO duplicated code for 4 cell types
        int region = attributeDataArray->GetValue(c);
        GEOSX_LOG_RANK_0_IF( region < 0,"Attribute value " << region << " not supported, please use value >=0 to describe regions");
        regions.insert(region);
        std::map< int, localIndex>::iterator it = regions_hex.find(region);
        if( it == regions_hex.end())
        {
          regions_hex[region] = 1;
        }
        else
        {
          it->second++;
        }
      }
    }
    else if( m_vtkMesh->GetCellType(c) == VTK_TETRA )
    {
      nbTet++;
      if (attributeDataArray != nullptr)
      {
        int region = attributeDataArray->GetValue(c);
        GEOSX_LOG_RANK_0_IF( region < 0,"Attribute value " << region << " not supported, please use value >=0 to describe regions");
        regions.insert(region);
        std::map< int, localIndex>::iterator it = regions_tetra.find(region);
        if( it == regions_tetra.end())
        {
          regions_tetra[region] = 1;
        }
        else
        {
          it->second++;
        }
      }
    }
    else if( m_vtkMesh->GetCellType(c) == VTK_WEDGE )
    {
      nbWedge++;
      if (attributeDataArray != nullptr)
      {
        int region = attributeDataArray->GetValue(c);
        GEOSX_LOG_RANK_0_IF( region < 0,"Attribute value " << region << " not supported, please use value >=0 to describe regions");
        regions.insert(region);
        std::map< int, localIndex>::iterator it = regions_wedges.find(region);
        if( it == regions_wedges.end())
        {
          regions_wedges[region] = 1;
        }
        else
        {
          it->second++;
        }
      }
    }
    else if( m_vtkMesh->GetCellType(c) == VTK_PYRAMID )
    {
      nbPyr++;
      if (attributeDataArray != nullptr)
      {
        int region = attributeDataArray->GetValue(c);
        GEOSX_LOG_RANK_0_IF( region < 0,"Attribute value " << region << " not supported, please use value >=0 to describe regions");
        regions.insert(region);
        std::map< int, localIndex>::iterator it = regions_pyramids.find(region);
        if( it == regions_pyramids.end())
        {
          regions_pyramids[region] = 1;
        }
        else
        {
          it->second++;
        }
      }
    }
    else if( m_vtkMesh->GetCellType(c) == VTK_TRIANGLE || m_vtkMesh->GetCellType(c) == VTK_QUAD )
    {
      if (attributeDataArray != nullptr)
      {
        surfaces.insert(attributeDataArray->GetValue(c));
      }
    }
  }


  // Communicate all the surfaces attributes to the ranks
  std::vector< int > surface_vector;
  surface_vector.reserve( surfaces.size() );
  for( int surface : surfaces)
  {
    surface_vector.push_back(surface);
  }
  int nbSurfaceId = surface_vector.size();
  int totalNbSurfaceId = 0;
  MpiWrapper::allReduce( &nbSurfaceId, &totalNbSurfaceId, 1, MPI_SUM,MPI_COMM_GEOSX );
  allSurfaces.resize(totalNbSurfaceId);
  MpiWrapper::allgather<int, int>(surface_vector.data(),  surface_vector.size(), allSurfaces.data(),surface_vector.size(),MPI_COMM_GEOSX);

  std::sort(allSurfaces.begin(), allSurfaces.end());
  auto last = std::unique( allSurfaces.begin(), allSurfaces.end());
  allSurfaces.erase(last, allSurfaces.end());
}

std::vector< vtkDataArray * > VTKMeshGenerator::FindArrayToBeImported()
{
  std::vector< vtkDataArray * > arrayToBeImported;
  vtkCellData * cellData = m_vtkMesh->GetCellData();
  for( int i = 0; i < cellData->GetNumberOfArrays(); i++)
  {
    vtkAbstractArray * curArray = cellData->GetAbstractArray(i);
    // We support only float and double type. Good assumption?
    int dataType = curArray->GetDataType();
    if( !curArray->IsNumeric() ) continue;
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

void VTKMeshGenerator::WriteCellBlocks( DomainPartition& domain, localIndex nbHex, localIndex nbTet, localIndex nbWedge, localIndex nbPyr,
                                        std::map<int,localIndex> & regions_hex, std::map<int,localIndex> & regions_tetra,
                                        std::map<int,localIndex> & regions_wedges, std::map<int,localIndex> & regions_pyramids,
                                        std::set< int > & regions, std::vector< vtkDataArray * > & arraysToBeImported)
{
  CellBlockManager & cellBlockManager = domain.getGroup< CellBlockManager >( keys::cellManager );

  if( regions.size() > 0 )
  {
    for (int region : regions)
    {
      std::map< int, localIndex>::iterator it;
      it = regions_hex.find(region);
      if( it != regions_hex.end())
      {
        WriteCellBlock("hexahedron", it->second, it->first, VTK_HEXAHEDRON, cellBlockManager, arraysToBeImported);
      }
      it = regions_tetra.find(region);
      if( it != regions_tetra.end())
      {

        WriteCellBlock("tetrahedron", it->second, it->first, VTK_TETRA, cellBlockManager, arraysToBeImported);
      }
      it = regions_wedges.find(region);
      if( it != regions_wedges.end())
      {
        WriteCellBlock("wedges", it->second, it->first, VTK_WEDGE, cellBlockManager, arraysToBeImported);
      }
      it = regions_pyramids.find(region);
      if( it != regions_pyramids.end())
        WriteCellBlock("pyramids", it->second, it->first, VTK_PYRAMID, cellBlockManager, arraysToBeImported);
    }
  }
  else
  {
      WriteCellBlock("hexahedron", nbHex, -1, VTK_HEXAHEDRON, cellBlockManager, arraysToBeImported);
      WriteCellBlock("tetrahedron", nbTet, -1, VTK_TETRA, cellBlockManager, arraysToBeImported);
      WriteCellBlock("wedges", nbWedge, -1, VTK_WEDGE, cellBlockManager, arraysToBeImported);
      WriteCellBlock("pyramids", nbPyr, -1, VTK_PYRAMID, cellBlockManager, arraysToBeImported);
  }
}


void VTKMeshGenerator::WriteSurfaces( NodeManager & nodeManager, std::vector<int> const & allSurfaces )
{
  Group & nodeSets = nodeManager.sets();
  vtkIntArray * attributeDataArray = GetAttributeDataArray();
  for( int surface : allSurfaces)
  {
    SortedArray< localIndex > & curNodeSet  = nodeSets.registerWrapper< SortedArray< localIndex > >( std::to_string( surface) ).reference();
    for( vtkIdType c = 0; c < m_vtkMesh->GetNumberOfCells(); c++)
    {
      vtkCell * currentCell = m_vtkMesh->GetCell(c);
      if( m_vtkMesh->GetCellType(c) == VTK_TRIANGLE ||  m_vtkMesh->GetCellType(c) == VTK_QUAD)
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

localIndex VTKMeshGenerator::GetNumberOfPoints( int cellType )
{
  if( cellType == VTK_HEXAHEDRON )
  {
    return 8;
  }
  else if( cellType == VTK_TETRA )
  {
    return 4;
  }
  else if( cellType == VTK_WEDGE )
  {
    return 6;
  }
  else if( cellType == VTK_PYRAMID )
  {
    return 5;
  }
  else
  {
    GEOSX_ERROR( cellType << " is not a recognized cell type to be used with the VTKMeshGenerator");
    return -1;
  }

}

void VTKMeshGenerator::WriteHexahedronVertices( CellBlock::NodeMapType & cellToVertex, int region_id,  arrayView1d< globalIndex > const & localToGlobal  )
{
  localIndex cellCount = 0;
  vtkIntArray * attributeDataArray =  GetAttributeDataArray();
  vtkIdTypeArray * globalCellIdsDataArray = GetCellGlobalIdDataArray();
  for( vtkIdType c = 0; c < m_vtkMesh->GetNumberOfCells(); c++)
  {
    vtkCell * currentCell = m_vtkMesh->GetCell(c);
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
      localToGlobal[cellCount] = globalCellIdsDataArray->GetValue(c);
      cellCount++;
    }
  }
}

void VTKMeshGenerator::WriteWedgeVertices( CellBlock::NodeMapType & cellToVertex, int region_id,  arrayView1d< globalIndex > const & localToGlobal )
{
  localIndex cellCount = 0;
  vtkIntArray * attributeDataArray =  GetAttributeDataArray();
  vtkIdTypeArray * globalCellIdsDataArray = GetCellGlobalIdDataArray();
  for( vtkIdType c = 0; c < m_vtkMesh->GetNumberOfCells(); c++)
  {
    vtkCell * currentCell = m_vtkMesh->GetCell(c);
    if( currentCell->GetCellType() == VTK_WEDGE )
    {
      if( attributeDataArray != nullptr && attributeDataArray->GetValue(c) != region_id) continue;
		  cellToVertex[cellCount][0] = currentCell->GetPointId(0);
		  cellToVertex[cellCount][1] = currentCell->GetPointId(3);
		  cellToVertex[cellCount][2] = currentCell->GetPointId(1);
		  cellToVertex[cellCount][3] = currentCell->GetPointId(4);
		  cellToVertex[cellCount][4] = currentCell->GetPointId(2);
		  cellToVertex[cellCount][4] = currentCell->GetPointId(5);
      localToGlobal[cellCount] = globalCellIdsDataArray->GetValue(c);
      cellCount++;
    }
  }
}
void VTKMeshGenerator::WriteCellBlock( string const & name, localIndex nbCells, int region_id, int cellType,
                                       CellBlockManager & cellBlockManager,
                                       std::vector< vtkDataArray * > const & arraysTobeImported )
{

  if( nbCells > 0)
  {
    string cellBlockName;
    assert(region_id >= -1 );
    if( region_id != -1 )
    {
      cellBlockName = std::to_string(region_id) + "_" + name;
    }
    else
    {
      cellBlockName = name;
    }

    GEOSX_LOG_RANK_0( "Writing " << nbCells << " cells to CellBlock " << cellBlockName);
    CellBlock & cellBlock = cellBlockManager.getGroup( keys::cellBlocks ).registerGroup< CellBlock >( cellBlockName );
    localIndex nbPointsPerCell = GetNumberOfPoints( cellType );
    cellBlock.setElementType( "C3D" + std::to_string(nbPointsPerCell));
    CellBlock::NodeMapType & cellToVertex = cellBlock.nodeList();
    cellBlock.resize( nbCells );
    arrayView1d< globalIndex > const & localToGlobal = cellBlock.localToGlobalMap();

    vtkIntArray * attributeDataArray =  GetAttributeDataArray();


    /// Writing properties
    for( vtkDataArray * array : arraysTobeImported )
    {
      int dimension = array->GetNumberOfComponents();
      if( dimension == 1 )
      {
        real64_array & property = cellBlock.addProperty< real64_array >( array->GetName() );
        localIndex cellCount = 0;
        for( vtkIdType c = 0; c < m_vtkMesh->GetNumberOfCells(); c++)
        {
          vtkCell * currentCell = m_vtkMesh->GetCell(c);
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
        for( vtkIdType c = 0; c < m_vtkMesh->GetNumberOfCells(); c++)
        {
          vtkCell * currentCell = m_vtkMesh->GetCell(c);
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
    cellToVertex.resize( nbCells, nbPointsPerCell );
    if( cellType == VTK_HEXAHEDRON ) // Special case for hexahedron because of the ordering
    {
      WriteHexahedronVertices( cellToVertex, region_id, localToGlobal );
      return;
    }
    if( cellType == VTK_WEDGE ) // Special case for hexahedron because of the ordering
    {
      WriteWedgeVertices( cellToVertex, region_id, localToGlobal );
      return;
    }
    vtkIdTypeArray * globalCellIdsDataArray = GetCellGlobalIdDataArray();
    localIndex cellCount = 0;
    for( vtkIdType c = 0; c < m_vtkMesh->GetNumberOfCells(); c++)
    {
      vtkCell * currentCell = m_vtkMesh->GetCell(c);
      if( currentCell->GetCellType() == cellType )
      {
        if( attributeDataArray != nullptr && attributeDataArray->GetValue(c) != region_id) continue;
        for( localIndex v = 0; v < nbPointsPerCell; v++)
        {
          cellToVertex[cellCount][v] = currentCell->GetPointId(v);
        }
        localToGlobal[cellCount] = globalCellIdsDataArray->GetValue(c);
        cellCount++;
      }
    }
  }
}
REGISTER_CATALOG_ENTRY( MeshGeneratorBase, VTKMeshGenerator, string const &, Group * const )
}
