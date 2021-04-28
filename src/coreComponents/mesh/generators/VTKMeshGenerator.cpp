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
#include "Mesh/MeshFactory.hpp"

#include "MeshDataWriters/MeshParts.hpp"

#include "mesh/MeshBody.hpp"

#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkCellData.h>
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
  string_array filePathTokenized = stringutilities::tokenize( m_filePath, "."); //TODO maybe code a method in Path to get the file extension?
  string extension = filePathTokenized[filePathTokenized.size() - 1];
  if( extension == "vtk")
  {
    vtkSmartPointer<vtkUnstructuredGridReader> vtkUgReader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
	vtkUgReader->SetFileName( m_filePath.c_str() );
    vtkUgReader->Update();
	m_vtkMesh = vtkUgReader->GetOutput();
  }
  else if( extension == "vtu")
  {
    vtkSmartPointer<vtkXMLUnstructuredGridReader> vtkUgReader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    GEOSX_LOG_RANK_0(m_filePath);
	vtkUgReader->SetFileName( m_filePath.c_str() );
    vtkUgReader->Update();
	m_vtkMesh = vtkUgReader->GetOutput();
  }
  else
  {
    GEOSX_ERROR( extension << " is not a recognized extension for using the VTK reader with GEOSX. Please use .vtk or .vtu" );
  }
}

Group * VTKMeshGenerator::createChild( string const & GEOSX_UNUSED_PARAM( childKey ), string const & GEOSX_UNUSED_PARAM( childName ) )
{
  return nullptr;
}

void VTKMeshGenerator::generateMesh( DomainPartition & domain )
{
  Group & meshBodies = domain.getGroup( string( "MeshBodies" ));
  MeshBody & meshBody = meshBodies.registerGroup< MeshBody >( this->getName() );

  //TODO for the moment we only consider on mesh level "Level0"
  MeshLevel & meshLevel0 = meshBody.registerGroup< MeshLevel >( string( "Level0" ));
  NodeManager & nodeManager = meshLevel0.getNodeManager();

  GEOSX_ERROR_IF( m_vtkMesh->GetNumberOfPoints() == 0, "Mesh is empty, aborting");
  GEOSX_LOG_RANK_0("Writing " << m_vtkMesh->GetNumberOfPoints() << " vertices");
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();
  nodeManager.resize( m_vtkMesh->GetNumberOfPoints() );
  Group & nodeSets = nodeManager.sets();
  SortedArray< localIndex > & allNodes  = nodeSets.registerWrapper< SortedArray< localIndex > >( string( "all" ) ).reference();
  real64 xMax[3] = { std::numeric_limits< real64 >::min() };
  real64 xMin[3] = { std::numeric_limits< real64 >::max() };
  for( vtkIdType v = 0; v < m_vtkMesh->GetNumberOfPoints(); v++)
  {
    double * point = m_vtkMesh->GetPoint( v );
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
  meshBody.setGlobalLengthScale( LvArray::tensorOps::l2Norm< 3 >( xMax ) );
  //TODO For the moment, we only consider ONE element region
  //TODO For the moment, we only deal with classical elements (tetra, hexa, wedges and pyramids)
  
  /// First pass to count the cells, and the regions
  localIndex nbHex =0;
  localIndex nbTet = 0;
  localIndex nbWedge = 0;
  localIndex nbPyr = 0;
  localIndex nbOthers = 0;
  // region id link to their number of cells
  std::set<int> regions;
  std::map<int,localIndex> regions_tetra;
  std::map<int,localIndex> regions_hex;
  std::map<int,localIndex> regions_wedges;
  std::map<int,localIndex> regions_pyramids;
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
    else
    {
      nbOthers++;
    }
  }

  GEOSX_ERROR_IF( nbHex + nbTet + nbWedge + nbPyr == 0, "Mesh contains no cell");
  GEOSX_LOG_RANK_0_IF( nbOthers > 0, nbOthers << " cells with unsupported types will be ignored and not written");
  GEOSX_LOG_RANK_0( regions.size() << " regions will be defined");
  CellBlockManager & cellBlockManager = domain.getGroup< CellBlockManager >( keys::cellManager );

  ///Write properties
  std::vector< vtkDataArray * > arrayToBeImported;
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
  if( regions.size() > 0 )
  {
    for (int region : regions)
    {
      std::map< int, localIndex>::iterator it;
      it = regions_hex.find(region);
      if( it != regions_hex.end())
      {
        WriteCellBlock("hexahedron", it->second, it->first, VTK_HEXAHEDRON, cellBlockManager, arrayToBeImported);
      }
      it = regions_tetra.find(region);
      if( it != regions_tetra.end())
      {

        WriteCellBlock("tetrahedron", it->second, it->first, VTK_TETRA, cellBlockManager, arrayToBeImported);
      }
      it = regions_wedges.find(region);
      if( it != regions_wedges.end())
      {
        WriteCellBlock("wedges", it->second, it->first, VTK_WEDGE, cellBlockManager, arrayToBeImported);
      }
      it = regions_pyramids.find(region);
      if( it != regions_pyramids.end())
        WriteCellBlock("pyramids", it->second, it->first, VTK_PYRAMID, cellBlockManager, arrayToBeImported);
    }
  }
  else
  {
      WriteCellBlock("hexahedron", nbHex, -1, VTK_HEXAHEDRON, cellBlockManager, arrayToBeImported);
      WriteCellBlock("tetrahedron", nbTet, -1, VTK_TETRA, cellBlockManager, arrayToBeImported);
      WriteCellBlock("wedges", nbWedge, -1, VTK_WEDGE, cellBlockManager, arrayToBeImported);
      WriteCellBlock("pyramids", nbPyr, -1, VTK_PYRAMID, cellBlockManager, arrayToBeImported);

  }
  std::cout << "C FINI" << std::endl;
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

void VTKMeshGenerator::WriteHexahedronVertices( CellBlock::NodeMapType & cellToVertex, int region_id )
{
  localIndex cellCount = 0;
  //TODO : duplicate code to get the attribute array
  int attributeIndex = -1;
  vtkCellData * cellData = m_vtkMesh->GetCellData();
  vtkAbstractArray * attributeArray = cellData->GetAbstractArray("attribute", attributeIndex);
  GEOSX_ERROR_IF(attributeIndex >=0 && attributeArray->GetDataType() != VTK_INT, "Array name \"attribute\" is reserved to store the regions " 
                << "used by GEOSX. The underlying type has to be an INT");
  vtkIntArray * attributeDataArray =  nullptr;
  if( attributeIndex >= 0)
  {
    attributeDataArray = vtkIntArray::SafeDownCast( attributeArray );
  }
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
      cellCount++;
    }
  }
}

void VTKMeshGenerator::WritePyramidVertices( CellBlock::NodeMapType & cellToVertex )
{
  localIndex cellCount = 0;
  for( vtkIdType c = 0; c < m_vtkMesh->GetNumberOfCells(); c++)
  {
    vtkCell * currentCell = m_vtkMesh->GetCell(c);
    if( currentCell->GetCellType() == VTK_PYRAMID )
    {
		  cellToVertex[cellCount][0] = currentCell->GetPointId(0);
		  cellToVertex[cellCount][1] = currentCell->GetPointId(1);
		  cellToVertex[cellCount][2] = currentCell->GetPointId(2);
		  cellToVertex[cellCount][3] = currentCell->GetPointId(3);
		  cellToVertex[cellCount][4] = currentCell->GetPointId(4);
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

    //TODO : duplicate code to get the attribute array
    int attributeIndex = -1;
    vtkCellData * cellData = m_vtkMesh->GetCellData();
    vtkAbstractArray * attributeArray = cellData->GetAbstractArray("attribute", attributeIndex);
    GEOSX_ERROR_IF(attributeIndex >=0 && attributeArray->GetDataType() != VTK_INT, "Array name \"attribute\" is reserved to store the regions " 
                  << "used by GEOSX. The underlying type has to be an INT");
    vtkIntArray * attributeDataArray =  nullptr;
    if( attributeIndex >= 0)
    {
      attributeDataArray = vtkIntArray::SafeDownCast( attributeArray );
    }

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
        for( vtkIdType c = 0; c < m_vtkMesh->GetNumberOfCells(); c++)
        {
          vtkCell * currentCell = m_vtkMesh->GetCell(c);
          localIndex cellCount = 0;
          if( currentCell->GetCellType() == cellType )
          {
            if( attributeDataArray != nullptr && attributeDataArray->GetValue(c) != region_id) continue;
            double * tuple = array->GetTuple(c);
            for(int i = 0; i < dimension; i++ )
            {
              property(cellCount,i) = tuple[i];
            }
          }
        }
      }
    }

    /// Writing connectivity
    cellToVertex.resize( nbCells, nbPointsPerCell );
    if( cellType == VTK_HEXAHEDRON ) // Special case for hexahedron because of the ordering
    {
      WriteHexahedronVertices( cellToVertex, region_id );
      return;
    }
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
        cellCount++;
      }
    }
  }
}
REGISTER_CATALOG_ENTRY( MeshGeneratorBase, VTKMeshGenerator, string const &, Group * const )
}
