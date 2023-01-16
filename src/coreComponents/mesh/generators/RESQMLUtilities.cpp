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
 * @file RESQMLUtilities.cpp
 */

#include "mesh/generators/RESQMLUtilities.hpp"

#include <vtkNew.h>
#include <vtkSmartPointer.h>
#include <vtkIdList.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>


#include "fesapi/resqml2/AbstractLocal3dCrs.h"
#include "fesapi/resqml2/CategoricalProperty.h"
#include "fesapi/resqml2/ContinuousProperty.h"
#include "fesapi/resqml2/DiscreteProperty.h"
#include "fesapi/resqml2/SubRepresentation.h"

#include <array>

namespace geosx
{

//----------------------------------------------------------------------------
void cellVtkTetra( vtkSmartPointer< vtkUnstructuredGrid > vtk_unstructuredGrid,
                   const RESQML2_NS::UnstructuredGridRepresentation *unstructuredGridRep,
                   ULONG64 const *cumulativeFaceCountPerCell,
                   unsigned char const *cellFaceNormalOutwardlyDirected,
                   ULONG64 cellIndex )
{
  unsigned int nodes[4];

  // Face 0
  ULONG64 const *nodeIndices = unstructuredGridRep->getNodeIndicesOfFaceOfCell( cellIndex, 0 );
  size_t cellFaceIndex = unstructuredGridRep->isFaceCountOfCellsConstant() || cellIndex == 0
                 ? cellIndex * 4
                 : cumulativeFaceCountPerCell[cellIndex - 1];
  if( cellFaceNormalOutwardlyDirected[cellFaceIndex] == 0 )
  { // The RESQML orientation of face 0 honors the VTK orientation of face 0 i.e. the face 0 normal defined using a right hand rule is
    // inwardly directed.
    nodes[0] = nodeIndices[0];
    nodes[1] = nodeIndices[1];
    nodes[2] = nodeIndices[2];
  }
  else
  { // The RESQML orientation of face 0 does not honor the VTK orientation of face 0
    nodes[0] = nodeIndices[2];
    nodes[1] = nodeIndices[1];
    nodes[2] = nodeIndices[0];
  }

  // Face 1
  nodeIndices = unstructuredGridRep->getNodeIndicesOfFaceOfCell( cellIndex, 1 );

  for( size_t index = 0; index < 3; ++index )
  {
    if( std::find( nodes, nodes + 3, nodeIndices[index] ) == nodes + 3 )
    {
      nodes[3] = nodeIndices[index];
      break;
    }
  }

  vtkNew< vtkIdList > tetra;
  for( vtkIdType pointId = 0; pointId < 4; ++pointId )
  {
    tetra->InsertId( pointId, nodes[pointId] );
  }
  vtk_unstructuredGrid->InsertNextCell( VTK_TETRA, tetra );
}

//----------------------------------------------------------------------------
void cellVtkWedgeOrPyramid( vtkSmartPointer< vtkUnstructuredGrid > vtk_unstructuredGrid,
                            const RESQML2_NS::UnstructuredGridRepresentation *unstructuredGridRep,
                            ULONG64 const *cumulativeFaceCountPerCell, unsigned char const *cellFaceNormalOutwardlyDirected,
                            ULONG64 cellIndex )
{
  // The global index of the first face of the polyhedron in the cellFaceNormalOutwardlyDirected array
  const size_t globalFirstFaceIndex = unstructuredGridRep->isFaceCountOfCellsConstant() || cellIndex == 0
    ? cellIndex * 5
    : cumulativeFaceCountPerCell[cellIndex - 1];

  std::vector< unsigned int > localFaceIndexWith4Nodes;
  for( unsigned int localFaceIndex = 0; localFaceIndex < 5; ++localFaceIndex )
  {
    if( unstructuredGridRep->getNodeCountOfFaceOfCell( cellIndex, localFaceIndex ) == 4 )
    {
      localFaceIndexWith4Nodes.push_back( localFaceIndex );
    }
  }
  if( localFaceIndexWith4Nodes.size() == 3 )
  { // VTK_WEDGE
    std::array< uint64_t, 6 > nodes;
    // Set the triangle base of the wedge
    unsigned int triangleIndex = 0;
    for(; triangleIndex < 5; ++triangleIndex )
    {
      const unsigned int localNodeCount = unstructuredGridRep->getNodeCountOfFaceOfCell( cellIndex, triangleIndex );
      if( localNodeCount == 3 )
      {
        uint64_t const *nodeIndices = unstructuredGridRep->getNodeIndicesOfFaceOfCell( cellIndex, triangleIndex );
        if( cellFaceNormalOutwardlyDirected[globalFirstFaceIndex + triangleIndex] == 0 )
        {
          for( size_t i = 0; i < 3; ++i )
          {
            nodes[i] = nodeIndices[2 - i];
          }
        }
        else
        {
          // The RESQML orientation of face 0 honors the VTK orientation of face 0 i.e. the face 0 normal defined using a right hand rule is
          // outwardly directed.
          for( size_t i = 0; i < 3; ++i )
          {
            nodes[i] = nodeIndices[i];
          }
        }
        ++triangleIndex;
        break;
      }
    }
    // Find the index of the vertex at the opposite triangle regarding the triangle base
    for( unsigned int localFaceIndex = 0; localFaceIndex < 5; ++localFaceIndex )
    {
      const unsigned int localNodeCount = unstructuredGridRep->getNodeCountOfFaceOfCell( cellIndex, localFaceIndex );
      if( localNodeCount == 4 )
      {
        uint64_t const *nodeIndices = unstructuredGridRep->getNodeIndicesOfFaceOfCell( cellIndex, localFaceIndex );
        if( nodeIndices[0] == nodes[0] )
        {
          nodes[3] = nodeIndices[1] == nodes[1] || nodeIndices[1] == nodes[2]
            ? nodeIndices[3] : nodeIndices[1];
          break;
        }
        else if( nodeIndices[1] == nodes[0] )
        {
          nodes[3] = nodeIndices[2] == nodes[1] || nodeIndices[2] == nodes[2]
            ? nodeIndices[0] : nodeIndices[2];
          break;
        }
        else if( nodeIndices[2] == nodes[0] )
        {
          nodes[3] = nodeIndices[3] == nodes[1] || nodeIndices[3] == nodes[2]
            ? nodeIndices[1] : nodeIndices[3];
          break;
        }
        else if( nodeIndices[3] == nodes[0] )
        {
          nodes[3] = nodeIndices[0] == nodes[1] || nodeIndices[0] == nodes[2]
            ? nodeIndices[2] : nodeIndices[0];
          break;
        }
      }
    }
    // Set the other triangle of the wedge
    for(; triangleIndex < 5; ++triangleIndex )
    {
      const unsigned int localNodeCount = unstructuredGridRep->getNodeCountOfFaceOfCell( cellIndex, triangleIndex );
      if( localNodeCount == 3 )
      {
        uint64_t const *nodeIndices = unstructuredGridRep->getNodeIndicesOfFaceOfCell( cellIndex, triangleIndex );
        if( nodeIndices[0] == nodes[3] )
        {
          if( cellFaceNormalOutwardlyDirected[globalFirstFaceIndex + triangleIndex] == 0 )
          {
            nodes[4] = nodeIndices[1];
            nodes[5] = nodeIndices[2];
          }
          else
          {
            nodes[4] = nodeIndices[2];
            nodes[5] = nodeIndices[1];
          }
        }
        else if( nodeIndices[1] == nodes[3] )
        {
          if( cellFaceNormalOutwardlyDirected[globalFirstFaceIndex + triangleIndex] == 0 )
          {
            nodes[4] = nodeIndices[2];
            nodes[5] = nodeIndices[0];
          }
          else
          {
            nodes[4] = nodeIndices[0];
            nodes[5] = nodeIndices[2];
          }
        }
        else if( nodeIndices[2] == nodes[3] )
        {
          if( cellFaceNormalOutwardlyDirected[globalFirstFaceIndex + triangleIndex] == 0 )
          {
            nodes[4] = nodeIndices[0];
            nodes[5] = nodeIndices[1];
          }
          else
          {
            nodes[4] = nodeIndices[1];
            nodes[5] = nodeIndices[0];
          }
        }
        break;
      }
    }

    vtkNew< vtkIdList > wedge;
    for( int nodesIndex = 0; nodesIndex < 6; ++nodesIndex )
    {
      wedge->InsertId( nodesIndex, nodes[nodesIndex] );
    }
    vtk_unstructuredGrid->InsertNextCell( VTK_WEDGE, wedge );
  }
  else if( localFaceIndexWith4Nodes.size() == 1 )
  { // VTK_PYRAMID
    ULONG64 nodes[5];

    ULONG64 const *nodeIndices = unstructuredGridRep->getNodeIndicesOfFaceOfCell( cellIndex, localFaceIndexWith4Nodes[0] );
    size_t cellFaceIndex = (unstructuredGridRep->isFaceCountOfCellsConstant() || cellIndex == 0
      ? cellIndex * 5
      : cumulativeFaceCountPerCell[cellIndex - 1]) +
                           localFaceIndexWith4Nodes[0];
    if( cellFaceNormalOutwardlyDirected[cellFaceIndex] == 0 )
    { // The RESQML orientation of the face honors the VTK orientation of face 0 i.e. the face 0 normal defined using a right hand rule is
      // inwardly directed.
      nodes[0] = nodeIndices[0];
      nodes[1] = nodeIndices[1];
      nodes[2] = nodeIndices[2];
      nodes[3] = nodeIndices[3];
    }
    else
    { // The RESQML orientation of the face does not honor the VTK orientation of face 0
      nodes[0] = nodeIndices[3];
      nodes[1] = nodeIndices[2];
      nodes[2] = nodeIndices[1];
      nodes[3] = nodeIndices[0];
    }

    // Face with 3 points
    nodeIndices = unstructuredGridRep->getNodeIndicesOfFaceOfCell( cellIndex, localFaceIndexWith4Nodes[0] == 0 ? 1 : 0 );

    for( size_t index = 0; index < 3; ++index )
    {
      if( std::find( nodes, nodes + 4, nodeIndices[index] ) == nodes + 4 )
      {
        nodes[4] = nodeIndices[index];
        break;
      }
    }

    vtkNew< vtkIdList > pyramid;
    for( int nodesIndex = 0; nodesIndex < 5; ++nodesIndex )
    {
      pyramid->InsertId( nodesIndex, nodes[nodesIndex] );
    }
    vtk_unstructuredGrid->InsertNextCell( VTK_PYRAMID, pyramid );
  }
  else
  {
    throw std::invalid_argument( "The cell index " + std::to_string( cellIndex ) + " is malformed : 5 faces but not a pyramid, not a wedge." );
  }
}

//----------------------------------------------------------------------------
bool cellVtkHexahedron( vtkSmartPointer< vtkUnstructuredGrid > vtk_unstructuredGrid,
                        const RESQML2_NS::UnstructuredGridRepresentation *unstructuredGridRep,
                        ULONG64 const *cumulativeFaceCountPerCell,
                        unsigned char const *cellFaceNormalOutwardlyDirected,
                        ULONG64 cellIndex )
{
  for( ULONG64 localFaceIndex = 0; localFaceIndex < 6; ++localFaceIndex )
  {
    if( unstructuredGridRep->getNodeCountOfFaceOfCell( cellIndex, localFaceIndex ) != 4 )
    {
      return false;
    }
  }

  ULONG64 nodes[8];

  ULONG64 const *nodeIndices = unstructuredGridRep->getNodeIndicesOfFaceOfCell( cellIndex, 0 );
  const size_t cellFaceIndex = unstructuredGridRep->isFaceCountOfCellsConstant() || cellIndex == 0
                   ? cellIndex * 6
                   : cumulativeFaceCountPerCell[cellIndex - 1];
  if( cellFaceNormalOutwardlyDirected[cellFaceIndex] == 0 )
  { // The RESQML orientation of the face honors the VTK orientation of face 0 i.e. the face 0 normal defined using a right hand rule is
    // inwardly directed.
    nodes[0] = nodeIndices[0];
    nodes[1] = nodeIndices[1];
    nodes[2] = nodeIndices[2];
    nodes[3] = nodeIndices[3];
  }
  else
  { // The RESQML orientation of the face does not honor the VTK orientation of face 0
    nodes[0] = nodeIndices[3];
    nodes[1] = nodeIndices[2];
    nodes[2] = nodeIndices[1];
    nodes[3] = nodeIndices[0];
  }

  // Find the opposite neighbors of the nodes already got
  bool alreadyTreated[4] = {false, false, false, false};
  for( unsigned int localFaceIndex = 1; localFaceIndex < 6; ++localFaceIndex )
  {
    nodeIndices = unstructuredGridRep->getNodeIndicesOfFaceOfCell( cellIndex, localFaceIndex );
    for( size_t index = 0; index < 4; ++index )
    {                                 // Loop on face nodes
      ULONG64 *itr = std::find( nodes, nodes + 4, nodeIndices[index] ); // Locate a node on face 0
      if( itr != nodes + 4 )
      {
        // A top neighbor node can be found
        const size_t topNeigborIdx = std::distance( nodes, itr );
        if( !alreadyTreated[topNeigborIdx] )
        {
          const size_t previousIndex = index == 0 ? 3 : index - 1;
          nodes[topNeigborIdx + 4] = std::find( nodes, nodes + 4, nodeIndices[previousIndex] ) != nodes + 4 // If previous index is also in
                                                                                                            // face 0
                           ? nodeIndices[index == 3 ? 0 : index + 1]            // Put next index
                           : nodeIndices[previousIndex];                  // Put previous index
          alreadyTreated[topNeigborIdx] = true;
        }
      }
    }
    if( localFaceIndex > 2 && std::find( alreadyTreated, alreadyTreated + 4, false ) == alreadyTreated + 4 )
    {
      // All top neighbor nodes have been found. No need to continue
      // A minimum of four faces is necessary in order to find all top neighbor nodes.
      break;
    }
  }

  vtkNew< vtkIdList > hexahedron;
  for( int nodesIndex = 0; nodesIndex < 8; ++nodesIndex )
  {
    hexahedron->InsertId( nodesIndex, nodes[nodesIndex] );
  }
  vtk_unstructuredGrid->InsertNextCell( VTK_HEXAHEDRON, hexahedron );
  return true;
}

//----------------------------------------------------------------------------
bool cellVtkPentagonalPrism( vtkSmartPointer< vtkUnstructuredGrid > vtk_unstructuredGrid,
                             const RESQML2_NS::UnstructuredGridRepresentation *unstructuredGridRep,
                             ULONG64 cellIndex )
{
  unsigned int faceTo5Nodes = 0;
  ULONG64 nodes[10];
  for( ULONG64 localFaceIndex = 0; localFaceIndex < 7; ++localFaceIndex )
  {
    const unsigned int localNodeCount = unstructuredGridRep->getNodeCountOfFaceOfCell( cellIndex, localFaceIndex );
    if( localNodeCount == 5 )
    {
      ULONG64 const *nodeIndices = unstructuredGridRep->getNodeIndicesOfFaceOfCell( cellIndex, localFaceIndex );
      for( unsigned int i = 0; i < localNodeCount; ++i )
      {
        nodes[faceTo5Nodes * 5 + i] = nodeIndices[i];
      }
      ++faceTo5Nodes;
    }
  }
  if( faceTo5Nodes == 2 )
  {
    vtkNew< vtkIdList > pentagonalPrism;
    for( int nodesIndex = 0; nodesIndex < 10; ++nodesIndex )
    {
      pentagonalPrism->InsertId( nodesIndex, nodes[nodesIndex] );
    }
    vtk_unstructuredGrid->InsertNextCell( VTK_PENTAGONAL_PRISM, pentagonalPrism );
    return true;
  }
  return false;
}

//----------------------------------------------------------------------------
bool cellVtkHexagonalPrism( vtkSmartPointer< vtkUnstructuredGrid > vtk_unstructuredGrid,
                            const RESQML2_NS::UnstructuredGridRepresentation *unstructuredGridRep,
                            ULONG64 cellIndex )
{
  unsigned int faceTo6Nodes = 0;
  ULONG64 nodes[12];
  for( ULONG64 localFaceIndex = 0; localFaceIndex < 8; ++localFaceIndex )
  {
    const unsigned int localNodeCount = unstructuredGridRep->getNodeCountOfFaceOfCell( cellIndex, localFaceIndex );
    if( localNodeCount == 6 )
    {
      ULONG64 const *nodeIndices = unstructuredGridRep->getNodeIndicesOfFaceOfCell( cellIndex, localFaceIndex );
      for( unsigned int i = 0; i < localNodeCount; ++i )
      {
        nodes[faceTo6Nodes * 6 + i] = nodeIndices[i];
      }
      ++faceTo6Nodes;
    }
  }
  if( faceTo6Nodes == 2 )
  {
    vtkNew< vtkIdList > hexagonalPrism;
    for( int nodesIndex = 0; nodesIndex < 12; ++nodesIndex )
    {
      hexagonalPrism->InsertId( nodesIndex, nodes[nodesIndex] );
    }
    vtk_unstructuredGrid->InsertNextCell( VTK_HEXAGONAL_PRISM, hexagonalPrism );
    return true;
  }
  return false;
}


int readContinuousProperty( RESQML2_NS::AbstractValuesProperty * valuesProperty, string name, vtkCellData * outDS )
{
  const unsigned int elementCountPerValue = valuesProperty->getElementCountPerValue();
  const unsigned int totalHDFElementcount = valuesProperty->getValuesCountOfPatch( 0 );

  double * valueIndices = new double[totalHDFElementcount]; // deleted by VTK cellData vtkSmartPointer

  valuesProperty->getDoubleValuesOfPatch( 0, valueIndices );

  //TODO find a better way to handle NaN data
  for( unsigned int x = 0; x < totalHDFElementcount; ++x )
  {
    if( std::isnan( valueIndices[x] ))
    {
      valueIndices[x] = 0.00000001;
    }
  }

  vtkNew< vtkDoubleArray > cellData;
  cellData->SetArray( valueIndices, totalHDFElementcount, 0, vtkAbstractArray::VTK_DATA_ARRAY_DELETE );
  cellData->SetName( name.c_str());
  cellData->SetNumberOfComponents( elementCountPerValue );

  outDS->AddArray( cellData );

  return 1;
}

int readDiscreteOrCategoricalProperty( RESQML2_NS::AbstractValuesProperty * valuesProperty, string name, vtkCellData * outDS )
{
  const unsigned int elementCountPerValue = valuesProperty->getElementCountPerValue();
  const unsigned int totalHDFElementcount = valuesProperty->getValuesCountOfPatch( 0 );

  int * valueIndices = new int[totalHDFElementcount]; // deleted by VTK cellData vtkSmartPointer

  valuesProperty->getIntValuesOfPatch( 0, valueIndices );

  //TODO handle NaN ?

  vtkNew< vtkIntArray > cellData;
  cellData->SetArray( valueIndices, totalHDFElementcount, 0, vtkAbstractArray::VTK_DATA_ARRAY_DELETE );
  cellData->SetName( name.c_str());
  cellData->SetNumberOfComponents( elementCountPerValue );

  outDS->AddArray( cellData );

  return 1;
}


vtkSmartPointer< vtkUnstructuredGrid >
loadUnstructuredGridRepresentation( RESQML2_NS::UnstructuredGridRepresentation *grid )
{
  auto vtk_unstructuredGrid = vtkSmartPointer< vtkUnstructuredGrid >::New();

  vtk_unstructuredGrid->Allocate( grid->getCellCount());


  uint64_t pointCount = grid->getXyzPointCountOfAllPatches();

  // POINTS
  double *allXyzPoints = new double[pointCount * 3]; // Will be deleted by VTK;
  grid->getXyzPointsOfAllPatches( allXyzPoints ); //getXyzPointsOfAllPatchesInGlobalCrs

  const size_t coordCount = pointCount * 3;

  if( grid->getLocalCrs( 0 ) != nullptr )
  {
    if( !grid->getLocalCrs( 0 )->isPartial())
    {
      const double zIndice = grid->getLocalCrs( 0 )->isDepthOriented() ? -1 : 1;

      for( size_t zCoordIndex = 2; zCoordIndex < coordCount; zCoordIndex += 3 )
      {
        allXyzPoints[zCoordIndex] *= zIndice;
      }
    }
  }

  vtkNew< vtkDoubleArray > vtkUnderlyingArray;
  vtkUnderlyingArray->SetNumberOfComponents( 3 );
  // Take ownership of the underlying C array
  vtkUnderlyingArray->SetArray( allXyzPoints, coordCount, vtkAbstractArray::VTK_DATA_ARRAY_DELETE );

  vtkNew< vtkPoints > vtkPts;
  vtkPts->SetData( vtkUnderlyingArray );

  vtk_unstructuredGrid->SetPoints( vtkPts );
  grid->loadGeometry();
  // CELLS
  const unsigned int nodeCount = grid->getNodeCount();
  std::unique_ptr< vtkIdType[] > pointIds( new vtkIdType[nodeCount] );
  for( unsigned int i = 0; i < nodeCount; ++i )
  {
    pointIds[i] = i;
  }
  const ULONG64 cellCount = grid->getCellCount();
  // This pointer is owned and managed by FESAPI
  ULONG64 const *cumulativeFaceCountPerCell = grid->isFaceCountOfCellsConstant()
                          ? nullptr
                          : grid->getCumulativeFaceCountPerCell();

  const uint64_t faceCount = cumulativeFaceCountPerCell == nullptr
    ? cellCount * grid->getConstantFaceCountOfCells()
    : cumulativeFaceCountPerCell[cellCount - 1];

  std::unique_ptr< unsigned char[] > cellFaceNormalOutwardlyDirected( new unsigned char[faceCount] );

  grid->getCellFaceIsRightHanded( cellFaceNormalOutwardlyDirected.get());

  auto * crs = grid->getLocalCrs( 0 );
  if( !crs->isPartial() && crs->isDepthOriented())
  {
    for( size_t i = 0; i < faceCount; ++i )
    {
      cellFaceNormalOutwardlyDirected[i] = !cellFaceNormalOutwardlyDirected[i];
    }
  }

  for( ULONG64 cellIndex = 0; cellIndex < cellCount; ++cellIndex )
  {
    bool isOptimizedCell = false;


    const ULONG64 localFaceCount = grid->getFaceCountOfCell( cellIndex );

    if( localFaceCount == 4 )
    { // VTK_TETRA
      cellVtkTetra( vtk_unstructuredGrid, grid, cumulativeFaceCountPerCell, cellFaceNormalOutwardlyDirected.get(), cellIndex );
      isOptimizedCell = true;
    }
    else if( localFaceCount == 5 )
    { // VTK_WEDGE or VTK_PYRAMID
      cellVtkWedgeOrPyramid( vtk_unstructuredGrid, grid, cumulativeFaceCountPerCell, cellFaceNormalOutwardlyDirected.get(), cellIndex );
      isOptimizedCell = true;
    }
    else if( localFaceCount == 6 )
    { // VTK_HEXAHEDRON
      isOptimizedCell = cellVtkHexahedron( vtk_unstructuredGrid, grid, cumulativeFaceCountPerCell, cellFaceNormalOutwardlyDirected.get(), cellIndex );
    }
    else if( localFaceCount == 7 )
    { // VTK_PENTAGONAL_PRISM
      isOptimizedCell = cellVtkPentagonalPrism( vtk_unstructuredGrid, grid, cellIndex );
    }
    else if( localFaceCount == 8 )
    { // VTK_HEXAGONAL_PRISM
      isOptimizedCell = cellVtkHexagonalPrism( vtk_unstructuredGrid, grid, cellIndex );
    }

    if( !isOptimizedCell )
    {
      vtkNew< vtkIdList > idList;

      // For polyhedron cell, a special ptIds input format is required : (numCellFaces, numFace0Pts, id1, id2, id3, numFace1Pts, id1, id2,
      // id3, ...)
      idList->InsertNextId( localFaceCount );
      for( ULONG64 localFaceIndex = 0; localFaceIndex < localFaceCount; ++localFaceIndex )
      {
        const unsigned int localNodeCount = grid->getNodeCountOfFaceOfCell( cellIndex, localFaceIndex );
        idList->InsertNextId( localNodeCount );
        ULONG64 const *nodeIndices = grid->getNodeIndicesOfFaceOfCell( cellIndex, localFaceIndex );
        for( unsigned int i = 0; i < localNodeCount; ++i )
        {
          idList->InsertNextId( nodeIndices[i] );
        }
      }

      vtk_unstructuredGrid->InsertNextCell( VTK_POLYHEDRON, idList );
    }
  }

  grid->unloadGeometry();

  return vtk_unstructuredGrid;
}

vtkSmartPointer< vtkDataSet >
loadProperty( vtkSmartPointer< vtkDataSet > dataset, RESQML2_NS::AbstractValuesProperty *valuesProperty, string fieldNameInGEOSX )
{
  const gsoap_eml2_3::resqml22__IndexableElement element = valuesProperty->getAttachmentKind();

  if( element != gsoap_eml2_3::resqml22__IndexableElement::cells )
  {
    GEOSX_ERROR( GEOSX_FMT( "Property indexable element must be cells. Error with {}", valuesProperty->getUuid() ) );
  }

  std::string typeProperty = valuesProperty->getXmlTag();
  if( typeProperty == RESQML2_NS::ContinuousProperty::XML_TAG )
  {
    readContinuousProperty( valuesProperty, fieldNameInGEOSX, dataset->GetCellData() );
  }
  else if( typeProperty == RESQML2_NS::DiscreteProperty::XML_TAG ||
           typeProperty == RESQML2_NS::CategoricalProperty::XML_TAG )
  {
    readDiscreteOrCategoricalProperty( valuesProperty, fieldNameInGEOSX, dataset->GetCellData() );
  }
  else
  {
    GEOSX_ERROR( GEOSX_FMT( "Property {} not supported...", valuesProperty->getUuid() ) );
  }

  return dataset;
}

vtkSmartPointer< vtkDataSet >
createRegions( vtkSmartPointer< vtkDataSet > dataset, std::vector< RESQML2_NS::SubRepresentation * > regions, string attributeName )
{
  int region_id = 0;

  vtkIntArray * arr = vtkIntArray::New();
  arr->SetName( attributeName.c_str());
  arr->SetNumberOfComponents( 1 );
  arr->SetNumberOfTuples( dataset->GetNumberOfCells());
  arr->FillValue( region_id );

  vtkDataArrayAccessor< vtkIntArray > attribute( arr );

  for( std::size_t i = 0; i < regions.size(); ++i )
  {
    region_id+=1;

    RESQML2_NS::SubRepresentation * region = regions[i];

    if( region->getElementKindOfPatch( 0, 0 ) != gsoap_eml2_3::resqml22__IndexableElement::faces )
    {
      uint64_t elementCountOfPatch = region->getElementCountOfPatch( 0 );
      std::unique_ptr< uint64_t[] > elementIndices( new uint64_t[elementCountOfPatch] );
      region->getElementIndicesOfPatch( 0, 0, elementIndices.get() );

      for( std::size_t j = 0; j < elementCountOfPatch; ++j )
      {
        attribute.Set( vtkIdType( elementIndices[j] ), 0, region_id );
      }
    }
  }
  dataset->GetCellData()->AddArray( arr );

  return dataset;
}

vtkUnstructuredGrid *
createSurfaces( vtkSmartPointer< vtkDataSet > dataset, std::vector< RESQML2_NS::SubRepresentation * > surfaces, string regionAttributeName )
{
  vtkUnstructuredGrid * grid = vtkUnstructuredGrid::SafeDownCast( dataset );

  vtkCellData * cell_data = grid->GetCellData();
  auto * region_ids = vtkIntArray::SafeDownCast( cell_data->GetAbstractArray( regionAttributeName.c_str()));
  auto * min_max = region_ids->GetRange();
  int region_id = min_max[1];

  for( std::size_t i = 0; i < surfaces.size(); ++i )
  {
    RESQML2_NS::SubRepresentation * surface = surfaces[i];

    // FACES
    auto *supportingGrid = dynamic_cast< RESQML2_NS::UnstructuredGridRepresentation * >(surface->getSupportingRepresentation( 0 ));

    const uint64_t gridFaceCount = supportingGrid->getFaceCount();
    std::unique_ptr< uint64_t[] > nodeCountOfFaces( new uint64_t[gridFaceCount] );
    supportingGrid->getCumulativeNodeCountPerFace( nodeCountOfFaces.get());

    std::unique_ptr< uint64_t[] > nodeIndices( new uint64_t[nodeCountOfFaces[gridFaceCount-1]] );
    supportingGrid->getNodeIndicesOfFaces( nodeIndices.get());

    const uint64_t subFaceCount = surface->getElementCountOfPatch( 0 );
    std::unique_ptr< uint64_t[] > elementIndices( new uint64_t[subFaceCount] );
    surface->getElementIndicesOfPatch( 0, 0, elementIndices.get());

    region_id += 1;
    for( uint64_t subFaceIndex =0; subFaceIndex < subFaceCount; ++subFaceIndex )
    {
      uint64_t faceIndex = elementIndices[subFaceIndex];
      auto first_indiceValue = faceIndex == 0 ? 0 : nodeCountOfFaces[faceIndex - 1];
      uint64_t nodeCount_OfFaceIndex = nodeCountOfFaces[faceIndex] - first_indiceValue;

      vtkIdType cell_type = VTK_POLYGON;
      if( nodeCount_OfFaceIndex == 3 )
      {
        cell_type = VTK_TRIANGLE;
      }
      else if( nodeCount_OfFaceIndex == 4 )
      {
        cell_type = VTK_QUAD;
      }

      vtkNew< vtkIdList > idList;
      for( uint64_t nodeIndex = 0; nodeIndex < nodeCount_OfFaceIndex; ++nodeIndex, ++first_indiceValue )
      {
        idList->InsertId( nodeIndex, nodeIndices[first_indiceValue] );
      }

      vtkIdType newCellId = grid->InsertNextCell( cell_type, idList );

      for( int arrayIdx = 0; arrayIdx < cell_data->GetNumberOfArrays(); ++arrayIdx )
      {
        auto * abArray = cell_data->GetAbstractArray( arrayIdx );
        if( abArray->GetName()==regionAttributeName )
        {
          vtkIntArray * ar = vtkArrayDownCast< vtkIntArray >( abArray );
          ar->InsertValue( newCellId, region_id );
        }
        else
        {
          int numComps = abArray->GetNumberOfComponents();
          vtkDoubleArray * ar = vtkArrayDownCast< vtkDoubleArray >( abArray );
          if( ar != nullptr )
          {
            for( int comp = 0; comp < numComps; ++comp )
            {
              ar->InsertTypedComponent( newCellId, comp, -9999. );
            }
          }
        }
      }
    }
  }

  return grid;
}

} // namespace geosx
