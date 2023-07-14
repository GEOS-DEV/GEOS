/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2020-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "FaceBlock.hpp"

namespace geos
{

localIndex FaceBlock::num2dElements() const
{
  return m_num2dElements;
}

localIndex FaceBlock::num2dFaces() const
{
  return m_num2dFaces;
}

ArrayOfArrays< localIndex > FaceBlock::get2dElemToNodes() const
{
  return m_2dElemToNodes;
}

ArrayOfArrays< localIndex > FaceBlock::get2dElemToEdges() const
{
  return m_2dElemToEdges;
}

ArrayOfArrays< localIndex > FaceBlock::get2dElemToFaces() const
{
  return m_2dElemToFaces;
}

ToCellRelation< ArrayOfArrays< localIndex > > FaceBlock::get2dElemToElems() const
{
  return m_2dElemToElems;
}

array1d< localIndex > FaceBlock::get2dFaceToEdge() const
{
  return m_2dFaceToEdge;
}

ArrayOfArrays< localIndex > FaceBlock::get2dFaceTo2dElems() const
{
  return m_2dFaceTo2dElems;
}

void FaceBlock::setNum2dElements( localIndex num2DElements )
{
  m_num2dElements = num2DElements;
}

void FaceBlock::setNum2dFaces( localIndex num2DFaces )
{
  m_num2dFaces = num2DFaces;
}

void FaceBlock::set2dElemToNodes( ArrayOfArrays< localIndex > && _2dElemToNodes )
{
  m_2dElemToNodes = _2dElemToNodes;
}

void FaceBlock::set2dElemToEdges( ArrayOfArrays< localIndex > && _2dElemToEdges )
{
  m_2dElemToEdges = _2dElemToEdges;
}

void FaceBlock::set2dElemToFaces( ArrayOfArrays< localIndex > && _2dElemToFaces )
{
  m_2dElemToFaces = _2dElemToFaces;
}

void FaceBlock::set2dFaceTo2dElems( ArrayOfArrays< localIndex > && _2dFaceTo2dElems )
{
  m_2dFaceTo2dElems = _2dFaceTo2dElems;
}

void FaceBlock::set2dFaceToEdge( array1d< localIndex > && _2dFaceToEdge )
{
  m_2dFaceToEdge = _2dFaceToEdge;
}

void FaceBlock::set2dElemToElems( ToCellRelation< ArrayOfArrays< localIndex > > && _2dElemToElems )
{
  m_2dElemToElems = _2dElemToElems;
}

ArrayOfArrays< globalIndex > FaceBlock::getCollocatedNodes() const
{
  return m_collocatedNodes;
}

void FaceBlock::setCollocatedNodes( ArrayOfArrays< globalIndex > && collocatedNodes )
{
  m_collocatedNodes = collocatedNodes;
}

ArrayOfArrays< globalIndex > FaceBlock::getCollocatedNodesOf2dElems() const
{
  return m_collocatedNodesOf2dElems;
}

void FaceBlock::setCollocatedNodesOf2dElems( ArrayOfArrays< geos::globalIndex > && collocatedNodes )
{
  m_collocatedNodesOf2dElems = collocatedNodes;
}

}
