#include "FaceBlock.hpp"

namespace geosx
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

array2d< localIndex > FaceBlock::get2dElemToFaces() const
{
  return m_2dElemToFaces;
}

ToCellRelation< array2d< localIndex > > FaceBlock::get2dElemToElems() const
{
  return ToCellRelation< array2d< localIndex > >( m_2dElemToCellBlock, m_2dElemToElems );
}

array1d< localIndex > FaceBlock::get2dFaceToEdge() const
{
  return m_2dFaceToEdge;
}

ArrayOfArrays< localIndex > FaceBlock::get2dFaceTo2dElems() const
{
  return m_2dFaceTo2dElems;
}

void FaceBlock::setNum2DElements( localIndex num2DElements )
{
  m_num2dElements = num2DElements;
}

void FaceBlock::setNum2DFaces( localIndex num2DFaces )
{
  m_num2dFaces = num2DFaces;
}

void FaceBlock::set2dElemToNodes( ArrayOfArrays< localIndex > _2dElemToNodes )
{
  m_2dElemToNodes = _2dElemToNodes;
}

void FaceBlock::set2dElemToEdges( ArrayOfArrays< localIndex > _2dElemToEdges )
{
  m_2dElemToEdges = _2dElemToEdges;
}

void FaceBlock::set2dElemToFaces( array2d< localIndex > _2dElemToFaces )
{
  m_2dElemToFaces = _2dElemToFaces;
}

void FaceBlock::set2dFaceTo2dElems( ArrayOfArrays< localIndex > _2dFaceTo2dElems )
{
  m_2dFaceTo2dElems = _2dFaceTo2dElems;
}

void FaceBlock::set2dFaceToEdge( array1d< localIndex > _2dFaceToEdge )
{
  m_2dFaceToEdge = _2dFaceToEdge;
}

void FaceBlock::set2dElemToElems( array2d< localIndex > _2dElemToElems )
{
  m_2dElemToElems = _2dElemToElems;
}

void FaceBlock::set2dElemToCellBlock( array2d< localIndex > _2dElemToCellBlock )
{
  m_2dElemToCellBlock = _2dElemToCellBlock;
}

}
