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

#ifndef GEOSX_FACEBLOCK_HPP
#define GEOSX_FACEBLOCK_HPP

#include "FaceBlockABC.hpp"

namespace geosx {

//class FaceBlock: public FaceBlockABC
//{
//public:
//  FaceBlock( string const & name, Group * const parent );
//
//  localIndex num2dElements() const override;
//
//  localIndex num2dFaces() const override;
//
//  ArrayOfArrays< localIndex > get2dElemToNodes() const override;
//
//  ArrayOfArrays< localIndex > get2dElemToEdges() const override;
//
//  array2d< localIndex > get2dElemToFaces() const override;
//
//  ToCellRelation< array2d< localIndex > > get2dElemToElems() const override;
//
//  array1d< localIndex > get2dFaceToEdge() const override;
//
//  ArrayOfArrays< localIndex > get2dFaceTo2dElems() const override;
//};

class FaceBlock: public FaceBlockABC
{
public:
  FaceBlock( string const & name,
             Group * const parent )
    :
    FaceBlockABC( name, parent )
  { }

  localIndex num2dElements() const override;

  localIndex num2dFaces() const override;

  ArrayOfArrays< localIndex > get2dElemToNodes() const override;

  ArrayOfArrays< localIndex > get2dElemToEdges() const override;

  array2d< localIndex > get2dElemToFaces() const override;

  ToCellRelation< array2d< localIndex > > get2dElemToElems() const override;

  array1d< localIndex > get2dFaceToEdge() const override;

  ArrayOfArrays< localIndex > get2dFaceTo2dElems() const override;

  void setNum2DElements( localIndex num2DElements );

  void setNum2DFaces( localIndex num2DFaces );

  void set2dElemToNodes( std::vector< std::vector< localIndex > > _2dElemToNodes );

  void set2dElemToEdges( std::vector< std::vector< localIndex > > _2dElemToEdges );

  void set2dElemToFaces( std::vector< std::vector< localIndex > > _2dElemToFaces );

  void set2dFaceTo2dElems( std::vector< std::vector< localIndex > > _2dFaceTo2dElems );

  void set2dFaceToEdge( std::vector< localIndex > _2dFaceToEdge );

  void set2dElemToElems( std::vector< std::vector< localIndex > > _2dElemToElems );

private:

  localIndex m_num2dElements;
  localIndex m_num2dFaces;

  std::vector< std::vector< localIndex > > m_2dElemToNodes;
  std::vector< std::vector< localIndex > > m_2dElemToEdges;
  std::vector< std::vector< localIndex > > m_2dElemToFaces;
  std::vector< std::vector< localIndex > > m_2dFaceTo2dElems;
  std::vector< localIndex > m_2dFaceToEdge;

  std::vector< std::vector< localIndex > > m_2dElemToElems;
};


}

#endif // include guard
