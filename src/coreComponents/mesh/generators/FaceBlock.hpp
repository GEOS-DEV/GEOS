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

#ifndef GEOS_FACEBLOCK_HPP
#define GEOS_FACEBLOCK_HPP

#include "FaceBlockABC.hpp"

namespace geos
{

/**
 * @brief Simple implementation of the @p FaceBlockABC contract.
 *
 * This class contains setters to store the externally computed mappings required by the @p FaceBlockABC contract.
 */
class FaceBlock : public FaceBlockABC
{
public:
  /**
   * @brief Constructor.
   * @param[in] name Name of this FaceBlock.
   * @param[in] parent Parent group.
   */
  FaceBlock( string const & name,
             Group * const parent )
    :
    FaceBlockABC( name, parent )
  { }

  localIndex num2dElements() const override;

  localIndex num2dFaces() const override;

  ArrayOfArrays< localIndex > get2dElemToNodes() const override;

  ArrayOfArrays< localIndex > get2dElemToEdges() const override;

  ArrayOfArrays< localIndex > get2dElemToFaces() const override;

  ToCellRelation< ArrayOfArrays< localIndex > > get2dElemToElems() const override;

  array1d< localIndex > get2dFaceToEdge() const override;

  ArrayOfArrays< localIndex > get2dFaceTo2dElems() const override;

  /**
   * @brief Defines the number of 2d elements.
   * @param num2DElements The input value.
   */
  void setNum2dElements( localIndex num2DElements );

  /**
   * @brief Defines the number of 2d faces.
   * @param num2DFaces The input value.
   */
  void setNum2dFaces( localIndex num2DFaces );

  /**
   * @brief Defines the 2d elements to nodes mapping.
   * @param _2dElemToNodes The input mapping.
   */
  void set2dElemToNodes( ArrayOfArrays< localIndex > && _2dElemToNodes );

  /**
   * @brief Defines the 2d elements to edges mapping.
   * @param _2dElemToEdges The input mapping.
   */
  void set2dElemToEdges( ArrayOfArrays< localIndex > && _2dElemToEdges );

  /**
   * @brief Defines the 2d elements to faces mapping.
   * @param _2dElemToFaces The input mapping.
   */
  void set2dElemToFaces( ArrayOfArrays< localIndex > && _2dElemToFaces );

  /**
   * @brief Defines the 2d faces to elements mapping.
   * @param _2dFaceTo2dElems The input mapping.
   */
  void set2dFaceTo2dElems( ArrayOfArrays< localIndex > && _2dFaceTo2dElems );

  /**
   * @brief Defines the 2d faces to edges mapping.
   * @param _2dFaceToEdge The input mapping.
   */
  void set2dFaceToEdge( array1d< localIndex > && _2dFaceToEdge );

  /**
   * @brief Defines the 2d elements to 3d elements mapping.
   * @param _2dElemToElems The input mapping.
   */
  void set2dElemToElems( ToCellRelation< ArrayOfArrays< localIndex > > && _2dElemToElems );

  void setLocalToGlobalMap( array1d< globalIndex > && l2g )
  { m_localToGlobalMap = l2g; }

  array1d< globalIndex > localToGlobalMap() const override
  { return m_localToGlobalMap; }

  ArrayOfArrays< globalIndex > getCollocatedNodes() const override;

  void setCollocatedNodes( ArrayOfArrays< globalIndex > && collocatedNodes );

private:

  localIndex m_num2dElements;
  localIndex m_num2dFaces;

  ArrayOfArrays< localIndex > m_2dElemToNodes;
  ArrayOfArrays< localIndex > m_2dElemToEdges;
  ArrayOfArrays< localIndex > m_2dElemToFaces;
  ToCellRelation< ArrayOfArrays< localIndex > > m_2dElemToElems;

  ArrayOfArrays< localIndex > m_2dFaceTo2dElems;
  array1d< localIndex > m_2dFaceToEdge;

  array1d< globalIndex > m_localToGlobalMap;

  ArrayOfArrays< globalIndex > m_collocatedNodes;
};


}

#endif // include guard
