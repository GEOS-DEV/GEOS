/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2020-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOSX_FACEBLOCKABC_HPP
#define GEOSX_FACEBLOCKABC_HPP

#include "CellBlockUtilities.hpp"
#include "dataRepository/Group.hpp"
#include "common/DataTypes.hpp"

namespace geosx
{

/**
 * @brief Block of 2d elements (geometrical surfaces in 3d).
 *
 * @details The @p FaceBlockABC represents a zone of 2d (@e i.e. surfacic) elements.
 * Unlike its volumic equivalent @p CellBlockABC, @p FaceBlockABC is @e not @e homogeneous
 * and may hold 2d elements of multiple types (triangles, quadrangles...).
 * @details In this class, we'll use the term @e 2d @e element for the elements of the @p FaceBlockABC,
 * which are geometrical surfaces (in 3d).
 * In the same way, we'll use the wording @e 2d @e face
 * to refer to the 2d boundaries of the @e 2d @e elements.
 * The @e 2d @e face are geometrical segments (in 3d).
 */
class FaceBlockABC : public dataRepository::Group
{
public:
  FaceBlockABC( string const & name,
                Group * const parent ):
    Group( name, parent )
  { }

  /**
   * @brief Get the number of 2d elements (geometrical surfaces in 3d).
   * @return Number of 2d elements in the current face block.
   */
  virtual localIndex num2dElements() const = 0;

  /**
   * @brief Get the number of 2d faces (geometrical segments in 3d).
   * @return Number of 2d faces in the current face block.
   * @deprecated
   */
  virtual localIndex num2dFaces() const = 0;

  /**
   * @brief Get the nodes of each 2d element (geometrical surfaces in 3d).
   * @return The mapping of first dimension @p num2dElements.
   * @note Element numbering is local to the FaceBlockABC. Node numbering is local to the rank.
   */
  virtual ArrayOfArrays< localIndex > get2dElemToNodes() const = 0;

  /**
   * @brief Get the 3d edges of each 2d elements (geometrical surfaces in 3d).
   * @return The mapping of first dimension @p num2dElements.
   * @note Element numbering is local to the FaceBlockABC. Edges numbering is local to the rank.
   */
  virtual ArrayOfArrays< localIndex > get2dElemToEdges() const = 0;

  /**
   * @brief Get the 3d faces (of the volumic mesh) that are aside each 2d element (geometrical surfaces in 3d) of the @p FaceBlockABC.
   * @return A mapping of size @p num2dFaces * 2 (each 2d element matches two 3d faces).
   * @note Mesh is supposed to be conformal, so each 2d element touches two faces (or one on boundaries).
   * In case of one unique touching face, second missing value with we @p -1.
   * @note Element numbering is local to the FaceBlockABC. Faces numbering is local to the rank.
   */
  virtual array2d< localIndex > get2dElemToFaces() const = 0;

  /**
   * @brief Get the 3d elements that are aside each 2d element (geometrical surfaces in 3d) of the @p FaceBlockABC.
   * @return The mapping of dimension @p num2dElements * 2.
   * @note Mesh is supposed to be conformal, so each 2d element touches two 3d elements (or one on boundaries).
   * In case of one unique touching face, second missing value with we @p -1.
   * @note 2d element numbering is local to the FaceBlockABC. 3d element numbering is local to the CellBlockABC.
   */
  virtual ToCellRelation< array2d< localIndex > > get2dElemToElems() const = 0;

  /**
   * @brief Get the 3d edge of each 2d faces (geometrical edges in 3d).
   * @return The mapping of size @p num2dFaces.
   * @note 2d face numbering is local to the FaceBlockABC. 3d edge numbering is local to the rank.
   *
   * The current @p FaceBlockABC has its own 2d faces (which are segments in terms of 3d geometry).
   * Those segments are edges in the volumic mesh.
   * The current mapping informs which edge (in the volumic mesh) is the geometrical equivalent of the 2d face of the @p FaceBlockABC.
   */
  virtual array1d< localIndex > get2dFaceToEdge() const = 0;

  /**
   * @brief Get the 2d element(s) (geometrical surfaces in 3d) connected to each 2d face (geometrical 3d segment).
   * @return The mapping of first dimension @p num2dFaces.
   * @note 2d face and 2d element numberings are both local to the @p FaceBlockABC.
   */
  virtual ArrayOfArrays< localIndex > get2dFaceTo2dElems() const = 0;
};

}

#endif // include guard
