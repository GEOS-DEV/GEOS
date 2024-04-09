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

#ifndef GEOS_EMBEDDEDSURFACEBLOCK_HPP
#define GEOS_EMBEDDEDSURFACEBLOCK_HPP

#include "EmbeddedSurfaceBlockABC.hpp"

namespace geos
{
    
/**
 * @brief Simple implementation of the @p EmbeddedSurfaceBlockABC contract
 * 
*/

class EmbeddedSurfaceBlock: public EmbeddedSurfaceBlockABC
{
public:
  /**
   * @brief Constructor.
   * @param[in] name Name of this EmbeddedSurfaceBlock.
   * @param[in] parent Parent group.
   */
  EmbeddedSurfaceBlock( string const & name,
             Group * const parent )
    :
    EmbeddedSurfaceBlockABC( name, parent )
  { }

    localIndex numEmbeddedSurfElem() const override;
    ArrayOfArrays<localIndex> getEmbeddedSurfElemToNodes() const override;
    ArrayOfArrays<localIndex> getEmbeddedSurfElemTo3dElem() const override;
    ArrayOfArrays<real64> getEmbeddedSurfElemNodes() const override;
    ArrayOfArrays<real64> getEmbeddedSurfElemNormalVectors() const override;
    ArrayOfArrays<real64> getEmbeddedSurfElemTangentialLengthVectors() const override;
    ArrayOfArrays<real64> getEmbeddedSurfElemTangentialWidthVectors() const override;
    

    /**
     * @brief Sets the number of embedded elements
     * @param _numEmbeddedSurfElem the input value
    */
    void setNumEmbeddedSurfElem(localIndex _numEmbeddedSurfElem);

    /**
     * @brief Sets the embedded elements to nodes mapping
     * @param _embeddedSurfElemToNodes the input mapping.
    */
    void setEmbeddedSurfElemToNodes(ArrayOfArrays<localIndex> && _embeddedSurfElemToNodes);

    /**
     * @brief Sets the embedded elements to 1 3d element mapping
     * @param _embeddedSurfElemTo3dElem the input mapping.
    */
    void setEmbeddedSurfElemTo3dElem(ArrayOfArrays<localIndex> && _embeddedSurfElemTo3dElem);

    /**
     * @brief Sets the embedded elements nodes coordinates
     * @param _embeddedSurfElemNodes the coordinates of embedded elements nodes.
    */
    void setEmbeddedSurfElemNodes(ArrayOfArrays<real64> && _embeddedSurfElemNodes);
    
    /**
     * @brief Sets the embedded elements normals 
     * @param _normals the coordinates of embedded elements normals.
    */
    void setEmbeddedSurfElemNormalVectors(ArrayOfArrays<real64> && _normals);

    /**
     * @brief Sets the embedded elements tangential length vectors (horizontal plan) 
     * @param _lengthVectors the coordinates of embedded elements tangential length vectors.
    */
    void setEmbeddedSurfElemTangentialLengthVectors(ArrayOfArrays<real64> && _lengthVectors);

    /**
     * @brief Sets the embedded elements tangential width vectors (vertical plan) 
     * @param _widthVectors the coordinates of embedded elements tangential width vectors.
    */
    void setEmbeddedSurfElemTangentialWidthVectors(ArrayOfArrays<real64> && _widthVectors);

private:
 
    localIndex m_numEmbeddedSurfaces;
    ArrayOfArrays< localIndex > m_embeddedSurfElemToNodes;
    ArrayOfArrays<localIndex> m_embeddedSurfElemTo3dElem;
    ArrayOfArrays< real64 > m_embeddedSurfElemNodes;
    ArrayOfArrays< real64 > m_embeddedSurfElemNormals;
    ArrayOfArrays< real64 > m_embeddedSurfElemLengthVectors;
    ArrayOfArrays< real64 > m_embeddedSurfElemWidthVectors;
};

}

#endif //inlcude guard