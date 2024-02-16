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

#ifndef GEOS_EMBEDDEDSURFACEBLOCKABC_HPP
#define GEOS_EMBEDDEDSURFACEBLOCKABC_HPP


#include "CellBlockUtilities.hpp"
#include "dataRepository/Group.hpp"
#include "common/DataTypes.hpp"


namespace geos
{

/**
 * @brief An embedded 2d element within a 3d element
 * 
 * @details The @p EmbeddedSurfaceBlockABC represents an array of 2d elements (@e surfacic) embedded within 
 * a 3d elements grid (2d element can intersect one or more 3d element). The 2d element is assumed to be a quad.  
 * @details In this class, we'll use the term @e 2d @e element for the elements of the @p EmbeddedSurfaceBlockABC,
 * which are geometrical quad surfaces (in 3d).
 * In the same way, we'll use the wording @e 2d @e face
 * to refer to the 2d boundaries of the @e 2d @e elements.
 * The @e 2d @e face are geometrical segments (in 3d).
 * 
*/
class EmbeddedSurfaceBlockABC : public dataRepository::Group
{
public:
 /**
  * @brief Constructor
  * @param name The name of this Group
  * @param parent The parent Group
 */
    
    EmbeddedSurfaceBlockABC(string const & name, 
                         Group * const parent):
      Group(name, parent)
    { }

   /**
    * @brief Get the number of embedded surface elements
    * @return Number of embedded surface elements
    * @details Return the number of embedded surface elements, each surface element
    * can intersect 1 or N 3d elements
    */ 
    virtual localIndex numEmbeddedSurfElem() const = 0;
    
    /**
    * @brief Get the indices of the nodes of all embedded surface elements
    * @return The mapping of first dimension @p numEmbeddedSurfElem.
    * Second dimension is 4 (quad), and represents the indices of each embedded surface element.
    *
    * @details each embedded surface element is supposed to have 4 nodes. Node indices for all the embedded
    * surface elements are given by getEmbeddedSurfaceElementsNodes. This method returns the mapping between
    * an embedded surface element and the 4 indices of its nodes.   
    */
    virtual ArrayOfArrays<localIndex> getEmbeddedSurfElemToNodes() const = 0;
    
    /**
    * @brief Get the indices of the parent 3d elements of each embedded surface element (1 or more)
    * @return The mapping of first dimension @p numEmbeddedSurfElem.
    * Second dimension 1 to number of 3d elements, and depends on how many 3d elements does the embedded surface element
    * intersect.
    *
    * @details each embedded surface element intersects 1 or more 3d elements. Indices of these 3d elements
    * are returned for each embedded surface element
    */
    virtual ToCellRelation<ArrayOfArrays< localIndex >> getEmbeddedSurfElemTo3dElem() const = 0;


    /**
    * @brief Get the x, y, and z coordinates of the embedded surface elements nodes. 
    * @return An array of x, y, z coordinates of all the embedded surface elements nodes.
    * first dimension is @p numEmbeddedSurfaceElements*4.
    * Second dimension is 3, and represents the x, y and z.
    */
    virtual ArrayOfArrays<real64> getEmbeddedSurfElemNodes() const = 0;
};
    
}

#endif // include guard