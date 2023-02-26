/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file MapMeshLevels.hpp
 */

#ifndef GEOSX_MESH_MESHLEVEL_HPP_
#define GEOSX_MESH_MESHLEVEL_HPP_

namespace geosx{

// GEOSX_HOST_DEVICE inline
// R1Tensor const coarseToFineStructuredElemMap(localIndex const coarseElemIndex,
//                                              localIndex const coarseNx, 
//                                              localIndex const coarseNy,
//                                              localIndex const coarseNz,
//                                              localIndex const fineRx,
//                                              localIndex const fineRy,
//                                              localIndex const fineRz)   
// {
//     localIndex X = coarseElemIndex/(coarseNy*coarseNz);
//     localIndex Y = (coarseElemIndex - u*(coarseNy*coarseNz)) / coarseNz;
//     localIndex Z = (coarseElemIndex - u*(coarseNy*coarseNz) - v * coarseNz);
//     //this only returns the triplet for the first element 
//     //to generate all, the equation is the following
//     //0<= u < fineRx
//     //0<= v < fineRy
//     //0<= w < fineRz
//     //fineIndex = (X*fineRx + u)*fineRy*coarseNy*fineRz*coarseNz 
//     // + (Y*fineRy + v)*fineRz*coarseNz + (Z*fineRz + w)
//     //vary u,v,w
//     R1Tensor fineRefElem = {X*fineRx, Y*fineRy, Z*fineRz};
//     return fineRefElem;
// } 

// GEOSX_HOST_DEVICE inline
// localIndex fineToCoarseStructuredElemMap(localIndex const fineElemIndex,
//                                                localIndex const coarseNx, 
//                                                localIndex const coarseNy,
//                                                localIndex const coarseNz,
//                                                localIndex const fineRx,
//                                                localIndex const fineRy,
//                                                localIndex const fineRz)   
// {
//     localIndex coarseElemIndex = 0;
//     localIndex x = fineElemIndex/(fineRy*coarseNy*fineRz*coarseNz);
//     localIndex y = (fineElemIndex - x * (fineRy*coarseNy*fineRz*coarseNz)) / (fineRz*coarseNz);
//     localIndex z = fineElemIndex - x * (fineRy*coarseNy*fineRz*coarseNz) - y * (fineRz*coarseNz);
//     coarseElemIndex = (x/fineRx) * coarseNy * coarseNz + (y/fineRy) * coarseNz + (z/fineRz); 
//     return coarseElemIndex;    
// } 

GEOSX_HOST_DEVICE inline
localIndex const coarseToFineStructuredNodeMap(localIndex const coarseNodeIndex
                                               localIndex const coarseNx, 
                                               localIndex const coarseNy,
                                               localIndex const coarseNz,
                                               localIndex const fineNx,
                                               localIndex const fineNx,
                                               localIndex const fineNx)   
{

    localIndex fineNodeIndex = 0;

    return fineNodeIndex;
}

GEOSX_HOST_DEVICE inline
localIndex const fineToCoarseStructuredNodeMap(localIndex const coarseNodeIndex
                                               localIndex const coarseNx, 
                                               localIndex const coarseNy,
                                               localIndex const coarseNz,
                                               localIndex const fineNx,
                                               localIndex const fineNx,
                                               localIndex const fineNx)   
{

    localIndex coarseNodeIndex = 0;

    return coarseNodeIndex;
}

} /* namespace geosx */

#endif /* GEOSX_MESH_MESHLEVEL_HPP_ */