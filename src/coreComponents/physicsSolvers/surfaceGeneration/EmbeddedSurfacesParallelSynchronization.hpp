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
 * @file EmebeddedSurfacesParallelSynchronization.hpp
 */

#ifndef SRC_COMPONENTS_SURFACEGENERATION_EMBEDDEDSURFACESPARALLELSYNCHRONIZATION_HPP_
#define SRC_COMPONENTS_SURFACEGENERATION_EMBEDDEDSURFACESPARALLELSYNCHRONIZATION_HPP_

#include "physicsSolvers/surfaceGeneration/EmbeddedSurfaceGenerator.hpp"

namespace geosx
{
class MeshLevel;
class NeighborCommunicator;
struct ModifiedObjectLists;

class EmebeddedSurfacesParallelSynchronization
{
public:

  EmebeddedSurfacesParallelSynchronization() = delete;

  ~EmebeddedSurfacesParallelSynchronization() = delete;

  EmebeddedSurfacesParallelSynchronization( EmebeddedSurfacesParallelSynchronization const & ) = delete;

  EmebeddedSurfacesParallelSynchronization( EmebeddedSurfacesParallelSynchronization && ) = delete;

  EmebeddedSurfacesParallelSynchronization & operator=( EmebeddedSurfacesParallelSynchronization const & ) = delete;

  EmebeddedSurfacesParallelSynchronization & operator=( EmebeddedSurfacesParallelSynchronization && ) = delete;

  static void synchronizeNewSurfaces( MeshLevel & mesh,
                                      std::vector< NeighborCommunicator > & neighbors,
                                      NewObjectLists & newObjects,
                                      int const mpiCommOrder );


  static void packNewObjectsToGhosts( NeighborCommunicator * const neighbor,
                                      int commID,
                                      MeshLevel & mesh,
                                      NewObjectLists & newObjects );

  static void unpackNewToGhosts( NeighborCommunicator * const neighbor,
                                 int commID,
                                 MeshLevel & mesh );


  static void synchronizeFracturedElements( MeshLevel & mesh,
                                            std::vector< NeighborCommunicator > & neighbors,
                                            string const fractureRegionName );

  static void packFracturedToGhosts( NeighborCommunicator * const neighbor,
                                     int commID,
                                     MeshLevel & mesh,
                                     string const fractureRegionName );

  static void unpackFracturedToGhosts( NeighborCommunicator * const neighbor,
                                       int commID,
                                       MeshLevel & mesh,
                                       string const fractureRegionName );

};

} /* namespace lvarray */

#endif /* SRC_COMPONENTS_SURFACEGENERATION_EMBEDDEDSURFACESPARALLELSYNCHRONIZATION_HPP_ */
