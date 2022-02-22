/**
 * @file ParallelTopologyChange.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SURFACEGENERATION_PARALLELTOPOLOGYCHANGE_HPP_
#define GEOSX_PHYSICSSOLVERS_SURFACEGENERATION_PARALLELTOPOLOGYCHANGE_HPP_

#include "physicsSolvers/surfaceGeneration/SurfaceGenerator.hpp"

namespace geosx
{
class MeshLevel;
class NeighborCommunicator;
struct ModifiedObjectLists;

namespace parallelTopologyChange
{

void synchronizeTopologyChange( MeshLevel * const mesh,
                                std::vector< NeighborCommunicator > & neighbors,
                                ModifiedObjectLists & modifiedObjects,
                                ModifiedObjectLists & receivedObjects,
                                int mpiCommOrder );

}

}

#endif /* GEOSX_PHYSICSSOLVERS_SURFACEGENERATION_PARALLELTOPOLOGYCHANGE_HPP_ */
