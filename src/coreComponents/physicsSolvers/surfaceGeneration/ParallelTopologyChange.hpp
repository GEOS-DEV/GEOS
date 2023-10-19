/**
 * @file ParallelTopologyChange.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_SURFACEGENERATION_PARALLELTOPOLOGYCHANGE_HPP_
#define GEOS_PHYSICSSOLVERS_SURFACEGENERATION_PARALLELTOPOLOGYCHANGE_HPP_

#include "physicsSolvers/surfaceGeneration/SurfaceGenerator.hpp"

namespace geos
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

#endif /* GEOS_PHYSICSSOLVERS_SURFACEGENERATION_PARALLELTOPOLOGYCHANGE_HPP_ */
