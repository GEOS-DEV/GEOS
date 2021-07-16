/**
 * @file ParallelTopologyChange.hpp
 */

#ifndef SRC_EXTERNALCOMPONENTS_HYDROFRACTURE_PARALLELTOPOLOGYCHANGE_HPP_
#define SRC_EXTERNALCOMPONENTS_HYDROFRACTURE_PARALLELTOPOLOGYCHANGE_HPP_

#include "physicsSolvers/surfaceGeneration/SurfaceGenerator.hpp"

namespace geosx
{
class MeshLevel;
class NeighborCommunicator;
struct ModifiedObjectLists;

class ParallelTopologyChange
{
public:
  ParallelTopologyChange();
  ~ParallelTopologyChange();

  static void synchronizeTopologyChange( MeshLevel * const mesh,
                                         std::vector< NeighborCommunicator > & neighbors,
                                         ModifiedObjectLists & modifiedObjects,
                                         ModifiedObjectLists & receivedObjects,
                                         int mpiCommOrder );

  static void packNewAndModifiedObjectsToOwningRanks( NeighborCommunicator * const neighbor,
                                                      MeshLevel * const meshLevel,
                                                      ModifiedObjectLists const & modifiedObjects,
                                                      int const commID );

  static localIndex
  unpackNewAndModifiedObjectsOnOwningRanks( NeighborCommunicator * const neighbor,
                                            MeshLevel * const mesh,
                                            int const commID,
                                            ModifiedObjectLists & receivedObjects );

  static void packNewModifiedObjectsToGhosts( NeighborCommunicator * const neighbor,
                                              int commID,
                                              MeshLevel * const mesh,
                                              ModifiedObjectLists & receivedObjects );

  static void unpackNewModToGhosts( NeighborCommunicator * const neighbor,
                                    int commID,
                                    MeshLevel * const mesh,
                                    ModifiedObjectLists & receivedObjects );


  static void updateConnectorsToFaceElems( std::set< localIndex > const & newFaceElements,
                                           FaceElementSubRegion const & faceElemSubRegion,
                                           EdgeManager & edgeManager );

};

} /* namespace lvarray */

#endif /* SRC_EXTERNALCOMPONENTS_HYDROFRACTURE_PARALLELTOPOLOGYCHANGE_HPP_ */
