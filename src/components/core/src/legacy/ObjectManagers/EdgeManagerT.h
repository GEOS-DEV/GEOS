/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistrubute it and/or modify it under
 * the terms of the GNU Lesser General Public Liscense (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file EdgeManagerT.h
 * @author settgast1
 * @date Jun 22, 2011
 */

#ifndef EDGEMANAGERT_H_
#define EDGEMANAGERT_H_

//#include "DataStructures/VectorFields/ObjectDataStructureBaseT.h"
#include "../../dataRepository/Group.hpp"


class FaceManagerT;
class NodeManager;

class EdgeManagerT : public ObjectDataStructureBaseT
{
public:
  EdgeManagerT();
  ~EdgeManagerT();

  void Initialize() {}

  void SetDomainBoundaryObjects( const ObjectDataStructureBaseT* const referenceObject = NULL);
  void SetIsExternal( const ObjectDataStructureBaseT* const referenceObject = NULL);
  void ExtractMapFromObjectForAssignGlobalObjectNumbers( const ObjectDataStructureBaseT& compositionObjectManager,
                                                         array<gArray1d>& objectToCompositionObject );


  void BuildEdges( const FaceManagerT& faceManager, const NodeManager& nodeManager );
  void BuildEdges( const ElementManagerT& elementManager, const NodeManager& nodeManager );

  template< typename T_indices >
  unsigned int PackEdges( const T_indices& sendedges,
                          const NodeManager& nodeManager,
                          const FaceManagerT& faceManager,
                          bufvector& buffer,
                          const bool packConnectivityToGlobal,
                          const bool packFields,
                          const bool packMaps,
                          const bool packSets  ) const;

  unsigned int UnpackEdges( const char*& buffer,
                            const NodeManager& nodeManager,
                            const FaceManagerT& faceManager,
                            lArray1d& edgeReceiveLocalIndices,
                            const bool unpackConnectivityToLocal,
                            const bool unpackFields,
                            const bool unpackMaps,
                            const bool unpackSets  );

  void ConnectivityFromGlobalToLocal( const lSet& indices,
                                      const std::map<globalIndex,localIndex>& nodeGlobalToLocal,
                                      const std::map<globalIndex,localIndex>& faceGlobalToLocal );

//  void UpdateEdgeExternalityFromSplit( const FaceManagerT& faceManager,
//                                     const lSet& newEdgeIndices,
//                                     const lSet& modifiedEdgeIndices );

  void EdgeCenter(const NodeManager& nodeManager, localIndex edge, R1Tensor& center) const;
  void EdgeVector(const NodeManager& nodeManager, localIndex edge, R1Tensor& vector) const;
  realT EdgeLength(const NodeManager& nodeManager, localIndex edge) const;

  void AddToEdgeToFaceMap( const FaceManagerT& faceManager,
                           const lArray1d& newFaceIndices );

  void SplitEdge( const localIndex indexToSplit,
                  const localIndex parentNodeIndex,
                  const localIndex childNodeIndex[2],
                  array<lSet>& nodesToEdges );

  bool hasNode( const localIndex edgeID, const localIndex nodeID ) const;

  localIndex FindEdgeFromNodeIDs(const localIndex nodeA, const localIndex nodeB, const NodeManager& nodeManager);

  void SetLayersFromDomainBoundary(const NodeManager& nodeManager);


  FixedOneToManyRelation& m_toNodesRelation;
  UnorderedVariableOneToManyRelation& m_toFacesRelation;

};

#endif /* EDGEMANAGERT_H_ */
