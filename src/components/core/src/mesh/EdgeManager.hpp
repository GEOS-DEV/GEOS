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
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
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

#include "InterObjectRelation.hpp"
#include "managers/ObjectManagerBase.hpp"


namespace geosx
{
class FaceManager;
class NodeManager;
class CellBlockManager;

class EdgeManager : public ObjectManagerBase
{
public:

  /**
    * @name Static Factory Catalog Functions
    */
   ///@{

   static const string CatalogName()
   { return "EdgeManager"; }

   virtual const string getCatalogName() const override final
   { return EdgeManager::CatalogName(); }


   ///@}


  EdgeManager( std::string const & name,
               ManagedGroup * const parent );
  ~EdgeManager();

//  void Initialize() {}

  void SetDomainBoundaryObjects( const ObjectDataStructureBaseT* const referenceObject = nullptr);
  void SetIsExternal( const ObjectDataStructureBaseT* const referenceObject = nullptr);
  void ExtractMapFromObjectForAssignGlobalObjectNumbers( const ObjectDataStructureBaseT& compositionObjectManager,
                                                         array<globalIndex_array>& objectToCompositionObject );


  void BuildEdges( FaceManager * const faceManager, NodeManager * const nodeManager );

  template< typename T_indices >
  unsigned int PackEdges( const T_indices& sendedges,
                          const NodeManager& nodeManager,
                          const FaceManager& faceManager,
//                          bufvector& buffer,
                          const bool packConnectivityToGlobal,
                          const bool packFields,
                          const bool packMaps,
                          const bool packSets  ) const;

  unsigned int UnpackEdges( const char*& buffer,
                            const NodeManager& nodeManager,
                            const FaceManager& faceManager,
                            localIndex_array& edgeReceiveLocalIndices,
                            const bool unpackConnectivityToLocal,
                            const bool unpackFields,
                            const bool unpackMaps,
                            const bool unpackSets  );

  void ConnectivityFromGlobalToLocal( const lSet& indices,
                                      const std::map<globalIndex,localIndex>& nodeGlobalToLocal,
                                      const std::map<globalIndex,localIndex>& faceGlobalToLocal );

//  void UpdateEdgeExternalityFromSplit( const FaceManager& faceManager,
//                                     const lSet& newEdgeIndices,
//                                     const lSet& modifiedEdgeIndices );

//  void EdgeCenter(const NodeManager& nodeManager, localIndex edge, R1Tensor&
// center)const;
//  void EdgeVector(const NodeManager& nodeManager, localIndex edge, R1Tensor&
// vector)const;
//  realT EdgeLength(const NodeManager& nodeManager, localIndex edge) const;

  void AddToEdgeToFaceMap( const FaceManager * faceManager,
                           const localIndex_array& newFaceIndices );

  void SplitEdge( const localIndex indexToSplit,
                  const localIndex parentNodeIndex,
                  const localIndex childNodeIndex[2],
                  array<lSet>& nodesToEdges );

  bool hasNode( const localIndex edgeID, const localIndex nodeID ) const;

//  localIndex FindEdgeFromNodeIDs(const localIndex nodeA, const localIndex
// nodeB, const NodeManager& nodeManager);

//  void SetLayersFromDomainBoundary(const NodeManager * const nodeManager);


//  FixedOneToManyRelation& m_toNodesRelation;
//  UnorderedVariableOneToManyRelation& m_toFacesRelation;

  struct viewKeyStruct : ObjectManagerBase::viewKeyStruct
  {
    static constexpr auto nodeListString              = "nodeList";
    static constexpr auto faceListString              = "faceList";
    static constexpr auto elementRegionListString     = "elemRegionList";
    static constexpr auto elementSubRegionListString  = "elemSubRegionList";
    static constexpr auto elementListString           = "elemList";

    dataRepository::ViewKey nodesList             = { nodeListString };
    dataRepository::ViewKey faceList              = { faceListString };
    dataRepository::ViewKey elementRegionList     = { elementRegionListString };
    dataRepository::ViewKey elementSubRegionList  = { elementSubRegionListString };
    dataRepository::ViewKey elementList           = { elementListString };

  } viewKeys;


  struct groupKeyStruct : ObjectManagerBase::groupKeyStruct
  {} groupKeys;


  FixedOneToManyRelation       & nodeList()       { return m_toNodesRelation; }
  FixedOneToManyRelation const & nodeList() const { return m_toNodesRelation; }

  UnorderedVariableOneToManyRelation       & faceList()       { return m_toFacesRelation; }
  UnorderedVariableOneToManyRelation const & faceList() const { return m_toFacesRelation; }


private:
  FixedOneToManyRelation m_toNodesRelation;
  UnorderedVariableOneToManyRelation m_toFacesRelation;

};
}
#endif /* EDGEMANAGERT_H_ */
