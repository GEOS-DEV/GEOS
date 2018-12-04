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
  ~EdgeManager() override;

//  void Initialize() {}

  void SetDomainBoundaryObjects( const ObjectDataStructureBaseT* const referenceObject = nullptr);
  void SetIsExternal( FaceManager const * const faceManager );
//  void ExtractMapFromObjectForAssignGlobalObjectNumbers( const ObjectDataStructureBaseT& compositionObjectManager,
//                                                         array<globalIndex_array>& objectToCompositionObject );

  void BuildEdges( FaceManager * const faceManager, NodeManager * const nodeManager );


  void ExtractMapFromObjectForAssignGlobalIndexNumbers( ObjectManagerBase const * const nodeManager,
                                                        array1d<globalIndex_array>& edgesToNodes ) override final;

  virtual localIndex PackUpDownMapsSize( arrayView1d<localIndex const> const & packList ) const override;
  virtual localIndex PackUpDownMaps( buffer_unit_type * & buffer,
                                     arrayView1d<localIndex const> const & packList ) const override;

  virtual localIndex UnpackUpDownMaps( buffer_unit_type const * & buffer,
                                       localIndex_array & packList ) override;

  virtual void FixUpDownMaps() override final;

  void ConnectivityFromGlobalToLocal( const set<localIndex>& indices,
                                      const std::map<globalIndex,localIndex>& nodeGlobalToLocal,
                                      const std::map<globalIndex,localIndex>& faceGlobalToLocal );

//  void UpdateEdgeExternalityFromSplit( const FaceManager& faceManager,
//                                     const set<localIndex>& newEdgeIndices,
//                                     const set<localIndex>& modifiedEdgeIndices );

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
                  array1d<set<localIndex>>& nodesToEdges );

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

  localIndex & nodeList( localIndex const edgeIndex, localIndex const nodeIndex )
  {
    return m_toNodesRelation(edgeIndex, nodeIndex);
  }
  localIndex nodeList( localIndex const edgeIndex, localIndex const nodeIndex ) const
  {
    return m_toNodesRelation(edgeIndex, nodeIndex);
  }


  UnorderedVariableOneToManyRelation       & faceList()       { return m_toFacesRelation; }
  UnorderedVariableOneToManyRelation const & faceList() const { return m_toFacesRelation; }


private:
  FixedOneToManyRelation m_toNodesRelation;
  UnorderedVariableOneToManyRelation m_toFacesRelation;

  map< localIndex, array1d<globalIndex> > m_unmappedGlobalIndicesInToNodes;
  map< localIndex, set<globalIndex> > m_unmappedGlobalIndicesInToFaces;


  template<bool DOPACK>
  localIndex PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                                    arrayView1d<localIndex const> const & packList ) const;

};
}
#endif /* EDGEMANAGERT_H_ */
