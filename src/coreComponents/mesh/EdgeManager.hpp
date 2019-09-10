/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
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
 * @file EdgeManager.hpp
 */

#ifndef EDGEMANAGERT_H_
#define EDGEMANAGERT_H_

#include "InterObjectRelation.hpp"
#include "managers/ObjectManagerBase.hpp"
#include "ToElementRelation.hpp"


namespace geosx
{
class FaceManager;
class NodeManager;
class CellBlockManager;

class EdgeManager : public ObjectManagerBase
{
public:

  using NodeMapType = FixedOneToManyRelation;
  using FaceMapType = UnorderedVariableOneToManyRelation;

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
               Group * const parent );
  ~EdgeManager() override;

//  void Initialize() {}

  void SetDomainBoundaryObjects( const ObjectDataStructureBaseT* const referenceObject = nullptr);
  void SetIsExternal( FaceManager const * const faceManager );
//  void ExtractMapFromObjectForAssignGlobalObjectNumbers( const ObjectDataStructureBaseT& compositionObjectManager,
//                                                         array<globalIndex_array>& objectToCompositionObject );

  void BuildEdges( FaceManager * const faceManager, NodeManager * const nodeManager );


  virtual void
  ExtractMapFromObjectForAssignGlobalIndexNumbers( ObjectManagerBase const * const  nodeManager,
                                                   std::vector< std::vector< globalIndex > >& faceToNodes ) override final;

  virtual localIndex PackUpDownMapsSize( arrayView1d<localIndex const> const & packList ) const override;
  virtual localIndex PackUpDownMaps( buffer_unit_type * & buffer,
                                     arrayView1d<localIndex const> const & packList ) const override;

  virtual localIndex UnpackUpDownMaps( buffer_unit_type const * & buffer,
                                       localIndex_array & packList,
                                       bool const overwriteUpMaps,
                                       bool const overwriteDownMaps ) override;

  void FixUpDownMaps( bool const clearIfUnmapped );

  void depopulateUpMaps( std::set<localIndex> const & receivedEdges,
                         array1d< array1d< localIndex > > const & facesToEdges );

  void ConnectivityFromGlobalToLocal( const set<localIndex>& indices,
                                      const map<globalIndex,localIndex>& nodeGlobalToLocal,
                                      const map<globalIndex,localIndex>& faceGlobalToLocal );

//  void UpdateEdgeExternalityFromSplit( const FaceManager& faceManager,
//                                     const set<localIndex>& newEdgeIndices,
//                                     const set<localIndex>& modifiedEdgeIndices );

  void AddToEdgeToFaceMap( const FaceManager * faceManager,
                           const localIndex_array& newFaceIndices );

  void SplitEdge( const localIndex indexToSplit,
                  const localIndex parentNodeIndex,
                  const localIndex childNodeIndex[2],
                  array1d<set<localIndex>>& nodesToEdges );

  bool hasNode( const localIndex edgeID, const localIndex nodeID ) const;

  void calculateCenter( localIndex const edgeIndex,
                        arraySlice1d<R1Tensor const> const & X,
                        R1Tensor & center ) const;

  void calculateLength( localIndex const edgeIndex,
                        arraySlice1d<R1Tensor const> const & X,
                        R1Tensor & center ) const;



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


    static constexpr auto edgesTofractureConnectorsEdgesString = "edgesToFractureConnectors";
    static constexpr auto fractureConnectorEdgesToEdgesString = "fractureConnectorsToEdges";
    static constexpr auto fractureConnectorsEdgesToFaceElementsIndexString = "fractureConnectorsToElementIndex";


    dataRepository::ViewKey nodesList             = { nodeListString };
    dataRepository::ViewKey faceList              = { faceListString };
    dataRepository::ViewKey elementRegionList     = { elementRegionListString };
    dataRepository::ViewKey elementSubRegionList  = { elementSubRegionListString };
    dataRepository::ViewKey elementList           = { elementListString };

  } viewKeys;


  struct groupKeyStruct : ObjectManagerBase::groupKeyStruct
  {} groupKeys;

  constexpr int maxEdgesPerNode() const { return 100; }

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


  // TODO These should be in their own subset of edges when we add that capability.
  /// maps from the edges to the fracture connectors index (edges that are fracture connectors)
  set< localIndex > m_recalculateFractureConnectorEdges;
  map< localIndex, localIndex > m_edgesToFractureConnectorsEdges;
  array1d<localIndex> m_fractureConnectorsEdgesToEdges;
  ArrayOfArrays<localIndex> m_fractureConnectorEdgesToFaceElements;


private:

  NodeMapType m_toNodesRelation;
  FaceMapType m_toFacesRelation;

  map< localIndex, array1d<globalIndex> > m_unmappedGlobalIndicesInToNodes;
  map< localIndex, set<globalIndex> > m_unmappedGlobalIndicesInToFaces;




  template<bool DOPACK>
  localIndex PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                                    arrayView1d<localIndex const> const & packList ) const;

};

inline void EdgeManager::calculateCenter( localIndex const edgeIndex,
                                          arraySlice1d<R1Tensor const> const & X,
                                          R1Tensor & center ) const
{
  center = X[m_toNodesRelation[edgeIndex][0]];
  center += X[m_toNodesRelation[edgeIndex][1]];
  center *= 0.5;
}

inline void EdgeManager::calculateLength( localIndex const edgeIndex,
                                          arraySlice1d<R1Tensor const> const & X,
                                          R1Tensor & center ) const
{
  center = X[m_toNodesRelation[edgeIndex][1]];
  center -= X[m_toNodesRelation[edgeIndex][0]];
}

}
#endif /* EDGEMANAGERT_H_ */
