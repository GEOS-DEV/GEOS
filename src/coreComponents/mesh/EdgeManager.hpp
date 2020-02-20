/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file EdgeManager.hpp
 */

#ifndef GEOSX_MESH_EDGEMANAGER_HPP_
#define GEOSX_MESH_EDGEMANAGER_HPP_

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
  using FaceMapType = InterObjectRelation< ArrayOfSets< localIndex > >;

  /**
    * @name Static Factory Catalog Functions
    */
   ///@{

   static const string CatalogName()
   { return "EdgeManager"; }

   virtual const string getCatalogName() const override final
   { return EdgeManager::CatalogName(); }

   ///@}

   static inline localIndex GetFaceMapOverallocation()
   { return 4; }

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

  void CompressRelationMaps( );

  void depopulateUpMaps( std::set<localIndex> const & receivedEdges,
                         ArrayOfArraysView< localIndex const > const & facesToEdges );

  void ConnectivityFromGlobalToLocal( const SortedArray<localIndex>& indices,
                                      const map<globalIndex,localIndex>& nodeGlobalToLocal,
                                      const map<globalIndex,localIndex>& faceGlobalToLocal );

//  void UpdateEdgeExternalityFromSplit( const FaceManager& faceManager,
//                                     const SortedArray<localIndex>& newEdgeIndices,
//                                     const SortedArray<localIndex>& modifiedEdgeIndices );

  void AddToEdgeToFaceMap( FaceManager const * const faceManager,
                           arrayView1d< localIndex const > const & newFaceIndices );

  void SplitEdge( const localIndex indexToSplit,
                  const localIndex parentNodeIndex,
                  const localIndex childNodeIndex[2],
                  array1d<SortedArray<localIndex>>& nodesToEdges );

  bool hasNode( const localIndex edgeID, const localIndex nodeID ) const;

  R1Tensor calculateCenter( localIndex const edgeIndex,
                            arrayView2d<real64 const, nodes::REFERENCE_POSITION_USD> const & X ) const;

  R1Tensor calculateLength( localIndex const edgeIndex,
                            arrayView2d<real64 const, nodes::REFERENCE_POSITION_USD> const & X ) const;



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

  NodeMapType       & nodeList()       { return m_toNodesRelation; }
  NodeMapType const & nodeList() const { return m_toNodesRelation; }

  localIndex & nodeList( localIndex const edgeIndex, localIndex const nodeIndex )
  {
    return m_toNodesRelation(edgeIndex, nodeIndex);
  }
  localIndex nodeList( localIndex const edgeIndex, localIndex const nodeIndex ) const
  {
    return m_toNodesRelation(edgeIndex, nodeIndex);
  }


  FaceMapType       & faceList()       { return m_toFacesRelation; }
  FaceMapType const & faceList() const { return m_toFacesRelation; }


  // TODO These should be in their own subset of edges when we add that capability.
  /// maps from the edges to the fracture connectors index (edges that are fracture connectors)
  SortedArray< localIndex > m_recalculateFractureConnectorEdges;
  map< localIndex, localIndex > m_edgesToFractureConnectorsEdges;
  array1d<localIndex> m_fractureConnectorsEdgesToEdges;
  ArrayOfArrays<localIndex> m_fractureConnectorEdgesToFaceElements;


private:
  NodeMapType m_toNodesRelation;
  FaceMapType m_toFacesRelation;

  map< localIndex, array1d<globalIndex> > m_unmappedGlobalIndicesInToNodes;
  map< localIndex, SortedArray<globalIndex> > m_unmappedGlobalIndicesInToFaces;


  template<bool DOPACK>
  localIndex PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                                    arrayView1d<localIndex const> const & packList ) const;

};

inline R1Tensor EdgeManager::calculateCenter( localIndex const edgeIndex,
                                              arrayView2d<real64 const, nodes::REFERENCE_POSITION_USD> const & X ) const
{
  R1Tensor center = X[m_toNodesRelation[edgeIndex][0]];
  center += X[m_toNodesRelation[edgeIndex][1]];
  center *= 0.5;
  return center;
}

inline R1Tensor EdgeManager::calculateLength( localIndex const edgeIndex,
                                              arrayView2d<real64 const, nodes::REFERENCE_POSITION_USD> const & X ) const
{
  R1Tensor length = X[m_toNodesRelation[edgeIndex][1]];
  length -= X[m_toNodesRelation[edgeIndex][0]];
  return length;
}

}
#endif /* GEOSX_MESH_EDGEMANAGER_HPP_ */
