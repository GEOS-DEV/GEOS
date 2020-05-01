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


/**
 * @class EdgeManager
 * @brief This class provides an interface to ObjectManagerBase in order to manage edge data.
 */

class EdgeManager : public ObjectManagerBase
{
public:

  using NodeMapType = FixedOneToManyRelation;
  using FaceMapType = InterObjectRelation< ArrayOfSets< localIndex > >;

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  
  /**
   * @return the string representing the edge manager name in the catalog
   */
  static const string CatalogName()
  { return "EdgeManager"; }

  
  /**
   * @brief Getter used to acces the edge manager catalog name.
   * @return the edge manager catalog name
   */
  virtual const string getCatalogName() const override
  { return EdgeManager::CatalogName(); }

  ///@}

  /// Oversize the Face mapping by this amount for each edge (hardcoded)
  static localIndex faceMapExtraSpacePerEdge()
  { return 4; }

  /**
   * @name Constructors/destructos
   */
  ///@{
  
  /**
   * @brief main constructor for EdgeManager Objects
   * @param[in] name the name of the EdgeManager object in the repository
   * @param[in] parent the parent group of the EdgeManager object defined
   */
  EdgeManager( std::string const & name,
               Group * const parent );

  
  /**
   * @brief default destructor
   */
  ~EdgeManager() override;
  ///@}
  
  
  /**
   * @brief Resize the EdgeManager object and all it members.
   * @param[in] newsize the new number of edges.
   */
  virtual void resize( localIndex const newsize ) override;

  
  /**
   * @brief Set the node of the domain boundary object.
   * @param[in] referenceObject the reference of the face manager.
   */
  void SetDomainBoundaryObjects( const ObjectDataStructureBaseT * const referenceObject = nullptr );

  
  /**
   * @brief Set external edges.
   * @param[in] faceManager the face manager to obtain external face
   */
  void SetIsExternal( FaceManager const * const faceManager );

  
  /**
   * @brief Build faces-to-edges and nodes-to-edges relation maps.
   * @param[in] faceManager manager of all faces in a mesh
   * @param[in] nodeManager manager of all nodes in a mesh
   */
  void BuildEdges( FaceManager * const faceManager, NodeManager * const nodeManager );

  
  /**
   * @brief Extract the map from nodeManager object to assign the global index numbers of edges.
   * @param[in] nodeManager the nodeManager object.
   * @param[out] globalEdgeNodes the globalEdgesNodes array to be assigned
   */
  virtual void
  ExtractMapFromObjectForAssignGlobalIndexNumbers( ObjectManagerBase const * const nodeManager,
                                                   std::vector< std::vector< globalIndex > > & faceToNodes ) override;

  
  /**
   * @brief Calculate the future size of a packed list.
   * @param[in] packList the list of edge indices to be packed
   * @return a local index representing the of the future paked list
   */
  virtual localIndex PackUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const override;

  
  /**
   * @brief Pack an array of edge indices in a buffer.
   * @param[inout] buffer a buffer to pack the edge index data into
   * @param[in] packList thes indices of edge to pack
   * @return a localIndex type value representing the size of the packed array
   */
  virtual localIndex PackUpDownMaps( buffer_unit_type * & buffer,
                                     arrayView1d< localIndex const > const & packList ) const override;

  
  /**
   * @brief Unpack a buffer to an array of indices that has been packed.
   * @param[in] buffer buffer containing the packed indices
   * @param[inout] packList the array of local index to be unpacked
   * @param[in] overwriteUpMaps boolean to state if the Up maps should be overwritten
   * @param[in] overwriteDownsMaps boolean to state if the Downs maps should be overwritten
   * @return a localIndex type value representing the size of the unpacked array
   */
  virtual localIndex UnpackUpDownMaps( buffer_unit_type const * & buffer,
                                       localIndex_array & packList,
                                       bool const overwriteUpMaps,
                                       bool const overwriteDownMaps ) override;

  
  /**
   * @brief Call FixUpDownMaps of the class ObjectManagerBase for nodes-to-edges and nodes-to-faces relation maps.
   * @param[in] clearIfUnmaped boolean to indicate if the unmaped indices should be removed or not
   */
  void FixUpDownMaps( bool const clearIfUnmapped );

  
  /**
   * @brief Compress all nodes-to-faces relation maps.
   */
  void compressRelationMaps();

  
  /**
   * @brief Clean up the faces-to-edges mapping with respect a new list of edges.
   * @param[in] receivedEdges the new list of edges indices
   * @param[in] facesToEdges map to go from faces to edges
   */
  void depopulateUpMaps( std::set< localIndex > const & receivedEdges,
                         ArrayOfArraysView< localIndex const > const & facesToEdges );

  
  /**
   * @brief Build the mapping faces-to-edges relation from the mapping global to local nodes.
   * @param[in] indices array of index of the global to local nodes
   * @param[in] nodeGlobalToLocal map of the global to local nodes
   * @param[in] faceGlobalToLocal GEOX UNUSED PARAMETER
   */
  void ConnectivityFromGlobalToLocal( const SortedArray< localIndex > & indices,
                                      const map< globalIndex, localIndex > & nodeGlobalToLocal,
                                      const map< globalIndex, localIndex > & faceGlobalToLocal );

  
  /**
   * @brief Add new faces to the faces to edges map.
   * @param[in] faceManager the face manager
   * @param[in] newFaceIndices the new face indices to be added
   */
  void AddToEdgeToFaceMap( FaceManager const * const faceManager,
                           arrayView1d< localIndex const > const & newFaceIndices );

  
  /**
   * @brief Split an edge (separate its two extremitiy nodes)
   * @param[in] indexToSplit Index of the edge to split
   * @param[in] parentNodeIndex index of the parent node
   * @param[in] childNodeIndex index of the child node
   * @param[in] nodesToEdges  array of nodes-to-edges list
   */
  void SplitEdge( const localIndex indexToSplit,
                  const localIndex parentNodeIndex,
                  const localIndex childNodeIndex[2],
                  array1d< SortedArray< localIndex > > & nodesToEdges );

  
  /**
   * @brief Check if there is a relation between a node and an edge.
   * @param[in] edgeID local index of the edge
   * @param[in] nodeID local index of the node
   * @return boolean : true if the node #nodeID is part of the edge #edgeID; false otherwise
   */
  bool hasNode( const localIndex edgeID, const localIndex nodeID ) const;

  
  /**
   * @brief Calculate the center of an edge given its index.
   * @param[in] edgeIndex index of the edge
   * @param[in] X array view of nodes associated with the edge
   * @return coordinate of the edge center
   */
  R1Tensor calculateCenter( localIndex const edgeIndex,
                            arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X ) const;

  
  /**
   * @brief Calculate the length of an edge.
   * @param[in] edgeIndex index of the edge 
   * @param[in] X array view of nodes associated with the edge 
   * @return length of the edge
   */
  R1Tensor calculateLength( localIndex const edgeIndex,
                            arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X ) const;


  /**
   * @name viewKeyStruct/groupKeyStruct
   */
  ///@{

  /**
   * @struct Container of keys needed to acces the data of the class member
   */
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

  ///}@
  
  /**
   * @brief Return the  maximum number of edges per node.
   * @return Maximum allowable number of edges connected to one node (hardcoded for now)
   */
  constexpr int maxEdgesPerNode() const { return 100; }

  
  /**
   * @name Getters for stored value.
   */
  ///@{
  
  /**
   * @brief Get a node list.
   * @return reference to the list of edge to node relation
   */
  NodeMapType & nodeList()       { return m_toNodesRelation; }

  
  /**
   * @brief Provide an accessor to an immutable node list.
   * @return reference to the list of edge to node relation
   */
  NodeMapType const & nodeList() const { return m_toNodesRelation; }

  
  /**
   * @brief Get a reference to a node list
   * @param[in] edgeIndex local index of the edge
   * @param[in] nodeIndex local index of the node
   * @return a reference of the index of the edge to node relation
   */
  localIndex & nodeList( localIndex const edgeIndex, localIndex const nodeIndex )
  { return m_toNodesRelation( edgeIndex, nodeIndex ); }

  
  /**
   * @brief Get a node list
   * @param[in] edgeIndex local index of the edge
   * @param[in] nodeIndex local index of the node
   * @return a reference of the index of the edge to node relation
   */
  localIndex nodeList( localIndex const edgeIndex, localIndex const nodeIndex ) const
  { return m_toNodesRelation( edgeIndex, nodeIndex ); }

  
  /**
   * @brief Get the nodes-to-face relation
   * @return A FaceMapType
   */
  FaceMapType & faceList()       { return m_toFacesRelation; }

  
  /**
   * @brief Get an immutable accessor to the nodes-to-face relation
   * @return An immutable FaceMapType
   */
  FaceMapType const & faceList() const { return m_toFacesRelation; }

  ///@}

  // TODO These should be in their own subset of edges when we add that capability.
  /// map from the edges to the fracture connectors index (edges that are fracture connectors)
  SortedArray< localIndex > m_recalculateFractureConnectorEdges;
  map< localIndex, localIndex > m_edgesToFractureConnectorsEdges;
  array1d< localIndex > m_fractureConnectorsEdgesToEdges;
  ArrayOfArrays< localIndex > m_fractureConnectorEdgesToFaceElements;


private:

  /// Map for the nodes-to-eges relation
  NodeMapType m_toNodesRelation;

  /// Map for the faces-to-eges relation
  FaceMapType m_toFacesRelation;

  /// Unmaped edge indices (those that are not in the nodes-to-edges relation map)
  map< localIndex, array1d< globalIndex > > m_unmappedGlobalIndicesInToNodes;

  /// Unmaped edge indices (those that are not in the faces-to-edges relation map)
  map< localIndex, SortedArray< globalIndex > > m_unmappedGlobalIndicesInToFaces;

  
  /**
   * @brief function to pack the upward and downward pointing maps.
   * @tparam DOPACK template argument to determine whether or not to pack the buffer. If false, the buffer is not
   *                packed and the function returns the size of the packing that would have occured if set to TRUE.
   * @param buffer the buffer to pack data into
   * @param packList the indices of nodes that should be packed.
   * @return size of data packed in terms of number of chars
   */
  template< bool DOPACK >
  localIndex PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                                    arrayView1d< localIndex const > const & packList ) const;

};

/**
 *@brief Calculate the coordinate of the center of an edge.
 *@param [in] edgeIndex index of the edge to calculate the center of
 *@param [in] X  array with node coordinates
 *@return Coordinates of the edge center
 */
inline R1Tensor EdgeManager::calculateCenter( localIndex const edgeIndex,
                                              arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X ) const
{
  R1Tensor center = X[m_toNodesRelation[edgeIndex][0]];
  center += X[m_toNodesRelation[edgeIndex][1]];
  center *= 0.5;
  return center;
}


/**
 *@brief Calculate the length of an edge.
 *@param [in] edgeIndex index of the edge to calculate the center of
 *@param [in] X array with node coordinates
 *@return Coordinates of the edge center
 */
inline R1Tensor EdgeManager::calculateLength( localIndex const edgeIndex,
                                              arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X ) const
{
  R1Tensor length = X[m_toNodesRelation[edgeIndex][1]];
  length -= X[m_toNodesRelation[edgeIndex][0]];
  return length;
}

}
#endif /* GEOSX_MESH_EDGEMANAGER_HPP_ */
