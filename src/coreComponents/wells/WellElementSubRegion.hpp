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

#ifndef GEOSX_WELLS_WELLELEMENTSUBREGION_HPP_
#define GEOSX_WELLS_WELLELEMENTSUBREGION_HPP_

#include "mesh/ElementSubRegionBase.hpp"
#include "mesh/InterObjectRelation.hpp"
#include "PerforationData.hpp"

namespace geosx
{

/**
 * @class WellElementSubRegion
 *
 * This class describes a collection of local well elements and perforations
 */  
class WellElementSubRegion : public ElementSubRegionBase
{
public:

  using NodeMapType = InterObjectRelation< array2d< localIndex, cells::NODE_MAP_PERMUTATION > >;
  using EdgeMapType = FixedOneToManyRelation; // unused but needed in MeshLevel::GenerateAdjacencyLists
  using FaceMapType = FixedOneToManyRelation; // unused but needed in MeshLevel::GenerateAdjacencyLists

  /**
   * @brief enumeration for values in segmentStatusList parameter of Generate()
   */
  enum WellElemStatus : unsigned
  {
    UNOWNED = 0,             // there are no perforations on this element
    REMOTE = 1,              // all perforations are remote
    LOCAL  = 2,              // all perforations are local
    SHARED = REMOTE | LOCAL  // both remote and local perforations
  };

  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */  
  WellElementSubRegion( string const & name, 
                        Group * const parent );

  /**
   * @brief default destructor
   */
  virtual ~WellElementSubRegion() override;

  /**
   * @return the name of this type in the catalog
   */  
  static const string CatalogName() { return "wellElementSubRegion"; }

  /**
   *
   * @return the name of this type in the catalog
   */
  virtual const string getCatalogName() const override { return WellElementSubRegion::CatalogName(); }

  virtual void CalculateElementGeometricQuantities( NodeManager const & GEOSX_UNUSED_PARAM( nodeManager ),
                                                    FaceManager const & GEOSX_UNUSED_PARAM( faceManager ) ) override 
  {}

  virtual void setupRelatedObjectsInRelations( MeshLevel const * const mesh ) override;

  /**
   * @return the element to edge map
   */
  FixedOneToManyRelation & edgeList() 
  { 
    return m_toEdgesRelation; 
  }

  /**
   * @return the element to edge map
   */
  FixedOneToManyRelation const & edgeList() const 
  { 
    return m_toEdgesRelation; 
  }

   /**
   * @return the element to face map
   */
  FixedOneToManyRelation & faceList() 
  { 
    return m_toFacesRelation; 
  }

  /**
   * @return the element to face map
   */
  FixedOneToManyRelation const & faceList() const 
  { 
    return m_toFacesRelation; 
  }

  /**
   * @return the element to node map
   */
  NodeMapType & nodeList() 
  { 
    return m_toNodesRelation; 
  }

  /**
   * @return the element to node map
   */
  NodeMapType const & nodeList() const 
  { 
    return m_toNodesRelation; 
  }

  /**
   * @brief Getter for the top element index
   * @return local index of well's top element or -1 if it is not on current rank
   */
  localIndex GetTopWellElementIndex() const
  {
    return m_topWellElementIndex;
  }

  /**
   * @brief Setter for the name of the WellControls object of this well
   * @param name the name of the WellControls object
   */
  void SetWellControlsName( string const & name ) 
  { 
    m_wellControlsName = name; 
  }

  /**
   * @brief Getter for the name of the WellControls object of this well
   * @return a string containing the name of the WellControls object
   */  
  string const & GetWellControlsName() const 
  { 
    return m_wellControlsName; 
  }

  /**
   * @brief Getter for the perforations
   * @return a pointer to the PerforationData object
   */
  PerforationData * GetPerforationData() 
  { 
    return &m_perforationData; 
  }

  /**
   * @brief Getter for the perforation data
   * @return a pointer to the const PerforationData object
   */
  PerforationData const * GetPerforationData() const 
  { 
    return &m_perforationData; 
  } 

  /**
   * @brief Setter fpr the MPI rank that owns this well (i.e. the top segment)
   * @param MPI rank of the owner process
   */
  void SetTopRank( int rank ) 
  { 
    m_topRank = rank; 
  }

  /**
   * @brief Getter for the MPI rank that owns this well (i.e. the top segment)
   * @return MPI rank of the owner process
   */
  int GetTopRank() const 
  { 
    return m_topRank; 
  }

  /**
   * @brief Check if well is owned by current rank
   * @return true if the well is owned by current rank, false otherwise
   */
  bool IsLocallyOwned() const;

  /**
   * @brief Build the local well elements from global well element data
   * @param[in] meshLevel the mesh object (single level only)
   * @param[in] wellGeometry the InternalWellGenerator containing the global well topology
   * @param[in] wellElemStatus list of well element status, as determined by perforations connected 
   *                           to local or remote mesh partitions. Status values are defined in
   *                           enum SegmentStatus. They are used to partition well elements.  
   * @param[in] nodeOffsetGlobal the offset of the first global well node ( = offset of last global mesh node + 1 )
   * @param[in] elemOffsetGlobal the offset of the first global well element ( = offset of last global mesh elem + 1 )
   */
  void Generate( MeshLevel                        & mesh, 
                 InternalWellGenerator      const & wellGeometry,
                 arrayView1d<integer>             & elemStatus,
                 globalIndex                        nodeOffsetGlobal,
                 globalIndex                        elemOffsetGlobal );

  /**
   * @brief For each perforation, find the reservoir element that contains the perforation
   * @param[in] meshLevel the mesh object (single level only)
   * @param[in] wellGeometry the InternalWellGenerator containing the global well topology
   */
  void ConnectPerforationsToMeshElements( MeshLevel                   & mesh,
                                          InternalWellGenerator const & wellGeometry );
  
  /*
   * @brief Reconstruct the (local) map nextWellElemId using nextWellElemIdGlobal after the ghost exchange
   */
  void ReconstructLocalConnectivity();

  virtual void ViewPackingExclusionList( SortedArray<localIndex> & exclusionList ) const override;

  virtual localIndex PackUpDownMapsSize( arrayView1d<localIndex const> const & packList ) const override;

  virtual localIndex PackUpDownMaps( buffer_unit_type * & buffer,
                                     arrayView1d<localIndex const> const & packList ) const override;

  virtual localIndex UnpackUpDownMaps( buffer_unit_type const * & buffer,
                                       localIndex_array & packList,
                                       bool const overwriteUpMaps,
                                       bool const overwriteDownMaps ) override;

  virtual void FixUpDownMaps( bool const clearIfUnmapped ) final override; 

  void DebugWellElementSubRegionsAfterSetupCommunications() const;

  struct viewKeyStruct : public ElementSubRegionBase::viewKeyStruct
  {
    static constexpr auto wellControlsString               = "wellControlsName";
    static constexpr auto wellNodeListString               = "nodeList";
    static constexpr auto nextWellElementIndexString       = "nextWellElementIndex";
    static constexpr auto nextWellElementIndexGlobalString = "nextWellElementIndexGlobal";
    static constexpr auto topWellElementIndexString        = "topWellElementIndex";
    static constexpr auto topRankString                    = "topRank";
  
    dataRepository::ViewKey wellControlsName     = { wellControlsString };
    dataRepository::ViewKey wellNodeList         = { wellNodeListString };
    dataRepository::ViewKey nextWellElementIndex = { nextWellElementIndexString };
    dataRepository::ViewKey nextWellElementIndexGlobal = { nextWellElementIndexGlobalString };
    dataRepository::ViewKey topWellElementIndex  = { topWellElementIndexString };
    dataRepository::ViewKey topRank            = { topRankString };

  } viewKeysWellElementSubRegion;
  
  struct groupKeyStruct : public ElementSubRegionBase::groupKeyStruct
  {
    static constexpr auto perforationDataString = "wellElementSubRegion";

    dataRepository::GroupKey perforationData = { perforationDataString };

  } groupKeysWellElementSubRegion;


private:

  /**
   * @brief Assign the unowned well elements ( = well elem without perforation ) that are 
            in the reservoir (and that can therefore be matched with a reservoir element) to an MPI rank
   * @param[in] meshLevel the mesh object (single level only)
   * @param[in] wellGeometry the InternalWellGenerator containing the global well topology
   * @param[in] unownedElems set of unowned well elems. 
   * @param[inout] localElems set of local well elems. It contains the perforated well elements
                              connected to local mesh elements before the call, and is filled 
                              with the newly assigned well elements in this function.
   * @param[inout] wellElemStatus list of current well element status. Status values are defined in
   *                           enum SegmentStatus. They are used to partition well elements.  
   */
  void AssignUnownedElementsInReservoir( MeshLevel                        & mesh,
                                         InternalWellGenerator      const & wellGeometry,
                                         SortedArray<globalIndex>           const & unownedElems,
                                         SortedArray<globalIndex>                 & localElems,
                                         arrayView1d<integer>             & elemStatusGlobal ) const;

  /**
   * @brief Check that all the well elements have been assigned to a single rank
   *        Also check that if two ranks are neighbors in the well, they are also neighbors in the mesh
   * @param[in] wellGeometry the InternalWellGenerator containing the global well topology
   * @param[inout] localElems set of local well elems. 
   * @param[inout] wellElemStatus list of current well element status. Status values are defined in
   *                           enum SegmentStatus. They are used to partition well elements.  
   */
  void CheckPartitioningValidity( InternalWellGenerator const & wellGeometry,
                                  SortedArray<globalIndex>            & localElems,
                                  arrayView1d<integer>        & elemStatusGlobal ) const;
  /**
   * @brief Now that the well elements are assigned, collect the nodes and tag the boundary nodes between ranks
            The function WellElementSubRegion::AssignUnownedElements must have been called before this function
   * @param[in] wellGeometry the InternalWellGenerator containing the global well topology
   * @param[in] localElems set of local well elems. At this point all the well elems have been assigned
   * @param[out] localNodes set of local well nodes (includes boundary nodes)
   * @param[out] boundaryNodes set of local well nodes that are at the boundary between this rank 
                               and another rank
   */ 
  void CollectLocalAndBoundaryNodes( InternalWellGenerator const & wellGeometry, 
                                     SortedArray<globalIndex>      const & localElems,
                                     SortedArray<globalIndex>            & localNodes,
                                     SortedArray<globalIndex>            & boundaryNodes ) const;

  /**
   * @brief Add the well nodes to the nodeManager (properly resized)
            The function WellElementSubRegion::CollectLocalAndBoundaryNodes must have been called before this function
   * @param[in] meshLevel the mesh object (single level only)
   * @param[in] wellGeometry the InternalWellGenerator containing the global well topology
   * @param[in] localNodes set of local well nodes (includes boundary nodes). At this point all the nodes have been collected
   * @param[in] boundaryNodes set of local well nodes that are at the boundary between this rank 
                               and another rank
   * @param[in] nodeOffsetGlobal the offset of the first global well node ( = offset of last global mesh node + 1 )
   */ 
  void UpdateNodeManagerSize( MeshLevel                    & mesh, 
                              InternalWellGenerator  const & wellGeometry,
                              SortedArray<globalIndex>       const & localNodes,
                              SortedArray<globalIndex>       const & boundaryNodes,
                              globalIndex                    nodeOffsetGlobal );

  /**
   * @brief Construct the subregion's local to global maps, as well as other local maps (toNodes, nextWellElemId, volume, etc)
            The function WellElementSubRegion::UpdateNodeManagerSize must have been called before this function
   * @param[in] meshLevel the mesh object (single level only)
   * @param[in] wellGeometry the InternalWellGenerator containing the global well topology
   * @param[in] localElems set of local well elems. At this point all the well elems have been assigned
   * @param[in] localNodes set of local well nodes (includes boundary nodes)
   * @param[in] nodeOffsetGlobal the offset of the first global well node ( = offset of last global mesh node + 1 )
   * @param[in] elemOffsetGlobal the offset of the first global well element ( = offset of last global mesh elem + 1 )
   */ 
  void ConstructSubRegionLocalElementMaps( MeshLevel                      & mesh, 
                                           InternalWellGenerator    const & wellGeometry,
                                           SortedArray<globalIndex> const & localElems,
                                           globalIndex                      nodeOffsetGlobal,
                                           globalIndex                      elemOffsetGlobal );

  /**
   * @brief This function is the equivalent of NodeManager::SetElementMaps for well elements.
            It constructs the toElementRegionList, toElementSubRegion, toElement maps
            The function WellElementSubRegion::ConstructSubRegionLocalElementMaps must have been called before this function
   * @param[in] meshLevel the mesh object (single level only)
   */ 
  void UpdateNodeManagerNodeToElementMap( MeshLevel & mesh );

  /**
   * @brief Search for the reservoir element that is the *closest* from the center of well element.
            Note that this reservoir element does not necessarily contain the center of the well element.
            This "init" reservoir element will be used in SearchLocalElements to find the reservoir element that
            contains the well element.
   * @param[in] meshLevel the mesh object (single level only)
   * @param[in] location the location of that we are trying to match with a reservoir element  
   * @param[inout] erInit the region index of the reservoir element from which we start the search
   * @param[inout] esrInit the subregion index of the reservoir element from which we start the search
   * @param[inout] eiInit the element index of the reservoir element from which we start the search
   */ 
  void InitializeLocalSearch( MeshLevel const & mesh,
                              R1Tensor  const & location,
                              localIndex      & erInit,
                              localIndex      & esrInit,
                              localIndex      & eiInit) const;

  /**
   * @brief Search for the reservoir element that contains the well element.
            To do that, loop over the reservoir elements that are in the neighborhood of (erInit,esrInit,eiInit)
   * @param[in] meshLevel the mesh object (single level only)
   * @param[in] location the location of that we are trying to match with a reservoir element  
   * @param[in] erInit the region index of the reservoir element from which we start the search
   * @param[in] esrInit the subregion index of the reservoir element from which we start the search
   * @param[in] eiInit the element index of the reservoir element from which we start the search
   * @param[inout] erMatched the region index of the reservoir element that contains "location", if any
   * @param[inout] esrMatched the subregion index of the reservoir element that contains "location", if any
   * @param[inout] eiMatched the element index of the reservoir element that contains "location", if any
   */ 
  bool SearchLocalElements( MeshLevel  const & mesh,
                            R1Tensor   const & location,
                            localIndex const & erInit,
                            localIndex const & esrInit,
                            localIndex const & eiInit,
                            localIndex       & erMatched,
                            localIndex       & esrMatched,
                            localIndex       & eiMatched ) const;

  /**
   * @brief Search the reservoir elements that can be accessed from the set "nodes".
            Stop if a reservoir element containing the perforation is found.
            If not, enlarge the set "nodes" 
   * @param[in] meshLevel the mesh object (single level only)
   * @param[in] location the location of that we are trying to match with a reservoir element  
   * @param[inout] nodes the nodes that have already been visited
   * @param[inout] elements the reservoir elements that have already been visited
   * @param[inout] erMatched the region index of the reservoir element that contains "location", if any
   * @param[inout] esrMatched the subregion index of the reservoir element that contains "location", if any
   * @param[inout] eiMatched the element index of the reservoir element that contains "location", if any
   */ 
  bool VisitNeighborElements( MeshLevel const          & mesh,
                              R1Tensor  const          & location,
                              SortedArray<localIndex>  & nodes,
                              SortedArray<globalIndex> & elements,
                              localIndex               & erMatched,
                              localIndex               & esrMatched,
                              localIndex               & eiMatched ) const;

  /**
   * @brief Collect the nodes of reservoir element ei
   * @param[in] subRegion the subRegion of reservoir element ei
   * @param[in] ei the index of the reservoir element
   * @param[inout] nodes the nodes that have already been visited
   */ 
  void CollectElementNodes( CellBlock const *         subRegion,
                            localIndex                ei,
                            SortedArray<localIndex> & nodes ) const;

  /**
   * @brief Check if "location" is contained in reservoir element ei
   * @param[in] subRegion the subRegion of reservoir element ei
   * @param[in] ei the index of the reservoir element
   * @return true if "location" is contained in reservoir element ei, false otherwise 
   */ 
  bool IsPointInsideElement( NodeManager const * const nodeManager,
                             R1Tensor    const & location,
                             CellBlock   const * subRegion,
                             localIndex          ei ) const;
  
  void DebugNodeManager( MeshLevel const & mesh ) const;

  void DebugWellElementSubRegions( arrayView1d<integer const> const & wellElemStatus, globalIndex elemOffsetGlobal ) const;

  template< bool DOPACK >
  localIndex PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                                    arrayView1d<localIndex const> const & packList ) const;

  map<localIndex, array1d<globalIndex> > m_unmappedGlobalIndicesInNodelist;

  /// name of the WellControls object for this well
  string m_wellControlsName;

  /// elements to nodes relation is one to one relation.
  NodeMapType  m_toNodesRelation;

  /// elements to edges relation
  EdgeMapType  m_toEdgesRelation; // unused but needed in MeshLevel::GenerateAdjacencyLists

  /// elements to faces relation
  FaceMapType  m_toFacesRelation; // unused but needed in MeshLevel::GenerateAdjacencyLists

  /// local indices of the next well element (used in solvers)
  array1d<localIndex> m_nextWellElementIndex;

  /// indices of the next well element (to reconstruct connectivity after ghost exchange)
  array1d<localIndex> m_nextWellElementIndexGlobal; 

  /// local index of well's top segment
  localIndex m_topWellElementIndex;

  /// perforations
  PerforationData m_perforationData;

  /// top rank
  integer m_topRank;

  /// depth of the local search to match perforation to res elements
  localIndex m_searchDepth;
  
};

} /* namespace geosx */

#endif /* GEOSX_WELLS_WELLELEMENTSUBREGION_HPP_ */

