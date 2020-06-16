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
 * @file NodeManager.hpp
 */


#ifndef GEOSX_MESH_NODEMANAGER_HPP_
#define GEOSX_MESH_NODEMANAGER_HPP_

#include "managers/ObjectManagerBase.hpp"
#include <string.h>
#include "CellBlockManager.hpp"
#include "ToElementRelation.hpp"


// *********************************************************************************************************************
// *********************************************************************************************************************
class SiloFile;

namespace geosx
{

class CellBlock;
class FaceManager;
class EdgeManager;
class ElementRegionManager;


/**
 * @class NodeManager
 * @brief The NodeManager class provides an interface to ObjectManagerBase in order to manage node data.
 *
 * The NodeManagerT class manages the node data using the
 * ObjectDataStructureBaseT as a data manager.
 * This means that each field is stored in an array where each array entry
 * corresponds to a node.
 */
class NodeManager : public ObjectManagerBase
{
public:

  //START_SPHINX_INCLUDE_01

  /// nodeToEdge map type
  using EdgeMapType = InterObjectRelation< ArrayOfSets< localIndex > >;

  /// nodeToFace map type
  using FaceMapType = InterObjectRelation< ArrayOfSets< localIndex > >;

  /// nodeToElement map type
  using ElemMapType = OrderedVariableToManyElementRelation;
  //END_SPHINX_INCLUDE_01

  /**
   * @brief Get maximum number of edges per node.
   * @return maximum number of edges per node
   */
  inline localIndex getEdgeMapOverallocation()
  { return 8; }

  /**
   * @brief Get maximum number of faces per node.
   * @return maximum number of faces per node
   */
  inline localIndex getFaceMapOverallocation()
  { return 8; }

  /**
   * @brief Get maximum number of elements per node.
   * @return maximum number of elements per node
   */
  inline localIndex getElemMapOverAllocation()
  { return 8; }

  /**
   * @brief main constructor for NodeManager Objects
   * @param name the name of this instantiation of NodeManager in the repository
   * @param parent the parent group of this instantiation of NodeManager
   */
  NodeManager( std::string const & name,
               dataRepository::Group * const parent );

  /**
   *  @brief default destructor
   */
  ~NodeManager() override;

  virtual void resize( localIndex const newsize ) override;

  /**
   * @brief name of the node manager in the object catalog.
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog
   */
  static string CatalogName()
  { return "NodeManager"; }

  /**
   * @brief virtual access to CatalogName().
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog
   */
  const string getCatalogName() const override final
  { return NodeManager::CatalogName(); }

  /**
   * @brief Set the nodeToEdge map.
   * @param edgeManager a constant pointer to the EdgeManager
   */
  void SetEdgeMaps( EdgeManager const * const edgeManager );

  /**
   * @brief Set the nodeToEdge map.
   * @param faceManager a constant pointer to the FaceeManager
   */
  void SetFaceMaps( FaceManager const * const faceManager );


  /**
   * @brief Set the nodeToEdge map.
   * @param elementRegionManager a constant pointer to the ElementRegionManager
   */
  void SetElementMaps( ElementRegionManager const * const elementRegionManager );

  /**
   * @brief Compress maps.
   */
  void CompressRelationMaps( );

//  void Initialize();

  virtual void ViewPackingExclusionList( SortedArray< localIndex > & exclusionList ) const override;

  virtual localIndex PackUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const override;

  virtual localIndex PackUpDownMaps( buffer_unit_type * & buffer,
                                     arrayView1d< localIndex const > const & packList ) const override;

  virtual localIndex UnpackUpDownMaps( buffer_unit_type const * & buffer,
                                       localIndex_array & packList,
                                       bool const overwriteUpMaps,
                                       bool const overwriteDownMaps ) override;


  /**
   * @brief Fix UpDownMaps.
   * @param clearIfUnmapped  a boolean variable....
   */
  void FixUpDownMaps( bool const clearIfUnmapped );

  /**
   * @brief Depopulate UpDownMaps.
   * @param receivedNodes set of received nodes
   * @param edgesToNodes edgeToNode map
   * @param facesToNodes faeToNode map
   * @param elemRegionManager a reference to the ElementRegionManager
   */
  void depopulateUpMaps( std::set< localIndex > const & receivedNodes,
                         array2d< localIndex > const & edgesToNodes,
                         ArrayOfArraysView< localIndex const > const & facesToNodes,
                         ElementRegionManager const & elemRegionManager );


  /**
   * @brief Struct containing the keys to all embedded surface element views.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : ObjectManagerBase::viewKeyStruct
  {
    /// String to access the reference position
    static constexpr auto referencePositionString       = "ReferencePosition";

    /// String to access the displacement
    static constexpr auto totalDisplacementString       = "TotalDisplacement";

    /// String to access the incremental displacement
    static constexpr auto incrementalDisplacementString = "IncrementalDisplacement";

    /// String to access the edge map
    static constexpr auto edgeListString                = "edgeList";

    /// String to access the face map
    static constexpr auto faceListString                = "faceList";

    /// String to access the element region map
    static constexpr auto elementRegionListString       = "elemRegionList";

    /// String to access the element subregion map
    static constexpr auto elementSubRegionListString    = "elemSubRegionList";

    /// String to access the element map
    static constexpr auto elementListString             = "elemList";

    /// String to access the reference position
    dataRepository::ViewKey referencePosition       = { referencePositionString };

    /// String to access the displacement
    dataRepository::ViewKey totalDisplacement       = { totalDisplacementString };

    /// String to access the incremental displacement
    dataRepository::ViewKey incrementalDisplacement = { incrementalDisplacementString };

    /// String to access the edge map
    dataRepository::ViewKey edgeList                = { edgeListString };

    /// String to access the face map
    dataRepository::ViewKey faceList                = { faceListString };

    /// String to access the element region map
    dataRepository::ViewKey elementRegionList       = { elementRegionListString };

    /// String to access the element subregion map
    dataRepository::ViewKey elementSubRegionList    = { elementSubRegionListString };

    /// String to access the element map
    dataRepository::ViewKey elementList             = { elementListString };

    /// String to access the velocity
    dataRepository::ViewKey velocity                = { dataRepository::keys::Velocity };

    /// String to access the acceleration
    dataRepository::ViewKey acceleration            = { dataRepository::keys::Acceleration };
  } viewKeys;


  /**
   * @brief Contains the groupkeys.
   * @struct groupKeyStruct
   */
  struct groupKeyStruct : ObjectManagerBase::groupKeyStruct
  {} groupKeys;


  /**
   * \defgroup getters for NodeManager fixed data
   * @{
   */

  /**
   * @brief Get the nodeToEdge map.
   * @return const reference to nodeToEdge map
   */
  EdgeMapType const & edgeList() const
  { return m_toEdgesRelation; }

  /**
   * @brief Get the nodeToEdge map.
   * @return reference to the nodeToEdge map
   */
  EdgeMapType & edgeList()
  { return m_toEdgesRelation; }

  /**
   * @brief Get the nodeToFace map.
   * @return reference to the nodeToFace map
   */
  FaceMapType & faceList()       { return m_toFacesRelation; }

  /**
   * @brief Get the nodeToFace map.
   * @return const reference to the nodeToFace map
   */
  FaceMapType const & faceList() const { return m_toFacesRelation; }

  /**
   * @brief Get the nodeToElement relation.
   * @return reference to the nodeToElement relation
   */
  OrderedVariableToManyElementRelation & toElementRelation() {return m_toElements;}
  /**
   * @brief Get the nodeToElement map.
   * @return const reference to the nodeToElement map
   */
  OrderedVariableToManyElementRelation const & toElementRelation() const {return m_toElements;}

  /**
   * @brief Get the nodeToRegion map.
   * @return reference to the nodeToRegion map
   */
  ArrayOfArrays< localIndex > & elementRegionList()       { return m_toElements.m_toElementRegion; }
  /**
   * @brief Get the nodeToRegion map.
   * @return const reference to the nodeToRegion map
   */
  ArrayOfArraysView< localIndex const > const & elementRegionList() const
  { return m_toElements.m_toElementRegion.toViewConst(); }

  /**
   * @brief Get the nodeToSubregion map.
   * @return reference to the nodeToSubregion map
   */
  ArrayOfArrays< localIndex > & elementSubRegionList()       { return m_toElements.m_toElementSubRegion; }

  /**
   * @brief Get the nodeToSubregion map.
   * @return const reference to the nodeToSubregion map
   */
  ArrayOfArraysView< localIndex const > const & elementSubRegionList() const
  { return m_toElements.m_toElementSubRegion.toViewConst(); }

  /**
   * @brief Get the nodeToElement map.
   * @return reference to the nodeToElement map
   */
  ArrayOfArrays< localIndex > & elementList()       { return m_toElements.m_toElementIndex; }

  /**
   * @brief Get the nodeToElement map.
   * @return const reference to the nodeToElement map
   */
  ArrayOfArraysView< localIndex const > const & elementList() const
  { return m_toElements.m_toElementIndex.toViewConst(); }

  //START_SPHINX_REFPOS_ACCESS
  /**
   * @brief Get the reference position array.
   * @return a reference to the reference position of the nodes
   */
  array2d< real64, nodes::REFERENCE_POSITION_PERM > & referencePosition()
  { return m_referencePosition; }

  /**
   * @brief Get an immutable arrayView of the reference position.
   * @return a constant reference to the reference position of the nodes
   */
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & referencePosition() const
  { return m_referencePosition; }
  //END_SPHINX_REFPOS_ACCESS

  /**
   * @brief Get the total displacement array if it exists, if not an error is thrown.
   * @return a reference to the total displacement
   */
  array2d< real64, nodes::TOTAL_DISPLACEMENT_PERM > & totalDisplacement()
  { return getReference< array2d< real64, nodes::TOTAL_DISPLACEMENT_PERM > >( viewKeys.totalDisplacement ); }

  /**
   * @brief Get an immutable arrayView of the total displacement if it exists, if not an error is thrown.
   * @return a constant reference to the total displacement
   */
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & totalDisplacement() const
  { return getReference< array2d< real64, nodes::TOTAL_DISPLACEMENT_PERM > >( viewKeys.totalDisplacement ); }

  /**
   * @brief Get the incremental displacement array if it exists, if not an error is thrown.
   * @return a reference to the incremental displacement array
   */
  array2d< real64, nodes::INCR_DISPLACEMENT_PERM > & incrementalDisplacement()
  { return getReference< array2d< real64, nodes::INCR_DISPLACEMENT_PERM > >( viewKeys.incrementalDisplacement ); }

  /**
   * @brief Get an immutable arrayView of the incremental displacement if it exists, if not an error is thrown.
   * @return a constant reference to the incremental displacement array
   */
  arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const & incrementalDisplacement() const
  { return getReference< array2d< real64, nodes::INCR_DISPLACEMENT_PERM > >( viewKeys.incrementalDisplacement ); }

  /**
   * @brief Get the velocity array if it exists, if not an error is thrown.
   * @return a  reference to the velocity vector
   */
  array2d< real64, nodes::VELOCITY_PERM > & velocity()
  { return getReference< array2d< real64, nodes::VELOCITY_PERM > >( viewKeys.velocity ); }

  /**
   * @brief Get an immutable arrayView of the velocity if it exists, if not an error is thrown.
   * @return a constant reference to the velocity vector
   */
  arrayView2d< real64 const, nodes::VELOCITY_USD > const & velocity() const
  { return getReference< array2d< real64, nodes::VELOCITY_PERM > >( viewKeys.velocity ); }

  /**
   * @brief Get the accleration array if it exists, if not an error is thrown.
   * @return a reference to the acceleration vector
   */
  array2d< real64, nodes::ACCELERATION_PERM > & acceleration()
  { return getReference< array2d< real64, nodes::ACCELERATION_PERM > >( viewKeys.acceleration ); }

  /**
   * @brief Get an immutable arrayView of the acceleration if it exists, if not an error is thrown.
   * @return a constant reference to the acceleration vector
   */
  arrayView2d< real64 const, nodes::ACCELERATION_USD > const & acceleration() const
  { return getReference< array2d< real64, nodes::ACCELERATION_PERM > >( viewKeys.acceleration ); }

private:
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


  //START_SPHINX_REFPOS
  /// reference position of the nodes
  array2d< real64, nodes::REFERENCE_POSITION_PERM > m_referencePosition;
  //END_SPHINX_REFPOS

  /// nodeToEdge relation
  EdgeMapType m_toEdgesRelation;

  /// nodeToFace relation
  FaceMapType m_toFacesRelation;

  /// nodeToElement relation
  ElemMapType m_toElements;

  /// Unmapped global indices to edge map
  map< localIndex, SortedArray< globalIndex > > m_unmappedGlobalIndicesInToEdges;

  /// Unmapped global indices to face map
  map< localIndex, SortedArray< globalIndex > > m_unmappedGlobalIndicesInToFaces;

  /// Unmapped global indices to element map
  map< localIndex, array1d< array1d< SortedArray< globalIndex > > > > m_unmappedGlobalIndicesInToElems;



  /// deleted constructor
  NodeManager() = delete;

  /// deleted copy constructor
  NodeManager( const NodeManager & init ) = delete;

  /// deleted assignement operator
  NodeManager & operator=( const NodeManager & ) = delete;

};
}


#endif // MESH_NODEMANAGER_HPP_
