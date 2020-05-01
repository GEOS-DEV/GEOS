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
  using EdgeMapType = InterObjectRelation< ArrayOfSets< localIndex > >;
  using FaceMapType = InterObjectRelation< ArrayOfSets< localIndex > >;
  using ElemMapType = OrderedVariableToManyElementRelation;
  //END_SPHINX_INCLUDE_01

  /**
   * @brief return default size of the value array in the node-to-edge mapping
   * @return default size of value array in the node-to-edge mapping
   */
  inline localIndex getEdgeMapOverallocation()
  { return 8; }
  /**
   * @brief return default size of the value in the node-to-face mapping
   * @return default size of value array in the node-to-face mapping
   */
  inline localIndex getFaceMapOverallocation()
  { return 8; }
  /**
   * @brief return default size of the value array in the node-to-element mapping
   * @return default size of value array in the node-to-element mapping
   */
  inline localIndex getElemMapOverAllocation()
  { return 8; }


/**
   * @name Constructors/destructor
   */
  ///@{

  /**
   * @brief Main constructor for NodeManager Objects.
   * @param [in] name the name of this instantiation of NodeManager in the repository
   * @param [in] parent the parent group of this instantiation of NodeManager
   */
  NodeManager( std::string const & name,
               dataRepository::Group * const parent );

  
  /**
   * @brief A default NodeManager destructor.
   */
  ~NodeManager() override;

  /**
  * @brief  deleted constructor
  */
  NodeManager() = delete;

  /**
   * @brief deleted copy constructor
   */
  NodeManager( const NodeManager & init ) = delete;

  /**
  * @brief deleted assignement operator
  */
  NodeManager & operator=( const NodeManager & ) = delete;

  ///@}
  
  /**
   * @brief Resizes the NodeManager, and all its member vectors that relate nodes to faces, to edges, and to elements.
   * @param [in] newsize the new number of nodes.
   */
  virtual void resize( localIndex const newsize ) override;

   /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  /**
   * @brief Returns the name of the node manager in the object catalog.
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string CatalogName()
  { return "NodeManager"; }

  
  /**
   * @brief Provides a virtual access to CatalogName().
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  const string getCatalogName() const override final
  { return NodeManager::CatalogName(); }
 ///@}
  
  /**
   * @brief Assigns an EgdeManager to a NodeManager, and performs the node-to-edge mapping.
   * @param [in] edgeManager the edgeManager to assign this NodeManager
   */
  void SetEdgeMaps( EdgeManager const * const edgeManager );

  
  /**
   * @brief Assigns an FaceManager to a NodeManager, and performs the node-to-face mapping.
   * @param [in] faceManager the faceManager to assign this NodeManager
   */
  void SetFaceMaps( FaceManager const * const faceManager );

  
  /**
   * @brief Assigns an ElementRegionManager to a NodeManager, and performs the node-to-element mapping in this region.
   * @param [in] elementRegionManager the ElementRegionManager to assign this NodeManager
   */
  void SetElementMaps( ElementRegionManager const * const elementRegionManager );

  
  /**
   * @brief Compress all NodeManager member arrays so that the values of each array are contiguous with no extra capacity in between.
   * @note The method used here on each arrays (compress) does not free any memory.
   */
  void CompressRelationMaps( );

 /**
   * @name Packing methods
   */
  ///@{

  /**
   * @brief Creates an array listing all excluded local indices values.
   * @param [inout] exclusionList Sorted array with excluded local indices
   */
  virtual void ViewPackingExclusionList( SortedArray< localIndex > & exclusionList ) const override;

  
  /**
   * @brief Calculates the size that a list would have if it were packed, but without actually packing it.
   * @param [in] packList the list of node indices that we wish to get the size of after packing
   * @return a localIndex value representing the size of packList if it were packed
   * @note This function does not perform any packing, it just evaluates and returns the possible packed size.
   */
  virtual localIndex PackUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const override;

  
  /**
   * @brief Packs an array of node indices into a buffer.
   * @param [inout] buffer buffer to pack the node index data into
   * @param [in] packList the indices of nodes that should be packed
   * @return a localIndex value representing the size of the packed data
   */
  virtual localIndex PackUpDownMaps( buffer_unit_type * & buffer,
                                     arrayView1d< localIndex const > const & packList ) const override;

  
  /**
   * @brief Unpacks a buffer to an array of node indices.
   * @param [in] buffer buffer with the packed data
   * @param [inout] packList an array of localIndex values that we wish to unpack to
   * @param [in] overwriteUpMaps boolean: true to overwrite the previous Up maps
   * @param [in] overwriteDownMaps boolean: true to overwrite the previous Down maps
   * @return a localIndex value representing the size of the unpacked list
   */
  virtual localIndex UnpackUpDownMaps( buffer_unit_type const * & buffer,
                                       localIndex_array & packList,
                                       bool const overwriteUpMaps,
                                       bool const overwriteDownMaps ) override;

  
  /**
   * @brief Calls FixUpDownMaps for nodes-to-edges and nodes-to-faces maps.
   * @param [in] clearIfUnmapped boolean: true to remove if it is not mapped
   */
  void FixUpDownMaps( bool const clearIfUnmapped );
  ///@}
  
  /**
   * @brief Clean up the mappings between nodes and edges, faces, elements based on a new (updated) list of nodes, in order to keep only relevant mappings.
   * @param [in] receivedNodes the new list of target node indices
   * @param [in] edgesToNodes map to go from edges to nodes
   * @param [in] facesToNodes map to go from faces to nodes
   * @param [in] elemRegionManager Element Region Manager
   */
  void depopulateUpMaps( std::set< localIndex > const & receivedNodes,
                         array2d< localIndex > const & edgesToNodes,
                         ArrayOfArraysView< localIndex const > const & facesToNodes,
                         ElementRegionManager const & elemRegionManager );

   /**
   * @name viewKeyStruct/groupKeyStruct
   */
  ///@{
  /**
   *  @struct Containing added view access key to be bound with class data member
   */
  struct viewKeyStruct : ObjectManagerBase::viewKeyStruct
  {
    static constexpr auto referencePositionString       = "ReferencePosition";
    static constexpr auto totalDisplacementString       = "TotalDisplacement";
    static constexpr auto incrementalDisplacementString = "IncrementalDisplacement";
    static constexpr auto edgeListString                = "edgeList";
    static constexpr auto faceListString                = "faceList";
    static constexpr auto elementRegionListString       = "elemRegionList";
    static constexpr auto elementSubRegionListString    = "elemSubRegionList";
    static constexpr auto elementListString             = "elemList";

    dataRepository::ViewKey referencePosition       = { referencePositionString };
    dataRepository::ViewKey totalDisplacement       = { totalDisplacementString };
    dataRepository::ViewKey incrementalDisplacement = { incrementalDisplacementString };
    dataRepository::ViewKey edgeList                = { edgeListString };
    dataRepository::ViewKey faceList                = { faceListString };
    dataRepository::ViewKey elementRegionList       = { elementRegionListString };
    dataRepository::ViewKey elementSubRegionList    = { elementSubRegionListString };
    dataRepository::ViewKey elementList             = { elementListString };
    dataRepository::ViewKey velocity                = { dataRepository::keys::Velocity };
    dataRepository::ViewKey acceleration            = { dataRepository::keys::Acceleration };
  } viewKeys;

 /**
   *  @struct Containing added group access key to be bound with class in group hierarchy
   */
  struct groupKeyStruct : ObjectManagerBase::groupKeyStruct
  {} groupKeys;
  ///@}

  /**
   * \defgroup Accessors for NodeManager fixed data
   * @{
   */

  
  /**
   * @brief Provides a const accessor to the nodes-to-edges relation.
   * @return const reference to  nodes-to-edges relation
   */
  EdgeMapType const & edgeList() const { return m_toEdgesRelation; }

  
  /**
   * @brief Gets the node-to-edge relation.
   * @return reference to nodes-to-edges relation
   */
  EdgeMapType & edgeList() { return m_toEdgesRelation; }


  /**
   * @brief Provides a const accessor to the nodes-to-faces relation.
   * @return const reference to nodes-to-faces relation
   */
  FaceMapType const & faceList() const { return m_toFacesRelation; }

  
  /**
   * @brief Gets the nodes-to-faces relation.
   * @return reference to nodes-to-faces relation
   */
  FaceMapType & faceList() { return m_toFacesRelation; }


  /**
   * @brief Gets the nodes-to-elements relation.
   * @return reference to nodes-to-elements relation
   */
  OrderedVariableToManyElementRelation & toElementRelation() {return m_toElements;}
  
  
  /**
   * @brief Provides a const accessor to the nodes-to-elements relation.
   * @return const reference to nodes-to-elements relation
   */
  OrderedVariableToManyElementRelation const & toElementRelation() const {return m_toElements;}

  
  /**
   * @brief Gets the nodes-to-elements-regions relation.
   * @return reference to nodes-to-elements-regions relation
   */
  ArrayOfArrays< localIndex > & elementRegionList() { return m_toElements.m_toElementRegion; }
  
  
  /**
   * @brief Provides an immutable arrayView to the nodes-to-elements-regions relation.
   * @return const reference to nodes-to-elements-regions relation
   */
  ArrayOfArraysView< localIndex const > const & elementRegionList() const { return m_toElements.m_toElementRegion.toViewConst(); }

  
  /**
   * @brief Gets the nodes-to-elements-subregions relation.
   * @return reference to nodes-to-elements-subregions relation
   */
  ArrayOfArrays< localIndex > & elementSubRegionList() { return m_toElements.m_toElementSubRegion; }
  
  
  /**
   * @brief Providesan immutable arrayView to the nodes-to-elements-subregions relation.
   * @return const reference to nodes-to-elements-subregions relation
   */
  ArrayOfArraysView< localIndex const > const & elementSubRegionList() const { return m_toElements.m_toElementSubRegion.toViewConst(); }

  
  /**
   * @brief Gets the nodes-to-elements indices.
   * @return reference to nodes-to-elements indices
   */
  ArrayOfArrays< localIndex > & elementList() { return m_toElements.m_toElementIndex; }
  
  
  /**
   * @brief Provides an immutable arrayView to the nodes-to-elements indices.
   * @return const reference to nodes-to-elements indices
   */
  ArrayOfArraysView< localIndex const > const & elementList() const { return m_toElements.m_toElementIndex.toViewConst(); }

  
  /**
   * @brief Gets the reference position array.
   * @return reference position array
   */
  array2d< real64, nodes::REFERENCE_POSITION_PERM > & referencePosition() { return m_referencePosition; }

  
  /**
   * @brief Provides an immutable arrayView of the reference position.
   * @return an immutable arrayView of the reference position.
   */
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & referencePosition() const { return m_referencePosition; }

  
  /**
   * @brief Get the total displacement array.
   * @return the total displacement array if it exists, or an error is thrown if it does not exist
   * @note An error is thrown if the total displacement does not exist
   */
  array2d< real64, nodes::TOTAL_DISPLACEMENT_PERM > & totalDisplacement()
  {
    return getReference< array2d< real64, nodes::TOTAL_DISPLACEMENT_PERM > >( viewKeys.totalDisplacement );
  }

  
  /**
   * @brief Provides an immutable arrayView to the total displacement array.
   * @return immutable arrayView of the total displacement array if it exists, or an error is thrown if it does not exist
   * @note An error is thrown if the total displacement does not exist
   */
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & totalDisplacement() const
  {
    return getReference< array2d< real64, nodes::TOTAL_DISPLACEMENT_PERM > >( viewKeys.totalDisplacement );
  }

  
  /**
   * @brief Get the incremental displacement array.
   * @return the incremental displacement array if it exists, or an error is thrown if it does not exist
   * @note An error is thrown if the incremental displacement does not exist
   */
  array2d< real64, nodes::INCR_DISPLACEMENT_PERM > & incrementalDisplacement()
  {
    return getReference< array2d< real64, nodes::INCR_DISPLACEMENT_PERM > >( viewKeys.incrementalDisplacement );
  }

  
  /**
   * @brief Provides an immutable arrayView to the incremental displacement array.
   * @return immutable arrayView of the incremental displacement array if it exists, or an error is thrown if it does not exist
   * @note An error is thrown if the total incremental does not exist
   */
  arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const & incrementalDisplacement() const
  {
    return getReference< array2d< real64, nodes::INCR_DISPLACEMENT_PERM > >( viewKeys.incrementalDisplacement );
  }


  /**
   * @brief Get the velocity array.
   * @return the velocity array if it exists, or an error is thrown if it does not exist
   * @note An error is thrown if the velocity array does not exist
   */
  array2d< real64, nodes::VELOCITY_PERM > & velocity()
  {
    return getReference< array2d< real64, nodes::VELOCITY_PERM > >( viewKeys.velocity );
  }

  
  /**
   * @brief Provides an immutable arrayView to the velocity array.
   * @return immutable arrayView of the velocity array if it exists, or an error is thrown if it does not exist
   * @note An error is thrown if the velocity array does not exist
   */
  arrayView2d< real64 const, nodes::VELOCITY_USD > const & velocity() const
  {
    return getReference< array2d< real64, nodes::VELOCITY_PERM > >( viewKeys.velocity );
  }

  
  /**
   * @brief Get the acceleration array.
   * @return the acceleration array if it exists, or an error is thrown if it does not exist
   * @note An error is thrown if the acceleration array does not exist
   */
  array2d< real64, nodes::ACCELERATION_PERM > & acceleration()
  {
    return getReference< array2d< real64, nodes::ACCELERATION_PERM > >( viewKeys.acceleration );
  }

  
  /**
   * @brief Provides an immutable arrayView to the acceleration array.
   * @return immutable arrayView of the acceleration array if it exists, or an error is thrown if it does not exist
   * @note An error is thrown if the acceleration array does not exist
   */
  arrayView2d< real64 const, nodes::ACCELERATION_USD > const & acceleration() const
  {
    return getReference< array2d< real64, nodes::ACCELERATION_PERM > >( viewKeys.acceleration );
  }

   ///@}
  
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

  /// Reference position of the nodes
  array2d< real64, nodes::REFERENCE_POSITION_PERM > m_referencePosition;

  /// nodes-to-edges relation
  EdgeMapType m_toEdgesRelation;

  /// nodes-to-faces relation
  FaceMapType m_toFacesRelation;

  /// nodes-to-element relation
  ElemMapType m_toElements;

  /// map of global  to local  indices for edges
  map< localIndex, SortedArray< globalIndex > > m_unmappedGlobalIndicesInToEdges;

  /// map of global  to local  indices for faces
  map< localIndex, SortedArray< globalIndex > > m_unmappedGlobalIndicesInToFaces;

  /// map of global  to local  indices for elements
  map< localIndex, array1d< array1d< SortedArray< globalIndex > > > > m_unmappedGlobalIndicesInToElems;


};
}

#endif // MESH_NODEMANAGER_HPP_
