/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file ParticleManager.hpp
 */

#ifndef GEOSX_MESH_PARTICLEMANAGER_HPP_
#define GEOSX_MESH_PARTICLEMANAGER_HPP_

#include "mesh/ObjectManagerBase.hpp"
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
 * @class ParticleManager
 * @brief The ParticleManager class provides an interface to ObjectManagerBase in order to manage particle data.
 *
 * The ParticleManagerT class manages the particle data using the
 * ObjectDataStructureBaseT as a data manager.
 * This means that each field is stored in an array where each array entry
 * corresponds to a particle.
 */
class ParticleManager : public ObjectManagerBase
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
   * @brief Main constructor for ParticleManager Objects.
   * @param [in] name the name of this instantiation of ParticleManager
   * @param [in] parent the parent group of this instantiation of ParticleManager
   */
  ParticleManager( string const & name,
               dataRepository::Group * const parent );

  /**
   * @brief The default ParticleManager destructor.
   */
  ~ParticleManager() override;

  /// @cond DO_NOT_DOCUMENT
  /**
   * @brief deleted constructor
   */
  ParticleManager() = delete;

  /**
   * @brief deleted copy constructor
   */
  ParticleManager( const ParticleManager & init ) = delete;

  /**
   * @brief deleted assignement operator
   */
  ParticleManager & operator=( const ParticleManager & ) = delete;
  /// @endcond

  ///@}

  /**
   * @brief Resize the ParticleManager, and all its member vectors that relate nodes to faces, to edges, and to elements.
   * @details the size of the ParticleManager is the number of nodes
   * @param[in] newsize the new size of the ParticleManager
   */
  virtual void resize( localIndex const newsize ) override;

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  /**
   * @brief Return the name of the node manager in the object catalog.
   * @return string that contains the ParticleManager catalog name
   */
  static string catalogName()
  { return "ParticleManager"; }

  /**
   * @brief Provide a virtual access to catalogName().
   * @return string that contains the ParticleManager catalog name
   */
  const string getCatalogName() const override final
  { return ParticleManager::catalogName(); }

  ///@}

  /**
   * @brief Link the EdgeManager \p edgeManager to the ParticleManager, and performs the node-to-edge mapping.
   * @param [in] edgeManager the edgeManager to assign this ParticleManager
   */
  void setEdgeMaps( EdgeManager const & edgeManager );

  /**
   * @brief Link the FaceManager \p faceManager to the ParticleManager, and performs the node-to-face mapping.
   * @param [in] faceManager the faceManager to assign this ParticleManager
   */
  void setFaceMaps( FaceManager const & faceManager );

  /**
   * @brief Assign the ElementRegionManager \p elementRegionManager to the ParticleManager, and performs the node-to-element mapping
   * @param [in] elementRegionManager the ElementRegionManager to assign this ParticleManager
   */
  void setElementMaps( ElementRegionManager const & elementRegionManager );

  /**
   * @brief Compress all ParticleManager member arrays so that the values of each array are contiguous with no extra capacity inbetween.
   * @note The method used here on each arrays (compress) does not free any memory.
   */
  void compressRelationMaps();

  /**
   * @name Packing methods
   */
  ///@{

  /**
   * @brief Creates an array listing all excluded local indices values.
   * @param [in,out] exclusionList Sorted array with excluded local indices
   */
  virtual void viewPackingExclusionList( SortedArray< localIndex > & exclusionList ) const override;

  /**
   * @brief Calculate the size that a list would have if it were packed, but without actually packing it.
   * @details Packed data are meant to be communicated to other MPI ranks
   * @param [in] packList the list of node indices that we wish to get the size of after packing
   * @return a localIndex value representing the size of packList if it were packed
   * @note This function does not perform any packing, it just evaluates and returns the possible packed size.
   */
  virtual localIndex packUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const override;

  /**
   * @brief Packs an array of node indices into a buffer.
   * @details Packed data are meant to be communicated to other MPI ranks
   * @param [in,out] buffer buffer to pack the node index data into
   * @param [in] packList the indices of nodes that should be packed
   * @return a localIndex value representing the size of the packed data
   */
  virtual localIndex packUpDownMaps( buffer_unit_type * & buffer,
                                     arrayView1d< localIndex const > const & packList ) const override;

  /**
   * @brief Unpack a buffer to an array of node indices.
   * @details Packed data are meant to be communicated to other MPI ranks
   * @param [in] buffer buffer with the packed data
   * @param [inout] packList an array of localIndex values that we wish to unpack to
   * @param [in] overwriteUpMaps boolean: true to overwrite the previous Up maps
   * @param [in] overwriteDownMaps boolean: true to overwrite the previous Down maps
   * @return a localIndex value representing the size of the unpacked list
   */
  virtual localIndex unpackUpDownMaps( buffer_unit_type const * & buffer,
                                       localIndex_array & packList,
                                       bool const overwriteUpMaps,
                                       bool const overwriteDownMaps ) override;

  /**
   * @brief Call fixUpDownMaps for nodes-to-edges and nodes-to-faces maps.
   * @details Packed data are meant to be communicated to other MPI ranks
   * @param [in] clearIfUnmapped boolean: true to remove if it is not mapped
   */
  void fixUpDownMaps( bool const clearIfUnmapped );

  ///@}

  /**
   * @brief Clean up the mappings between nodes and edges, faces, elements based on a new (updated) list of nodes, in order to keep only
   * relevant mappings.
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
   *  @brief contains the added view access keys to be bound with class data member.
   *  @struct viewKeyStruct
   */
  struct viewKeyStruct : ObjectManagerBase::viewKeyStruct
  {
    /// @return String to access the reference position
    static constexpr char const * referencePositionString() { return "ReferencePosition"; }

    /// @return String to access the location of the nodes
    static constexpr char const * EmbSurfNodesPositionString() { return "EmbSurfNodesPosition"; }

    /// @return String to access the displacement
    static constexpr char const * totalDisplacementString() { return "TotalDisplacement"; }

    /// @return String to access the incremental displacement
    static constexpr char const * incrementalDisplacementString() { return "IncrementalDisplacement"; }

    /// @return String to access the edge map
    static constexpr char const * edgeListString() { return "edgeList"; }

    /// @return String to access the face map
    static constexpr char const * faceListString() { return "faceList"; }

    /// @return String to access the element region map
    static constexpr char const * elementRegionListString() { return "elemRegionList"; }

    /// @return String to access the element subregion map
    static constexpr char const * elementSubRegionListString() { return "elemSubRegionList"; }

    /// @return String to access the element map
    static constexpr char const * elementListString() { return "elemList"; }

    /// Accessor to reference position
    dataRepository::ViewKey referencePosition       = { referencePositionString() };

    /// Accessor to displacement
    dataRepository::ViewKey totalDisplacement       = { totalDisplacementString() };

    /// Accessor to incremental displacement
    dataRepository::ViewKey incrementalDisplacement = { incrementalDisplacementString() };

    /// Accessor to edge map
    dataRepository::ViewKey edgeList                = { edgeListString() };

    /// Accessor to face map
    dataRepository::ViewKey faceList                = { faceListString() };

    /// Accessor to element region map
    dataRepository::ViewKey elementRegionList       = { elementRegionListString() };

    /// Accessor to element subregion map
    dataRepository::ViewKey elementSubRegionList    = { elementSubRegionListString() };

    /// Accessor to element map
    dataRepository::ViewKey elementList             = { elementListString() };

    /// Accessor to velocity
    dataRepository::ViewKey velocity                = { dataRepository::keys::Velocity };

    /// Accessor to acceleration
    dataRepository::ViewKey acceleration            = { dataRepository::keys::Acceleration };
  }
  /// viewKeys
  viewKeys;

  ///@}

  /**
   * \defgroup Accessors for ParticleManager fixed data
   * @{
   */

  /**
   * @brief Provide an immutable accessor to the nodes-to-edges relation.
   * @return const reference to  nodes-to-edges relation
   */
  EdgeMapType const & edgeList() const { return m_toEdgesRelation; }

  /**
   * @brief Get a mutable accessor to the node-to-edges relation.
   * @return reference to nodes-to-edges relation
   */
  EdgeMapType & edgeList() { return m_toEdgesRelation; }

  /**
   * @brief Provide a const accessor to the nodes-to-faces relation.
   * @return const reference to nodes-to-faces relation
   */
  FaceMapType const & faceList() const { return m_toFacesRelation; }

  /**
   * @brief Get the nodes-to-faces relation.
   * @return reference to nodes-to-faces relation
   */
  FaceMapType & faceList() { return m_toFacesRelation; }

  /**
   * @brief Get the nodes-to-elements relation.
   * @return reference to nodes-to-elements relation
   */
  ElemMapType & toElementRelation() {return m_toElements;}

  /**
   * @brief Provide a const accessor to the nodes-to-elements relation.
   * @details The returned ElemMapType gives access, for one node
   * to the element index, the element sub region, and the element region
   * in relation with a node
   * @return const reference to nodes-to-elements relation
   */
  ElemMapType const & toElementRelation() const {return m_toElements;}

  /**
   * @brief Get the mutable nodes-to-elements-regions relation.
   * @return reference to nodes-to-elements-regions relation
   */
  ArrayOfArrays< localIndex > & elementRegionList() { return m_toElements.m_toElementRegion; }

  /**
   * @brief Provide an immutable arrayView to the nodes-to-elements-regions relation.
   * @return const reference to nodes-to-elements-regions relation
   */
  ArrayOfArraysView< localIndex const > elementRegionList() const { return m_toElements.m_toElementRegion.toViewConst(); }

  /**
   * @brief Get the mutable nodes-to-elements-subregions relation.
   * @return reference to nodes-to-elements-subregions relation
   */
  ArrayOfArrays< localIndex > & elementSubRegionList() { return m_toElements.m_toElementSubRegion; }

  /**
   * @brief Provide an immutable arrayView to the nodes-to-elements-subregions relation.
   * @return const reference to nodes-to-elements-subregions relation
   */
  ArrayOfArraysView< localIndex const > elementSubRegionList() const { return m_toElements.m_toElementSubRegion.toViewConst(); }

  /**
   * @brief Get the mutable nodes-to-elements indices.
   * @return reference to nodes-to-elements indices
   */
  ArrayOfArrays< localIndex > & elementList() { return m_toElements.m_toElementIndex; }

  /**
   * @brief Provide an immutable arrayView to the nodes-to-elements indices.
   * @return const reference to nodes-to-elements indices
   */

  ArrayOfArraysView< localIndex const > elementList() const
  { return m_toElements.m_toElementIndex.toViewConst(); }

  //START_SPHINX_REFPOS_ACCESS
  /**
   * @brief Get the mutable reference position array. This table will contain all the node coordinates.
   * @return reference position array
   */
  array2d< real64, nodes::REFERENCE_POSITION_PERM > & referencePosition() { return m_referencePosition; }

  /**
   * @brief Provide an immutable arrayView of the reference position. This table will contain all the node coordinates.
   * @return an immutable arrayView of the reference position.
   */

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > referencePosition() const
  { return m_referencePosition; }
  //END_SPHINX_REFPOS_ACCESS

  /**
   * @brief Get a mutable total displacement array.
   * @return the total displacement array if it exists, or an error is thrown if it does not exist
   * @note An error is thrown if the total displacement does not exist
   */
  array2d< real64, nodes::TOTAL_DISPLACEMENT_PERM > & totalDisplacement()
  { return getReference< array2d< real64, nodes::TOTAL_DISPLACEMENT_PERM > >( viewKeys.totalDisplacement ); }

  /**
   * @brief Provide an immutable arrayView to the total displacement array.
   * @return immutable arrayView of the total displacement array if it exists, or an error is thrown if it does not exist
   * @note An error is thrown if the total displacement does not exist
   */
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > totalDisplacement() const
  {return getReference< array2d< real64, nodes::TOTAL_DISPLACEMENT_PERM > >( viewKeys.totalDisplacement ); }

  /**
   * @brief Get a mutable incremental displacement array.
   * @return the incremental displacement array if it exists, or an error is thrown if it does not exist
   * @note An error is thrown if the incremental displacement does not exist
   */
  array2d< real64, nodes::INCR_DISPLACEMENT_PERM > & incrementalDisplacement()
  { return getReference< array2d< real64, nodes::INCR_DISPLACEMENT_PERM > >( viewKeys.incrementalDisplacement ); }

  /**
   * @brief Provide an immutable arrayView to the incremental displacement array.
   * @return immutable arrayView of the incremental displacement array if it exists, or an error is thrown if it does not exist
   * @note An error is thrown if the total incremental does not exist
   */
  arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > incrementalDisplacement() const
  { return getReference< array2d< real64, nodes::INCR_DISPLACEMENT_PERM > >( viewKeys.incrementalDisplacement ); }

  /**
   * @brief Get a mutable velocity array.
   * @return the velocity array if it exists, or an error is thrown if it does not exist
   * @note An error is thrown if the velocity array does not exist
   */
  array2d< real64, nodes::VELOCITY_PERM > & velocity()
  { return getReference< array2d< real64, nodes::VELOCITY_PERM > >( viewKeys.velocity ); }

  /**
   * @brief Provide an immutable arrayView to the velocity array.
   * @return immutable arrayView of the velocity array if it exists, or an error is thrown if it does not exist
   * @note An error is thrown if the velocity array does not exist
   */
  arrayView2d< real64 const, nodes::VELOCITY_USD > velocity() const
  { return getReference< array2d< real64, nodes::VELOCITY_PERM > >( viewKeys.velocity ); }

  /**
   * @brief Get a mutable acceleration array.
   * @return the acceleration array if it exists, or an error is thrown if it does not exist
   * @note An error is thrown if the acceleration array does not exist
   */
  array2d< real64, nodes::ACCELERATION_PERM > & acceleration()
  { return getReference< array2d< real64, nodes::ACCELERATION_PERM > >( viewKeys.acceleration ); }

  /**
   * @brief Provide an immutable arrayView to the acceleration array.
   * @return immutable arrayView of the acceleration array if it exists, or an error is thrown if it does not exist
   * @note An error is thrown if the acceleration array does not exist
   */
  arrayView2d< real64 const, nodes::ACCELERATION_USD > acceleration() const
  { return getReference< array2d< real64, nodes::ACCELERATION_PERM > >( viewKeys.acceleration ); }

  ///@}

private:

  /**
   * @brief Pack the upward and downward pointing maps into a buffer.
   * @tparam DOPACK template argument to determine whether or not to pack the buffer. If false, the buffer is not
   *                packed and the function returns the size of the packing that would have occured if set to TRUE.
   * @param buffer the buffer to pack data into
   * @param packList the indices of nodes that should be packed.
   * @return size of data packed in terms of number of chars
   */
  template< bool DOPACK >
  localIndex packUpDownMapsPrivate( buffer_unit_type * & buffer,
                                    arrayView1d< localIndex const > const & packList ) const;



  //START_SPHINX_REFPOS
  /// reference position of the nodes
  array2d< real64, nodes::REFERENCE_POSITION_PERM > m_referencePosition;
  //END_SPHINX_REFPOS

  /// nodes-to-edges relation
  EdgeMapType m_toEdgesRelation;

  /// nodes-to-faces relation
  FaceMapType m_toFacesRelation;

  /// nodes-to-element relation
  ElemMapType m_toElements;

  /// map of global to local indices for edges
  map< localIndex, SortedArray< globalIndex > > m_unmappedGlobalIndicesInToEdges;

  /// map of global to local indices for faces
  map< localIndex, SortedArray< globalIndex > > m_unmappedGlobalIndicesInToFaces;

  /// map of global to local indices for elements
  map< localIndex, array1d< array1d< SortedArray< globalIndex > > > > m_unmappedGlobalIndicesInToElems;

};
}

#endif // MESH_PARTICLEMANAGER_HPP_
