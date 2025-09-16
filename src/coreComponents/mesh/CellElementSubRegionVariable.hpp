/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


#ifndef GEOS_MESH_CELLELEMENTSUBREGIONVARIABLE_HPP_
#define GEOS_MESH_CELLELEMENTSUBREGIONVARIABLE_HPP_

#include "mesh/generators/CellBlockABC.hpp"
#include "mesh/NodeManager.hpp"
#include "mesh/FaceManager.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"
#include "ElementSubRegionBase.hpp"


namespace geos
{

class MeshLevel;

/**
 * @class CellElementSubRegionVariable
 * Class specializing the element subregion for a cell element.
 * This is the class used in the physics solvers to represent a collection of mesh cell elements
 */
class CellElementSubRegionVariable : public ElementSubRegionBase
{
public:

  /// Alias for the type of the element-to-node map
  using NodeMapType = InterObjectRelation< ArrayOfArrays< localIndex > >;
  /// Alias for the type of the element-to-edge map
  using EdgeMapType = InterObjectRelation< ArrayOfArrays< localIndex > >;
  /// Alias for the type of the element-to-face map
  using FaceMapType = InterObjectRelation< ArrayOfArrays< localIndex > >;
  /// Type of map between cell blocks and embedded elements
  using EmbSurfMapType = InterObjectRelation< ArrayOfArrays< localIndex > >;

  /**
   * @brief Const getter for the catalog name.
   * @return the name of this type in the catalog
   */
  static string catalogName()
  { return "CellElementSubRegionVariable"; }

  /**
   * @copydoc catalogName()
   */
  virtual string getCatalogName() const override final
  { return catalogName(); }

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor for this class.
   * @param[in] name the name of this object manager
   * @param[in] parent the parent Group
   */
  CellElementSubRegionVariable( string const & name, Group * const parent );

  ///@}

  /**
   * @name Helpers for CellElementSubRegionVariable construction
   */
  ///@{

  /**
   * @brief Fill the CellElementSubRegionVariable by copying those of the source CellBlock
   * @param cellBlock the CellBlock which properties (connectivity info) will be copied.
   */
  void copyFromCellBlock( CellBlockABC const & cellBlock );

  ///@}


  /**
   * @name Overriding packing / Unpacking functions
   */
  ///@{

  virtual localIndex packUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const override;

  virtual localIndex packUpDownMaps( buffer_unit_type * & buffer,
                                     arrayView1d< localIndex const > const & packList ) const override;

  virtual localIndex unpackUpDownMaps( buffer_unit_type const * & buffer,
                                       array1d< localIndex > & packList,
                                       bool const overwriteUpMaps,
                                       bool const overwriteDownMaps ) override;


  ///@}

  /**
   * @brief struct to serve as a container for variable strings and keys
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : public ElementSubRegionBase::viewKeyStruct
  {
  }
  /// viewKey struct for the CellElementSubRegionVariable class
  m_CellBlockSubRegionViewKeys;

  virtual viewKeyStruct & viewKeys() override { return m_CellBlockSubRegionViewKeys; }
  virtual viewKeyStruct const & viewKeys() const override { return m_CellBlockSubRegionViewKeys; }

  /**
   * @brief Get the local indices of the nodes in a face of the element.
   * @param[in] elementIndex The local index of the target element.
   * @param[in] localFaceIndex The local index of the target face in the element (this will be [0, numFacesInElement[)
   * @param[out] nodeIndices Memory to which node indices for the face will be written, must have sufficient size.
   * @return tne number of values written into @p nodeIndices
   * @deprecated This method will be removed soon.
   */
  localIndex getFaceNodes( localIndex const elementIndex,
                           localIndex const localFaceIndex,
                           Span< localIndex > const nodeIndices ) const;

  /**
   * @brief Get the element-to-node map.
   * @return a reference to the element-to-node map
   */
  NodeMapType & nodeList() { return m_toNodesRelation; }

  /**
   * @copydoc nodeList()
   */
  NodeMapType const & nodeList() const { return m_toNodesRelation; }

  /**
   * @brief Get the element-to-edge map.
   * @return a reference to element-to-edge map
   */
  EdgeMapType & edgeList() { return m_toEdgesRelation; }

  /**
   * @copydoc edgeList()
   */
  EdgeMapType const & edgeList() const { return m_toEdgesRelation; }

  /**
   * @brief Get the element-to-face map.
   * @return a reference to the element to face map
   */
  FaceMapType & faceList() { return m_toFacesRelation; }

  /**
   * @copydoc faceList()
   */
  FaceMapType const & faceList() const { return m_toFacesRelation; }



  /**
   * @brief Compute the center of each element in the subregion.
   * @param[in] X an arrayView of (const) node positions
   */
  void calculateElementCenters( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X ) const
  {
    ElementSubRegionBase::calculateElementCenters( m_toNodesRelation, X );
  }

  void calculateElementGeometricQuantities( NodeManager const & nodeManager,
                                            FaceManager const & faceManager ) override;

private:

  /// Element-to-node relation
  NodeMapType m_toNodesRelation;

  /// Element-to-edge relation
  EdgeMapType m_toEdgesRelation;

  /// Element-to-face relation
  FaceMapType m_toFacesRelation;

  /// Name of the properties registered from an external mesh
  string_array m_externalPropertyNames;

  /**
   * @brief Compute the volume of the k-th element in the subregion.
   * @param[in] k the index of the element in the subregion
   * @param[in] X an arrayView of (const) node positions
   */
  void calculateCellVolumesKernel( localIndex const k,
                                   arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X ) const;

  void calculateElementCenterAndVolume( localIndex const k,
                                        arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X ) const;

  /// The array of jacobian determinantes.
  array2d< real64 > m_detJ;

  /// Map of unmapped global indices in the element-to-node map
  map< localIndex, array1d< globalIndex > > m_unmappedGlobalIndicesInNodelist;

  /// Map of unmapped global indices in the element-to-face map
  map< localIndex, array1d< globalIndex > > m_unmappedGlobalIndicesInEdgelist;

  /// Map of unmapped global indices in the element-to-face map
  map< localIndex, array1d< globalIndex > > m_unmappedGlobalIndicesInFacelist;


  /**
   * @brief Pack element-to-node and element-to-face maps
   * @tparam DO_PACKING the flag for the bufferOps::Pack function
   * @param buffer the buffer used in the bufferOps::Pack function
   * @param packList the packList used in the bufferOps::Pack function
   * @return the pack size
   */
  template< bool DO_PACKING >
  localIndex packUpDownMapsImpl( buffer_unit_type * & buffer,
                                 arrayView1d< localIndex const > const & packList ) const;

  /**
   * @brief Links the managers to their mappings.
   * @param[in] mesh Holds the node/edge/face managers.
   *
   * Defines links the element to nodes, edges and faces to their respective node/edge/face managers.
   */
  void setupRelatedObjectsInRelations( MeshLevel const & mesh ) override;

  template< bool DO_PACKING >
  localIndex packFracturedElementsImpl( buffer_unit_type * & buffer,
                                        arrayView1d< localIndex const > const & packList,
                                        arrayView1d< globalIndex const > const & embeddedSurfacesLocalToGlobal ) const;
};

} /* namespace geos */

#endif /* GEOS_MESH_CELLELEMENTSUBREGIONVARIABLE_HPP_ */
