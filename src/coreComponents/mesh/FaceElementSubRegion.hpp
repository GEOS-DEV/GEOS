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
 * @file FaceElementSubRegion.hpp
 */

#ifndef GEOSX_MESH_FACEELEMENTSUBREGION_HPP_
#define GEOSX_MESH_FACEELEMENTSUBREGION_HPP_

#include "SurfaceElementSubRegion.hpp"

namespace geosx
{

/**
 * @class FaceElementSubRegion
 *
 * The FaceElementSubRegion class contains the functionality to support the concept of Elements that are comprised of
 * a face, or a pair of faces. This class which derives from ElementSubRegionBase, has specific connectivity maps and
 * and methods to support the specific geometry of an element comprised of a reduced dimensionality face element (i.e.
 * face area and aperture = volume)
 */
class FaceElementSubRegion : public SurfaceElementSubRegion
{
public:

  /// Face element to faces map type
  using FaceMapType = InterObjectRelation< array2d< localIndex > >;

  /**
   * @name Static factory catalog functions
   */
  ///@{

  /**
   * @brief Get catalog name.
   * @return the catalog name
   */
  static string catalogName()
  { return "FaceElementSubRegion"; }

  /**
   * @brief Get catalog name.
   * @return the catalog name
   */
  virtual string getCatalogName() const override
  {
    return catalogName();
  }

  ///@}

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor.
   * @param name the group name
   * @param parent the parent group
   */
  FaceElementSubRegion( string const & name,
                        dataRepository::Group * const parent );

  ///@}
  void copyFromCellBlock( CellBlockABC const & cellBlock, CellBlockManagerABC const & cellBlockManager );

  /**
   * @name Geometry computation / Connectivity
   */
  ///@{

  virtual void calculateElementGeometricQuantities( NodeManager const & nodeManager,
                                                    FaceManager const & faceManager ) override;
  /**
   * @brief Function to compute the geometric quantities of a specific face element.
   * @param k index of the face element
   * @param faceArea surface area of the face
   */
  void calculateSingleElementGeometricQuantities( localIndex const k,
                                                  arrayView1d< real64 const > const & faceArea );

  virtual localIndex packUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const override;

  virtual localIndex packUpDownMaps( buffer_unit_type * & buffer,
                                     arrayView1d< localIndex const > const & packList ) const override;

  virtual localIndex unpackUpDownMaps( buffer_unit_type const * & buffer,
                                       array1d< localIndex > & packList,
                                       bool const overwriteUpMaps,
                                       bool const overwriteDownMaps ) override;

  virtual void fixUpDownMaps( bool const clearIfUnmapped ) override;

  ///@}

  /**
   * @brief Function to set the ghostRank for a list of FaceElements and set them to the value of their bounding faces.
   * @param faceManager The face manager group
   * @param indices The list of indices to set value of ghostRank
   */
  void inheritGhostRankFromParentFace( FaceManager const & faceManager,
                                       std::set< localIndex > const & indices );

  /**
   * @brief Struct containing the keys to all face element views.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : SurfaceElementSubRegion::viewKeyStruct
  {
    /// @return String key for the derivatives of the shape functions with respect to the reference configuration
    static constexpr char const * dNdXString() { return "dNdX"; }
    /// @return String key for the derivative of the jacobian.
    static constexpr char const * detJString() { return "detJ"; }
    /// @return String key to the map of edge local indices to the fracture connector local indices.
    static constexpr char const * edgesTofractureConnectorsEdgesString() { return "edgesToFractureConnectors"; }
    /// @return String key to the map of fracture connector local indices to edge local indices.
    static constexpr char const * fractureConnectorEdgesToEdgesString() { return "fractureConnectorsToEdges"; }
    /// @return String key to the map of fracture connector local indices face element local indices.
    static constexpr char const * fractureConnectorsEdgesToFaceElementsIndexString() { return "fractureConnectorsToElementIndex"; }

#if GEOSX_USE_SEPARATION_COEFFICIENT
    /// Separation coefficient string.
    constexpr static char const * separationCoeffString() { return "separationCoeff"; }
    /// dSepCoeffdAper string.
    constexpr static char const * dSeparationCoeffdAperString() { return "dSeparationCoeffdAper"; }
#endif
  };

  virtual void setupRelatedObjectsInRelations( MeshLevel const & mesh ) override;


  /**
   * @name Relation getters
   * @brief Getter functions for the various inter-object relations
   */
  ///@{

  /**
   * @brief Get the face element to faces map.
   * @return the face element to edges map
   */
  FaceMapType const & faceList() const
  {
    return m_toFacesRelation;
  }

  /**
   * @copydoc faceList() const
   */
  FaceMapType & faceList()
  {
    return m_toFacesRelation;
  }
  ///@}


  /**
   * @name Properties Getters
   * @brief Getters to face element properties.
   */
  ///@{

  /**
   * @brief Get the number of nodes per face element
   * @return the number of nodes per face element
   */
  //virtual localIndex numNodesPerElement( localIndex const k ) const override { return m_toNodesRelation[k].size(); }

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
  /**
   * @brief Get separation coefficient.
   * @return the separation coefficient
   */
  arrayView1d< real64 > getSeparationCoefficient() { return m_separationCoefficient; }
  /**
   * @copydoc getSeparationCoefficient()
   */
  arrayView1d< real64 const > getSeparationCoefficient() const { return m_separationCoefficient; }
#endif

  ///@}

  /// Unmapped face elements to edges map
  map< localIndex, array1d< globalIndex > > m_unmappedGlobalIndicesInToEdges;

  /// Unmapped face elements to faces map
  map< localIndex, array1d< globalIndex > > m_unmappedGlobalIndicesInToFaces;

  /// List of the new face elements that have been generated
  SortedArray< localIndex > m_newFaceElements;

  /// map from the edges to the fracture connectors index (edges that are fracture connectors)
  SortedArray< localIndex > m_recalculateFractureConnectorEdges;

  /// A map of edge local indices to the fracture connector local indices.
  map< localIndex, localIndex > m_edgesToFractureConnectorsEdges;

  /// A map of fracture connector local indices to edge local indices.
  array1d< localIndex > m_fractureConnectorsEdgesToEdges;

  /// A map of fracture connector local indices face element local indices.
  ArrayOfArrays< localIndex > m_fractureConnectorEdgesToFaceElements;

  /**
   * @brief @return The array of shape function derivatives.
   */
  array4d< real64 > & dNdX()
  { return m_dNdX; }

  /**
   * @brief @return The array of shape function derivatives.
   */
  arrayView4d< real64 const > dNdX() const
  { return m_dNdX.toViewConst(); }

  /**
   * @brief @return The array of jacobian determinantes.
   */
  array2d< real64 > & detJ()
  { return m_detJ; }

  /**
   * @brief @return The array of jacobian determinantes.
   */
  arrayView2d< real64 const > detJ() const
  { return m_detJ; }

private:

  /**
   * @brief Pack element-to-node and element-to-face maps
   * @tparam the flag for the bufferOps::Pack function
   * @param buffer the buffer used in the bufferOps::Pack function
   * @param packList the packList used in the bufferOps::Pack function
   * @return the pack size
   */
  template< bool DO_PACKING >
  localIndex packUpDownMapsImpl( buffer_unit_type * & buffer,
                                 arrayView1d< localIndex const > const & packList ) const;

  /// The array of shape function derivaties.
  array4d< real64 > m_dNdX;

  /// The array of jacobian determinantes.
  array2d< real64 > m_detJ;

  /// Element-to-face relation
  FaceMapType m_toFacesRelation;

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
  /// Separation coefficient
  array1d< real64 > m_separationCoefficient;
#endif

};

} /* namespace geosx */

#endif /* GEOSX_MESH_FACEELEMENTSUBREGION_HPP_ */
