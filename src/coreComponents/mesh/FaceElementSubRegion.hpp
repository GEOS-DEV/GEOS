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
 * @file FaceElementSubRegion.hpp
 */

#ifndef GEOSX_MESH_FACEELEMENTSUBREGION_HPP_
#define GEOSX_MESH_FACEELEMENTSUBREGION_HPP_

#include "ElementSubRegionBase.hpp"
#include "InterObjectRelation.hpp"
#include "ToElementRelation.hpp"

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
class FaceElementSubRegion : public ElementSubRegionBase
{
public:

  /// Face element to nodes map type
  using NodeMapType = InterObjectRelation< array1d< array1d< localIndex > > >;

  /// Face element to edges map type
  using EdgeMapType = InterObjectRelation< array1d< array1d< localIndex > > >;

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
  static const string CatalogName()
  { return "FaceElementSubRegion"; }

  /**
   * @brief Get catalog name.
   * @return the catalog name
   */
  virtual const string getCatalogName() const override
  {
    return FaceElementSubRegion::CatalogName();
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

  virtual ~FaceElementSubRegion() override;

  ///@}

  /**
   * @name Geometry computation / Connectivity
   */
  ///@{

  virtual void CalculateElementGeometricQuantities( NodeManager const & nodeManager,
                                                    FaceManager const & faceManager ) override;
  /**
   * @brief Function to compute the geometric quantities of a specific face element.
   * @param index index of the face element
   * @param faceArea surface area of the face
   */
  void CalculateElementGeometricQuantities( localIndex const index,
                                            arrayView1d< real64 const > const & faceArea );

  ///@}

  /**
   * @name Overriding packing / Unpacking functions
   */
  ///@{

  virtual localIndex PackUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const override;

  virtual localIndex PackUpDownMaps( buffer_unit_type * & buffer,
                                     arrayView1d< localIndex const > const & packList ) const override;

  virtual localIndex UnpackUpDownMaps( buffer_unit_type const * & buffer,
                                       localIndex_array & packList,
                                       bool const overwriteUpMaps,
                                       bool const overwriteDownMaps ) override;

  virtual void FixUpDownMaps( bool const clearIfUnmapped ) override;

  virtual void ViewPackingExclusionList( SortedArray< localIndex > & exclusionList ) const override;

  ///@}

  /**
   * @brief Function to set the ghostRank for a list of FaceElements and set them to the value of their bounding faces.
   * @param faceManager The face manager group
   * @param indices The list of indices to set value of ghostRank
   */
  void inheritGhostRankFromParentFace( FaceManager const * const faceManager,
                                       std::set< localIndex > const & indices );

  /**
   * @brief Struct containing the keys to all face element views.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : ElementSubRegionBase::viewKeyStruct
  {
    /// Face element aperture string.
    static constexpr auto elementApertureString        = "elementAperture";

    /// Face element area string.
    static constexpr auto elementAreaString            = "elementArea";

    /// Face element to cell regions map string.
    static constexpr auto faceElementsToCellRegionsString    = "fractureElementsToCellRegions";

    /// Face element to cell subregions map string.
    static constexpr auto faceElementsToCellSubRegionsString    = "fractureElementsToCellSubRegions";

    /// Face element to cell indices map string.
    static constexpr auto faceElementsToCellIndexString    = "fractureElementsToCellIndices";

    /// Mass creation string.
    constexpr static auto creationMassString = "creationMass";

#if GEOSX_USE_SEPARATION_COEFFICIENT

    /// Separation coefficient string.
    constexpr static auto separationCoeffString = "separationCoeff";

    /// dSepCoeffdAper string.
    constexpr static auto dSeparationCoeffdAperString = "dSeparationCoeffdAper";
#endif

  };

  virtual void setupRelatedObjectsInRelations( MeshLevel const * const mesh ) override;


  /**
   * @name Relation getters
   * @brief Getter functions for the various inter-object relations
   */
  ///@{

  /**
   * @brief Get the face element to nodes map.
   * @return the face element to node map
   */
  NodeMapType const & nodeList() const
  {
    return m_toNodesRelation;
  }

  /**
   * @copydoc nodeList() const
   */
  NodeMapType & nodeList()
  {
    return m_toNodesRelation;
  }

  /**
   * @brief Get the face element to edges map.
   * @return The face element to edge map
   */
  EdgeMapType const & edgeList() const
  {
    return m_toEdgesRelation;
  }

  /**
   * @copydoc edgeList() const
   */
  EdgeMapType & edgeList()
  {
    return m_toEdgesRelation;
  }

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

  /**
   * @brief Get face element aperture.
   * @return the aperture of the face elements
   */
  arrayView1d< real64 > const & getElementAperture()       { return m_elementAperture; }

  /**
   * @copydoc getElementAperture()
   */
  arrayView1d< real64 const > const & getElementAperture() const { return m_elementAperture; }

  /**
   * @brief Get face element surface area.
   * @return the surface area of the face element
   */
  arrayView1d< real64 > const & getElementArea()       { return m_elementArea; }

  /**
   * @copydoc getElementArea()
   */
  arrayView1d< real64 const > const & getElementArea() const { return m_elementArea; }

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
  /**
   * @brief Get separation coefficient.
   * @return the separation coefficient
   */
  arrayView1d< real64 > const & getSeparationCoefficient()       { return m_separationCoefficient; }
  /**
   * @copydoc getSeparationCoefficient()
   */
  arrayView1d< real64 const > const & getSeparationCoefficient() const { return m_separationCoefficient; }
#endif

  ///@}

  /// Unmapped face elements to nodes map
  map< localIndex, array1d< globalIndex > > m_unmappedGlobalIndicesInToNodes;

  /// Unmapped face elements to edges map
  map< localIndex, array1d< globalIndex > > m_unmappedGlobalIndicesInToEdges;

  /// Unmapped face elements to faces map
  map< localIndex, array1d< globalIndex > > m_unmappedGlobalIndicesInToFaces;

  /// Map between the face elements and the cells
  FixedToManyElementRelation m_faceElementsToCells;

  /// List of the new face elements that have been generated
  SortedArray< localIndex > m_newFaceElements;

private:

  /**
   * @brief Pack element-to-node and element-to-face maps
   * @tparam the flag for the bufferOps::Pack function
   * @param buffer the buffer used in the bufferOps::Pack function
   * @param packList the packList used in the bufferOps::Pack function
   * @return the pack size
   */
  template< bool DOPACK >
  localIndex PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                                    arrayView1d< localIndex const > const & packList ) const;

  /// Element-to-node relation
  NodeMapType m_toNodesRelation;

  /// Element-to-edge relation
  EdgeMapType m_toEdgesRelation;

  /// Element-to-face relation
  FaceMapType m_toFacesRelation;

  /// Member level field for the element center
  array1d< real64 > m_elementAperture;

  /// Member level field for the element center
  array1d< real64 > m_elementArea;

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
  /// Separation coefficient
  array1d< real64 > m_separationCoefficient;
#endif

};

} /* namespace geosx */

#endif /* GEOSX_MESH_FACEELEMENTSUBREGION_HPP_ */
