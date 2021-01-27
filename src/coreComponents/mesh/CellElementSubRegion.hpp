/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


#ifndef GEOSX_MESH_CELLELEMENTSUBREGION_HPP_
#define GEOSX_MESH_CELLELEMENTSUBREGION_HPP_

#include "CellBlock.hpp"

namespace geosx
{

/**
 * @class CellElementSubRegion
 * Class deriving from CellBlock further specializing the element subregion
 * for a cell element. This is the class used in the physics solvers to
 * represent a collection of mesh cell elements
 */
class CellElementSubRegion : public CellBlock
{
public:

  /// Type of map between cell blocks and embedded elements
  using EmbSurfMapType = InterObjectRelation< ArrayOfArrays< localIndex > >;

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor for this class.
   * @param[in] name the name of this object manager
   * @param[in] parent the parent Group
   */
  CellElementSubRegion( std::string const & name, Group * const parent );

  /**
   * @brief Destructor.
   */
  virtual ~CellElementSubRegion() override;

  ///@}

  /**
   * @name Helpers for CellElementSubRegion construction
   */
  ///@{

  /**
   * @brief Fill the CellElementSubRegion by copying those of the source CellBlock
   * @param source the CellBlock whose properties (connectivity info) will be copied
   */
  void copyFromCellBlock( CellBlock * source );

  /**
   * @brief Fill the CellElementSubRegion by querying a target set into the faceManager
   * @param[in] faceManager a pointer to the faceManager
   * @param[in] setName a reference to string containing the name of the set
   */
  void constructSubRegionFromFaceSet( FaceManager const * const faceManager,
                                      std::string const & setName );

  ///@}

  /**
   * @brief Add fractured element to list and relative entries to the map.
   * @param cellElemIndex cell element index
   * @param embSurfIndex embedded surface element index
   */
  void addFracturedElement( localIndex const cellElemIndex,
                            localIndex const embSurfIndex );

  /**
   * @name Overriding packing / Unpacking functions
   */
  ///@{

  virtual void viewPackingExclusionList( SortedArray< localIndex > & exclusionList ) const override;

  virtual localIndex packUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const override;

  virtual localIndex packUpDownMaps( buffer_unit_type * & buffer,
                                     arrayView1d< localIndex const > const & packList ) const override;

  virtual localIndex unpackUpDownMaps( buffer_unit_type const * & buffer,
                                       array1d< localIndex > & packList,
                                       bool const overwriteUpMaps,
                                       bool const overwriteDownMaps ) override;

  virtual void fixUpDownMaps( bool const clearIfUnmapped ) final override;

  ///@}

  /**
   * @name Miscellaneous
   */
  ///@{

  /**
   * @brief Helper function to apply a lambda function over all constructive groups
   * @tparam LAMBDA the type of the lambda function
   * @param lambda the lambda function
   */
  template< typename LAMBDA >
  void forMaterials( LAMBDA lambda )
  {

    for( auto & constitutiveGroup : m_constitutiveGrouping )
    {
      lambda( constitutiveGroup );
    }
  }

  ///@}

  /**
   * @brief struct to serve as a container for variable strings and keys
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : public CellBlock::viewKeyStruct
  {
    /// String key for the constitutive point volume fraction
    static constexpr auto constitutivePointVolumeFraction = "ConstitutivePointVolumeFraction";
    /// String key for the derivatives of the shape functions with respect to the reference configuration
    static constexpr auto dNdXString = "dNdX";
    /// String key for the derivative of the jacobian.
    static constexpr auto detJString = "detJ";
    /// String key for the constitutive grouping
    static constexpr auto constitutiveGroupingString = "ConstitutiveGrouping";
    /// String key for the constitutive map
    static constexpr auto constitutiveMapString = "ConstitutiveMap";
    /// String key to embSurfMap
    static constexpr auto toEmbSurfString = "ToEmbeddedSurfaces";

    /// ViewKey for the constitutive grouping
    dataRepository::ViewKey constitutiveGrouping  = { constitutiveGroupingString };
    /// ViewKey for the constitutive map
    dataRepository::ViewKey constitutiveMap       = { constitutiveMapString };
  }
  /// viewKey struct for the CellElementSubRegion class
  m_CellBlockSubRegionViewKeys;

  virtual viewKeyStruct & viewKeys() override { return m_CellBlockSubRegionViewKeys; }
  virtual viewKeyStruct const & viewKeys() const override { return m_CellBlockSubRegionViewKeys; }

  /**
   * @brief @return The array of shape function derivatives.
   */
  array4d< real64 > & dNdX()
  { return m_dNdX; }

  /**
   * @brief @return The array of shape function derivatives.
   */
  arrayView4d< real64 const > dNdX() const
  { return m_dNdX; }

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

  /**
   * @brief @return The sorted array of fractured elements.
   */
  SortedArray< localIndex > & fracturedElementsList()
  { return m_fracturedCells; }

  /**
   * @brief @return The sorted array view of fractured elements.
   */
  SortedArrayView< localIndex const > const fracturedElementsList() const
  { return m_fracturedCells.toViewConst(); }

  /**
   * @brief @return The map to the embedded surfaces
   */
  EmbSurfMapType & embeddedSurfacesList() { return m_toEmbeddedSurfaces; }

  /**
   * @brief @return The map to the embedded surfaces
   */
  EmbSurfMapType const & embeddedSurfacesList() const { return m_toEmbeddedSurfaces; }

  /// Map used for constitutive grouping
  map< string, localIndex_array > m_constitutiveGrouping;

  /// Array of constitutive point volume fraction
  array3d< real64 > m_constitutivePointVolumeFraction;

private:

  /// The array of shape function derivaties.
  array4d< real64 > m_dNdX;

  /// The array of jacobian determinantes.
  array2d< real64 > m_detJ;

  /// Map of unmapped global indices in the element-to-node map
  map< localIndex, array1d< globalIndex > > m_unmappedGlobalIndicesInNodelist;

  /// Map of unmapped global indices in the element-to-face map
  map< localIndex, array1d< globalIndex > > m_unmappedGlobalIndicesInFacelist;

  /// List of fractured elements
  SortedArray< localIndex > m_fracturedCells;

  /// Map from Cell Elements to Embedded Surfaces
  EmbSurfMapType m_toEmbeddedSurfaces;

  /**
   * @brief Pack element-to-node and element-to-face maps
   * @tparam the flag for the bufferOps::Pack function
   * @param buffer the buffer used in the bufferOps::Pack function
   * @param packList the packList used in the bufferOps::Pack function
   * @return the pack size
   */
  template< bool DOPACK >
  localIndex packUpDownMapsPrivate( buffer_unit_type * & buffer,
                                    arrayView1d< localIndex const > const & packList ) const;

};

} /* namespace geosx */

#endif /* GEOSX_MESH_CELLELEMENTSUBREGION_HPP_ */
