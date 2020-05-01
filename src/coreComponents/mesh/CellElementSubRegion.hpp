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

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor for this class.
   * @param[in] name the name of this object manager
   * @param[in] parent the parent Group
   */
  CellElementSubRegion( string const & name, Group * const parent );

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
  void CopyFromCellBlock( CellBlock * source );

  /**
   * @brief Fill the CellElementSubRegion by querying a target set into the faceManager
   * @param[in] faceManager a pointer to the faceManager
   * @param[in] setName a reference to string containing the name of the set
   */
  void ConstructSubRegionFromFaceSet( FaceManager const * const faceManager,
                                      string const & setName );

  ///@}

  /**
   * @name Overridden packing / Unpacking functions
   */
  ///@{

  virtual void ViewPackingExclusionList( SortedArray< localIndex > & exclusionList ) const override;

  virtual localIndex PackUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const override;

  virtual localIndex PackUpDownMaps( buffer_unit_type * & buffer,
                                     arrayView1d< localIndex const > const & packList ) const override;

  virtual localIndex UnpackUpDownMaps( buffer_unit_type const * & buffer,
                                       localIndex_array & packList,
                                       bool const overwriteUpMaps,
                                       bool const overwriteDownMaps ) override;

  virtual void FixUpDownMaps( bool const clearIfUnmapped ) final override;

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
    /// The string key for the constitutive point volume fraction
    static constexpr auto constitutivePointVolumeFraction = "ConstitutivePointVolumeFraction";
    /// The string key for the derivatives of the shape functions with respect to the reference configuration
    static constexpr auto dNdXString = "dNdX";

    /// The string key for the constitutive grouping
    static constexpr auto constitutiveGroupingString = "ConstitutiveGrouping";
    /// The string key for the constitutive map
    static constexpr auto constitutiveMapString = "ConstitutiveMap";

    /// The ViewKey for the constitutive grouping
    dataRepository::ViewKey constitutiveGrouping  = { constitutiveGroupingString };
    /// The ViewKey for the constitutive map
    dataRepository::ViewKey constitutiveMap       = { constitutiveMapString };
    /// The ViewKey for the derivatives of the shape functions with respect to the reference configuration
    dataRepository::ViewKey dNdX                  = { dNdXString };

  }
  /// xxx
  m_CellBlockSubRegionViewKeys;

  virtual viewKeyStruct & viewKeys() override { return m_CellBlockSubRegionViewKeys; }
  virtual viewKeyStruct const & viewKeys() const override { return m_CellBlockSubRegionViewKeys; }

  /// The map used for constitutive grouping
  map< string, localIndex_array > m_constitutiveGrouping;

  /// The array of constitutive point volume fraction
  array3d< real64 > m_constitutivePointVolumeFraction;

  /// The derivatives of the shape functions wrt the reference configuration
  array3d< R1Tensor > m_dNdX;


private:

  /// The map of unmapped global indices in the element-to-node map
  map< localIndex, array1d< globalIndex > > m_unmappedGlobalIndicesInNodelist;

  /// The map of unmapped global indices in the element-to-face map
  map< localIndex, array1d< globalIndex > > m_unmappedGlobalIndicesInFacelist;

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

};

} /* namespace geosx */

#endif /* GEOSX_MESH_CELLELEMENTSUBREGION_HPP_ */
