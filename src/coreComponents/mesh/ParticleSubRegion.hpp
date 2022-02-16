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


#ifndef GEOSX_MESH_PARTICLEELEMENTSUBREGION_HPP_
#define GEOSX_MESH_PARTICLEELEMENTSUBREGION_HPP_

#include "mesh/generators/ParticleBlockABC.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"
#include "ParticleSubRegionBase.hpp"


namespace geosx
{

class MeshLevel;

/**
 * @class ParticleSubRegion
 * Class specializing the particle subregion for a cell particle.
 * This is the class used in the physics solvers to represent a collection of mesh cell particles
 */
class ParticleSubRegion : public ParticleSubRegionBase
{
public:

  /**
   * @brief Const getter for the catalog name.
   * @return the name of this type in the catalog
   */
  static const string catalogName()
  { return "ParticleSubRegion"; }

  /**
   * @copydoc catalogName()
   */
  virtual const string getCatalogName() const override final
  { return ParticleSubRegion::catalogName(); }

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor for this class.
   * @param[in] name the name of this object manager
   * @param[in] parent the parent Group
   */
  ParticleSubRegion( string const & name, Group * const parent );

  /**
   * @brief Destructor.
   */
  virtual ~ParticleSubRegion() override;

  ///@}

  /**
   * @name Helpers for ParticleSubRegion construction
   */
  ///@{

  /**
   * @brief Fill the ParticleSubRegion by copying those of the source ParticleBlock
   * @param particleBlock the ParticleBlock which properties (connectivity info) will be copied.
   */
  void copyFromParticleBlock( ParticleBlockABC & particleBlock );

  ///@}

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
  struct viewKeyStruct : public ParticleSubRegionBase::viewKeyStruct
  {
    /// @return String key for the constitutive point volume fraction
    static constexpr char const * constitutivePointVolumeFractionString() { return "ConstitutivePointVolumeFraction"; }
    /// @return String key for the derivatives of the shape functions with respect to the reference configuration
    static constexpr char const * dNdXString() { return "dNdX"; }
    /// @return String key for the derivative of the jacobian.
    static constexpr char const * detJString() { return "detJ"; }
    /// @return String key for the constitutive grouping
    static constexpr char const * constitutiveGroupingString() { return "ConstitutiveGrouping"; }
    /// @return String key for the constitutive map
    static constexpr char const * constitutiveMapString() { return "ConstitutiveMap"; }

    /// ViewKey for the constitutive grouping
    dataRepository::ViewKey constitutiveGrouping  = { constitutiveGroupingString() };
    /// ViewKey for the constitutive map
    dataRepository::ViewKey constitutiveMap       = { constitutiveMapString() };
  }
  /// viewKey struct for the ParticleSubRegion class
  m_ParticleBlockSubRegionViewKeys;

  virtual viewKeyStruct & viewKeys() override { return m_ParticleBlockSubRegionViewKeys; }
  virtual viewKeyStruct const & viewKeys() const override { return m_ParticleBlockSubRegionViewKeys; }

  /**
   * @brief Get the local indices of the nodes in a face of the particle.
   * @param[in] particleIndex The local index of the target particle.
   * @param[in] localFaceIndex The local index of the target face in the particle (this will be [0, numFacesInParticle[)
   * @param[out] nodeIndices A reference to the array of node indices of the face. Gets resized at the proper size.
   * @deprecated This method will be removed soon.
   */
  void getFaceNodes( localIndex const particleIndex,
                     localIndex const localFaceIndex,
                     array1d< localIndex > & nodeIndices ) const;





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

private:

  /// Map used for constitutive grouping
  map< string, localIndex_array > m_constitutiveGrouping;

  /// Array of constitutive point volume fraction
  array3d< real64 > m_constitutivePointVolumeFraction;

  /// Name of the properties registered from an external mesh
  string_array m_externalPropertyNames;

  /// The array of shape function derivaties.
  array4d< real64 > m_dNdX;

  /// The array of jacobian determinantes.
  array2d< real64 > m_detJ;

  /**
   * @brief Pack particle-to-node and particle-to-face maps
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
