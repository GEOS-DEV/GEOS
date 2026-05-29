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

#ifndef GEOS_MESH_COHESIVEZONEREGIONBASE_HPP
#define GEOS_MESH_COHESIVEZONEREGIONBASE_HPP

#include "mesh/ObjectManagerBase.hpp"

namespace geos
{

namespace constitutive
{
class CohesiveZoneBase;
}


/**
 * @class CohesiveZoneRegionBase
 * @brief The CohesiveZoneRegionBase is the base class to manage the cohesive zone data stored at the node level.
 *
 * The CohesiveZoneRegionBase is the base class for the CohesiveZoneRegion class. It may be depreciated at
 * some point since no other classes are currently derived from CohesiveZoneRegionBase.
 */
class CohesiveZoneRegionBase : public ObjectManagerBase
{
public:

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Deleted default constructor.
   */
  CohesiveZoneRegionBase() = delete;

  /**
   * @brief Main constructor.
   * @param name the name of the particle region
   * @param parent the pointer to the parent group
   */
  CohesiveZoneRegionBase( string const & name, Group * const parent );


  /**
   * @brief Copy constructor.
   * @param init the particle region to be copied
   */
  CohesiveZoneRegionBase( const CohesiveZoneRegionBase & init );

  /**
   * @brief Default destructor.
   */
  virtual ~CohesiveZoneRegionBase() override;

  ///@}

  /**
   * @name Generation of the mesh region
   */
  ///@{

  /**
   * @brief Generate mesh.
   * @param blocks blocks where the mesh is generated
   */
  virtual void generateMesh( Group & blocks ) final
  {
    GEOS_UNUSED_VAR( blocks );
    GEOS_ERROR( "CohesiveZoneRegionBase::GenerateMesh() is non-op." );
  }

  ///@}

  /**
   * @name Getters / Setters
   */
  ///@{

  void setTag( localIndex tag ) { m_tag = tag; }

  localIndex getTag() const { return m_tag; }

  localIndex getFieldA() const { return m_fieldA; }

  localIndex getFieldB() const { return m_fieldB; }

  /**
   * @brief Get a pointer to the constitutive model.
   * @tparam T The type of the constitutive model.
   * @param name The name of the constitutive model.
   * @return A pointer to the constitutive model.
   */
  template< typename T = constitutive::CohesiveZoneBase >
  T const & getConstitutiveModel() const
  { return this->getGroup< T >( m_constitutiveModelName ); }

  /**
   * @copydoc getConstitutiveModel() const
   */
  template< typename T = constitutive::CohesiveZoneBase >
  T & getConstitutiveModel()
  { return this->getGroup< T >( m_constitutiveModelName ); }

  void setConstitutiveModelName( string const & constitutiveModelName ) { m_constitutiveModelName = constitutiveModelName; }
  string const & getConstitutiveModelName() const { return m_constitutiveModelName; }

  void setInitialized( int const & initialized ) { m_initialized = initialized; }
  int isInitalized() const { return m_initialized; }

  /**
   * @brief Get the global ID of each cohesive zone node.
   * @return an arrayView1d of const node global ID
   */
  SortedArrayView< globalIndex const > getGlobalID() const
  { return m_globalID.toViewConst(); }

  void setGlobalID(SortedArrayView< globalIndex const > const & globalID )
  {
    m_globalID.insert( globalID.begin(), globalID.end() );
  }

  void setCZVolumeNormalization( int const & czVolumeNormalization ) { m_czVolumeNormalization = czVolumeNormalization; }
  void setComputeParticleSurfaceNormalsAndPositions( int const & computeParticleSurfaceNormalsAndPositions ) { m_computeParticleSurfaceNormalsAndPositions = computeParticleSurfaceNormalsAndPositions; }
  void setNormalsAndPositionsMethod( int const & normalsAndPositionsMethod ) { m_normalsAndPositionsMethod = normalsAndPositionsMethod; }
  void setCZSurfaceDisplacementUpdate( int const & czSurfaceDisplacementUpdate ) { m_czSurfaceDisplacementUpdate = czSurfaceDisplacementUpdate; }

  int getCZVolumeNormalization() const { return m_czVolumeNormalization; }
  int getComputeNormalsAndPositions() const { return m_computeParticleSurfaceNormalsAndPositions; }
  int getNormalsAndPositionsMethod() const { return m_normalsAndPositionsMethod; }
  int getCZSurfaceDisplacementUpdate() const { return m_czSurfaceDisplacementUpdate; }

  /**
   * @brief Get the reference partitioning surface normal of each cohesive zone node.
   * @return an arrayView2d of const node reference partitioning surface normal
   */
  arrayView2d< real64 const > getReferencePartitioningSurfaceNormal() const
  { return m_referencePartitioningSurfaceNormal; }

  /**
   * @copydoc getReferencePartitioningSurfaceNormal() const
   */
  arrayView2d< real64 > getReferencePartitioningSurfaceNormal()
  { return m_referencePartitioningSurfaceNormal; }

  /**
   * @brief Get the reference surface normal of each cohesive zone node.
   * @return an arrayView2d of const node reference surface normal
   */
  arrayView3d< real64 const > getReferenceSurfaceNormal() const
  { return m_referenceSurfaceNormal; }

  /**
   * @copydoc getReferenceSurfaceNormal() const
   */
  arrayView3d< real64 > getReferenceSurfaceNormal()
  { return m_referenceSurfaceNormal; }

  /**
   * @brief Get the reference area of each cohesive zone node.
   * @return an arrayView2d of const node reference area
   */
  arrayView2d< real64 const > getReferenceArea() const
  { return m_referenceArea; }

  /**
   * @copydoc getReferencePosition() const
   */
  arrayView2d< real64 > getReferenceArea()
  { return m_referenceArea; }

  /**
   * @brief Get the reference position of each cohesive zone node.
   * @return an arrayView2d of const node reference position
   */
  arrayView2d< real64 const > getReferencePosition() const
  { return m_referencePosition; }

  /**
   * @copydoc getReferencePosition() const
   */
  arrayView2d< real64 > getReferencePosition()
  { return m_referencePosition; }

  ///@}

  /**
   * @brief Struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : public ObjectManagerBase::viewKeyStruct
  {
    /// @return String key for the member level field for the cohesive zone tag.
    static constexpr char const * tagString() { return "tag"; }

    /// @return String key for the member level field for the cohesive zone node global ID.
    static constexpr char const * globalIDString() { return "globalID"; }

    /// @return String key for the member level field for the cohesive zone node reference partitioning surface normals.
    static constexpr char const * referencePartitioningSurfaceNormalString() { return "referencePartitioningSurfaceNormal"; }

    /// @return String key for the member level field for the cohesive zone node reference surface normals.
    static constexpr char const * referenceSurfaceNormalString() { return "referenceSurfaceNormal"; }

    /// @return String key for the member level field for the cohesive zone node reference areas.
    static constexpr char const * referenceAreaString() { return "referenceArea"; }

    /// @return String key for the cohesive-zone surface displacement update method.
    static constexpr char const * czSurfaceDisplacementUpdateString() { return "czSurfaceDisplacementUpdate"; }

    /// @return String key for the member level field for the cohesive zone node reference position.
    static constexpr char const * referencePositionString() { return "referencePosition"; }
  };

private:

  CohesiveZoneRegionBase & operator=( const CohesiveZoneRegionBase & rhs );

  int m_initialized;
  int m_czVolumeNormalization;
  int m_computeParticleSurfaceNormalsAndPositions;
  int m_normalsAndPositionsMethod; // Should be the enum from SolidMechanicsMPM, but currently circular dependences that needs to be resolved
  int m_czSurfaceDisplacementUpdate; // Should be the enum from SolidMechanicsMPM, but currently circular dependency prevents that

  string m_constitutiveModelName;

  localIndex m_tag;

  // Indices of fields for either side of the cohesive zone
  localIndex m_fieldA;
  localIndex m_fieldB;

  // Reference fields
  SortedArray< globalIndex > m_globalID;
  array2d< real64 > m_referencePosition;
  array2d< real64 > m_referencePartitioningSurfaceNormal;
  array3d< real64 > m_referenceSurfaceNormal;
  array2d< real64 > m_referenceArea;
};

}



#endif /* GEOS_MESH_COHESIVEZONEREGIONBASE_HPP */
