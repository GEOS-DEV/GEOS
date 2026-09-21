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
 * @brief Manages sparse cohesive physical-node, field-slot, and field-pair data.
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

  /**
   * @brief Get the compact field-slot indices on the two sides of each cohesive pair.
   * @return An array with dimensions @c {numberOfPairs, 2}.
   */
  arrayView2d< localIndex const > getFieldSlots() const
  { return m_fieldSlots; }

  /** @copydoc getFieldSlots() const */
  arrayView2d< localIndex > getFieldSlots()
  { return m_fieldSlots; }

  /** @brief Get the compact physical-node index for every field slot. */
  arrayView1d< localIndex const > getFieldSlotNode() const
  { return m_fieldSlotNode; }

  /** @copydoc getFieldSlotNode() const */
  arrayView1d< localIndex > getFieldSlotNode()
  { return m_fieldSlotNode; }

  /** @brief Get the global velocity-field index for every compact field slot. */
  arrayView1d< localIndex const > getFieldSlotVelocityField() const
  { return m_fieldSlotVelocityField; }

  /** @copydoc getFieldSlotVelocityField() const */
  arrayView1d< localIndex > getFieldSlotVelocityField()
  { return m_fieldSlotVelocityField; }

  /**
   * @brief Resize the independently-sized physical-node and field-slot arrays.
   * @param numNodes Number of unique physical cohesive nodes.
   * @param numFieldSlots Number of active node/velocity-field combinations.
   */
  void resizeTopology( localIndex const numNodes, localIndex const numFieldSlots )
  {
    m_globalID.resize( numNodes );
    m_referencePosition.resize( numNodes, 3 );
    m_referencePartitioningSurfaceNormal.resize( numNodes, 3 );
    m_fieldSlotNode.resize( numFieldSlots );
    m_fieldSlotVelocityField.resize( numFieldSlots );
    m_referenceSurfaceNormal.resize( numFieldSlots, 3 );
    m_referenceArea.resize( numFieldSlots );
  }

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
  int isInitialized() const { return m_initialized; }

  void setEnabled( int const & enabled ) { m_enabled = enabled; }
  int isEnabled() const { return m_enabled; }

  /**
   * @brief Get the sorted global IDs of unique physical cohesive nodes.
   * @return an array view of physical grid-node global IDs
   */
  arrayView1d< globalIndex const > getGlobalID() const
  { return m_globalID.toViewConst(); }

  /** @copydoc getGlobalID() const */
  arrayView1d< globalIndex > getGlobalID()
  { return m_globalID.toView(); }

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
   * @brief Get the reference surface normal of each compact field slot.
   * @return an array view with dimensions @c {numberOfFieldSlots, 3}
   */
  arrayView2d< real64 const > getReferenceSurfaceNormal() const
  { return m_referenceSurfaceNormal; }

  /**
   * @copydoc getReferenceSurfaceNormal() const
   */
  arrayView2d< real64 > getReferenceSurfaceNormal()
  { return m_referenceSurfaceNormal; }

  /**
   * @brief Get the reference area of each compact field slot.
   * @return an array view with one scalar per field slot
   */
  arrayView1d< real64 const > getReferenceArea() const
  { return m_referenceArea; }

  /** @copydoc getReferenceArea() const */
  arrayView1d< real64 > getReferenceArea()
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

    /// @return String key for the compact field slots forming each cohesive pair.
    static constexpr char const * fieldSlotsString() { return "fieldSlots"; }

    /// @return String key for the physical-node index of each compact field slot.
    static constexpr char const * fieldSlotNodeString() { return "fieldSlotNode"; }

    /// @return String key for the global velocity-field index of each compact field slot.
    static constexpr char const * fieldSlotVelocityFieldString() { return "fieldSlotVelocityField"; }

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

  int m_enabled;
  int m_initialized;
  
  int m_czVolumeNormalization;
  int m_computeParticleSurfaceNormalsAndPositions;
  int m_normalsAndPositionsMethod; // Should be the enum from SolidMechanicsMPM, but currently circular dependences that needs to be resolved
  int m_czSurfaceDisplacementUpdate; // Stored as int to avoid a mesh-to-solver dependency cycle.

  string m_constitutiveModelName;

  localIndex m_tag;

  // Sparse topology. Each cohesive constitutive point is a binary pair, while
  // field slots deduplicate the kinematic and surface data shared at junctions.
  array2d< localIndex > m_fieldSlots;
  array1d< localIndex > m_fieldSlotNode;
  array1d< localIndex > m_fieldSlotVelocityField;

  // Physical-node reference fields (one entry per unique physical node).
  array1d< globalIndex > m_globalID;
  array2d< real64 > m_referencePosition;
  array2d< real64 > m_referencePartitioningSurfaceNormal;

  // Side geometry (one entry per active node/velocity-field slot).
  array2d< real64 > m_referenceSurfaceNormal;
  array1d< real64 > m_referenceArea;
};

}



#endif /* GEOS_MESH_COHESIVEZONEREGIONBASE_HPP */
