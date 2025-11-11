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

/**
 * @file CohesiveZoneRegion.hpp
 *
 */

#ifndef GEOS_MESH_COHESIVEREGION_HPP_
#define GEOS_MESH_COHESIVEREGION_HPP_

#include "CohesiveZoneRegionBase.hpp"

namespace geos
{
/**
 * @class CohesiveZoneRegion
 *
 * The CohesiveZoneRegion class contains the functionality to support the concept of a CohesiveZoneRegion in the element
 * hierarchy. CohesiveZoneRegion derives from CohesiveZoneRegionBase and has an entry in the ObjectManagerBase catalog.
 *
 *
 */
class CohesiveZoneRegion : public CohesiveZoneRegionBase
{
public:

  /**
   * @name Constructor / Destructor
   */
  ///@{


  /**
   * @brief Constructor.
   * @param name the name of the object in the data hierarchy.
   * @param parent a pointer to the parent group in the data hierarchy.
   */
  CohesiveZoneRegion( string const & name, Group * const parent );

  /**
   * @brief Deleted default constructor.
   */
  CohesiveZoneRegion() = delete;

  /**
   * @brief Destructor.
   */
  virtual ~CohesiveZoneRegion() override;

  /**
   * @name Static factory catalog functions
   */
  ///@{

  /**
   * @brief The key name for the FaceElementRegion in the object catalog.
   * @return A string containing the key name.
   */
  static string catalogName()
  { return "CohesiveZoneRegion"; }

  /**
   * @copydoc catalogName()
   */
  virtual string getCatalogName() const override final
  { return CohesiveZoneRegion::catalogName(); }

  ///@}

  /**
   * @name Getters
   */
  ///@{

  /**
   * @brief Get the global ID of each cohesive zone node.
   * @return an arrayView1d of const node global ID
   */
  arrayView1d< globalIndex const > getGlobalID() const
  { return m_globalID; }

  /**
   * @copydoc getGlobalID() const
   */
  arrayView1d< globalIndex > getGlobalID()
  { return m_globalID; }

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
  arrayView2d< real64 const > getReferenceSurfaceNormal() const
  { return m_referenceSurfaceNormal; }

  /**
   * @copydoc getReferenceSurfaceNormal() const
   */
  arrayView2d< real64 > getReferenceSurfaceNormal()
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

  /**
   * @brief Get the max normal displacement of each cohesive zone node.
   * @return an arrayView2d of const node max normal displacement
   */
  arrayView2d< real64 const > getMaxNormalDisplacement() const
  { return m_maxNormalDisplacement; }

  /**
   * @copydoc getMaxNormalDisplacement() const
   */
  arrayView2d< real64 > getMaxNormalDisplacement()
  { return m_maxNormalDisplacement; }

  /**
   * @brief Get the max tangential displacement of each cohesive zone node.
   * @return an arrayView2d of const node max tangential displacement
   */
  arrayView2d< real64 const > getMaxTangentialDisplacement() const
  { return m_maxTangentialDisplacement; }

  /**
   * @copydoc getMaxTangentialDisplacement() const
   */
  arrayView2d< real64 > getMaxTangentialDisplacement()
  { return m_maxTangentialDisplacement; }

  /**
   * @brief Get the damage of each cohesive zone node.
   * @return an arrayView2d of const node damage
   */
  arrayView2d< real64 const > getDamage() const
  { return m_damage; }

  /**
   * @copydoc getDamage() const
   */
  arrayView2d< real64 > getDamage()
  { return m_damage; }

  ///@}


  /**
   * @brief A struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : public CohesiveZoneRegionBase::viewKeyStruct
  {
    /// @return String key for the member level field for the cohesive zone node global ID.
    static constexpr char const * globalIDString() { return "globalID"; }

    /// @return String key for the member level field for the cohesive zone node reference partitioning surface normals.
    static constexpr char const * referencePartitioningSurfaceNormalString() { return "referencePartitioningSurfaceNormal"; }

    /// @return String key for the member level field for the cohesive zone node reference surface normals.
    static constexpr char const * referenceSurfaceNormalString() { return "referenceSurfaceNormal"; }

    /// @return String key for the member level field for the cohesive zone node reference areas.
    static constexpr char const * referenceAreaString() { return "referenceArea"; }

    /// @return String key for the member level field for the cohesive zone node reference position.
    static constexpr char const * referencePositionString() { return "referencePosition"; }

    /// @return String key for the member level field for the cohesive zone node max normal displacement.
    static constexpr char const * maxNormalDisplacementString() { return "maxNormalDisplacement"; }

    /// @return String key for the member level field for the cohesive zone node max tangential displacement.
    static constexpr char const * maxTangentialDisplacementString() { return "maxTangentialDisplacement"; }

    /// @return String key for the member level field for the cohesive zone node damages.
    static constexpr char const * damageString() { return "damage"; }
  };

private:

  // Cohesive zone parameters

  array1d< globalIndex > m_globalID;

  array2d< real64 > m_referencePartitioningSurfaceNormal;
  array3d< real64 > m_referenceSurfaceNormal;
  array2d< real64 > m_referenceArea;
  array2d< real64 > m_referencePosition;
  array1d< real64 > m_maxNormalDisplacement;
  array1d< real64 > m_maxTangentialDisplacement;
  array1d< real64 > m_damage;

};

} /* namespace geos */

#endif /* GEOS_MESH_COHESIVEZONEREGION_HPP_ */
