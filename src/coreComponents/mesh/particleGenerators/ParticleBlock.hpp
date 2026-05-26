/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOS_MESH_PARTICLEBLOCK_HPP_
#define GEOS_MESH_PARTICLEBLOCK_HPP_

#include "ParticleBlockABC.hpp"
#include "mesh/ParticleType.hpp"

#include "dataRepository/Group.hpp"

namespace geos
{

/**
 * This implementation of ParticleBlockABC mainly use the cell patterns/shapes
 * to build all the particle to nodes, faces and edges mappings.
 */
class ParticleBlock : public ParticleBlockABC
{
public:

  /**
   * @name Constructor / Destructor
   */
  ///@{

  ParticleBlock() = delete;

  /**
   * @brief Constructor for this class.
   * @param[in] name the name of this object manager
   * @param[in] parent the parent Group
   */
  ParticleBlock( string const & name, Group * const parent );

  /**
   * @brief Copy constructor.
   * @param[in] init the source to copy
   */
  ParticleBlock( const ParticleBlock & init ) = delete;

  /**
   * @brief Default destructor.
   */
  ~ParticleBlock();

  ///@}

  /**
   * @name Geometry computation / Connectivity
   */
  ///@{

  ///@}
  /**
   * @name Getters / Setters
   */
  ///@{

  /**
   * @brief Defines the underlying particle type (hex, tet...)
   * @param[in] particleType the particle type
   *
   * @note Allocates the values of the particle to nodes, edges, faces accordingly.
   */
  void setParticleType( ParticleType particleType );

  ParticleType getParticleType() const override
  { return m_particleType; }

  array1d< globalIndex > getParticleID() const override
  { return m_particleID; }

  /**
   * @brief Sets the global IDs of particles in this subregion.
   * @param particleID The input list of global IDs
   */
  void setParticleID( array1d< globalIndex > const particleID )
  { m_particleID = particleID; }

  array2d< real64 > getParticleCenter() const override
  { return m_particleCenter; }

  /**
   * @brief Set the list of particle center locations in this subregion.
   * @param particleCenter The input list of particle center coordinates
   */
  void setParticleCenter( array2d< real64 > const particleCenter )
  { m_particleCenter = particleCenter; }

  array2d< real64 > getParticleVelocity() const override
  { return m_particleVelocity; }

  /**
   * @brief Set the list of particle velocities in this subregion.
   * @param particleVelocity The input list of velocities
   */
  void setParticleVelocity( array2d< real64 > const particleVelocity )
  { m_particleVelocity = particleVelocity; }

  array3d< real64 > getParticleMaterialDirection() const override
  { return m_particleMaterialDirection; }

  /**
   * @brief Set the list of material directions in this subregion.
   * @param particleMaterialDirection The input list of directions
   */
  void setParticleMaterialDirection( array3d< real64 > const particleMaterialDirection )
  { m_particleMaterialDirection = particleMaterialDirection; }

  array1d< int > getParticleGroup() const override
  { return m_particleGroup; }

  /**
   * @brief Set the list of particle group numbers (for contact) in this subregion.
   * @param particleGroup The input list of contact group numbers.
   */
  void setParticleGroup( array1d< int > const particleGroup )
  { m_particleGroup = particleGroup; } // TODO: Rename to particleContactGroup

  array1d< integer > getParticleSurfaceFlag() const override
  { return m_particleSurfaceFlag; }

  /**
   * @brief Set the list of particle surface flags in this subregion.
   * @param particleSurfaceFlag The input list of surface flags.
   */
  void setParticleSurfaceFlag( array1d< integer > const particleSurfaceFlag )
  { m_particleSurfaceFlag = particleSurfaceFlag; }

  array1d< real64 > getParticleDamage() const override
  { return m_particleDamage; }

  /**
   * @brief Set the list of particle damage values in this subregion.
   * @param particleDamage The input list of damage values
   */
  void setParticleDamage( array1d< real64 > const particleDamage )
  { m_particleDamage = particleDamage; }

  array1d< real64 > getParticlePorosity() const override
  { return m_particlePorosity; }

  /**
   * @brief Set the list of particle porosity values in this subregion.
   * @param particlePorosity The input list of porosity values
   */
  void setParticlePorosity( array1d< real64 > const particlePorosity )
  { m_particlePorosity = particlePorosity; }

  array1d< real64 > getParticleTemperature() const override
  { return m_particleTemperature; }

  /**
   * @brief Set the list of particle temperature values in this subregion.
   * @param particleTemperature The input list of temperature values
   */
  void setParticleTemperature( array1d< real64 > const particleTemperature )
  { m_particleTemperature = particleTemperature; }

  array1d< real64 > getParticleTemperatureRate() const override
  { return m_particleTemperatureRate; }

  /**
   * @brief Set the list of particle temperature rate values in this subregion.
   * @param particleTemperatureRate The input list of temperature rate values
   */
  void setParticleTemperatureRate( array1d< real64 > const particleTemperatureRate )
  { m_particleTemperatureRate = particleTemperatureRate; }

  array1d< real64 > getParticleStrengthScale() const override
  { return m_particleStrengthScale; }

  /**
   * @brief Set the list of particle strength scale values in this subregion.
   * @param particleStrengthScale The input list of strength scale values
   */
  void setParticleStrengthScale( array1d< real64 > const particleStrengthScale )
  { m_particleStrengthScale = particleStrengthScale; }

  array1d< real64 > getParticleVolume() const override
  { return m_particleVolume; }

  /**
   * @brief Set the list of particle volumes in this subregion.
   * @param particleVolume The input list of volumes
   */
  void setParticleVolume( array1d< real64 > const particleVolume )
  { m_particleVolume = particleVolume; }

  array3d< real64 > getParticleRVectors() const override
  { return m_particleRVectors; }

  /**
   * @brief Set the list of particle r-vectors in this subregion.
   * @param particleRVectors The input list of r-vectors
   */
  void setParticleRVectors( array3d< real64 > const particleRVectors )
  { m_particleRVectors = particleRVectors; }

  bool hasRVectors() const override
  { return m_hasRVectors; }

  array2d< real64 > getParticleSurfaceNormal() const override
  { return m_particleSurfaceNormal; }

  /**
   * @brief Set the list of surface normals in this subregion.
   * @param particleSurfaceNormal The input list of normals
   */
  void setParticleSurfaceNormal( array2d< real64 > const particleSurfaceNormal )
  { m_particleSurfaceNormal = particleSurfaceNormal; }

  array2d< real64 > getParticleSurfacePosition() const override
  { return m_particleSurfacePosition; }

  /**
   * @brief Set the list of surface positions in this subregion.
   * @param particleSurfacePosition The input list of surface positions
   */
  void setParticleSurfacePosition( array2d< real64 > const particleSurfacePosition )
  { m_particleSurfacePosition = particleSurfacePosition; }

  array2d< real64 > getParticleSurfaceTraction() const override
  { return m_particleSurfaceTraction; }

  /**
   * @brief Set the list of surface tractions in this subregion.
   * @param particleSurfaceTraction The input list of surface positions
   */
  void setParticleSurfaceTraction( array2d< real64 > const particleSurfaceTraction )
  { m_particleSurfaceTraction = particleSurfaceTraction; }

  array1d< localIndex > getParticleCZTag() const override
  { return m_particleCZTag; }

  /**
   * @brief Set the list of cohesive zone tags in this subregion.
   * @param particleCZTag The input list of cohesive zone tags
   */
  void setParticleCZTag( array1d< localIndex > const particleCZTag )
  { m_particleCZTag = particleCZTag; }

  localIndex numParticles() const override
  { return size(); }

  /**
   * @brief Get local to global map, non-const version.
   * @return The mapping relationship as a array.
   *
   * @deprecated This accessor is meant to be used like a setter even though it's a bit like having public attribute...
   * Use a real setter instead.
   */
  arrayView1d< globalIndex > localToGlobalMap()
  { return m_localToGlobalMap; }

  array1d< globalIndex > localToGlobalMap() const override
  { return m_localToGlobalMap; }

  /**
   * @brief Resize the cell block to hold @p numParticles
   * @param numParticles The new number of particles.
   */
  void resize( dataRepository::indexType const numParticles ) override final;

  ///@}

  /**
   * @name Properties
   */
  ///@{

  ///@}

private:

  /// Contains the global index of each object.
  array1d< globalIndex > m_localToGlobalMap;

  /// Name of the properties registered from an external mesh
  string_array m_externalPropertyNames;

  /// Member level field for the particle global ID.
  array1d< globalIndex > m_particleID;

  /// Member level field for the particle contact group.
  array1d< int > m_particleGroup;

  /// Member level field for the particle surface flag.
  array1d< integer > m_particleSurfaceFlag;

  /// Member level field for the particle damage.
  array1d< real64 > m_particleDamage;

  /// Member level field for the particle porosity.
  array1d< real64 > m_particlePorosity;

  /// Member level field for the particle temperature.
  array1d< real64 > m_particleTemperature;

  /// Member level field for the particle temperature.
  array1d< real64 > m_particleTemperatureRate;

  /// Member level field for the particle strength scale.
  array1d< real64 > m_particleStrengthScale;

  /// Member level field for the particle center.
  array2d< real64 > m_particleCenter;

  /// Member level field for the particle velocity.
  array2d< real64 > m_particleVelocity;

  /// Member level field for the particle material direction.
  array3d< real64 > m_particleMaterialDirection;

  /// Member level field for the particle volume.
  array1d< real64 > m_particleVolume;

  /// Type of particles in this subregion.
  ParticleType m_particleType;

  /// Bool flag for whether particles in this block have r-vectors defining its domain
  bool m_hasRVectors;

  /// R-vectors
  array3d< real64 > m_particleRVectors;

  /// Member level field for the particle surface normal.
  array2d< real64 > m_particleSurfaceNormal;

  /// Member level field for the particle surface position.
  array2d< real64 > m_particleSurfacePosition;

  /// Member level field for the particle surface traction.
  array2d< real64 > m_particleSurfaceTraction;

  /// Member level field for the particle cohesive zone tag.
  array1d< localIndex > m_particleCZTag;

  std::list< dataRepository::WrapperBase * > getExternalProperties() override
  {
    std::list< dataRepository::WrapperBase * > result;
    for( string const & externalPropertyName : m_externalPropertyNames )
    {
      result.push_back( &this->getWrapperBase( externalPropertyName ) );
    }
    return result;
  }
};

}

#endif /* GEOS_MESH_CELLBLOCK_HPP_ */
