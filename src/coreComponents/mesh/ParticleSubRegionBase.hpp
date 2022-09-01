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
 * @file ParticleSubRegionBase.hpp
 */

#ifndef GEOSX_MESH_PARTICLESUBREGIONBASE_HPP_
#define GEOSX_MESH_PARTICLESUBREGIONBASE_HPP_

#include "mesh/ParticleType.hpp"
#include "mesh/ObjectManagerBase.hpp"

namespace geosx
{

class ParticleManager;
class MeshLevel;

namespace constitutive
{
class ConstitutiveBase;
}

/**
 * @class ParticleSubRegionBase
 * Abstract class for a collection of mesh particles that
 * will be derived and specialized for particles
 */
class ParticleSubRegionBase : public ObjectManagerBase
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
  ParticleSubRegionBase( string const & name, dataRepository::Group * const parent );

  /**
   * @brief Destructor.
   */
  ~ParticleSubRegionBase();

  ///@}

  /**
   * @name Getters / Setters
   */
  ///@{

  /**
   * @brief Get whether particle has r-vectors.
   */
  bool hasRVectors() const
  { return m_hasRVectors; }

  /**
   * @brief Set whether particle has r-vectors.
   * @param[in] Boolean indicating whether particle has r-vectors.
   */
  void setHasRVectors( bool hasRVectors )
  { m_hasRVectors = hasRVectors; }

  /**
   * @brief Get the global ID of each particle in this subregion.
   * @return an arrayView1d of const particle global IDs
   */
  arrayView1d< int const > getParticleID() const
  { return m_particleID; }

  /**
   * @copydoc getParticleID() const
   */
  arrayView1d< int > getParticleID()
  { return m_particleID; }

    /**
   * @brief Get the contact group of each particle in this subregion.
   * @return an arrayView1d of const particle contact groups
   */
  arrayView1d< int const > getParticleGroup() const
  { return m_particleGroup; }

  /**
   * @copydoc getParticleGroup() const
   */
  arrayView1d< int > getParticleGroup()
  { return m_particleGroup; }

  /**
   * @brief Get the ghost rank of each particle in this subregion.
   * @return an arrayView1d of const particle ghost ranks
   */
  arrayView1d< int const > getParticleRank() const
  { return m_particleRank; }

  /**
   * @copydoc getParticleRank() const
   */
  arrayView1d< int > getParticleRank()
  { return m_particleRank; }

  /**
   * @brief Get the center of each particle in this subregion.
   * @return an arrayView1d of const particle centers
   */
  arrayView2d< real64 const > getParticleCenter() const
  { return m_particleCenter; }

  /**
   * @copydoc getParticleCenter() const
   */
  arrayView2d< real64 > getParticleCenter()
  { return m_particleCenter; }

  /**
   * @brief Set the center of each particle in this subregion.
   */
  void setParticleCenter( array2d< real64 > particleCenter )
  { m_particleCenter = particleCenter; }

  /**
   * @brief Get the velocity of each particle in this subregion.
   * @return an arrayView1d of const particle velocities
   */
  arrayView2d< real64 const > getParticleVelocity() const
  { return m_particleVelocity; }

  /**
   * @copydoc getParticleVelocity() const
   */
  arrayView2d< real64 > getParticleVelocity()
  { return m_particleVelocity; }

  /**
   * @brief Get the current volume of each particle in this subregion.
   * @return an arrayView1d of const current particle volumes
   */
  arrayView1d< real64 const > getParticleVolume() const
  { return m_particleVolume; }

  /**
   * @copydoc getParticleVolume() const
   */
  arrayView1d< real64 > getParticleVolume()
  { return m_particleVolume; }

  /**
   * @brief Get the initial volume of each particle in this subregion.
   * @return an arrayView1d of const initial particle volumes
   */
  arrayView1d< real64 const > getParticleInitialVolume() const
  { return m_particleInitialVolume; }

  /**
   * @copydoc getParticleInitialVolume() const
   */
  arrayView1d< real64 > getParticleInitialVolume()
  { return m_particleInitialVolume; }

  /**
   * @brief Get the mass of each particle in this subregion.
   * @return an arrayView1d of const particle masses
   */
  arrayView1d< real64 const > getParticleMass() const
  { return m_particleMass; }

  /**
   * @copydoc getParticleMass() const
   */
  arrayView1d< real64 > getParticleMass()
  { return m_particleMass; }

  /**
   * @brief Get the deformation gradient of each particle in this subregion.
   * @return an arrayView3d of const particle deformation gradients
   */
  arrayView3d< real64 const > getParticleDeformationGradient() const
  { return m_particleDeformationGradient; }

  /**
   * @copydoc getParticleDeformationGradient() const
   */
  arrayView3d< real64 > getParticleDeformationGradient()
  { return m_particleDeformationGradient; }

  /**
   * @brief Get the r-vectors of each particle in this subregion.
   * @return an arrayView3d of const particle r-vectors
   */
  arrayView3d< real64 const > getParticleRVectors() const
  { return m_particleRVectors; }

  /**
   * @copydoc getParticleRVectors() const
   */
  arrayView3d< real64 > getParticleRVectors()
  { return m_particleRVectors; }


  /**
   * @brief Get the group in which the constitutive models of this subregion are registered.
   * @return a pointer to the const group in which the constitutive models are registered
   */
  dataRepository::Group const & getConstitutiveModels() const
  { return m_constitutiveModels; }

  /**
   * @copydoc getConstitutiveModels() const
   */
  dataRepository::Group & getConstitutiveModels()
  { return m_constitutiveModels; }

  /**
   * @brief Get a pointer to the constitutive model.
   * @tparam T The type of the constitutive model.
   * @param name The name of the constitutive model.
   * @return A pointer to the constitutive model.
   */
  template< typename T = constitutive::ConstitutiveBase >
  T const & getConstitutiveModel( string const & name ) const
  { return m_constitutiveModels.getGroup< T >( name ); }

  /**
   * @copydoc getConstitutiveModel( string const & ) const
   */
  template< typename T = constitutive::ConstitutiveBase >
  T & getConstitutiveModel( string const & name )
  { return m_constitutiveModels.getGroup< T >( name ); }


  /**
   * @brief Get the type of particle in this subregion.
   * @return the type of particle in this subregion
   */
  ParticleType getParticleType() const
  { return m_particleType; }

  /**
   * @brief Set the type of particle in this subregion.
   * @param[in] particleType the particle type
   */
  virtual void setParticleType( ParticleType const particleType )
  { m_particleType = particleType; }

  ///@}

  /**
   * @brief A struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : ObjectManagerBase::viewKeyStruct
  {
    /// @return String key for the member level field for the particle global ID.
    static constexpr char const * particleIDString() { return "particleID"; }

    /// @return String key for the member level field for the particle contact group.
    static constexpr char const * particleGroupString() { return "particleGroup"; }

    /// @return String key for the member level field for the particle ghost rank.
    static constexpr char const * particleRankString() { return "particleRank"; }

    /// @return String key for the member level field for the particle center.
    static constexpr char const * particleCenterString() { return "particleCenter"; }

    /// @return String key for the member level field for the particle velocity.
    static constexpr char const * particleVelocityString() { return "particleVelocity"; }

    /// @return String key for the member level field for the current particle volume.
    static constexpr char const * particleVolumeString() { return "particleVolume"; }

    /// @return String key for the member level field for the initial particle volume.
    static constexpr char const * particleInitialVolumeString() { return "particleInitialVolume"; }

    /// @return String key for the member level field for the particle volume.
    static constexpr char const * particleMassString() { return "particleMass"; }

    /// @return String key for the member level field for the particle volume.
    static constexpr char const * particleDeformationGradientString() { return "particleDeformationGradient"; }

    /// @return String key for the member level field for the particle volume.
    static constexpr char const * particleRVectorsString() { return "particleRVectors"; }

    /// @return String key for the member level field for the particle volume.
    static constexpr char const * particleInitialRVectorsString() { return "particleInitialRVectors"; }
  };

  /**
   * @brief struct to serve as a container for group strings and keys
   * @struct groupKeyStruct
   */
  struct groupKeyStruct : public ObjectManagerBase::groupKeyStruct
  {
    /// @return String key for the group in which the constitutive models of this subregion are registered.
    static constexpr auto constitutiveModelsString() { return "ConstitutiveModels"; }
  };

  unsigned int particlePack( buffer_type & buffer,
                             arrayView1d< localIndex > const & localIndices,
                             bool doPack ) const;

  void particleUnpack( buffer_type & buffer,
                       int const & startingIndex,
                       int const & numberOfIncomingParticles );

  void erase(localIndex pp);

  void eraseVector(array2d< real64 > & vector, localIndex index);

  void eraseTensor(array3d< real64 > & tensor, localIndex index);

private:
  /// Group in which the constitutive models of this subregion are registered
  dataRepository::Group m_constitutiveModels;

protected:
  /// Boolean indicating whether the particle subregion contains particles needing r-vectors defining their domain extent.
  bool m_hasRVectors;

  /// Member level field for particle ghost ranks
  array1d< int > m_particleRank;

  /// Member level field for the particle global ID.
  array1d< int > m_particleID;

    /// Member level field for the particle contact group.
  array1d< int > m_particleGroup;

  /// Member level field for the particle center.
  array2d< real64 > m_particleCenter;

  /// Member level field for the particle velocity.
  array2d< real64 > m_particleVelocity;

  /// Member level field for the current particle volume.
  array1d< real64 > m_particleVolume;

  /// Member level field for the initial particle volume.
  array1d< real64 > m_particleInitialVolume;

  /// Member level field for the particle mass.
  array1d< real64 > m_particleMass;

  /// Member level field for the particle deformation gradient.
  array3d< real64 > m_particleDeformationGradient;

  /// Type of particles in this subregion.
  ParticleType m_particleType;

  /// current half-R-vectors (center to face)
  array3d< real64 > m_particleRVectors;

  /// initial half-R-vectors (center to face)
  array3d< real64 > m_particleInitialRVectors;

};


} /* namespace geosx */

#endif /* GEOSX_MESH_ELEMENTSUBREGIONBASE_HPP_ */
