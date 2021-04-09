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

/**
 * @file ParticleElementRegion.hpp
 *
 */

#ifndef GEOSX_MESH_CELLELEMENTREGION_HPP_
#define GEOSX_MESH_CELLELEMENTREGION_HPP_

#include "ElementRegionBase.hpp"

namespace geosx
{
class EdgeManager;
class EmbeddedSurfaceGenerator;

/**
 * @class ParticleElementRegion
 *
 * The ParticleElementRegion class contains the functionality to support the concept of a ParticleElementRegion in the element
 * hierarchy. ParticleElementRegion derives from ElementRegionBase and has an entry in the ObjectManagerBase catalog.
 *
 *
 */
class ParticleElementRegion : public ElementRegionBase
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
  ParticleElementRegion( string const & name, Group * const parent );

  /**
   * @brief Deleted default constructor.
   */
  ParticleElementRegion() = delete;

  /**
   * @brief Destructor.
   */
  virtual ~ParticleElementRegion() override;

  /**
   * @name Static factory catalog functions
   */
  ///@{

  /**
   * @brief The key name for the FaceElementRegion in the object catalog.
   * @return A string containing the key name.
   */
  static const string catalogName()
  { return "ParticleElementRegion"; }

  /**
   * @copydoc catalogName()
   */
  virtual const string getCatalogName() const override final
  { return ParticleElementRegion::catalogName(); }

  ///@}

  /**
   * @name Generation of the particle element mesh region
   */
  ///@{

  /**
   * @brief Add a particleBlockRegion name to the list.
   * @param particleBlockName string containing the particle block region name.
   */
  void addParticleBlockName( string const & particleBlockName )
  {
    m_particleBlockNames.emplace_back( particleBlockName );
  }

  virtual void generateMesh( Group & particleBlocks ) override;

  /**
   * @brief Generate the aggregates.
   * @param faceManager a pointer to the FaceManager
   * @param nodeManager a pointer to the NodeManager
   */
  void generateAggregates( FaceManager const & faceManager, NodeManager const & nodeManager );

  ///@}

  /**
   * @brief A struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : public ElementRegionBase::viewKeyStruct
  {
    /// @return String key for the coarsening ratio
    static constexpr char const * coarseningRatioString() {return "coarseningRatio"; }

    /// @return String key for the particle block names
    static constexpr char const * sourceParticleBlockNamesString() {return "particleBlocks"; }
  };


private:

  // Particle block names
  string_array m_particleBlockNames;

  // Coarsening ratio
  real64 m_coarseningRatio;

};

} /* namespace geosx */

#endif /* GEOSX_MESH_CELLELEMENTREGION_HPP_ */
