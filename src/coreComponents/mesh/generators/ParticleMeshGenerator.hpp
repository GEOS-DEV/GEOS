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
 * @file ParticleMeshGenerator.hpp
 */

#ifndef GEOSX_MESH_GENERATORS_PARTICLEMESHGENERATOR_HPP
#define GEOSX_MESH_GENERATORS_PARTICLEMESHGENERATOR_HPP

#include "codingUtilities/EnumStrings.hpp"
#include "mesh/generators/MeshGeneratorBase.hpp"

namespace geosx
{

class NodeManager; // Probably don't need this?
class ParticleManager; // This needs to be created
class SpatialPartition; // Need this to initially partition particles

/**
 * @class ParticleMeshGenerator
 * @brief The ParticleMeshGenerator class is a class handling import of particle data from an externel particle file.
 */
class ParticleMeshGenerator : public MeshGeneratorBase
{
public:

  /**
   * @brief Main constructor for ParticleMeshGenerator.
   * @param[in] name of the ParticleMeshGenerator
   * @param[in] parent point to the parent Group of the ParticleMeshGenerator
   */
  ParticleMeshGenerator( const string & name, Group * const parent );

  virtual ~ParticleMeshGenerator() override = default;

  /**
   * @brief Return the name of the ParticleMeshGenerator in object Catalog.
   * @return string that contains the key name to ParticleMeshGenerator in the Catalog
   */
  static string catalogName() { return "ParticleMesh"; }

  /**
   * @brief Create a new geometric object (box, plane, etc) as a child of this group.
   * @param childKey the catalog key of the new geometric object to create
   * @param childName the name of the new geometric object in the repository
   * @return the group child
   */
  virtual Group * createChild( string const & childKey, string const & childName ) override;

  virtual void generateMesh( DomainPartition & domain ) override;

  /**
   * @return Whether or not a Cartesian mesh is being generated.
   */
  virtual inline bool isCartesian() const
  {
    return true;
  }


protected:

  ///@cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    constexpr static char const * filePathString() { return "file"; }
    constexpr static char const * particleBlockNamesString() { return "particleBlockNames"; }
    constexpr static char const * positionToleranceString() { return "positionTolerance"; }
  };
  /// @endcond

  void postProcessInput() override;

  /// Mesh number of dimension
  int m_dim;

  /// Minimum extent of particle coordinates
  real64 m_min[3];

  /// Maximum extent of particle coordinates
  real64 m_max[3];

  /// Array of particle coordinates
  array1d< real64 > m_pCoord[3];

private:

  /// Path to the particle file
  Path m_filePath;

  /// String array of particle region names
  array1d< string > m_regionNames;

public:

};

} /* namespace geosx */

#endif /* GEOSX_MESH_GENERATORS_PARTICLEMESHGENERATOR_HPP */
