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
  static string catalogName()
  { return "ParticleSubRegion"; }

  /**
   * @copydoc catalogName()
   */
  virtual string getCatalogName() const override final
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
   * @brief Update the local-to-global and global-to-local maps
   */
  void updateMaps();

  /**
   * @name Helpers for ParticleSubRegion construction
   */
  ///@{

  /**
   * @brief Set the ghost rank of particles in this subregion
   * @param rank the mpi rank to which all particles in this subregion (in this partition) will have their ghost rank set to
   */
  void setParticleRank( int rank, int np );

  /**
   * @brief Get the maximum number of nodes a particle could map to
   */
  int numNodesMappedTo();

  /**
   * @brief Fill the ParticleSubRegion by copying those of the source ParticleBlock
   * @param particleBlock the ParticleBlock which properties (connectivity info) will be copied.
   */
  void copyFromParticleBlock( ParticleBlockABC & particleBlock );

  ///@}

  void flagOutOfRangeParticles( std::array< real64, 3 > const & xGlobalMin,
                                std::array< real64, 3 > const & xGlobalMax,
                                std::array< real64, 3 > const & hEl,
                                arrayView1d< int > const isBad );

  /**
   * @brief This function modifies m_particleRVectors using the particle deformation gradient. Used by solvers.
   * @param p the local particle index
   * @param F the particle deformation gradient
   * @param intialRVectors the particle's initial R-Vectors
   */
  void computeRVectors( int const p,
                        arraySlice2d< real64 > const F,
                        arraySlice2d< real64 > const initialRVectors );

  void cpdiDomainScaling( real64 lCrit,
                          int m_planeStrain );

  void getMappedNodes( int const p,
                       std::array< real64, 3 > const & xMin,
                       std::array< real64, 3 > const & hx,
                       array3d< int > const & ijkMap,
                       arrayView1d< localIndex > const nodeIDs );

  void getAllWeights( int const p,
                      std::array< real64, 3 > const & xMin,
                      std::array< real64, 3 > const & hx,
                      array3d< int > const & ijkMap,
                      arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const g_X,
                      arrayView1d< int > const nodeIDs,
                      arrayView1d< real64 > const weights,
                      arrayView2d< real64 > const gradWeights );

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

private:

  /// Map used for constitutive grouping
  map< string, localIndex_array > m_constitutiveGrouping;

  /// Name of the properties registered from an external mesh
  string_array m_externalPropertyNames;

};

} /* namespace geosx */

#endif /* GEOSX_MESH_CELLELEMENTSUBREGION_HPP_ */
