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
 * @file MeshBody.hpp
 */

#ifndef GEOSX_MESH_MESHBODY_HPP_
#define GEOSX_MESH_MESHBODY_HPP_

#include "MeshLevel.hpp"


namespace geosx
{

class MeshLevel;

/**
 * @class MeshBody
 * @brief The class is used to manage mesh body
 */
class MeshBody : public dataRepository::Group
{
public:

  /**
   * @brief Constructor for MeshBody object
   * @param [in] name the name of this instantiation of MeshBody
   * @param [in] parent the parent group of this instantiation of MeshBody
   */
  MeshBody( string const & name,
            Group * const parent );

  /**
   * @brief Create a new mesh level
   * @param [in] newLevel index of the new mesh level
   * @return reference to the created MeshLevel
   */
  MeshLevel & createMeshLevel( localIndex const newLevel );

  /**
   * @brief Get the meshLevels group
   * @return reference to the meshLevels group.
   */
  Group & getMeshLevels() { return m_meshLevels; }
  /**
   * @copydoc getMeshLevels()
   */
  Group const & getMeshLevels() const { return m_meshLevels; }

  /**
   * @brief Get mesh level
   * @param [in] level index of the mesh level
   * @return reference to MeshLevel
   */
  MeshLevel & getMeshLevel( string const & level )
  { return m_meshLevels.getGroup< MeshLevel >( level ); }

  /**
   * @brief Get mesh level
   * @param [in] level index of the mesh level
   * @return reference to const MeshLevel
   */
  MeshLevel const & getMeshLevel( string const & level ) const
  { return m_meshLevels.getGroup< MeshLevel >( level ); }

  /**
   * @brief Get mesh level
   * @param [in] level index of the mesh level
   * @return pointer to MeshLevel
   */
  MeshLevel & getMeshLevel( localIndex const level )
  { return getMeshLevel( intToMeshLevelString( level ) ); }

  /**
   * @brief Get mesh level
   * @param [in] level index of the mesh level
   * @return pointer to const MeshLevel
   */
  MeshLevel const & getMeshLevel( localIndex const level ) const
  { return getMeshLevel( intToMeshLevelString( level ) ); }

  /**
   * @brief Apply the given functor to all meshLevels on this meshBody.
   * @tparam FUNCTION the type of functor to call
   * @param[in] function  the functor to call
   */
  template< typename FUNCTION >
  void forMeshLevels( FUNCTION && function ) const
  {
    m_meshLevels.forSubGroups< MeshLevel >( std::forward< FUNCTION >( function ) );
  }

  /**
   * @copydoc forMeshLevels(FUNCTION &&) const
   */
  template< typename FUNCTION >
  void forMeshLevels( FUNCTION && function )
  {
    m_meshLevels.forSubGroups< MeshLevel >( std::forward< FUNCTION >( function ) );
  }

  /**
   * @brief Set mesh length scale used to define an absolute length tolerance
   * @param [in] scale length scale
   */
  void setGlobalLengthScale( real64 scale );

  /**
   * @brief Get mesh length scale
   * @return value of mesh length scale
   */
  real64 getGlobalLengthScale() const
  {
    return m_globalLengthScale;
  }

  /**
   * @brief Data repository keys
   */
  struct viewKeysStruct
  {} viewKeys; ///< viewKeys

  /**
   * @brief Group keys
   */
  struct groupStructKeys
  {
    /// @return The key/string used to register/access the Group that contains the MeshLevel objects.
    static constexpr char const * meshLevelsString() { return "meshLevels"; }
  } groupKeys; ///< groupKeys

private:
  Group & m_meshLevels;

  /// Mesh length scale used to define an absolute length tolerance
  /// The default value can be set to another value
  real64 m_globalLengthScale { 0. };


  static string intToMeshLevelString( localIndex const meshLevel );

};

} /* namespace geosx */

#endif /* GEOSX_MESH_MESHBODY_HPP_ */
