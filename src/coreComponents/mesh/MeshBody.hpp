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
   * @brief Create a new mesh level.
   * @param[in] name The name of the new MeshLevel.
   * @return A reference to the new MeshLevel.
   */
  MeshLevel & createMeshLevel( string const & name );

  /**
   * @brief Create a new mesh level from a source MeshLevel.
   * @param sourceLevelName The name of the source MeshLevel.
   * @param newLevelName The name of the new MeshLevel.
   * @param order The order of the new MeshLevel.
   * @return A reference to the new MeshLevel.
   */
  MeshLevel & createMeshLevel( string const & sourceLevelName,
                               string const & newLevelName,
			       CellBlockManager & cellBlockManager,
                               int const order );

  /**
   * @brief Creates a mesh level in which the member pointers are set to the
   *        allocations from another MeshLevel.
   * @param sourceLevelName The MeshLevel to be "copied"
   * @param newLevelName The name of the new shallow copy.
   * @return A reference to the new shallow MeshLevel.
   */
  MeshLevel & createShallowMeshLevel( string const & sourceLevelName,
                                      string const & newLevelName );

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
   * @brief Get a reference to a MeshLevel.
   * @tparam T The type of the lookup key. Function is only implemented for
   *  string and const char * key types.
   * @param[in] level The lookup key of the MeshLevel
   * @return const reference to the MeshLevel
   */
  template< typename T, std::enable_if_t< std::is_same< T, string >::value ||
                                          std::is_same< T, const char * >::value, bool > = false >
  MeshLevel & getMeshLevel( T const & level ) const
  { return m_meshLevels.getGroup< MeshLevel >( level ); }

  /**
   * @brief Get a reference to a MeshLevel.
   * @tparam T The type of the lookup key. Function is only implemented for
   *  string and const char * key types.
   * @param[in] level The lookup key of the MeshLevel
   * @return Reference to the MeshLevel
   */
  template< typename T, std::enable_if_t< std::is_same< T, string >::value ||
                                          std::is_same< T, const char * >::value, bool > = false >
  MeshLevel & getMeshLevel( T const & level )
  { return m_meshLevels.getGroup< MeshLevel >( level ); }


  /**
   * @brief Convenience function to access the baseDiscretization.
   * @return Reference to the base discretization.
   */
  MeshLevel & getBaseDiscretization()
  {
    return getMeshLevel( groupStructKeys::baseDiscretizationString() );
  }

  /**
   * @brief Convenience function to access the baseDiscretization.
   * @return Reference to the base discretization.
   */
  MeshLevel const & getBaseDiscretization() const
  {
    return getMeshLevel( groupStructKeys::baseDiscretizationString() );
  }

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

    /// @return The key/string used to register/access the Group that contains the base discretization.
    static constexpr char const * baseDiscretizationString() { return "Level0"; }
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
