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
class FiniteElementDiscretization;


using DiscretizationContainer = map< std::pair< string, string >, FiniteElementDiscretization const & >;
using DiscretizationEntry = DiscretizationContainer::value_type;

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
   * @brief Destructor
   */
  virtual ~MeshBody();

  /**
   * @brief Create a new mesh level
   * @param [in] newLevel index of the new mesh level
   * @return pointer to the created MeshLevel
   */
  MeshLevel & createMeshLevel( localIndex const newLevel );

  MeshLevel & createMeshLevel( string const & name  );

  MeshLevel & createMeshLevel( string const & sourceLevelName,
                               string const & newLevelName,
                               int const order,
                               arrayView1d<string const> const & regions  );


  /**
   * @brief Get mesh level
   * @param [in] level index of the mesh level
   * @return pointer to MeshLevel
   */
  template< typename KEY >
  MeshLevel & getMeshLevel( KEY const & key )
  { return this->getGroup< MeshLevel >( key ); }

  /**
   * @brief Get mesh level
   * @param [in] level index of the mesh level
   * @return pointer to const MeshLevel
   */
  template< typename KEY >
  MeshLevel const & getMeshLevel( KEY const & key ) const
  { return this->getGroup< MeshLevel >( key ); }

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
  {
    /// The key for MeshLevel
    dataRepository::ViewKey meshLevels                = { "meshLevels" };
  } viewKeys; ///< viewKeys

  /**
   * @brief Group keys
   */
  struct groupStructKeys
  {} groupKeys; ///< groupKeys

private:
  /// Mesh length scale used to define an absolute length tolerance
  /// The default value can be set to another value
  real64 m_globalLengthScale { 0. };


};

} /* namespace geosx */

#endif /* GEOSX_MESH_MESHBODY_HPP_ */
