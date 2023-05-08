/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file MeshObjectManager.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_MULTISCALE_MULTISCALEOBJECTMANAGER_HPP
#define GEOSX_LINEARALGEBRA_MULTISCALE_MULTISCALEOBJECTMANAGER_HPP

#include "mesh/ObjectManagerBase.hpp"
#include "mesh/InterObjectRelation.hpp"

namespace geos
{
namespace multiscale
{

/**
 * @brief Mesh object manager used in multiscale preconditioners to keep
 *        a simplified (node/cell only) representation of (part of) the mesh.
 *
 * The same manager class is used for both cells and nodes. The term "dual" is
 * used to refer to the counterpart object type (e.g. cells are dual objects to
 * nodes and vice versa). Compared to geos::MeshLevel, notions of regions and
 * subregions are removed, and everything is presented in an unstructured way
 * (i.e. using ArrayOfSets for mesh maps). New contiguous global indices that
 * are consistent with dof numbering in the original system matrix are assigned
 * to mesh objects. This approach is intended to simplify further mesh processing
 * (coarsening and construction of multiscale basis support regions).
 */
class MeshObjectManager : public ObjectManagerBase
{
public:

  /// Alias for relation map type
  using MapType = InterObjectRelation< ArrayOfSets< localIndex > >;

  /// Alias for relation map const view type
  using MapViewConst = decltype( std::declval< MapType >().base().toViewConst() );

  /**
   * @brief Constructor.
   * @param name the name of this instance
   * @param parent the parent group
   */
  MeshObjectManager( string const & name, dataRepository::Group * parent );

  string getCatalogName() const final
  {
    return "MultiscaleObjectManager"; // doesn't matter really
  }

  /**
   * @brief @return number of locally owned objects
   */
  localIndex numOwnedObjects() const { return m_numOwnedObjects; }

  /**
   * @brief Set the number of locally owned objects.
   * @param n the number of locally owned objects
   */
  void setNumOwnedObjects( localIndex const n );

  /**
   * @brief @return the dual object adjacency map
   */
  MapType & toDualRelation() { return m_toDualRelation; }

  /**
   * @brief @return the dual object adjacency map (const view)
   */
  MapViewConst toDualRelation() const { return m_toDualRelation.base().toViewConst(); }

  /**
   * @brief contains the added view access keys to be bound with class data member.
   */
  struct viewKeyStruct : ObjectManagerBase::viewKeyStruct
  {
    /// @return String to access the reference position
    static constexpr char const * dualObjectString() { return "dualObject"; }
  };

private:

  localIndex m_numOwnedObjects; ///< Number of locally owned (non-ghosted) objects

  MapType m_toDualRelation; ///< Map for cell-to-node or node-to-cell relations

};

} // namespace multiscale
} // namespace geos

#endif //GEOSX_LINEARALGEBRA_MULTISCALE_MULTISCALEOBJECTMANAGER_HPP
