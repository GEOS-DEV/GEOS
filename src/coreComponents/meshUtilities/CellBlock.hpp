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
 * @file CellBlock.hpp
 */

#ifndef GEOSX_MESH_CELLBLOCK_HPP_
#define GEOSX_MESH_CELLBLOCK_HPP_

#include "mesh/FaceManager.hpp"
#include "dataRepository/Group.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

namespace geosx
{

/**
 * @class CellBlock
 * Class deriving from ElementSubRegionBase specializing the element subregion
 * for a cell element. In particular, this class adds connectivity maps (namely,
 * element-to-node map, element-to-face map, and element-to-edge map) and
 * and methods to compute the geometry of an element (cell center and volume)
 */
class CellBlock : public dataRepository::Group
{
public:

  /// Alias for the type of the element-to-node map
  using NodeMapType = InterObjectRelation< array2d< localIndex, cells::NODE_MAP_PERMUTATION > >;
  /// Alias for the type of the element-to-edge map
  using EdgeMapType = FixedOneToManyRelation;
  /// Alias for the type of the element-to-face map
  using FaceMapType = FixedOneToManyRelation;

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Deleted default constructor.
   */
  CellBlock() = delete;

  /**
   * @brief Constructor for this class.
   * @param[in] name the name of this object manager
   * @param[in] parent the parent Group
   */
  CellBlock( string const & name, Group * const parent );

  /**
   * @brief Copy constructor.
   * @param[in] init the source to copy
   */
  CellBlock( const CellBlock & init ) = delete;

  /**
   * @brief Destructor.
   */
  ~CellBlock() override = default;

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

  void setElementType( string const & elementType );

  /**
   * @brief Get the type of element in this subregion.
   * @return a string specifying the type of element in this subregion
   *
   * See class FiniteElementBase for possible element type.
   */
  string getElementTypeString() const
  { return m_elementTypeString; }

  /**
   * @brief Get the number of nodes per element.
   * @return number of nodes per element
   */
  localIndex const & numNodesPerElement() const { return m_numNodesPerElement; }

  /**
   * @brief Get the number of faces per element.
   * @return number of faces per element
   */
  localIndex const & numFacesPerElement() const { return m_numFacesPerElement; }

  /**
   * @brief Get the element-to-node map.
   * @return a reference to the element-to-node map
   */
  NodeMapType & nodeList() { return m_toNodesRelation; }

  /**
   * @copydoc nodeList()
   */
  NodeMapType const & nodeList() const { return m_toNodesRelation; }

  /**
   * @brief Get the element-to-edge map.
   * @return a reference to element-to-edge map
   */
  FixedOneToManyRelation & edgeList() { return m_toEdgesRelation; }

  /**
   * @copydoc edgeList()
   */
  FixedOneToManyRelation const & edgeList() const { return m_toEdgesRelation; }

  /**
   * @brief Get the element-to-face map.
   * @return a reference to the element to face map
   */
  FixedOneToManyRelation & faceList() { return m_toFacesRelation; }

  /**
   * @copydoc faceList()
   */
  FixedOneToManyRelation const & faceList() const { return m_toFacesRelation; }

  /**
   * @brief Get local to global map.
   * @return The mapping relationship as a array.
   */
  arrayView1d< globalIndex > localToGlobalMap()
  { return m_localToGlobalMap; }

  /**
   * @brief Get local to global map, const version.
   * @return The mapping relationship as a array.
   */
  arrayView1d< globalIndex const > localToGlobalMap() const
  { return m_localToGlobalMap; }

  void resize( dataRepository::indexType const newSize ) final;

  ///@}

  /**
   * @name Properties
   */
  ///@{

  /**
   * @brief Add a property to the CellBlock.
   * @tparam T type of the property
   * @param[in] propertyName the name of the property
   * @return a non-const reference to the property
   */
  template< typename T >
  T & addProperty( string const & propertyName )
  {
    m_externalPropertyNames.emplace_back( propertyName );
    return this->registerWrapper< T >( propertyName ).reference();
  }

  /**
   * @brief Helper function to apply a lambda function over all the external properties of the subregion
   * @tparam LAMBDA the type of the lambda function
   * @param lambda lambda function that is applied to the wrappers of external properties
   */
  template< typename LAMBDA >
  void forExternalProperties( LAMBDA && lambda )
  {
    for( auto & externalPropertyName : m_externalPropertyNames )
    {
      lambda( this->getWrapperBase( externalPropertyName ) );
    }
  }

  ///@}

private:

  /// Number of nodes per element in this subregion.
  localIndex m_numNodesPerElement;

  /// Number of edges per element in this subregion.
  localIndex m_numEdgesPerElement;

  /// Number of faces per element in this subregion.
  localIndex m_numFacesPerElement;

  /// Element-to-node relation
  NodeMapType m_toNodesRelation;

  /// Element-to-edge relation
  EdgeMapType m_toEdgesRelation;

  /// Element-to-face relation
  FaceMapType m_toFacesRelation;

  /// Contains the global index of each object.
  array1d< globalIndex > m_localToGlobalMap;

  /// Name of the properties registered from an external mesh
  string_array m_externalPropertyNames;

  /// Type of element in this subregion.
  string m_elementTypeString;
};

}

#endif /* GEOSX_MESH_CELLBLOCK_HPP_ */
