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

#include "dataRepository/Group.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "mesh/generators/CellBlockABC.hpp"

namespace geosx
{

/**
 * @class CellBlock
 * Class deriving from ElementSubRegionBase specializing the element subregion
 * for a cell element. In particular, this class adds connectivity maps (namely,
 * element-to-node map, element-to-face map, and element-to-edge map) and
 * and methods to compute the geometry of an element (cell center and volume)
 */
class CellBlock : public CellBlockABC
{
public:

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

  string getElementTypeString() const override
  { return m_elementTypeString; }

  localIndex numNodesPerElement() const override { return m_numNodesPerElement; }

  localIndex numFacesPerElement() const override { return m_numFacesPerElement; }

  localIndex numElements() const override { return size(); }

  std::vector< localIndex > getFaceNodes( localIndex iElement,
                                          localIndex iFace ) const override;

  /**
   * @copydoc nodeList()
   */
  NodeMapType & getElemToNode() { return m_toNodesRelation; }

  NodeMapType const & getElemToNode() const override { return m_toNodesRelation; }

  array2d< localIndex > const & getElemToFaces() const override { return m_toFacesRelation; }

  void setElementToFaces( localIndex iFace, localIndex j, localIndex curFaceID );

  /**
   * @brief Get local to global map, non-const version.
   * @return The mapping relationship as a array.
   */
  arrayView1d< globalIndex > localToGlobalMap()
  { return m_localToGlobalMap; }

  arrayView1d< globalIndex const > localToGlobalMap() const override
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

  /// Element-to-node relation
  array2d< localIndex > m_toFacesRelation;

  /// Contains the global index of each object.
  array1d< globalIndex > m_localToGlobalMap;

  /// Name of the properties registered from an external mesh
  string_array m_externalPropertyNames;

  /// Type of element in this subregion.
  string m_elementTypeString;

  std::list< dataRepository::WrapperBase * > getExternalProperties() override
  {
    std::list< dataRepository::WrapperBase * > result;
    for( auto & externalPropertyName : m_externalPropertyNames )
    {
      result.push_back( &this->getWrapperBase( externalPropertyName ) );
    }
    return result;
  }
};

}

#endif /* GEOSX_MESH_CELLBLOCK_HPP_ */
