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
 * @file ParticleBlock.hpp
 */

#ifndef GEOSX_MESH_PARTICLEBLOCK_HPP_
#define GEOSX_MESH_PARTICLEBLOCK_HPP_

#include "ElementSubRegionBase.hpp"
#include "FaceManager.hpp"
#include "utilities/ComputationalGeometry.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geosx
{

/**
 * @class ParticleBlock
 * Class deriving from ElementSubRegionBase specializing the particle subregion
 * for a particle. In particular, this class adds connectivity maps (namely,
 * element-to-node map, element-to-face map, and element-to-edge map) and
 * and methods to compute the geometry of an element (cell center and volume)
 */
class ParticleBlock : public ElementSubRegionBase
{
public:

  /// Alias for the type of the element-to-node map
  using NodeMapType = InterObjectRelation< array2d< localIndex, cells::NODE_MAP_PERMUTATION > >;
  /// Alias for the type of the element-to-edge map
  using EdgeMapType = FixedOneToManyRelation;
  /// Alias for the type of the element-to-face map
  using FaceMapType = FixedOneToManyRelation;

  /**
   * @name Static factory catalog functions
   */
  ///@{

  /**
   * @brief Const getter for the catalog name.
   * @return the name of this type in the catalog
   */
  static const string catalogName()
  { return "ParticleBlock"; }

  /**
   * @copydoc catalogName()
   */
  virtual const string getCatalogName() const override final
  { return ParticleBlock::catalogName(); }

  ///@}

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Deleted default constructor.
   */
  ParticleBlock() = delete;

  /**
   * @brief Constructor for this class.
   * @param[in] name the name of this object manager
   * @param[in] parent the parent Group
   */
  ParticleBlock( string const & name, Group * const parent );

  /**
   * @brief Copy constructor.
   * @param[in] init the source to copy
   */
  ParticleBlock( const ParticleBlock & init ) = delete;

  /**
   * @brief Destructor.
   */
  virtual ~ParticleBlock() override;

  ///@}

  /**
   * @name Geometry computation / Connectivity
   */
  ///@{

  virtual void calculateElementGeometricQuantities( NodeManager const & nodeManager,
                                                    FaceManager const & faceManager ) override;

  virtual void setupRelatedObjectsInRelations( MeshLevel const & mesh ) override;

  ///@}
  /**
   * @name Getters / Setters
   */
  ///@{

  virtual void setElementType( ElementType const elementType ) override;

  /**
   * @brief Get the number of the nodes in a face of the element.
   * @param elementIndex the local index of the target element
   * @param localFaceIndex the local index of the target face in the element  (this will be 0-numFacesInElement)
   * @return the number of nodes of this face
   */
  localIndex getNumFaceNodes( localIndex const elementIndex,
                              localIndex const localFaceIndex ) const;

  /**
   * @brief Get the local indices of the nodes in a face of the element.
   * @param elementIndex the local index of the target element
   * @param localFaceIndex the local index of the target face in the element  (this will be 0-numFacesInElement)
   * @param nodeIndices a pointer to the node indices of the face
   * @return the number of nodes in the face
   */
  localIndex getFaceNodes( localIndex const elementIndex,
                           localIndex const localFaceIndex,
                           localIndex * const nodeIndices ) const;


  /**
   * @brief Get the local indices of the nodes in a face of the element.
   * @param elementIndex the local index of the target element
   * @param localFaceIndex the local index of the target face in the element  (this will be 0-numFacesInElement)
   * @param nodeIndices a reference to the array of node indices of the face
   */
  void getFaceNodes( localIndex const elementIndex,
                     localIndex const localFaceIndex,
                     localIndex_array & nodeIndices ) const;

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
   * @brief Get the local index of the a-th node of the k-th element.
   * @param[in] k the index of the element
   * @param[in] a the index of the node in the element
   * @return a reference to the local index of the node
   */
  localIndex & nodeList( localIndex const k, localIndex a ) { return m_toNodesRelation( k, a ); }

  /**
   * @copydoc nodeList( localIndex const k, localIndex a )
   */
  localIndex const & nodeList( localIndex const k, localIndex a ) const { return m_toNodesRelation( k, a ); }

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

  ///@}

  /**
   * @name Properties
   */
  ///@{

  /**
   * @brief Add a property to the ParticleBlock.
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

protected:

  /// Element-to-node relation
  NodeMapType m_toNodesRelation;

  /// Element-to-edge relation
  EdgeMapType m_toEdgesRelation;

  /// Element-to-face relation
  FaceMapType m_toFacesRelation;

private:
  /// Name of the properties registered from an external mesh
  string_array m_externalPropertyNames;

  /**
   * @brief Compute the volume of the k-th element in the subregion.
   * @param[in] k the index of the element in the subregion
   * @param[in] X an arrayView of (const) node positions
   */
  inline void calculateCellVolumesKernel( localIndex const k,
                                          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X ) const
  {
    LvArray::tensorOps::fill< 3 >( m_elementCenter[ k ], 0 );

    real64 Xlocal[10][3];

    for( localIndex a = 0; a < m_numNodesPerElement; ++a )
    {
      LvArray::tensorOps::copy< 3 >( Xlocal[ a ], X[ m_toNodesRelation( k, a ) ] );
      LvArray::tensorOps::add< 3 >( m_elementCenter[ k ], X[ m_toNodesRelation( k, a ) ] );
    }
    LvArray::tensorOps::scale< 3 >( m_elementCenter[ k ], 1.0 / m_numNodesPerElement );

    switch( m_elementType )
    {
      case ElementType::Hexahedron:
      {
        m_elementVolume[k] = computationalGeometry::hexVolume( Xlocal );
        break;
      }
      case ElementType::Tetrahedron:
      {
        m_elementVolume[k] = computationalGeometry::tetVolume( Xlocal );
        break;
      }
      case ElementType::Prism:
      {
        m_elementVolume[k] = computationalGeometry::wedgeVolume( Xlocal );
        break;
      }
      case ElementType::Pyramid:
      {
        m_elementVolume[k] = computationalGeometry::pyramidVolume( Xlocal );
        break;
      }
      default:
      {
        GEOSX_ERROR( "Volume calculation not supported for element type: " << m_elementType );
      }
    }
  }
};

}

#endif /* GEOSX_MESH_PARTICLEBLOCK_HPP_ */
