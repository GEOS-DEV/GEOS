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
 * @file EmbeddedSurfaceSubRegion.hpp
 */

#ifndef EMBEDDEDSURFACESUBREGION_HPP_
#define EMBEDDEDSURFACESUBREGION_HPP_

#include "ElementSubRegionBase.hpp"
#include "InterObjectRelation.hpp"
#include "ToElementRelation.hpp"
#include "EdgeManager.hpp"
#include "CellElementSubRegion.hpp"
#include "meshUtilities/SimpleGeometricObjects/BoundedPlane.hpp"

namespace geosx
{

/**
 * @class EmbeddedSurfaceSubRegion
 *
 * The EmbeddedSurfaceSubRegion class contains the functionality to support the concept of Embedded
 *  Surface Elements.
 */
class EmbeddedSurfaceSubRegion : public ElementSubRegionBase
{
public:

  using NodeMapType = InterObjectRelation< array2d< localIndex, cells::NODE_MAP_PERMUTATION > >;
  using FaceMapType = FixedOneToManyRelation;
  using EdgeMapType = InterObjectRelation< ArrayOfArrays< localIndex > >;

  static const string CatalogName()
  { return "EmbeddedSurfaceSubRegion"; }

  virtual const string getCatalogName() const override
  {
    return EmbeddedSurfaceSubRegion::CatalogName();
  }

  EmbeddedSurfaceSubRegion( string const & name,
                            dataRepository::Group * const parent );

  virtual ~EmbeddedSurfaceSubRegion() override;

  virtual void CalculateElementGeometricQuantities( NodeManager const & nodeManager,
                                                    FaceManager const & facemanager ) override;

  void CalculateElementGeometricQuantities( localIndex const index );

  void AddNewEmbeddedSurface( localIndex const cellIndex,
                              R1Tensor normalVector );

  bool AddNewEmbeddedSurface( localIndex const cellIndex,
                              localIndex const regionIndex,
                              localIndex const subRegionIndex,
                              NodeManager const & nodeManager,
                              EdgeManager const & edgeManager,
                              FixedOneToManyRelation const & cellToEdges,
                              BoundedPlane const * plane );

  void CalculateElementGeometricQuantities( array1d< R1Tensor > const intersectionPoints,
                                            localIndex k );

  /**
   * @brief function to set the ghostRank for a list of FaceElements and set them to the value of their bounding faces.
   * @param[in] faceManager The face group.
   * @param[in] indices The list of indices to set value of ghostRank.
   */

  struct viewKeyStruct : ElementSubRegionBase::viewKeyStruct
  {
    static constexpr auto elementApertureString        = "elementAperture";
    static constexpr auto elementAreaString            = "elementArea";
    static constexpr auto cellListString               = "fractureElementsToCellIndices";
    static constexpr auto regionListString             = "fractureElementsToRegionIndex";
    static constexpr auto subregionListString          = "fractureElementsToSubRegionIndex";
    static constexpr auto normalVectorString           = "normalVector";
  };

  virtual void setupRelatedObjectsInRelations( MeshLevel const * const mesh ) override;

  virtual string GetElementTypeString() const override { return "Embedded"; }


  /**
   * @name Relation Accessors
   * @brief Accessor function for the various inter-object relations
   */
  ///@{
  NodeMapType const & nodeList() const
  {
    return m_toNodesRelation;
  }

  NodeMapType & nodeList()
  {
    return m_toNodesRelation;
  }

  ///@}

  /**
   * @return number of nodes per element
   */
  //virtual localIndex numNodesPerElement( localIndex const k ) const override { return m_toNodesRelation[k].size(); }

  arrayView1d< real64 > const & getElementAperture()       { return m_elementAperture; }
  arrayView1d< real64 const > const & getElementAperture() const { return m_elementAperture; }

  arrayView1d< real64 > const & getElementArea()       { return m_elementArea; }
  arrayView1d< real64 const > const & getElementArea() const { return m_elementArea; }

  arrayView1d< localIndex > const & getSurfaceToRegionList()       { return m_embeddedSurfaceToRegion; }
  arrayView1d< localIndex const > const & getSurfaceToRegionList() const { return m_embeddedSurfaceToRegion; }

  arrayView1d< localIndex > const & getSurfaceToSubRegionList()       { return m_embeddedSurfaceToSubRegion; }
  arrayView1d< localIndex const > const & getSurfaceToSubRegionList() const { return m_embeddedSurfaceToSubRegion; }

  arrayView1d< localIndex > const & getSurfaceToCellList()       { return m_embeddedSurfaceToCell; }
  arrayView1d< localIndex const > const & getSurfaceToCellList() const { return m_embeddedSurfaceToCell; }

  array1d< R1Tensor > & getNormalVector()       { return m_normalVector; }
  array1d< R1Tensor > const & getNormalVector() const { return m_normalVector; }

  R1Tensor & getNormalVector( localIndex k )       { return m_normalVector[k];}
  R1Tensor const & getNormalVector( localIndex k ) const { return m_normalVector[k];}

  array1d< R1Tensor > & getTangentVector1()       { return m_tangentVector1; }
  array1d< R1Tensor > const & getTangentVector1() const { return m_tangentVector1; }

  R1Tensor & getTangentVector1( localIndex k )       { return m_tangentVector1[k];}
  R1Tensor const & getTangentVector1( localIndex k ) const { return m_tangentVector1[k];}

  array1d< R1Tensor > & getTangentVector2()       { return m_tangentVector2; }
  array1d< R1Tensor > const & getTangentVector2() const { return m_tangentVector2; }

  R1Tensor & getTangentVector2( localIndex k )       { return m_tangentVector2[k];}
  R1Tensor const & getTangentVector2( localIndex k ) const { return m_tangentVector2[k];}

private:
  /// normal vector to the embedded surface element
  array1d< R1Tensor > m_normalVector;
  // tangential direction 1
  array1d< R1Tensor > m_tangentVector1;
  // tangential direction 2
  array1d< R1Tensor > m_tangentVector2;

  /// list of regions
  array1d< localIndex > m_embeddedSurfaceToRegion;
  /// list of subregions
  array1d< localIndex > m_embeddedSurfaceToSubRegion;
  /// list of elements cut by the embedded surface elem
  array1d< localIndex > m_embeddedSurfaceToCell;

  /// list of nodes
  NodeMapType m_toNodesRelation;    // Not used for now. Will need for Flow?

  /// list of edges (if necessary)
  EdgeMapType m_toEdgesRelation;    // Not used for now. Will need for Flow?

  /// The member level field for the element center
  array1d< real64 > m_elementAperture;

  /// The member level field for the element center
  array1d< real64 > m_elementArea;
};

} /* namespace geosx */

#endif /* EMBEDDEDSURFACESUBREGION_HPP_ */
