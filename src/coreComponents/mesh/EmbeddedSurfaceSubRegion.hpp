/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

    using NodeMapType = InterObjectRelation<array1d<array1d<localIndex>>>;
    using EdgeMapType = InterObjectRelation<array1d<array1d<localIndex>>>;
    using FaceMapType = InterObjectRelation<array2d<localIndex>>;

    static const string CatalogName()
    { return "EmbeddedSurfaceSubRegion"; }

    virtual const string getCatalogName() const override
    {
      return EmbeddedSurfaceSubRegion::CatalogName();
    }

    EmbeddedSurfaceSubRegion( string const & name,
                       dataRepository::Group * const parent );

    virtual ~EmbeddedSurfaceSubRegion() override;

    virtual R1Tensor const & calculateElementCenter(localIndex k,
                                                    const NodeManager& GEOSX_UNUSED_ARG( nodeManager ),
                                                    const bool GEOSX_UNUSED_ARG( useReferencePos ) = true) const override
    {
      return m_elementCenter[k];
    }

    virtual void CalculateElementGeometricQuantities( NodeManager const & nodeManager,
                                                      FaceManager const & facemanager ) override;

    void CalculateElementGeometricQuantities( localIndex const index);

    void AddNewEmbeddedSurface(localIndex const cellIndex,
                               R1Tensor normalVector);


    void CalculateElementGeometricQuantities(NodeManager const & nodeManager,
                                             EdgeManager const & edgeManager,
                                             FixedOneToManyRelation const & cellToEdges,
                                             R1Tensor origin);

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
      static constexpr auto normalVectorString           = "normalVector";

      //
      //static constexpr auto faceElementsToCellRegionsString    = "fractureElementsToCellRegions";
      //static constexpr auto faceElementsToCellSubRegionsString    = "fractureElementsToCellSubRegions";
      //static constexpr auto faceElementsToCellIndexString    = "fractureElementsToCellIndices";
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

    arrayView1d< real64 > const &       getElementAperture()       { return m_elementAperture; }
    arrayView1d< real64 const > const & getElementAperture() const { return m_elementAperture; }

    arrayView1d< real64 > const &       getElementArea()       { return m_elementArea; }
    arrayView1d< real64 const > const & getElementArea() const { return m_elementArea; }

    arrayView1d< localIndex > const &       getSurfaceToCellList()       { return m_embeddedSurfaceToCell; }
    arrayView1d< localIndex const > const & getSurfaceToCellList() const { return m_embeddedSurfaceToCell; }

    array1d< R1Tensor >  &       getNormalVector()       { return m_normalVector; }
    array1d< R1Tensor > const &  getNormalVector() const { return m_normalVector; }

private:
    /// normal vector to the embedded surface element
    array1d < R1Tensor > m_normalVector;

    /// list of elements cut by the embedded surface el
    array1d< localIndex > m_embeddedSurfaceToCell;

    /// list of nodes
    NodeMapType  m_toNodesRelation; // ?

    /// list of edges (if necessary)
    EdgeMapType  m_toEdgesRelation; // ?

    /// The member level field for the element center
    array1d< real64 > m_elementAperture;

    /// The member level field for the element center
    array1d< real64 > m_elementArea;
};

} /* namespace geosx */

#endif /* EMBEDDEDSURFACESUBREGION_HPP_ */
