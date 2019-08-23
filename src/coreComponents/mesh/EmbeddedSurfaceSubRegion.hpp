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

#ifndef FACECELLSUBREGION_HPP_
#define FACECELLSUBREGION_HPP_

#include "ElementSubRegionBase.hpp"
#include "InterObjectRelation.hpp"
#include "ToElementRelation.hpp"

namespace geosx
{

/**
 * @class EmbeddedSurfaceSubRegion
 *
 * The EmbeddedSurfaceSubRegion class contains the functionality to support the concept of Elements that are comprised of
 * a face, or a pair of faces. This class which derives from ElementSubRegionBase, has specific connectivity maps and
 * and methods to support the specific geometry of an element comprised of a reduced dimensionality face element (i.e.
 * face area and aperture = volume)
 */
class EmbeddedSurfaceSubRegion : public ElementSubRegionBase
{
public:

  using NodeMapType=OrderedVariableOneToManyRelation;
  using EdgeMapType=OrderedVariableOneToManyRelation;
  using FaceMapType=FixedOneToManyRelation;

  static const string CatalogName()
  { return "EmbeddedSurfaceSubRegion"; }

  virtual const string getCatalogName() const override
  {
    return EmbeddedSurfaceSubRegion::CatalogName();
  }

  EmbeddedSurfaceSubRegion( string const & name,
                     dataRepository::ManagedGroup * const parent );
  virtual ~EmbeddedSurfaceSubRegion() override;

  virtual R1Tensor const & calculateElementCenter( localIndex k,
                                     const NodeManager& nodeManager,
                                     const bool useReferencePos = true) const override;

  virtual void CalculateElementGeometricQuantities( NodeManager const & nodeManager,
                                                    FaceManager const & facemanager ) override;

  void CalculateElementGeometricQuantities( localIndex const index,
                                            arrayView1d<real64 const> const & faceArea );

  virtual localIndex PackUpDownMapsSize( arrayView1d<localIndex const> const & packList ) const override;
  virtual localIndex PackUpDownMaps( buffer_unit_type * & buffer,
                                     arrayView1d<localIndex const> const & packList ) const override;

  virtual localIndex UnpackUpDownMaps( buffer_unit_type const * & buffer,
                                       localIndex_array & packList,
                                       bool const overwriteUpMaps,
                                       bool const overwriteDownMaps ) override;

  virtual void FixUpDownMaps( bool const clearIfUnmapped ) override;

  virtual void ViewPackingExclusionList( set<localIndex> & exclusionList ) const override;

  /**
   * @brief function to set the ghostRank for a list of FaceElements and set them to the value of their bounding faces.
   * @param[in] faceManager The face group.
   * @param[in] indices The list of indices to set value of ghostRank.
   */
  void inheritGhostRankFromParentFace( FaceManager const * const faceManager,
                                       std::set<localIndex> const & indices );

  // Matteo: need to check what this is for
  struct viewKeyStruct : ElementSubRegionBase::viewKeyStruct
  {
    static constexpr auto elementApertureString        = "elementAperture";
    static constexpr auto elementAreaString            = "elementArea";
  };

  virtual void setupRelatedObjectsInRelations( MeshLevel const * const mesh ) override;

  virtual string GetElementTypeString() const override { return "C3D8"; }

  /**
   * @name Relation Accessors
   * @brief Accessor function for the various inter-object relations
   */

  arrayView1d< real64 > const &       getElementAperture()       { return m_elementAperture; }
  arrayView1d< real64 const > const & getElementAperture() const { return m_elementAperture; }

  arrayView1d< real64 > const &       getElementArea()       { return m_elementArea; }
  arrayView1d< real64 const > const & getElementArea() const { return m_elementArea; }

  arrayView1d< real64 > const &       getSurfaceToCellList()       { return m_embeddedSurfaceToCell; }
  arrayView1d< real64 const > const & getSurfaceToCellList() const { return m_embeddedSurfaceToCell; }

  // FixedToManyElementRelation m_faceElementsToCells; Matteo: what is this ??

  set< localIndex > m_newFaceElements; // Matteo: I am assuming this is for propagation?

private:
  template<bool DOPACK>
  localIndex PackUpDownMapsPrivate( buffer_unit_type * & buffer,
                                    arrayView1d<localIndex const> const & packList ) const;

  /// The member connection to cell elements
  array1d< globalIndex > m_embeddedSurfaceToCell;

  /// The member level field for the element center
  array1d< real64 > m_elementAperture;

  /// The member level field for the element center
  array1d< real64 > m_elementArea;

};

} /* namespace geosx */

#endif /* FACECELLSUBREGION_HPP_ */
