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

#ifndef FACECELLSUBREGION_HPP_
#define FACECELLSUBREGION_HPP_

#include "CellBase.hpp"
#include "InterObjectRelation.hpp"

namespace geosx
{

class FaceCellSubRegion : public CellBase
{
public:

  using ToNodeMap=OrderedVariableOneToManyRelation;
  using ToEdgesMap=OrderedVariableOneToManyRelation;
  using ToFacesMap=FixedOneToManyRelation;

  static const string CatalogName()
  { return "FaceCell"; }

  virtual const string getCatalogName() const override
  {
    return FaceCellSubRegion::CatalogName();
  }

  FaceCellSubRegion( string const & name,
                     dataRepository::ManagedGroup * const parent );
  virtual ~FaceCellSubRegion() override;

  virtual R1Tensor const & calculateElementCenter( localIndex k,
                                     const NodeManager& nodeManager,
                                     const bool useReferencePos = true) const override;

  virtual void CalculateCellVolumes( array1d<localIndex> const & indices,
                                     array1d<R1Tensor> const & X ) override
  {
    CellBase::CalculateCellVolumes<FaceCellSubRegion>( *this,
                                                       indices,
                                                       X );
  }

  inline void CalculateCellVolumesKernel( localIndex const k,
                                          array1d<R1Tensor> const & X )
  {
    m_elementArea[k] = 1;
    m_elementCenter[k] = 1;
    m_elementVolume[k] = 1;
  }

  struct viewKeyStruct : ObjectManagerBase::viewKeyStruct
  {
    static constexpr auto nodeListString               = "nodeList";
    static constexpr auto edgeListString               = "edgeList";
    static constexpr auto faceListString               = "faceList";
    static constexpr auto elementApertureString        = "elementAperture";
    static constexpr auto elementAreaString            = "elementArea";
    static constexpr auto elementCenterString          = "elementCenter";
    static constexpr auto elementVolumeString          = "elementVolume";
  };

  virtual void setupRelatedObjectsInRelations( MeshLevel const * const mesh ) override;


  OrderedVariableOneToManyRelation const & nodeList() const
  {
    return m_toNodesRelation;
  }

  OrderedVariableOneToManyRelation & nodeList()
  {
    return m_toNodesRelation;
  }

  virtual arraySlice1d<localIndex const> nodeList( localIndex const k ) const override
  {
    return m_toNodesRelation[k];
  }

  virtual arraySlice1d<localIndex> nodeList( localIndex const k ) override
  {
    return m_toNodesRelation[k];
  }

private:

  /// The elements to nodes relation
  OrderedVariableOneToManyRelation  m_toNodesRelation;

  /// The elements to edges relation
  OrderedVariableOneToManyRelation  m_toEdgesRelation;

  /// The elements to faces relation
  FixedOneToManyRelation  m_toFacesRelation;

  /// The member level field for the element center
  array1d< R1Tensor > m_elementAperture;

  /// The member level field for the element center
  array1d< R1Tensor > m_elementArea;
};

} /* namespace geosx */

#endif /* FACECELLSUBREGION_HPP_ */
