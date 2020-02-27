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

#ifndef GEOSX_MESH_AGGREGATEELEMENTSUBREGION_HPP_
#define GEOSX_MESH_AGGREGATEELEMENTSUBREGION_HPP_

#include "ElementSubRegionBase.hpp"
#include "InterObjectRelation.hpp"



namespace geosx
{

class AggregateElementSubRegion : public ElementSubRegionBase
{
public:

  using NodeMapType=FixedOneToManyRelation;

  static const string CatalogName()
  { return "AggregateCell"; }

  virtual const string getCatalogName() const override
  {
    return AggregateElementSubRegion::CatalogName();
  }

  template< typename LAMBDA >
  void forFineCellsInAggregate( localIndex aggregateIndex, LAMBDA lambda )
  {
    for(localIndex fineCell = m_nbFineCellsPerCoarseCell[aggregateIndex]; 
        fineCell < m_nbFineCellsPerCoarseCell[aggregateIndex+1]; fineCell++)
    {
      lambda(m_fineToCoarse[fineCell]);
    }
  }

  localIndex GetNbCellsPerAggregate( localIndex aggregateIndex ) const
  {
    return m_nbFineCellsPerCoarseCell[aggregateIndex + 1] - m_nbFineCellsPerCoarseCell[aggregateIndex];
  }

  AggregateElementSubRegion( string const & name,
                             dataRepository::Group * const parent );

  virtual ~AggregateElementSubRegion() override;
 
  void CreateFromFineToCoarseMap( localIndex nbAggregates,
                                  array1d< localIndex > const & fineToCoarse,
                                  array1d< R1Tensor > const & barycenters);

  const array1d< localIndex >& GetFineToCoarseMap()
  {
    return m_fineToCoarse;
  }

  virtual void CalculateElementGeometricQuantities( NodeManager const & GEOSX_UNUSED_PARAM( nodeManager ),
                                                    FaceManager const & GEOSX_UNUSED_PARAM( faceManager ) ) override
  {
      //TODO ?
  }

  virtual void setupRelatedObjectsInRelations( MeshLevel const * const GEOSX_UNUSED_PARAM( mesh ) ) override
  {
    //TODO ?
  }

  struct viewKeyStruct : ObjectManagerBase::viewKeyStruct
  {
    static constexpr auto elementVolumeString          = "elementVolume";
    static constexpr auto fineElementsListString       = "fineElements";
  };

private:
  /// The elements to nodes relation is one to one relation.
  NodeMapType  m_toNodesRelation;

  /// Relation between fine and coarse elements ordered by aggregates
  array1d< localIndex > m_fineToCoarse;

  /// Number of fine cells per aggregate
  array1d< localIndex > m_nbFineCellsPerCoarseCell;
};
}

#endif
