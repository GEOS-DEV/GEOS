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

#ifndef AGGREGATECELLSUBREGION_HPP_
#define AGGREGATECELLSUBREGION_HPP_

#include "ElementSubRegionBase.hpp"
#include "InterObjectRelation.hpp"

#include "metis.h"


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
                          dataRepository::ManagedGroup * const parent );

  virtual ~AggregateElementSubRegion() override;
 
  void CreateFromFineToCoarseMap( localIndex nbAggregates,
                                  const array1d< idx_t >& fineToCoarse,
                                  const array1d< R1Tensor >& barycenters);

  const array1d< localIndex >& GetFineToCoarseMap()
  {
    return m_fineToCoarse;
  }
  
  virtual R1Tensor const & calculateElementCenter( localIndex k,
                                     const NodeManager& nodeManager,
                                     const bool useReferencePos = true) const override
  {
    return m_elementCenter[k];
  }

  virtual void CalculateCellVolumes( array1d<localIndex> const & indices,
                                     array1d<R1Tensor> const & X ) override
  {
      //TODO ?
  }

  virtual void setupRelatedObjectsInRelations( MeshLevel const * const mesh ) override
  {
    //TODO ?
  }

  struct viewKeyStruct : ObjectManagerBase::viewKeyStruct
  {
    static constexpr auto elementVolumeString          = "elementVolume";
    static constexpr auto fineElementsListString       = "fineElements";
  };

  /*!
   * @brief returns the element to node relations.
   * @details The aggregates are elements composed of 1 node.
   * @param[in] k the index of the element.
   */
  virtual arraySlice1dRval<localIndex const> nodeList( localIndex const k ) const override
  { 
    return m_toNodesRelation[k];
  }

  /*!
   * @brief returns the element to node relations.
   * @details The aggregates are elements composed of 1 node.
   * @param[in] k the index of the element.
   */
  virtual arraySlice1dRval<localIndex> nodeList( localIndex const k ) override
  {
    return m_toNodesRelation[k];
  }

private:
  /// The elements to nodes relation is one to one relation.
  NodeMapType  m_toNodesRelation;

  /// Relation between fine and coarse elements ordered by aggregates
  array1d< localIndex > m_fineToCoarse;

  /// Number of fine cells per aggregate
  array1d< localIndex > m_nbFineCellsPerCoarseCell;

  /// Barycenters of aggregates
  array1d< R1Tensor > m_aggregateBarycenters;
};
}

#endif
