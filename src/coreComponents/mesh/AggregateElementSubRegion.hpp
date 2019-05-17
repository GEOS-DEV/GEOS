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
  void forGlobalFineCellsInAggregate( localIndex aggregateIndex, LAMBDA lambda )
  {
    for(globalIndex fineCell :  m_globalFineByAggregates[aggregateIndex])
    {
      lambda(fineCell);
    }
  }

  localIndex GetNbCellsPerAggregate( localIndex aggregateIndex ) const
  {
    return m_globalFineByAggregates[aggregateIndex].size();
  }

  AggregateElementSubRegion( string const & name,
                             dataRepository::ManagedGroup * const parent );

  virtual ~AggregateElementSubRegion() override;
 
  void CreateFromFineToCoarseMap( localIndex nbAggregates,
                                  array1d< localIndex > const & fineToCoarse,
                                  array1d< R1Tensor > const & barycenters,
                                  array1d< real64 > const & volumes );

  virtual R1Tensor const & calculateElementCenter( localIndex k,
                                                   NodeManager const & nodeManager,
                                                   const bool useReferencePos = true) const override
  {
    return m_elementCenter[k];
  }

  virtual void CalculateCellVolumes( array1d<localIndex> const & indices,
                                     array1d<R1Tensor> const & X ) override
  {
  }

  virtual void setupRelatedObjectsInRelations( MeshLevel const * const mesh ) override
  {
  }

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

  /// Fine cell global index within each aggregate
  array1d< array1d < localIndex > > m_globalFineByAggregates;
};
}

#endif
