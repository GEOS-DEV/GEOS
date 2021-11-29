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

#ifndef GEOSX_MESH_AGGREGATEELEMENTSUBREGION_HPP_
#define GEOSX_MESH_AGGREGATEELEMENTSUBREGION_HPP_

#include "ElementSubRegionBase.hpp"
#include "InterObjectRelation.hpp"



namespace geosx
{


/**
 * @class AggregateElementSubRegion
 * @brief The AggregateElementSubRegion class provides an interface to aggregate fine cells into coarse cells.
 */

class AggregateElementSubRegion : public ElementSubRegionBase
{
public:

  /// AggregateToNode map type
  using NodeMapType=FixedOneToManyRelation;

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  /**
   * @brief Return the name of the aggregate element sub-region in the object catalog.
   * @return string that contains the AggregateElementSubRegion catalog name
   */
  static const string catalogName()
  { return "AggregateCell"; }

  /**
   * @brief Provide a virtual access to catalogName().
   * @return string that contains the AggregateElementSubRegion catalog name
   */
  virtual const string getCatalogName() const override
  {
    return AggregateElementSubRegion::catalogName();
  }
  ///@}

  /**
   * @brief Find all the fine cell in a given aggregate coarse cell.
   * @tparam LAMBDA template argument
   * @param[in] aggregateIndex index of the aggregate
   * @param[in,out] lambda all the fine cells in the aggregate
   */
  template< typename LAMBDA >
  void forFineCellsInAggregate( localIndex aggregateIndex, LAMBDA lambda )
  {
    for( localIndex fineCell = m_nbFineCellsPerCoarseCell[aggregateIndex];
         fineCell < m_nbFineCellsPerCoarseCell[aggregateIndex+1]; fineCell++ )
    {
      lambda( m_fineToCoarse[fineCell] );
    }
  }


  /**
   * @brief Gives the number of fine cells of an aggregate coarse cell.
   * @param[in] aggregateIndex index of the aggregate coarse cell
   * @return the number of fine cell in the aggregate
   */
  localIndex getNbCellsPerAggregate( localIndex aggregateIndex ) const
  {
    return m_nbFineCellsPerCoarseCell[aggregateIndex + 1] - m_nbFineCellsPerCoarseCell[aggregateIndex];
  }

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor for this class.
   * @param[in] name the name of this object manager
   * @param[in] parent the parent Group
   */
  AggregateElementSubRegion( string const & name,
                             dataRepository::Group * const parent );
  /**
   * @brief Destructor.
   */
  virtual ~AggregateElementSubRegion() override;
  ///@}

  /**
   * @brief Construct the relation map between fine and coarse cells and order by aggregates.
   * @param[in] nbAggregates number of aggregate
   * @param[in] fineToCoarse index array of fine cells to be aggregated to form coarse cells
   * @param[in] barycenters coordinates of the elements center
   */
  void createFromFineToCoarseMap( localIndex nbAggregates,
                                  arrayView1d< localIndex const > const & fineToCoarse,
                                  arrayView2d< real64 const > const & barycenters );

  /**
   * @brief Accessor to the relation array between fine and coarse elements.
   * @return the relation array between fine and coarse elements ordered by aggregates
   */
  const array1d< localIndex > & getFineToCoarseMap()
  {
    return m_fineToCoarse;
  }

  virtual void calculateElementGeometricQuantities( NodeManager const & nodeManager,
                                                    FaceManager const & faceManager ) override
  {
    GEOSX_UNUSED_VAR( nodeManager );
    GEOSX_UNUSED_VAR( faceManager );
    //TODO ?
  }

  virtual void setupRelatedObjectsInRelations( MeshLevel const & mesh ) override
  {
    GEOSX_UNUSED_VAR( mesh );
    //TODO ?
  }

  /**
   * @name viewKeyStruct
   */
  ///@{

  /**
   *  @struct viewKeyStruct
   *  @brief Contains added view access key to be bound with class data member.
   */
  struct viewKeyStruct : ObjectManagerBase::viewKeyStruct
  {
    /// @cond DO_NOT_DOCUMENT
    static constexpr char const * elementVolumeString() { return "elementVolume"; }
    static constexpr char const * fineElementsListString() { return "fineElements"; }
    /// @endcond
  };
  ///@}

private:
  /// The elements to nodes relation is one to one relation.
  NodeMapType m_toNodesRelation;

  /// Relation between fine and coarse elements ordered by aggregates
  array1d< localIndex > m_fineToCoarse;

  /// Number of fine cells per aggregate
  array1d< localIndex > m_nbFineCellsPerCoarseCell;
};
}

#endif
