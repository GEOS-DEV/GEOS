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
 * @file CellElementStencilTPFA.hpp
 */

#ifndef GEOSX_FINITEVOLUME_STENCILBASE_HPP_
#define GEOSX_FINITEVOLUME_STENCILBASE_HPP_

#include "common/DataTypes.hpp"
#include "codingUtilities/Utilities.hpp"

namespace geosx
{

/**
 * @class StencilBase
 *
 * Provides management of the interior stencil points when using Two-Point flux approximation.
 */
template< typename LEAFCLASSTRAITS, typename LEAFCLASS >
class StencilBase
{
public:

  StencilBase():
    m_elementRegionIndices(),
    m_elementSubRegionIndices(),
    m_elementIndices(),
    m_weights(),
    m_connectorIndices()
  {}

  virtual ~StencilBase() = default;
  StencilBase( StencilBase const & ) = default;
  StencilBase( StencilBase && ) = default;
  /**
   * @brief reserve the size of the stencil
   * @param[in] size the size of the stencil to reserve
   */
  virtual void reserve( localIndex const size );

  /**
   * @brief Add an entry to the stencil
   * @param[in] numPts The number of points in the stencil entry
   * @param[in] elementRegionIndices The element region indices for each point in the stencil entry
   * @param[in] elementSubRegionIndices The element sub-region indices for each point in the stencil entry
   * @param[in] elementIndices The element indices for each point in the stencil entry
   * @param[in] weights The weights each point in the stencil entry
   * @param[in] connectorIndex The index of the connector element that the stencil acts across
   */
  virtual void add( localIndex const numPts,
                    localIndex const * const elementRegionIndices,
                    localIndex const * const elementSubRegionIndices,
                    localIndex const * const elementIndices,
                    real64 const * const weights,
                    localIndex const connectorIndex ) = 0;

  /**
   * @brief Zero weights for a stencil entry
   * @param[in] connectorIndex The index of the connector element that the stencil acts across for which the weights are
   *                           to be zero.
   * @return True if a valid connectorIndex was found, and had its corresponding weights set to zero.
   */
  virtual bool zero( localIndex const connectorIndex );

  /**
   * @brief Give the number of stencil entries
   * @return The number of stencil entries
   */
  virtual localIndex size() const = 0;


  /**
   * @brief Const access to the element regions indices
   * @return A view to const
   */
  typename LEAFCLASSTRAITS::IndexContainerViewConstType const & getElementRegionIndices() const { return m_elementRegionIndices; }

  /**
   * @brief Const access to the element subregions indices
   * @return A view to const
   */
  typename LEAFCLASSTRAITS::IndexContainerViewConstType const & getElementSubRegionIndices() const { return m_elementSubRegionIndices; }

  /**
   * @brief Const access to the element indices
   * @return A view to const
   */
  typename LEAFCLASSTRAITS::IndexContainerViewConstType const & getElementIndices() const { return m_elementIndices; }

  /**
   * @brief Const access to the stencil weights
   * @return A view to const
   */
  typename LEAFCLASSTRAITS::WeightContainerViewConstType const & getWeights() const { return m_weights; }

protected:
  /// The container for the element region indices for each point in each stencil
  typename LEAFCLASSTRAITS::IndexContainerType m_elementRegionIndices;

  /// The container for the element sub region indices for each point in each stencil
  typename LEAFCLASSTRAITS::IndexContainerType m_elementSubRegionIndices;

  /// The container for the element indices for each point in each stencil
  typename LEAFCLASSTRAITS::IndexContainerType m_elementIndices;

  /// The container for the weights for each point in each stencil
  typename LEAFCLASSTRAITS::WeightContainerType m_weights;

  /// The map that provides the stencil index given the index of the underlying connector object.
  map< localIndex, localIndex > m_connectorIndices;

};



template< typename LEAFCLASSTRAITS, typename LEAFCLASS >
void StencilBase< LEAFCLASSTRAITS, LEAFCLASS >::reserve( localIndex const size )
{
  m_elementRegionIndices.reserve( size * 2 );
  m_elementSubRegionIndices.reserve( size * 2 );
  m_elementIndices.reserve( size * 2 );
  m_weights.reserve( size * 2 );
}


template< typename LEAFCLASSTRAITS, typename LEAFCLASS >
bool StencilBase< LEAFCLASSTRAITS, LEAFCLASS >::zero( localIndex const connectorIndex )
{
  return
    executeOnMapValue( m_connectorIndices, connectorIndex, [&]( localIndex const connectionListIndex )
  {
    for( localIndex i = 0; i < static_cast< LEAFCLASS * >(this)->stencilSize( connectorIndex ); ++i )
    {
      m_weights[connectionListIndex][i] = 0;
    }
  } );
}
} /* namespace geosx */

#endif /* GEOSX_FINITEVOLUME_STENCILBASE_HPP_ */
