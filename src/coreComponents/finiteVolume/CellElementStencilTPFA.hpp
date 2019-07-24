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
 * @file CellElementStencilTPFA.hpp
 */

#ifndef SRC_CORECOMPONENTS_FINITEVOLUME_CELLELEMENTSTENCILTPFA_HPP_
#define SRC_CORECOMPONENTS_FINITEVOLUME_CELLELEMENTSTENCILTPFA_HPP_

#include "common/DataTypes.hpp"

namespace geosx
{

/**
 * @class CellElementStencilTPFA
 *
 * Provides management of the interior stencil points when using Two-Point flux approximation.
 */
class CellElementStencilTPFA
{
public:

  /// The array type that will be used to store the indices of the stencil contributors
  using INDEX_TYPE = array2d<localIndex>;

  /// The array type that is used to store the weights of the stencil contributors
  using WEIGHT_TYPE = array2d<real64>;

  /// The array view type for the stencil indices
  using INDEX_VIEW_TYPE = arrayView2d<localIndex>;
  /// The array view to const type for the stencil indices
  using INDEX_VIEW_CONST_TYPE = arrayView2d<localIndex const>;

  /// The array view type for the stencil weights
  using WEIGHT_VIEW_TYPE = arrayView2d<real64>;
  /// The array view to const type for the stencil weights
  using WEIGHT_VIEW_CONST_TYPE = arrayView2d<real64 const>;

  /// Number of points the flux is between (always 2 for TPFA)
  static localIndex constexpr NUM_POINT_IN_FLUX = 2;

  /// Maximum number of points in a stencil (this is 2 for TPFA)
  static localIndex constexpr MAX_STENCIL_SIZE = 2;

  /// default constructor
  CellElementStencilTPFA();

  /**
   * @brief reserve the size of the stencil
   * @param[in] size the size of the stencil to reserve
   */
  void reserve( localIndex const size );

  /**
   * @brief Add an entry to the stencil
   * @param[in] numPts The number of points in the stencil entry
   * @param[in] elementRegionIndices The element region indices for each point in the stencil entry
   * @param[in] elementSubRegionIndices The element sub-region indices for each point in the stencil entry
   * @param[in] elementIndices The element indices for each point in the stencil entry
   * @param[in] weights The weights each point in the stencil entry
   * @param[in] connectorIndex The index of the connector element that the stencil acts across
   */
  void add( localIndex const numPts,
            localIndex  const * const elementRegionIndices,
            localIndex  const * const elementSubRegionIndices,
            localIndex  const * const elementIndices,
            real64 const * const weights,
            localIndex const connectorIndex );

  /**
   * @brief Zero weights for a stencil entry
   * @param[in] connectorIndex The index of the connector element that the stencil acts across for which the weights are
   *                           to be zero.
   * @return True if a valid connectorIndex was found, and had its corresponding weights set to zero.
   */
  bool zero( localIndex const connectorIndex );

  /**
   * @brief Give the number of stencil entries
   * @return The number of stencil entries
   */
  localIndex size() const { return m_elementRegionIndices.size(0); }

  /**
   * @brief Gives the number of points in a stencil entry.
   * @param[in] index of the stencil entry for which to query the size
   * @return the size of a stencil entry
   */
  constexpr localIndex stencilSize( localIndex index ) const { return MAX_STENCIL_SIZE; }

  /**
   * @brief Const access to the element regions indices
   * @return A view to const
   */
  INDEX_VIEW_CONST_TYPE const &  getElementRegionIndices() const { return m_elementRegionIndices; }

  /**
   * @brief Const access to the element subregions indices
   * @return A view to const
   */
  INDEX_VIEW_CONST_TYPE const &  getElementSubRegionIndices() const { return m_elementSubRegionIndices; }

  /**
   * @brief Const access to the element indices
   * @return A view to const
   */
  INDEX_VIEW_CONST_TYPE const &  getElementIndices() const { return m_elementIndices; }

  /**
   * @brief Const access to the stencil weights
   * @return A view to const
   */
  WEIGHT_VIEW_CONST_TYPE const & getWeights() const { return m_weights; }


private:
  /// The container for the element region indices for each point in each stencil
  INDEX_TYPE  m_elementRegionIndices;

  /// The container for the element sub region indices for each point in each stencil
  INDEX_TYPE  m_elementSubRegionIndices;

  /// The container for the element indices for each point in each stencil
  INDEX_TYPE  m_elementIndices;

  /// The container for the weights for each point in each stencil
  WEIGHT_TYPE m_weights;

  /// The map that provides the stencil index given the index of the underlying connector object.
  map<localIndex, localIndex> m_connectorIndices;

};

} /* namespace geosx */

#endif /* SRC_CORECOMPONENTS_FINITEVOLUME_CELLELEMENTSTENCILTPFA_HPP_ */
