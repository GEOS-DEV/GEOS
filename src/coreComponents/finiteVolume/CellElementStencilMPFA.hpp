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

/**
 * @file CellElementStencilMPFA.hpp
 */

#ifndef GEOSX_FINITEVOLUME_CELLELEMENTSTENCILMPFA_HPP_
#define GEOSX_FINITEVOLUME_CELLELEMENTSTENCILMPFA_HPP_

#include "StencilBase.hpp"

namespace geos
{

/**
 * @brief Describes properties of CellElementStencilMPFA.
 *
 * This type of stencil supports connecting exactly two elements in each
 * flux with a larger number of points involved in the flux computation.
 */
using CellElementStencilMPFATraits = StencilTraits< ArrayOfArrays, 2, 18, 1 >;

/**
 * @brief Provides management of the interior stencil points when using a Multi-point flux approximation.
 */
class CellElementStencilMPFA final : public StencilBase< CellElementStencilMPFATraits, CellElementStencilMPFA >
{
public:

  virtual void reserve( localIndex const size ) override;

  virtual void add( localIndex const numPts,
                    localIndex const * const elementRegionIndices,
                    localIndex const * const elementSubRegionIndices,
                    localIndex const * const elementIndices,
                    real64 const * const weights,
                    localIndex const connectorIndex ) override final;

  /**
   * @brief Return the stencil size.
   * @return the stencil size
   */
  virtual localIndex size() const override
  { return m_elementRegionIndices.size(); }

  /**
   * @brief Give the number of points in a stencil entry.
   * @param[in] index of the stencil entry for which to query the size
   * @return the size of a stencil entry
   */
  localIndex stencilSize( localIndex index ) const
  { return m_elementRegionIndices.sizeOfArray( index ); }

  /**
   * @brief Give the number of points between which the flux is.
   * @param[in] index of the stencil entry for which to query the size
   * @return the number of points.
   */
  constexpr localIndex numPointsInFlux( localIndex index ) const
  {
    GEOS_UNUSED_VAR( index );
    return maxNumPointsInFlux;
  }

};

} /* namespace geos */

#endif /* GEOSX_FINITEVOLUME_CELLELEMENTSTENCILMPFA_HPP_ */
