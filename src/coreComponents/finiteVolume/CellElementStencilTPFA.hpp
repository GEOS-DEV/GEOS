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

#include "StencilBase.hpp"

namespace geosx
{

/**
 * @struct CellElementStencilTPFA_Traits
 * Struct to predeclare the types and consexpr values of CellElementStencilTPFA so that they may be used in
 * StencilBase.
 */
struct CellElementStencilTPFA_Traits
{
  /// The array type that will be used to store the indices of the stencil contributors
  using IndexContainerType = array2d<localIndex>;

  /// The array view type for the stencil indices
  using IndexContainerViewType = arrayView2d<localIndex>;

  /// The array view to const type for the stencil indices
  using IndexContainerViewConstType = arrayView2d<localIndex const>;

  /// The array type that is used to store the weights of the stencil contributors
  using WeightContainerType = array2d<real64>;

  /// The array view type for the stencil weights
  using WeightContainerViewType = arrayView2d<real64>;

  /// The array view to const type for the stencil weights
  using WeightContainerViewConstType = arrayView2d<real64 const>;


  /// Number of points the flux is between (always 2 for TPFA)
  static constexpr localIndex NUM_POINT_IN_FLUX = 2;

  /// Maximum number of points in a stencil (this is 2 for TPFA)
  static constexpr localIndex MAX_STENCIL_SIZE = 2;
};

/**
 * @class CellElementStencilTPFA
 *
 * Provides management of the interior stencil points when using Two-Point flux approximation.
 */
class CellElementStencilTPFA : public StencilBase<CellElementStencilTPFA_Traits, CellElementStencilTPFA>,
                               public CellElementStencilTPFA_Traits
{
public:

  /// default constructor
  CellElementStencilTPFA();

  virtual void add( localIndex const numPts,
                    localIndex  const * const elementRegionIndices,
                    localIndex  const * const elementSubRegionIndices,
                    localIndex  const * const elementIndices,
                    real64 const * const weights,
                    localIndex const connectorIndex ) override final;

  virtual localIndex size() const override final
  { return m_elementRegionIndices.size(0); }

  /**
   * @brief Gives the number of points in a stencil entry.
   * @param[in] index of the stencil entry for which to query the size
   * @return the size of a stencil entry
   */
  constexpr localIndex stencilSize( localIndex index ) const
  { return MAX_STENCIL_SIZE; }

};

} /* namespace geosx */

#endif /* SRC_CORECOMPONENTS_FINITEVOLUME_CELLELEMENTSTENCILTPFA_HPP_ */
