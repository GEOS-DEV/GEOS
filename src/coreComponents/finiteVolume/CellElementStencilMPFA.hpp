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
 * @file CellElementStencilMPFA.hpp
 */

#ifndef SRC_CORECOMPONENTS_FINITEVOLUME_CELLELEMENTSTENCILMPFA_HPP_
#define SRC_CORECOMPONENTS_FINITEVOLUME_CELLELEMENTSTENCILMPFA_HPP_

#include "StencilBase.hpp"

namespace geosx
{


/**
 * @struct CellElementStencilMPFA_Traits
 * Struct to predeclare the types and consexpr values of CellElementStencilMPFA so that they may be used in
 * StencilBase.
 */
struct CellElementStencilMPFA_Traits
{
  /// The array type that will be used to store the indices of the stencil contributors
  using IndexContainerType = ArrayOfArrays<localIndex>;

  /// The array view type for the stencil indices
  using IndexContainerViewType = ArrayOfArraysView<localIndex>;

  /// The array view to const type for the stencil indices
  using IndexContainerViewConstType = ArrayOfArraysView<localIndex const>;

  /// The array type that is used to store the weights of the stencil contributors
  using WeightContainerType = ArrayOfArrays<real64>;

  /// The array view type for the stencil weights
  using WeightContainerViewType = ArrayOfArraysView<real64>;

  /// The array view to const type for the stencil weights
  using WeightContainerViewConstType = ArrayOfArraysView<real64 const>;

  /// Number of points the flux is between (always 2)
  static localIndex constexpr NUM_POINT_IN_FLUX = 2;

  /// Maximum number of points in a stencil
  static localIndex constexpr MAX_STENCIL_SIZE = 18;

};

/**
 * @class CellElementStencilMPFA
 *
 * Provides management of the interior stencil points when using a Multi-point flux approximation.
 */
class CellElementStencilMPFA : public StencilBase<CellElementStencilMPFA_Traits,CellElementStencilMPFA>,
                               public CellElementStencilMPFA_Traits
{
public:

  /// default constructor
  CellElementStencilMPFA();

  virtual void reserve( localIndex const size ) override final;

  virtual void add( localIndex const numPts,
                    localIndex  const * const elementRegionIndices,
                    localIndex  const * const elementSubRegionIndices,
                    localIndex  const * const elementIndices,
                    real64 const * const weights,
                    localIndex const connectorIndex ) override final;

  virtual localIndex size() const override final
  { return m_elementRegionIndices.size(); }

  /**
   * @brief Gives the number of points in a stencil entry.
   * @param[in] index of the stencil entry for which to query the size
   * @return the size of a stencil entry
   */
  localIndex stencilSize( localIndex index ) const
  { return m_elementRegionIndices.sizeOfArray(index); }

};

} /* namespace geosx */

#endif /* SRC_CORECOMPONENTS_FINITEVOLUME_CELLELEMENTSTENCILMPFA_HPP_ */
