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
 * @file BlockOperatorWrapper.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_UTILITIES_BLOCKOPERATORWRAPPER_HPP_
#define GEOSX_LINEARALGEBRA_UTILITIES_BLOCKOPERATORWRAPPER_HPP_

#include "linearAlgebra/utilities/BlockOperatorView.hpp"
#include "linearAlgebra/common.hpp"

namespace geosx
{

template< typename VECTOR, typename OPERATOR >
class BlockOperatorWrapper : public BlockOperatorView< VECTOR, OPERATOR >
{
public:

  using Base = BlockOperatorView< VECTOR, OPERATOR >;
  using Vector = typename Base::Vector;

  /**
   * @brief Create a vector wrapper of @p nBlocks blocks.
   * @param nBlocks number of blocks
   */
  explicit BlockOperatorWrapper( localIndex const nRows, localIndex const nCols )
    : Base( nRows, nCols )
  {}

  /**
   * @brief Deleted copy constructor
   * @param rhs the block operator to copy
   */
  BlockOperatorWrapper( BlockOperatorWrapper const & rhs ) = delete;

  /**
   * @brief Deleted move constructor
   * @param rhs the block operator to move from
   */
  BlockOperatorWrapper( BlockOperatorWrapper && rhs ) = delete;

  /**
   * @brief Destructor
   */
  virtual ~BlockOperatorWrapper() override = default;

  /**
   * @brief Set block (@p blockRowIndex, @p blockColIndex) using @p matrix.
   */
  void set( localIndex const blockRowIndex,
            localIndex const blockColIndex,
            OPERATOR & op )
  {
    this->setPointer( blockRowIndex, blockColIndex, &op );
  }
};

} // namespace geosx

#endif //GEOSX_LINEARALGEBRA_UTILITIES_BLOCKOPERATORWRAPPER_HPP_
