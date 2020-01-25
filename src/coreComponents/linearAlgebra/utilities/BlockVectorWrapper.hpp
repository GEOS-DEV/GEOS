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
 * @file BlockVectorWrapper.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_UTILITIES_BLOCKVECTORWRAPPER_HPP_
#define GEOSX_LINEARALGEBRA_UTILITIES_BLOCKVECTORWRAPPER_HPP_

#include "linearAlgebra/utilities/BlockVectorView.hpp"

namespace geosx
{

/**
 * @brief An 'thin' BlockVectorView that allows setting individual vector pointers.
 * @tparam VECTOR type of underlying vector
 */
template< typename VECTOR >
class BlockVectorWrapper : public BlockVectorView< VECTOR >
{
public:

  /**
   * @brief Create a vector wrapper of @p nBlocks blocks.
   * @param nBlocks number of blocks
   */
  explicit BlockVectorWrapper( localIndex const nBlocks );

  /**
   * @brief Deleted copy constructor
   * @param rhs the block vector to copy
   */
  BlockVectorWrapper( BlockVectorWrapper< VECTOR > const & rhs ) = delete;

  /**
   * @brief Deleted move constructor
   * @param rhs the block vector to move from
   */
  BlockVectorWrapper( BlockVectorWrapper< VECTOR > && rhs ) = delete;

  /**
   * @brief Copy assignment
   * @param rhs the block vector to copy
   * @return reference to current object
   */
  BlockVectorWrapper< VECTOR > & operator=( BlockVectorWrapper< VECTOR > const & rhs ) = delete;

  /**
   * @brief Move assignment
   * @param rhs the block vector to move from
   * @return reference to current object
   */
  BlockVectorWrapper< VECTOR > & operator=( BlockVectorWrapper< VECTOR > && rhs ) = delete;

  /**
   * @brief Destructor
   */
  virtual ~BlockVectorWrapper() override;

  /**
   * @brief Create a shallow copy by copying pointers.
   * @param x Vector to copy.
   */
  void shallowCopy( BlockVectorView<VECTOR> const & x );

  /**
   * @brief Set block <tt>blockIndex</tt> using <tt>vector</tt>.
   * @param blockIndex Index of the block to return.
   * @param vector Input vector to put in the block <tt>blockIndex</tt>.
   */
  void set( localIndex const blockIndex, VECTOR & x );

private:

  using BlockVectorView< VECTOR >::m_vectors;
};

} //namespace geosx

#endif //GEOSX_LINEARALGEBRA_UTILITIES_BLOCKVECTORWRAPPER_HPP_

