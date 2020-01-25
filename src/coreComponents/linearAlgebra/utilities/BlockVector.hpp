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
 * @file BlockVector.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_UTILITIES_BLOCKVECTOR_HPP_
#define GEOSX_LINEARALGEBRA_UTILITIES_BLOCKVECTOR_HPP_

#include "linearAlgebra/utilities/BlockVectorView.hpp"

namespace geosx
{

template< typename VECTOR >
class BlockVector : public BlockVectorView< VECTOR >
{
public:

  /**
   * @brief Vector of <tt>nBlocks</tt> blocks.
   *
   * Create a block vector of size <tt>nBlocks</tt>.
   *
   * @param nBlocks Number of blocks.
   *
   */
  explicit BlockVector( localIndex const nBlocks );

  /**
   * @brief Copy constructor that performs a deep copy of each sub-vector
   * @param rhs the block vector to copy
   */
  BlockVector( BlockVector< VECTOR > const & rhs );

  /**
   * @brief Move constructor
   * @param rhs the block vector to move from
   */
  BlockVector( BlockVector< VECTOR > && rhs ) noexcept;

  /**
   * @brief Conversion from a compatible view with a deep copy of each sub-vector
   * @param rhs the block vector view to copy from
   *
   * @note declared explicit to avoid unintended deep copying
   */
  explicit BlockVector( BlockVectorView< VECTOR > const & rhs );

  /**
   * @brief Copy assignment
   * @param rhs the block vector to copy
   * @return reference to current object
   */
  BlockVector< VECTOR > & operator=( BlockVector< VECTOR > const & rhs );

  /**
   * @brief Move assignment
   * @param rhs the block vector to move from
   * @return reference to current object
   */
  BlockVector< VECTOR > & operator=( BlockVector< VECTOR > && rhs ) noexcept;

  /**
   * @brief Desctructor.
   */
  virtual ~BlockVector() override;

private:

  using BlockVectorView< VECTOR >::m_vectors;

  void setPointers();

  array1d< VECTOR > m_vectorStorage;

};

} //namespace geosx

#endif //GEOSX_LINEARALGEBRA_UTILITIES_BLOCKVECTOR_HPP_
