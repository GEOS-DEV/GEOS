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
 * @file BlockVectorView.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_UTILITIES_BLOCKVECTORVIEW_HPP_
#define GEOSX_LINEARALGEBRA_UTILITIES_BLOCKVECTORVIEW_HPP_

#include "common/DataTypes.hpp"

namespace geosx
{

/**
 * @brief This class creates and provides basic support for block vectors objects.
 * @tparam VECTOR type of vectors
 */
template< typename VECTOR >
class BlockVectorView
{
public:

  using Vector = VECTOR;

  /**
   * @brief Destructor.
   */
  virtual ~BlockVectorView() = default;

  ///@}
  //! @name Setters
  ///@{

  /**
   * @brief Update vector <tt>y</tt> as <tt>y</tt> = <tt>x</tt>.
   * @param src Vector to copy
   */
  void copy( BlockVectorView< VECTOR > const & src );

  ///@}
  //! @name Linear Algebra Operations
  ///@{

  /**
   * @brief Scale the block vector with <tt>factor</tt>.
   * @param factor Multiplication factor.
   *
   */
  void scale( real64 const factor );

  /**
   * @brief Set the vector to zero
   */
  void zero();

  /**
   * @brief Set vector elements to random entries.
   */
  void rand( unsigned const seed = 1984 );

  /**
   * @brief Dot product.
   * @param blockVec Block vector.
   * @param result Result of the dot product.
   *
   */
  real64 dot( BlockVectorView<VECTOR> const &x ) const;

  /**
   * @brief 2-norm of the block vector.
   * @param result 2-norm of the block vector.
   *
   */
  real64 norm2() const;

  /**
   * @brief Inf-norm of the block vector.
   * @param result Inf-norm of the block vector.
   */
  real64 normInf() const;

  /**
   * @brief Update vector <tt>y</tt> as <tt>y = alpha*x + y</tt>.
   * @note The naming convention follows the logic of the BLAS library.
   * @param alpha Scaling factor for added vector.
   * @param x Vector to add.
   */
  void axpy( real64 const alpha,
             BlockVectorView<VECTOR> const &x );

  /**
   * @brief Update vector <tt>y</tt> as <tt>y = alpha*x + beta*y</tt>.
   *
   * @note The naming convention follows the logic of the BLAS library.
   *
   * @param alpha Scaling factor for added vector.
   * @param x Vector to add.
   * @param beta Scaling factor for self vector.
   */
  void axpby( real64 const alpha,
              BlockVectorView<VECTOR> const &x,
              real64 const beta );

  ///@}
  //! @name Accessors
  ///@{

  /**
   * @brief Get block size.
   */
  localIndex blockSize() const;

  /**
   * @brief Get global size.
   */
  globalIndex globalSize() const;

  /**
   * @brief Print the block vector
   * @param os the stream to print to
   */
  void print( std::ostream & os = std::cout ) const;

  /**
   * @brief Get a reference to the vector corresponding to block <tt>blockRowIndex</tt>.
   *
   * This enables convenient calls to individual block operations, e.g.
   * <tt>block_vector.block(1).norm2()</tt>
   *
   * @param blockIndex Index of the block to return.
   *
   */
  VECTOR & block( localIndex const blockIndex ) const;
  ///@}

protected:

  /**
   * @name Constructor/assignment Methods
   *
   * @note Protected to prevent users from creating/copying/moving base class objects.
   *       Derived classes should still use them to size/copy/move pointer array.
   */
  ///@{

  /**
   * @brief Create a vector of @p nBlocks blocks.
   * @param nBlocks Number of blocks.
   */
  explicit BlockVectorView( localIndex const nBlocks );

  /**
   * @brief Copy constructor
   */
  BlockVectorView( BlockVectorView< VECTOR > const & x );

  /**
   * @brief Move constructor
   */
  BlockVectorView( BlockVectorView< VECTOR > && x );

  /**
   * @brief Copy assignment
   */
  BlockVectorView<VECTOR> & operator=( BlockVectorView< VECTOR > const & x );

  /**
   * @brief Move assignment
   */
  BlockVectorView<VECTOR> & operator=( BlockVectorView< VECTOR > && x );

  ///@}

  /// Resizable array of pointers to GEOSX vectors.
  array1d< VECTOR * > m_vectors;
};

/**
 * @brief Stream insertion operator
 * @param os the output stream
 * @param vec the vector to be printed
 * @return reference to the output stream
 */
template<typename VECTOR>
std::ostream & operator<<( std::ostream & os, BlockVectorView<VECTOR> const & vec )
{
  vec.print( os );
  return os;
}

} // end geosx namespace


#endif /* GEOSX_LINEARALGEBRA_UTILITIES_BLOCKVECTORVIEW_HPP_ */
