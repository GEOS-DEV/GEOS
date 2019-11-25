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
 * \class BlockVectorView
 * \brief This class creates and provides basic support for block
 *        vectors objects (templated on the LA interface).
 */
template< typename LAI >
class BlockVectorView
{

  using ParallelVector = typename LAI::ParallelVector;

public:

  //! @name Constructor/Destructor Methods
  //@{
  /**
   * @brief Empty vector constructor.
   *
   * Create an empty block vector.
   */
  BlockVectorView();

  /**
   * @brief Vector of <tt>nBlocks</tt> blocks.
   *
   * Create a block vector of size <tt>nBlocks</tt>.
   *
   * \param nBlocks Number of blocks.
   *
   */
  BlockVectorView( localIndex const nBlocks );

  /* DISABLED: view usage can lead to non-intuitive behavior
   *
   * @brief Shallow copy constructor.
   *
   * Create a shallow copy of the block vector <tt>x</tt>.  Because
   * the block vector is just a view and doesn't store the actual
   * vectors, we are merely copying pointers to the original ParallelVectors.
   *
   * \param x Block vector to copy.
   *
   */
  //BlockVectorView( BlockVectorView<LAI> const &x );

  /**
   * @brief Destructor.
   */
  ~BlockVectorView() = default;

  //@}
  //! @name Setters
  //@{

  /**
   * @brief Create a shallow copy <tt>y = x</tt>.
   *
   * @note Since the BlockVector does not store any actual data, we
   *       just copy the pointers to the original ParallelVectors
   *       that make up the block.
   *
   * \param x Vector to copy.
   */
  void shallowCopy( BlockVectorView<LAI> const &x );

  /**
   * @brief Set block <tt>blockIndex</tt> using <tt>vector</tt>.
   *
   * \param blockIndex Index of the block to return.
   * \param vector Input vector to put in the block <tt>blockIndex</tt>.
   *
   */
  void set( localIndex const blockIndex,
            ParallelVector &x );

  //@}
  //! @name Linear Algebra Operations
  //@{

  /**
   * @brief Scale the block vector with <tt>factor</tt>.
   *
   * \param factor Multiplication factor.
   *
   */
  void scale( real64 const factor );

  /**
   * @brief Dot product.
   *
   * \param blockVec Block vector.
   * \param result Result of the dot product.
   *
   */
  real64 dot( BlockVectorView<LAI> const &x ) const;

  /**
   * @brief 2-norm of the block vector.
   *
   * \param result 2-norm of the block vector.
   *
   */
  real64 norm2() const;

  /**
   * @brief Inf-norm of the block vector.
   *
   * \param result Inf-norm of the block vector.
   */
  real64 normInf() const;

  /**
   * @brief Update vector <tt>y</tt> as <tt>y = alpha*x + y</tt>.
   *
   * @note The naming convention follows the logic of the BLAS library.
   *
   * \param alpha Scaling factor for added vector.
   * \param x Vector to add.
   */
  void axpy( real64 const alpha,
             BlockVectorView<LAI> const &x );

  /**
   * @brief Update vector <tt>y</tt> as <tt>y = alpha*x + beta*y</tt>.
   *
   * @note The naming convention follows the logic of the BLAS library.
   *
   * \param alpha Scaling factor for added vector.
   * \param x Vector to add.
   * \param beta Scaling factor for self vector.
   */
  void axpby( real64 const alpha,
              BlockVectorView<LAI> const &x,
              real64 const beta );

  //@}
  //! @name Accessors
  //@{

  /**
   * @brief Get block size.
   */
  localIndex blockSize() const;

  /**
   * @brief Get global size.
   */
  globalIndex globalSize() const;

  /**
   * @brief Get a reference to the vector corresponding to block <tt>blockRowIndex</tt>.
   *
   * This enables convenient calls to individual block operations, e.g.
   * <tt>block_vector.block(1).norm2()</tt>
   *
   * \param blockIndex Index of the block to return.
   *
   */
  ParallelVector & block( localIndex const blockIndex ) const;
  //@}

private:

  // Resizable array of pointers to GEOSX vectors.
  array1d<ParallelVector *> m_vectors;
};

} // end geosx namespace


#endif /* GEOSX_LINEARALGEBRA_UTILITIES_BLOCKVECTORVIEW_HPP_ */
