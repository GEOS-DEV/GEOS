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

#include "linearAlgebra/common.hpp"

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
   * @brief Deleted copy assignment
   */
  BlockVectorView & operator=( BlockVectorView const & x ) = delete;

  /**
   * @brief Deleted move assignment
   */
  BlockVectorView & operator=( BlockVectorView && x ) noexcept = delete;

  /**
   * @brief Destructor.
   */
  virtual ~BlockVectorView() = default;

  ///@}

  /// @name Setters
  ///@{

  /**
   * @brief Update vector <tt>y</tt> as <tt>y</tt> = <tt>x</tt>.
   * @param src Vector to copy
   */
  void copy( BlockVectorView const & src );

  ///@}

  /// @name Linear Algebra Operations
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
  real64 dot( BlockVectorView const & x ) const;

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
             BlockVectorView const & x );

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
              BlockVectorView const & x,
              real64 const beta );

  ///@}

  /// @name Accessors
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
   * @brief Get global size.
   */
  localIndex localSize() const;

  /**
   * @brief Print the block vector
   * @param os the stream to print to
   */
  void print( std::ostream & os = std::cout ) const;

  /**
   * @brief Get a reference to the vector corresponding to block @p blockRowIndex.
   * @param blockIndex Index of the block to return.
   */
  VECTOR const & block( localIndex const blockIndex ) const;

  /**
   * @copydoc block(localIndex const) const
   */
  VECTOR & block( localIndex const blockIndex );

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
  BlockVectorView( BlockVectorView const & x ) = default;

  /**
   * @brief Move constructor
   */
  BlockVectorView( BlockVectorView && x ) = default;

  ///@}

  /**
   * @brief Set pointer to a vector
   * @param i index of vector
   * @param vec new pointer to vector
   */
  void setPointer( localIndex i, VECTOR * vec )
  {
    m_vectors( i ) = vec;
  }

private:

  /// Resizable array of pointers to GEOSX vectors.
  array1d< VECTOR * > m_vectors;
};

template< typename VECTOR >
BlockVectorView< VECTOR >::BlockVectorView( localIndex const nBlocks )
  : m_vectors( nBlocks )
{ }

template< typename VECTOR >
void BlockVectorView< VECTOR >::copy( BlockVectorView const & src )
{
  GEOSX_LAI_ASSERT_EQ( blockSize(), src.blockSize() );
  for( localIndex i = 0; i < blockSize(); i++ )
  {
    block( i ).copy( src.block( i ) );
  }
}

template< typename VECTOR >
void BlockVectorView< VECTOR >::scale( real64 const factor )
{
  for( localIndex i = 0; i < blockSize(); i++ )
  {
    block( i ).scale( factor );
  }
}

template< typename VECTOR >
void BlockVectorView< VECTOR >::zero()
{
  for( localIndex i = 0; i < blockSize(); i++ )
  {
    block( i ).zero();
  }
}

template< typename VECTOR >
void BlockVectorView< VECTOR >::rand( unsigned const seed )
{
  for( localIndex i = 0; i < blockSize(); i++ )
  {
    block( i ).rand( seed );
  }
}

template< typename VECTOR >
real64 BlockVectorView< VECTOR >::dot( BlockVectorView const & src ) const
{
  GEOSX_LAI_ASSERT_EQ( blockSize(), src.blockSize() );
  real64 accum = 0;
  for( localIndex i = 0; i < blockSize(); i++ )
  {
    accum += block( i ).dot( src.block( i ) );
  }
  return accum;
}

template< typename VECTOR >
real64 BlockVectorView< VECTOR >::norm2() const
{
  real64 accum = 0;
  for( localIndex i = 0; i < blockSize(); i++ )
  {
    real64 const temp = block( i ).norm2();
    accum = accum + temp * temp;
  }
  return std::sqrt( accum );
}

template< typename VECTOR >
real64 BlockVectorView< VECTOR >::normInf() const
{
  real64 result = 0.0;
  for( localIndex i = 0; i < blockSize(); i++ )
  {
    result = std::fmax( result, block( i ).normInf() );
  }
  return result;
}

template< typename VECTOR >
void BlockVectorView< VECTOR >::axpy( real64 const alpha,
                                      BlockVectorView const & x )
{
  GEOSX_LAI_ASSERT_EQ( blockSize(), x.blockSize() );
  for( localIndex i = 0; i < blockSize(); i++ )
  {
    block( i ).axpy( alpha, x.block( i ) );
  }
}

template< typename VECTOR >
void BlockVectorView< VECTOR >::axpby( real64 const alpha,
                                       BlockVectorView const & x,
                                       real64 const beta )
{
  GEOSX_LAI_ASSERT_EQ( blockSize(), x.blockSize() );
  for( localIndex i = 0; i < blockSize(); i++ )
  {
    block( i ).axpby( alpha, x.block( i ), beta );
  }
}

template< typename VECTOR >
localIndex BlockVectorView< VECTOR >::blockSize() const
{
  return m_vectors.size();
}

template< typename VECTOR >
globalIndex BlockVectorView< VECTOR >::globalSize() const
{
  globalIndex size = 0;
  for( localIndex i = 0; i < blockSize(); i++ )
  {
    size += block( i ).globalSize();
  }
  return size;
}

template< typename VECTOR >
localIndex BlockVectorView< VECTOR >::localSize() const
{
  localIndex size = 0;
  for( localIndex i = 0; i < blockSize(); i++ )
  {
    size += block( i ).localSize();
  }
  return size;
}

template< typename VECTOR >
VECTOR const & BlockVectorView< VECTOR >::block( localIndex blockIndex ) const
{
  GEOSX_LAI_ASSERT( m_vectors[blockIndex] != nullptr );
  return *m_vectors[blockIndex];
}

template< typename VECTOR >
VECTOR & BlockVectorView< VECTOR >::block( localIndex blockIndex )
{
  GEOSX_LAI_ASSERT( m_vectors[blockIndex] != nullptr );
  return *m_vectors[blockIndex];
}

template< typename VECTOR >
void BlockVectorView< VECTOR >::print( std::ostream & os ) const
{
  for( localIndex i = 0; i < m_vectors.size(); i++ )
  {
    os << "Block " << i << " of " << blockSize() << ":" << std::endl;
    os << "=============" << std::endl;
    block( i ).print( os );
  }
}

/**
 * @brief Stream insertion operator
 * @param os the output stream
 * @param vec the vector to be printed
 * @return reference to the output stream
 */
template< typename VECTOR >
std::ostream & operator<<( std::ostream & os,
                           BlockVectorView< VECTOR > const & vec )
{
  vec.print( os );
  return os;
}

} // end geosx namespace


#endif /* GEOSX_LINEARALGEBRA_UTILITIES_BLOCKVECTORVIEW_HPP_ */
