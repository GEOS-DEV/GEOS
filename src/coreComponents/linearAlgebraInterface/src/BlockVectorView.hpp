/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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
 * @file BlockVectorView.hpp
 *
 *  Created on: Aug 24, 2018
 *      Author: Matthias Cremon
 */

#ifndef SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLOCKVECTORVIEW_HPP_
#define SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLOCKVECTORVIEW_HPP_

#include "TrilinosInterface.hpp"
//#include "HypreInterface.hpp"

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
  BlockVectorView( typename LAI::laiLID const nBlocks );

  /**
   * @brief Copy constructor.
   *
   * Create a copy of the block vector <tt>in_blockVec</tt>.
   *
   * \param in_blockVec Block vector to copy.
   *
   */
  BlockVectorView( BlockVectorView<LAI> const &in_blockVec );

  /**
   * @brief Virtual destructor.
   */
  virtual ~BlockVectorView() = default;
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
   * @brief Scale the vector corresponding to block (<tt>blockIndex</tt>)
   * with <tt>factor</tt>.
   *
   * \param blockIndex Index of the block to scale.
   * \param factor Multiplication factor.
   *
   */
  void scale( typename LAI::laiLID const blockIndex,
              real64 const factor );

  /**
   * @brief Dot product.
   *
   * \param blockVec Block vector.
   * \param result Result of the dot product.
   *
   */
  void dot( BlockVectorView<LAI> const &blockVec,
            real64 &result ) const;

  /**
   * @brief 2-norm of the block vector.
   *
   * \param result 2-norm of the block vector.
   *
   */
  void norm2( real64 &result ) const;

  /**
   * @brief 2-norm of the vector corresponding to block (<tt>blockIndex</tt>).
   *
   * \param blockIndex Index of the block to scale.
   * \param result 2-norm of the <tt>blockIndex</tt>-th vector.
   *
   */
  void norm2( typename LAI::laiLID const blockIndex,
              real64 &result ) const;

  /**
   * @brief Inf-norm of the block vector.
   *
   * \param result Inf-norm of the block vector.
   */
  void normInf( real64 &result ) const;

  /**
   * @brief Inf-norm of the vector corresponding to block (<tt>blockIndex</tt>).
   *
   * \param blockIndex Index of the block to scale.
   * \param result 2-norm of the <tt>blockIndex</tt>-th vector.
   *
   */
  void normInf( typename LAI::laiLID blockIndex,
                real64 &result ) const;

  /**
   * @brief Update vector <tt>y</tt> as <tt>y = x</tt>.
   *
   * @note The naming convention follows the logic of the BLAS library.
   *
   * \param x Vector to copy.
   */
  void copy( BlockVectorView<LAI> const &x );

  /**
   * @brief Update vector <tt>y(blockIndex)</tt> as <tt>y(blockIndex) = x</tt>.
   *
   * @note The naming convention follows the logic of the BLAS library.
   *
   * \param x Vector to copy.
   */
  void copy( typename LAI::laiLID blockIndex,
             ParallelVector const &x );

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
   * @brief Update the vector corresponding to block (<tt>blockIndex</tt>)
   *  as <tt>y(blockIndex) = alpha*x + beta*y(blockIndex)</tt>.
   *
   * @note The naming convention follows the logic of the BLAS library.
   *
   * \param alpha Scaling factor for added vector.
   * \param x Vector to add.
   *
   */
  void axpy( typename LAI::laiLID const blockIndex,
              real64 const alpha,
              ParallelVector const &x );

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

  /**
   * @brief Update the vector corresponding to block (<tt>blockIndex</tt>)
   *  as <tt>y(blockIndex) = alpha*x + beta*y(blockIndex)</tt>.
   *
   * @note The naming convention follows the logic of the BLAS library.
   *
   * \param alpha Scaling factor for added vector.
   * \param x Vector to add.
   * \param beta Scaling factor for self vector.
   *
   */
  void axpby( typename LAI::laiLID const blockIndex,
              real64 const alpha,
              ParallelVector const &x,
              real64 const beta );

  //@}

  //! @name Accessors/Setters
  //@{

  /**
   * @brief Get block size.
   */
  typename LAI::laiLID blockSize() const;

  /**
   * @brief Get global size.
   */
  typename LAI::laiGID globalSize() const;

  /**
   * @brief Get the vector corresponding to block (<tt>blockRowIndex</tt>,<tt>blockColIndex</tt>).
   *
   * \param blockIndex Index of the block to return.
   *
   */
  ParallelVector * getBlock( typename LAI::laiLID const blockIndex ) const;

  /**
   * @brief Get the vector corresponding to block <tt>name</tt>.
   *
   * \param blockName Name of the block to return.
   *
   */
  ParallelVector * getBlock( std::string const blockName ) const;

  /**
   * @brief Set block <tt>blockIndex</tt> using <tt>vector</tt>.
   *
   * \param blockIndex Index of the block to return.
   * \param vector Input vector to put in the block <tt>blockIndex</tt>.
   *
   */
  void setBlock( typename LAI::laiLID const blockIndex,
                 ParallelVector &vector );

  //@}

private:

  // Resizable array of pointers to GEOSX vectors.
  array1d<ParallelVector *> m_vectors;

};

// Empty constructor (inlined)
template< typename LAI >
inline
BlockVectorView<LAI>::BlockVectorView()
{}

// Constructor with a size
template< typename LAI >
BlockVectorView<LAI>::BlockVectorView( typename LAI::laiLID const nBlocks )
{
  m_vectors.resize( nBlocks );
}

// Copy constructor
template< typename LAI >
BlockVectorView<LAI>::BlockVectorView( BlockVectorView<LAI> const &in_blockVec )
{
  m_vectors.resize( in_blockVec.blockSize() );
  for( typename LAI::laiLID i = 0 ; i < m_vectors.size() ; i++ )
  {
    //ParallelVector tempVec( *in_blockVec.getBlock( i ) );
    m_vectors[i] = new ParallelVector ( *in_blockVec.getBlock( i ));
  }
}

// Scale a block.
template< typename LAI >
void BlockVectorView<LAI>::scale( real64 const factor )
{
  for( typename LAI::laiLID i = 0 ; i < m_vectors.size() ; i++ )
    m_vectors[i]->scale( factor );
}

// Scale a block.
template< typename LAI >
void BlockVectorView<LAI>::scale( typename LAI::laiLID const blockIndex,
                                  real64 const factor )
{
  m_vectors[blockIndex]->scale( factor );
}

// Dot product.
template< typename LAI >
void BlockVectorView<LAI>::dot( BlockVectorView<LAI> const &blockVec,
                                real64 &result ) const
{
  real64 accum = 0;
  for( typename LAI::laiLID i = 0 ; i < m_vectors.size() ; i++ )
  {
    real64 temp = 0;
    m_vectors[i]->dot( *blockVec.getBlock( i ), temp );
    accum = accum + temp;
  }
  result = accum;
}

// Compute the 2 norm of one block vector.
template< typename LAI >
void BlockVectorView<LAI>::norm2( typename LAI::laiLID const blockIndex,
                                  real64 &result ) const
{
  m_vectors[blockIndex]->norm2( result );
}

// Compute the 2 norm of the block vector.
template< typename LAI >
void BlockVectorView<LAI>::norm2( real64 &result ) const
{
  real64 accum = 0;
  for( typename LAI::laiLID i = 0 ; i < m_vectors.size() ; i++ )
  {
    real64 temp;
    m_vectors[i]->norm2( temp );
    accum = accum + temp*temp;
  }
  result = std::sqrt( accum );
}

// Compute the 2 norm of one block vector.
template< typename LAI >
void BlockVectorView<LAI>::normInf( typename LAI::laiLID const blockIndex,
                                    real64 &result ) const
{
  m_vectors[blockIndex]->normInf( result );
}

// Compute the 2 norm of the block vector.
template< typename LAI >
void BlockVectorView<LAI>::normInf( real64 &result ) const
{
  real64 temp;
  result = 0;
  temp = 0;
  for( typename LAI::laiLID i = 0 ; i < m_vectors.size() ; i++ )
  {
    m_vectors[i]->normInf( temp );
    if( temp > result )
      result = temp;
  }
}

// Copy the vector.
template< typename LAI >
void BlockVectorView<LAI>::copy( BlockVectorView<LAI> const &x )
{
  for( typename LAI::laiLID i = 0 ; i < m_vectors.size() ; i++ )
    m_vectors[i]->copy( *x.getBlock( i ) );
}

// Copy a ParallelVector in a block.
template< typename LAI >
void BlockVectorView<LAI>::copy( typename LAI::laiLID blockIndex,
                                 ParallelVector const &x )
{
  m_vectors[blockIndex]->copy( x );
}

// Update the vector.
template< typename LAI >
void BlockVectorView<LAI>::axpy( real64 const alpha,
                                 BlockVectorView<LAI> const &x )
{
  for( typename LAI::laiLID i = 0 ; i < m_vectors.size() ; i++ )
    m_vectors[i]->axpby( alpha, *x.getBlock( i ), 1. );
}

// Update a block.
template< typename LAI >
void BlockVectorView<LAI>::axpy( typename LAI::laiLID blockIndex,
                                 real64 const alpha,
                                 ParallelVector const &x )
{
  m_vectors[blockIndex]->axpby( alpha, x, 1. );
}

// Update the vector.
template< typename LAI >
void BlockVectorView<LAI>::axpby( real64 const alpha,
                                  BlockVectorView<LAI> const &x,
                                  real64 const beta )
{
  for( typename LAI::laiLID i = 0 ; i < m_vectors.size() ; i++ )
    m_vectors[i]->axpby( alpha, *x.getBlock( i ), beta );
}

// Update a block.
template< typename LAI >
void BlockVectorView<LAI>::axpby( typename LAI::laiLID blockIndex,
                                  real64 const alpha,
                                  ParallelVector const &x,
                                  real64 const beta )
{
  m_vectors[blockIndex]->axpby( alpha, x, beta );
}

// Setter for block.
template< typename LAI >
void BlockVectorView<LAI>::setBlock( typename LAI::laiLID const blockIndex,
                                     ParallelVector &vector )
{
  m_vectors[blockIndex] = &vector;
}

// Accessor for block size (number of blocks)
template< typename LAI >
typename LAI::laiLID BlockVectorView<LAI>::blockSize() const
{
  return static_cast<typename LAI::laiLID>(m_vectors.size());
}

// Accessor for size (total number of  elements)
template< typename LAI >
typename LAI::laiGID BlockVectorView<LAI>::globalSize(  ) const
{
  typename LAI::laiGID size = 0;
  for( typename LAI::laiLID i = 0 ; i < m_vectors.size() ; i++ )
    size = m_vectors[i]->globalSize() + size;
  return size;
}

// Accessor for block.
template< typename LAI >
typename LAI::ParallelVector * BlockVectorView<LAI>::getBlock( typename LAI::laiLID blockIndex ) const
{
  return m_vectors[blockIndex];
}

}


#endif /* SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLOCKMATRIXVIEW_HPP_ */
