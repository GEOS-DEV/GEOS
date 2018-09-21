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
 * @file BlockMatrixView.hpp
 *
 *  Created on: Aug 24, 2018
 *      Author: Matthias Cremon
 */

#ifndef SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLOCKMATRIXVIEW_HPP_
#define SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLOCKMATRIXVIEW_HPP_

#include "TrilinosInterface.hpp"
//#include "HypreInterface.hpp"
#include "BlockVectorView.hpp"

namespace geosx
{

/**
 * \class BlockMatrixView
 * \brief This class creates and provides basic support for block
 *        matrices objects (templated on the LA interface).
 */

template< typename LAI >
class BlockMatrixView
{

  using ParallelMatrix = typename LAI::ParallelMatrix;
  using ParallelVector = typename LAI::ParallelVector;

public:

  //! @name Constructor/Destructor Methods
  //@{
  /**
   * @brief Empty matrix constructor.
   *
   * Create an empty block matrix.
   */
  BlockMatrixView();

  /**
   * @brief Matrix of (<tt>nRows</tt>,<tt>nCols</tt>) blocks.
   *
   * Create a block matrix of size (<tt>nRows</tt>,<tt>nCols</tt>).
   */
  BlockMatrixView( typename LAI::laiLID const nRows,
                   typename LAI::laiLID const nCols );

  /**
   * @brief Virtual destructor.
   */
  virtual ~BlockMatrixView() = default;
  //@}

  //! @name Linear Algebra Methods
  //@{
  /**
   * @brief Apply the block matrix to a block vector.
   *
   * Computes the matrix-vector product <tt>Ax = b</tt>.
   *
   * \param x Input vector.
   * \param b Output vector.
   *
   */
  void multiply( BlockVectorView<LAI> const &x,
                 BlockVectorView<LAI> &b ) const;

  /**
   * @brief Compute residual <tt>r = b - Ax</tt>.
   *
   * \param x Input solution.
   * \param b Input right hand size.
   * \param r Output residual.
   *
   */
  void residual( BlockVectorView<LAI> const &x,
                 BlockVectorView<LAI> const &b,
                 BlockVectorView<LAI> &r ) const;


  /**
   * @brief Compute generalized matrix-vector product <tt>y = alpha*Ax + beta*y</tt>.
   *
   * \param alpha Scaling factor for Ax.
   * \param x Input vector.
   * \param beta Scaling factor for y.
   * \param y Updated vector.
   *
   */
  void gemv( real64 alpha,
             BlockVectorView<LAI> const &x,
             real64 beta,
             BlockVectorView<LAI> &y ) const;

  /**
   * @brief Scale block (<tt>blockRowIndex</tt>,<tt>blockColIndex</tt>) using <tt>factor</tt>.
   *
   * \param blockRowIndex Row index.
   * \param blockColIndex Column index.
   * \param factor Scaling factor.
   *
   */
  void scale( typename LAI::laiLID const blockRowIndex,
              typename LAI::laiLID const blockColIndex,
              real64 const factor );

  /**
   * @brief Scale matrix using <tt>factor</tt>.
   */
  void scale( real64 const factor );

  /**
   * @brief Clear row and multiply the diagonal entry by <tt>factor</tt>.
   */
  void clearRow( typename LAI::laiGID const rowIndex,
                 real64 const factor );
  //@}

  //! @name Accessors/Setters
  //@{
  /**
   * @brief Get the matrix corresponding to block (<tt>blockRowIndex</tt>,<tt>blockColIndex</tt>).
   */
  ParallelMatrix * getBlock( typename LAI::laiLID const blockRowIndex,
                             typename LAI::laiLID const blockColIndex ) const;

  /**
   * @brief Get the matrix corresponding to block <tt>name</tt>.
   */
  ParallelMatrix * getBlock( std::string const blockName ) const;

  /**
   * @brief Set block (<tt>i</tt>,<tt>j</tt>) using <tt>matrix</tt>.
   */
  void setBlock( typename LAI::laiLID const blockRowIndex,
                 typename LAI::laiLID const blockColIndex,
                 ParallelMatrix &matrix );

  //@}

private:

  array2d<ParallelMatrix *> m_matrices;

};

// Empty constructor (inlined)
template< typename LAI >
inline
BlockMatrixView<LAI>::BlockMatrixView()
{}

// Constructor with a size (inlined)
template< typename LAI >
inline
BlockMatrixView<LAI>::BlockMatrixView( typename LAI::laiLID const nRows,
                                       typename LAI::laiLID const nCols )
{
  m_matrices.resize( nRows, nCols );
}

// Apply the block matrix to a block vector.
template< typename LAI >
void BlockMatrixView<LAI>::multiply( BlockVectorView<LAI> const &x,
                                     BlockVectorView<LAI> &b ) const
{
  for( typename LAI::laiLID row = 0 ; row < m_matrices.size( 0 ) ; row++ )
  {
    b.scale( row, 0. );
    ParallelVector temp( *b.getBlock( row ) );
    for( typename LAI::laiLID col = 0 ; col < m_matrices.size( 1 ) ; col++ )
    {
      if( m_matrices[row][col] != nullptr )
      {
        m_matrices[row][col]->multiply( *x.getBlock( col ), temp );
        b.axpby( row, 1.0, temp, 1.0 );
      }
    }
  }
}

// Set to residual form.
template< typename LAI >
void BlockMatrixView<LAI>::residual( BlockVectorView<LAI> const &x,
                                     BlockVectorView<LAI> const &b,
                                     BlockVectorView<LAI> &r ) const
{
  BlockVectorView<LAI> Ax( b );
  Ax.scale( 0. );
  for( typename LAI::laiLID row = 0 ; row < m_matrices.size( 0 ) ; row++ )
  {
    ParallelVector temp( *b.getBlock( row ) );
    temp.scale( 0. );
    for( typename LAI::laiLID col = 0 ; col < m_matrices.size( 1 ) ; col++ )
    {
      if( m_matrices[row][col] != nullptr )
      {
        m_matrices[row][col]->multiply( *x.getBlock( col ), temp );
        Ax.axpby( row, 1.0, temp, 1.0 );
      }
    }
  }
  // r = b
  r.axpby( 1.0, b, 0. );
  // r = r - Ax
  r.axpby( -1.0, Ax, 1.0 );
}

// Generalized matrix-vector product.
template< typename LAI >
void BlockMatrixView<LAI>::gemv( real64 alpha,
                                 BlockVectorView<LAI> const &x,
                                 real64 beta,
                                 BlockVectorView<LAI> &y ) const
{
  BlockVectorView<LAI> Ax( y );
  Ax.scale( 0. );
  for( typename LAI::laiLID row = 0 ; row < m_matrices.size( 0 ) ; row++ )
  {
    ParallelVector temp( *y.getBlock( row ) );
    temp.scale( 0. );
    for( typename LAI::laiLID col = 0 ; col < m_matrices.size( 1 ) ; col++ )
    {
      if( m_matrices[row][col] != nullptr )
      {
        m_matrices[row][col]->multiply( *x.getBlock( col ), temp );
        Ax.axpby( row, 1.0, temp, 1.0 );
      }
    }
  }
  // y = alpha*Ax + beta*y
  y.axpby( alpha, Ax, beta );
}

// Scale the matrix with the factor <tt>factor</tt>.
template< typename LAI >
void BlockMatrixView<LAI>::scale( real64 const factor )
{
  for( typename LAI::laiLID row = 0 ; row < m_matrices.size( 0 ) ; row++ )
  {
    for( typename LAI::laiLID col = 0 ; col < m_matrices.size( 1 ) ; col++ )
    {
      if( m_matrices[row][col] != nullptr )
      {
        m_matrices[row][col]->scale( factor );
      }
    }
  }
}

// Scale the block (<tt>blockRowIndex</tt>,<tt>blockRowIndex</tt>) with the factor <tt>factor</tt>.
template< typename LAI >
void BlockMatrixView<LAI>::scale( typename LAI::laiLID const blockRowIndex,
                                  typename LAI::laiLID const blockColIndex,
                                  real64 const factor )
{
  m_matrices[blockRowIndex][blockColIndex]->scale( factor );
}

// Clear row and multiply the diagonal entry by <tt>factor</tt>.
template< typename LAI >
void BlockMatrixView<LAI>::clearRow( typename LAI::laiGID const rowIndex,
                                     real64 const factor )
{}

// Accessor for block.
template< typename LAI >
typename LAI::ParallelMatrix * BlockMatrixView<LAI>::getBlock( typename LAI::laiLID const blockRowIndex,
                                                               typename LAI::laiLID const blockColIndex ) const
{
  return m_matrices[blockRowIndex][blockColIndex];
}

// Setter for block.
template< typename LAI >
void BlockMatrixView<LAI>::setBlock( typename LAI::laiLID const blockRowIndex,
                                     typename LAI::laiLID const blockColIndex,
                                     typename LAI::ParallelMatrix &matrix )
{
  m_matrices[blockRowIndex][blockColIndex] = &matrix;
}

}


#endif /* SRC_COMPONENTS_LINEARALGEBRAINTERFACE_SRC_BLOCKMATRIXVIEW_HPP_ */
