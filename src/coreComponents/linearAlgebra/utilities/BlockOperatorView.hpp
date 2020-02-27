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
 * @file BlockOperatorView.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_UTILITIES_BLOCKOPERATORVIEW_HPP_
#define GEOSX_LINEARALGEBRA_UTILITIES_BLOCKOPERATORVIEW_HPP_

#include "linearAlgebra/interfaces/LinearOperator.hpp"
#include "linearAlgebra/utilities/BlockVectorView.hpp"

namespace geosx
{

/**
 * @brief This class creates and provides basic support for block operator objects.
 * @tparam Vector type of vector that sub-blocks of this view can operate on
 */

template< typename VECTOR, typename OPERATOR = LinearOperator<VECTOR> >
class BlockOperatorView : public LinearOperator< BlockVectorView<VECTOR> >
{

public:

  /// the type of vector this linear operator operates on
  using Vector   = BlockVectorView<VECTOR>;

  /// the underlying operator type for each block
  using Operator = OPERATOR;

  /**
   * @brief Destructor.
   */
  virtual ~BlockOperatorView() override = default;

  /**
   * @brief Deleted copy assignment
   */
  BlockOperatorView & operator=( BlockOperatorView const & rhs ) = delete;

  /**
   * @brief Deleted move assignment
   */
  BlockOperatorView & operator=( BlockOperatorView && rhs ) = delete;
  
  /**
   * @brief Apply the block matrix to a block vector.
   *
   * Computes the matrix-vector product <tt>Ax = b</tt>.
   *
   * @param x Input vector.
   * @param b Output vector.
   *
   */
  virtual void multiply( BlockVectorView<VECTOR> const & x,
                         BlockVectorView<VECTOR> & b ) const override;

  /**
   * @brief
   * @return number of block rows
   */
  localIndex numBlockRows() const
  {
    return m_operators.size(0);
  }

  /**
   * @brief
   * @return number of block columns
   */
  localIndex numBlockCols() const
  {
    return m_operators.size(1);
  }

  /**
   * @brief Get the matrix corresponding to block (@p blockRowIndex, @p blockColIndex).
   */
   Operator const & block( localIndex const blockRowIndex, localIndex const blockColIndex ) const
  {
    return *m_operators( blockRowIndex, blockColIndex );
  }

  /**
   * @copydoc block( localIndex const, localIndex const )
   */
  Operator & block( localIndex const blockRowIndex, localIndex const blockColIndex )
  {
    return *m_operators( blockRowIndex, blockColIndex );
  }

protected:

  /**
   * @brief Create a matrix of (@P nRows, @p nCols) blocks.
   */
  BlockOperatorView( localIndex const nRows, localIndex const nCols )
    : m_operators( nRows, nCols )
  {}

  /**
   * @brief Copy constructor
   */
  BlockOperatorView( BlockOperatorView< VECTOR, OPERATOR > const & x ) = default;

  /**
   * @brief Move constructor
   */
  BlockOperatorView( BlockOperatorView< VECTOR, OPERATOR > && x ) noexcept = default;

  /// Array of pointers to blocks
  array2d< Operator * > m_operators;

};

template< typename VECTOR, typename OPERATOR >
void BlockOperatorView< VECTOR, OPERATOR >::multiply( BlockVectorView<VECTOR> const & x,
                                                      BlockVectorView<VECTOR> & b ) const
{
  for( localIndex row = 0; row < m_operators.size( 0 ); row++ )
  {
    b.block( row ).zero();
    VECTOR temp( b.block( row ) );
    for( localIndex col = 0; col < m_operators.size( 1 ); col++ )
    {
      if( m_operators[row][col] != nullptr )
      {
        m_operators[row][col]->multiply( x.block( col ), temp );
        b.block( row ).axpy( 1.0, temp );
      }
    }
  }
}

}// end geosx namespace


#endif /*GEOSX_LINEARALGEBRA_UTILITIES_BLOCKOPERATORVIEW_HPP_*/
