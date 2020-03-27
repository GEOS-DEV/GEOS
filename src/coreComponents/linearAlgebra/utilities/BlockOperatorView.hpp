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

template< typename VECTOR, typename OPERATOR = LinearOperator< VECTOR > >
class BlockOperatorView : public LinearOperator< BlockVectorView< VECTOR > >
{

public:

  /// Base type
  using Base = LinearOperator< BlockVectorView< VECTOR > >;

  /// the type of vector this linear operator operates on
  using Vector = typename Base::Vector;

  /**
   * @name Constructors/destructors
   */
  ///@{

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

  ///@}

  /**
   * @name LinearOperator interface
   */
  ///@{

  virtual void apply( BlockVectorView< VECTOR > const & x,
                      BlockVectorView< VECTOR > & b ) const override;

  virtual globalIndex numGlobalRows() const override;

  virtual globalIndex numGlobalCols() const override;

  virtual localIndex numLocalRows() const;

  virtual localIndex numLocalCols() const;

  ///@}

  /**
   * @name Getters
   */
  ///@{

  /**
   * @brief Get number of block rows
   * @return number of block rows
   */
  localIndex numBlockRows() const
  {
    return m_operators.size( 0 );
  }

  /**
   * @brief Get number of block columns
   * @return number of block columns
   */
  localIndex numBlockCols() const
  {
    return m_operators.size( 1 );
  }

  /**
   * @brief Get the matrix corresponding to block (@p blockRowIndex, @p blockColIndex).
   */
  OPERATOR const & block( localIndex const blockRowIndex, localIndex const blockColIndex ) const
  {
    return *m_operators( blockRowIndex, blockColIndex );
  }

  /**
   * @copydoc block( localIndex const, localIndex const )
   */
  OPERATOR & block( localIndex const blockRowIndex, localIndex const blockColIndex )
  {
    return *m_operators( blockRowIndex, blockColIndex );
  }

  ///@}

protected:

  /**
   * @brief Create a matrix of (@P nRows, @p nCols) blocks.
   */
  BlockOperatorView( localIndex const nRows, localIndex const nCols )
    : m_operators( nRows, nCols )
  {
    GEOSX_LAI_ASSERT_GT( nRows, 0 );
    GEOSX_LAI_ASSERT_GT( nCols, 0 );
  }

  /**
   * @brief Copy constructor
   */
  BlockOperatorView( BlockOperatorView< VECTOR, OPERATOR > const & x ) = default;

  /**
   * @brief Move constructor
   */
  BlockOperatorView( BlockOperatorView< VECTOR, OPERATOR > && x ) = default;

  /**
   * @brief Set/replace a pointer to a block
   * @param blockRowIndex row index of the block
   * @param blockColIndex column index of the block
   * @param op the new pointer
   */
  void setPointer( localIndex const blockRowIndex, localIndex const blockColIndex, OPERATOR * op )
  {
    GEOSX_LAI_ASSERT_GE( blockRowIndex, 0 );
    GEOSX_LAI_ASSERT_GT( numBlockRows(), blockRowIndex );
    GEOSX_LAI_ASSERT_GE( blockColIndex, 0 );
    GEOSX_LAI_ASSERT_GT( numBlockCols(), blockColIndex );
    m_operators( blockRowIndex, blockColIndex ) = op;
  }

private:

  /// Array of pointers to blocks
  array2d< OPERATOR * > m_operators;

};

template< typename VECTOR, typename OPERATOR >
void BlockOperatorView< VECTOR, OPERATOR >::apply( BlockVectorView< VECTOR > const & x,
                                                   BlockVectorView< VECTOR > & b ) const
{
  for( localIndex i = 0; i < m_operators.size( 0 ); i++ )
  {
    b.block( i ).zero();
    VECTOR temp( b.block( i ) );
    for( localIndex j = 0; j < m_operators.size( 1 ); j++ )
    {
      if( m_operators( i, j ) != nullptr )
      {
        m_operators( i, j )->apply( x.block( j ), temp );
        b.block( i ).axpy( 1.0, temp );
      }
    }
  }
}

template< typename VECTOR, typename OPERATOR >
globalIndex BlockOperatorView< VECTOR, OPERATOR >::numGlobalRows() const
{
  globalIndex numRows = 0;
  for( localIndex i = 0; i < numBlockRows(); ++i )
  {
    for( localIndex j = 0; j < numBlockCols(); ++j )
    {
      if( m_operators( i, j ) != nullptr )
      {
        numRows += block( i, j ).numGlobalRows();
        break;
      }
    }
  }
  return numRows;
}

template< typename VECTOR, typename OPERATOR >
globalIndex BlockOperatorView< VECTOR, OPERATOR >::numGlobalCols() const
{
  globalIndex numCols = 0;
  for( localIndex j = 0; j < numBlockCols(); j++ )
  {
    for( localIndex i = 0; i < numBlockRows(); ++i )
    {
      if( m_operators( i, j ) != nullptr )
      {
        numCols += block( i, j ).numGlobalCols();
        break;
      }
    }
  }
  return numCols;
}

template< typename VECTOR, typename OPERATOR >
localIndex BlockOperatorView< VECTOR, OPERATOR >::numLocalRows() const
{
  localIndex numRows = 0;
  for( localIndex i = 0; i < numBlockRows(); ++i )
  {
    for( localIndex j = 0; j < numBlockCols(); ++j )
    {
      if( m_operators( i, j ) != nullptr )
      {
        numRows += block( i, j ).numLocalRows();
        break;
      }
    }
  }
  return numRows;
}

template< typename VECTOR, typename OPERATOR >
localIndex BlockOperatorView< VECTOR, OPERATOR >::numLocalCols() const
{
  localIndex numCols = 0;
  for( localIndex j = 0; j < numBlockCols(); j++ )
  {
    for( localIndex i = 0; i < numBlockRows(); ++i )
    {
      if( m_operators( i, j ) != nullptr )
      {
        numCols += block( i, j ).numLocalCols();
        break;
      }
    }
  }
  return numCols;
}

}// end geosx namespace


#endif /*GEOSX_LINEARALGEBRA_UTILITIES_BLOCKOPERATORVIEW_HPP_*/
