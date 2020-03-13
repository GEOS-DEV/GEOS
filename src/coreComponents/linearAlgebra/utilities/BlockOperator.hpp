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
 * @file BlockOperator.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_UTILITIES_BLOCKOPERATOR_HPP_
#define GEOSX_LINEARALGEBRA_UTILITIES_BLOCKOPERATOR_HPP_

#include "linearAlgebra/utilities/BlockOperatorView.hpp"

namespace geosx
{

template< typename VECTOR, typename OPERATOR >
class BlockOperator : public BlockOperatorView< VECTOR, OPERATOR >
{
public:

  /// The base class
  using Base = BlockOperatorView< VECTOR, OPERATOR >;

  /// The type of vector this linear operator operates on
  using Vector = typename Base::Vector;

  /**
   * @brief Create an operator with (@p nRows, @p nCols) blocks.
   */
  BlockOperator( localIndex const nRows, localIndex const nCols );

  /**
   * @brief Copy constructor
   * @param rhs the block operator to copy from
   */
  BlockOperator( BlockOperator const & rhs );

  /**
   * @brief Move constructor
   * @param rhs the block operator to move from
   */
  BlockOperator( BlockOperator && rhs );

  /**
   * @brief Destructor.
   */
  virtual ~BlockOperator() override = default;

private:

  void setPointers();

  /// Actual storage for blocks
  array2d< OPERATOR > m_operatorStorage;
};

template< typename VECTOR, typename OPERATOR >
BlockOperator< VECTOR, OPERATOR >::BlockOperator( localIndex const nRows, localIndex const nCols )
  : Base( nRows, nCols ),
  m_operatorStorage( nRows, nCols )
{
  setPointers();
}

template< typename VECTOR, typename OPERATOR >
void BlockOperator< VECTOR, OPERATOR >::setPointers()
{
  GEOSX_LAI_ASSERT_EQ( this->numBlockRows(), m_operatorStorage.size( 0 ) );
  GEOSX_LAI_ASSERT_EQ( this->numBlockCols(), m_operatorStorage.size( 1 ) );
  for( localIndex i = 0; i < m_operatorStorage.size( 0 ); ++i )
  {
    for( localIndex j = 0; j < m_operatorStorage.size( 1 ); ++j )
    {
      this->setPointer( i, j, &m_operatorStorage( i, j ) );
    }
  }
}

template< typename VECTOR, typename OPERATOR >
BlockOperator< VECTOR, OPERATOR >::BlockOperator( BlockOperator const & rhs )
  : Base( rhs ),
  m_operatorStorage( rhs.m_operatorStorage )
{
  setPointers();
}

template< typename VECTOR, typename OPERATOR >
BlockOperator< VECTOR, OPERATOR >::BlockOperator( BlockOperator && rhs )
  : Base( std::move( rhs ) ),
  m_operatorStorage( std::move( rhs.m_operatorStorage ) )
{
  setPointers();
}

} // namespace geosx

#endif //GEOSX_LINEARALGEBRA_UTILITIES_BLOCKOPERATOR_HPP_
