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

  using Base = BlockVectorView< VECTOR >;

  /**
   * @brief Create a vector of @p nBlocks blocks.
   * @param nBlocks Number of blocks.
   */
  explicit BlockVector( localIndex const nBlocks );

  /**
   * @brief Copy constructor that performs a deep copy of each sub-vector
   * @param rhs the block vector to copy
   */
  BlockVector( BlockVector const & rhs );

  /**
   * @brief Move constructor
   * @param rhs the block vector to move from
   */
  BlockVector( BlockVector && rhs );

  /**
   * @brief Conversion constructor from a compatible view with a deep copy of each sub-vector
   * @param rhs the block vector view to copy from
   * @note declared explicit to avoid unintended deep copying
   */
  explicit BlockVector( BlockVectorView< VECTOR > const & rhs );

  /**
   * @brief Destructor.
   */
  virtual ~BlockVector() override = default;

private:

  void setPointers();

  /// storage for actual vectors
  array1d< VECTOR > m_vectorStorage;
};

template< typename VECTOR >
BlockVector< VECTOR >::BlockVector( localIndex const nBlocks )
  : Base( nBlocks ),
  m_vectorStorage( nBlocks )
{
  setPointers();
}

template< typename VECTOR >
BlockVector< VECTOR >::BlockVector( BlockVector< VECTOR > const & rhs )
  : Base( rhs ),
  m_vectorStorage( rhs.m_vectorStorage )
{
  setPointers();
}

template< typename VECTOR >
BlockVector< VECTOR >::BlockVector( BlockVector && rhs )
  : Base( std::move( rhs ) ),
  m_vectorStorage( std::move( rhs.m_vectorStorage ) )
{
  setPointers();
}

template< typename VECTOR >
BlockVector< VECTOR >::BlockVector( BlockVectorView< VECTOR > const & rhs )
  : Base( rhs.blockSize() )
{
  for( localIndex i = 0; i < rhs.blockSize(); ++i )
  {
    m_vectorStorage.push_back( rhs.block( i ) );
  }
  setPointers();
}

template< typename VECTOR >
void BlockVector< VECTOR >::setPointers()
{
  GEOSX_LAI_ASSERT_EQ( this->blockSize(), m_vectorStorage.size() );
  for( localIndex i = 0; i < m_vectorStorage.size(); ++i )
  {
    this->setPointer( i, &m_vectorStorage[i] );
  }
}

} //namespace geosx

#endif //GEOSX_LINEARALGEBRA_UTILITIES_BLOCKVECTOR_HPP_
