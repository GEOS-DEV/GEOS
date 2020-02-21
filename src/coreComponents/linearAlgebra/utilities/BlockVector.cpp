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
 * @file BlockVector.cpp
 */

#include "BlockVector.hpp"

#include "linearAlgebra/interfaces/InterfaceTypes.hpp"

namespace geosx
{

template< typename VECTOR >
BlockVector< VECTOR >::BlockVector( localIndex const nBlocks )
  : BlockVectorView< VECTOR >( nBlocks )
{
  m_vectorStorage.resize( nBlocks );
  setPointers();
}

template< typename VECTOR >
BlockVector< VECTOR >::BlockVector( BlockVector< VECTOR > const & rhs )
: BlockVectorView< VECTOR >( rhs ),
  m_vectorStorage( rhs.m_vectorStorage )
{
  setPointers();
}

template< typename VECTOR >
BlockVector< VECTOR >::BlockVector( BlockVector< VECTOR > && rhs ) noexcept
: BlockVectorView< VECTOR >( std::move(rhs) ),
  m_vectorStorage( std::move( rhs.m_vectorStorage ) )
{
  setPointers();
}

template< typename VECTOR >
BlockVector< VECTOR >::BlockVector( BlockVectorView< VECTOR > const & rhs )
: BlockVectorView< VECTOR >( rhs.blockSize() )
{
  for( localIndex i = 0; i < rhs.blockSize(); ++i )
  {
    m_vectorStorage.push_back( rhs.block( i ) );
  }
  setPointers();
}

template< typename VECTOR >
BlockVector< VECTOR >::~BlockVector()
{

}

template< typename VECTOR >
void BlockVector< VECTOR >::setPointers()
{
  GEOSX_ASSERT( m_vectors.size() == m_vectorStorage.size() );
  for( localIndex i = 0; i < m_vectorStorage.size(); ++i )
  {
    m_vectors[i] = &m_vectorStorage[i];
  }
}

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOSX_USE_TRILINOS
template class BlockVector<TrilinosInterface::ParallelVector>;
#endif

#ifdef GEOSX_USE_HYPRE
template class BlockVector<HypreInterface::ParallelVector>;
#endif

#ifdef GEOSX_USE_PETSC
template class BlockVector<PetscInterface::ParallelVector>;
#endif

} //namespace geosx

