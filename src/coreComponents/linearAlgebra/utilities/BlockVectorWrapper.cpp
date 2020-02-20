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
 * @file BlockVectorWrapper.cpp
 */

#include "BlockVectorWrapper.hpp"

#include "linearAlgebra/interfaces/InterfaceTypes.hpp"

namespace geosx
{

template< typename VECTOR >
BlockVectorWrapper< VECTOR >::BlockVectorWrapper( localIndex const nBlocks )
  : BlockVectorView<VECTOR>( nBlocks )
{

}

template< typename VECTOR >
BlockVectorWrapper< VECTOR >::~BlockVectorWrapper( ) = default;

template< typename VECTOR >
void BlockVectorWrapper< VECTOR >::shallowCopy( BlockVectorView< VECTOR > const & src )
{
  m_vectors.resize( src.blockSize()  );
  for( localIndex i = 0 ; i < m_vectors.size() ; i++ )
  {
    m_vectors[i] = &( src.block( i ) );
  }
}

template< typename VECTOR >
void BlockVectorWrapper< VECTOR >::set( localIndex const blockIndex,
                                        VECTOR & vector )
{
  m_vectors[blockIndex] = &vector;
}

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOSX_USE_TRILINOS
template class BlockVectorWrapper<TrilinosInterface::ParallelVector>;
#endif

#ifdef GEOSX_USE_HYPRE
template class BlockVectorWrapper<HypreInterface::ParallelVector>;
#endif

#ifdef GEOSX_USE_PETSC
template class BlockVectorWrapper<PetscInterface::ParallelVector>;
#endif

} //namespace geosx