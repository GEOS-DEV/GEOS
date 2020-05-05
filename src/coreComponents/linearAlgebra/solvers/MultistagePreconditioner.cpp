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
 * @file MultistagePreconditioner.cpp
 */

#include "MultistagePreconditioner.hpp"
#include "linearAlgebra/DofManager.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"

namespace geosx
{

template< typename LAI >
MultistagePreconditioner< LAI >::MultistagePreconditioner() = default;

template< typename LAI >
MultistagePreconditioner< LAI >::~MultistagePreconditioner() = default;

template< typename LAI >
void MultistagePreconditioner< LAI >::compute( Matrix const & mat,
                                               DofManager const & dofManager )
{
  Base::compute( mat, dofManager );
  for( localIndex k = 0; k < numStages(); ++k )
  {
    m_stages[k]->compute( mat, dofManager );
  }
  m_residual.createWithLocalSize( mat.numLocalRows(), mat.getComm() );
  m_solution.createWithLocalSize( mat.numLocalRows(), mat.getComm() );
}

template< typename LAI >
void MultistagePreconditioner< LAI >::apply( Vector const & src,
                                             Vector & dst ) const
{
  GEOSX_LAI_ASSERT_EQ( this->numGlobalRows(), dst.globalSize() );
  GEOSX_LAI_ASSERT_EQ( this->numGlobalCols(), src.globalSize() );

  m_residual.copy( src );
  dst.zero();
  for( localIndex k = 0; k < numStages(); ++k )
  {
    for( localIndex s = 0; s < m_numSteps[k]; ++s )
    {
      m_stages[k]->apply( m_residual, m_solution );
      dst.axpy( 1.0, m_solution );
      if( k < numStages() - 1 || s < m_numSteps[k] - 1 )
      {
        this->matrix().gemv( -1.0, m_solution, 1.0, m_residual );
      }
    }
  }
}

template< typename LAI >
void MultistagePreconditioner< LAI >::addStage( std::unique_ptr< PreconditionerBase< LAI > > stage,
                                                localIndex numSteps )
{
  m_stages.emplace_back( std::move( stage ) );
  m_numSteps.push_back( numSteps );
}

template< typename LAI >
void MultistagePreconditioner< LAI >::clear()
{
  m_stages.clear();
  m_numSteps.clear();
}

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOSX_USE_TRILINOS
template class MultistagePreconditioner< TrilinosInterface >;
#endif

#ifdef GEOSX_USE_HYPRE
template class MultistagePreconditioner< HypreInterface >;
#endif

#ifdef GEOSX_USE_PETSC
template class MultistagePreconditioner< PetscInterface >;
#endif

} //namespace geosx
