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
 * @file MsrsbLevelBuilderBase.cpp
 */

#include "MsrsbLevelBuilderBase.hpp"

#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/multiscale/msrsb/MsrsbUtils.hpp"

namespace geos
{
namespace multiscale
{

template< typename LAI >
MsrsbLevelBuilderBase< LAI >::MsrsbLevelBuilderBase( string name, LinearSolverParameters params )
  : LevelBuilderBase< LAI >( std::move( name ), std::move( params ) ),
  m_dofManager( m_name )
{}

template< typename LAI >
void MsrsbLevelBuilderBase< LAI >::compute( Matrix const & fineMatrix )
{
  GEOS_MARK_FUNCTION;

  if( !m_restriction )
  {
    // Setting up finest level: just compute smoothers
    GEOS_MARK_SCOPE( setup smoother );
    m_preSmoother->setup( fineMatrix );
    m_postSmoother->setup( fineMatrix );
    return;
  }

  // Recompute coarse operator only if prolongation has changed sufficiently
  if( updateProlongation( fineMatrix ) || !m_matrix.ready() )
  {
    {
      GEOS_MARK_SCOPE( compute RAP );
      GEOS_LOG_RANK_0_IF( m_params.multiscale.debugLevel >= 2, GEOS_FMT( "[MsRSB] {}: computing RAP", m_name ) );
      m_matrix = msrsb::computeRAP( m_params.multiscale, fineMatrix, m_prolongation, *m_restriction );
    }
    if( m_params.multiscale.debugLevel >= 4 )
    {
      m_matrix.write( m_name + ".mtx", LAIOutputFormat::MATRIX_MARKET );
    }
    {
      GEOS_MARK_SCOPE( setup smoother );
      m_preSmoother->setup( m_matrix );
      m_postSmoother->setup( m_matrix );
    }
  }
}

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOSX_USE_TRILINOS
template class MsrsbLevelBuilderBase< TrilinosInterface >;
#endif

#ifdef GEOSX_USE_HYPRE
template class MsrsbLevelBuilderBase< HypreInterface >;
#endif

#ifdef GEOSX_USE_PETSC
template class MsrsbLevelBuilderBase< PetscInterface >;
#endif

} // namespace multiscale
} // namespace geos
