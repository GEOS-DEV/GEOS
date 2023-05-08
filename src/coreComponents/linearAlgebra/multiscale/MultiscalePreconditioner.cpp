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
 * @file MultiscalePreconditioner.cpp
 */

#include "MultiscalePreconditioner.hpp"

#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/multiscale/LevelBuilderBase.hpp"

namespace geos
{

template< typename LAI >
MultiscalePreconditioner< LAI >::MultiscalePreconditioner( LinearSolverParameters params,
                                                           DomainPartition & domain )
  : Base(),
  m_params( std::move( params ) ),
  m_domain( domain )
{}

template< typename LAI >
MultiscalePreconditioner< LAI >::~MultiscalePreconditioner() = default;

template< typename LAI >
void MultiscalePreconditioner< LAI >::printLevelInfo() const
{
  constexpr char const lineFormat[] = "{:>2}  {:>10}  {:>12}  {:>7.4f}  {:>7.2f}  {:>7.2f}\n";
  constexpr char const headFormat[] = "{:>2}  {:>10}  {:>12}  {:>7}  {:>7}  {:>7}\n";

  std::ostringstream os;
  string const header = GEOS_FMT( headFormat, "L", "rows", "entries", "sparse", "nnz/row", "ratio" );
  os << "\nOperators:\n" << header;
  os << string( header.length() - 1, '=' ) << "\n";

  globalIndex totalNumRows = 0;
  globalIndex totalNumNonzeros = 0;
  globalIndex totalMemory = 0;

  for( size_t levelIndex = 0; levelIndex < m_levels.size(); ++levelIndex )
  {
    LevelData const & level = m_levels[levelIndex];
    globalIndex const nrow = level.matrix->numGlobalRows();
    globalIndex const nnz = level.matrix->numGlobalNonzeros();
    globalIndex const prevNrow = levelIndex > 0 ? m_levels[levelIndex-1].matrix->numGlobalRows() : nrow;
    os << GEOS_FMT( lineFormat, levelIndex, nrow, nnz, real64( nnz ) / ( nrow * nrow ), real64( nnz ) / nrow, real64( prevNrow ) / nrow );

    totalNumRows += nrow;
    totalNumNonzeros += nnz;
    totalMemory += nnz;
    if( levelIndex > 0 )
    {
      totalMemory += level.builder->prolongation().numGlobalNonzeros() + level.builder->restriction().numGlobalNonzeros();
    }
  }

  Matrix const & fineMat = *m_levels[0].matrix;
  constexpr char const compFormat[] = "  {:>8} = {:>6.4f}\n";
  os << "\nComplexities:\n";
  os << GEOS_FMT( compFormat, "grid", real64( totalNumRows ) / fineMat.numGlobalRows() );
  os << GEOS_FMT( compFormat, "operator", real64( totalNumNonzeros ) / fineMat.numGlobalNonzeros() );
  os << GEOS_FMT( compFormat, "memory", real64( totalMemory ) / fineMat.numGlobalNonzeros() );

  GEOS_LOG_RANK_0( os.str() );
}

template< typename LAI >
void MultiscalePreconditioner< LAI >::logMessage( integer const minLevel, string const & msg ) const
{
  GEOS_LOG_RANK_0_IF( m_params.logLevel >= minLevel,
                      GEOS_FMT( "[Multiscale] {}: {}", m_params.multiscale.label, msg ) );
}

template< typename LAI >
void MultiscalePreconditioner< LAI >::computeLevel( integer const level ) const
{
  CALI_CXX_MARK_SCOPE( GEOS_FMT( "Level {}", level ).c_str() );
  logMessage( 3, GEOS_FMT( "computing operators for level {}", level ) );
  m_levels[level].builder->compute( *m_levels[std::max( level - 1, 0 )].matrix );
}

template< typename LAI >
void MultiscalePreconditioner< LAI >::createLevels( Matrix const & mat )
{
  m_levels.clear();

  auto const levelName = [&]( integer const level ) { return GEOS_FMT( "{}_L{}", m_params.multiscale.label, level ); };

  // create fine level
  {
    integer const level = 0;
    m_levels.emplace_back();
    LevelData & curr = m_levels[level];
    curr.builder = multiscale::LevelBuilderBase< LAI >::create( levelName( level ), m_params );
    curr.builder->initializeFineLevel( m_domain, *mat.dofManager(), mat.comm() );
    curr.matrix = &mat;
    computeLevel( level );
  }

  // create coarse levels
  for( integer level = 1; level < m_params.multiscale.maxLevels; ++level )
  {
    m_levels.emplace_back();
    LevelData & curr = m_levels[level];
    LevelData & prev = m_levels[level - 1];
    curr.builder = multiscale::LevelBuilderBase< LAI >::create( levelName( level ), m_params );
    curr.builder->initializeCoarseLevel( *prev.builder, *prev.matrix );
    computeLevel( level );
    curr.matrix = &curr.builder->matrix();

    if( curr.matrix->numGlobalRows() == prev.matrix->numGlobalRows() )
    {
      // Coarsening no longer effective, stop do not produce more levels
      m_levels.pop_back();
      break;
    }
    if( curr.matrix->numGlobalRows() <= m_params.multiscale.coarsening.maxCoarseDof )
    {
      // Prevent further coarsening globally by truncating level hierarchy
      break;
    }
  }

  // create coarse solver
  m_coarse_solver = m_levels.back().builder->makeCoarseSolver();

  // create vectors
  for( LevelData & level : m_levels )
  {
    level.rhs.create( level.matrix->numLocalRows(), level.matrix->comm() );
    level.sol.create( level.matrix->numLocalRows(), level.matrix->comm() );
    level.tmp.create( level.matrix->numLocalRows(), level.matrix->comm() );
  }
}

template< typename LAI >
void MultiscalePreconditioner< LAI >::setup( Matrix const & mat )
{
  GEOS_MARK_FUNCTION;

  Base::setup( mat );

  if( !m_initialized )
  {
    logMessage( 3, "creating level hierarchy" );
    createLevels( mat );
    m_initialized = true;
  }
  else
  {
    m_levels[0].matrix = &mat;
    for( integer level = 0; level < static_cast< integer >( m_levels.size() ); ++level )
    {
      computeLevel( level );
    }
  }

  // setup coarse solver
  {
    GEOS_MARK_SCOPE( coarse solver );
    logMessage( 3, "setting up coarse solver" );
    m_coarse_solver->setup( *m_levels.back().matrix );
  }

  logMessage( 1, "statistics:" );
  if( m_params.logLevel >= 3 )
  {
    printLevelInfo();
  }
}

template< typename VECTOR >
class NormTracker
{
public:

  using Vector = VECTOR;

  NormTracker( Vector const & vec, char const * const label, int const level, bool const active )
    : m_vec( vec ),
    m_label( label ),
    m_level( level ),
    m_active( active )
  {
    if( m_active )
    {
      m_normInit = m_vec.norm2();
    }
  }

  ~NormTracker()
  {
    if( m_active )
    {
      real64 const normFinal = m_vec.norm2();
      GEOS_LOG_RANK_0( GEOS_FMT( "\tLevel {} {}: {:e} -> {:e} | x {:>9.6f}",
                                 m_level, m_label, m_normInit, normFinal, normFinal / m_normInit ) );
    }
  }

private:

  Vector const & m_vec;
  char const * m_label;
  int m_level;
  bool m_active;
  real64 m_normInit{};
};

template< typename LAI >
void MultiscalePreconditioner< LAI >::apply( Vector const & src,
                                             Vector & dst ) const
{
  GEOS_MARK_FUNCTION;

  // TODO: remove hardcoded V-cycle, abstract into a separate component

  auto const applySmoother = [&]( LevelData const & level, Operator const * smoother )
  {
    GEOS_MARK_SCOPE( smoother );
    for( integer s = 0; s < m_params.multiscale.smoother.numSweeps; ++s )
    {
      smoother->apply( level.rhs, level.tmp );
      level.sol.axpy( 1.0, level.tmp );
      level.matrix->residual( level.tmp, level.rhs, level.rhs );
    }
  };

  bool const trackRnorm = m_params.logLevel >= 4;

  m_levels[0].rhs.copy( src );
  int const numLevels = LvArray::integerConversion< int >( m_levels.size() );

  // down phase
  {
    GEOS_MARK_SCOPE( v-cycle down phase );
    for( int levelIndex = 0; levelIndex < numLevels - 1; ++levelIndex )
    {
      CALI_CXX_MARK_SCOPE( GEOS_FMT( "Level {}", levelIndex ).c_str() );
      LevelData const & fine = m_levels[levelIndex];
      LevelData const & coarse = m_levels[levelIndex + 1];

      fine.sol.zero();
      NormTracker< Vector > tracker( fine.rhs, " pre-smoothing", levelIndex, trackRnorm );
      applySmoother( fine, fine.builder->presmoother() );

      {
        GEOS_MARK_SCOPE( restriction );
        coarse.builder->restriction().apply( fine.rhs, coarse.rhs );
      }
    }
  }

  // coarse level solve
  {
    GEOS_MARK_SCOPE( coarse solve );
    LevelData const & coarse = m_levels.back();
    NormTracker< Vector > tracker( coarse.rhs, "  coarse solve", m_levels.size() - 1, trackRnorm );
    m_coarse_solver->apply( coarse.rhs, coarse.sol );
    if( trackRnorm )
    {
      coarse.matrix->residual( coarse.sol, coarse.rhs, coarse.rhs );
    }
  }

  // up phase
  {
    GEOS_MARK_SCOPE( v-cycle up phase );
    for( int levelIndex = numLevels - 2; levelIndex >= 0; --levelIndex )
    {
      CALI_CXX_MARK_SCOPE( GEOS_FMT( "Level {}", levelIndex ).c_str() );
      LevelData const & fine = m_levels[levelIndex];
      LevelData const & coarse = m_levels[levelIndex + 1];

      {
        GEOS_MARK_SCOPE( prolongation );
        coarse.builder->prolongation().apply( coarse.sol, fine.tmp );
      }

      {
        GEOS_MARK_SCOPE( residual );
        fine.sol.axpy( 1.0, fine.tmp );
        fine.matrix->residual( fine.tmp, fine.rhs, fine.rhs );
      }

      NormTracker< Vector > tracker( fine.rhs, "post-smoothing", levelIndex, trackRnorm );
      applySmoother( fine, fine.builder->postsmoother() );
    }
  }

  dst.copy( m_levels[0].sol );
}

template< typename LAI >
void MultiscalePreconditioner< LAI >::clear()
{
  m_levels.clear();
  m_coarse_solver.reset();
  m_initialized = false;
  Base::clear();
}

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOSX_USE_TRILINOS
template class MultiscalePreconditioner< TrilinosInterface >;
#endif

#ifdef GEOSX_USE_HYPRE
template class MultiscalePreconditioner< HypreInterface >;
#endif

#ifdef GEOSX_USE_PETSC
template class MultiscalePreconditioner< PetscInterface >;
#endif

} // namespace geos
