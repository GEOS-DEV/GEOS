/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file KrylovSolver.cpp
 */

#include "KrylovSolver.hpp"
#include "linearAlgebra/solvers/BicgstabSolver.hpp"
#include "linearAlgebra/solvers/CgSolver.hpp"
#include "linearAlgebra/solvers/GmresSolver.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"

namespace geosx
{

template< typename VECTOR >
KrylovSolver< VECTOR >::KrylovSolver( LinearSolverParameters params,
                                      LinearOperator< Vector > const & matrix,
                                      LinearOperator< Vector > const & precond )
  : Base(),
  m_params( std::move( params ) ),
  m_operator( matrix ),
  m_precond( precond )
{
  GEOSX_ERROR_IF_LE_MSG( m_params.krylov.maxIterations, 0, "Krylov solver: max number of iteration must be positive." );
  GEOSX_LAI_ASSERT_EQ( m_operator.numLocalRows(), m_operator.numLocalCols() );
  GEOSX_LAI_ASSERT_EQ( m_operator.numLocalRows(), m_precond.numLocalRows() );
  GEOSX_LAI_ASSERT_EQ( m_operator.numLocalCols(), m_precond.numLocalCols() );
}

template< typename VECTOR >
std::unique_ptr< KrylovSolver< VECTOR > >
KrylovSolver< VECTOR >::create( LinearSolverParameters const & parameters,
                                LinearOperator< VECTOR > const & matrix,
                                LinearOperator< VECTOR > const & precond )
{
  switch( parameters.solverType )
  {
    case LinearSolverParameters::SolverType::cg:
    {
      return std::make_unique< CgSolver< Vector > >( parameters,
                                                     matrix,
                                                     precond );
    }
    case LinearSolverParameters::SolverType::bicgstab:
    {
      return std::make_unique< BicgstabSolver< Vector > >( parameters,
                                                           matrix,
                                                           precond );
    }
    case LinearSolverParameters::SolverType::gmres:
    {
      return std::make_unique< GmresSolver< Vector > >( parameters,
                                                        matrix,
                                                        precond );
    }
    default:
    {
      GEOSX_ERROR( "Unsupported linear solver type: " << parameters.solverType );
    }
  }
  return {};
}

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOSX_USE_TRILINOS
template class KrylovSolver< TrilinosInterface::ParallelVector >;
template class KrylovSolver< BlockVectorView< TrilinosInterface::ParallelVector > >;
#endif

#ifdef GEOSX_USE_HYPRE
template class KrylovSolver< HypreInterface::ParallelVector >;
template class KrylovSolver< BlockVectorView< HypreInterface::ParallelVector > >;
#endif

#ifdef GEOSX_USE_PETSC
template class KrylovSolver< PetscInterface::ParallelVector >;
template class KrylovSolver< BlockVectorView< PetscInterface::ParallelVector > >;
#endif

}
