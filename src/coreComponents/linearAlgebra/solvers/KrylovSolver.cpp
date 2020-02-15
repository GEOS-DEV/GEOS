/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#include "KrylovSolver.hpp"

#include "linearAlgebra/interfaces/InterfaceTypes.hpp"

namespace geosx
{

template< typename VECTOR >
KrylovSolver< VECTOR >::KrylovSolver( LinearOperator<Vector> const & A,
                                      LinearOperator<Vector> const & M,
                                      real64 const tolerance,
                                      localIndex const maxIterations,
                                      integer const verbosity )
: m_operator( A ),
  m_precond( M ),
  m_tolerance( tolerance ),
  m_maxIterations( maxIterations ),
  m_verbosity( verbosity )
{

}

template< typename VECTOR >
KrylovSolver< VECTOR >::~KrylovSolver() = default;

template< typename VECTOR >
void KrylovSolver< VECTOR >::multiply( Vector const & src,
                                       Vector & dst ) const
{
  solve( src, dst );
}

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOSX_USE_TRILINOS
template class KrylovSolver<TrilinosInterface::ParallelVector>;
template class KrylovSolver<BlockVectorView<TrilinosInterface::ParallelVector>>;
#endif

#ifdef GEOSX_USE_HYPRE
template class KrylovSolver<HypreInterface::ParallelVector>;
template class KrylovSolver<BlockVectorView<HypreInterface::ParallelVector>>;
#endif

#ifdef GEOSX_USE_PETSC
template class KrylovSolver<PetscInterface::ParallelVector>;
template class KrylovSolver<BlockVectorView<PetscInterface::ParallelVector>>;
#endif

} //namespace geosx
