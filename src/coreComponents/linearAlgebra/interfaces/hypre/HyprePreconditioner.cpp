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
 * @file HyprePreconditioner.cpp
 */

#include "HyprePreconditioner.hpp"

#include "linearAlgebra/interfaces/hypre/HypreUtils.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"

#include <_hypre_utilities.h>
#include <_hypre_parcsr_ls.h>
#include <_hypre_IJ_mv.h>
#include <krylov.h>

namespace geosx
{

HyprePreconditioner::HyprePreconditioner( LinearSolverParameters params )
  : Base{},
  m_parameters( std::move( params ) ),
  m_precond{},
  m_functions( std::make_unique< HyprePrecFuncs >() )
{
  if( m_precond == nullptr )
  {
    if( m_parameters.preconditionerType == "none" )
    {
      m_functions->setup = (HYPRE_PtrToParSolverFcn) hypre_ParKrylovIdentitySetup;
      m_functions->apply = (HYPRE_PtrToParSolverFcn) hypre_ParKrylovIdentity;
    }
    else if( m_parameters.preconditionerType == "jacobi" )
    {
      m_functions->setup = (HYPRE_PtrToParSolverFcn) HYPRE_ParCSRDiagScaleSetup;
      m_functions->apply = (HYPRE_PtrToParSolverFcn) HYPRE_ParCSRDiagScale;
    }
    if( m_parameters.preconditionerType == "amg" )
    {
      createAMG();
    }
    else if( m_parameters.preconditionerType == "iluk" )
    {
      createILU();
    }
    else if( m_parameters.preconditionerType == "ilut" )
    {
      createILUT();
    }
    else
    {
      GEOSX_ERROR( "Unsupported preconditioner type: " << m_parameters.preconditionerType );
    }
  }
}

HyprePreconditioner::~HyprePreconditioner()
{
  clear();
}

namespace
{

HYPRE_Int getHypreAMGCycleType( string const & type )
{
  static std::map< string, HYPRE_Int > const typeMap =
  {
    { "V", 1 },
    { "W", 2 },
  };

  GEOSX_LAI_ASSERT_MSG( typeMap.count( type ) > 0, "Unsupported Hypre AMG cycle option: " << type );
  return typeMap.at( type );
}

HYPRE_Int getHypreAMGRelaxationType( string const & type )
{
  static std::map< string, HYPRE_Int > const typeMap =
  {
    { "jacobi", 0 },
    { "hybridForwardGaussSeidel", 3 },
    { "hybridBackwardGaussSeidel", 4 },
    { "hybridSymmetricGaussSeidel", 6 },
    { "gaussSeidel", 6 },
    { "L1hybridSymmetricGaussSeidel", 8 },
    { "chebyshev", 16 },
    { "L1jacobi", 18 },
  };

  GEOSX_LAI_ASSERT_MSG( typeMap.count( type ) > 0, "Unsupported Hypre AMG relaxation option: " << type );
  return typeMap.at( type );
}

}

void HyprePreconditioner::createAMG()
{
  GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGCreate( &m_precond ) );


  // Hypre's parameters to use BoomerAMG as a preconditioner
  GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetTol( m_precond, 0.0 ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxIter( m_precond, 1 ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetPrintLevel( m_precond, toHYPRE_Int( m_parameters.logLevel ) ) );;

  // Set maximum number of multigrid levels (default 25)
  GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxLevels( m_precond, toHYPRE_Int( m_parameters.amg.maxLevels ) ) );

  // Set type of cycle (1: V-cycle (default); 2: W-cycle)
  GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCycleType( m_precond, getHypreAMGCycleType( m_parameters.amg.cycleType ) ) );

  // Set smoother to be used (other options available, see hypre's documentation)
  // (default "gaussSeidel", i.e. local symmetric Gauss-Seidel)
  if( m_parameters.amg.smootherType.substr( 0, 3 ) == "ilu" )
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetSmoothType( m_precond, 5 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetType( m_precond, 0 ) );
    if( m_parameters.amg.smootherType == "ilu1" )
    {
      GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetLevelOfFill( m_precond, 1 ) );
    }
  }
  else
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxType( m_precond, getHypreAMGRelaxationType( m_parameters.amg.smootherType ) ) );

    if( m_parameters.amg.smootherType == "chebyshev" )
    {
      // Set order for Chebyshev smoother valid options 1, 2 (default), 3, 4)
      if( ( 0 < m_parameters.amg.numSweeps ) && ( m_parameters.amg.numSweeps < 5 ) )
      {
        GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetChebyOrder( m_precond, toHYPRE_Int( m_parameters.amg.numSweeps ) ) );
      }
    }
  }

  // Set coarsest level solver
  // (by default for coarsest grid size above 5,000 superlu_dist is used)
  GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetDSLUThreshold( m_precond, 5000 ) );
  if( m_parameters.amg.coarseType == "direct" )
  {
    GEOSX_LAI_CHECK_ERROR( hypre_BoomerAMGSetCycleRelaxType( m_precond, 9, 3 ) );
  }

  // Set the number of sweeps
  if( m_parameters.amg.preOrPostSmoothing == "both" )
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumSweeps( m_precond, toHYPRE_Int( m_parameters.amg.numSweeps ) ) );
  }
  else if( m_parameters.amg.preOrPostSmoothing == "pre" )
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCycleNumSweeps( m_precond, toHYPRE_Int( m_parameters.amg.numSweeps ), 1 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCycleNumSweeps( m_precond, 0, 2 ) );
  }
  else if( m_parameters.amg.preOrPostSmoothing == "post" )
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCycleNumSweeps( m_precond, 0, 1 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCycleNumSweeps( m_precond, toHYPRE_Int( m_parameters.amg.numSweeps ), 2 ) );
  }

  // Set strength of connection
  if( m_parameters.amg.threshold > 0.0 )
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetStrongThreshold( m_precond, m_parameters.amg.threshold ) );
  }

  m_functions->setup = HYPRE_BoomerAMGSetup;
  m_functions->apply = HYPRE_BoomerAMGSolve;
  m_functions->destroy = HYPRE_BoomerAMGDestroy;
}

void HyprePreconditioner::createILU()
{
  GEOSX_LAI_CHECK_ERROR( HYPRE_ILUCreate( &m_precond ) );

  // Hypre's parameters to use ParCSR ILU as a preconditioner
  GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetMaxIter( m_precond, 1 ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetTol( m_precond, 0.0 ) );

  if( m_parameters.ilu.fill >= 0 )
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetLevelOfFill( m_precond, toHYPRE_Int( m_parameters.ilu.fill ) ) );
  }

  m_functions->setup = HYPRE_ILUSetup;
  m_functions->apply = HYPRE_ILUSolve;
  m_functions->destroy = HYPRE_ILUDestroy;
}

void HyprePreconditioner::createILUT()
{
  GEOSX_LAI_CHECK_ERROR( HYPRE_ILUCreate( &m_precond ) );

  // Hypre's parameters to use ParCSR ILU as a preconditioner
  GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetMaxIter( m_precond, 1 ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetTol( m_precond, 0.0 ) );

  if( m_parameters.ilu.fill >= 0 )
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetLevelOfFill( m_precond, toHYPRE_Int( m_parameters.ilu.fill ) ) );
  }
  if( m_parameters.ilu.threshold >= 0 )
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetDropThreshold( m_precond,
                                                      m_parameters.ilu.threshold ) );
  }

  m_functions->setup = HYPRE_ILUSetup;
  m_functions->apply = HYPRE_ILUSolve;
  m_functions->destroy = HYPRE_ILUDestroy;
}

void HyprePreconditioner::compute( Matrix const & mat )
{
  PreconditionerBase::compute( mat );
  GEOSX_LAI_CHECK_ERROR( m_functions->setup( m_precond, mat.unwrapped(), nullptr, nullptr ) );
}

void HyprePreconditioner::apply( Vector const & src,
                                 Vector & dst ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( src.ready() );
  GEOSX_LAI_ASSERT( dst.ready() );
  GEOSX_LAI_ASSERT_EQ( src.globalSize(), this->numGlobalCols() );
  GEOSX_LAI_ASSERT_EQ( dst.globalSize(), this->numGlobalRows() );

  m_functions->apply( m_precond, this->matrix().unwrapped(), src.unwrapped(), dst.unwrapped() );
}

void HyprePreconditioner::clear()
{
  PreconditionerBase::clear();
  if( m_precond != nullptr && m_functions && m_functions->destroy != nullptr )
  {
    m_functions->destroy( m_precond );
  }
  m_functions.reset();
}

HYPRE_Solver const & HyprePreconditioner::unwrapped() const
{
  return m_precond;
}

HyprePrecFuncs const & HyprePreconditioner::unwrappedFuncs() const
{
  GEOSX_LAI_ASSERT( m_functions );
  return *m_functions;
}

}
