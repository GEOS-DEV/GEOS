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
 * @file PreconditionerAMG.cpp
 */

#include "TrilinosPreconditioner.hpp"

#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"

#include <Epetra_Operator.h>
#include <Epetra_FEVector.h>
#include <Epetra_FECrsMatrix.h>
#include <ml_MultiLevelPreconditioner.h>
#include <Ifpack.h>

namespace geosx
{

TrilinosPreconditioner::TrilinosPreconditioner( LinearSolverParameters params )
  : Base{},
  m_parameters( std::move( params ) ),
  m_precond{}
{}

TrilinosPreconditioner::~TrilinosPreconditioner() = default;


namespace
{

string const & getMLCycleType( string const & value )
{
  static std::map< string, string > const OPTION_MAP =
  {
    { "V", "MGV" },
    { "W", "MGW" },
  };

  GEOSX_LAI_ASSERT_MSG( OPTION_MAP.count( value ) > 0, "Unsupported Trilinos/ML cycle option: " << value );
  return OPTION_MAP.at( value );
}

string const & getMLSmootherType( string const & value )
{
  static std::map< string, string > const OPTION_MAP =
  {
    { "jacobi", "Jacobi" },
    { "blockJacobi", "block Jacobi" },
    { "gaussSeidel", "Gauss-Seidel" },
    { "blockGaussSeidel", "block Gauss-Seidel" },
    { "chebyshev", "Chebyshev" },
    { "icc", "IC" },
    { "ilu", "ILU" },
    { "ilut", "ILUT" },
  };

  GEOSX_LAI_ASSERT_MSG( OPTION_MAP.count( value ) > 0, "Unsupported Trilinos/ML smoother option: " << value );
  return OPTION_MAP.at( value );
}

string const & getMLCoarseType( string const & value )
{
  static std::map< string, string > const OPTION_MAP =
  {
    { "jacobi", "Jacobi" },
    { "gaussSeidel", "Gauss-Seidel" },
    { "blockGaussSeidel", "block Gauss-Seidel" },
    { "chebyshev", "Chebyshev" },
    { "direct", "Amesos-KLU"},
  };

  GEOSX_LAI_ASSERT_MSG( OPTION_MAP.count( value ) > 0, "Unsupported Trilinos/ML coarse solver option: " << value );
  return OPTION_MAP.at( value );
}

std::unique_ptr< Epetra_Operator >
CreateMLOperator( LinearSolverParameters const & params, Epetra_RowMatrix const & matrix )
{
  Teuchos::ParameterList list;
  GEOSX_LAI_CHECK_ERROR( ML_Epetra::SetDefaults( params.isSymmetric ? "SA" : "NSSA", list ) );

  list.set( "ML output", params.logLevel );
  list.set( "max levels", params.amg.maxLevels );
  list.set( "aggregation: type", "Uncoupled" );
  list.set( "aggregation: threshold", params.amg.threshold );
  list.set( "PDE equations", params.dofsPerNode );
  list.set( "smoother: sweeps", params.amg.numSweeps );
  list.set( "prec type", getMLCycleType( params.amg.cycleType ) );
  list.set( "smoother: type", getMLSmootherType( params.amg.smootherType ) );
  list.set( "smoother: pre or post", params.amg.preOrPostSmoothing );
  list.set( "coarse: type", getMLCoarseType( params.amg.coarseType ) );

  std::unique_ptr< Epetra_Operator > precond =
    std::make_unique< ML_Epetra::MultiLevelPreconditioner >( matrix, list );

  return precond;
}

Ifpack::EPrecType getIfpackPrecondType( LinearSolverParameters::PreconditionerType const & type )
{
  static std::map< LinearSolverParameters::PreconditionerType, Ifpack::EPrecType > const TYPE_MAP =
  {
    { LinearSolverParameters::PreconditionerType::iluk, Ifpack::ILU },
    { LinearSolverParameters::PreconditionerType::ilut, Ifpack::ILUT },
    { LinearSolverParameters::PreconditionerType::icc, Ifpack::IC },
    { LinearSolverParameters::PreconditionerType::ict, Ifpack::ICT },
    { LinearSolverParameters::PreconditionerType::jacobi, Ifpack::POINT_RELAXATION },
    { LinearSolverParameters::PreconditionerType::gs, Ifpack::POINT_RELAXATION },
    { LinearSolverParameters::PreconditionerType::sgs, Ifpack::POINT_RELAXATION },
  };

  GEOSX_LAI_ASSERT_MSG( TYPE_MAP.count( type ) > 0, "Unsupported Trilinos/Ifpack preconditioner option: " << type );
  return TYPE_MAP.at( type );
}
string getIfpackRelaxationType( LinearSolverParameters::PreconditionerType const & type )
{
  static std::map< LinearSolverParameters::PreconditionerType, string > const TYPE_MAP =
  {
    { LinearSolverParameters::PreconditionerType::jacobi, "Jacobi" },
    { LinearSolverParameters::PreconditionerType::gs, "Gauss-Seidel" },
    { LinearSolverParameters::PreconditionerType::sgs, "symmetric Gauss-Seidel" },
  };

  GEOSX_LAI_ASSERT_MSG( TYPE_MAP.count( type ) > 0, "Unsupported Trilinos/Ifpack preconditioner option: " << type );
  return TYPE_MAP.at( type );
}

std::unique_ptr< Epetra_Operator >
CreateIfpackOperator( LinearSolverParameters const & params, Epetra_RowMatrix const & matrix )
{
  Ifpack::EPrecType const type = getIfpackPrecondType( params.preconditionerType );
  std::unique_ptr< Ifpack_Preconditioner > precond( Ifpack::Create( type,
                                                                    const_cast< Epetra_RowMatrix * >( &matrix ),
                                                                    params.dd.overlap ) );
  Teuchos::ParameterList list;
  if( type == Ifpack::POINT_RELAXATION )
  {
    list.set( "relaxation: type", getIfpackRelaxationType( params.preconditionerType ) );
  }
  else
  {
    list.set( "fact: level-of-fill", params.ilu.fill );
    list.set( "fact: ilut level-of-fill", std::max( real64( params.ilu.fill ), 1.0 ) );
    list.set( "fact: ict level-of-fill", std::max( real64( params.ilu.fill ), 1.0 ) );
  }

  GEOSX_LAI_CHECK_ERROR( precond->SetParameters( list ) );
  GEOSX_LAI_CHECK_ERROR( precond->Compute() );

  return precond;
}

} // namespace


void TrilinosPreconditioner::compute( Matrix const & mat )
{
  Base::compute( mat );

  // To be able to use Trilinos solvers we need to disable floating point exceptions
  // In particular, we want to avoid a ML FPE crashes on Mac
  LvArray::system::FloatingPointExceptionGuard guard;

  switch( m_parameters.preconditionerType )
  {
    case LinearSolverParameters::PreconditionerType::amg:
    {
      m_precond = CreateMLOperator( m_parameters, mat.unwrapped() );
      break;
    }
    case LinearSolverParameters::PreconditionerType::jacobi:
    case LinearSolverParameters::PreconditionerType::gs:
    case LinearSolverParameters::PreconditionerType::sgs:
    case LinearSolverParameters::PreconditionerType::iluk:
    case LinearSolverParameters::PreconditionerType::ilut:
    case LinearSolverParameters::PreconditionerType::icc:
    case LinearSolverParameters::PreconditionerType::ict:
    {
      m_precond = CreateIfpackOperator( m_parameters, mat.unwrapped() );
      break;
    }
    case LinearSolverParameters::PreconditionerType::none:
    {
      m_precond.reset();
      break;
    }
    default:
    {
      GEOSX_ERROR( "Preconditioner type not supported in Trilinos interface: " << m_parameters.preconditionerType );
    }
  }
}

void TrilinosPreconditioner::apply( Vector const & src,
                                    Vector & dst ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( src.ready() );
  GEOSX_LAI_ASSERT( dst.ready() );
  GEOSX_LAI_ASSERT_EQ( src.globalSize(), this->numGlobalCols() );
  GEOSX_LAI_ASSERT_EQ( dst.globalSize(), this->numGlobalRows() );

  if( !m_precond )
  {
    dst.copy( src );
  }
  else
  {
    GEOSX_LAI_ASSERT( m_precond );
    m_precond->ApplyInverse( src.unwrapped(), dst.unwrapped() );
  }
}

void TrilinosPreconditioner::clear()
{
  PreconditionerBase::clear();
  m_precond.reset();
}

Epetra_Operator const & TrilinosPreconditioner::unwrapped() const
{
  GEOSX_LAI_ASSERT( ready() );
  return *m_precond;
}

Epetra_Operator & TrilinosPreconditioner::unwrapped()
{
  GEOSX_LAI_ASSERT( ready() );
  return *m_precond;
}

}
