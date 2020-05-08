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
 * @file TrilinosTpetraPreconditioner.cpp
 */

#include "TrilinosTpetraPreconditioner.hpp"

#include <Tpetra_Operator.hpp>
#include <Ifpack2_Factory.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>

namespace geosx
{

TrilinosTpetraPreconditioner::TrilinosTpetraPreconditioner( LinearSolverParameters params )
  : Base{},
  m_parameters( std::move( params ) ),
  m_precond{}
{ }

TrilinosTpetraPreconditioner::~TrilinosTpetraPreconditioner() = default;

namespace
{

string const & getMueLuCoarseType( LinearSolverParameters::PreconditionerType const & value )
{
  static std::map< LinearSolverParameters::PreconditionerType, string > const optionMap =
  {
    { LinearSolverParameters::PreconditionerType::direct, "KLU2" },
  };

  GEOSX_LAI_ASSERT_MSG( optionMap.count( value ) > 0, "Unsupported Trilinos/MueLu coarse solver option: " << value );
  return optionMap.at( value );
}

string const & getMueLuVerbosity( integer const value )
{
  static std::array< string, 5 > const optionMap =
  {
    "none",
    "low",
    "medium",
    "high",
    "extreme"
  };

  return optionMap[ std::max( std::min( value, integer( optionMap.size() - 1 ) ), 0 ) ];
}

string getIfpack2PrecondType( LinearSolverParameters::PreconditionerType const & type )
{
  static std::map< LinearSolverParameters::PreconditionerType, string > const typeMap =
  {
    { LinearSolverParameters::PreconditionerType::iluk, "RILUK" },
    { LinearSolverParameters::PreconditionerType::ilut, "ILUT" },
    { LinearSolverParameters::PreconditionerType::jacobi, "RELAXATION" },
    { LinearSolverParameters::PreconditionerType::gs, "RELAXATION" },
    { LinearSolverParameters::PreconditionerType::sgs, "RELAXATION" },
    { LinearSolverParameters::PreconditionerType::direct, "AMESOS2" }
  };

  GEOSX_LAI_ASSERT_MSG( typeMap.count( type ) > 0, "Unsupported Trilinos/Ifpack2 preconditioner option: " << type );
  return typeMap.at( type );
}

string getIfpack2RelaxationType( LinearSolverParameters::PreconditionerType const & type )
{
  static std::map< LinearSolverParameters::PreconditionerType, string > const typeMap =
  {
    { LinearSolverParameters::PreconditionerType::jacobi, "Jacobi" },
    { LinearSolverParameters::PreconditionerType::gs, "Gauss-Seidel" },
    { LinearSolverParameters::PreconditionerType::sgs, "Symmetric Gauss-Seidel" },
  };

  GEOSX_LAI_ASSERT_MSG( typeMap.count( type ) > 0, "Unsupported Trilinos/Ifpack2 preconditioner option: " << type );
  return typeMap.at( type );
}

std::unique_ptr< TrilinosTpetraPreconditioner::Tpetra_Operator >
CreateMueLuOperator( LinearSolverParameters const & params,
                     TrilinosTpetraPreconditioner::Tpetra_RowMatrix const & matrix )
{
  using Tpetra_Operator  = TrilinosTpetraPreconditioner::Tpetra_Operator;
  using Tpetra_RowMatrix = TrilinosTpetraPreconditioner::Tpetra_RowMatrix;

  Teuchos::ParameterList list;
  list.set( "verbosity", getMueLuVerbosity( params.logLevel ) );
  list.set( "problem: symmetric", params.isSymmetric );
  list.set( "max levels", params.amg.maxLevels );
  list.set( "cycle type", params.amg.cycleType );
  list.set( "number of equations", params.dofsPerNode );
  list.set( "smoother: pre or post", params.amg.preOrPostSmoothing );
  list.set( "smoother: type", getIfpack2PrecondType( params.amg.smootherType ) );
  list.set( "coarse: type", getMueLuCoarseType( params.amg.coarseType ) );
  list.set( "smoother: overlap", params.dd.overlap );
  list.set( "aggregation: drop tol", params.amg.threshold );

  Teuchos::ParameterList smootherList;
  smootherList.set( "relaxation: type", getIfpack2RelaxationType( params.amg.smootherType ) );
  smootherList.set( "relaxation: sweeps", params.amg.numSweeps );
  smootherList.set( "relaxation: use l1", true );
  list.set( "smoother: params", smootherList );

  // Const cast unfortunately required, see e.g. https://github.com/trilinos/Trilinos/issues/6390
  // TODO: upon upgrading Trilinos release, check if this has been fixed
  Teuchos::RCP< Tpetra_Operator > matRCP( const_cast< Tpetra_RowMatrix * >( &matrix ), false );

  std::unique_ptr< Tpetra_Operator > precond( MueLu::CreateTpetraPreconditioner( matRCP, list ).release().get() );

  return precond;
}

std::unique_ptr< TrilinosTpetraPreconditioner::Tpetra_Operator >
CreateIfpackOperator( LinearSolverParameters const & params,
                      TrilinosTpetraPreconditioner::Tpetra_RowMatrix const & matrix )
{
  using Tpetra_Operator = TrilinosTpetraPreconditioner::Tpetra_Operator;
  using Ifpack2_Preconditioner = Ifpack2::Preconditioner< Tpetra_Operator::scalar_type,
                                                          Tpetra_Operator::local_ordinal_type,
                                                          Tpetra_Operator::global_ordinal_type,
                                                          Tpetra_Operator::node_type >;

  Teuchos::RCP< Ifpack2_Preconditioner > precond = Ifpack2::Factory::create( "SCHWARZ", Teuchos::rcp( &matrix, false ) );

  string const innerPrecType = getIfpack2PrecondType( params.preconditionerType );

  Teuchos::ParameterList innerParams;

  if( innerPrecType == "RELAXATION" )
  {
    innerParams.set( "relaxation: type", getIfpack2RelaxationType( params.preconditionerType ) );
  }
  else if( innerPrecType == "RILUK" )
  {
    innerParams.set( "fact: iluk level-of-fill", params.ilu.fill );
    innerParams.set( "fact: drop tolerance", params.ilu.threshold );
  }
  else if( innerPrecType == "ILUT" )
  {
    innerParams.set( "fact: ilut level-of-fill", std::fmax( params.ilu.fill, 1.0 ) );
    innerParams.set( "fact: drop tolerance", params.ilu.threshold );
  }

  Teuchos::ParameterList outerParams;
  outerParams.set( "schwarz: overlap level", params.dd.overlap );
  outerParams.set( "schwarz: subdomain solver name", innerPrecType );
  outerParams.set( "schwarz: subdomain solver parameters", innerParams );
  precond->setParameters( outerParams );

  precond->compute();

  return std::unique_ptr< Tpetra_Operator >( precond.release().get() );
}

std::unique_ptr< TrilinosTpetraPreconditioner::Tpetra_Operator >
CreateIfpackDirect( LinearSolverParameters const & GEOSX_UNUSED_PARAM( params ),
                    TrilinosTpetraPreconditioner::Tpetra_RowMatrix const & matrix )
{
  using Tpetra_Operator = TrilinosTpetraPreconditioner::Tpetra_Operator;
  using Ifpack2_Preconditioner = Ifpack2::Preconditioner< Tpetra_Operator::scalar_type,
                                                          Tpetra_Operator::local_ordinal_type,
                                                          Tpetra_Operator::global_ordinal_type,
                                                          Tpetra_Operator::node_type >;

  Teuchos::RCP< Ifpack2_Preconditioner > precond = Ifpack2::Factory::create( "AMESOS2", Teuchos::rcp( &matrix, false ) );
  precond->compute();

  return std::unique_ptr< Tpetra_Operator >( precond.release().get() );
}

} // namespace


void TrilinosTpetraPreconditioner::compute( Matrix const & mat )
{
  Base::compute( mat );

  switch( m_parameters.preconditionerType )
  {
    case LinearSolverParameters::PreconditionerType::amg:
    {
      m_precond = CreateMueLuOperator( m_parameters, mat.unwrapped() );
      break;
    }
    case LinearSolverParameters::PreconditionerType::iluk:
    case LinearSolverParameters::PreconditionerType::ilut:
    {
      m_precond = CreateIfpackOperator( m_parameters, mat.unwrapped() );
      break;
    }
    case LinearSolverParameters::PreconditionerType::direct:
    {
      m_precond = CreateIfpackDirect( m_parameters, mat.unwrapped() );
      break;
    }
    case LinearSolverParameters::PreconditionerType::none:
    {
      m_precond.reset();
      break;
    }
    default:
    {
      GEOSX_ERROR( "Unsupported preconditioner type: " << m_parameters.preconditionerType );
    }
  }
}

void TrilinosTpetraPreconditioner::apply( Vector const & src,
                                          Vector & dst ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( src.ready() );
  GEOSX_LAI_ASSERT( dst.ready() );
  GEOSX_LAI_ASSERT_EQ( src.globalSize(), this->numGlobalCols() );
  GEOSX_LAI_ASSERT_EQ( dst.globalSize(), this->numGlobalRows() );

  if( m_parameters.preconditionerType == LinearSolverParameters::PreconditionerType::none )
  {
    dst.copy( src );
  }
  else
  {
    GEOSX_LAI_ASSERT( m_precond );
    m_precond->apply( src.unwrapped(), dst.unwrapped() );
  }
}

void TrilinosTpetraPreconditioner::clear()
{
  PreconditionerBase::clear();
  m_precond.reset();
}

TrilinosTpetraPreconditioner::Tpetra_Operator const & TrilinosTpetraPreconditioner::unwrapped() const
{
  GEOSX_LAI_ASSERT( ready() );
  return *m_precond;
}

TrilinosTpetraPreconditioner::Tpetra_Operator & TrilinosTpetraPreconditioner::unwrapped()
{
  GEOSX_LAI_ASSERT( ready() );
  return *m_precond;
}

} // namespace geosx
