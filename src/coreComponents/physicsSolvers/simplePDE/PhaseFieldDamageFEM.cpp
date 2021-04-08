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
 * @file PhaseFieldDamageFEM.hpp
 */

#include "PhaseFieldDamageFEM.hpp"
#include "PhaseFieldDamageFEMKernels.hpp"
#include "PoroPhaseFieldDamageFEMKernels.hpp"
#include <math.h>
#include <vector>

#include "common/TimingMacros.hpp"
#include "dataRepository/Group.hpp"
#include "mpiCommunications/CommunicationTools.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"

#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "constitutive/ConstitutiveBase.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "constitutive/solid/Damage.hpp"
#include "constitutive/solid/SolidBase.hpp"
#include "finiteElement/FiniteElementDiscretization.hpp"
#include "finiteElement/FiniteElementDiscretizationManager.hpp"
#include "finiteElement/Kinematics.h"
#include "managers/NumericalMethodsManager.hpp"

#include "managers/DomainPartition.hpp"

namespace geosx
{

namespace dataRepository
{
namespace keys
{}
} // namespace dataRepository

using namespace dataRepository;
using namespace constitutive;

PhaseFieldDamageFEM::PhaseFieldDamageFEM( const string & name,
                                          Group * const parent ):
  SolverBase( name, parent ),
  m_fieldName( "primaryField" ),
  m_solidModelNames()
{

  registerWrapper< string >( PhaseFieldDamageFEMViewKeys.timeIntegrationOption.key() ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "option for default time integration method" );

  registerWrapper< string >( PhaseFieldDamageFEMViewKeys.fieldVarName.key(), &m_fieldName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "name of field variable" );

  registerWrapper( viewKeyStruct::localDissipationOptionString(), &m_localDissipationOption ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Type of local dissipation function. Can be Linear or Quadratic" );

  registerWrapper( viewKeyStruct::solidModelNamesString(), &m_solidModelNames ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "name of solid constitutive model" );
}

PhaseFieldDamageFEM::~PhaseFieldDamageFEM()
{
  // TODO Auto-generated destructor stub
}

void PhaseFieldDamageFEM::registerDataOnMesh( Group & meshBodies )
{
  meshBodies.forSubGroups< MeshBody >( [&] ( MeshBody & meshBody )
  {
    MeshLevel & meshLevel = meshBody.getMeshLevel( 0 );

    NodeManager & nodes = meshLevel.getNodeManager();

    nodes.registerWrapper< real64_array >( m_fieldName )
      .setApplyDefaultValue( 0.0 )
      .setPlotLevel( PlotLevel::LEVEL_0 )
      .setDescription( "Primary field variable" );

    ElementRegionManager & elemManager = meshLevel.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion >( [ &]( CellElementSubRegion & subRegion )
    {
      subRegion.registerWrapper( viewKeyStruct::coeffNameString(), &m_coeff ).
        setApplyDefaultValue( 0.0 ).
        setPlotLevel( PlotLevel::LEVEL_0 ).
        setDescription( "field variable representing the diffusion coefficient" );
    } );
  } );
}

void PhaseFieldDamageFEM::postProcessInput()
{
  SolverBase::postProcessInput();

  string tiOption = this->getReference< string >(
    PhaseFieldDamageFEMViewKeys.timeIntegrationOption );

  if( tiOption == "SteadyState" )
  {
    this->m_timeIntegrationOption = timeIntegrationOption::SteadyState;
  }
  else if( tiOption == "ImplicitTransient" )
  {
    this->m_timeIntegrationOption = timeIntegrationOption::ImplicitTransient;
  }
  else if( tiOption == "ExplicitTransient" )
  {
    this->m_timeIntegrationOption = timeIntegrationOption::ExplicitTransient;
  }
  else
  {
    GEOSX_ERROR( "invalid time integration option" );
  }

  if( m_localDissipationOption != "Linear" && m_localDissipationOption != "Quadratic" )
  {
    GEOSX_ERROR( "invalid local dissipation option - must be Linear or Quadratic" );
  }

  // Set basic parameters for solver
  // m_linearSolverParameters.logLevel = 0;
  // m_linearSolverParameters.solverType = "gmres";
  // m_linearSolverParameters.krylov.tolerance = 1e-8;
  // m_linearSolverParameters.krylov.maxIterations = 250;
  // m_linearSolverParameters.krylov.maxRestart = 250;
  // m_linearSolverParameters.preconditionerType = "amg";
  // m_linearSolverParameters.amg.smootherType = "gaussSeidel";
  // m_linearSolverParameters.amg.coarseType = "direct";
}

real64 PhaseFieldDamageFEM::solverStep( real64 const & time_n,
                                        real64 const & dt,
                                        const int cycleNumber,
                                        DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  real64 dtReturn = dt;
  if( m_timeIntegrationOption == timeIntegrationOption::ExplicitTransient )
  {
    dtReturn = explicitStep( time_n, dt, cycleNumber, domain );
  }
  else if( m_timeIntegrationOption ==
           timeIntegrationOption::ImplicitTransient ||
           m_timeIntegrationOption == timeIntegrationOption::SteadyState )
  {
    this->setupSystem( domain, m_dofManager, m_localMatrix, m_localRhs, m_localSolution, false );

    dtReturn = this->nonlinearImplicitStep( time_n,
                                            dt,
                                            cycleNumber,
                                            domain );
  }
  return dtReturn;
}

real64 PhaseFieldDamageFEM::explicitStep(
  real64 const & GEOSX_UNUSED_PARAM( time_n ),
  real64 const & dt,
  const int GEOSX_UNUSED_PARAM( cycleNumber ),
  DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{
  return dt;
}

void PhaseFieldDamageFEM::setupSystem( DomainPartition & domain,
                                       DofManager & dofManager,
                                       CRSMatrix< real64, globalIndex > & localMatrix,
                                       array1d< real64 > & localRhs,
                                       array1d< real64 > & localSolution,
                                       bool const setSparsity )
{
  GEOSX_MARK_FUNCTION;
  SolverBase::setupSystem( domain, dofManager, localMatrix, localRhs, localSolution, setSparsity );
}

void PhaseFieldDamageFEM::implicitStepComplete(
  real64 const & GEOSX_UNUSED_PARAM( time_n ),
  real64 const & GEOSX_UNUSED_PARAM( dt ),
  DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{}

void PhaseFieldDamageFEM::setupDofs(
  DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
  DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;
  dofManager.addField( m_fieldName, DofManager::Location::Node );

  dofManager.addCoupling( m_fieldName,
                          m_fieldName,
                          DofManager::Connector::Elem );

}

void PhaseFieldDamageFEM::assembleSystem( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                          real64 const GEOSX_UNUSED_PARAM( dt ),
                                          DomainPartition & domain,
                                          DofManager const & dofManager,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  NodeManager & nodeManager = mesh.getNodeManager();

  arrayView1d< globalIndex const > const & dofIndex = nodeManager.getReference< array1d< globalIndex > >( dofManager.getKey( m_fieldName ) );

  // Initialize all entries to zero
#if 1 // Andre...this is the new code
  localMatrix.setValues< parallelDevicePolicy< 32 > >( 0 );
  localRhs.setValues< parallelDevicePolicy< 32 > >( 0 );
  if (m_pressureTerm == 0) {
    finiteElement::
      regionBasedKernelApplication< serialPolicy,
				    constitutive::DamageBase,
				    CellElementSubRegion,
				    PhaseFieldDamageKernel >( mesh,
                                                            targetRegionNames(),
                                                            this->getDiscretizationName(),
                                                            m_solidModelNames,
                                                            dofIndex,
                                                            dofManager.rankOffset(),
                                                            localMatrix,
                                                            localRhs,
                                                            m_fieldName,
                                                            m_localDissipationOption=="Linear" ? 1 : 2 );
  }
  else {
    finiteElement::
      regionBasedKernelApplication< serialPolicy,
				    constitutive::DamageBase,
				    CellElementSubRegion,
				    PoroPhaseFieldDamageKernel >( mesh,
                                                            targetRegionNames(),
                                                            this->getDiscretizationName(),
                                                            m_solidModelNames,
                                                            dofIndex,
                                                            dofManager.rankOffset(),
                                                            localMatrix,
                                                            localRhs,
                                                            m_fieldName,
                                                            m_localDissipationOption=="Linear" ? 1 : 2 );
    
  }
    
}

void PhaseFieldDamageFEM::applySystemSolution( DofManager const & dofManager,
                                               arrayView1d< real64 const > const & localSolution,
                                               real64 const scalingFactor,
                                               DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  dofManager.addVectorToField( localSolution, m_fieldName, m_fieldName, scalingFactor );

  // Syncronize ghost nodes
  std::map< string, string_array > fieldNames;
  fieldNames["node"].emplace_back( m_fieldName );

  getGlobalState().getCommunicationTools().synchronizeFields( fieldNames,
                                                              mesh,
                                                              domain.getNeighbors(),
                                                              false );
}

void PhaseFieldDamageFEM::applyBoundaryConditions(
  real64 const time_n,
  real64 const dt, DomainPartition & domain,
  DofManager const & dofManager,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;
  applyDirichletBCImplicit( time_n + dt, dofManager, domain, localMatrix, localRhs );

  if( getLogLevel() == 2 )
  {
    GEOSX_LOG_RANK_0( "After PhaseFieldDamageFEM::applyBoundaryConditions" );
    GEOSX_LOG_RANK_0( "\nJacobian:\n" );
    std::cout << localMatrix.toViewConst();
    GEOSX_LOG_RANK_0( "\nResidual:\n" );
    std::cout << localRhs;
  }
//
//  if( getLogLevel() >= 3 )
//  {
//    NonlinearSolverParameters & solverParams = getNonlinearSolverParameters();
//    integer newtonIter = solverParams.m_numNewtonIterations;
//
//    string filename_mat = "matrix_bc_" + std::to_string( time_n ) + "_" +
//                          std::to_string( newtonIter ) + ".mtx";
//    matrix.write( filename_mat );
//
//    string filename_rhs = "rhs_bc_" + std::to_string( time_n ) + "_" +
//                          std::to_string( newtonIter ) + ".mtx";
//    rhs.write( filename_rhs );
//
//    GEOSX_LOG_RANK_0( "After PhaseFieldDamageFEM::applyBoundaryConditions" );
//    GEOSX_LOG_RANK_0( "Jacobian: written to " << filename_mat );
//    GEOSX_LOG_RANK_0( "Residual: written to " << filename_rhs );
//  }
}

real64
PhaseFieldDamageFEM::calculateResidualNorm( DomainPartition const & domain,
                                            DofManager const & dofManager,
                                            arrayView1d< real64 const > const & localRhs )
{
  GEOSX_MARK_FUNCTION;
  const MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  const NodeManager & nodeManager = mesh.getNodeManager();
  const arrayView1d< const integer > & ghostRank = nodeManager.ghostRank();

  const arrayView1d< const globalIndex > &
  dofNumber = nodeManager.getReference< array1d< globalIndex > >( dofManager.getKey( m_fieldName ) );
  const globalIndex rankOffset = dofManager.rankOffset();

  RAJA::ReduceSum< parallelDeviceReduce, real64 > localSum( 0.0 );

  forAll< parallelDevicePolicy<> >( nodeManager.size(),
                                    [localRhs, localSum, dofNumber, rankOffset, ghostRank] GEOSX_HOST_DEVICE ( localIndex const k )
  {
    if( ghostRank[k] < 0 )
    {
      localIndex const localRow = LvArray::integerConversion< localIndex >( dofNumber[k] - rankOffset );
      localSum += localRhs[localRow] * localRhs[localRow];
    }
  } );

  const real64 localResidualNorm[2] =
  { localSum.get(), 1 };

  // globalResidualNorm[0]: the sum of all the local sum(rhs^2).
  // globalResidualNorm[1]: max of max force of each rank. Basically max force globally
  real64 globalResidualNorm[2] = {0, 0};

  const int rank = MpiWrapper::commRank( MPI_COMM_GEOSX );
  const int size = MpiWrapper::commSize( MPI_COMM_GEOSX );
  array1d< real64 > globalValues( size * 2 );

  // Everything is done on rank 0
  MpiWrapper::gather( localResidualNorm,
                      2,
                      globalValues.data(),
                      2,
                      0,
                      MPI_COMM_GEOSX );

  if( rank==0 )
  {
    for( int r=0; r<size; ++r )
    {
      // sum/max across all ranks
      globalResidualNorm[0] += globalValues[r*2];
      globalResidualNorm[1] = std::max( globalResidualNorm[1], globalValues[r*2+1] );
    }
  }

  MpiWrapper::bcast( globalResidualNorm, 2, 0, MPI_COMM_GEOSX );


  const real64 residual = sqrt( globalResidualNorm[0] ) / ( globalResidualNorm[1] );

  return residual;
}

void PhaseFieldDamageFEM::solveSystem( DofManager const & dofManager,
                                       ParallelMatrix & matrix,
                                       ParallelVector & rhs,
                                       ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;
  rhs.scale( -1.0 ); // TODO decide if we want this here
  solution.zero();

//  GEOSX_LOG_RANK_0( "Before PhaseFieldDamageFEM::SolveSystem" );
//  std::cout << matrix<<std::endl;
//  std::cout<< rhs << std::endl;

  SolverBase::solveSystem( dofManager, matrix, rhs, solution );

  if( getLogLevel() == 2 )
  {
    GEOSX_LOG_RANK_0( "After PhaseFieldDamageFEM::SolveSystem" );
    GEOSX_LOG_RANK_0( "\nSolution\n" );
    std::cout << solution;
  }
}

void PhaseFieldDamageFEM::applyDirichletBCImplicit( real64 const time,
                                                    DofManager const & dofManager,
                                                    DomainPartition & domain,
                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                    arrayView1d< real64 > const & localRhs )

{
  FieldSpecificationManager const & fsManager = getGlobalState().getFieldSpecificationManager();
  fsManager.apply( time,
                   domain,
                   "nodeManager",
                   m_fieldName,
                   [&]( FieldSpecificationBase const & bc, string const &,
                        SortedArrayView< localIndex const > const & targetSet,
                        Group & targetGroup,
                        string const GEOSX_UNUSED_PARAM( fieldName ) ) -> void
  {
    bc.applyBoundaryConditionToSystem< FieldSpecificationEqual,
                                       parallelDevicePolicy< 32 > >( targetSet,
                                                                     time,
                                                                     targetGroup,
                                                                     m_fieldName,
                                                                     dofManager.getKey( m_fieldName ),
                                                                     dofManager.rankOffset(),
                                                                     localMatrix,
                                                                     localRhs );
  } );

  fsManager.applyFieldValue< serialPolicy >( time, domain, "ElementRegions", viewKeyStruct::coeffNameString() );
}

REGISTER_CATALOG_ENTRY( SolverBase, PhaseFieldDamageFEM, string const &,
                        Group * const )
} // namespace geosx
