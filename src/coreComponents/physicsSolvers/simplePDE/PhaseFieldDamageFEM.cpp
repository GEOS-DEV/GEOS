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

// this should be part of the input file

double myFunc2( double, double, double )
{
  return 0;
  // return pow(x, 2) + pow(y, 2) + pow(z, 2) + 6;
//  return x * (1 - x) * y * (1 - y) * z * (1 - z) -
//         2 * (x - 1) * x * (y - 1) * y - 2 * (x - 1) * x * (z - 1) * z -
//         2 * (y - 1) * y * (z - 1) * z;
}

///////////////////////////////////

namespace geosx
{

namespace dataRepository
{
namespace keys
{}
} // namespace dataRepository

using namespace dataRepository;
using namespace constitutive;

PhaseFieldDamageFEM::PhaseFieldDamageFEM( const std::string & name,
                                          Group * const parent ):
  SolverBase( name, parent ),
  m_fieldName( "primaryField" ),
  m_solidModelNames()
{

  registerWrapper< string >( PhaseFieldDamageFEMViewKeys.timeIntegrationOption.Key() )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "option for default time integration method" );

  registerWrapper< string >( PhaseFieldDamageFEMViewKeys.fieldVarName.Key(), &m_fieldName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "name of field variable" );

  registerWrapper( viewKeyStruct::localDissipationOption, &m_localDissipationOption )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Type of local dissipation function. Can be Linear or Quadratic" );

  registerWrapper( viewKeyStruct::lengthScale, &m_lengthScale )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "lenght scale l in the phase-field equation" );

  registerWrapper( viewKeyStruct::criticalFractureEnergy, &m_criticalFractureEnergy )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "critical fracture energy" );

  registerWrapper( viewKeyStruct::solidModelNamesString, &m_solidModelNames )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "name of solid constitutive model" );
}

PhaseFieldDamageFEM::~PhaseFieldDamageFEM()
{
  // TODO Auto-generated destructor stub
}

void PhaseFieldDamageFEM::RegisterDataOnMesh( Group * const MeshBodies )
{
  for( auto & mesh : MeshBodies->GetSubGroups() )
  {

    MeshLevel *meshLevel = Group::group_cast< MeshBody * >( mesh.second )->getMeshLevel( 0 );

    NodeManager * const nodes = meshLevel->getNodeManager();

    nodes->registerWrapper< real64_array >( m_fieldName )
      ->setApplyDefaultValue( 0.0 )
      ->setPlotLevel( PlotLevel::LEVEL_0 )
      ->setDescription( "Primary field variable" );

    ElementRegionManager * const elemManager = meshLevel->getElemManager();

    elemManager->forElementSubRegions< CellElementSubRegion >( [ &]( CellElementSubRegion & subRegion )
    {
      subRegion.registerWrapper( viewKeyStruct::coeffName, &m_coeff )->
        setApplyDefaultValue( 0.0 )->
        setPlotLevel( PlotLevel::LEVEL_0 )->
        setDescription( "field variable representing the diffusion coefficient" );
    } );

  }
}

void PhaseFieldDamageFEM::PostProcessInput()
{
  SolverBase::PostProcessInput();

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

real64 PhaseFieldDamageFEM::SolverStep( real64 const & time_n,
                                        real64 const & dt,
                                        const int cycleNumber,
                                        DomainPartition & domain )
{
  real64 dtReturn = dt;
  if( m_timeIntegrationOption == timeIntegrationOption::ExplicitTransient )
  {
    dtReturn = ExplicitStep( time_n, dt, cycleNumber, domain );
  }
  else if( m_timeIntegrationOption ==
           timeIntegrationOption::ImplicitTransient ||
           m_timeIntegrationOption == timeIntegrationOption::SteadyState )
  {
    this->SetupSystem( domain, m_dofManager, m_localMatrix, m_localRhs, m_localSolution, false );

    dtReturn = this->NonlinearImplicitStep( time_n,
                                            dt,
                                            cycleNumber,
                                            domain );
  }
  return dtReturn;
}

real64 PhaseFieldDamageFEM::ExplicitStep(
  real64 const & GEOSX_UNUSED_PARAM( time_n ),
  real64 const & dt,
  const int GEOSX_UNUSED_PARAM( cycleNumber ),
  DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{
  return dt;
}

void PhaseFieldDamageFEM::SetupSystem( DomainPartition & domain,
                                       DofManager & dofManager,
                                       CRSMatrix< real64, globalIndex > & localMatrix,
                                       array1d< real64 > & localRhs,
                                       array1d< real64 > & localSolution,
                                       bool const setSparisty )
{
  GEOSX_MARK_FUNCTION;
  SolverBase::SetupSystem( domain, dofManager, localMatrix, localRhs, localSolution, setSparisty );
}

void PhaseFieldDamageFEM::ImplicitStepComplete(
  real64 const & GEOSX_UNUSED_PARAM( time_n ),
  real64 const & GEOSX_UNUSED_PARAM( dt ),
  DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{}

void PhaseFieldDamageFEM::SetupDofs(
  DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
  DofManager & dofManager ) const
{
  dofManager.addField( m_fieldName, DofManager::Location::Node );

  dofManager.addCoupling( m_fieldName,
                          m_fieldName,
                          DofManager::Connector::Elem );

}

void PhaseFieldDamageFEM::AssembleSystem( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                          real64 const GEOSX_UNUSED_PARAM( dt ),
                                          DomainPartition & domain,
                                          DofManager const & dofManager,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs )
{
  MeshLevel * const mesh = domain.getMeshBody( 0 )->getMeshLevel( 0 );
  NodeManager * const nodeManager = mesh->getNodeManager();

  arrayView1d< globalIndex const > const & dofIndex = nodeManager->getReference< array1d< globalIndex > >( dofManager.getKey( m_fieldName ) );

  // Initialize all entries to zero
  localMatrix.setValues< parallelDevicePolicy< 32 > >( 0 );
  localRhs.setValues< parallelDevicePolicy< 32 > >( 0 );

  finiteElement::
    regionBasedKernelApplication< serialPolicy,
                                  constitutive::SolidBase,
                                  CellElementSubRegion,
                                  PhaseFieldDamageKernel >( *mesh,
                                                            targetRegionNames(),
                                                            this->getDiscretizationName(),
                                                            m_solidModelNames,
                                                            dofIndex,
                                                            dofManager.rankOffset(),
                                                            localMatrix,
                                                            localRhs,
                                                            m_fieldName,
                                                            m_criticalFractureEnergy,
                                                            m_lengthScale,
                                                            m_localDissipationOption=="Linear" ? 1 : 2 );


}

void PhaseFieldDamageFEM::ApplySystemSolution( DofManager const & dofManager,
                                               arrayView1d< real64 const > const & localSolution,
                                               real64 const scalingFactor,
                                               DomainPartition & domain )
{
  MeshLevel * const mesh = domain.getMeshBody( 0 )->getMeshLevel( 0 );

  dofManager.addVectorToField( localSolution,
                               m_fieldName,
                               m_fieldName,
                               scalingFactor );

  // Syncronize ghost nodes
  std::map< string, string_array > fieldNames;
  fieldNames["node"].emplace_back( m_fieldName );

  CommunicationTools::SynchronizeFields( fieldNames,
                                         mesh,
                                         domain.getNeighbors() );
}

void PhaseFieldDamageFEM::ApplyBoundaryConditions(
  real64 const time_n,
  real64 const dt, DomainPartition & domain,
  DofManager const & dofManager,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs )
{
  ApplyDirichletBC_implicit( time_n + dt, dofManager, domain, localMatrix, localRhs );

  if( getLogLevel() == 2 )
  {
    GEOSX_LOG_RANK_0( "After PhaseFieldDamageFEM::ApplyBoundaryConditions" );
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
//    GEOSX_LOG_RANK_0( "After PhaseFieldDamageFEM::ApplyBoundaryConditions" );
//    GEOSX_LOG_RANK_0( "Jacobian: written to " << filename_mat );
//    GEOSX_LOG_RANK_0( "Residual: written to " << filename_rhs );
//  }
}

real64
PhaseFieldDamageFEM::CalculateResidualNorm( DomainPartition const & domain,
                                            DofManager const & dofManager,
                                            arrayView1d< real64 const > const & localRhs )
{
  const MeshLevel & mesh = *( domain.getMeshBody( 0 )->getMeshLevel( 0 ) );
  const NodeManager & nodeManager = *mesh.getNodeManager();
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
      for( localIndex dim = 0; dim < 3; ++dim )
      {
        localSum += localRhs[localRow + dim] * localRhs[localRow + dim];
      }
    }
  } );

  const real64 localResidualNorm[2] =
  { localSum.get(), 1 };

  // globalResidualNorm[0]: the sum of all the local sum(rhs^2).
  // globalResidualNorm[1]: max of max force of each rank. Basically max force globally
  real64 globalResidualNorm[2] = {0, 0};

  const int rank = MpiWrapper::Comm_rank( MPI_COMM_GEOSX );
  const int size = MpiWrapper::Comm_size( MPI_COMM_GEOSX );
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

void PhaseFieldDamageFEM::SolveSystem( DofManager const & dofManager,
                                       ParallelMatrix & matrix,
                                       ParallelVector & rhs,
                                       ParallelVector & solution )
{
  rhs.scale( -1.0 ); // TODO decide if we want this here
  solution.zero();

//  GEOSX_LOG_RANK_0( "Before PhaseFieldDamageFEM::SolveSystem" );
//  std::cout << matrix<<std::endl;
//  std::cout<< rhs << std::endl;

  SolverBase::SolveSystem( dofManager, matrix, rhs, solution );

  if( getLogLevel() == 2 )
  {
    GEOSX_LOG_RANK_0( "After PhaseFieldDamageFEM::SolveSystem" );
    GEOSX_LOG_RANK_0( "\nSolution\n" );
    std::cout << solution;
  }
}

void PhaseFieldDamageFEM::ApplyDirichletBC_implicit( real64 const time,
                                                     DofManager const & dofManager,
                                                     DomainPartition & domain,
                                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                     arrayView1d< real64 > const & localRhs )

{
  FieldSpecificationManager const & fsManager = FieldSpecificationManager::get();
  fsManager.Apply( time,
                   &domain,
                   "nodeManager",
                   m_fieldName,
                   [&]( FieldSpecificationBase const * const bc, string const &,
                        SortedArrayView< localIndex const > const & targetSet,
                        Group * const targetGroup,
                        string const GEOSX_UNUSED_PARAM( fieldName ) ) -> void
  {
    bc->ApplyBoundaryConditionToSystem< FieldSpecificationEqual,
                                        parallelDevicePolicy< 32 > >( targetSet,
                                                                      time,
                                                                      targetGroup,
                                                                      m_fieldName,
                                                                      dofManager.getKey( m_fieldName ),
                                                                      dofManager.rankOffset(),
                                                                      localMatrix,
                                                                      localRhs );
  } );

  fsManager.ApplyFieldValue< serialPolicy >( time, &domain, "ElementRegions", viewKeyStruct::coeffName );
}

REGISTER_CATALOG_ENTRY( SolverBase, PhaseFieldDamageFEM, std::string const &,
                        Group * const )
} // namespace geosx
