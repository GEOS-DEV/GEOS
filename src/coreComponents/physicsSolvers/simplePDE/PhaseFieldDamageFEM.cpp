/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
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
#include "PhaseFieldDamagePressureFEMKernels.hpp"
#include <math.h>
#include <vector>

#include "common/TimingMacros.hpp"
#include "dataRepository/Group.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/mpiCommunications/NeighborCommunicator.hpp"

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
#include "discretizationMethods/NumericalMethodsManager.hpp"

#include "mesh/DomainPartition.hpp"

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
  m_damageName( "damage" ),
  m_pressureEffectsFlag( 0 )
{

  registerWrapper< string >( PhaseFieldDamageFEMViewKeys.timeIntegrationOption.key() ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "option for default time integration method" );

  // registerWrapper< string >( PhaseFieldDamageFEMViewKeys.fieldVarName.key(), &m_fieldName ).
  //   setInputFlag( InputFlags::REQUIRED ).
  //   setDescription( "name of field variable" );

  registerWrapper( viewKeyStruct::localDissipationOptionString(), &m_localDissipationOption ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Type of local dissipation function. Can be Linear or Quadratic" );

  registerWrapper( viewKeyStruct::irreversibilityFlagString(), &m_irreversibilityFlag ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "The flag to indicate whether to apply the irreversibility constraint" );

  registerWrapper( viewKeyStruct::damageUpperBoundString(), &m_damageUpperBound ).
    setApplyDefaultValue( 1.5 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "The upper bound of the damage" );
}

PhaseFieldDamageFEM::~PhaseFieldDamageFEM()
{
  // TODO Auto-generated destructor stub
}

void PhaseFieldDamageFEM::registerDataOnMesh( Group & meshBodies )
{
  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & regionNames )
  {
    NodeManager & nodes = mesh.getNodeManager();

    nodes.registerWrapper< real64_array >( "damage" )
      .setApplyDefaultValue( 0.0 )
      .setPlotLevel( PlotLevel::LEVEL_0 )
      .setDescription( "Damage variable" );

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [ &]( localIndex const, CellElementSubRegion & subRegion )
    {
    #if 0
      subRegion.registerWrapper( viewKeyStruct::coeffNameString(), &m_coeff ).
        setApplyDefaultValue( 0.0 ).
        setPlotLevel( PlotLevel::LEVEL_0 ).
        setDescription( "field variable representing the diffusion coefficient" );
    #endif


      subRegion.registerWrapper< string >( viewKeyStruct::solidModelNamesString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRestartFlags( RestartFlags::NO_WRITE ).
        setSizedFromParent( 0 );

      string & solidMaterialName = subRegion.getReference< string >( viewKeyStruct::solidModelNamesString() );
      solidMaterialName = SolverBase::getConstitutiveName< SolidBase >( subRegion );
      GEOSX_ERROR_IF( solidMaterialName.empty(), GEOSX_FMT( "SolidBase model not found on subregion {}", subRegion.getName() ) );

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
    this->setupSystem( domain, m_dofManager, m_localMatrix, m_rhs, m_solution, false );

    dtReturn = this->nonlinearImplicitStep( time_n,
                                            dt,
                                            cycleNumber,
                                            domain );
  }
  return dtReturn;
}

real64 PhaseFieldDamageFEM::explicitStep( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                          real64 const & dt,
                                          const int GEOSX_UNUSED_PARAM( cycleNumber ),
                                          DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{
  return dt;
}

void PhaseFieldDamageFEM::implicitStepComplete( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                                real64 const & GEOSX_UNUSED_PARAM( dt ),
                                                DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{}

void PhaseFieldDamageFEM::setupDofs( DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
                                     DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;
  dofManager.addField( m_damageName,
                       FieldLocation::Node,
                       1,
                       getMeshTargets() );

  dofManager.addCoupling( m_damageName,
                          m_damageName,
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
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    arrayView1d< globalIndex const > const & dofIndex = nodeManager.getReference< array1d< globalIndex > >( dofManager.getKey( m_damageName ) );

    // Initialize all entries to zero
    localMatrix.zero();
    localRhs.zero();

    if( m_pressureEffectsFlag == 1 )
    {
      //use pressurized phase-field kernels
      PhaseFieldDamagePressureKernelFactory kernelFactory( dofIndex,
                                                           dofManager.rankOffset(),
                                                           localMatrix,
                                                           localRhs,
                                                           m_damageName,
                                                           m_localDissipationOption=="Linear" ? 1 : 2 );

      finiteElement::
        regionBasedKernelApplication< parallelDevicePolicy<>,
                                      constitutive::DamageBase,
                                      CellElementSubRegion >( mesh,
                                                              regionNames,
                                                              this->getDiscretizationName(),
                                                              viewKeyStruct::solidModelNamesString(),
                                                              kernelFactory );
    }
    else //use standard phase-field kernels
    {
      PhaseFieldDamageKernelFactory kernelFactory( dofIndex,
                                                   dofManager.rankOffset(),
                                                   localMatrix,
                                                   localRhs,
                                                   m_damageName,
                                                   m_localDissipationOption=="Linear" ? 1 : 2 );

      finiteElement::
        regionBasedKernelApplication< parallelDevicePolicy<>,
                                      constitutive::DamageBase,
                                      CellElementSubRegion >( mesh,
                                                              regionNames,
                                                              this->getDiscretizationName(),
                                                              viewKeyStruct::solidModelNamesString(),
                                                              kernelFactory );
    }


  } );
}

void PhaseFieldDamageFEM::applySystemSolution( DofManager const & dofManager,
                                               arrayView1d< real64 const > const & localSolution,
                                               real64 const scalingFactor,
                                               DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  dofManager.addVectorToField( localSolution, m_damageName, m_damageName, scalingFactor );

  // Syncronize ghost nodes
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addFields( FieldLocation::Node, { m_damageName } );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                         mesh,
                                                         domain.getNeighbors(),
                                                         false );
  } );

}

void PhaseFieldDamageFEM::updateState( DomainPartition & domain )
{
  GEOSX_UNUSED_VAR( domain );
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

  // Apply the crack irreversibility constraint
  if( m_irreversibilityFlag )
  {
    applyIrreversibilityConstraint( dofManager, domain, localMatrix, localRhs );
  }

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

  RAJA::ReduceSum< parallelDeviceReduce, real64 > localSum( 0.0 );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & )
  {
    const NodeManager & nodeManager = mesh.getNodeManager();
    const arrayView1d< const integer > & ghostRank = nodeManager.ghostRank();

    const arrayView1d< const globalIndex > &
    dofNumber = nodeManager.getReference< array1d< globalIndex > >( dofManager.getKey( m_damageName ) );
    const globalIndex rankOffset = dofManager.rankOffset();


    forAll< parallelDevicePolicy<> >( nodeManager.size(),
                                      [localRhs, localSum, dofNumber, rankOffset, ghostRank] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      if( ghostRank[k] < 0 )
      {
        localIndex const localRow = LvArray::integerConversion< localIndex >( dofNumber[k] - rankOffset );
        localSum += localRhs[localRow] * localRhs[localRow];
      }
    } );
  } );
  const real64 localResidualNorm[2] = { localSum.get(), 1};

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

void PhaseFieldDamageFEM::applyDirichletBCImplicit( real64 const time,
                                                    DofManager const & dofManager,
                                                    DomainPartition & domain,
                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                    arrayView1d< real64 > const & localRhs )

{
  FieldSpecificationManager const & fsManager = FieldSpecificationManager::getInstance();
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & )
  {
    fsManager.template apply< NodeManager >( time,
                                             mesh,
                                             m_damageName,
                                             [&]( FieldSpecificationBase const & bc,
                                                  string const &,
                                                  SortedArrayView< localIndex const > const & targetSet,
                                                  NodeManager & targetGroup,
                                                  string const GEOSX_UNUSED_PARAM( damageName ) ) -> void
    {
      bc.applyBoundaryConditionToSystem< FieldSpecificationEqual,
                                         parallelDevicePolicy< 32 > >( targetSet,
                                                                       time,
                                                                       targetGroup,
                                                                       m_damageName,
                                                                       dofManager.getKey( m_damageName ),
                                                                       dofManager.rankOffset(),
                                                                       localMatrix,
                                                                       localRhs );
    } );

    //fsManager.applyFieldValue< serialPolicy >( time, mesh, viewKeyStruct::coeffNameString() );
  } );
}

void PhaseFieldDamageFEM::applyIrreversibilityConstraint( DofManager const & dofManager,
                                                          DomainPartition & domain,
                                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                          arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    arrayView1d< globalIndex const > const dofIndex = nodeManager.getReference< array1d< globalIndex > >( dofManager.getKey( m_fieldName ) );

    arrayView1d< real64 const > const nodalDamage = nodeManager.getReference< array1d< real64 > >( m_fieldName );

    globalIndex const rankOffSet = dofManager.rankOffset();

    real64 const damangeUpperBound = m_damageUpperBound;

    forAll< parallelDevicePolicy<> >( nodeManager.size(), [=] GEOSX_HOST_DEVICE ( localIndex const nodeIndex )
    {
      localIndex const dof = dofIndex[nodeIndex];

      if( dof > -1 )
      {
        real64 const damageAtNode = nodalDamage[nodeIndex];

        if( damageAtNode >= damangeUpperBound )
        {

          // Specify the contribution to rhs
          real64 rhsContribution;

          FieldSpecificationEqual::SpecifyFieldValue( dof,
                                                      rankOffSet,
                                                      localMatrix,
                                                      rhsContribution,
                                                      damangeUpperBound,
                                                      damageAtNode );

          globalIndex const localRow = dof - rankOffSet;

          if( localRow >= 0 && localRow < localRhs.size() )
          {
            localRhs[ localRow ] = rhsContribution;
          }
        }
      }
    } );

  } );
}


REGISTER_CATALOG_ENTRY( SolverBase, PhaseFieldDamageFEM, string const &,
                        Group * const )
} // namespace geosx
