/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file PhaseFieldDamageFEM.cpp
 */

#include "PhaseFieldDamageFEM.hpp"
#include "PhaseFieldDamageFEMKernels.hpp"
#include "PhaseFieldPressurizedDamageFEMKernels.hpp"

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "constitutive/solid/Damage.hpp"
#include "constitutive/solid/SolidBase.hpp"
#include "fieldSpecification/FieldSpecificationImpl.hpp"
#include "finiteElement/FiniteElementDiscretization.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/LogLevelsInfo.hpp"

namespace geos
{

using namespace dataRepository;
using namespace constitutive;

PhaseFieldDamageFEM::PhaseFieldDamageFEM( const string & name,
                                          Group * const parent ):
  PhysicsSolverBase( name, parent ),
  m_fieldName( "primaryField" ),
  m_fracturePressureTermFlag( 0 )
{

  registerWrapper< TimeIntegrationOption >( PhaseFieldDamageFEMViewKeys.timeIntegrationOption.key(), &m_timeIntegrationOption ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "option for default time integration method" );

  registerWrapper< string >( PhaseFieldDamageFEMViewKeys.fieldVarName.key(), &m_fieldName ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRef ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "name of field variable" );

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

  registerWrapper( viewKeyStruct::fracturePressureTermFlagString(), &m_fracturePressureTermFlag ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "The flag to indicate whether to add the fracture pressure contribution" );

  addLogLevel< logInfo::ResidualNorm >();
}

PhaseFieldDamageFEM::~PhaseFieldDamageFEM() = default;

void PhaseFieldDamageFEM::registerDataOnMesh( Group & meshBodies )
{
  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    string_array const & regionNames )
  {
    NodeManager & nodes = mesh.getNodeManager();

    nodes.registerWrapper< real64_array >( m_fieldName )
      .setApplyDefaultValue( 0.0 )
      .setPlotLevel( PlotLevel::LEVEL_0 )
      .setDescription( "Primary field variable" );

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [ &]( localIndex const, CellElementSubRegion & subRegion )
    {
      subRegion.registerWrapper( viewKeyStruct::coeffNameString(), &m_coeff ).
        setApplyDefaultValue( 0.0 ).
        setPlotLevel( PlotLevel::LEVEL_0 ).
        setDescription( "field variable representing the diffusion coefficient" );

      // TODO this should be in setConstitutiveNames
      setConstitutiveName< SolidBase >( subRegion, viewKeyStruct::solidModelNamesString(), "solid" );
    } );
  } );
}

void PhaseFieldDamageFEM::initializePreSubGroups()
{
  PhysicsSolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const, CellElementSubRegion & subRegion )
    {
      string const & solidModelName = subRegion.getReference< string >( viewKeyStruct::solidModelNamesString() );
      SolidBase & solidModel = subRegion.getConstitutiveModel< SolidBase >( solidModelName );

      ConstitutivePassThru< DamageBase >::execute( solidModel, [&]( auto & damageModel )
      {
        GEOS_ERROR_IF( damageModel.getExtDrivingForceFlag() && m_localDissipationOption != LocalDissipation::Linear,
                       GEOS_FMT( "{}: the external driving force (extDrivingForceFlag) is only supported "
                                 "with the Linear local dissipation function.", getName() ) );
      } );
    } );
  } );
}

real64 PhaseFieldDamageFEM::solverStep( real64 const & time_n,
                                        real64 const & dt,
                                        const int cycleNumber,
                                        DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;
  real64 dtReturn = dt;
  if( m_timeIntegrationOption == TimeIntegrationOption::ExplicitTransient )
  {
    dtReturn = explicitStep( time_n, dt, cycleNumber, domain );
  }
  else if( m_timeIntegrationOption == TimeIntegrationOption::ImplicitTransient ||
           m_timeIntegrationOption == TimeIntegrationOption::SteadyState )
  {
    this->setupSystem( domain, m_dofManager, m_localMatrix, m_rhs, m_solution, false );

    dtReturn = this->nonlinearImplicitStep( time_n,
                                            dt,
                                            cycleNumber,
                                            domain );
  }
  return dtReturn;
}

real64 PhaseFieldDamageFEM::explicitStep( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                          real64 const & dt,
                                          const int GEOS_UNUSED_PARAM( cycleNumber ),
                                          DomainPartition & GEOS_UNUSED_PARAM( domain ) )
{
  return dt;
}

void PhaseFieldDamageFEM::implicitStepComplete( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                                real64 const & GEOS_UNUSED_PARAM( dt ),
                                                DomainPartition & GEOS_UNUSED_PARAM( domain ) )
{}

void PhaseFieldDamageFEM::setupDofs( DomainPartition const & GEOS_UNUSED_PARAM( domain ),
                                     DofManager & dofManager ) const
{
  GEOS_MARK_FUNCTION;
  dofManager.addField( m_fieldName,
                       FieldLocation::Node,
                       1,
                       getMeshTargets() );

  dofManager.addCoupling( m_fieldName,
                          m_fieldName,
                          DofManager::Connector::Elem );

}

void PhaseFieldDamageFEM::assembleSystem( real64 const GEOS_UNUSED_PARAM( time_n ),
                                          real64 const dt,
                                          DomainPartition & domain,
                                          DofManager const & dofManager,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  localMatrix.zero();
  localRhs.zero();

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    arrayView1d< globalIndex const > const & dofIndex =
      nodeManager.getReference< array1d< globalIndex > >( dofManager.getKey( m_fieldName ) );

    auto const launchKernels = [&]( auto kernelFactory )
    {
      finiteElement::
        regionBasedKernelApplication< parallelDevicePolicy<>,
                                      constitutive::DamageBase,
                                      CellElementSubRegion >( mesh,
                                                              regionNames,
                                                              this->getDiscretizationName(),
                                                              viewKeyStruct::solidModelNamesString(),
                                                              kernelFactory );
    };

    if( m_fracturePressureTermFlag )
    {
      launchKernels( PhaseFieldPressurizedDamageKernelFactory( dofIndex,
                                                               dofManager.rankOffset(),
                                                               localMatrix,
                                                               localRhs,
                                                               dt,
                                                               m_fieldName,
                                                               m_localDissipationOption ) );
    }
    else
    {
      launchKernels( PhaseFieldDamageKernelFactory( dofIndex,
                                                    dofManager.rankOffset(),
                                                    localMatrix,
                                                    localRhs,
                                                    dt,
                                                    m_fieldName,
                                                    m_localDissipationOption ) );
    }
  } );
}

void PhaseFieldDamageFEM::applySystemSolution( DofManager const & dofManager,
                                               arrayView1d< real64 const > const & localSolution,
                                               real64 const scalingFactor,
                                               real64 const dt,
                                               DomainPartition & domain )
{
  GEOS_UNUSED_VAR( dt );
  GEOS_MARK_FUNCTION;

  dofManager.addVectorToField( localSolution, m_fieldName, m_fieldName, scalingFactor );

  // Syncronize ghost nodes
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {
    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addFields( FieldLocation::Node, { m_fieldName } );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                         mesh,
                                                         domain.getNeighbors(),
                                                         false );
  } );

}

void PhaseFieldDamageFEM::updateState( DomainPartition & domain )
{
  GEOS_UNUSED_VAR( domain );
}

void PhaseFieldDamageFEM::applyBoundaryConditions(
  real64 const time_n,
  real64 const dt, DomainPartition & domain,
  DofManager const & dofManager,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;
  applyDirichletBCImplicit( time_n + dt, dofManager, domain, localMatrix, localRhs );

  // Apply the crack irreversibility constraint
  if( m_irreversibilityFlag )
  {
    applyIrreversibilityConstraint( dofManager, domain, localMatrix, localRhs );
  }

  if( getLogLevel() == 2 )
  {
    GEOS_LOG_RANK_0( "After PhaseFieldDamageFEM::applyBoundaryConditions" );
    GEOS_LOG_RANK_0( "\nJacobian:\n" );
    std::cout << localMatrix.toViewConst();
    GEOS_LOG_RANK_0( "\nResidual:\n" );
    std::cout << localRhs;
  }
}

real64
PhaseFieldDamageFEM::calculateResidualNorm( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                            real64 const & GEOS_UNUSED_PARAM( dt ),
                                            DomainPartition const & domain,
                                            DofManager const & dofManager,
                                            arrayView1d< real64 const > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  RAJA::ReduceSum< parallelDeviceReduce, real64 > localSum( 0.0 );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                string_array const & )
  {
    const NodeManager & nodeManager = mesh.getNodeManager();
    const arrayView1d< const integer > & ghostRank = nodeManager.ghostRank();

    const arrayView1d< const globalIndex > &
    dofNumber = nodeManager.getReference< array1d< globalIndex > >( dofManager.getKey( m_fieldName ) );
    const globalIndex rankOffset = dofManager.rankOffset();


    forAll< parallelDevicePolicy<> >( nodeManager.size(),
                                      [localRhs, localSum, dofNumber, rankOffset, ghostRank] GEOS_HOST_DEVICE ( localIndex const k )
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

  const int rank = MpiWrapper::commRank( MPI_COMM_GEOS );
  const int size = MpiWrapper::commSize( MPI_COMM_GEOS );
  array1d< real64 > globalValues( size * 2 );

  // Everything is done on rank 0
  MpiWrapper::gather( localResidualNorm,
                      2,
                      globalValues.data(),
                      2,
                      0,
                      MPI_COMM_GEOS );

  if( rank==0 )
  {
    for( int r=0; r<size; ++r )
    {
      // sum/max across all ranks
      globalResidualNorm[0] += globalValues[r*2];
      globalResidualNorm[1] = std::max( globalResidualNorm[1], globalValues[r*2+1] );
    }
  }

  MpiWrapper::bcast( globalResidualNorm, 2, 0, MPI_COMM_GEOS );


  const real64 residual = sqrt( globalResidualNorm[0] ) / ( globalResidualNorm[1] );

  GEOS_LOG_LEVEL_RANK_0_NLR( logInfo::ResidualNorm,
                             GEOS_FMT( "        ( R{} ) = ( {:4.2e} )", coupledSolverAttributePrefix(), residual ))


  getConvergenceStats().setResidualValue( GEOS_FMT( "R{}", coupledSolverAttributePrefix()), residual );

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
                                                               string_array const & )
  {
    fsManager.template apply< NodeManager >( time,
                                             mesh,
                                             m_fieldName,
                                             [&]( FieldSpecification const & bc,
                                                  string const &,
                                                  SortedArrayView< localIndex const > const & targetSet,
                                                  NodeManager & targetGroup,
                                                  string const GEOS_UNUSED_PARAM( fieldName ) ) -> void
    {
      FieldSpecificationImpl::applyBoundaryConditionToSystem< FieldSpecificationEqual,
                                                              parallelDevicePolicy< > >( bc,
                                                                                         targetSet,
                                                                                         time,
                                                                                         targetGroup,
                                                                                         m_fieldName,
                                                                                         dofManager.getKey( m_fieldName ),
                                                                                         dofManager.rankOffset(),
                                                                                         localMatrix,
                                                                                         localRhs );
    } );

    fsManager.applyFieldValue< serialPolicy >( time, mesh, viewKeyStruct::coeffNameString() );
  } );
}

void PhaseFieldDamageFEM::applyIrreversibilityConstraint( DofManager const & dofManager,
                                                          DomainPartition & domain,
                                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                          arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    arrayView1d< globalIndex const > const dofIndex = nodeManager.getReference< array1d< globalIndex > >( dofManager.getKey( m_fieldName ) );

    arrayView1d< real64 const > const nodalDamage = nodeManager.getReference< array1d< real64 > >( m_fieldName );

    globalIndex const rankOffSet = dofManager.rankOffset();

    real64 const damageUpperBound = m_damageUpperBound;

    forAll< parallelDevicePolicy<> >( nodeManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const nodeIndex )
    {
      localIndex const dof = dofIndex[nodeIndex];

      if( dof > -1 )
      {
        real64 const damageAtNode = nodalDamage[nodeIndex];

        if( damageAtNode >= damageUpperBound )
        {

          // Specify the contribution to rhs
          real64 rhsContribution;

          FieldSpecificationEqual::SpecifyFieldValue( dof,
                                                      rankOffSet,
                                                      localMatrix,
                                                      rhsContribution,
                                                      1.0,
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

void PhaseFieldDamageFEM::saveSequentialIterationState( DomainPartition & GEOS_UNUSED_PARAM( domain ) )
{
  // nothing to save yet
}

REGISTER_CATALOG_ENTRY( PhysicsSolverBase, PhaseFieldDamageFEM, string const &,
                        Group * const )
} // namespace geos
