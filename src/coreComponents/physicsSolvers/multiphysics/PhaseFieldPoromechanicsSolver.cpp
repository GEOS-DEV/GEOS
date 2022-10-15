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
 * @file PhaseFieldPoromechanicsSolver.cpp
 *
 */

#include "PhaseFieldPoromechanicsSolver.hpp"

#include "constitutive/ConstitutiveManager.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "fieldSpecification/TractionBoundaryCondition.hpp"
#include "finiteElement/Kinematics.h"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/simplePDE/PhaseFieldDamageFEM.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/multiphysics/SinglePhasePoromechanicsSolver.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

PhaseFieldPoromechanicsSolver::PhaseFieldPoromechanicsSolver( const string & name,
                                                              Group * const parent ):
  SolverBase( name, parent ),
  m_poromechanicsSolverName(),
  m_damageSolverName(),
  m_couplingTypeOption( CouplingTypeOption::FixedStress )

{
  registerWrapper( viewKeyStruct::poromechanicsSolverNameString(), &m_poromechanicsSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription(
    "Name of the poromechanics solver to use in the PhaseFieldPoromechanics solver" );

  registerWrapper( viewKeyStruct::damageSolverNameString(), &m_damageSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription(
    "Name of the damage mechanics solver to use in the PhaseFieldPoromechanics solver" );

  registerWrapper( viewKeyStruct::couplingTypeOptionString(), &m_couplingTypeOption ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Coupling option. Valid options:\n* " + EnumStrings< CouplingTypeOption >::concat( "\n* " ) );

  registerWrapper( viewKeyStruct::subcyclingOptionString(), &m_subcyclingOption ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "turn on subcycling on each load step" );

}

void PhaseFieldPoromechanicsSolver::registerDataOnMesh( Group & meshBodies )
{
  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & meshLevel,
                                                    arrayView1d< string const > const & )
  {
    ElementRegionManager & elemManager = meshLevel.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion,
                                      FaceElementSubRegion >( [&] ( auto & elementSubRegion )
    {
      elementSubRegion.template registerWrapper< array1d< real64 > >( viewKeyStruct::totalMeanStressString() ).
        setDescription( "Total Mean Stress" );
      elementSubRegion.template registerWrapper< array1d< real64 > >( viewKeyStruct::oldTotalMeanStressString() ).
        setDescription( "Total Mean Stress" );
    } );
  } );
}

void PhaseFieldPoromechanicsSolver::implicitStepSetup( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                                       real64 const & GEOSX_UNUSED_PARAM( dt ),
                                                       DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const totalMeanStress =
      elemManager.constructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::totalMeanStressString() );

    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > oldTotalMeanStress =
      elemManager.constructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::oldTotalMeanStressString() );

    //***** loop over all elements and initialize the derivative arrays *****
    forAllElemsInMesh( mesh, [ &]( localIndex const er, localIndex const esr, localIndex const k )
    {
      oldTotalMeanStress[er][esr][k] = totalMeanStress[er][esr][k];
    } );
  } );
}

void PhaseFieldPoromechanicsSolver::implicitStepComplete( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                                          real64 const & GEOSX_UNUSED_PARAM( dt ),
                                                          DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{}

void PhaseFieldPoromechanicsSolver::postProcessInput()
{
  if( m_couplingTypeOption == CouplingTypeOption::FixedStress )
  {
    // For this coupled solver the minimum number of Newton Iter should be 0 for both flow and solid solver otherwise it
    // will never converge.
    SinglePhasePoromechanicsSolver &
    poromechanicsSolver = this->getParent().getGroup< SinglePhasePoromechanicsSolver >( m_poromechanicsSolverName );
    integer & minNewtonIterPoro = poromechanicsSolver.getNonlinearSolverParameters().m_minIterNewton;

    PhaseFieldDamageFEM &
    damageSolver = this->getParent().getGroup< PhaseFieldDamageFEM >( m_damageSolverName );
    integer & minNewtonIterDamage = damageSolver.getNonlinearSolverParameters().m_minIterNewton;

    minNewtonIterPoro = 0;
    minNewtonIterDamage = 0;
  }
}

void PhaseFieldPoromechanicsSolver::initializePostInitialConditionsPreSubGroups()
{}

PhaseFieldPoromechanicsSolver::~PhaseFieldPoromechanicsSolver()
{
  // TODO Auto-generated destructor stub
}

void PhaseFieldPoromechanicsSolver::resetStateToBeginningOfStep( DomainPartition & domain )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const totalMeanStress =
      elemManager.constructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::totalMeanStressString() );

    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > oldTotalMeanStress =
      elemManager.constructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::oldTotalMeanStressString() );

    //***** loop over all elements and initialize the derivative arrays *****
    forAllElemsInMesh( mesh, [ &]( localIndex const er,
                                   localIndex const esr,
                                   localIndex const k )
    {
      totalMeanStress[er][esr][k] = oldTotalMeanStress[er][esr][k];
    } );
  } );
}

real64 PhaseFieldPoromechanicsSolver::solverStep( real64 const & time_n,
                                                  real64 const & dt,
                                                  int const cycleNumber,
                                                  DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  real64 dtReturn = dt;
  if( m_couplingTypeOption == CouplingTypeOption::FixedStress )
  {
    dtReturn = splitOperatorStep( time_n, dt, cycleNumber, domain );
  }
  else if( m_couplingTypeOption == CouplingTypeOption::TightlyCoupled )
  {
    GEOSX_ERROR( "CouplingTypeOption::FullyImplicit not yet implemented" );
  }
  return dtReturn;
}

real64 PhaseFieldPoromechanicsSolver::splitOperatorStep( real64 const & time_n,
                                                         real64 const & dt,
                                                         integer const cycleNumber,
                                                         DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  real64 dtReturn = dt;
  real64 dtReturnTemporary;

  SinglePhasePoromechanicsSolver &
  poromechanicsSolver = this->getParent().getGroup< SinglePhasePoromechanicsSolver >( m_poromechanicsSolverName );

  PhaseFieldDamageFEM &
  damageSolver = this->getParent().getGroup< PhaseFieldDamageFEM >( m_damageSolverName );

  damageSolver.setupSystem( domain,
                            damageSolver.getDofManager(),
                            damageSolver.getLocalMatrix(),
                            damageSolver.getSystemRhs(),
                            damageSolver.getSystemSolution(),
                            true );

  poromechanicsSolver.setupSystem( domain,
                                   poromechanicsSolver.getDofManager(),
                                   poromechanicsSolver.getLocalMatrix(),
                                   poromechanicsSolver.getSystemRhs(),
                                   poromechanicsSolver.getSystemSolution() );

  damageSolver.implicitStepSetup( time_n, dt, domain );

  poromechanicsSolver.implicitStepSetup( time_n, dt, domain );

  this->implicitStepSetup( time_n, dt, domain );

  NonlinearSolverParameters & solverParams = getNonlinearSolverParameters();
  integer & iter = solverParams.m_numNewtonIterations;
  iter = 0;
  bool isConverged = false;
  while( iter < solverParams.m_maxIterNewton )
  {
    if( iter == 0 )
    {
      // reset the states of all slave solvers if any of them has been reset
      damageSolver.resetStateToBeginningOfStep( domain );
      poromechanicsSolver.resetStateToBeginningOfStep( domain );
      resetStateToBeginningOfStep( domain );
    }

    GEOSX_LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1 << ", PoromechanicsSolver: " );

    applyDamageOnTractionBC( domain ); 

    dtReturnTemporary = poromechanicsSolver.nonlinearImplicitStep( time_n,
                                                                   dtReturn,
                                                                   cycleNumber,
                                                                   domain );

    if( dtReturnTemporary < dtReturn )
    {
      iter = 0;
      dtReturn = dtReturnTemporary;
      continue;
    }

    if( poromechanicsSolver.getNonlinearSolverParameters().m_numNewtonIterations == 0 && iter > 0 )
    {
      GEOSX_LOG_LEVEL_RANK_0( 1, "***** The iterative coupling has converged in " << iter << " iterations! *****\n" );
      isConverged = true;
      break;
    }
    else if( m_subcyclingOption == 0 && iter > 0 )
    {
      GEOSX_LOG_LEVEL_RANK_0( 1, "***** Single Pass solver, no subcycling *****\n" );
      isConverged = true;
      break;
    }

    GEOSX_LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1 << ", DamageSolver: " );

    dtReturnTemporary = damageSolver.nonlinearImplicitStep( time_n,
                                                            dtReturn,
                                                            cycleNumber,
                                                            domain );

    mapDamageAndGradientToQuadrature( domain );

    //std::cout << "Here: " << dtReturnTemporary << std::endl;

    if( dtReturnTemporary < dtReturn )
    {
      iter = 0;
      dtReturn = dtReturnTemporary;
      continue;
    }

    if( m_subcyclingOption == 0 )
    {
      GEOSX_LOG_LEVEL_RANK_0( 1, "***** Single Pass solver, no subcycling *****\n" );
      isConverged = true;
      break;
    }

    ++iter;
  }

  GEOSX_ERROR_IF( !isConverged, "PhaseFieldPoromechanicsSolver::SplitOperatorStep() did not converge" );

  damageSolver.implicitStepComplete( time_n, dt, domain );
  poromechanicsSolver.implicitStepComplete( time_n, dt, domain );
  this->implicitStepComplete( time_n, dt, domain );

  return dtReturn;
}

void PhaseFieldPoromechanicsSolver::mapDamageAndGradientToQuadrature( DomainPartition & domain )
{

  GEOSX_MARK_FUNCTION;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const xNodes = nodeManager.referencePosition();

    PhaseFieldDamageFEM const &
    damageSolver = this->getParent().getGroup< PhaseFieldDamageFEM >( m_damageSolverName );

    string const & damageFieldName = damageSolver.getFieldName();

    //should get reference to damage field here.
    arrayView1d< real64 const > const nodalDamage = nodeManager.getReference< array1d< real64 > >( damageFieldName );

    ElementRegionManager & elemManager = mesh.getElemManager();

    // begin region loop
    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [this, xNodes, nodalDamage]
                                                                ( localIndex const,
                                                                CellElementSubRegion & elementSubRegion )
    {
      string const & solidModelName = elementSubRegion.getReference< string >( SolidMechanicsLagrangianFEM::viewKeyStruct::solidMaterialNamesString());
      constitutive::SolidBase &
      solidModel = elementSubRegion.getConstitutiveModel< constitutive::SolidBase >( solidModelName );

      ConstitutivePassThru< DamageBase >::execute( solidModel, [this, &elementSubRegion, xNodes, nodalDamage]( auto & damageModel )
      {
        using CONSTITUTIVE_TYPE = TYPEOFREF( damageModel );
        typename CONSTITUTIVE_TYPE::KernelWrapper constitutiveUpdate = damageModel.createKernelUpdates();

        arrayView2d< real64 > const damageFieldOnMaterial = constitutiveUpdate.m_newDamage;
        arrayView3d< real64 > const damageGradOnMaterial = constitutiveUpdate.m_damageGrad;
        arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemNodes = elementSubRegion.nodeList();

        auto const & elemsToNodes = elementSubRegion.nodeList().toViewConst();

        finiteElement::FiniteElementBase const &
        fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( m_discretizationName );

        finiteElement::dispatch3D( fe, [xNodes, nodalDamage, &elementSubRegion, damageFieldOnMaterial, damageGradOnMaterial, elemNodes, elemsToNodes]( auto & finiteElement )
        {
          forAll< serialPolicy >( elementSubRegion.size(), [xNodes, nodalDamage, damageFieldOnMaterial, damageGradOnMaterial, elemNodes, elemsToNodes] ( localIndex const k )
          {
            using FE_TYPE = TYPEOFREF( finiteElement );
            constexpr localIndex numNodesPerElement = FE_TYPE::numNodes;
            constexpr localIndex n_q_points = FE_TYPE::numQuadraturePoints;

            real64 xLocal[ numNodesPerElement ][ 3 ];
            real64 nodalDamageLocal[ numNodesPerElement ];

            for( localIndex a = 0; a < numNodesPerElement; ++a )
            {
              localIndex const localNodeIndex = elemsToNodes( k, a );

              for( int dim=0; dim < 3; ++dim )
              {
                xLocal[a][dim] = xNodes[ localNodeIndex ][dim];
              }

              nodalDamageLocal[ a ] = nodalDamage[ localNodeIndex ];
            }

            for( localIndex q = 0; q < n_q_points; ++q )
            {
              real64 N[ numNodesPerElement ];
              FE_TYPE::calcN( q, N );

              real64 dNdX[ numNodesPerElement ][ 3 ];

              real64 const detJ = FE_TYPE::calcGradN( q, xLocal, dNdX );

              GEOSX_UNUSED_VAR( detJ );

              real64 qDamage = 0.0;
              real64 qDamageGrad[3] = {0, 0, 0};
              FE_TYPE::valueAndGradient( N, dNdX, nodalDamageLocal, qDamage, qDamageGrad );

              damageFieldOnMaterial( k, q ) = qDamage;

              for( int dim=0; dim < 3; ++dim )
              {
                damageGradOnMaterial[k][q][dim] = qDamageGrad[dim];
              }
            }
          } );
        } );
      } );
    } );
  } );
}


void PhaseFieldPoromechanicsSolver::applyDamageOnTractionBC( DomainPartition & domain )
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  GEOSX_MARK_FUNCTION;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    PhaseFieldDamageFEM const &
    damageSolver = this->getParent().getGroup< PhaseFieldDamageFEM >( m_damageSolverName );

    string const & damageFieldName = damageSolver.getFieldName();

    // Get an array of nodal damage values
    arrayView1d< real64 const > const nodalDamage = nodeManager.getReference< array1d< real64 > >( damageFieldName );

    fsManager.forSubGroups< TractionBoundaryCondition >( [&] ( TractionBoundaryCondition & fs )
    {
      string_array const targetPath = stringutilities::tokenize( fs.getObjectPath(), "/" );

      dataRepository::Group * targetGroup = &mesh;

      dataRepository::Group * const elemRegionSubGroup = targetGroup->getGroupPointer( ElementRegionManager::groupKeyStruct::elementRegionsGroup() );

      if( elemRegionSubGroup != nullptr )
      {
        targetGroup = elemRegionSubGroup;
      }

      dataRepository::Group * const elemSubRegionSubGroup = targetGroup->getGroupPointer( ElementRegionBase::viewKeyStruct::elementSubRegions() );
      if( elemSubRegionSubGroup != nullptr )
      {
        targetGroup = elemSubRegionSubGroup;
      }

      targetGroup = &targetGroup->getGroup( targetPath[0] );

      Group & target = *targetGroup;

      dataRepository::Group const & setGroup = target.getGroup( ObjectManagerBase::groupKeyStruct::setsString() );
      string_array setNames = fs.getSetNames();
      for( auto & setName : setNames )
      {
        if( setGroup.hasWrapper( setName ) )
        {
          SortedArrayView< localIndex const > const & targetSet = setGroup.getReference< SortedArray< localIndex > >( setName );

          fs.reinitScaleSet( faceManager,
                             targetSet,
                             nodalDamage );
        }
      }
    } );
  } );

}

REGISTER_CATALOG_ENTRY( SolverBase, PhaseFieldPoromechanicsSolver, string const &, Group * const )

} /* namespace geosx */
