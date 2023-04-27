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
 * @file PhaseFieldFractureSolver.cpp
 *
 */

#include "PhaseFieldFractureSolver.hpp"

#include "constitutive/ConstitutiveManager.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "finiteElement/Kinematics.h"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

PhaseFieldFractureSolver::PhaseFieldFractureSolver( const string & name,
                                                    Group * const parent ):
  Base( name, parent )
{
  registerWrapper( viewKeyStruct::pressureEffectsString(), &m_pressureEffects ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( 0 ).
    setDescription( "consider background pressure effects (matrix and fracture pressures)" );
}

void PhaseFieldFractureSolver::registerDataOnMesh( Group & meshBodies )
{
  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & meshLevel,
                                                    arrayView1d< string const > const & )
  {
    ElementRegionManager & elemManager = meshLevel.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion,
                                      FaceElementSubRegion >( [&] ( auto & elementSubRegion )
    {
      if( m_pressureEffects )
      {
        //register fake pressures for testing purposes
        elementSubRegion.template registerWrapper< array1d< real64 > >( "hardCodedPMatrixName" ).
          setPlotLevel( PlotLevel::LEVEL_1 ).
          setDescription( "matrix pressure field for testing purposes only" );
        elementSubRegion.template registerWrapper< array1d< real64 > >( "hardCodedPFractureName" ).
          setPlotLevel( PlotLevel::LEVEL_1 ).
          setDescription( "fracture pressure field for testing purposes only" );
      }
    } );

    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const matrixPressure =
      elemManager.constructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( "hardCodedPMatrixName" );

    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > fracturePressure =
      elemManager.constructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( "hardCodedPFractureName" );

    // if( m_pressureEffects == 1 )
    // {
    //   //***** loop over all elements and write manufactured pressure fields *****
    //   forAllElemsInMesh( meshLevel, [ &]( localIndex const er,
    //                                       localIndex const esr,
    //                                       localIndex const k )
    //   {
    //     //make a non-trivial field to test
    //     matrixPressure[er][esr][k] = k*1.0e-2;
    //     fracturePressure[er][esr][k] = k*1.0e0;
    //   } );
    // }

  } );
}

PhaseFieldFractureSolver::~PhaseFieldFractureSolver()
{
  // TODO Auto-generated destructor stub
}

void PhaseFieldFractureSolver::postProcessInput()
{
  Base::postProcessInput();

  if( m_pressureEffects == 1 ) //this will add background pressure effects to solid and damage solvers - must be done before assemble
                               // routine
  {
    this->solidMechanicsSolver()->setPressureEffects();
    this->damageSolver()->setPressureEffects();
    //imposeFakeBackgroundPressures( domain ); //DOMAIN IS COMING OUT EMPTY
  }

  GEOSX_WARNING_IF( getNonlinearSolverParameters().m_couplingType == NonlinearSolverParameters::CouplingType::FullyImplicit,
                    "FullyImplicit coupling not implemented for this solver. A sequential coupling approach will be used." );
  getNonlinearSolverParameters().m_couplingType = NonlinearSolverParameters::CouplingType::Sequential;
}

// void PhaseFieldFractureSolver::imposeFakeBackgroundPressures( DomainPartition & domain )
// {
//   forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
//                                                                 MeshLevel & mesh,
//                                                                 arrayView1d< string const > const & )
//   {
//     ElementRegionManager & elemManager = mesh.getElemManager();
//     std::cout<<"entering meshBody loop\n";
//     ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const matrixPressure =
//       elemManager.constructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( "hardCodedPMatrixName" );

//     ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > fracturePressure =
//       elemManager.constructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( "hardCodedPFractureName" );

//     //***** loop over all elements and write manufactured pressure fields *****
//     forAllElemsInMesh( mesh, [ &]( localIndex const er,
//                                    localIndex const esr,
//                                    localIndex const k )
//     {
//       //make a non-trivial field to test
//       matrixPressure[er][esr][k] = k*1.0e-2;
//       fracturePressure[er][esr][k] = k*1.0e0;
//     } );
//   } );
// }

// real64 PhaseFieldFractureSolver::splitOperatorStep( real64 const & time_n,
//                                                     real64 const & dt,
//                                                     integer const cycleNumber,
//                                                     DomainPartition & domain )
// {
//   GEOSX_MARK_FUNCTION;
//   real64 dtReturn = dt;
//   real64 dtReturnTemporary;

//   SolidMechanicsLagrangianFEM &
//   solidSolver = this->getParent().getGroup< SolidMechanicsLagrangianFEM >( m_solidSolverName );

//   PhaseFieldDamageFEM &
//   damageSolver = this->getParent().getGroup< PhaseFieldDamageFEM >( m_damageSolverName );

//   if( m_pressureEffects == 1 ) //this will add background pressure effects to solid and damage solvers - must be done before assemble
//                                // routine
//   {
//     solidSolver.setPressureEffects();
//     damageSolver.setPressureEffects();
//     imposeFakeBackgroundPressures( domain );
//   }

//   damageSolver.setupSystem( domain,
//                             damageSolver.getDofManager(),
//                             damageSolver.getLocalMatrix(),
//                             damageSolver.getSystemRhs(),
//                             damageSolver.getSystemSolution(),
//                             true );

//   solidSolver.setupSystem( domain,
//                            solidSolver.getDofManager(),
//                            solidSolver.getLocalMatrix(),
//                            solidSolver.getSystemRhs(),
//                            solidSolver.getSystemSolution() );

//   damageSolver.implicitStepSetup( time_n, dt, domain );

//   solidSolver.implicitStepSetup( time_n, dt, domain );

//   this->implicitStepSetup( time_n, dt, domain );

//   NonlinearSolverParameters & solverParams = getNonlinearSolverParameters();
//   integer & iter = solverParams.m_numNewtonIterations;
//   iter = 0;
//   bool isConverged = false;
//   while( iter < solverParams.m_maxIterNewton )
//   {
//     if( iter == 0 )
//     {
//       // reset the states of all slave solvers if any of them has been reset
//       damageSolver.resetStateToBeginningOfStep( domain );
//       solidSolver.resetStateToBeginningOfStep( domain );
//       resetStateToBeginningOfStep( domain );
//     }

//     GEOSX_LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1 << ", MechanicsSolver: " );

//     dtReturnTemporary = solidSolver.nonlinearImplicitStep( time_n,
//                                                            dtReturn,
//                                                            cycleNumber,
//                                                            domain );

//     if( dtReturnTemporary < dtReturn )
//     {
//       iter = 0;
//       dtReturn = dtReturnTemporary;
//       continue;
//     }

//     if( solidSolver.getNonlinearSolverParameters().m_numNewtonIterations == 0 && iter > 0 )
//     {
//       GEOSX_LOG_LEVEL_RANK_0( 1, "***** The iterative coupling has converged in " << iter << " iterations! *****\n" );
//       isConverged = true;
//       break;
//     }
//     else if( m_subcyclingOption == 0 && iter > 0 )
//     {
//       GEOSX_LOG_LEVEL_RANK_0( 1, "***** Single Pass solver, no subcycling *****\n" );
//       isConverged = true;
//       break;
//     }

//     GEOSX_LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1 << ", DamageSolver: " );

//     dtReturnTemporary = damageSolver.nonlinearImplicitStep( time_n,
//                                                             dtReturn,
//                                                             cycleNumber,
//                                                             domain );

//     if( m_pressureEffects )
//     {
//       mapDamageAndGradientToQuadrature( domain );
//     }
//     else
//     {
//       mapDamageToQuadrature( domain );
//     }
//     //std::cout << "Here: " << dtReturnTemporary << std::endl;

//     if( dtReturnTemporary < dtReturn )
//     {
//       iter = 0;
//       dtReturn = dtReturnTemporary;
//       continue;
//     }

//     if( m_subcyclingOption == 0 )
//     {
//       GEOSX_LOG_LEVEL_RANK_0( 1, "***** Single Pass solver, no subcycling *****\n" );
//       isConverged = true;
//       break;
//     }

//     ++iter;
//   }

//   GEOSX_ERROR_IF( !isConverged, "PhaseFieldFractureSolver::SplitOperatorStep() did not converge" );

//   damageSolver.implicitStepComplete( time_n, dt, domain );
//   solidSolver.implicitStepComplete( time_n, dt, domain );
//   this->implicitStepComplete( time_n, dt, domain );

//   return dtReturn;
// }

// void PhaseFieldFractureSolver::mapSolutionBetweenSolvers( DomainPartition & domain, integer const solverType )
// {

//   GEOSX_MARK_FUNCTION;
//   if( solverType ==  static_cast< integer >( SolverType::Damage ) )
//   {
//     forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
//                                                                   MeshLevel & mesh,
//                                                                   arrayView1d< string const > const & regionNames )
//     {
//       NodeManager & nodeManager = mesh.getNodeManager();

//       string const & damageFieldName = damageSolver()->getFieldName();

//       string const & discretizationName = damageSolver()->getDiscretizationName();

//       //should get reference to damage field here.
//       arrayView1d< real64 const > const nodalDamage = nodeManager.getReference< array1d< real64 > >( damageFieldName );

//       ElementRegionManager & elemManager = mesh.getElemManager();

//       // begin region loop
//       elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [discretizationName, nodalDamage]
//                                                                   ( localIndex const,
//                                                                   CellElementSubRegion & elementSubRegion )
//       {
//         string const & solidModelName = elementSubRegion.getReference< string >(
// SolidMechanicsLagrangianFEM::viewKeyStruct::solidMaterialNamesString());
//         constitutive::SolidBase &
//         solidModel = elementSubRegion.getConstitutiveModel< constitutive::SolidBase >( solidModelName );

//         ConstitutivePassThru< DamageBase >::execute( solidModel, [&elementSubRegion, discretizationName, nodalDamage]( auto & damageModel
// )
//         {
//           using CONSTITUTIVE_TYPE = TYPEOFREF( damageModel );
//           typename CONSTITUTIVE_TYPE::KernelWrapper constitutiveUpdate = damageModel.createKernelUpdates();

//           arrayView2d< real64 > const damageFieldOnMaterial = constitutiveUpdate.m_newDamage;
//           arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemToNodes = elementSubRegion.nodeList();

//           finiteElement::FiniteElementBase const &
//           fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( discretizationName );

//           finiteElement::FiniteElementDispatchHandler< ALL_FE_TYPES >::dispatch3D( fe, [=, &elementSubRegion] ( auto & finiteElement )
//           {
//             using FE_TYPE = TYPEOFREF( finiteElement );

//             DamageInterpolationKernel< FE_TYPE > interpolationKernel( elementSubRegion );

//             interpolationKernel.interpolateDamage( elemToNodes, nodalDamage, damageFieldOnMaterial );
//           } );
//         } );
//       } );
//     } );
//   }
// }

void PhaseFieldFractureSolver::mapSolutionBetweenSolvers( DomainPartition & domain, integer const solverType )
//void PhaseFieldFractureSolver::mapDamageAndGradientToQuadrature( DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  if( solverType ==  static_cast< integer >( SolverType::Damage ) )
  {
    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                  MeshLevel & mesh,
                                                                  arrayView1d< string const > const & regionNames )
    {
      NodeManager & nodeManager = mesh.getNodeManager();

      arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const xNodes = nodeManager.referencePosition();

      string const & damageFieldName = damageSolver()->getFieldName();
      string const & discretizationName = damageSolver()->getDiscretizationName();

      //should get reference to damage field here.
      arrayView1d< real64 const > const nodalDamage = nodeManager.getReference< array1d< real64 > >( damageFieldName );

      ElementRegionManager & elemManager = mesh.getElemManager();

      // ConstitutiveManager & constitutiveManager = domain.getGroup< ConstitutiveManager >( keys::ConstitutiveManager );

      // ElementRegionManager::ConstitutiveRelationAccessor< ConstitutiveBase >
      // constitutiveRelations = elemManager.constructFullConstitutiveAccessor< ConstitutiveBase >( constitutiveManager );
      // begin region loop
      elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [discretizationName, xNodes, nodalDamage]
                                                                  ( localIndex const,
                                                                  CellElementSubRegion & elementSubRegion )
      {
        string const & solidModelName = elementSubRegion.getReference< string >( SolidMechanicsLagrangianFEM::viewKeyStruct::solidMaterialNamesString());
        constitutive::SolidBase &
        solidModel = elementSubRegion.getConstitutiveModel< constitutive::SolidBase >( solidModelName );
//         ConstitutivePassThru< DamageBase >::execute( solidModel, [&elementSubRegion, discretizationName, nodalDamage]( auto & damageModel
// )

        ConstitutivePassThru< DamageBase >::execute( solidModel, [&elementSubRegion, xNodes, discretizationName, nodalDamage]( auto & damageModel )
        {
          using CONSTITUTIVE_TYPE = TYPEOFREF( damageModel );
          typename CONSTITUTIVE_TYPE::KernelWrapper constitutiveUpdate = damageModel.createKernelUpdates();

          arrayView2d< real64 > const damageFieldOnMaterial = constitutiveUpdate.m_newDamage;
          arrayView3d< real64 > const damageGradOnMaterial = constitutiveUpdate.m_damageGrad;
          //arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemNodes = elementSubRegion.nodeList();

          auto const & elemsToNodes = elementSubRegion.nodeList().toViewConst();

          finiteElement::FiniteElementBase const &
          fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( discretizationName );
          finiteElement::FiniteElementDispatchHandler< ALL_FE_TYPES >::dispatch3D( fe, [=, &elementSubRegion] ( auto & finiteElement )
                                                                                  //  finiteElement::FiniteElementDispatchHandler<
                                                                                  //  ALL_FE_TYPES >::dispatch3D( fe, [xNodes, nodalDamage,
                                                                                  //  &elementSubRegion, damageFieldOnMaterial,
                                                                                  //  damageGradOnMaterial, elemNodes, elemsToNodes]( auto &
                                                                                  //  finiteElement )
          {
            using FE_TYPE = TYPEOFREF( finiteElement );

            DamageInterpolationKernel< FE_TYPE > interpolationKernel( elementSubRegion );

            interpolationKernel.interpolateDamage( elemsToNodes, xNodes, nodalDamage, damageFieldOnMaterial, damageGradOnMaterial );

            // forAll< serialPolicy >( elementSubRegion.size(), [xNodes, nodalDamage, damageFieldOnMaterial, damageGradOnMaterial,
            // elemNodes, elemsToNodes] ( localIndex const k )
            // {
            //   using FE_TYPE = TYPEOFREF( finiteElement );
            //   constexpr localIndex numNodesPerElement = FE_TYPE::numNodes;
            //   constexpr localIndex n_q_points = FE_TYPE::numQuadraturePoints;

            //   real64 xLocal[ numNodesPerElement ][ 3 ];
            //   real64 nodalDamageLocal[ numNodesPerElement ];

            //   for( localIndex a = 0; a < numNodesPerElement; ++a )
            //   {
            //     localIndex const localNodeIndex = elemsToNodes( k, a );

            //     for( int dim=0; dim < 3; ++dim )
            //     {
            //       xLocal[a][dim] = xNodes[ localNodeIndex ][dim];
            //     }

            //     nodalDamageLocal[ a ] = nodalDamage[ localNodeIndex ];
            //   }

            //   for( localIndex q = 0; q < n_q_points; ++q )
            //   {
            //     real64 N[ numNodesPerElement ];
            //     FE_TYPE::calcN( q, N );

            //     real64 dNdX[ numNodesPerElement ][ 3 ];

            //     real64 const detJ = FE_TYPE::calcGradN( q, xLocal, dNdX );

            //     GEOSX_UNUSED_VAR( detJ );

            //     real64 qDamage = 0.0;
            //     real64 qDamageGrad[3] = {0, 0, 0};
            //     FE_TYPE::valueAndGradient( N, dNdX, nodalDamageLocal, qDamage, qDamageGrad );

            //     damageFieldOnMaterial( k, q ) = qDamage;

            //     for( int dim=0; dim < 3; ++dim )
            //     {
            //       damageGradOnMaterial[k][q][dim] = qDamageGrad[dim];
            //     }
            //   }
            // } );
          } );
        } );
      } );
    } );
  }
}

REGISTER_CATALOG_ENTRY( SolverBase, PhaseFieldFractureSolver, string const &, Group * const )

} /* namespace geosx */
