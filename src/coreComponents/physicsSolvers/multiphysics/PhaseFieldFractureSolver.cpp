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

    // if( m_subdomainFlag == 1 )
    // {
    //   this->solidMechanicsSolver()->setSubdomainFlag();
    //   this->damageSolver()->setSubdomainFlag();
    //   this->solidMechanicsSolver()->setSubdomainElemSet(m_subdomainElems);
    //   this->damageSolver()->setSubdomainElemSet(m_subdomainElems);      
    // }
    //imposeFakeBackgroundPressures( domain ); //DOMAIN IS COMING OUT EMPTY
  }

  GEOSX_WARNING_IF( getNonlinearSolverParameters().m_couplingType == NonlinearSolverParameters::CouplingType::FullyImplicit,
                    "FullyImplicit coupling not implemented for this solver. A sequential coupling approach will be used." );
  getNonlinearSolverParameters().m_couplingType = NonlinearSolverParameters::CouplingType::Sequential;
}

void PhaseFieldFractureSolver::mapSolutionBetweenSolvers( DomainPartition & domain, integer const solverType )
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
