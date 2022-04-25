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
 * @file MultiphasePoroelasticSolver.cpp
 *
 */

#include "MultiphasePoromechanicsSolver.hpp"

#include "common/DataLayouts.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/solid/PorousSolid.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "finiteElement/Kinematics.h"
#include "linearAlgebra/solvers/BlockPreconditioner.hpp"
#include "linearAlgebra/solvers/SeparateComponentPreconditioner.hpp"
#include "mesh/DomainPartition.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"
#include "constitutive/fluid/MultiFluidBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

#include "physicsSolvers/multiphysics/MultiphasePoromechanicsKernel.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

MultiphasePoromechanicsSolver::MultiphasePoromechanicsSolver( const string & name,
                                                              Group * const parent ):
  SolverBase( name, parent ),
  m_solidSolverName(),
  m_flowSolverName(),
  m_useStab( 0 )

{
  registerWrapper( viewKeyStruct::solidSolverNameString(), &m_solidSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the solid mechanics solver to use in the poroelastic solver" );

  registerWrapper( viewKeyStruct::fluidSolverNameString(), &m_flowSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the fluid mechanics solver to use in the poroelastic solver" );

  registerWrapper( viewKeyStruct::useStabFlagString(), &m_useStab ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Use pressure jump stabilization in flux" );

  m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::multiphasePoromechanics;
  m_linearSolverParameters.get().mgr.separateComponents = true;
  m_linearSolverParameters.get().mgr.displacementFieldName = keys::TotalDisplacement;
  m_linearSolverParameters.get().dofsPerNode = 3;
}

void MultiphasePoromechanicsSolver::setupDofs( DomainPartition const & domain,
                                               DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->setupDofs( domain, dofManager );
  m_flowSolver->setupDofs( domain, dofManager );

  dofManager.addCoupling( keys::TotalDisplacement,
                          CompositionalMultiphaseBase::viewKeyStruct::elemDofFieldString(),
                          DofManager::Connector::Elem );
}

void MultiphasePoromechanicsSolver::registerDataOnMesh( Group & meshBodies )
{
  SolverBase::registerDataOnMesh( meshBodies );

  forMeshTargets( meshBodies, [&] ( string const &,
                                    MeshLevel & mesh,
                                    arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                              [&]( localIndex const,
                                                                   ElementSubRegionBase & subRegion )
    {
      subRegion.registerWrapper< string >( viewKeyStruct::porousMaterialNamesString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRestartFlags( RestartFlags::NO_WRITE ).
        setSizedFromParent( 0 );

        subRegion.registerWrapper< array1d <integer> >( viewKeyStruct::elementMacroIDString() ).
        setApplyDefaultValue( -1 ).
        setPlotLevel( PlotLevel::LEVEL_1 );
    } );
  } );
}

void MultiphasePoromechanicsSolver::initializePreSubGroups()
{
  SolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elementRegionManager = mesh.getElemManager();
    elementRegionManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                                       [&]( localIndex const,
                                                                            ElementSubRegionBase & subRegion )
    {
      string & porousName = subRegion.getReference< string >( viewKeyStruct::porousMaterialNamesString() );
      porousName = getConstitutiveName< CoupledSolidBase >( subRegion );
      GEOSX_ERROR_IF( porousName.empty(), GEOSX_FMT( "Solid model not found on subregion {}", subRegion.getName() ) );
    } );
  } );
}

void MultiphasePoromechanicsSolver::initializePostInitialConditionsPreSubGroups()
{

  SolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                    MeshLevel & mesh,
                                    arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

      GEOSX_UNUSED_VAR(regionNames);

      NodeManager const & nodeManager = mesh.getNodeManager();

      ElementRegionManager::ElementViewAccessor< arrayView1d< integer > > elemMacroID =
            elemManager.constructViewAccessor< array1d< integer >, arrayView1d< integer > >( viewKeyStruct::elementMacroIDString() );

      array1d< integer > nodeVisited( nodeManager.size() );
      nodeVisited.setValues< serialPolicy>( 0 );
      arrayView1d< integer > const nodeVisitedView = nodeVisited.toView();

      arrayView1d<integer const> const bdryNodes = nodeManager.getDomainBoundaryIndicator();

      ArrayOfArraysView< localIndex const > elemRegionList = nodeManager.elementRegionList();
      ArrayOfArraysView< localIndex const > elemSubRegionList = nodeManager.elementSubRegionList();
      ArrayOfArraysView< localIndex const > elemList = nodeManager.elementList();

      integer currentID = 0;

      forAll< serialPolicy >( nodeManager.size(), [&] GEOSX_HOST_DEVICE ( localIndex const a )
      {

        if (bdryNodes[a] == 1) 
        {
          nodeVisitedView[a] = 1;
        }

        if (nodeVisitedView[a] != 1) 
        {

          for( localIndex k = 0; k < elemRegionList[a].size(); ++k )
          {
            // collect the element number
            localIndex const er = elemRegionList[a][k];
            localIndex const esr = elemSubRegionList[a][k];
            localIndex const ei = elemList[a][k];

            elemMacroID[er][esr][ei] = currentID;

            // get the elemToNodes maps 
            ElementRegionBase const & region = elemManager.getRegion( er );
            CellElementSubRegion const & subRegion = region.getSubRegion< CellElementSubRegion, localIndex >( esr );
            arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = subRegion.nodeList(); 

            // get the nodes connected to this element
            for( localIndex l = 0; l < elemsToNodes[ei].size(); ++l )
            {
              localIndex const iNode = elemsToNodes[ei][l]; // ensure this is a global index compatible with vector
              nodeVisitedView[iNode] = 1;
            }     
          }

          ++currentID;
        }

      } );


  // part 2 - assign any unassigned cell
  // loop through all elements and check if each has been assigned a macroelement (ie, when you index into ID it should be >0)
  // If not, check elements that share a face and add it to one of their macroelements
  // If no neighbors in a macroelement, reconsider this element at the end of the loop

  // This may not be strictly necessary, per our discussions. Any clumps of elements that were not visited have ID of -1, and therefore will
  // also be treated as macro elements. The main issue may be if there are a lot of single element macro elements in the mesh, but we 
  // will have to see how it works. 


 } );

}

void MultiphasePoromechanicsSolver::setupSystem( DomainPartition & domain,
                                                 DofManager & dofManager,
                                                 CRSMatrix< real64, globalIndex > & localMatrix,
                                                 ParallelVector & rhs,
                                                 ParallelVector & solution,
                                                 bool const setSparsity )
{
//  if( m_precond )
//  {
//    m_precond->clear();
//  }

  // setup monolithic coupled system
  SolverBase::setupSystem( domain, dofManager, localMatrix, rhs, solution, setSparsity );

//  if( !m_precond && m_linearSolverParameters.get().solverType != LinearSolverParameters::SolverType::direct )
//  {
//    createPreconditioner();
//  }
}

void MultiphasePoromechanicsSolver::implicitStepSetup( real64 const & time_n,
                                                       real64 const & dt,
                                                       DomainPartition & domain )
{
  m_flowSolver->implicitStepSetup( time_n, dt, domain );
  m_solidSolver->implicitStepSetup( time_n, dt, domain );
}

void MultiphasePoromechanicsSolver::implicitStepComplete( real64 const & time_n,
                                                          real64 const & dt,
                                                          DomainPartition & domain )
{
  m_solidSolver->implicitStepComplete( time_n, dt, domain );
  m_flowSolver->implicitStepComplete( time_n, dt, domain );

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                              [&]( localIndex const,
                                                                   ElementSubRegionBase & subRegion )
    {
      string const porousMaterialName = subRegion.getReference< string >( viewKeyStruct::porousMaterialNamesString() );
      ConstitutiveBase const & porousMaterial = subRegion.getConstitutiveModel< ConstitutiveBase >( porousMaterialName );
      porousMaterial.saveConvergedState();
    } );
  } );
}

void MultiphasePoromechanicsSolver::postProcessInput()
{
  SolverBase::postProcessInput();

  m_flowSolver = &this->getParent().getGroup< CompositionalMultiphaseBase >( m_flowSolverName );
  m_solidSolver = &this->getParent().getGroup< SolidMechanicsLagrangianFEM >( m_solidSolverName );
}

MultiphasePoromechanicsSolver::~MultiphasePoromechanicsSolver()
{
  // TODO Auto-generated destructor stub
}

void MultiphasePoromechanicsSolver::resetStateToBeginningOfStep( DomainPartition & domain )
{
  m_flowSolver->resetStateToBeginningOfStep( domain );
  m_solidSolver->resetStateToBeginningOfStep( domain );
}

real64 MultiphasePoromechanicsSolver::solverStep( real64 const & time_n,
                                                  real64 const & dt,
                                                  int const cycleNumber,
                                                  DomainPartition & domain )
{
  real64 dt_return = dt;

  setupSystem( domain,
               m_dofManager,
               m_localMatrix,
               m_rhs,
               m_solution );

  implicitStepSetup( time_n, dt, domain );

  dt_return = nonlinearImplicitStep( time_n, dt, cycleNumber, domain );

  implicitStepComplete( time_n, dt_return, domain );

  return dt_return;
}

void MultiphasePoromechanicsSolver::assembleSystem( real64 const time_n,
                                                    real64 const dt,
                                                    DomainPartition & domain,
                                                    DofManager const & dofManager,
                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                    arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  GEOSX_UNUSED_VAR( time_n );

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & regionNames )
  {

    NodeManager const & nodeManager = mesh.getNodeManager();

    string const displacementDofKey = dofManager.getKey( dataRepository::keys::TotalDisplacement );
    arrayView1d< globalIndex const > const & displacementDofNumber = nodeManager.getReference< globalIndex_array >( displacementDofKey );

    string const flowDofKey = dofManager.getKey( CompositionalMultiphaseBase::viewKeyStruct::elemDofFieldString() );

    localIndex const numComponents = m_flowSolver->numFluidComponents();
    localIndex const numPhases = m_flowSolver->numFluidPhases();

    real64 const gravityVectorData[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );

    poromechanicsKernels::MultiphaseKernelFactory kernelFactory( displacementDofNumber,
                                                                 flowDofKey,
                                                                 dofManager.rankOffset(),
                                                                 gravityVectorData,
                                                                 numComponents,
                                                                 numPhases,
                                                                 FlowSolverBase::viewKeyStruct::fluidNamesString(),
                                                                 localMatrix,
                                                                 localRhs );

    // Cell-based contributions
    m_solidSolver->getMaxForce() =
      finiteElement::
        regionBasedKernelApplication< parallelDevicePolicy< 32 >,
                                      constitutive::PorousSolidBase,
                                      CellElementSubRegion >( mesh,
                                                              regionNames,
                                                              this->getDiscretizationName(),
                                                              viewKeyStruct::porousMaterialNamesString(),
                                                              kernelFactory );
  } );

  // Face-based contributions
  if( m_useStab )
  {
    m_flowSolver->assembleStabilizedFluxTerms( dt,
                                               domain,
                                               dofManager,
                                               localMatrix,
                                               localRhs );
  }
  else
  {
    m_flowSolver->assembleFluxTerms( dt,
                                     domain,
                                     dofManager,
                                     localMatrix,
                                     localRhs );
  }
}

void MultiphasePoromechanicsSolver::applyBoundaryConditions( real64 const time_n,
                                                             real64 const dt,
                                                             DomainPartition & domain,
                                                             DofManager const & dofManager,
                                                             CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                             arrayView1d< real64 > const & localRhs )
{
  m_solidSolver->applyBoundaryConditions( time_n, dt,
                                          domain,
                                          dofManager,
                                          localMatrix,
                                          localRhs );

  m_flowSolver->applyBoundaryConditions( time_n, dt,
                                         domain,
                                         dofManager,
                                         localMatrix,
                                         localRhs );
}

real64 MultiphasePoromechanicsSolver::calculateResidualNorm( DomainPartition const & domain,
                                                             DofManager const & dofManager,
                                                             arrayView1d< real64 const > const & localRhs )
{
  // compute norm of momentum balance residual equations
  real64 const momementumResidualNorm = m_solidSolver->calculateResidualNorm( domain, dofManager, localRhs );

  // compute norm of mass balance residual equations
  real64 const massResidualNorm = m_flowSolver->calculateResidualNorm( domain, dofManager, localRhs );

  GEOSX_LOG_LEVEL_RANK_0( 1, GEOSX_FMT( "    ( Rsolid, Rfluid ) = ( {:4.2e}, {:4.2e} )", momementumResidualNorm, massResidualNorm ) );

  return sqrt( momementumResidualNorm * momementumResidualNorm + massResidualNorm * massResidualNorm );
}

void MultiphasePoromechanicsSolver::solveSystem( DofManager const & dofManager,
                                                 ParallelMatrix & matrix,
                                                 ParallelVector & rhs,
                                                 ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  solution.zero();
  SolverBase::solveSystem( dofManager, matrix, rhs, solution );
}

void MultiphasePoromechanicsSolver::applySystemSolution( DofManager const & dofManager,
                                                         arrayView1d< real64 const > const & localSolution,
                                                         real64 const scalingFactor,
                                                         DomainPartition & domain )
{
  // update displacement field
  m_solidSolver->applySystemSolution( dofManager, localSolution, scalingFactor, domain );
  // update pressure field
  m_flowSolver->applySystemSolution( dofManager, localSolution, -scalingFactor, domain );
}

void MultiphasePoromechanicsSolver::updateState( DomainPartition & domain )
{
  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();
    elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
    {
      m_flowSolver->updateFluidState( subRegion );
    } );

  } );
}

REGISTER_CATALOG_ENTRY( SolverBase, MultiphasePoromechanicsSolver, string const &, Group * const )

} /* namespace geosx */
