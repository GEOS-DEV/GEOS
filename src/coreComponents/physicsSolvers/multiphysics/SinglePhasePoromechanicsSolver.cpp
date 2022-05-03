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
 * @file SinglePhasePoromechanicsSolver.cpp
 */

#include "SinglePhasePoromechanicsSolver.hpp"

#include "constitutive/solid/PorousSolid.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "linearAlgebra/multiscale/MultiscalePreconditioner.hpp"
#include "linearAlgebra/solvers/BlockPreconditioner.hpp"
#include "linearAlgebra/solvers/SeparateComponentPreconditioner.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/multiphysics/SinglePhasePoromechanicsKernel.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

SinglePhasePoromechanicsSolver::SinglePhasePoromechanicsSolver( const string & name,
                                                                Group * const parent )
  : Base( name, parent ),
    m_systemScaling( 0 )
{
  registerWrapper( viewKeyStruct::linearSystemScalingString(), &m_systemScaling ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_systemScaling ).
    setDescription( "Whether block system scaling should be performed" );

  LinearSolverParameters & linParams = m_linearSolverParameters.get();
  linParams.mgr.strategy = LinearSolverParameters::MGR::StrategyType::singlePhasePoromechanics;
  linParams.mgr.separateComponents = true;
  linParams.mgr.displacementFieldName = keys::TotalDisplacement;
  linParams.dofsPerNode = 3;
  linParams.multiscale.label = "poro";
}

void SinglePhasePoromechanicsSolver::registerDataOnMesh( Group & meshBodies )
{
  SolverBase::registerDataOnMesh( meshBodies );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
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
    } );
  } );
}

void SinglePhasePoromechanicsSolver::setupCoupling( DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
                                                    DofManager & dofManager ) const
{
  dofManager.addCoupling( keys::TotalDisplacement,
                          SinglePhaseBase::viewKeyStruct::elemDofFieldString(),
                          DofManager::Connector::Elem );
}

void SinglePhasePoromechanicsSolver::initializePreSubGroups()
{
  SolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
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

void SinglePhasePoromechanicsSolver::setupSystem( DomainPartition & domain,
                                                  DofManager & dofManager,
                                                  CRSMatrix< real64, globalIndex > & localMatrix,
                                                  ParallelVector & rhs,
                                                  ParallelVector & solution,
                                                  bool const setSparsity )
{
  // setup monolithic coupled system
  SolverBase::setupSystem( domain, dofManager, localMatrix, rhs, solution, setSparsity );

  if( !m_precond && m_linearSolverParameters.get().solverType != LinearSolverParameters::SolverType::direct )
  {
    m_precond = createPreconditioner( domain );
  }
}

void SinglePhasePoromechanicsSolver::initializePostInitialConditionsPreSubGroups()
{
  using StrategyType = LinearSolverParameters::MGR::StrategyType;
  LinearSolverParameters & linParams = m_linearSolverParameters.get();
  if( flowSolver()->getLinearSolverParameters().mgr.strategy == StrategyType::singlePhaseHybridFVM )
  {
    linParams.mgr.strategy = StrategyType::hybridSinglePhasePoromechanics;
  }
  linParams.block.subParams.emplace_back( &solidMechanicsSolver()->getLinearSolverParameters() );
  linParams.block.subParams.emplace_back( &flowSolver()->getLinearSolverParameters() );
}

real64 SinglePhasePoromechanicsSolver::solverStep( real64 const & time_n,
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

void SinglePhasePoromechanicsSolver::assembleSystem( real64 const time_n,
                                                     real64 const dt,
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
    NodeManager const & nodeManager = mesh.getNodeManager();

    string const dofKey = dofManager.getKey( dataRepository::keys::TotalDisplacement );
    arrayView1d< globalIndex const > const & dispDofNumber = nodeManager.getReference< globalIndex_array >( dofKey );

    string const pDofKey = dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );

    real64 const gravityVectorData[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );

    poromechanicsKernels::SinglePhaseKernelFactory kernelFactory( dispDofNumber,
                                                                  pDofKey,
                                                                  dofManager.rankOffset(),
                                                                  localMatrix,
                                                                  localRhs,
                                                                  gravityVectorData,
                                                                  FlowSolverBase::viewKeyStruct::fluidNamesString() );

    // Cell-based contributions
    solidMechanicsSolver()->getMaxForce() =
      finiteElement::
        regionBasedKernelApplication< parallelDevicePolicy< 32 >,
                                      constitutive::PorousSolidBase,
                                      CellElementSubRegion >( mesh,
                                                              regionNames,
                                                              solidMechanicsSolver()->getDiscretizationName(),
                                                              viewKeyStruct::porousMaterialNamesString(),
                                                              kernelFactory );

  } );

  flowSolver()->assemblePoroelasticFluxTerms( time_n, dt,
                                              domain,
                                              dofManager,
                                              localMatrix,
                                              localRhs,
                                              " " );
}

std::unique_ptr< PreconditionerBase< LAInterface> >
SinglePhasePoromechanicsSolver::createPreconditioner( DomainPartition & domain ) const
{
  LinearSolverParameters const & linParams = m_linearSolverParameters.get();
  switch( linParams.preconditionerType )
  {
    case LinearSolverParameters::PreconditionerType::block:
    {
      auto precond = std::make_unique< BlockPreconditioner< LAInterface > >( linParams.block );

      precond->setupBlock( 0,
                           { { keys::TotalDisplacement, { 3, true } } },
                           solidMechanicsSolver()->createPreconditioner( domain ) );
      precond->setupBlock( 1,
                           { { SinglePhaseBase::viewKeyStruct::elemDofFieldString(), { 1, true } } },
                           flowSolver()->createPreconditioner( domain ) );

      return precond;
    }
    case LinearSolverParameters::PreconditionerType::multiscale:
    {
      return std::make_unique< MultiscalePreconditioner< LAInterface > >( linParams, domain );
    }
    default:
    {
      return SolverBase::createPreconditioner( domain );
    }
  }
}

void SinglePhasePoromechanicsSolver::solveLinearSystem( DofManager const & dofManager,
                                                        ParallelMatrix & matrix,
                                                        ParallelVector & rhs,
                                                        ParallelVector & solution )
{
  if( m_systemScaling )
  {
    // Only compute this once and reuse for the entire simulation
    if( !m_scalingVector.created() )
    {
      // TODO: currently only handles displacement and cell pressure blocks, ignores face pressure in HybridFVM
      DofManager::SubComponent const fields[2] =
      {
        { keys::TotalDisplacement, DofManager::CompMask{ 3, true } },
        { SinglePhaseBase::viewKeyStruct::elemDofFieldString(), DofManager::CompMask{ 1, true } }
      };

      real64 norms[2];
      for( integer i = 0; i < 2; ++i )
      {
        ParallelMatrix P, A;
        dofManager.makeRestrictor( { fields[i] }, matrix.comm(), true, P );
        matrix.multiplyPtAP( P, A );
        norms[i] = A.normFrobenius();
      }
      real64 const scale[2] = { std::min( norms[1] / norms[0], 1.0 ), std::min( norms[0] / norms[1], 1.0 ) };

      m_scalingVector.create( rhs.localSize(), rhs.comm() );
      m_scalingVector.set( 1.0 );

      localIndex offset = 0;
      arrayView1d< real64 > const values = m_scalingVector.open();
      for( integer i = 0; i < 2; ++i )
      {
        localIndex const numDof = dofManager.numLocalDofs( fields[i].fieldName );
        forAll< parallelDevicePolicy<> >( numDof, [=] GEOSX_HOST_DEVICE ( localIndex const k )
        {
          values[offset + k] = scale[i];
        } );
        offset += numDof;
      }
      m_scalingVector.close();
    }

    matrix.leftScale( m_scalingVector );
    rhs.pointwiseProduct( m_scalingVector, rhs );
  }

  SolverBase::solveLinearSystem( dofManager, matrix, rhs, solution );
}

void SinglePhasePoromechanicsSolver::updateState( DomainPartition & domain )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   CellElementSubRegion & subRegion )
    {
      flowSolver()->updateFluidState( subRegion );

    } );
  } );
}

REGISTER_CATALOG_ENTRY( SolverBase, SinglePhasePoromechanicsSolver, string const &, Group * const )

} /* namespace geosx */
