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
 * @file SinglePhasePoromechanics.cpp
 */

#define GEOSX_DISPATCH_VEM /// enables VEM in FiniteElementDispatch

#include "SinglePhasePoromechanics.hpp"

#include "constitutive/solid/PorousSolid.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "linearAlgebra/solvers/BlockPreconditioner.hpp"
#include "linearAlgebra/solvers/SeparateComponentPreconditioner.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/SinglePhasePoromechanics.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/ThermalSinglePhasePoromechanics.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/solidMechanics/kernels/ImplicitSmallStrainQuasiStatic.hpp"

namespace geos
{

using namespace constitutive;
using namespace dataRepository;
using namespace fields;

SinglePhasePoromechanics::SinglePhasePoromechanics( const string & name,
                                                    Group * const parent )
  : Base( name, parent ),
  m_isThermal( 0 ),
  m_systemScaling( 0 )
{
  this->registerWrapper( viewKeyStruct::isThermalString(), &m_isThermal ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag indicating whether the problem is thermal or not. Set isThermal=\"1\" to enable the thermal coupling" );

  this->registerWrapper( viewKeyStruct::performStressInitializationString(), &m_performStressInitialization ).
    setApplyDefaultValue( false ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Flag to indicate that the solver is going to perform stress initialization" );

  registerWrapper( viewKeyStruct::linearSystemScalingString(), &m_systemScaling ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( m_systemScaling ).
    setDescription( "Whether block system scaling should be performed" );

  LinearSolverParameters & linParams = m_linearSolverParameters.get();
  linParams.mgr.strategy = LinearSolverParameters::MGR::StrategyType::singlePhasePoromechanics;
  linParams.mgr.separateComponents = true;
  linParams.mgr.displacementFieldName = solidMechanics::totalDisplacement::key();
  linParams.dofsPerNode = 3;
}

void SinglePhasePoromechanics::registerDataOnMesh( Group & meshBodies )
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

      if( getNonlinearSolverParameters().m_couplingType == NonlinearSolverParameters::CouplingType::Sequential )
      {
        // register the bulk density for use in the solid mechanics solver
        // ideally we would resize it here as well, but the solid model name is not available yet (see below)
        subRegion.registerField< fields::poromechanics::bulkDensity >( getName() );
      }

    } );
  } );
}

void SinglePhasePoromechanics::setupCoupling( DomainPartition const & GEOS_UNUSED_PARAM( domain ),
                                              DofManager & dofManager ) const
{
  dofManager.addCoupling( solidMechanics::totalDisplacement::key(),
                          SinglePhaseBase::viewKeyStruct::elemDofFieldString(),
                          DofManager::Connector::Elem );
}

void SinglePhasePoromechanics::initializePreSubGroups()
{
  SolverBase::initializePreSubGroups();

  if( getNonlinearSolverParameters().m_couplingType == NonlinearSolverParameters::CouplingType::Sequential )
  {
    solidMechanicsSolver()->turnOnFixedStressThermoPoromechanicsFlag();
  }

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
      GEOS_THROW_IF( porousName.empty(),
                     GEOS_FMT( "{} {} : Solid model not found on subregion {}", catalogName(), getName(), subRegion.getName() ),
                     InputError );

      if( subRegion.hasField< fields::poromechanics::bulkDensity >() )
      {
        // get the solid model to know the number of quadrature points and resize the bulk density
        CoupledSolidBase const & solid = getConstitutiveModel< CoupledSolidBase >( subRegion, porousName );
        subRegion.getField< fields::poromechanics::bulkDensity >().resizeDimension< 1 >( solid.getDensity().size( 1 ) );
      }

    } );
  } );

}

void SinglePhasePoromechanics::setupSystem( DomainPartition & domain,
                                            DofManager & dofManager,
                                            CRSMatrix< real64, globalIndex > & localMatrix,
                                            ParallelVector & rhs,
                                            ParallelVector & solution,
                                            bool const setSparsity )
{
  if( m_precond )
  {
    m_precond->clear();
  }

  // setup monolithic coupled system
  SolverBase::setupSystem( domain, dofManager, localMatrix, rhs, solution, setSparsity );

  if( !m_precond && m_linearSolverParameters.get().solverType != LinearSolverParameters::SolverType::direct )
  {
    m_precond = createPreconditioner( domain );
  }
}

void SinglePhasePoromechanics::initializePostInitialConditionsPreSubGroups()
{

  SolverBase::initializePostInitialConditionsPreSubGroups();

  integer & isFlowThermal = flowSolver()->getReference< integer >( FlowSolverBase::viewKeyStruct::isThermalString() );
  GEOS_LOG_RANK_0_IF( m_isThermal && !isFlowThermal,
                      GEOS_FMT( "{} {}: The attribute `{}` of the flow solver `{}` is set to 1 since the poromechanics solver is thermal",
                                catalogName(), getName(), FlowSolverBase::viewKeyStruct::isThermalString(), flowSolver()->getName() ) );
  isFlowThermal = m_isThermal;

  using StrategyType = LinearSolverParameters::MGR::StrategyType;
  LinearSolverParameters & linParams = m_linearSolverParameters.get();
  if( m_isThermal )
  {
    linParams.mgr.strategy = StrategyType::thermalSinglePhasePoromechanics;
  }
  else
  {
    if( flowSolver()->getLinearSolverParameters().mgr.strategy == StrategyType::singlePhaseHybridFVM )
    {
      linParams.mgr.strategy = StrategyType::hybridSinglePhasePoromechanics;
    }
  }

  // Populate sub-block solver parameters for block preconditioner
  linParams.block.resize( 2 );
  linParams.block.subParams[toUnderlying( SolverType::SolidMechanics )] = &solidMechanicsSolver()->getLinearSolverParameters();
  linParams.block.subParams[toUnderlying( SolverType::Flow )] = &flowSolver()->getLinearSolverParameters();

  // For classical fixed-stress scheme the order must be mechanics-flow
  linParams.block.order[toUnderlying( SolverType::SolidMechanics )] = 0;
  linParams.block.order[toUnderlying( SolverType::Flow )] = 1;
}

void SinglePhasePoromechanics::assembleSystem( real64 const time_n,
                                               real64 const dt,
                                               DomainPartition & domain,
                                               DofManager const & dofManager,
                                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                               arrayView1d< real64 > const & localRhs )
{

  GEOS_MARK_FUNCTION;

  real64 poromechanicsMaxForce = 0.0;
  real64 mechanicsMaxForce = 0.0;

  // step 1: apply the full poromechanics coupling on the target regions on the poromechanics solver

  set< string > poromechanicsRegionNames;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    poromechanicsRegionNames.insert( regionNames.begin(), regionNames.end() );

    string const flowDofKey = dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );

    if( m_isThermal )
    {
      poromechanicsMaxForce =
        assemblyLaunch< constitutive::PorousSolid< ElasticIsotropic >, // TODO: change once there is a cmake solution
                        thermalPoromechanicsKernels::ThermalSinglePhasePoromechanicsKernelFactory >( mesh,
                                                                                                     dofManager,
                                                                                                     regionNames,
                                                                                                     viewKeyStruct::porousMaterialNamesString(),
                                                                                                     localMatrix,
                                                                                                     localRhs,
                                                                                                     flowDofKey,
                                                                                                     FlowSolverBase::viewKeyStruct::fluidNamesString() );
    }
    else
    {
      poromechanicsMaxForce =
        assemblyLaunch< constitutive::PorousSolidBase,
                        poromechanicsKernels::SinglePhasePoromechanicsKernelFactory >( mesh,
                                                                                       dofManager,
                                                                                       regionNames,
                                                                                       viewKeyStruct::porousMaterialNamesString(),
                                                                                       localMatrix,
                                                                                       localRhs,
                                                                                       flowDofKey,
                                                                                       FlowSolverBase::viewKeyStruct::fluidNamesString() );
    }
  } );

  // step 2: apply mechanics solver on its target regions not included in the poromechanics solver target regions

  solidMechanicsSolver()->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                                        MeshLevel & mesh,
                                                                                        arrayView1d< string const > const & regionNames )
  {
    // collect the target region of the mechanics solver not included in the poromechanics target regions
    array1d< string > filteredRegionNames;
    filteredRegionNames.reserve( regionNames.size() );
    for( string const & regionName : regionNames )
    {
      // if the mechanics target region is not included in the poromechanics target region, save the string
      if( poromechanicsRegionNames.count( regionName ) == 0 )
      {
        filteredRegionNames.emplace_back( regionName );
      }
    }

    // if the array is empty, the mechanics and poromechanics solver target regions overlap perfectly, there is nothing to do
    if( filteredRegionNames.empty() )
    {
      return;
    }

    mechanicsMaxForce =
      assemblyLaunch< constitutive::SolidBase,
                      solidMechanicsLagrangianFEMKernels::QuasiStaticFactory >( mesh,
                                                                                dofManager,
                                                                                filteredRegionNames.toViewConst(),
                                                                                SolidMechanicsLagrangianFEM::viewKeyStruct::solidMaterialNamesString(),
                                                                                localMatrix,
                                                                                localRhs );
  } );

  solidMechanicsSolver()->getMaxForce() = LvArray::math::max( mechanicsMaxForce, poromechanicsMaxForce );


  // tell the flow solver that this is a stress initialization step
  flowSolver()->keepFlowVariablesConstantDuringInitStep( m_performStressInitialization );

  // step 3: compute the fluxes (face-based contributions)

  if( m_isThermal )
  {
    flowSolver()->assembleFluxTerms( time_n, dt,
                                     domain,
                                     dofManager,
                                     localMatrix,
                                     localRhs );
  }
  else
  {
    flowSolver()->assemblePoroelasticFluxTerms( time_n, dt,
                                                domain,
                                                dofManager,
                                                localMatrix,
                                                localRhs,
                                                " " );
  }

}

std::unique_ptr< PreconditionerBase< LAInterface > >
SinglePhasePoromechanics::createPreconditioner( DomainPartition & domain ) const
{
  LinearSolverParameters const & linParams = m_linearSolverParameters.get();
  switch( linParams.preconditionerType )
  {
    case LinearSolverParameters::PreconditionerType::block:
    {
      auto precond = std::make_unique< BlockPreconditioner< LAInterface > >( linParams.block );

      precond->setupBlock( linParams.block.order[toUnderlying( SolverType::SolidMechanics )],
                           { { solidMechanics::totalDisplacement::key(), { 3, true } } },
                           solidMechanicsSolver()->createPreconditioner( domain ) );
      precond->setupBlock( linParams.block.order[toUnderlying( SolverType::Flow )],
                           { { SinglePhaseBase::viewKeyStruct::elemDofFieldString(), { 1, true } } },
                           flowSolver()->createPreconditioner( domain ) );

      return precond;
    }
    default:
    {
      return SolverBase::createPreconditioner( domain );
    }
  }
}

void SinglePhasePoromechanics::solveLinearSystem( DofManager const & dofManager,
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
      DofManager::SubComponent fields[2];
      fields[toUnderlying( SolverType::SolidMechanics )] = { solidMechanics::totalDisplacement::key(), DofManager::CompMask{ 3, true } };
      fields[toUnderlying( SolverType::Flow )] = { SinglePhaseBase::viewKeyStruct::elemDofFieldString(), DofManager::CompMask{ 1, true } };

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
        forAll< parallelDevicePolicy<> >( numDof, [=] GEOS_HOST_DEVICE ( localIndex const k )
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

void SinglePhasePoromechanics::updateState( DomainPartition & domain )
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
      if( m_isThermal )
      {
        flowSolver()->updateSolidInternalEnergyModel( subRegion );
      }
    } );
  } );
}

void SinglePhasePoromechanics::mapSolutionBetweenSolvers( DomainPartition & domain, integer const solverType )
{
  GEOS_MARK_FUNCTION;
  if( solverType == static_cast< integer >( SolverType::SolidMechanics ) )
  {
    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                 MeshLevel & mesh,
                                                                 arrayView1d< string const > const & regionNames )
    {
      mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                            auto & subRegion )
      {
        // update the porosity after a change in displacement (after mechanics solve)
        // or a change in pressure/temperature (after a flow solve)
        flowSolver()->updatePorosityAndPermeability( subRegion );

        // update the bulk density
        // in fact, this is only needed after a flow solve, but we don't have a mechanism to know where we are in the outer loop
        // TODO: ideally, we would not recompute the bulk density, but a more general "rhs" containing the body force and the
        // pressure/temperature terms
        updateBulkDensity( subRegion );
      } );
    } );
  }
}

void SinglePhasePoromechanics::updateBulkDensity( ElementSubRegionBase & subRegion )
{
  // get the fluid model (to access fluid density)
  string const fluidName = subRegion.getReference< string >( FlowSolverBase::viewKeyStruct::fluidNamesString() );
  SingleFluidBase const & fluid = getConstitutiveModel< SingleFluidBase >( subRegion, fluidName );

  // get the solid model (to access porosity and solid density)
  string const solidName = subRegion.getReference< string >( viewKeyStruct::porousMaterialNamesString() );
  CoupledSolidBase const & solid = getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );

  // update the bulk density
  poromechanicsKernels::
    SinglePhaseBulkDensityKernelFactory::
    createAndLaunch< parallelDevicePolicy<> >( fluid,
                                               solid,
                                               subRegion );
}

REGISTER_CATALOG_ENTRY( SolverBase, SinglePhasePoromechanics, string const &, Group * const )

} /* namespace geos */
