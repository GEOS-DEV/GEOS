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
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "linearAlgebra/solvers/BlockPreconditioner.hpp"
#include "linearAlgebra/solvers/SeparateComponentPreconditioner.hpp"
#include "mesh/utilities/AverageOverQuadraturePointsKernel.hpp"
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
  m_isThermal( 0 )
{
  this->registerWrapper( viewKeyStruct::isThermalString(), &m_isThermal ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag indicating whether the problem is thermal or not. Set isThermal=\"1\" to enable the thermal coupling" );

  this->registerWrapper( viewKeyStruct::performStressInitializationString(), &m_performStressInitialization ).
    setApplyDefaultValue( false ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Flag to indicate that the solver is going to perform stress initialization" );

  m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::singlePhasePoromechanics;
  m_linearSolverParameters.get().mgr.separateComponents = true;
  m_linearSolverParameters.get().mgr.displacementFieldName = solidMechanics::totalDisplacement::key();
  m_linearSolverParameters.get().dofsPerNode = 3;
}

void SinglePhasePoromechanics::registerDataOnMesh( Group & meshBodies )
{
  SolverBase::registerDataOnMesh( meshBodies );

  if( getNonlinearSolverParameters().m_couplingType == NonlinearSolverParameters::CouplingType::Sequential )
  {
    // to let the solid mechanics solver that there is a pressure and temperature RHS in the mechanics solve
    solidMechanicsSolver()->enableFixedStressPoromechanicsUpdate();
    // to let the flow solver that saving pressure_k and temperature_k is necessary (for the fixed-stress porosity terms)
    flowSolver()->enableFixedStressPoromechanicsUpdate();
  }

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

      // This is need by the way the surface generator currently does things.  
      subRegion.registerWrapper< string >( CoupledSolidBase::viewKeyStruct::porosityModelNameString() ).
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

void SinglePhasePoromechanics::setupDofs( DomainPartition const & domain,
                                          DofManager & dofManager ) const
{
  // note that the order of operations matters a lot here (for instance for the MGR labels)
  // we must set up dofs for solid mechanics first, and then for flow
  // that's the reason why this function is here and not in CoupledSolvers.hpp
  solidMechanicsSolver()->setupDofs( domain, dofManager );
  flowSolver()->setupDofs( domain, dofManager );

  setupCoupling( domain, dofManager );
}


void SinglePhasePoromechanics::initializePreSubGroups()
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
      GEOS_THROW_IF( porousName.empty(),
                     GEOS_FMT( "{} {} : Solid model not found on subregion {}",
                               catalogName(), getDataContext().toString(), subRegion.getName() ),
                     InputError );
      string & porosityModelName = subRegion.getReference< string >( CoupledSolidBase::viewKeyStruct::porosityModelNameString() );
      porosityModelName = getConstitutiveName< PorosityBase >( subRegion );
      GEOS_THROW_IF( porosityModelName.empty(),
                     GEOS_FMT( "{} {} : Porosity model not found on subregion {}",
                               catalogName(), getDataContext().toString(), subRegion.getName() ),
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

void SinglePhasePoromechanics::implicitStepSetup( real64 const & time_n,
                                                  real64 const & dt,
                                                  DomainPartition & domain )
{
  flowSolver()->keepFlowVariablesConstantDuringInitStep( m_performStressInitialization );
  Base::implicitStepSetup( time_n, dt, domain );
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
    createPreconditioner();
  }
}

void SinglePhasePoromechanics::initializePostInitialConditionsPreSubGroups()
{
  SolverBase::initializePostInitialConditionsPreSubGroups();

  arrayView1d< string const > const & poromechanicsTargetRegionNames =
    this->getReference< array1d< string > >( SolverBase::viewKeyStruct::targetRegionsString() );
  arrayView1d< string const > const & flowTargetRegionNames =
    flowSolver()->getReference< array1d< string > >( SolverBase::viewKeyStruct::targetRegionsString() );
  for( integer i = 0; i < poromechanicsTargetRegionNames.size(); ++i )
  {
    GEOS_THROW_IF( std::find( flowTargetRegionNames.begin(), flowTargetRegionNames.end(), poromechanicsTargetRegionNames[i] )
                   == flowTargetRegionNames.end(),
                   GEOS_FMT( "{} {}: region `{}` must be a target region of `{}`",
                             catalogName(), getName(), poromechanicsTargetRegionNames[i], flowSolver()->getName() ),
                   InputError );
  }

  integer & isFlowThermal = flowSolver()->getReference< integer >( FlowSolverBase::viewKeyStruct::isThermalString() );
  GEOS_LOG_RANK_0_IF( m_isThermal && !isFlowThermal,
                      GEOS_FMT( "{} {}: The attribute `{}` of the flow solver `{}` is set to 1 since the poromechanics solver is thermal",
                                catalogName(), getDataContext().toString(),
                                FlowSolverBase::viewKeyStruct::isThermalString(), flowSolver()->getName() ) );
  isFlowThermal = m_isThermal;

  if( m_isThermal )
  {
    m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::thermalSinglePhasePoromechanics;
  }
  else
  {
    if( flowSolver()->getLinearSolverParameters().mgr.strategy == LinearSolverParameters::MGR::StrategyType::singlePhaseHybridFVM )
    {
      m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::hybridSinglePhasePoromechanics;
    }
  }
}

void SinglePhasePoromechanics::assembleSystem( real64 const time_n,
                                               real64 const dt,
                                               DomainPartition & domain,
                                               DofManager const & dofManager,
                                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                               arrayView1d< real64 > const & localRhs )
{

  GEOS_MARK_FUNCTION;

  // Steps 1 and 2: compute element-based terms (mechanics and local flow terms)
  assembleElementBasedTerms( time_n,
                             dt,
                             domain,
                             dofManager,
                             localMatrix,
                             localRhs );

  // Step 3: compute the fluxes (face-based contributions)
  flowSolver()->assembleFluxTerms( time_n, dt,
                                   domain,
                                   dofManager,
                                   localMatrix,
                                   localRhs );

}

void SinglePhasePoromechanics::assembleElementBasedTerms( real64 const time_n,
                                                          real64 const dt,
                                                          DomainPartition & domain,
                                                          DofManager const & dofManager,
                                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                          arrayView1d< real64 > const & localRhs )
{
  GEOS_UNUSED_VAR( time_n );
  GEOS_UNUSED_VAR( dt );

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
                                                                                                     dt,
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
                                                                                       dt,
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
                                                                                localRhs,
                                                                                dt );
  } );

  solidMechanicsSolver()->getMaxForce() = LvArray::math::max( mechanicsMaxForce, poromechanicsMaxForce );
}

void SinglePhasePoromechanics::createPreconditioner()
{
  if( m_linearSolverParameters.get().preconditionerType == LinearSolverParameters::PreconditionerType::block )
  {
    auto precond = std::make_unique< BlockPreconditioner< LAInterface > >( BlockShapeOption::UpperTriangular,
                                                                           SchurComplementOption::RowsumDiagonalProbing,
                                                                           BlockScalingOption::FrobeniusNorm );

    auto mechPrecond = LAInterface::createPreconditioner( solidMechanicsSolver()->getLinearSolverParameters() );
    precond->setupBlock( 0,
                         { { solidMechanics::totalDisplacement::key(), { 3, true } } },
                         std::make_unique< SeparateComponentPreconditioner< LAInterface > >( 3, std::move( mechPrecond ) ) );

    auto flowPrecond = LAInterface::createPreconditioner( flowSolver()->getLinearSolverParameters() );
    precond->setupBlock( 1,
                         { { flow::pressure::key(), { 1, true } } },
                         std::move( flowPrecond ) );

    m_precond = std::move( precond );
  }
  else
  {
    //TODO: Revisit this part such that is coherent across physics solver
    //m_precond = LAInterface::createPreconditioner( m_linearSolverParameters.get() );
  }
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

  /// After the flow solver
  if( solverType == static_cast< integer >( SolverType::Flow ) )
  {
    // save pressure and temperature at the end of this iteration
    flowSolver()->saveIterationState( domain );

    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                 MeshLevel & mesh,
                                                                 arrayView1d< string const > const & regionNames )
    {

      mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                            auto & subRegion )
      {
        // update the bulk density
        // TODO: ideally, we would not recompute the bulk density, but a more general "rhs" containing the body force and the
        // pressure/temperature terms
        updateBulkDensity( subRegion );

      } );
    } );
  }

  /// After the solid mechanics solver
  if( solverType == static_cast< integer >( SolverType::SolidMechanics ) )
  {
    // compute the average of the mean stress increment over quadrature points
    averageMeanStressIncrement( domain );

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

void SinglePhasePoromechanics::averageMeanStressIncrement( DomainPartition & domain )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          auto & subRegion )
    {
      // get the solid model (to access stress increment)
      string const solidName = subRegion.template getReference< string >( viewKeyStruct::porousMaterialNamesString() );
      CoupledSolidBase & solid = getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );

      arrayView2d< real64 const > const meanStressIncrement_k = solid.getMeanEffectiveStressIncrement_k();
      arrayView1d< real64 > const averageMeanStressIncrement_k = solid.getAverageMeanEffectiveStressIncrement_k();

      finiteElement::FiniteElementBase & subRegionFE =
        subRegion.template getReference< finiteElement::FiniteElementBase >( solidMechanicsSolver()->getDiscretizationName() );

      // determine the finite element type
      finiteElement::FiniteElementDispatchHandler< BASE_FE_TYPES >::
      dispatch3D( subRegionFE, [&] ( auto const finiteElement )
      {
        using FE_TYPE = decltype( finiteElement );

        // call the factory and launch the kernel
        AverageOverQuadraturePoints1DKernelFactory::
          createAndLaunch< CellElementSubRegion,
                           FE_TYPE,
                           parallelDevicePolicy<> >( mesh.getNodeManager(),
                                                     mesh.getEdgeManager(),
                                                     mesh.getFaceManager(),
                                                     subRegion,
                                                     finiteElement,
                                                     meanStressIncrement_k,
                                                     averageMeanStressIncrement_k );
      } );
    } );
  } );
}

REGISTER_CATALOG_ENTRY( SolverBase, SinglePhasePoromechanics, string const &, Group * const )

} /* namespace geos */
