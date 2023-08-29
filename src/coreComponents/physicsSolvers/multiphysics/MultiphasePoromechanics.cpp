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
 * @file MultiphasePoromechanics.cpp
 */

#define GEOSX_DISPATCH_VEM /// enables VEM in FiniteElementDispatch

#include "MultiphasePoromechanics.hpp"

#include "constitutive/fluid/MultiFluidBase.hpp"
#include "constitutive/solid/PorousSolid.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/MultiphasePoromechanics.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/ThermalMultiphasePoromechanics.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/solidMechanics/kernels/ImplicitSmallStrainQuasiStatic.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace fields;

MultiphasePoromechanics::MultiphasePoromechanics( const string & name,
                                                  Group * const parent )
  : Base( name, parent ),
  m_isThermal( 0 )
{
  registerWrapper( viewKeyStruct::stabilizationTypeString(), &m_stabilizationType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Stabilization type. Options are:\n" +
                    toString( StabilizationType::None ) + " - Add no stabilization to mass equation,\n" +
                    toString( StabilizationType::Global ) + " - Add stabilization to all faces,\n" +
                    toString( StabilizationType::Local ) + " - Add stabilization only to interiors of macro elements." );

  registerWrapper( viewKeyStruct::stabilizationRegionNamesString(), &m_stabilizationRegionNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Regions where stabilization is applied." );

  registerWrapper( viewKeyStruct::stabilizationMultiplierString(), &m_stabilizationMultiplier ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Constant multiplier of stabilization strength." );

  registerWrapper( viewKeyStruct::isThermalString(), &m_isThermal ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag indicating whether the problem is thermal or not. Set isThermal=\"1\" to enable the thermal coupling" );

  registerWrapper( viewKeyStruct::performStressInitializationString(), &m_performStressInitialization ).
    setApplyDefaultValue( false ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Flag to indicate that the solver is going to perform stress initialization" );

  m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::multiphasePoromechanics;
  m_linearSolverParameters.get().mgr.separateComponents = true;
  m_linearSolverParameters.get().mgr.displacementFieldName = solidMechanics::totalDisplacement::key();
  m_linearSolverParameters.get().dofsPerNode = 3;
}

void MultiphasePoromechanics::registerDataOnMesh( Group & meshBodies )
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

      if( m_stabilizationType == StabilizationType::Global ||
          m_stabilizationType == StabilizationType::Local )
      {
        subRegion.registerField< fields::flow::macroElementIndex >( getName() );
        subRegion.registerField< fields::flow::elementStabConstant >( getName() );
      }
    } );
  } );
}

void MultiphasePoromechanics::setupCoupling( DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
                                             DofManager & dofManager ) const
{
  dofManager.addCoupling( solidMechanics::totalDisplacement::key(),
                          CompositionalMultiphaseBase::viewKeyStruct::elemDofFieldString(),
                          DofManager::Connector::Elem );
}

void MultiphasePoromechanics::assembleSystem( real64 const GEOSX_UNUSED_PARAM( time ),
                                              real64 const dt,
                                              DomainPartition & domain,
                                              DofManager const & dofManager,
                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                              arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  real64 poromechanicsMaxForce = 0.0;
  real64 mechanicsMaxForce = 0.0;

  // step 1: apply the full poromechanics coupling on the target regions on the poromechanics solver

  set< string > poromechanicsRegionNames;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    poromechanicsRegionNames.insert( regionNames.begin(), regionNames.end() );

    string const flowDofKey = dofManager.getKey( CompositionalMultiphaseBase::viewKeyStruct::elemDofFieldString() );

    if( m_isThermal )
    {
      poromechanicsMaxForce =
        assemblyLaunch< constitutive::PorousSolid< ElasticIsotropic >, // TODO: change once there is a cmake solution
                        thermalPoromechanicsKernels::ThermalMultiphasePoromechanicsKernelFactory >( mesh,
                                                                                                    dofManager,
                                                                                                    regionNames,
                                                                                                    viewKeyStruct::porousMaterialNamesString(),
                                                                                                    localMatrix,
                                                                                                    localRhs,
                                                                                                    flowDofKey,
                                                                                                    flowSolver()->numFluidComponents(),
                                                                                                    flowSolver()->numFluidPhases(),
                                                                                                    FlowSolverBase::viewKeyStruct::fluidNamesString() );
    }
    else
    {
      poromechanicsMaxForce =
        assemblyLaunch< constitutive::PorousSolidBase,
                        poromechanicsKernels::MultiphasePoromechanicsKernelFactory >( mesh,
                                                                                      dofManager,
                                                                                      regionNames,
                                                                                      viewKeyStruct::porousMaterialNamesString(),
                                                                                      localMatrix,
                                                                                      localRhs,
                                                                                      flowDofKey,
                                                                                      flowSolver()->numFluidComponents(),
                                                                                      flowSolver()->numFluidPhases(),
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

  if( m_stabilizationType == StabilizationType::Global ||
      m_stabilizationType == StabilizationType::Local )
  {
    updateStabilizationParameters( domain );
    flowSolver()->assembleStabilizedFluxTerms( dt,
                                               domain,
                                               dofManager,
                                               localMatrix,
                                               localRhs );
  }
  else
  {
    flowSolver()->assembleFluxTerms( dt,
                                     domain,
                                     dofManager,
                                     localMatrix,
                                     localRhs );
  }
}

void MultiphasePoromechanics::updateState( DomainPartition & domain )
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

void MultiphasePoromechanics::initializePostInitialConditionsPreSubGroups()
{
  SolverBase::initializePostInitialConditionsPreSubGroups();

  integer & isFlowThermal = flowSolver()->getReference< integer >( FlowSolverBase::viewKeyStruct::isThermalString() );
  GEOSX_LOG_RANK_0_IF( m_isThermal && !isFlowThermal,
                       GEOSX_FMT( "{} {}: The attribute `{}` of the flow solver `{}` is set to 1 since the poromechanics solver is thermal",
                                  catalogName(), getName(), FlowSolverBase::viewKeyStruct::isThermalString(), flowSolver()->getName() ) );
  isFlowThermal = m_isThermal;

  if( m_isThermal )
  {
    m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::thermalMultiphasePoromechanics;
  }
}

void MultiphasePoromechanics::initializePreSubGroups()
{
  SolverBase::initializePreSubGroups();

  if( getNonlinearSolverParameters().m_couplingType == NonlinearSolverParameters::CouplingType::Sequential )
  {
    solidMechanicsSolver()->turnOnFixedStressThermoPoroElasticityFlag();
  }

  GEOSX_THROW_IF( m_stabilizationType == StabilizationType::Local,
                  catalogName() << " " << getName() << ": Local stabilization has been disabled temporarily",
                  InputError );

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

void MultiphasePoromechanics::resetStateToBeginningOfStep( DomainPartition & domain )
{
  Base::resetStateToBeginningOfStep( domain ); 

  flowSolver()->keepFlowVariablesConstantDuringInitStep( m_performStressInitialization );
}  

void MultiphasePoromechanics::updateStabilizationParameters( DomainPartition & domain ) const
{
  // Step 1: we loop over the regions where stabilization is active and collect their name

  set< string > regionFilter;
  for( string const & regionName : m_stabilizationRegionNames )
  {
    regionFilter.insert( regionName );
  }

  // Step 2: loop over the target regions of the solver, and tag the elements belonging to stabilization regions
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & targetRegionNames )
  {
    // keep only the target regions that are in the filter
    array1d< string > filteredTargetRegionNames;
    filteredTargetRegionNames.reserve( targetRegionNames.size() );

    for( string const & targetRegionName : targetRegionNames )
    {
      if( regionFilter.count( targetRegionName ) )
      {
        filteredTargetRegionNames.emplace_back( targetRegionName );
      }
    }

    // loop over the elements and update the stabilization constant
    mesh.getElemManager().forElementSubRegions( filteredTargetRegionNames.toViewConst(), [&]( localIndex const,
                                                                                              ElementSubRegionBase & subRegion )

    {
      arrayView1d< integer > const macroElementIndex = subRegion.getField< fields::flow::macroElementIndex >();
      arrayView1d< real64 > const elementStabConstant = subRegion.getField< fields::flow::elementStabConstant >();

      geosx::constitutive::CoupledSolidBase const & porousSolid =
        getConstitutiveModel< geosx::constitutive::CoupledSolidBase >( subRegion, subRegion.getReference< string >( viewKeyStruct::porousMaterialNamesString() ) );

      arrayView1d< real64 const > const bulkModulus = porousSolid.getBulkModulus();
      arrayView1d< real64 const > const shearModulus = porousSolid.getShearModulus();
      arrayView1d< real64 const > const biotCoefficient = porousSolid.getBiotCoefficient();

      real64 const stabilizationMultiplier = m_stabilizationMultiplier;

      forAll< parallelDevicePolicy<> >( subRegion.size(), [bulkModulus,
                                                           shearModulus,
                                                           biotCoefficient,
                                                           stabilizationMultiplier,
                                                           macroElementIndex,
                                                           elementStabConstant] GEOSX_HOST_DEVICE ( localIndex const ei )
      {
        real64 const bM = bulkModulus[ei];
        real64 const sM = shearModulus[ei];
        real64 const bC = biotCoefficient[ei];

        macroElementIndex[ei] = 1;
        elementStabConstant[ei] = stabilizationMultiplier * 9.0 * (bC * bC) / (32.0 * (10.0 * sM / 3.0 + bM));

      } );
    } );
  } );
}

void MultiphasePoromechanics::mapSolutionBetweenSolvers( DomainPartition & domain, integer const solverType )
{
  GEOSX_MARK_FUNCTION;
  if( solverType == static_cast< integer >( SolverType::SolidMechanics ) )
  {
    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                 MeshLevel & mesh,
                                                                 arrayView1d< string const > const & regionNames )
    {
      mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                            auto & subRegion )
      {
        flowSolver()->updatePorosityAndPermeability( subRegion );
      } );
    } );
  }
}



REGISTER_CATALOG_ENTRY( SolverBase, MultiphasePoromechanics, string const &, Group * const )

} /* namespace geosx */
