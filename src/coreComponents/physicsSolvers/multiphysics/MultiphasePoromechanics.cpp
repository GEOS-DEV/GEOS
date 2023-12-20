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

#include "constitutive/fluid/multifluid/MultiFluidBase.hpp"
#include "constitutive/solid/PorousSolid.hpp"
#include "mesh/utilities/AverageOverQuadraturePointsKernel.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/multiphysics/CompositionalMultiphaseReservoirAndWells.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/MultiphasePoromechanics.hpp"
#include "physicsSolvers/multiphysics/poromechanicsKernels/ThermalMultiphasePoromechanics.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/solidMechanics/kernels/ImplicitSmallStrainQuasiStatic.hpp"

namespace geos
{

using namespace dataRepository;
using namespace constitutive;
using namespace fields;
using namespace stabilization;

namespace
{

// This is meant to be specialized to work, see below
template< typename FLOW_SOLVER > class
  MultiphaseCatalogNames {};
// Class specialization for a FLOW_SOLVER set to CompositionalMultiphaseBase
template<> class MultiphaseCatalogNames< CompositionalMultiphaseBase >
{
public:
  static string name() { return "MultiphasePoromechanics"; }
};
// Class specialization for a FLOW_SOLVER set to CompositionalMultiphaseReservoirAndWells
template<> class MultiphaseCatalogNames< CompositionalMultiphaseReservoirAndWells< CompositionalMultiphaseBase > >
{
public:
  static string name() { return CompositionalMultiphaseReservoirAndWells< CompositionalMultiphaseBase >::catalogName() + "Poromechanics"; }
};
}

// provide a definition for catalogName()
template< typename FLOW_SOLVER >
string
MultiphasePoromechanics< FLOW_SOLVER >::
catalogName()
{
  return MultiphaseCatalogNames< FLOW_SOLVER >::name();
}

template< typename FLOW_SOLVER >
MultiphasePoromechanics< FLOW_SOLVER >::MultiphasePoromechanics( const string & name,
                                                                 Group * const parent )
  : Base( name, parent ),
  m_isThermal( 0 )
{
  this->registerWrapper( viewKeyStruct::stabilizationTypeString(), &m_stabilizationType ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Stabilization type. Options are:\n" +
                    toString( StabilizationType::None ) + " - Add no stabilization to mass equation,\n" +
                    toString( StabilizationType::Global ) + " - Add stabilization to all faces,\n" +
                    toString( StabilizationType::Local ) + " - Add stabilization only to interiors of macro elements." );

  this->registerWrapper( viewKeyStruct::stabilizationRegionNamesString(), &m_stabilizationRegionNames ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Regions where stabilization is applied." );

  this->registerWrapper( viewKeyStruct::stabilizationMultiplierString(), &m_stabilizationMultiplier ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Constant multiplier of stabilization strength." );

  this->registerWrapper( viewKeyStruct::isThermalString(), &m_isThermal ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag indicating whether the problem is thermal or not. Set isThermal=\"1\" to enable the thermal coupling" );

  this->registerWrapper( viewKeyStruct::performStressInitializationString(), &m_performStressInitialization ).
    setApplyDefaultValue( false ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "Flag to indicate that the solver is going to perform stress initialization" );

  LinearSolverParameters & linearSolverParameters = this->m_linearSolverParameters.get();
  linearSolverParameters.mgr.strategy = LinearSolverParameters::MGR::StrategyType::multiphasePoromechanics;
  linearSolverParameters.mgr.separateComponents = true;
  linearSolverParameters.mgr.displacementFieldName = solidMechanics::totalDisplacement::key();
  linearSolverParameters.dofsPerNode = 3;
}

template< typename FLOW_SOLVER >
void MultiphasePoromechanics< FLOW_SOLVER >::postProcessInput()
{
  Base::postProcessInput();

  GEOS_ERROR_IF( flowSolver()->catalogName() == "CompositionalMultiphaseReservoir" &&
                 this->getNonlinearSolverParameters().couplingType() != NonlinearSolverParameters::CouplingType::Sequential,
                 GEOS_FMT( "{}: {} solver is only designed to work for {} = {}",
                           this->getDataContext(), catalogName(), NonlinearSolverParameters::viewKeysStruct::couplingTypeString(),
                           EnumStrings< NonlinearSolverParameters::CouplingType >::toString( NonlinearSolverParameters::CouplingType::Sequential )
                           ));
}

template< typename FLOW_SOLVER >
void MultiphasePoromechanics< FLOW_SOLVER >::registerDataOnMesh( Group & meshBodies )
{
  SolverBase::registerDataOnMesh( meshBodies );

  if( this->getNonlinearSolverParameters().m_couplingType == NonlinearSolverParameters::CouplingType::Sequential )
  {
    // to let the solid mechanics solver that there is a pressure and temperature RHS in the mechanics solve
    solidMechanicsSolver()->enableFixedStressPoromechanicsUpdate();
    // to let the flow solver that saving pressure_k and temperature_k is necessary (for the fixed-stress porosity terms)
    flowSolver()->enableFixedStressPoromechanicsUpdate();
  }

  this->template forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
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
        subRegion.registerField< fields::flow::macroElementIndex >( this->getName() );
        subRegion.registerField< fields::flow::elementStabConstant >( this->getName() );
      }

      if( this->getNonlinearSolverParameters().m_couplingType == NonlinearSolverParameters::CouplingType::Sequential )
      {
        // register the bulk density for use in the solid mechanics solver
        // ideally we would resize it here as well, but the solid model name is not available yet (see below)
        subRegion.registerField< fields::poromechanics::bulkDensity >( this->getName() );
      }
    } );
  } );
}

template< typename FLOW_SOLVER >
void MultiphasePoromechanics< FLOW_SOLVER >::setupCoupling( DomainPartition const & GEOS_UNUSED_PARAM( domain ),
                                                            DofManager & dofManager ) const
{
  dofManager.addCoupling( solidMechanics::totalDisplacement::key(),
                          CompositionalMultiphaseBase::viewKeyStruct::elemDofFieldString(),
                          DofManager::Connector::Elem );
}

template< typename FLOW_SOLVER >
void MultiphasePoromechanics< FLOW_SOLVER >::setupDofs( DomainPartition const & domain,
                                                        DofManager & dofManager ) const
{
  // note that the order of operations matters a lot here (for instance for the MGR labels)
  // we must set up dofs for solid mechanics first, and then for flow
  // that's the reason why this function is here and not in CoupledSolvers.hpp
  solidMechanicsSolver()->setupDofs( domain, dofManager );
  flowSolver()->setupDofs( domain, dofManager );

  setupCoupling( domain, dofManager );
}


template< typename FLOW_SOLVER >
void MultiphasePoromechanics< FLOW_SOLVER >::assembleSystem( real64 const GEOS_UNUSED_PARAM( time ),
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

  this->template forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
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
                                                                                                    dt,
                                                                                                    flowDofKey,
                                                                                                    flowSolver()->numFluidComponents(),
                                                                                                    flowSolver()->numFluidPhases(),
                                                                                                    flowSolver()->useTotalMassEquation(),
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
                                                                                      dt,
                                                                                      flowDofKey,
                                                                                      flowSolver()->numFluidComponents(),
                                                                                      flowSolver()->numFluidPhases(),
                                                                                      flowSolver()->useTotalMassEquation(),
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

template< typename FLOW_SOLVER >
void MultiphasePoromechanics< FLOW_SOLVER >::updateState( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  real64 maxDeltaPhaseVolFrac = 0.0;
  this->template forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                               MeshLevel & mesh,
                                                                               arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();
    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   CellElementSubRegion & subRegion )
    {
      real64 const deltaPhaseVolFrac = flowSolver()->updateFluidState( subRegion );
      maxDeltaPhaseVolFrac = LvArray::math::max( maxDeltaPhaseVolFrac, deltaPhaseVolFrac );
      if( m_isThermal )
      {
        flowSolver()->updateSolidInternalEnergyModel( subRegion );
      }
    } );
  } );

  maxDeltaPhaseVolFrac = MpiWrapper::max( maxDeltaPhaseVolFrac );

  GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "        {}: Max phase volume fraction change = {}",
                                      this->getName(), GEOS_FMT( "{:.{}f}", maxDeltaPhaseVolFrac, 4 ) ) );
}

template< typename FLOW_SOLVER >
void MultiphasePoromechanics< FLOW_SOLVER >::initializePostInitialConditionsPreSubGroups()
{
  SolverBase::initializePostInitialConditionsPreSubGroups();

  arrayView1d< string const > const & poromechanicsTargetRegionNames =
    this->template getReference< array1d< string > >( SolverBase::viewKeyStruct::targetRegionsString() );
  arrayView1d< string const > const & solidMechanicsTargetRegionNames =
    solidMechanicsSolver()->template getReference< array1d< string > >( SolverBase::viewKeyStruct::targetRegionsString() );
  arrayView1d< string const > const & flowTargetRegionNames =
    flowSolver()->template getReference< array1d< string > >( SolverBase::viewKeyStruct::targetRegionsString() );
  for( integer i = 0; i < poromechanicsTargetRegionNames.size(); ++i )
  {
    GEOS_THROW_IF( std::find( solidMechanicsTargetRegionNames.begin(), solidMechanicsTargetRegionNames.end(),
                              poromechanicsTargetRegionNames[i] )
                   == solidMechanicsTargetRegionNames.end(),
                   GEOS_FMT( "{} {}: region {} must be a target region of {}",
                             getCatalogName(), this->getDataContext(), poromechanicsTargetRegionNames[i],
                             solidMechanicsSolver()->getDataContext() ),
                   InputError );
    GEOS_THROW_IF( std::find( flowTargetRegionNames.begin(), flowTargetRegionNames.end(), poromechanicsTargetRegionNames[i] )
                   == flowTargetRegionNames.end(),
                   GEOS_FMT( "{} {}: region `{}` must be a target region of `{}`",
                             getCatalogName(), this->getDataContext(), poromechanicsTargetRegionNames[i], flowSolver()->getDataContext() ),
                   InputError );
  }

  integer & isFlowThermal = flowSolver()->isThermal();
  GEOS_WARNING_IF( m_isThermal && !isFlowThermal,
                   GEOS_FMT( "{} {}: The attribute `{}` of the flow solver `{}` is set to 1 since the poromechanics solver is thermal",
                             getCatalogName(), this->getName(), FlowSolverBase::viewKeyStruct::isThermalString(), flowSolver()->getName() ) );
  isFlowThermal = m_isThermal;

  if( m_isThermal )
  {
    this->m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::thermalMultiphasePoromechanics;
  }
}

template< typename FLOW_SOLVER >
void MultiphasePoromechanics< FLOW_SOLVER >::initializePreSubGroups()
{
  SolverBase::initializePreSubGroups();

  GEOS_THROW_IF( m_stabilizationType == StabilizationType::Local,
                 this->getWrapperDataContext( viewKeyStruct::stabilizationTypeString() ) <<
                 ": Local stabilization has been disabled temporarily",
                 InputError );

  DomainPartition & domain = this->template getGroupByPath< DomainPartition >( "/Problem/domain" );

  this->template forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                               MeshLevel & mesh,
                                                                               arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elementRegionManager = mesh.getElemManager();
    elementRegionManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                                       [&]( localIndex const,
                                                                            ElementSubRegionBase & subRegion )
    {
      string & porousName = subRegion.getReference< string >( viewKeyStruct::porousMaterialNamesString() );
      porousName = this->template getConstitutiveName< CoupledSolidBase >( subRegion );
      GEOS_ERROR_IF( porousName.empty(), GEOS_FMT( "{}: Solid model not found on subregion {}",
                                                   this->getDataContext(), subRegion.getName() ) );

      if( subRegion.hasField< fields::poromechanics::bulkDensity >() )
      {
        // get the solid model to know the number of quadrature points and resize the bulk density
        CoupledSolidBase const & solid = this->template getConstitutiveModel< CoupledSolidBase >( subRegion, porousName );
        subRegion.getField< fields::poromechanics::bulkDensity >().resizeDimension< 1 >( solid.getDensity().size( 1 ) );
      }
    } );
  } );
}

template< typename FLOW_SOLVER >
void MultiphasePoromechanics< FLOW_SOLVER >::implicitStepSetup( real64 const & time_n,
                                                                real64 const & dt,
                                                                DomainPartition & domain )
{
  flowSolver()->keepFlowVariablesConstantDuringInitStep( m_performStressInitialization );
  Base::implicitStepSetup( time_n, dt, domain );
}

template< typename FLOW_SOLVER >
void MultiphasePoromechanics< FLOW_SOLVER >::updateStabilizationParameters( DomainPartition & domain ) const
{
  // Step 1: we loop over the regions where stabilization is active and collect their name

  set< string > regionFilter;
  for( string const & regionName : m_stabilizationRegionNames )
  {
    regionFilter.insert( regionName );
  }

  // Step 2: loop over the target regions of the solver, and tag the elements belonging to stabilization regions
  this->template forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
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

      geos::constitutive::CoupledSolidBase const & porousSolid =
        this->template getConstitutiveModel< geos::constitutive::CoupledSolidBase >( subRegion, subRegion.getReference< string >( viewKeyStruct::porousMaterialNamesString() ) );

      arrayView1d< real64 const > const bulkModulus = porousSolid.getBulkModulus();
      arrayView1d< real64 const > const shearModulus = porousSolid.getShearModulus();
      arrayView1d< real64 const > const biotCoefficient = porousSolid.getBiotCoefficient();

      real64 const stabilizationMultiplier = m_stabilizationMultiplier;

      forAll< parallelDevicePolicy<> >( subRegion.size(), [bulkModulus,
                                                           shearModulus,
                                                           biotCoefficient,
                                                           stabilizationMultiplier,
                                                           macroElementIndex,
                                                           elementStabConstant] GEOS_HOST_DEVICE ( localIndex const ei )
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

template< typename FLOW_SOLVER >
void MultiphasePoromechanics< FLOW_SOLVER >::mapSolutionBetweenSolvers( DomainPartition & domain, integer const solverType )
{
  GEOS_MARK_FUNCTION;

  /// After the flow solver
  if( solverType == static_cast< integer >( Base::SolverType::Flow ) )
  {
    this->template forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
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
  if( solverType == static_cast< integer >( Base::SolverType::SolidMechanics ) )
  {
    // compute the average of the mean total stress increment over quadrature points
    averageMeanTotalStressIncrement( domain );

    this->template forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
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
  // call base method (needed to perform nonlinear acceleration)
  Base::mapSolutionBetweenSolvers( domain, solverType );
}

template< typename FLOW_SOLVER >
void MultiphasePoromechanics< FLOW_SOLVER >::updateBulkDensity( ElementSubRegionBase & subRegion )
{
  // get the fluid model (to access fluid density)
  string const fluidName = subRegion.getReference< string >( FlowSolverBase::viewKeyStruct::fluidNamesString() );
  MultiFluidBase const & fluid = this->template getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );

  // get the solid model (to access porosity and solid density)
  string const solidName = subRegion.getReference< string >( viewKeyStruct::porousMaterialNamesString() );
  CoupledSolidBase const & solid = this->template getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );

  // update the bulk density
  poromechanicsKernels::
    MultiphaseBulkDensityKernelFactory::
    createAndLaunch< parallelDevicePolicy<> >( flowSolver()->numFluidPhases(),
                                               fluid,
                                               solid,
                                               subRegion );
}

template< typename FLOW_SOLVER >
void MultiphasePoromechanics< FLOW_SOLVER >::averageMeanTotalStressIncrement( DomainPartition & domain )
{
  this->template forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                              MeshLevel & mesh,
                                                                              arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          auto & subRegion )
    {
      // get the solid model (to access stress increment)
      string const solidName = subRegion.template getReference< string >( viewKeyStruct::porousMaterialNamesString() );
      CoupledSolidBase & solid = this->template getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );

      arrayView2d< real64 const > const meanTotalStressIncrement_k = solid.getMeanTotalStressIncrement_k();
      arrayView1d< real64 > const averageMeanTotalStressIncrement_k = solid.getAverageMeanTotalStressIncrement_k();

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
                                                     meanTotalStressIncrement_k,
                                                     averageMeanTotalStressIncrement_k );
      } );
    } );
  } );
}

template class MultiphasePoromechanics< CompositionalMultiphaseBase >;
template class MultiphasePoromechanics< CompositionalMultiphaseReservoirAndWells< CompositionalMultiphaseBase > >;

namespace
{
typedef MultiphasePoromechanics< CompositionalMultiphaseReservoirAndWells< CompositionalMultiphaseBase > > MultiphaseReservoirPoromechanics;
REGISTER_CATALOG_ENTRY( SolverBase, MultiphaseReservoirPoromechanics, string const &, Group * const )
typedef MultiphasePoromechanics< CompositionalMultiphaseBase > MultiphasePoromechanics;
REGISTER_CATALOG_ENTRY( SolverBase, MultiphasePoromechanics, string const &, Group * const )
}

} /* namespace geos */
