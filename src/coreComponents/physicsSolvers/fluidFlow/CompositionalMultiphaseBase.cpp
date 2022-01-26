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
 * @file CompositionalMultiphaseBase.cpp
 */

#include "CompositionalMultiphaseBase.hpp"

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/capillaryPressure/CapillaryPressureExtrinsicData.hpp"
#include "constitutive/capillaryPressure/capillaryPressureSelector.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "constitutive/fluid/MultiFluidExtrinsicData.hpp"
#include "constitutive/fluid/multiFluidSelector.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityExtrinsicData.hpp"
#include "constitutive/relativePermeability/relativePermeabilitySelector.hpp"
#include "constitutive/thermalConductivity/thermalConductivitySelector.hpp"
#include "constitutive/permeability/PermeabilityExtrinsicData.hpp"
#include "fieldSpecification/AquiferBoundaryCondition.hpp"
#include "fieldSpecification/EquilibriumInitialCondition.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fieldSpecification/SourceFluxBoundaryCondition.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseKernels.hpp"

#if defined( __INTEL_COMPILER )
#pragma GCC optimize "O0"
#endif

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace CompositionalMultiphaseBaseKernels;

CompositionalMultiphaseBase::CompositionalMultiphaseBase( const string & name,
                                                          Group * const parent )
  :
  FlowSolverBase( name, parent ),
  m_numPhases( 0 ),
  m_numComponents( 0 ),
  m_computeCFLNumbers( 0 ),
  m_capPressureFlag( 0 ),
  m_thermalFlag( 0 ),
  m_maxCompFracChange( 1.0 ),
  m_minScalingFactor( 0.01 ),
  m_allowCompDensChopping( 1 )
{
//START_SPHINX_INCLUDE_00
  this->registerWrapper( viewKeyStruct::inputTemperatureString(), &m_inputTemperature ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Temperature" );
//END_SPHINX_INCLUDE_00
  this->registerWrapper( viewKeyStruct::useMassFlagString(), &m_useMass ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Use mass formulation instead of molar" );

  this->registerWrapper( viewKeyStruct::computeCFLNumbersString(), &m_computeCFLNumbers ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag indicating whether CFL numbers are computed or not" );

  this->registerWrapper( viewKeyStruct::relPermNamesString(), &m_relPermModelNames ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Name of the relative permeability constitutive model to use" );

  this->registerWrapper( viewKeyStruct::capPressureNamesString(), &m_capPressureModelNames ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of the capillary pressure constitutive model to use" );

  this->registerWrapper( viewKeyStruct::thermalConductivityNamesString(), &m_thermalConductivityModelNames ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::FALSE ). // input disabled temporarily
    setDescription( "Name of the thermal conductivity constitutive model to use" );

  this->registerWrapper( viewKeyStruct::maxCompFracChangeString(), &m_maxCompFracChange ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1.0 ).
    setDescription( "Maximum (absolute) change in a component fraction between two Newton iterations" );

  this->registerWrapper( viewKeyStruct::allowLocalCompDensChoppingString(), &m_allowCompDensChopping ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1 ).
    setDescription( "Flag indicating whether local (cell-wise) chopping of negative compositions is allowed" );

}

void CompositionalMultiphaseBase::postProcessInput()
{
  FlowSolverBase::postProcessInput();
  checkModelNames( m_relPermModelNames, viewKeyStruct::relPermNamesString() );
  m_capPressureFlag = checkModelNames( m_capPressureModelNames, viewKeyStruct::capPressureNamesString(), true );
  m_thermalFlag = checkModelNames( m_thermalConductivityModelNames, viewKeyStruct::thermalConductivityNamesString(), true );

  GEOSX_ERROR_IF_GT_MSG( m_maxCompFracChange, 1.0,
                         "The maximum absolute change in component fraction must smaller or equal to 1.0" );
  GEOSX_ERROR_IF_LT_MSG( m_maxCompFracChange, 0.0,
                         "The maximum absolute change in component fraction must larger or equal to 0.0" );
}

void CompositionalMultiphaseBase::registerDataOnMesh( Group & meshBodies )
{
  using namespace extrinsicMeshData::flow;

  FlowSolverBase::registerDataOnMesh( meshBodies );

  DomainPartition const & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  ConstitutiveManager const & cm = domain.getConstitutiveManager();

  // 1. Set key dimensions of the problem
  // Empty check needed to avoid accessing m_fluidModelNames when running in schema generation mode.
  if( !m_fluidModelNames.empty() )
  {
    MultiFluidBase const & fluid0 = cm.getConstitutiveRelation< MultiFluidBase >( m_fluidModelNames[0] );
    m_numPhases = fluid0.numFluidPhases();
    m_numComponents = fluid0.numFluidComponents();
  }
  m_numDofPerCell = m_numComponents + 1;

  // 2. Register and resize all fields as necessary
  meshBodies.forSubGroups< MeshBody >( [&]( MeshBody & meshBody )
  {
    MeshLevel & mesh = meshBody.getMeshLevel( 0 );

    forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion )
    {
      MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, m_fluidModelNames[targetIndex] );

      subRegion.registerExtrinsicData< pressure >( getName() );
      subRegion.registerExtrinsicData< initialPressure >( getName() );
      subRegion.registerExtrinsicData< deltaPressure >( getName() );

      subRegion.registerExtrinsicData< bcPressure >( getName() );

      subRegion.registerExtrinsicData< temperature >( getName() );

      // The resizing of the arrays needs to happen here, before the call to initializePreSubGroups,
      // to make sure that the dimensions are properly set before the timeHistoryOutput starts its initialization.

      subRegion.registerExtrinsicData< globalCompDensity >( getName() ).
        setDimLabels( 1, fluid.componentNames() ).
        reference().resizeDimension< 1 >( m_numComponents );
      subRegion.registerExtrinsicData< deltaGlobalCompDensity >( getName() ).
        reference().resizeDimension< 1 >( m_numComponents );

      subRegion.registerExtrinsicData< globalCompFraction >( getName() ).
        setDimLabels( 1, fluid.componentNames() ).
        reference().resizeDimension< 1 >( m_numComponents );
      subRegion.registerExtrinsicData< dGlobalCompFraction_dGlobalCompDensity >( getName() ).
        reference().resizeDimension< 1, 2 >( m_numComponents, m_numComponents );

      subRegion.registerExtrinsicData< phaseVolumeFraction >( getName() ).
        setDimLabels( 1, fluid.phaseNames() ).
        reference().resizeDimension< 1 >( m_numPhases );
      subRegion.registerExtrinsicData< dPhaseVolumeFraction_dPressure >( getName() ).
        reference().resizeDimension< 1 >( m_numPhases );
      subRegion.registerExtrinsicData< dPhaseVolumeFraction_dGlobalCompDensity >( getName() ).
        reference().resizeDimension< 1, 2 >( m_numPhases, m_numComponents );

      subRegion.registerExtrinsicData< phaseMobility >( getName() ).
        setDimLabels( 1, fluid.phaseNames() ).
        reference().resizeDimension< 1 >( m_numPhases );
      subRegion.registerExtrinsicData< dPhaseMobility_dPressure >( getName() ).
        reference().resizeDimension< 1 >( m_numPhases );
      subRegion.registerExtrinsicData< dPhaseMobility_dGlobalCompDensity >( getName() ).
        reference().resizeDimension< 1, 2 >( m_numPhases, m_numComponents );

      if( m_computeCFLNumbers )
      {
        subRegion.registerExtrinsicData< phaseOutflux >( getName() ).
          reference().resizeDimension< 1 >( m_numPhases );
        subRegion.registerExtrinsicData< componentOutflux >( getName() ).
          reference().resizeDimension< 1 >( m_numComponents );
        subRegion.registerExtrinsicData< phaseCFLNumber >( getName() );
        subRegion.registerExtrinsicData< componentCFLNumber >( getName() );
      }

      subRegion.registerExtrinsicData< phaseVolumeFractionOld >( getName() ).
        reference().resizeDimension< 1 >( m_numPhases );
      subRegion.registerExtrinsicData< totalDensityOld >( getName() );
      subRegion.registerExtrinsicData< phaseDensityOld >( getName() ).
        reference().resizeDimension< 1 >( m_numPhases );
      subRegion.registerExtrinsicData< phaseMobilityOld >( getName() ).
        reference().resizeDimension< 1 >( m_numPhases );
      subRegion.registerExtrinsicData< phaseComponentFractionOld >( getName() ).
        reference().resizeDimension< 1, 2 >( m_numPhases, m_numComponents );
    } );

    FaceManager & faceManager = mesh.getFaceManager();
    {
      faceManager.registerExtrinsicData< facePressure >( getName() );
    }

  } );
}

namespace
{

template< typename MODEL1_TYPE, typename MODEL2_TYPE >
void compareMultiphaseModels( MODEL1_TYPE const & lhs, MODEL2_TYPE const & rhs )
{
  GEOSX_THROW_IF_NE_MSG( lhs.numFluidPhases(), rhs.numFluidPhases(),
                         GEOSX_FMT( "Mismatch in number of phases between constitutive models {} and {}", lhs.getName(), rhs.getName() ),
                         InputError );

  for( localIndex ip = 0; ip < lhs.numFluidPhases(); ++ip )
  {
    GEOSX_THROW_IF_NE_MSG( lhs.phaseNames()[ip], rhs.phaseNames()[ip],
                           GEOSX_FMT( "Mismatch in phase names between constitutive models {} and {}", lhs.getName(), rhs.getName() ),
                           InputError );
  }
}

template< typename MODEL1_TYPE, typename MODEL2_TYPE >
void compareMulticomponentModels( MODEL1_TYPE const & lhs, MODEL2_TYPE const & rhs )
{
  GEOSX_THROW_IF_NE_MSG( lhs.numFluidComponents(), rhs.numFluidComponents(),
                         GEOSX_FMT( "Mismatch in number of components between constitutive models {} and {}", lhs.getName(), rhs.getName() ),
                         InputError );

  for( localIndex ic = 0; ic < lhs.numFluidComponents(); ++ic )
  {
    GEOSX_THROW_IF_NE_MSG( lhs.componentNames()[ic], rhs.componentNames()[ic],
                           GEOSX_FMT( "Mismatch in component names between constitutive models {} and {}", lhs.getName(), rhs.getName() ),
                           InputError );
  }
}

}

void CompositionalMultiphaseBase::validateConstitutiveModels( constitutive::ConstitutiveManager const & cm ) const
{
  MultiFluidBase const & fluid0 = cm.getConstitutiveRelation< MultiFluidBase >( m_fluidModelNames[0] );

  for( localIndex i = 1; i < m_fluidModelNames.size(); ++i )
  {
    MultiFluidBase const & fluid = cm.getConstitutiveRelation< MultiFluidBase >( m_fluidModelNames[i] );
    compareMultiphaseModels( fluid, fluid0 );
    compareMulticomponentModels( fluid, fluid0 );
  }

  RelativePermeabilityBase const & relPerm0 = cm.getConstitutiveRelation< RelativePermeabilityBase >( m_relPermModelNames[0] );
  compareMultiphaseModels( relPerm0, fluid0 );

  for( localIndex i = 1; i < m_relPermModelNames.size(); ++i )
  {
    RelativePermeabilityBase const & relPerm = cm.getConstitutiveRelation< RelativePermeabilityBase >( m_relPermModelNames[i] );
    compareMultiphaseModels( relPerm, relPerm0 );
  }

  if( m_capPressureFlag )
  {
    CapillaryPressureBase const & capPres0 = cm.getConstitutiveRelation< CapillaryPressureBase >( m_capPressureModelNames[0] );
    compareMultiphaseModels( capPres0, fluid0 );

    for( localIndex i = 1; i < m_capPressureModelNames.size(); ++i )
    {
      CapillaryPressureBase const & capPres = cm.getConstitutiveRelation< CapillaryPressureBase >( m_capPressureModelNames[i] );
      compareMultiphaseModels( capPres, capPres0 );
    }
  }

  if( m_thermalFlag )
  {
    ThermalConductivityBase const & conductivity0 = cm.getConstitutiveRelation< ThermalConductivityBase >( m_thermalConductivityModelNames[0] );
    compareMultiphaseModels( conductivity0, fluid0 );

    for( localIndex i = 1; i < m_thermalConductivityModelNames.size(); ++i )
    {
      ThermalConductivityBase const & conductivity = cm.getConstitutiveRelation< ThermalConductivityBase >( m_thermalConductivityModelNames[i] );
      compareMultiphaseModels( conductivity, conductivity0 );
    }
  }

}

void CompositionalMultiphaseBase::initializeAquiferBC( ConstitutiveManager const & cm ) const
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  MultiFluidBase const & fluid0 = cm.getConstitutiveRelation< MultiFluidBase >( m_fluidModelNames[0] );

  fsManager.forSubGroups< AquiferBoundaryCondition >( [&] ( AquiferBoundaryCondition & bc )
  {
    // set the gravity vector (needed later for the potential diff calculations)
    bc.setGravityVector( gravityVector() );

    // set the water phase index in the Aquifer boundary condition
    // note: if the water phase is not found, the fluid model is going to throw an error
    integer const waterPhaseIndex = fluid0.getWaterPhaseIndex();
    bc.setWaterPhaseIndex( waterPhaseIndex );
  } );
}

void CompositionalMultiphaseBase::validateAquiferBC( ConstitutiveManager const & cm ) const
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  MultiFluidBase const & fluid0 = cm.getConstitutiveRelation< MultiFluidBase >( m_fluidModelNames[0] );

  fsManager.forSubGroups< AquiferBoundaryCondition >( [&] ( AquiferBoundaryCondition const & bc )
  {
    arrayView1d< real64 const > const & aquiferWaterPhaseCompFrac = bc.getWaterPhaseComponentFraction();
    arrayView1d< string const > const & aquiferWaterPhaseCompNames = bc.getWaterPhaseComponentNames();

    GEOSX_ERROR_IF_NE_MSG( fluid0.numFluidComponents(), aquiferWaterPhaseCompFrac.size(),
                           "Mismatch in number of components between constitutive model "
                           << fluid0.getName() << " and the water phase composition in aquifer " << bc.getName() );

    for( localIndex ic = 0; ic < fluid0.numFluidComponents(); ++ic )
    {
      GEOSX_ERROR_IF_NE_MSG( fluid0.componentNames()[ic], aquiferWaterPhaseCompNames[ic],
                             "Mismatch in component names between constitutive model "
                             << fluid0.getName() << " and the water phase components in aquifer " << bc.getName() );
    }
  } );
}

void CompositionalMultiphaseBase::initializePreSubGroups()
{
  FlowSolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  ConstitutiveManager const & cm = domain.getConstitutiveManager();

  // 1. Validate various models against each other (must have same phases and components)
  validateConstitutiveModels( cm );

  // 2. Validate constitutive models in regions
  for( auto & mesh : domain.getMeshBodies().getSubGroups() )
  {
    MeshLevel & meshLevel = dynamicCast< MeshBody * >( mesh.second )->getMeshLevel( 0 );

    validateModelMapping< MultiFluidBase >( meshLevel.getElemManager(), m_fluidModelNames );
    validateModelMapping< RelativePermeabilityBase >( meshLevel.getElemManager(), m_relPermModelNames );
    if( m_capPressureFlag )
    {
      validateModelMapping< CapillaryPressureBase >( meshLevel.getElemManager(), m_capPressureModelNames );
    }

    // 3. Set the value of temperature
    forTargetSubRegions( meshLevel, [&]( localIndex const, ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 > const temp = subRegion.getExtrinsicData< extrinsicMeshData::flow::temperature >();
      temp.setValues< parallelHostPolicy >( m_inputTemperature );
    } );

  }

  // 3. Initialize and validate the aquifer boundary condition
  initializeAquiferBC( cm );
  validateAquiferBC( cm );
}

void CompositionalMultiphaseBase::updateComponentFraction( ObjectManagerBase & dataGroup ) const
{
  GEOSX_MARK_FUNCTION;

  ComponentFractionKernelFactory::
    createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                               dataGroup );

}

void CompositionalMultiphaseBase::updatePhaseVolumeFraction( ObjectManagerBase & dataGroup,
                                                             localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( dataGroup, m_fluidModelNames[targetIndex] );

  PhaseVolumeFractionKernelFactory::
    createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                               m_numPhases,
                                               dataGroup,
                                               fluid );

}

void CompositionalMultiphaseBase::updateFluidModel( ObjectManagerBase & dataGroup, localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  arrayView1d< real64 const > const pres = dataGroup.getExtrinsicData< extrinsicMeshData::flow::pressure >();
  arrayView1d< real64 const > const dPres = dataGroup.getExtrinsicData< extrinsicMeshData::flow::deltaPressure >();
  arrayView1d< real64 const > const temp = dataGroup.getExtrinsicData< extrinsicMeshData::flow::temperature >();
  arrayView2d< real64 const, compflow::USD_COMP > const compFrac =
    dataGroup.getExtrinsicData< extrinsicMeshData::flow::globalCompFraction >();

  MultiFluidBase & fluid = getConstitutiveModel< MultiFluidBase >( dataGroup, m_fluidModelNames[targetIndex] );

  constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
  {
    using FluidType = TYPEOFREF( castedFluid );
    using ExecPolicy = typename FluidType::exec_policy;
    typename FluidType::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

    FluidUpdateKernel::launch< ExecPolicy >( dataGroup.size(),
                                             fluidWrapper,
                                             pres,
                                             dPres,
                                             temp,
                                             compFrac );
  } );
}

void CompositionalMultiphaseBase::updateRelPermModel( ObjectManagerBase & dataGroup, localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  arrayView2d< real64 const, compflow::USD_PHASE > const phaseVolFrac =
    dataGroup.getExtrinsicData< extrinsicMeshData::flow::phaseVolumeFraction >();

  RelativePermeabilityBase & relPerm =
    getConstitutiveModel< RelativePermeabilityBase >( dataGroup, m_relPermModelNames[targetIndex] );

  constitutive::constitutiveUpdatePassThru( relPerm, [&] ( auto & castedRelPerm )
  {
    typename TYPEOFREF( castedRelPerm ) ::KernelWrapper relPermWrapper = castedRelPerm.createKernelWrapper();

    RelativePermeabilityUpdateKernel::launch< parallelDevicePolicy<> >( dataGroup.size(),
                                                                        relPermWrapper,
                                                                        phaseVolFrac );
  } );
}

void CompositionalMultiphaseBase::updateCapPressureModel( ObjectManagerBase & dataGroup, localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  if( m_capPressureFlag )
  {
    arrayView2d< real64 const, compflow::USD_PHASE > const phaseVolFrac =
      dataGroup.getExtrinsicData< extrinsicMeshData::flow::phaseVolumeFraction >();

    CapillaryPressureBase & capPressure =
      getConstitutiveModel< CapillaryPressureBase >( dataGroup, m_capPressureModelNames[targetIndex] );

    constitutive::constitutiveUpdatePassThru( capPressure, [&] ( auto & castedCapPres )
    {
      typename TYPEOFREF( castedCapPres ) ::KernelWrapper capPresWrapper = castedCapPres.createKernelWrapper();

      CapillaryPressureUpdateKernel::launch< parallelDevicePolicy<> >( dataGroup.size(),
                                                                       capPresWrapper,
                                                                       phaseVolFrac );
    } );
  }
}

void CompositionalMultiphaseBase::updateFluidState( ObjectManagerBase & subRegion, localIndex targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  updateComponentFraction( subRegion );
  updateFluidModel( subRegion, targetIndex );
  updatePhaseVolumeFraction( subRegion, targetIndex );
  updateRelPermModel( subRegion, targetIndex );
  updatePhaseMobility( subRegion, targetIndex );
  updateCapPressureModel( subRegion, targetIndex );
  // note: for now, thermal conductivity is treated explicitly, so no update here
}

void CompositionalMultiphaseBase::initializeFluidState( MeshLevel & mesh )
{
  GEOSX_MARK_FUNCTION;

  integer const numComp = m_numComponents;
  integer const numPhase = m_numPhases;

  // 1. Compute hydrostatic equilibrium in the regions for which corresponding field specification tag has been specified
  computeHydrostaticEquilibrium();

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion )
  {
    // 2. Assume global component fractions have been prescribed.
    // Initialize constitutive state to get fluid density.
    updateFluidModel( subRegion, targetIndex );

    // 3. Back-calculate global component densities from fractions and total fluid density
    // in order to initialize the primary solution variables
    MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidModelNames()[targetIndex] );
    arrayView2d< real64 const, multifluid::USD_FLUID > const totalDens = fluid.totalDensity();

    arrayView2d< real64 const, compflow::USD_COMP > const compFrac =
      subRegion.getExtrinsicData< extrinsicMeshData::flow::globalCompFraction >();
    arrayView2d< real64, compflow::USD_COMP > const compDens =
      subRegion.getExtrinsicData< extrinsicMeshData::flow::globalCompDensity >();

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      for( localIndex ic = 0; ic < numComp; ++ic )
      {
        compDens[ei][ic] = totalDens[ei][0] * compFrac[ei][ic];
      }
    } );

  } );

  // for some reason CUDA does not want the host_device lambda to be defined inside the generic lambda
  // I need the exact type of the subRegion for updateSolidflowProperties to work well.
  forTargetSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( mesh, [&]( localIndex const targetIndex,
                                                                                   auto & subRegion )
  {
    // 4. Update dependent state quantities
    // Note that the order used below is important, as some constitutive models are initialized with values from other constitutive models

    // 4.1 First, we update the porosity and permeability, and save the porosity into the "old porosity"
    updatePorosityAndPermeability( subRegion, targetIndex );
    CoupledSolidBase const & porousMaterial = getConstitutiveModel< CoupledSolidBase >( subRegion, m_solidModelNames[targetIndex] );
    porousMaterial.initializeState();

    // 4.2 Then, we initialize the capillary pressure model (which can depend on porosity and permeability)
    // note: this **must** be called after the porosity update, and **before** calling updateCapPressureModel
    if( m_capPressureFlag )
    {
      // initialized porosity
      arrayView2d< real64 const > const porosity = porousMaterial.getPorosity();

      PermeabilityBase const & permeabilityMaterial =
        getConstitutiveModel< PermeabilityBase >( subRegion, m_permeabilityModelNames[targetIndex] );
      // initialized permeability
      arrayView3d< real64 const > const permeability = permeabilityMaterial.permeability();

      CapillaryPressureBase const & capPressureMaterial =
        getConstitutiveModel< CapillaryPressureBase >( subRegion, m_capPressureModelNames[targetIndex] );
      capPressureMaterial.initializeRockState( porosity, permeability );
    }

    // 4.3 Then, we call the remaining constitutive models to perform the updates
    updatePhaseVolumeFraction( subRegion, targetIndex );
    updateRelPermModel( subRegion, targetIndex );
    updatePhaseMobility( subRegion, targetIndex );
    updateCapPressureModel( subRegion, targetIndex );
    // thermal conductivity is explicitly, so no update here

    // 4.4 Finally, we initialize the thermal conductivity (which can depend on porosity and phase volume fraction)
    if( m_thermalFlag )
    {
      // initialized porosity
      arrayView2d< real64 const > const porosity = porousMaterial.getPorosity();

      // initialized phase volume fraction
      arrayView2d< real64 const, compflow::USD_PHASE > const phaseVolFrac =
        subRegion.template getExtrinsicData< extrinsicMeshData::flow::phaseVolumeFraction >();

      ThermalConductivityBase const & conductivityMaterial =
        getConstitutiveModel< ThermalConductivityBase >( subRegion, m_thermalConductivityModelNames[targetIndex] );
      conductivityMaterial.initializeRockFluidState( porosity, phaseVolFrac );
    }

  } );

  // 5. Save initial pressure and total mass density (needed by the poromechanics solvers)
  //    Specifically, the initial pressure is used to compute a deltaPressure = currentPres - initPres in the total stress
  //    And the initial total mass density is used to compute a deltaBodyForce
  forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion )
  {
    arrayView1d< real64 const > const pres = subRegion.getExtrinsicData< extrinsicMeshData::flow::pressure >();
    arrayView1d< real64 > const initPres = subRegion.getExtrinsicData< extrinsicMeshData::flow::initialPressure >();
    arrayView2d< real64 const, compflow::USD_PHASE > const phaseVolFrac =
      subRegion.getExtrinsicData< extrinsicMeshData::flow::phaseVolumeFraction >();


    MultiFluidBase & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidModelNames()[targetIndex] );
    arrayView3d< real64 const, multifluid::USD_PHASE > const phaseMassDens = fluid.phaseMassDensity();
    arrayView2d< real64, multifluid::USD_FLUID > const initTotalMassDens =
      fluid.getReference< extrinsicMeshData::multifluid::initialTotalMassDensity::type >( extrinsicMeshData::multifluid::initialTotalMassDensity::key() );

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      initPres[ei] = pres[ei];
      initTotalMassDens[ei][0] = 0.0;
      for( localIndex ip = 0; ip < numPhase; ++ip )
      {
        initTotalMassDens[ei][0] += phaseVolFrac[ei][ip] * phaseMassDens[ei][0][ip];
      }
    } );
  } );
}

void CompositionalMultiphaseBase::computeHydrostaticEquilibrium()
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  localIndex const numComps = m_numComponents;
  localIndex const numPhases = m_numPhases;

  real64 const gravVector[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );

  // Step 1: count individual equilibriums (there may be multiple ones)

  std::map< string, localIndex > equilNameToEquilId;
  localIndex equilCounter = 0;

  fsManager.forSubGroups< EquilibriumInitialCondition >( [&] ( EquilibriumInitialCondition const & bc )
  {

    // collect all the equilibrium names to idx
    equilNameToEquilId[bc.getName()] = equilCounter;
    equilCounter++;

    // check that the gravity vector is aligned with the z-axis
    GEOSX_THROW_IF( !isZero( gravVector[0] ) || !isZero( gravVector[1] ),
                    catalogName() << " " << getName() <<
                    ": the gravity vector specified in this simulation (" << gravVector[0] << " " << gravVector[1] << " " << gravVector[2] <<
                    ") is not aligned with the z-axis. \n"
                    "This is incompatible with the " << EquilibriumInitialCondition::catalogName() << " called " << bc.getName() <<
                    "used in this simulation. To proceed, you can either: \n" <<
                    "   - Use a gravityVector aligned with the z-axis, such as (0.0,0.0,-9.81)\n" <<
                    "   - Remove the hydrostatic equilibrium initial condition from the XML file",
                    InputError );

  } );

  if( equilCounter == 0 )
  {
    return;
  }

  // Step 2: find the min elevation and the max elevation in the targetSets

  array1d< real64 > globalMaxElevation( equilNameToEquilId.size() );
  array1d< real64 > globalMinElevation( equilNameToEquilId.size() );
  findMinMaxElevationInEquilibriumTarget( domain,
                                          equilNameToEquilId,
                                          globalMaxElevation,
                                          globalMinElevation );

  // Step 3: for each equil, compute a fine table with hydrostatic pressure vs elevation if the region is a target region

  // first compute the region filter
  std::set< string > regionFilter;
  for( string const & regionName : targetRegionNames() )
  {
    regionFilter.insert( regionName );
  }

  fsManager.apply< EquilibriumInitialCondition >( 0.0,
                                                  domain,
                                                  "ElementRegions",
                                                  EquilibriumInitialCondition::catalogName(),
                                                  [&] ( EquilibriumInitialCondition const & fs,
                                                        string const &,
                                                        SortedArrayView< localIndex const > const & targetSet,
                                                        Group & subRegion,
                                                        string const & )
  {
    // Step 3.1: retrieve the data necessary to construct the pressure table in this subregion

    integer const maxNumEquilIterations = fs.getMaxNumEquilibrationIterations();
    real64 const equilTolerance = fs.getEquilibrationTolerance();
    real64 const datumElevation = fs.getDatumElevation();
    real64 const datumPressure = fs.getDatumPressure();
    string const initPhaseName = fs.getInitPhaseName(); // will go away when GOC/WOC are implemented

    localIndex const equilIndex = equilNameToEquilId.at( fs.getName() );
    real64 const minElevation = LvArray::math::min( globalMinElevation[equilIndex], datumElevation );
    real64 const maxElevation = LvArray::math::max( globalMaxElevation[equilIndex], datumElevation );
    real64 const elevationIncrement = LvArray::math::min( fs.getElevationIncrement(), maxElevation - minElevation );
    localIndex const numPointsInTable = std::ceil( (maxElevation - minElevation) / elevationIncrement ) + 1;

    real64 const eps = 0.1 * (maxElevation - minElevation); // we add a small buffer to only log in the pathological cases
    GEOSX_LOG_RANK_0_IF( ( (datumElevation > globalMaxElevation[equilIndex]+eps)  || (datumElevation < globalMinElevation[equilIndex]-eps) ),
                         CompositionalMultiphaseBase::catalogName() << " " << getName()
                                                                    << ": By looking at the elevation of the cell centers in this model, GEOSX found that "
                                                                    << "the min elevation is " << globalMinElevation[equilIndex] << " and the max elevation is " << globalMaxElevation[equilIndex] <<
                         "\n"
                                                                    << "But, a datum elevation of " << datumElevation << " was specified in the input file to equilibrate the model.\n "
                                                                    << "The simulation is going to proceed with this out-of-bound datum elevation, but the initial condition may be inaccurate." );

    array1d< array1d< real64 > > elevationValues;
    array1d< real64 > pressureValues;
    elevationValues.resize( 1 );
    elevationValues[0].resize( numPointsInTable );
    pressureValues.resize( numPointsInTable );

    // Step 3.2: retrieve the user-defined tables (temperature and comp fraction)

    FunctionManager & functionManager = FunctionManager::getInstance();

    array1d< TableFunction::KernelWrapper > compFracTableWrappers;
    arrayView1d< string const > compFracTableNames = fs.getComponentFractionVsElevationTableNames();
    for( localIndex ic = 0; ic < numComps; ++ic )
    {
      TableFunction const & compFracTable = functionManager.getGroup< TableFunction >( compFracTableNames[ic] );
      compFracTableWrappers.emplace_back( compFracTable.createKernelWrapper() );
    }

    string const tempTableName = fs.getTemperatureVsElevationTableName();
    TableFunction const & tempTable = functionManager.getGroup< TableFunction >( tempTableName );
    TableFunction::KernelWrapper tempTableWrapper = tempTable.createKernelWrapper();

    // Step 3.3: retrieve the fluid model to compute densities
    // we end up with the same issue as in applyDirichletBC: there is not a clean way to retrieve the fluid info

    Group const & region = subRegion.getParent().getParent();
    auto itRegionFilter = regionFilter.find( region.getName() );
    if( itRegionFilter == regionFilter.end() )
    {
      return; // the region is not in target, there is nothing to do
    }
    string const & fluidName = m_fluidModelNames[ targetRegionIndex( region.getName() ) ];
    MultiFluidBase & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );

    arrayView1d< string const > componentNames = fs.getComponentNames();
    GEOSX_THROW_IF( fluid.componentNames().size() != componentNames.size(),
                    "Mismatch in number of components between constitutive model "
                    << fluid.getName() << " and the Equilibrium initial condition " << fs.getName(),
                    InputError );
    for( localIndex ic = 0; ic < fluid.numFluidComponents(); ++ic )
    {
      GEOSX_THROW_IF( fluid.componentNames()[ic] != componentNames[ic],
                      "Mismatch in component names between constitutive model "
                      << fluid.getName() << " and the Equilibrium initial condition " << fs.getName(),
                      InputError );
    }

    // Note: for now, we assume that the reservoir is in a single-phase state at initialization
    arrayView1d< string const > phaseNames = fluid.phaseNames();
    auto const itPhaseNames = std::find( std::begin( phaseNames ), std::end( phaseNames ), initPhaseName );
    GEOSX_THROW_IF( itPhaseNames == std::end( phaseNames ),
                    CompositionalMultiphaseBase::catalogName() << " " << getName() << ": phase name " << initPhaseName
                                                               << " not found in the phases of " << fluid.getName(),
                    InputError );
    integer const ipInit = std::distance( std::begin( phaseNames ), itPhaseNames );

    // Step 3.4: compute the hydrostatic pressure values

    constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
    {
      using FluidType = TYPEOFREF( castedFluid );
      typename FluidType::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

      // note: inside this kernel, serialPolicy is used, and elevation/pressure values don't go to the GPU
      HydrostaticPressureKernel::ReturnType const returnValue =
        HydrostaticPressureKernel::launch( numPointsInTable,
                                           numComps,
                                           numPhases,
                                           ipInit,
                                           maxNumEquilIterations,
                                           equilTolerance,
                                           gravVector,
                                           minElevation,
                                           elevationIncrement,
                                           datumElevation,
                                           datumPressure,
                                           fluidWrapper,
                                           compFracTableWrappers.toViewConst(),
                                           tempTableWrapper,
                                           elevationValues.toNestedView(),
                                           pressureValues.toView() );

      GEOSX_THROW_IF( returnValue == HydrostaticPressureKernel::ReturnType::FAILED_TO_CONVERGE,
                      CompositionalMultiphaseBase::catalogName() << " " << getName()
                                                                 << ": hydrostatic pressure initialization failed to converge in region " << region.getName() << "! \n"
                                                                 << "Try to loosen the equilibration tolerance, or increase the number of equilibration iterations. \n"
                                                                 << "If nothing works, something may be wrong in the fluid model, see <Constitutive> ",
                      std::runtime_error );

      GEOSX_LOG_RANK_0_IF( returnValue == HydrostaticPressureKernel::ReturnType::DETECTED_MULTIPHASE_FLOW,
                           CompositionalMultiphaseBase::catalogName() << " " << getName()
                                                                      << ": currently, GEOSX assumes that there is only one mobile phase when computing the hydrostatic pressure. \n"
                                                                      << "We detected multiple phases using the provided datum pressure, temperature, and component fractions. \n"
                                                                      << "Please make sure that only one phase is mobile at the beginning of the simulation. \n"
                                                                      << "If this is not the case, the problem will not be at equilibrium when the simulation starts" );

    } );

    // Step 3.5: create hydrostatic pressure table

    string const tableName = fs.getName() + "_" + subRegion.getName() + "_" + phaseNames[ipInit] + "_table";
    TableFunction * const presTable = dynamicCast< TableFunction * >( functionManager.createChild( TableFunction::catalogName(), tableName ) );
    presTable->setTableCoordinates( elevationValues );
    presTable->setTableValues( pressureValues );
    presTable->setInterpolationMethod( TableFunction::InterpolationType::Linear );
    TableFunction::KernelWrapper presTableWrapper = presTable->createKernelWrapper();

    // Step 4: assign pressure, temperature, and component fraction as a function of elevation
    // TODO: this last step should probably be delayed to wait for the creation of FaceElements
    // TODO: this last step should be modified to account for GOC and WOC
    arrayView2d< real64 const > const elemCenter =
      subRegion.getReference< array2d< real64 > >( ElementSubRegionBase::viewKeyStruct::elementCenterString() );

    arrayView1d< real64 > const pres = subRegion.getReference< array1d< real64 > >( extrinsicMeshData::flow::pressure::key() );
    arrayView1d< real64 > const temp = subRegion.getReference< array1d< real64 > >( extrinsicMeshData::flow::temperature::key() );
    arrayView2d< real64, compflow::USD_COMP > const compFrac =
      subRegion.getReference< array2d< real64, compflow::LAYOUT_COMP > >( extrinsicMeshData::flow::globalCompFraction::key() );
    arrayView1d< TableFunction::KernelWrapper const > compFracTableWrappersViewConst =
      compFracTableWrappers.toViewConst();

    forAll< parallelDevicePolicy<> >( targetSet.size(), [targetSet,
                                                         elemCenter,
                                                         presTableWrapper,
                                                         tempTableWrapper,
                                                         compFracTableWrappersViewConst,
                                                         numComps,
                                                         pres,
                                                         temp,
                                                         compFrac] GEOSX_HOST_DEVICE ( localIndex const i )
    {
      localIndex const k = targetSet[i];
      real64 const elevation = elemCenter[k][2];

      pres[k] = presTableWrapper.compute( &elevation );
      temp[k] = tempTableWrapper.compute( &elevation );
      for( localIndex ic = 0; ic < numComps; ++ic )
      {
        compFrac[k][ic] = compFracTableWrappersViewConst[ic].compute( &elevation );
      }
    } );
  } );
}

void CompositionalMultiphaseBase::initializePostInitialConditionsPreSubGroups()
{
  GEOSX_MARK_FUNCTION;

  FlowSolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( extrinsicMeshData::flow::pressure::key() );
  fieldNames["elems"].emplace_back( extrinsicMeshData::flow::globalCompDensity::key() );

  CommunicationTools::getInstance().synchronizeFields( fieldNames, mesh, domain.getNeighbors(), false );

  // set mass fraction flag on fluid models
  forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion )
  {
    MultiFluidBase & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, m_fluidModelNames[targetIndex] );
    fluid.setMassFlag( m_useMass );
  } );

  // Initialize primary variables from applied initial conditions
  initializeFluidState( mesh );
}

real64 CompositionalMultiphaseBase::solverStep( real64 const & time_n,
                                                real64 const & dt,
                                                integer const cycleNumber,
                                                DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  // Only build the sparsity pattern once
  // TODO: this should be triggered by a topology change indicator
  static bool systemSetupDone = false;
  if( !systemSetupDone )
  {
    setupSystem( domain, m_dofManager, m_localMatrix, m_rhs, m_solution );
    systemSetupDone = true;
  }

  implicitStepSetup( time_n, dt, domain );

  // currently the only method is implicit time integration
  real64 const dt_return = nonlinearImplicitStep( time_n, dt, cycleNumber, domain );

  // final step for completion of timestep. typically secondary variable updates and cleanup.
  implicitStepComplete( time_n, dt_return, domain );

  return dt_return;
}

void CompositionalMultiphaseBase::backupFields( MeshLevel & mesh ) const
{
  GEOSX_MARK_FUNCTION;

  integer const numComp = m_numComponents;
  integer const numPhase = m_numPhases;

  // backup some fields used in time derivative approximation
  forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion )
  {
    arrayView1d< integer const > const elemGhostRank = subRegion.ghostRank();

    arrayView2d< real64 const, compflow::USD_PHASE > const phaseVolFrac =
      subRegion.getExtrinsicData< extrinsicMeshData::flow::phaseVolumeFraction >();
    arrayView2d< real64 const, compflow::USD_PHASE > const & phaseMob =
      subRegion.getExtrinsicData< extrinsicMeshData::flow::phaseMobility >();

    MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidModelNames()[targetIndex] );
    arrayView2d< real64 const, multifluid::USD_FLUID > const totalDens = fluid.totalDensity();
    arrayView3d< real64 const, multifluid::USD_PHASE > const phaseDens = fluid.phaseDensity();
    arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const phaseCompFrac = fluid.phaseCompFraction();

    arrayView1d< real64 > const totalDensOld =
      subRegion.getExtrinsicData< extrinsicMeshData::flow::totalDensityOld >();

    arrayView2d< real64, compflow::USD_PHASE > const phaseDensOld =
      subRegion.getExtrinsicData< extrinsicMeshData::flow::phaseDensityOld >();
    arrayView2d< real64, compflow::USD_PHASE > const phaseVolFracOld =
      subRegion.getExtrinsicData< extrinsicMeshData::flow::phaseVolumeFractionOld >();
    arrayView2d< real64, compflow::USD_PHASE > const phaseMobOld =
      subRegion.getExtrinsicData< extrinsicMeshData::flow::phaseMobilityOld >();
    arrayView3d< real64, compflow::USD_PHASE_COMP > const phaseCompFracOld =
      subRegion.getExtrinsicData< extrinsicMeshData::flow::phaseComponentFractionOld >();

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      if( elemGhostRank[ei] >= 0 )
        return;

      for( localIndex ip = 0; ip < numPhase; ++ip )
      {
        phaseDensOld[ei][ip] = phaseDens[ei][0][ip];
        phaseVolFracOld[ei][ip] = phaseVolFrac[ei][ip];
        phaseMobOld[ei][ip] = phaseMob[ei][ip];

        for( localIndex ic = 0; ic < numComp; ++ic )
        {
          phaseCompFracOld[ei][ip][ic] = phaseCompFrac[ei][0][ip][ic];
        }
      }
      totalDensOld[ei] = totalDens[ei][0];
    } );
  } );
}

void
CompositionalMultiphaseBase::implicitStepSetup( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                                real64 const & GEOSX_UNUSED_PARAM( dt ),
                                                DomainPartition & domain )
{
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  // set deltas to zero and recompute dependent quantities
  resetStateToBeginningOfStep( domain );

  // backup fields used in time derivative approximation
  backupFields( mesh );
}

void CompositionalMultiphaseBase::assembleSystem( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                                  real64 const dt,
                                                  DomainPartition & domain,
                                                  DofManager const & dofManager,
                                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                  arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  assembleAccumulationAndVolumeBalanceTerms( domain,
                                             dofManager,
                                             localMatrix,
                                             localRhs );

  assembleFluxTerms( dt,
                     domain,
                     dofManager,
                     localMatrix,
                     localRhs );


}

void CompositionalMultiphaseBase::assembleAccumulationAndVolumeBalanceTerms( DomainPartition & domain,
                                                                             DofManager const & dofManager,
                                                                             CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                             arrayView1d< real64 > const & localRhs ) const
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  forTargetSubRegions( mesh,
                       [&]( localIndex const targetIndex,
                            ElementSubRegionBase const & subRegion )
  {
    string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
    MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidModelNames()[targetIndex] );
    CoupledSolidBase const & solid = getConstitutiveModel< CoupledSolidBase >( subRegion, m_solidModelNames[targetIndex] );

    ElementBasedAssemblyKernelFactory::
      createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                 m_numPhases,
                                                 dofManager.rankOffset(),
                                                 dofKey,
                                                 subRegion,
                                                 fluid,
                                                 solid,
                                                 localMatrix,
                                                 localRhs );

  } );
}

void CompositionalMultiphaseBase::applyBoundaryConditions( real64 const time_n,
                                                           real64 const dt,
                                                           DomainPartition & domain,
                                                           DofManager const & dofManager,
                                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                           arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  // apply pressure boundary conditions.
  applyDirichletBC( time_n, dt, dofManager, domain, localMatrix.toViewConstSizes(), localRhs.toView() );

  // apply flux boundary conditions
  applySourceFluxBC( time_n, dt, dofManager, domain, localMatrix.toViewConstSizes(), localRhs.toView() );

  // apply aquifer boundary conditions
  applyAquiferBC( time_n, dt, dofManager, domain, localMatrix.toViewConstSizes(), localRhs.toView() );
}

namespace internal
{
string const bcLogMessage = string( "CompositionalMultiphaseBase {}: at time {}s, " )
                            + string( "the <{}> boundary condition '{}' is applied to the element set '{}' in subRegion '{}'. " )
                            + string( "\nThe scale of this boundary condition is {} and multiplies the value of the provided function (if any). " )
                            + string( "\nThe total number of target elements (including ghost elements) is {}. " )
                            + string( "\nNote that if this number is equal to zero for all subRegions, the boundary condition will not be applied on this element set." );
}

void CompositionalMultiphaseBase::applySourceFluxBC( real64 const time,
                                                     real64 const dt,
                                                     DofManager const & dofManager,
                                                     DomainPartition & domain,
                                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                     arrayView1d< real64 > const & localRhs ) const
{
  GEOSX_MARK_FUNCTION;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  fsManager.apply( time + dt,
                   domain,
                   "ElementRegions",
                   FieldSpecificationBase::viewKeyStruct::fluxBoundaryConditionString(),
                   [&]( FieldSpecificationBase const & fs,
                        string const & setName,
                        SortedArrayView< localIndex const > const & targetSet,
                        Group & subRegion,
                        string const & )
  {
    if( fs.getLogLevel() >= 1 && m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
    {
      globalIndex const numTargetElems = MpiWrapper::sum< globalIndex >( targetSet.size() );
      GEOSX_LOG_RANK_0( GEOSX_FMT( geosx::internal::bcLogMessage,
                                   getName(), time+dt, SourceFluxBoundaryCondition::catalogName(),
                                   fs.getName(), setName, subRegion.getName(), fs.getScale(), numTargetElems ) );
    }

    arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const ghostRank =
      subRegion.getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString() );

    // Step 1: get the values of the source boundary condition that need to be added to the rhs
    // We don't use FieldSpecificationBase::applyConditionToSystem here because we want to account for the row permutation used in the
    // compositional solvers

    array1d< globalIndex > dofArray( targetSet.size() );
    array1d< real64 > rhsContributionArray( targetSet.size() );
    arrayView1d< real64 > rhsContributionArrayView = rhsContributionArray.toView();
    localIndex const rankOffset = dofManager.rankOffset();

    // note that the dofArray will not be used after this step (simpler to use dofNumber instead)
    fs.computeRhsContribution< FieldSpecificationAdd,
                               parallelDevicePolicy<> >( targetSet.toViewConst(),
                                                         time + dt,
                                                         dt,
                                                         subRegion,
                                                         dofNumber,
                                                         rankOffset,
                                                         localMatrix,
                                                         dofArray.toView(),
                                                         rhsContributionArrayView,
                                                         [] GEOSX_HOST_DEVICE ( localIndex const )
    {
      return 0.0;
    } );

    // Step 2: we are ready to add the right-hand side contributions, taking into account our equation layout

    integer const fluidComponentId = fs.getComponent();
    integer const numFluidComponents = m_numComponents;
    forAll< parallelDevicePolicy<> >( targetSet.size(), [targetSet,
                                                         rankOffset,
                                                         ghostRank,
                                                         fluidComponentId,
                                                         numFluidComponents,
                                                         dofNumber,
                                                         rhsContributionArrayView,
                                                         localRhs] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      // we need to filter out ghosts here, because targetSet may contain them
      localIndex const ei = targetSet[a];
      if( ghostRank[ei] >= 0 )
      {
        return;
      }

      // for all "fluid components", we add the value to the total mass balance equation
      globalIndex const totalMassBalanceRow = dofNumber[ei] - rankOffset;
      localRhs[totalMassBalanceRow] += rhsContributionArrayView[a];

      // for all "fluid components" except the last one, we add the value to the component mass balance equation (shifted appropriately)
      if( fluidComponentId < numFluidComponents - 1 )
      {
        globalIndex const compMassBalanceRow = totalMassBalanceRow + fluidComponentId + 1; // component mass bal equations are shifted
        localRhs[compMassBalanceRow] += rhsContributionArrayView[a];
      }
    } );

  } );
}

namespace
{

bool validateDirichletBC( DomainPartition & domain,
                          integer const numComp,
                          real64 const time )
{
  constexpr integer MAX_NC = MultiFluidBase::MAX_NUM_COMPONENTS;
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  map< string, map< string, map< string, ComponentMask< MAX_NC > > > > bcStatusMap; // map to check consistent application of BC
  bool bcConsistent = true;

  // 1. Check pressure Dirichlet BCs
  fsManager.apply( time,
                   domain,
                   "ElementRegions",
                   extrinsicMeshData::flow::pressure::key(),
                   [&]( FieldSpecificationBase const &,
                        string const & setName,
                        SortedArrayView< localIndex const > const &,
                        Group & subRegion,
                        string const & )
  {
    // 1.0. Check whether pressure has already been applied to this set
    string const & subRegionName = subRegion.getName();
    string const & regionName = subRegion.getParent().getParent().getName();

    auto & subRegionSetMap = bcStatusMap[regionName][subRegionName];
    if( subRegionSetMap.count( setName ) > 0 )
    {
      bcConsistent = false;
      GEOSX_WARNING( GEOSX_FMT( "Conflicting pressure boundary conditions on set {}/{}/{}", regionName, subRegionName, setName ) );
    }
    subRegionSetMap[setName].setNumComp( numComp );
  } );

  // 2. Check composition BC (global component fraction)
  fsManager.apply( time,
                   domain,
                   "ElementRegions",
                   extrinsicMeshData::flow::globalCompFraction::key(),
                   [&] ( FieldSpecificationBase const & fs,
                         string const & setName,
                         SortedArrayView< localIndex const > const &,
                         Group & subRegion,
                         string const & )
  {
    // 2.0. Check pressure and record composition bc application
    string const & subRegionName = subRegion.getName();
    string const & regionName = subRegion.getParent().getParent().getName();
    integer const comp = fs.getComponent();

    auto & subRegionSetMap = bcStatusMap[regionName][subRegionName];
    if( subRegionSetMap.count( setName ) == 0 )
    {
      bcConsistent = false;
      GEOSX_WARNING( GEOSX_FMT( "Pressure boundary condition not prescribed on set {}/{}/{}", regionName, subRegionName, setName ) );
    }
    if( comp < 0 || comp >= numComp )
    {
      bcConsistent = false;
      GEOSX_WARNING( GEOSX_FMT( "Invalid component index [{}] in composition boundary condition {}", comp, fs.getName() ) );
      return; // can't check next part with invalid component id
    }

    ComponentMask< MAX_NC > & compMask = subRegionSetMap[setName];
    if( compMask[comp] )
    {
      bcConsistent = false;
      GEOSX_WARNING( GEOSX_FMT( "Conflicting composition[{}] boundary conditions on set {}/{}/{}", comp, regionName, subRegionName, setName ) );
    }
    compMask.set( comp );
  } );

  // 2.3 Check consistency between composition BC applied to sets
  for( auto const & regionEntry : bcStatusMap )
  {
    for( auto const & subRegionEntry : regionEntry.second )
    {
      for( auto const & setEntry : subRegionEntry.second )
      {
        ComponentMask< MAX_NC > const & compMask = setEntry.second;
        for( integer ic = 0; ic < numComp; ++ic )
        {
          if( !compMask[ic] )
          {
            bcConsistent = false;
            GEOSX_WARNING( GEOSX_FMT( "Boundary condition not applied to composition[{}] on set {}/{}/{}",
                                      ic, regionEntry.first, subRegionEntry.first, setEntry.first ) );
          }
        }
      }
    }
  }

  return bcConsistent;
}

}


void CompositionalMultiphaseBase::applyDirichletBC( real64 const time,
                                                    real64 const dt,
                                                    DofManager const & dofManager,
                                                    DomainPartition & domain,
                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                    arrayView1d< real64 > const & localRhs ) const
{
  GEOSX_MARK_FUNCTION;

  // Only validate BC at the beginning of Newton loop
  if( m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
  {
    bool const bcConsistent = validateDirichletBC( domain, m_numComponents, time + dt );
    GEOSX_ERROR_IF( !bcConsistent, GEOSX_FMT( "CompositionalMultiphaseBase {}: inconsistent boundary conditions", getName() ) );
  }

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  // 1. Apply pressure Dirichlet BCs, store in a separate field
  fsManager.apply( time + dt,
                   domain,
                   "ElementRegions",
                   extrinsicMeshData::flow::pressure::key(),
                   [&]( FieldSpecificationBase const & fs,
                        string const & setName,
                        SortedArrayView< localIndex const > const & targetSet,
                        Group & subRegion,
                        string const & )
  {
    if( fs.getLogLevel() >= 1 && m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
    {
      globalIndex const numTargetElems = MpiWrapper::sum< globalIndex >( targetSet.size() );
      GEOSX_LOG_RANK_0( GEOSX_FMT( geosx::internal::bcLogMessage,
                                   getName(), time+dt, FieldSpecificationBase::catalogName(),
                                   fs.getName(), setName, subRegion.getName(), fs.getScale(), numTargetElems ) );
    }

    fs.applyFieldValue< FieldSpecificationEqual, parallelDevicePolicy<> >( targetSet,
                                                                           time + dt,
                                                                           subRegion,
                                                                           extrinsicMeshData::flow::bcPressure::key() );
  } );

  // 2. Apply composition BC (global component fraction) and store them for constitutive call
  fsManager.apply( time + dt,
                   domain,
                   "ElementRegions",
                   extrinsicMeshData::flow::globalCompFraction::key(),
                   [&] ( FieldSpecificationBase const & fs,
                         string const &,
                         SortedArrayView< localIndex const > const & targetSet,
                         Group & subRegion,
                         string const & )
  {
    fs.applyFieldValue< FieldSpecificationEqual, parallelDevicePolicy<> >( targetSet,
                                                                           time + dt,
                                                                           subRegion,
                                                                           extrinsicMeshData::flow::globalCompFraction::key() );
  } );

  globalIndex const rankOffset = dofManager.rankOffset();
  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  // 3. Call constitutive update, back-calculate target global component densities and apply to the system
  fsManager.apply( time + dt,
                   domain,
                   "ElementRegions",
                   extrinsicMeshData::flow::pressure::key(),
                   [&] ( FieldSpecificationBase const &,
                         string const &,
                         SortedArrayView< localIndex const > const & targetSet,
                         Group & subRegion,
                         string const & )
  {
    // TODO: hack! Find a better way to get the fluid
    Group const & region = subRegion.getParent().getParent();
    string const & fluidName = m_fluidModelNames[ targetRegionIndex( region.getName() ) ];
    MultiFluidBase & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );

    arrayView1d< real64 const > const bcPres =
      subRegion.getReference< array1d< real64 > >( extrinsicMeshData::flow::bcPressure::key() );
    arrayView1d< real64 const > const temp =
      subRegion.getReference< array1d< real64 > >( extrinsicMeshData::flow::temperature::key() );
    arrayView2d< real64 const, compflow::USD_COMP > const compFrac =
      subRegion.getReference< array2d< real64, compflow::LAYOUT_COMP > >( extrinsicMeshData::flow::globalCompFraction::key() );

    constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
    {
      using FluidType = TYPEOFREF( castedFluid );
      using ExecPolicy = typename FluidType::exec_policy;
      typename FluidType::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

      FluidUpdateKernel::launch< ExecPolicy >( targetSet,
                                               fluidWrapper,
                                               bcPres,
                                               temp,
                                               compFrac );
    } );

    arrayView1d< integer const > const ghostRank =
      subRegion.getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString() );
    arrayView1d< globalIndex const > const dofNumber =
      subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< real64 const > const pres =
      subRegion.getReference< array1d< real64 > >( extrinsicMeshData::flow::pressure::key() );
    arrayView1d< real64 const > const dPres =
      subRegion.getReference< array1d< real64 > >( extrinsicMeshData::flow::deltaPressure::key() );
    arrayView2d< real64 const, compflow::USD_COMP > const compDens =
      subRegion.getReference< array2d< real64, compflow::LAYOUT_COMP > >( extrinsicMeshData::flow::globalCompDensity::key() );
    arrayView2d< real64 const, compflow::USD_COMP > const dCompDens =
      subRegion.getReference< array2d< real64, compflow::LAYOUT_COMP > >( extrinsicMeshData::flow::deltaGlobalCompDensity::key() );
    arrayView2d< real64 const, multifluid::USD_FLUID > const totalDens = fluid.totalDensity();

    integer const numComp = m_numComponents;
    forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      localIndex const ei = targetSet[a];
      if( ghostRank[ei] >= 0 )
      {
        return;
      }

      globalIndex const dofIndex = dofNumber[ei];
      localIndex const localRow = dofIndex - rankOffset;
      real64 rhsValue;

      // 3.1. Apply pressure value to the matrix/rhs
      FieldSpecificationEqual::SpecifyFieldValue( dofIndex,
                                                  rankOffset,
                                                  localMatrix,
                                                  rhsValue,
                                                  bcPres[ei],
                                                  pres[ei] + dPres[ei] );
      localRhs[localRow] = rhsValue;

      // 3.2. For each component, apply target global density value
      for( localIndex ic = 0; ic < numComp; ++ic )
      {
        FieldSpecificationEqual::SpecifyFieldValue( dofIndex + ic + 1,
                                                    rankOffset,
                                                    localMatrix,
                                                    rhsValue,
                                                    totalDens[ei][0] * compFrac[ei][ic],
                                                    compDens[ei][ic] + dCompDens[ei][ic] );
        localRhs[localRow + ic + 1] = rhsValue;
      }
    } );
  } );

}

void CompositionalMultiphaseBase::solveSystem( DofManager const & dofManager,
                                               ParallelMatrix & matrix,
                                               ParallelVector & rhs,
                                               ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  rhs.scale( -1.0 );
  solution.zero();

  SolverBase::solveSystem( dofManager, matrix, rhs, solution );
}

void CompositionalMultiphaseBase::chopNegativeDensities( DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  integer const numComp = m_numComponents;
  forTargetSubRegions( mesh, [&]( localIndex const, ElementSubRegionBase & subRegion )
  {
    arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

    arrayView2d< real64 const, compflow::USD_COMP > const compDens =
      subRegion.getExtrinsicData< extrinsicMeshData::flow::globalCompDensity >();
    arrayView2d< real64, compflow::USD_COMP > const dCompDens =
      subRegion.getExtrinsicData< extrinsicMeshData::flow::deltaGlobalCompDensity >();

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      if( ghostRank[ei] < 0 )
      {
        for( localIndex ic = 0; ic < numComp; ++ic )
        {
          real64 const newDens = compDens[ei][ic] + dCompDens[ei][ic];
          if( newDens < minDensForDivision )
          {
            dCompDens[ei][ic] = -compDens[ei][ic] + minDensForDivision;
          }
        }
      }
    } );
  } );
}

void CompositionalMultiphaseBase::resetStateToBeginningOfStep( DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  forTargetSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( mesh, [&]( localIndex const targetIndex, auto & subRegion )
  {
    arrayView1d< real64 > const & dPres =
      subRegion.template getReference< array1d< real64 > >( extrinsicMeshData::flow::deltaPressure::key() );
    arrayView2d< real64, compflow::USD_COMP > const & dCompDens =
      subRegion.template getExtrinsicData< extrinsicMeshData::flow::deltaGlobalCompDensity >();

    dPres.zero();
    dCompDens.zero();

    // update porosity and permeability
    updatePorosityAndPermeability( subRegion, targetIndex );
    // update all fluid properties
    updateFluidState( subRegion, targetIndex );
  } );
}

void CompositionalMultiphaseBase::implicitStepComplete( real64 const & time,
                                                        real64 const & dt,
                                                        DomainPartition & domain )
{
  integer const numComp = m_numComponents;

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  // Step 1: save the aquifer converged state
  // note: we have to save the aquifer state **before** updating the pressure,
  // otherwise the aquifer flux is saved with the wrong pressure time level
  saveAquiferConvergedState( time, dt, domain );

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion )
  {
    arrayView1d< real64 const > const dPres =
      subRegion.getExtrinsicData< extrinsicMeshData::flow::deltaPressure >();
    arrayView2d< real64 const, compflow::USD_COMP > const dCompDens =
      subRegion.getExtrinsicData< extrinsicMeshData::flow::deltaGlobalCompDensity >();

    arrayView1d< real64 > const pres =
      subRegion.getExtrinsicData< extrinsicMeshData::flow::pressure >();
    arrayView2d< real64, compflow::USD_COMP > const compDens =
      subRegion.getExtrinsicData< extrinsicMeshData::flow::globalCompDensity >();

    // Step 2: increment the primary variables with the accumulated Newton updates
    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      pres[ei] += dPres[ei];
      for( localIndex ic = 0; ic < numComp; ++ic )
      {
        compDens[ei][ic] += dCompDens[ei][ic];
      }
    } );

    // Step 3: save the converged solid state (porosity, solid internal energy, etc)
    CoupledSolidBase const & porousMaterial = getConstitutiveModel< CoupledSolidBase >( subRegion, m_solidModelNames[targetIndex] );
    porousMaterial.saveConvergedState();

    // Step 4: if capillary pressure is supported, send the converged porosity and permeability to the capillary pressure model
    // note: this is needed when the capillary pressure depends on porosity and permeability (Leverett J-function for instance)
    if( m_capPressureFlag )
    {
      arrayView2d< real64 const > const porosity = porousMaterial.getPorosity();

      PermeabilityBase const & permeabilityMaterial =
        getConstitutiveModel< PermeabilityBase >( subRegion, m_permeabilityModelNames[targetIndex] );
      arrayView3d< real64 const > const permeability = permeabilityMaterial.permeability();

      CapillaryPressureBase const & capPressureMaterial =
        getConstitutiveModel< CapillaryPressureBase >( subRegion, m_capPressureModelNames[targetIndex] );
      capPressureMaterial.saveConvergedRockState( porosity, permeability );
    }

    // Step 5: if the thermal option is on, send the converged porosity and phase volume fraction to the thermal conductivity model
    // note: this is needed because the phaseVolFrac-weighted thermal conductivity treats phaseVolumeFraction explicitly for now
    if( m_thermalFlag )
    {
      arrayView2d< real64 const > const porosity = porousMaterial.getPorosity();

      arrayView2d< real64 const, compflow::USD_PHASE > const phaseVolFrac =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::phaseVolumeFraction >();

      ThermalConductivityBase const & thermalConductivityMaterial =
        getConstitutiveModel< ThermalConductivityBase >( subRegion, m_thermalConductivityModelNames[targetIndex] );
      thermalConductivityMaterial.saveConvergedRockFluidState( porosity, phaseVolFrac );
    }

  } );
}

void CompositionalMultiphaseBase::updateState( DomainPartition & domain )
{
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  forTargetSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( mesh, [&]( localIndex const targetIndex, auto & subRegion )
  {
    // update porosity and permeability
    updatePorosityAndPermeability( subRegion, targetIndex );
    // update all fluid properties
    updateFluidState( subRegion, targetIndex );
  } );
}

} // namespace geosx
