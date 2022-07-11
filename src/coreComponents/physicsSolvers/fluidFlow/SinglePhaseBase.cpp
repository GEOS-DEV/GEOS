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
 * @file SinglePhaseBase.cpp
 */

#include "SinglePhaseBase.hpp"


#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "constitutive/fluid/SingleFluidExtrinsicData.hpp"
#include "constitutive/fluid/singleFluidSelector.hpp"
#include "constitutive/permeability/PermeabilityExtrinsicData.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "constitutive/solid/SolidInternalEnergy.hpp"
#include "constitutive/thermalConductivity/singlePhaseThermalConductivitySelector.hpp"
#include "fieldSpecification/AquiferBoundaryCondition.hpp"
#include "fieldSpecification/EquilibriumInitialCondition.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fieldSpecification/SourceFluxBoundaryCondition.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "functions/TableFunction.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/ThermalSinglePhaseBaseKernels.hpp"


namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace singlePhaseBaseKernels;

SinglePhaseBase::SinglePhaseBase( const string & name,
                                  Group * const parent ):
  FlowSolverBase( name, parent )
{
  this->registerWrapper( viewKeyStruct::inputTemperatureString(), &m_inputTemperature ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Temperature" );
}


void SinglePhaseBase::registerDataOnMesh( Group & meshBodies )
{
  using namespace extrinsicMeshData::flow;

  FlowSolverBase::registerDataOnMesh( meshBodies );

  m_numDofPerCell = m_isThermal ? 2 : 1;

  forMeshTargets( meshBodies, [&] ( string const &,
                                    MeshLevel & mesh,
                                    arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                              [&]( localIndex const,
                                                                   ElementSubRegionBase & subRegion )
    {
      subRegion.registerExtrinsicData< pressure_n >( getName() );
      subRegion.registerExtrinsicData< initialPressure >( getName() );
      subRegion.registerExtrinsicData< pressure >( getName() );

      subRegion.registerExtrinsicData< bcPressure >( getName() ); // needed for the application of boundary conditions

      subRegion.registerExtrinsicData< deltaVolume >( getName() );

      subRegion.registerExtrinsicData< temperature >( getName() );
      subRegion.registerExtrinsicData< temperature_n >( getName() );

      subRegion.registerExtrinsicData< bcTemperature >( getName() ); // needed for the application of boundary conditions

      subRegion.registerExtrinsicData< mobility >( getName() );
      subRegion.registerExtrinsicData< dMobility_dPressure >( getName() );

      if( m_isThermal )
      {
        subRegion.registerExtrinsicData< dMobility_dTemperature >( getName() );
      }

    } );

    FaceManager & faceManager = mesh.getFaceManager();
    {
      faceManager.registerExtrinsicData< facePressure >( getName() );
    }
  } );
}

void SinglePhaseBase::setConstitutiveNamesCallSuper( ElementSubRegionBase & subRegion ) const
{
  FlowSolverBase::setConstitutiveNamesCallSuper( subRegion );
}

void SinglePhaseBase::setConstitutiveNames( ElementSubRegionBase & subRegion ) const
{
  string & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  fluidName = getConstitutiveName< SingleFluidBase >( subRegion );
  GEOSX_ERROR_IF( fluidName.empty(), GEOSX_FMT( "Fluid model not found on subregion {}", subRegion.getName() ) );

  if( m_isThermal )
  {
    string & thermalConductivityName = subRegion.registerWrapper< string >( viewKeyStruct::thermalConductivityNamesString() ).
                                         setPlotLevel( PlotLevel::NOPLOT ).
                                         setRestartFlags( RestartFlags::NO_WRITE ).
                                         setSizedFromParent( 0 ).
                                         setDescription( "Name of the thermal conductivity constitutive model to use" ).
                                         reference();

    thermalConductivityName = getConstitutiveName< SinglePhaseThermalConductivityBase >( subRegion );
    GEOSX_THROW_IF( thermalConductivityName.empty(),
                    GEOSX_FMT( "Thermal conductivity model not found on subregion {}", subRegion.getName() ),
                    InputError );
  }
}


void SinglePhaseBase::initializeAquiferBC() const
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  fsManager.forSubGroups< AquiferBoundaryCondition >( [&] ( AquiferBoundaryCondition & bc )
  {
    // set the gravity vector (needed later for the potential diff calculations)
    bc.setGravityVector( gravityVector() );
  } );
}

void SinglePhaseBase::validateConstitutiveModels( DomainPartition & domain ) const
{
  GEOSX_MARK_FUNCTION;

  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                               MeshLevel & mesh,
                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const,
                                                                  ElementSubRegionBase & subRegion )
    {
      string & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      fluidName = getConstitutiveName< SingleFluidBase >( subRegion );
      GEOSX_THROW_IF( fluidName.empty(),
                      GEOSX_FMT( "Fluid model not found on subregion {}", subRegion.getName() ),
                      InputError );

      SingleFluidBase const & fluid = getConstitutiveModel< SingleFluidBase >( subRegion, fluidName );

      constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
      {
        string const fluidModelName = castedFluid.catalogName();
        GEOSX_THROW_IF( m_isThermal && (fluidModelName != "ThermalCompressibleSinglePhaseFluid"),
                        GEOSX_FMT( "SingleFluidBase {}: the thermal option is enabled in the solver, but the fluid model `{}` is not for thermal fluid",
                                   getName(), fluid.getName() ),
                        InputError );
        GEOSX_THROW_IF( !m_isThermal && (fluidModelName == "ThermalCompressibleSinglePhaseFluid"),
                        GEOSX_FMT( "SingleFluidBase {}: the fluid model is for thermal fluid `{}`, but the solver option is incompatible with the fluid model",
                                   getName(), fluid.getName() ),
                        InputError );
      } );
    } );
  } );
}

SinglePhaseBase::FluidPropViews SinglePhaseBase::getFluidProperties( ConstitutiveBase const & fluid ) const
{
  SingleFluidBase const & singleFluid = dynamicCast< SingleFluidBase const & >( fluid );
  return { singleFluid.density(),
           singleFluid.dDensity_dPressure(),
           singleFluid.viscosity(),
           singleFluid.dViscosity_dPressure(),
           singleFluid.defaultDensity(),
           singleFluid.defaultViscosity() };
}

SinglePhaseBase::ThermalFluidPropViews SinglePhaseBase::getThermalFluidProperties( ConstitutiveBase const & fluid ) const
{
  SingleFluidBase const & singleFluid = dynamicCast< SingleFluidBase const & >( fluid );
  return { singleFluid.dDensity_dTemperature(),
           singleFluid.dViscosity_dTemperature() };
}

void SinglePhaseBase::initializePreSubGroups()
{
  FlowSolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  // 1. Validate various models against each other (must have same phases and components)
  validateConstitutiveModels( domain );

  // 2. Set the value of temperature
  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                               MeshLevel & mesh,
                                               arrayView1d< string const > const & regionNames )

  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 > const temp = subRegion.getExtrinsicData< extrinsicMeshData::flow::temperature >();
      temp.setValues< parallelHostPolicy >( m_inputTemperature );
    } );
  } );

  // 3. Initialize the aquifer boundary condition
  initializeAquiferBC();
}

void SinglePhaseBase::updateFluidModel( ObjectManagerBase & dataGroup ) const
{
  GEOSX_MARK_FUNCTION;

  arrayView1d< real64 const > const pres = dataGroup.getExtrinsicData< extrinsicMeshData::flow::pressure >();
  arrayView1d< real64 const > const temp = dataGroup.getExtrinsicData< extrinsicMeshData::flow::temperature >();

  SingleFluidBase & fluid =
    getConstitutiveModel< SingleFluidBase >( dataGroup, dataGroup.getReference< string >( viewKeyStruct::fluidNamesString() ) );

  constitutiveUpdatePassThru( fluid, [&]( auto & castedFluid )
  {
    typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();
    thermalSinglePhaseBaseKernels::FluidUpdateKernel::launch( fluidWrapper, pres, temp );
  } );
}

void SinglePhaseBase::updateSolidInternalEnergyModel( ObjectManagerBase & dataGroup ) const
{
  arrayView1d< real64 const > const temp = dataGroup.getExtrinsicData< extrinsicMeshData::flow::temperature >();

  string const & solidInternalEnergyName = dataGroup.getReference< string >( viewKeyStruct::solidInternalEnergyNamesString() );
  SolidInternalEnergy & solidInternalEnergy = getConstitutiveModel< SolidInternalEnergy >( dataGroup, solidInternalEnergyName );

  SolidInternalEnergy::KernelWrapper solidInternalEnergyWrapper = solidInternalEnergy.createKernelUpdates();

  thermalSinglePhaseBaseKernels::SolidInternalEnergyUpdateKernel::launch< parallelDevicePolicy<> >( dataGroup.size(), solidInternalEnergyWrapper, temp );
}

void SinglePhaseBase::updateFluidState( ObjectManagerBase & subRegion ) const
{
  updateFluidModel( subRegion );
  updateMobility( subRegion );
}

void SinglePhaseBase::updateMobility( ObjectManagerBase & dataGroup ) const
{
  GEOSX_MARK_FUNCTION;

  // output

  arrayView1d< real64 > const mob =
    dataGroup.getExtrinsicData< extrinsicMeshData::flow::mobility >();

  arrayView1d< real64 > const dMob_dPres =
    dataGroup.getExtrinsicData< extrinsicMeshData::flow::dMobility_dPressure >();

  // input

  SingleFluidBase & fluid =
    getConstitutiveModel< SingleFluidBase >( dataGroup, dataGroup.getReference< string >( viewKeyStruct::fluidNamesString() ) );
  FluidPropViews fluidProps = getFluidProperties( fluid );

  if( m_isThermal )
  {
    arrayView1d< real64 > const dMob_dTemp =
      dataGroup.getExtrinsicData< extrinsicMeshData::flow::dMobility_dTemperature >();

    ThermalFluidPropViews thermalFluidProps = getThermalFluidProperties( fluid );

    thermalSinglePhaseBaseKernels::MobilityKernel::launch< parallelDevicePolicy<> >( dataGroup.size(),
                                                                                     fluidProps.dens,
                                                                                     fluidProps.dDens_dPres,
                                                                                     thermalFluidProps.dDens_dTemp,
                                                                                     fluidProps.visc,
                                                                                     fluidProps.dVisc_dPres,
                                                                                     thermalFluidProps.dVisc_dTemp,
                                                                                     mob,
                                                                                     dMob_dPres,
                                                                                     dMob_dTemp );
  }
  else
  {
    singlePhaseBaseKernels::MobilityKernel::launch< parallelDevicePolicy<> >( dataGroup.size(),
                                                                              fluidProps.dens,
                                                                              fluidProps.dDens_dPres,
                                                                              fluidProps.visc,
                                                                              fluidProps.dVisc_dPres,
                                                                              mob,
                                                                              dMob_dPres );
  }
}

void SinglePhaseBase::initializePostInitialConditionsPreSubGroups()
{
  GEOSX_MARK_FUNCTION;

  FlowSolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );


  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                               MeshLevel & mesh,
                                               arrayView1d< string const > const & regionNames )
  {
    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addElementFields( { extrinsicMeshData::flow::pressure::key() },
                                     regionNames );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), false );

    // Moved the following part from ImplicitStepSetup to here since it only needs to be initialized once
    // They will be updated in applySystemSolution and ImplicitStepComplete, respectively
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                                                   auto & subRegion )
    {
      // Compute hydrostatic equilibrium in the regions for which corresponding field specification tag has been specified
      computeHydrostaticEquilibrium();

      // 1. update porosity, permeability, and density/viscosity
      SingleFluidBase const & fluid = getConstitutiveModel< SingleFluidBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::fluidNamesString() ) );

      updatePorosityAndPermeability( subRegion );
      updateFluidState( subRegion );

      // 2. save the initial density (for use in the single-phase poromechanics solver to compute the deltaBodyForce)
      fluid.initializeState();

      // 3. save the initial/old porosity
      CoupledSolidBase const & porousSolid = getConstitutiveModel< CoupledSolidBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::solidNamesString() ) );

      porousSolid.initializeState();

      // 4. initialize the rock thermal quantities: conductivity and solid internal energy
      if( m_isThermal )
      {
        // initialized porosity
        arrayView2d< real64 const > const porosity = porousSolid.getPorosity();

        string const & thermalConductivityName = subRegion.template getReference< string >( viewKeyStruct::thermalConductivityNamesString() );
        SinglePhaseThermalConductivityBase const & conductivityMaterial =
          getConstitutiveModel< SinglePhaseThermalConductivityBase >( subRegion, thermalConductivityName );
        conductivityMaterial.initializeRockFluidState( porosity );
        // note that there is nothing to update here because thermal conductivity is explicit for now

        updateSolidInternalEnergyModel( subRegion );
        string const & solidInternalEnergyName = subRegion.template getReference< string >( viewKeyStruct::solidInternalEnergyNamesString() );
        SolidInternalEnergy const & solidInternalEnergyMaterial =
          getConstitutiveModel< SolidInternalEnergy >( subRegion, solidInternalEnergyName );
        solidInternalEnergyMaterial.saveConvergedState();

      }
    } );

    mesh.getElemManager().forElementRegions< SurfaceElementRegion >( regionNames,
                                                                     [&]( localIndex const,
                                                                          SurfaceElementRegion & region )
    {
      region.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
      {
        ConstitutiveBase & fluid = getConstitutiveModel( subRegion, subRegion.getReference< string >( viewKeyStruct::fluidNamesString() )  );
        real64 const defaultDensity = getFluidProperties( fluid ).defaultDensity;

        subRegion.getWrapper< real64_array >( extrinsicMeshData::flow::hydraulicAperture::key() ).
          setApplyDefaultValue( region.getDefaultAperture() );

        subRegion.getWrapper< real64_array >( FaceElementSubRegion::viewKeyStruct::creationMassString() ).
          setApplyDefaultValue( defaultDensity * region.getDefaultAperture() );
      } );
    } );

    // Save initial pressure field (needed by the poromechanics solvers to compute the deltaPressure needed by the total stress)
    mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const,
                                                                  ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 const > const pres = subRegion.getExtrinsicData< extrinsicMeshData::flow::pressure >();
      arrayView1d< real64 > const presInit   = subRegion.getExtrinsicData< extrinsicMeshData::flow::initialPressure >();
      presInit.setValues< parallelDevicePolicy<> >( pres );
    } );

  } );
}

void SinglePhaseBase::computeHydrostaticEquilibrium()
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  real64 const gravVector[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );

  // Step 1: count individual equilibriums (there may be multiple ones)

  std::map< string, localIndex > equilNameToEquilId;
  localIndex equilCounter = 0;

  fsManager.forSubGroups< EquilibriumInitialCondition >( [&] ( EquilibriumInitialCondition const & bc )
  {
    // collect all the equil name to idx
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
  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel &,
                                                arrayView1d< string const > const & regionNames )
  {
    for( string const & regionName : regionNames )
    {
      regionFilter.insert( regionName );
    }
  } );

  // then start the actual table construction
  fsManager.apply< EquilibriumInitialCondition >( 0.0,
                                                  domain.getMeshBody( 0 ).getMeshLevel( 0 ),
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

    localIndex const equilIndex = equilNameToEquilId.at( fs.getName() );
    real64 const minElevation = LvArray::math::min( globalMinElevation[equilIndex], datumElevation );
    real64 const maxElevation = LvArray::math::max( globalMaxElevation[equilIndex], datumElevation );
    real64 const elevationIncrement = LvArray::math::min( fs.getElevationIncrement(), maxElevation - minElevation );
    localIndex const numPointsInTable = std::ceil( (maxElevation - minElevation) / elevationIncrement ) + 1;

    real64 const eps = 0.1 * (maxElevation - minElevation); // we add a small buffer to only log in the pathological cases
    GEOSX_LOG_RANK_0_IF( ( (datumElevation > globalMaxElevation[equilIndex]+eps)  || (datumElevation < globalMinElevation[equilIndex]-eps) ),
                         SinglePhaseBase::catalogName() << " " << getName()
                                                        << ": By looking at the elevation of the cell centers in this model, GEOSX found that "
                                                        << "the min elevation is " << globalMinElevation[equilIndex] << " and the max elevation is " << globalMaxElevation[equilIndex] << "\n"
                                                        << "But, a datum elevation of " << datumElevation << " was specified in the input file to equilibrate the model.\n "
                                                        << "The simulation is going to proceed with this out-of-bound datum elevation, but the initial condition may be inaccurate." );

    array1d< array1d< real64 > > elevationValues;
    array1d< real64 > pressureValues;
    elevationValues.resize( 1 );
    elevationValues[0].resize( numPointsInTable );
    pressureValues.resize( numPointsInTable );

    // Step 3.2: retrieve the fluid model to compute densities
    // we end up with the same issue as in applyDirichletBC: there is not a clean way to retrieve the fluid info

    // filter out region not in target
    Group const & region = subRegion.getParent().getParent();
    auto it = regionFilter.find( region.getName() );
    if( it == regionFilter.end() )
    {
      return; // the region is not in target, there is nothing to do
    }

    string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString());

    // filter out the proppant fluid constitutive models
    ConstitutiveBase & fluid = getConstitutiveModel( subRegion, fluidName );
    if( !dynamicCast< SingleFluidBase * >( &fluid ) )
    {
      return;
    }
    SingleFluidBase & singleFluid = dynamicCast< SingleFluidBase & >( fluid );

    // Step 3.3: compute the hydrostatic pressure values

    constitutiveUpdatePassThru( singleFluid, [&] ( auto & castedFluid )
    {
      using FluidType = TYPEOFREF( castedFluid );
      typename FluidType::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

      // note: inside this kernel, serialPolicy is used, and elevation/pressure values don't go to the GPU
      bool const equilHasConverged =
        HydrostaticPressureKernel::launch( numPointsInTable,
                                           maxNumEquilIterations,
                                           equilTolerance,
                                           gravVector,
                                           minElevation,
                                           elevationIncrement,
                                           datumElevation,
                                           datumPressure,
                                           fluidWrapper,
                                           elevationValues.toNestedView(),
                                           pressureValues.toView() );

      GEOSX_THROW_IF( !equilHasConverged,
                      SinglePhaseBase::catalogName() << " " << getName()
                                                     << ": hydrostatic pressure initialization failed to converge in region " << region.getName() << "!",
                      std::runtime_error );
    } );

    // Step 3.4: create hydrostatic pressure table

    FunctionManager & functionManager = FunctionManager::getInstance();

    string const tableName = fs.getName() + "_" + subRegion.getName() + "_table";
    TableFunction * const presTable = dynamicCast< TableFunction * >( functionManager.createChild( TableFunction::catalogName(), tableName ) );
    presTable->setTableCoordinates( elevationValues );
    presTable->setTableValues( pressureValues );
    presTable->setInterpolationMethod( TableFunction::InterpolationType::Linear );
    TableFunction::KernelWrapper presTableWrapper = presTable->createKernelWrapper();

    // Step 4: assign pressure as a function of elevation
    // TODO: this last step should probably be delayed to wait for the creation of FaceElements
    arrayView2d< real64 const > const elemCenter =
      subRegion.getReference< array2d< real64 > >( ElementSubRegionBase::viewKeyStruct::elementCenterString() );

    arrayView1d< real64 > const pres =
      subRegion.getReference< array1d< real64 > >( extrinsicMeshData::flow::pressure::key() );

    forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const i )
    {
      localIndex const k = targetSet[i];
      real64 const elevation = elemCenter[k][2];

      pres[k] = presTableWrapper.compute( &elevation );
    } );
  } );
}


real64 SinglePhaseBase::solverStep( real64 const & time_n,
                                    real64 const & dt,
                                    const int cycleNumber,
                                    DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  real64 dt_return;

  // setup dof numbers and linear system
  setupSystem( domain, m_dofManager, m_localMatrix, m_rhs, m_solution );

  implicitStepSetup( time_n, dt, domain );

  // currently the only method is implicit time integration
  dt_return = nonlinearImplicitStep( time_n, dt, cycleNumber, domain );

  // final step for completion of timestep. typically secondary variable updates and cleanup.
  implicitStepComplete( time_n, dt_return, domain );

  return dt_return;
}

void SinglePhaseBase::setupSystem( DomainPartition & domain,
                                   DofManager & dofManager,
                                   CRSMatrix< real64, globalIndex > & localMatrix,
                                   ParallelVector & rhs,
                                   ParallelVector & solution,
                                   bool const setSparsity )
{
  GEOSX_MARK_FUNCTION;

  SolverBase::setupSystem( domain,
                           dofManager,
                           localMatrix,
                           rhs,
                           solution,
                           setSparsity );
}

void SinglePhaseBase::implicitStepSetup( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                         real64 const & GEOSX_UNUSED_PARAM( dt ),
                                         DomainPartition & domain )
{
  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                               MeshLevel & mesh,
                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                                                   auto & subRegion )
    {
      arrayView1d< real64 const > const & pres = subRegion.template getExtrinsicData< extrinsicMeshData::flow::pressure >();
      arrayView1d< real64 > const & pres_n = subRegion.template getExtrinsicData< extrinsicMeshData::flow::pressure_n >();
      pres_n.setValues< parallelDevicePolicy<> >( pres );

      if( m_isThermal )
      {
        arrayView1d< real64 const > const & temp = subRegion.template getExtrinsicData< extrinsicMeshData::flow::temperature >();
        arrayView1d< real64 > const & temp_n = subRegion.template getExtrinsicData< extrinsicMeshData::flow::temperature_n >();
        temp_n.setValues< parallelDevicePolicy<> >( temp );
      }

      arrayView1d< real64 > const & dVol = subRegion.template getExtrinsicData< extrinsicMeshData::flow::deltaVolume >();
      dVol.zero();

      // This should fix NaN density in newly created fracture elements
      updatePorosityAndPermeability( subRegion );
      updateFluidState( subRegion );
      // for thermal simulations, update solid internal energy
      if( m_isThermal )
      {
        updateSolidInternalEnergyModel( subRegion );
      }

    } );

    mesh.getElemManager().forElementSubRegions< FaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          FaceElementSubRegion & subRegion )
    {
      arrayView1d< real64 const > const aper = subRegion.getExtrinsicData< extrinsicMeshData::flow::hydraulicAperture >();
      arrayView1d< real64 > const aper0 = subRegion.getExtrinsicData< extrinsicMeshData::flow::aperture0 >();
      aper0.setValues< parallelDevicePolicy<> >( aper );

      // Needed coz faceElems don't exist when initializing.
      CoupledSolidBase const & porousSolid = getConstitutiveModel< CoupledSolidBase >( subRegion, subRegion.getReference< string >( viewKeyStruct::solidNamesString() ) );
      porousSolid.saveConvergedState();

      updatePorosityAndPermeability( subRegion );
      updateFluidState( subRegion );

      // This call is required by the proppant solver, but should not be here
      SingleFluidBase const & fluid =
        getConstitutiveModel< SingleFluidBase >( subRegion, subRegion.getReference< string >( viewKeyStruct::fluidNamesString() ) );
      fluid.saveConvergedState();

    } );
  } );



}

void SinglePhaseBase::implicitStepComplete( real64 const & time,
                                            real64 const & dt,
                                            DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  // note: we have to save the aquifer state **before** updating the pressure,
  // otherwise the aquifer flux is saved with the wrong pressure time level
  saveAquiferConvergedState( time, dt, domain );

  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                               MeshLevel & mesh,
                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const,
                                                                  ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 const > const dVol = subRegion.getExtrinsicData< extrinsicMeshData::flow::deltaVolume >();
      arrayView1d< real64 > const vol = subRegion.getReference< array1d< real64 > >( CellElementSubRegion::viewKeyStruct::elementVolumeString() );

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
      {
        vol[ei] += dVol[ei];
      } );

      SingleFluidBase const & fluid =
        getConstitutiveModel< SingleFluidBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::fluidNamesString() ) );
      fluid.saveConvergedState();

      CoupledSolidBase const & porousSolid =
        getConstitutiveModel< CoupledSolidBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::solidNamesString() ) );
      porousSolid.saveConvergedState();

      if( m_isThermal )
      {
        arrayView2d< real64 const > const porosity = porousSolid.getPorosity();

        SinglePhaseThermalConductivityBase const & conductivityMaterial =
          getConstitutiveModel< SinglePhaseThermalConductivityBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::thermalConductivityNamesString() ) );

        conductivityMaterial.saveConvergedRockFluidState( porosity );
      }

    } );

    mesh.getElemManager().forElementSubRegions< FaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          FaceElementSubRegion & subRegion )
    {
      arrayView1d< integer const > const elemGhostRank = subRegion.ghostRank();
      arrayView1d< real64 const > const volume = subRegion.getElementVolume();
      arrayView1d< real64 > const creationMass = subRegion.getReference< real64_array >( FaceElementSubRegion::viewKeyStruct::creationMassString() );

      SingleFluidBase const & fluid =
        getConstitutiveModel< SingleFluidBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::fluidNamesString() ) );
      arrayView2d< real64 const > const density_n = fluid.density_n();

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
      {
        if( elemGhostRank[ei] < 0 )
        {
          if( volume[ei] * density_n[ei][0] > 1.1 * creationMass[ei] )
          {
            creationMass[ei] *= 0.75;
            if( creationMass[ei]<1.0e-20 )
            {
              creationMass[ei] = 0.0;
            }
          }
        }
      } );
    } );

  } );
}


void SinglePhaseBase::assembleSystem( real64 const time_n,
                                      real64 const dt,
                                      DomainPartition & domain,
                                      DofManager const & dofManager,
                                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                      arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  assembleAccumulationTerms( domain,
                             dofManager,
                             localMatrix,
                             localRhs );

  assembleFluxTerms( time_n,
                     dt,
                     domain,
                     dofManager,
                     localMatrix,
                     localRhs );
}

void SinglePhaseBase::assembleAccumulationTerms( DomainPartition & domain,
                                                 DofManager const & dofManager,
                                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                 arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                               MeshLevel & mesh,
                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion,
                                                SurfaceElementSubRegion >( regionNames,
                                                                           [&]( localIndex const,
                                                                                auto & subRegion )
    {
      SingleFluidBase const & fluid =
        getConstitutiveModel< SingleFluidBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::fluidNamesString() ) );
      //START_SPHINX_INCLUDE_COUPLEDSOLID
      CoupledSolidBase const & solid =
        getConstitutiveModel< CoupledSolidBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::solidNamesString() ) );
      //END_SPHINX_INCLUDE_COUPLEDSOLID

      string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

      if( m_isThermal )
      {
        thermalSinglePhaseBaseKernels::
          ElementBasedAssemblyKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( dofManager.rankOffset(),
                                                     dofKey,
                                                     subRegion,
                                                     fluid,
                                                     solid,
                                                     localMatrix,
                                                     localRhs );
      }
      else
      {
        singlePhaseBaseKernels::
          ElementBasedAssemblyKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( dofManager.rankOffset(),
                                                     dofKey,
                                                     subRegion,
                                                     fluid,
                                                     solid,
                                                     localMatrix,
                                                     localRhs );
      }
    } );
  } );
}

void SinglePhaseBase::applyBoundaryConditions( real64 time_n,
                                               real64 dt,
                                               DomainPartition & domain,
                                               DofManager const & dofManager,
                                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                               arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  applySourceFluxBC( time_n, dt, domain, dofManager, localMatrix, localRhs );
  applyDirichletBC( time_n, dt, domain, dofManager, localMatrix, localRhs );
  applyAquiferBC( time_n, dt, domain, dofManager, localMatrix, localRhs );
}

namespace internal
{
string const bcLogMessage = string( "SinglePhaseBase {}: at time {}s, " )
                            + string( "the <{}> boundary condition '{}' is applied to the element set '{}' in subRegion '{}'. " )
                            + string( "\nThe scale of this boundary condition is {} and multiplies the value of the provided function (if any). " )
                            + string( "\nThe total number of target elements (including ghost elements) is {}. " )
                            + string( "\nNote that if this number is equal to zero for all subRegions, the boundary condition will not be applied on this element set." );
}

void SinglePhaseBase::applyDirichletBC( real64 const time_n,
                                        real64 const dt,
                                        DomainPartition & domain,
                                        DofManager const & dofManager,
                                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                        arrayView1d< real64 > const & localRhs ) const
{
  GEOSX_MARK_FUNCTION;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & )
  {
    // 1. Apply pressure Dirichlet BCs
    fsManager.apply( time_n + dt,
                     mesh,
                     "ElementRegions",
                     extrinsicMeshData::flow::pressure::key(),
                     [&]( FieldSpecificationBase const & fs,
                          string const & setName,
                          SortedArrayView< localIndex const > const & lset,
                          Group & subRegion,
                          string const & )
    {
      if( fs.getLogLevel() >= 1 && m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
      {
        globalIndex const numTargetElems = MpiWrapper::sum< globalIndex >( lset.size() );
        GEOSX_LOG_RANK_0( GEOSX_FMT( geosx::internal::bcLogMessage,
                                     getName(), time_n+dt, FieldSpecificationBase::catalogName(),
                                     fs.getName(), setName, subRegion.getName(), fs.getScale(), numTargetElems ) );
      }

      // Specify the bc value of pressure
      fs.applyFieldValue< FieldSpecificationEqual,
                          parallelDevicePolicy<> >( lset,
                                                    time_n + dt,
                                                    subRegion,
                                                    extrinsicMeshData::flow::bcPressure::key() );
    } );

    // 2. Apply temperature Dirichlet BCs
    if( m_isThermal )
    {
      fsManager.apply( time_n + dt,
                       mesh,
                       "ElementRegions",
                       extrinsicMeshData::flow::temperature::key(),
                       [&]( FieldSpecificationBase const & fs,
                            string const &,
                            SortedArrayView< localIndex const > const & lset,
                            Group & subRegion,
                            string const & )
      {
        // Specify the bc value of temperature
        fs.applyFieldValue< FieldSpecificationEqual,
                            parallelDevicePolicy<> >( lset,
                                                      time_n + dt,
                                                      subRegion,
                                                      extrinsicMeshData::flow::bcTemperature::key() );
      } );
    }

    globalIndex const rankOffset = dofManager.rankOffset();

    // 3. Apply the pressure bc to the system
    fsManager.apply( time_n + dt,
                     mesh,
                     "ElementRegions",
                     extrinsicMeshData::flow::pressure::key(),
                     [&]( FieldSpecificationBase const & GEOSX_UNUSED_PARAM( fs ),
                          string const &,
                          SortedArrayView< localIndex const > const & lset,
                          Group & subRegion,
                          string const & )
    {
      arrayView1d< real64 const > const bcPres =
        subRegion.getReference< array1d< real64 > >( extrinsicMeshData::flow::bcPressure::key() );

      arrayView1d< integer const > const ghostRank =
        subRegion.getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString() );
      arrayView1d< globalIndex const > const dofNumber =
        subRegion.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< real64 const > const pres =
        subRegion.getReference< array1d< real64 > >( extrinsicMeshData::flow::pressure::key() );

      forAll< parallelDevicePolicy<> >( lset.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
      {
        localIndex const ei = lset[a];
        if( ghostRank[ei] >= 0 )
        {
          return;
        }

        globalIndex const dofIndex = dofNumber[ei];
        localIndex const localRow = dofIndex - rankOffset;
        real64 rhsValue;

        // Apply pressure value to the matrix/rhs
        FieldSpecificationEqual::SpecifyFieldValue( dofIndex,
                                                    rankOffset,
                                                    localMatrix,
                                                    rhsValue,
                                                    bcPres[ei],
                                                    pres[ei] );
        localRhs[localRow] = rhsValue;
      } );
    } );

    // 3. Apply the temperature bc to the system
    if( m_isThermal )
    {
      fsManager.apply( time_n + dt,
                       mesh,
                       "ElementRegions",
                       extrinsicMeshData::flow::temperature::key(),
                       [&]( FieldSpecificationBase const & GEOSX_UNUSED_PARAM( fs ),
                            string const &,
                            SortedArrayView< localIndex const > const & lset,
                            Group & subRegion,
                            string const & )
      {
        arrayView1d< real64 const > const bcTemp =
          subRegion.getReference< array1d< real64 > >( extrinsicMeshData::flow::bcTemperature::key() );

        arrayView1d< integer const > const ghostRank =
          subRegion.getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString() );
        arrayView1d< globalIndex const > const dofNumber =
          subRegion.getReference< array1d< globalIndex > >( dofKey );
        arrayView1d< real64 const > const temp =
          subRegion.getReference< array1d< real64 > >( extrinsicMeshData::flow::temperature::key() );

        forAll< parallelDevicePolicy<> >( lset.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
        {
          localIndex const ei = lset[a];
          if( ghostRank[ei] >= 0 )
          {
            return;
          }

          globalIndex const dofIndex = dofNumber[ei];
          localIndex const localRow = dofIndex - rankOffset;
          real64 rhsValue;

          // Apply temperature value to the matrix/rhs
          FieldSpecificationEqual::SpecifyFieldValue( dofIndex + 1,
                                                      rankOffset,
                                                      localMatrix,
                                                      rhsValue,
                                                      bcTemp[ei],
                                                      temp[ei] );
          localRhs[localRow + 1] = rhsValue;
        } );
      } );
    }
  } );
}

void SinglePhaseBase::applySourceFluxBC( real64 const time_n,
                                         real64 const dt,
                                         DomainPartition & domain,
                                         DofManager const & dofManager,
                                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                         arrayView1d< real64 > const & localRhs ) const
{
  GEOSX_MARK_FUNCTION;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  // Step 1: count individual source flux boundary conditions

  std::map< string, localIndex > bcNameToBcId;
  localIndex bcCounter = 0;

  fsManager.forSubGroups< SourceFluxBoundaryCondition >( [&] ( SourceFluxBoundaryCondition const & bc )
  {
    // collect all the bc names to idx
    bcNameToBcId[bc.getName()] = bcCounter;
    bcCounter++;
  } );

  if( bcCounter == 0 )
  {
    return;
  }

  // Step 2: count the set size for each source flux (each source flux may have multiple target sets)

  array1d< globalIndex > bcAllSetsSize( bcNameToBcId.size() );

  computeSourceFluxSizeScalingFactor( time_n,
                                      dt,
                                      domain,
                                      bcNameToBcId,
                                      bcAllSetsSize.toView() );

  // Step 3: we are ready to impose the boundary condition, normalized by the set size

  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                               MeshLevel & mesh,
                                               arrayView1d< string const > const & )
  {
    fsManager.apply( time_n + dt,
                     mesh,
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
                                     getName(), time_n+dt, SourceFluxBoundaryCondition::catalogName(),
                                     fs.getName(), setName, subRegion.getName(), fs.getScale(), numTargetElems ) );
      }

      if( targetSet.size() == 0 )
      {
        return;
      }

      arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< integer const > const ghostRank =
        subRegion.getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString() );

      // Step 3.1: get the values of the source boundary condition that need to be added to the rhs

      array1d< globalIndex > dofArray( targetSet.size() );
      array1d< real64 > rhsContributionArray( targetSet.size() );
      arrayView1d< real64 > rhsContributionArrayView = rhsContributionArray.toView();
      localIndex const rankOffset = dofManager.rankOffset();

      // note that the dofArray will not be used after this step (simpler to use dofNumber instead)
      fs.computeRhsContribution< FieldSpecificationAdd,
                                 parallelDevicePolicy<> >( targetSet.toViewConst(),
                                                           time_n + dt,
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

      // Step 3.2: we are ready to add the right-hand side contributions, taking into account our equation layout

      // get the normalizer
      real64 const sizeScalingFactor = bcAllSetsSize[bcNameToBcId.at( fs.getName())];

      forAll< parallelDevicePolicy<> >( targetSet.size(), [sizeScalingFactor,
                                                           targetSet,
                                                           rankOffset,
                                                           ghostRank,
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

        // add the value to the mass balance equation
        globalIndex const rowIndex = dofNumber[ei] - rankOffset;
        localRhs[rowIndex] += rhsContributionArrayView[a] / sizeScalingFactor; // scale the contribution by the sizeScalingFactor here!
      } );
    } );
  } );
}

void SinglePhaseBase::updateState( DomainPartition & domain )
{

// set mass fraction flag on fluid models
  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                               MeshLevel & mesh,
                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                                                   auto & subRegion )
    {
      updatePorosityAndPermeability( subRegion );
      updateFluidState( subRegion );

      if( m_isThermal )
      {
        updateSolidInternalEnergyModel( subRegion );
      }
    } );
  } );
}

void SinglePhaseBase::solveLinearSystem( DofManager const & dofManager,
                                         ParallelMatrix & matrix,
                                         ParallelVector & rhs,
                                         ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  rhs.scale( -1.0 );
  solution.zero();

  SolverBase::solveLinearSystem( dofManager, matrix, rhs, solution );
}

void SinglePhaseBase::resetStateToBeginningOfStep( DomainPartition & domain )
{
  // set mass fraction flag on fluid models
  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                               MeshLevel & mesh,
                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                                                   auto & subRegion )
    {
      arrayView1d< real64 > const pres = subRegion.template getExtrinsicData< extrinsicMeshData::flow::pressure >();
      arrayView1d< real64 const > const pres_n = subRegion.template getExtrinsicData< extrinsicMeshData::flow::pressure_n >();
      pres.setValues< parallelDevicePolicy<> >( pres_n );

      if( m_isThermal )
      {
        arrayView1d< real64 > const temp = subRegion.template getExtrinsicData< extrinsicMeshData::flow::temperature >();
        arrayView1d< real64 const > const temp_n = subRegion.template getExtrinsicData< extrinsicMeshData::flow::temperature_n >();
        temp.setValues< parallelDevicePolicy<> >( temp_n );
      }

      updatePorosityAndPermeability( subRegion );
      updateFluidState( subRegion );

      if( m_isThermal )
      {
        updateSolidInternalEnergyModel( subRegion );
      }
    } );
  } );
}

} /* namespace geosx */
