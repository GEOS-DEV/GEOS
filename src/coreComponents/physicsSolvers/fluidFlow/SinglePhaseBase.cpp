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
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidFields.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidSelector.hpp"
#include "constitutive/permeability/PermeabilityFields.hpp"
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
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseFields.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/ThermalSinglePhaseBaseKernels.hpp"


namespace geos
{

using namespace dataRepository;
using namespace constitutive;
using namespace singlePhaseBaseKernels;

SinglePhaseBase::SinglePhaseBase( const string & name,
                                  Group * const parent ):
  FlowSolverBase( name, parent ),
  m_keepFlowVariablesConstantDuringInitStep( 0 )
{
  this->registerWrapper( viewKeyStruct::inputTemperatureString(), &m_inputTemperature ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Temperature" );
}


void SinglePhaseBase::registerDataOnMesh( Group & meshBodies )
{
  using namespace fields::flow;

  FlowSolverBase::registerDataOnMesh( meshBodies );

  m_numDofPerCell = m_isThermal ? 2 : 1;

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                              [&]( localIndex const,
                                                                   ElementSubRegionBase & subRegion )
    {
      subRegion.registerField< pressure >( getName() );
      subRegion.registerField< pressure_n >( getName() );
      subRegion.registerField< initialPressure >( getName() );
      subRegion.registerField< deltaPressure >( getName() ); // for reporting/stats purposes
      subRegion.registerField< bcPressure >( getName() ); // needed for the application of boundary conditions
      if( m_isFixedStressPoromechanicsUpdate )
      {
        subRegion.registerField< pressure_k >( getName() ); // needed for the fixed-stress porosity update
      }

      subRegion.registerField< deltaVolume >( getName() );

      subRegion.registerField< temperature >( getName() );
      subRegion.registerField< temperature_n >( getName() );
      subRegion.registerField< initialTemperature >( getName() );
      subRegion.registerField< bcTemperature >( getName() ); // needed for the application of boundary conditions
      if( m_isFixedStressPoromechanicsUpdate )
      {
        subRegion.registerField< temperature_k >( getName() ); // needed for the fixed-stress porosity update
      }

      subRegion.registerField< mobility >( getName() );
      subRegion.registerField< dMobility_dPressure >( getName() );

      if( m_isThermal )
      {
        subRegion.registerField< dMobility_dTemperature >( getName() );
      }

    } );

    FaceManager & faceManager = mesh.getFaceManager();
    {
      faceManager.registerField< facePressure >( getName() );

      if( m_isThermal )
      {
        faceManager.registerField< faceTemperature >( getName() );
      }
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
  GEOS_ERROR_IF( fluidName.empty(), GEOS_FMT( "Fluid model not found on subregion {}", subRegion.getName() ) );

  if( m_isThermal )
  {
    string & thermalConductivityName = subRegion.registerWrapper< string >( viewKeyStruct::thermalConductivityNamesString() ).
                                         setPlotLevel( PlotLevel::NOPLOT ).
                                         setRestartFlags( RestartFlags::NO_WRITE ).
                                         setSizedFromParent( 0 ).
                                         setDescription( "Name of the thermal conductivity constitutive model to use" ).
                                         reference();

    thermalConductivityName = getConstitutiveName< SinglePhaseThermalConductivityBase >( subRegion );
    GEOS_THROW_IF( thermalConductivityName.empty(),
                   GEOS_FMT( "Thermal conductivity model not found on subregion {}", subRegion.getName() ),
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
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const,
                                                                  ElementSubRegionBase & subRegion )
    {
      string & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      fluidName = getConstitutiveName< SingleFluidBase >( subRegion );
      GEOS_THROW_IF( fluidName.empty(),
                     GEOS_FMT( "Fluid model not found on subregion {}", subRegion.getName() ),
                     InputError );

      SingleFluidBase const & fluid = getConstitutiveModel< SingleFluidBase >( subRegion, fluidName );

      constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
      {
        string const fluidModelName = castedFluid.catalogName();
        GEOS_THROW_IF( m_isThermal && (fluidModelName != "ThermalCompressibleSinglePhaseFluid"),
                       GEOS_FMT( "SingleFluidBase {}: the thermal option is enabled in the solver, but the fluid model `{}` is not for thermal fluid",
                                 getName(), fluid.getName() ),
                       InputError );
        GEOS_THROW_IF( !m_isThermal && (fluidModelName == "ThermalCompressibleSinglePhaseFluid"),
                       GEOS_FMT( "SingleFluidBase {}: the fluid model is for thermal fluid `{}`, but the solver option is incompatible with the fluid model",
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
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )

  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 > const temp = subRegion.getField< fields::flow::temperature >();
      temp.setValues< parallelHostPolicy >( m_inputTemperature );
    } );
  } );

  // 3. Initialize the aquifer boundary condition
  initializeAquiferBC();
}

void SinglePhaseBase::updateFluidModel( ObjectManagerBase & dataGroup ) const
{
  GEOS_MARK_FUNCTION;

  arrayView1d< real64 const > const pres = dataGroup.getField< fields::flow::pressure >();
  arrayView1d< real64 const > const temp = dataGroup.getField< fields::flow::temperature >();

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
  arrayView1d< real64 const > const temp = dataGroup.getField< fields::flow::temperature >();

  string const & solidInternalEnergyName = dataGroup.getReference< string >( viewKeyStruct::solidInternalEnergyNamesString() );
  SolidInternalEnergy & solidInternalEnergy = getConstitutiveModel< SolidInternalEnergy >( dataGroup, solidInternalEnergyName );

  SolidInternalEnergy::KernelWrapper solidInternalEnergyWrapper = solidInternalEnergy.createKernelUpdates();

  thermalSinglePhaseBaseKernels::SolidInternalEnergyUpdateKernel::launch< parallelDevicePolicy<> >( dataGroup.size(), solidInternalEnergyWrapper, temp );
}

void SinglePhaseBase::updateThermalConductivity( ElementSubRegionBase & subRegion ) const
{
  CoupledSolidBase const & porousSolid =
    getConstitutiveModel< CoupledSolidBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::solidNamesString() ) );

  arrayView2d< real64 const > const porosity = porousSolid.getPorosity();

  string const & thermalConductivityName = subRegion.template getReference< string >( viewKeyStruct::thermalConductivityNamesString() );
  SinglePhaseThermalConductivityBase const & conductivityMaterial =
    getConstitutiveModel< SinglePhaseThermalConductivityBase >( subRegion, thermalConductivityName );
  conductivityMaterial.update( porosity );
}

void SinglePhaseBase::updateFluidState( ObjectManagerBase & subRegion ) const
{
  updateFluidModel( subRegion );
  updateMobility( subRegion );
}

void SinglePhaseBase::updateMobility( ObjectManagerBase & dataGroup ) const
{
  GEOS_MARK_FUNCTION;

  // output

  arrayView1d< real64 > const mob = dataGroup.getField< fields::flow::mobility >();
  arrayView1d< real64 > const dMob_dPres = dataGroup.getField< fields::flow::dMobility_dPressure >();

  // input

  SingleFluidBase & fluid =
    getConstitutiveModel< SingleFluidBase >( dataGroup, dataGroup.getReference< string >( viewKeyStruct::fluidNamesString() ) );
  FluidPropViews fluidProps = getFluidProperties( fluid );

  if( m_isThermal )
  {
    arrayView1d< real64 > const dMob_dTemp =
      dataGroup.getField< fields::flow::dMobility_dTemperature >();

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
  GEOS_MARK_FUNCTION;

  FlowSolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );


  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addElementFields( { fields::flow::pressure::key() },
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
      // In addition, to avoid multiplying permeability/porosity bay netToGross in the assembly kernel, we do it once and for all here
      arrayView1d< real64 const > const netToGross = subRegion.template getField< fields::flow::netToGross >();
      CoupledSolidBase const & porousSolid =
        getConstitutiveModel< CoupledSolidBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::solidNamesString() ) );
      PermeabilityBase const & permeabilityModel =
        getConstitutiveModel< PermeabilityBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::permeabilityNamesString() ) );
      permeabilityModel.scaleHorizontalPermeability( netToGross );
      porousSolid.scaleReferencePorosity( netToGross );
      saveConvergedState( subRegion ); // necessary for a meaningful porosity update in sequential schemes
      updatePorosityAndPermeability( subRegion );

      SingleFluidBase const & fluid =
        getConstitutiveModel< SingleFluidBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::fluidNamesString() ) );
      updateFluidState( subRegion );

      // 2. save the initial density (for use in the single-phase poromechanics solver to compute the deltaBodyForce)
      fluid.initializeState();

      // 3. save the initial/old porosity
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

        subRegion.getWrapper< real64_array >( fields::flow::hydraulicAperture::key() ).
          setApplyDefaultValue( region.getDefaultAperture() );

        subRegion.getWrapper< real64_array >( FaceElementSubRegion::viewKeyStruct::creationMassString() ).
          setApplyDefaultValue( defaultDensity * region.getDefaultAperture() );
      } );
    } );

    // Save initial pressure field
    mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const,
                                                                  ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 const > const pres = subRegion.getField< fields::flow::pressure >();
      arrayView1d< real64 > const initPres = subRegion.getField< fields::flow::initialPressure >();
      arrayView1d< real64 const > const & temp = subRegion.template getField< fields::flow::temperature >();
      arrayView1d< real64 > const initTemp = subRegion.template getField< fields::flow::initialTemperature >();
      initPres.setValues< parallelDevicePolicy<> >( pres );
      initTemp.setValues< parallelDevicePolicy<> >( temp );

      if (m_isFixedStressPoromechanicsUpdate)
      {
        arrayView1d< real64 > const pres_k = subRegion.getField < fields::flow::pressure_k >();
        pres_k.setValues< parallelDevicePolicy<> >( pres );
        arrayView1d< real64 > const temp_k = subRegion.getField < fields::flow::temperature_k >();
        temp_k.setValues< parallelDevicePolicy<> >( temp );

        std::cout << pres_k[0] << "\t" << pres[0] << "\t" << initPres[0] << std::endl;
      }

    } );
  } );

  // report to the user if some pore volumes are very small
  // note: this function is here because: 1) porosity has been initialized and 2) NTG has been applied
  validatePoreVolumes( domain );
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
    GEOS_THROW_IF( !isZero( gravVector[0] ) || !isZero( gravVector[1] ),
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
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel &,
                                                                arrayView1d< string const > const & regionNames )
  {
    for( string const & regionName : regionNames )
    {
      regionFilter.insert( regionName );
    }
  } );

  // then start the actual table construction
  fsManager.apply< ElementSubRegionBase,
                   EquilibriumInitialCondition >( 0.0,
                                                  domain.getMeshBody( 0 ).getBaseDiscretization(),
                                                  EquilibriumInitialCondition::catalogName(),
                                                  [&] ( EquilibriumInitialCondition const & fs,
                                                        string const &,
                                                        SortedArrayView< localIndex const > const & targetSet,
                                                        ElementSubRegionBase & subRegion,
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
    localIndex const numPointsInTable = ( elevationIncrement > 0 ) ? std::ceil( (maxElevation - minElevation) / elevationIncrement ) + 1 : 1;

    real64 const eps = 0.1 * (maxElevation - minElevation); // we add a small buffer to only log in the pathological cases
    GEOS_LOG_RANK_0_IF( ( (datumElevation > globalMaxElevation[equilIndex]+eps)  || (datumElevation < globalMinElevation[equilIndex]-eps) ),
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

      GEOS_THROW_IF( !equilHasConverged,
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
    arrayView2d< real64 const > const elemCenter = subRegion.getElementCenter();
    arrayView1d< real64 > const pres = subRegion.getField< fields::flow::pressure >();

    RAJA::ReduceMin< parallelDeviceReduce, real64 > minPressure( LvArray::NumericLimits< real64 >::max );

    forAll< parallelDevicePolicy< > >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const i )
    {
      localIndex const k = targetSet[i];
      real64 const elevation = elemCenter[k][2];
      pres[k] = presTableWrapper.compute( &elevation );
      minPressure.min( pres[k] );
    } );

    // For single phase flow, just issue a warning, because the simulation can proceed with a negative pressure
    GEOS_WARNING_IF( minPressure.get() <= 0.0,
                     GEOS_FMT( "A negative pressure of {} Pa was found during hydrostatic initialization in region/subRegion {}/{}",
                               minPressure.get(), region.getName(), subRegion.getName() ) );
  } );
}

void SinglePhaseBase::implicitStepSetup( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                         real64 const & GEOS_UNUSED_PARAM( dt ),
                                         DomainPartition & domain )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                                                   auto & subRegion )
    {
      arrayView1d< real64 const > const & pres = subRegion.template getField< fields::flow::pressure >();
      arrayView1d< real64 const > const & initPres = subRegion.template getField< fields::flow::initialPressure >();
      arrayView1d< real64 > const & deltaPres = subRegion.template getField< fields::flow::deltaPressure >();

      singlePhaseBaseKernels::StatisticsKernel::
        saveDeltaPressure( subRegion.size(), pres, initPres, deltaPres );
      saveConvergedState( subRegion );

      arrayView1d< real64 > const & dVol = subRegion.template getField< fields::flow::deltaVolume >();
      dVol.zero();

      // This should fix NaN density in newly created fracture elements
      updatePorosityAndPermeability( subRegion );
      updateFluidState( subRegion );
      // for thermal simulations, update solid internal energy
      if( m_isThermal )
      {
        updateSolidInternalEnergyModel( subRegion );
        updateThermalConductivity( subRegion );

      }

    } );

    mesh.getElemManager().forElementSubRegions< FaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          FaceElementSubRegion & subRegion )
    {
      arrayView1d< real64 const > const aper = subRegion.getField< fields::flow::hydraulicAperture >();
      arrayView1d< real64 > const aper0 = subRegion.getField< fields::flow::aperture0 >();
      aper0.setValues< parallelDevicePolicy<> >( aper );

      // Needed coz faceElems don't exist when initializing.
      CoupledSolidBase const & porousSolid = getConstitutiveModel< CoupledSolidBase >( subRegion, subRegion.getReference< string >( viewKeyStruct::solidNamesString() ) );
      porousSolid.saveConvergedState();

      saveConvergedState( subRegion ); // necessary for a meaningful porosity update in sequential schemes
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
  GEOS_MARK_FUNCTION;

  // note: we have to save the aquifer state **before** updating the pressure,
  // otherwise the aquifer flux is saved with the wrong pressure time level
  saveAquiferConvergedState( time, dt, domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const,
                                                                  ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 const > const dVol = subRegion.getField< fields::flow::deltaVolume >();
      arrayView1d< real64 > const vol = subRegion.getReference< array1d< real64 > >( CellElementSubRegion::viewKeyStruct::elementVolumeString() );

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
      {
        vol[ei] += dVol[ei];
      } );

      SingleFluidBase const & fluid =
        getConstitutiveModel< SingleFluidBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::fluidNamesString() ) );
      fluid.saveConvergedState();

      CoupledSolidBase const & porousSolid =
        getConstitutiveModel< CoupledSolidBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::solidNamesString() ) );
      if( m_keepFlowVariablesConstantDuringInitStep )
      {
        porousSolid.ignoreConvergedState(); // newPorosity <- porosity_n
      }
      else
      {
        porousSolid.saveConvergedState(); // porosity_n <- porosity
      }

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

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
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
  GEOS_MARK_FUNCTION;

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
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion,
                                                SurfaceElementSubRegion >( regionNames,
                                                                           [&]( localIndex const,
                                                                                auto & subRegion )
    {
      accumulationAssemblyLaunch( dofManager, subRegion, localMatrix, localRhs );
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
  GEOS_MARK_FUNCTION;

  if( m_keepFlowVariablesConstantDuringInitStep )
  {
    // this function is going to force the current flow state to be constant during the time step
    // this is used when the poromechanics solver is performing the stress initialization
    // TODO: in the future, a dedicated poromechanics kernel should eliminate the flow vars to construct a reduced system
    //       which will remove the need for this brittle passing aroung of flag
    keepFlowVariablesConstantDuringInitStep( time_n, dt, dofManager, domain, localMatrix.toViewConstSizes(), localRhs.toView() );
  }
  else
  {
    applySourceFluxBC( time_n, dt, domain, dofManager, localMatrix, localRhs );
    applyDirichletBC( time_n, dt, domain, dofManager, localMatrix, localRhs );
    applyAquiferBC( time_n, dt, domain, dofManager, localMatrix, localRhs );
  }
}

namespace
{

char const bcLogMessage[] =
  "SinglePhaseBase {}: at time {}s, "
  "the <{}> boundary condition '{}' is applied to the element set '{}' in subRegion '{}'. "
  "\nThe scale of this boundary condition is {} and multiplies the value of the provided function (if any). "
  "\nThe total number of target elements (including ghost elements) is {}. "
  "\nNote that if this number is equal to zero for all subRegions, the boundary condition will not be applied on this element set.";

void applyAndSpecifyFieldValue( real64 const & time_n,
                                real64 const & dt,
                                MeshLevel & mesh,
                                globalIndex const rankOffset,
                                string const dofKey,
                                bool const isFirstNonlinearIteration,
                                string const solverName,
                                integer const idof,
                                string const fieldKey,
                                string const boundaryFieldKey,
                                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                arrayView1d< real64 > const & localRhs )
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  fsManager.apply< ElementSubRegionBase >( time_n + dt,
                                           mesh,
                                           fieldKey,
                                           [&]( FieldSpecificationBase const & fs,
                                                string const & setName,
                                                SortedArrayView< localIndex const > const & lset,
                                                ElementSubRegionBase & subRegion,
                                                string const & )
  {
    if( fs.getLogLevel() >= 1 && isFirstNonlinearIteration )
    {
      globalIndex const numTargetElems = MpiWrapper::sum< globalIndex >( lset.size() );
      GEOS_LOG_RANK_0( GEOS_FMT( bcLogMessage,
                                 solverName, time_n+dt, FieldSpecificationBase::catalogName(),
                                 fs.getName(), setName, subRegion.getName(), fs.getScale(), numTargetElems ) );
    }

    // Specify the bc value of the field
    fs.applyFieldValue< FieldSpecificationEqual,
                        parallelDevicePolicy<> >( lset,
                                                  time_n + dt,
                                                  subRegion,
                                                  boundaryFieldKey );

    arrayView1d< integer const > const ghostRank = subRegion.ghostRank();
    arrayView1d< globalIndex const > const dofNumber =
      subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< real64 const > const bcField =
      subRegion.getReference< array1d< real64 > >( boundaryFieldKey );
    arrayView1d< real64 const > const field =
      subRegion.getReference< array1d< real64 > >( fieldKey );

    forAll< parallelDevicePolicy<> >( lset.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      localIndex const ei = lset[a];
      if( ghostRank[ei] >= 0 )
      {
        return;
      }

      globalIndex const dofIndex = dofNumber[ei];
      localIndex const localRow = dofIndex - rankOffset;
      real64 rhsValue;

      // Apply field value to the matrix/rhs
      FieldSpecificationEqual::SpecifyFieldValue( dofIndex + idof,
                                                  rankOffset,
                                                  localMatrix,
                                                  rhsValue,
                                                  bcField[ei],
                                                  field[ei] );
      localRhs[localRow + idof] = rhsValue;
    } );
  } );
}

}

void SinglePhaseBase::applyDirichletBC( real64 const time_n,
                                        real64 const dt,
                                        DomainPartition & domain,
                                        DofManager const & dofManager,
                                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                        arrayView1d< real64 > const & localRhs ) const
{
  GEOS_MARK_FUNCTION;

  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
  globalIndex const rankOffset = dofManager.rankOffset();
  bool const isFirstNonlinearIteration = ( m_nonlinearSolverParameters.m_numNewtonIterations == 0 );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    applyAndSpecifyFieldValue( time_n, dt, mesh, rankOffset, dofKey, isFirstNonlinearIteration, getName(),
                               0, fields::flow::pressure::key(), fields::flow::bcPressure::key(),
                               localMatrix, localRhs );
    if( m_isThermal )
    {
      applyAndSpecifyFieldValue( time_n, dt, mesh, rankOffset, dofKey, isFirstNonlinearIteration, getName(),
                                 1, fields::flow::temperature::key(), fields::flow::bcTemperature::key(),
                                 localMatrix, localRhs );
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
  GEOS_MARK_FUNCTION;

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

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & )
  {
    integer const isThermal = m_isThermal;

    fsManager.apply< ElementSubRegionBase,
                     SourceFluxBoundaryCondition >( time_n + dt,
                                                    mesh,
                                                    SourceFluxBoundaryCondition::catalogName(),
                                                    [&, isThermal]( SourceFluxBoundaryCondition const & fs,
                                                                    string const & setName,
                                                                    SortedArrayView< localIndex const > const & targetSet,
                                                                    ElementSubRegionBase & subRegion,
                                                                    string const & )
    {
      if( fs.getLogLevel() >= 1 && m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
      {
        globalIndex const numTargetElems = MpiWrapper::sum< globalIndex >( targetSet.size() );
        GEOS_LOG_RANK_0( GEOS_FMT( bcLogMessage,
                                   getName(), time_n+dt, SourceFluxBoundaryCondition::catalogName(),
                                   fs.getName(), setName, subRegion.getName(), fs.getScale(), numTargetElems ) );

        if( isThermal )
        {
          char const msg[] = "SinglePhaseBase {} with isThermal = 1. At time {}s, "
                             "the <{}> source flux boundary condition '{}' will be applied with the following behavior"
                             "\n - negative value (injection): the mass balance equation is modified to considered the additional source term"
                             "\n - positive value (production): both the mass balance and the energy balance equations are modified to considered the additional source term. " \
                             "\n For the energy balance equation, the mass flux is multipied by the enthalpy in the cell from which the fluid is being produced.";
          GEOS_LOG_RANK_0( GEOS_FMT( msg,
                                     getName(), time_n+dt, SourceFluxBoundaryCondition::catalogName(), fs.getName() ) );
        }
      }

      if( targetSet.size() == 0 )
      {
        return;
      }

      arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

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
                                                           [] GEOS_HOST_DEVICE ( localIndex const )
      {
        return 0.0;
      } );

      // Step 3.2: we are ready to add the right-hand side contributions, taking into account our equation layout

      // get the normalizer
      real64 const sizeScalingFactor = bcAllSetsSize[bcNameToBcId.at( fs.getName())];

      if( isThermal )
      {
        SingleFluidBase const & fluid =
          getConstitutiveModel< SingleFluidBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::fluidNamesString() ) );

        arrayView2d< real64 const > const enthalpy = fluid.enthalpy();
        arrayView2d< real64 const > const dEnthalpy_dTemperature = fluid.dEnthalpy_dTemperature();
        arrayView2d< real64 const > const dEnthalpy_dPressure    = fluid.dEnthalpy_dPressure();
        forAll< parallelDevicePolicy<> >( targetSet.size(), [sizeScalingFactor,
                                                             targetSet,
                                                             rankOffset,
                                                             ghostRank,
                                                             dofNumber,
                                                             enthalpy,
                                                             dEnthalpy_dTemperature,
                                                             dEnthalpy_dPressure,
                                                             rhsContributionArrayView,
                                                             localRhs,
                                                             localMatrix] GEOS_HOST_DEVICE ( localIndex const a )
        {
          // we need to filter out ghosts here, because targetSet may contain them
          localIndex const ei = targetSet[a];
          if( ghostRank[ei] >= 0 )
          {
            return;
          }

          // add the value to the mass balance equation
          globalIndex const massRowIndex   = dofNumber[ei] - rankOffset;
          globalIndex const energyRowIndex = massRowIndex + 1;
          real64 const rhsValue = rhsContributionArrayView[a] / sizeScalingFactor; // scale the contribution by the sizeScalingFactor here!
          localRhs[massRowIndex] += rhsValue;
          //add the value to the energey balance equation if the flux is positive (i.e., it's a producer)
          if( rhsContributionArrayView[a] > 0.0 )
          {
            globalIndex const pressureDofIndex    = dofNumber[ei] - rankOffset;
            globalIndex const temperatureDofIndex = pressureDofIndex + 1;

            localRhs[energyRowIndex] += enthalpy[ei][0] * rhsValue;

            globalIndex dofIndices[2]{pressureDofIndex, temperatureDofIndex};
            real64 jacobian[2]{rhsValue * dEnthalpy_dPressure[ei][0], rhsValue * dEnthalpy_dTemperature[ei][0]};

            localMatrix.template addToRow< serialAtomic >( energyRowIndex,
                                                           dofIndices,
                                                           jacobian,
                                                           2 );
          }
        } );
      }
      else
      {
        forAll< parallelDevicePolicy<> >( targetSet.size(), [sizeScalingFactor,
                                                             targetSet,
                                                             rankOffset,
                                                             ghostRank,
                                                             dofNumber,
                                                             rhsContributionArrayView,
                                                             localRhs] GEOS_HOST_DEVICE ( localIndex const a )
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
      }
    } );
  } );
}

void SinglePhaseBase::keepFlowVariablesConstantDuringInitStep( real64 const time,
                                                               real64 const dt,
                                                               DofManager const & dofManager,
                                                               DomainPartition & domain,
                                                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                               arrayView1d< real64 > const & localRhs ) const
{
  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( time, dt );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase const & subRegion )
    {
      globalIndex const rankOffset = dofManager.rankOffset();
      string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();
      arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );

      arrayView1d< real64 const > const pres = subRegion.getField< fields::flow::pressure >();
      arrayView1d< real64 const > const temp = subRegion.getField< fields::flow::temperature >();

      integer const isThermal = m_isThermal;
      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
      {
        if( ghostRank[ei] >= 0 )
        {
          return;
        }

        globalIndex const dofIndex = dofNumber[ei];
        localIndex const localRow = dofIndex - rankOffset;
        real64 rhsValue;

        // 4.1. Apply pressure value to the matrix/rhs
        FieldSpecificationEqual::SpecifyFieldValue( dofIndex,
                                                    rankOffset,
                                                    localMatrix,
                                                    rhsValue,
                                                    pres[ei], // freeze the current pressure value
                                                    pres[ei] );
        localRhs[localRow] = rhsValue;

        // 4.2. Apply temperature value to the matrix/rhs
        if( isThermal )
        {
          FieldSpecificationEqual::SpecifyFieldValue( dofIndex + 1,
                                                      rankOffset,
                                                      localMatrix,
                                                      rhsValue,
                                                      temp[ei], // freeze the current temperature value
                                                      temp[ei] );
          localRhs[localRow + 1] = rhsValue;
        }
      } );
    } );
  } );
}


void SinglePhaseBase::updateState( DomainPartition & domain )
{

// set mass fraction flag on fluid models
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
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

void SinglePhaseBase::resetStateToBeginningOfStep( DomainPartition & domain )
{
  // set mass fraction flag on fluid models
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                                                   auto & subRegion )
    {
      arrayView1d< real64 > const pres = subRegion.template getField< fields::flow::pressure >();
      arrayView1d< real64 const > const pres_n = subRegion.template getField< fields::flow::pressure_n >();
      pres.setValues< parallelDevicePolicy<> >( pres_n );

      if( m_isThermal )
      {
        arrayView1d< real64 > const temp = subRegion.template getField< fields::flow::temperature >();
        arrayView1d< real64 const > const temp_n = subRegion.template getField< fields::flow::temperature_n >();
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

} /* namespace geos */
