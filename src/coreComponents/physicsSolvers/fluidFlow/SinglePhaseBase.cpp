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
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "constitutive/solid/SolidBase.hpp"
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

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace SinglePhaseBaseKernels;

SinglePhaseBase::SinglePhaseBase( const string & name,
                                  Group * const parent ):
  FlowSolverBase( name, parent )
{
  m_numDofPerCell = 1;
}


void SinglePhaseBase::registerDataOnMesh( Group & meshBodies )
{
  using namespace extrinsicMeshData::flow;

  FlowSolverBase::registerDataOnMesh( meshBodies );

  meshBodies.forSubGroups< MeshBody >( [&] ( MeshBody & meshBody )
  {
    MeshLevel & meshLevel = meshBody.getMeshLevel( 0 );

    ElementRegionManager & elemManager = meshLevel.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
    {
      subRegion.registerExtrinsicData< pressure >( getName() );
      subRegion.registerExtrinsicData< initialPressure >( getName() );
      subRegion.registerExtrinsicData< deltaPressure >( getName() );

      subRegion.registerExtrinsicData< deltaVolume >( getName() );

      subRegion.registerExtrinsicData< mobility >( getName() );
      subRegion.registerExtrinsicData< dMobility_dPressure >( getName() );

      subRegion.registerExtrinsicData< densityOld >( getName() );

      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::totalCompressibilityString() ).
              setRestartFlags( RestartFlags::NO_WRITE );

            subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::referencePressureString() ).
              setRestartFlags( RestartFlags::NO_WRITE );

            subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::fluidMassString() ).
            		setPlotLevel( PlotLevel::LEVEL_0 );
    } );

    elemManager.forElementSubRegions< FaceElementSubRegion, EmbeddedSurfaceSubRegion >( [&] ( auto & subRegion )
    {
      subRegion.template registerExtrinsicData< extrinsicMeshData::flow::pressure >( getName() );
      subRegion.template registerExtrinsicData< extrinsicMeshData::flow::initialPressure >( getName() );
      subRegion.template registerExtrinsicData< extrinsicMeshData::flow::deltaPressure >( getName() );

      subRegion.template registerExtrinsicData< extrinsicMeshData::flow::deltaVolume >( getName() );

      subRegion.template registerExtrinsicData< extrinsicMeshData::flow::mobility >( getName() );
      subRegion.template registerExtrinsicData< extrinsicMeshData::flow::dMobility_dPressure >( getName() );

      subRegion.template registerExtrinsicData< extrinsicMeshData::flow::densityOld >( getName() );
    } );

    FaceManager & faceManager = meshLevel.getFaceManager();
    {
      faceManager.registerExtrinsicData< facePressure >( getName() );
    }
  } );
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


void SinglePhaseBase::validateFluidModels( DomainPartition const & domain ) const
{
  for( auto & mesh : domain.getMeshBodies().getSubGroups() )
  {
    MeshLevel const & meshLevel = dynamicCast< MeshBody const * >( mesh.second )->getMeshLevel( 0 );
    validateModelMapping< SingleFluidBase >( meshLevel.getElemManager(), m_fluidModelNames );
  }
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

void SinglePhaseBase::initializePreSubGroups()
{
  FlowSolverBase::initializePreSubGroups();

  validateFluidModels( this->getGroupByPath< DomainPartition >( "/Problem/domain" ) );

  initializeAquiferBC();
}


void SinglePhaseBase::updateEOSExplicit( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                         real64 const & GEOSX_UNUSED_PARAM( dt ),
                                         DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{
	/*GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex,
                                  ElementSubRegionBase & subRegion )
  {
    arrayView1d< real64 const > const vol = subRegion.getReference< array1d< real64 > >( CellBlock::viewKeyStruct::elementVolumeString );
    arrayView1d< real64 const > const poro = subRegion.getReference< array1d< real64 > >( viewKeyStruct::porosityString );
    arrayView1d< real64 const > const mass = subRegion.getReference< array1d< real64 > >( viewKeyStruct::fluidMassString );
    arrayView1d< real64 > const pres = subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );

    CompressibleSinglePhaseFluid & fluid = GetConstitutiveModel< CompressibleSinglePhaseFluid >( subRegion, m_fluidModelNames[targetIndex] );
    arrayView2d< real64 > const dens = fluid.getReference< array2d< real64 > >( SingleFluidBase::viewKeyStruct::densityString );
    real64 pressureCap = 1e8;
    real64 referenceDensity = fluid.referenceDensity();
    real64 referencePressure = fluid.referencePressure();
    real64 compressibility = fluid.compressibility();

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      dens[ei][0] = mass[ei] / ( vol[ei] * poro[ei] );

      // TODO: this formulation of explicit solver may require applying density boundary condition before updating the pressure based on the density
      // allow negative pressure
      pres[ei] =  pres[ei] <= 0.5 * pressureCap ? ( 1 - referenceDensity / dens[ei][0] ) / compressibility + referencePressure
    		  : 0.5 * pressureCap + ( pres[ei] - 0.5 * pressureCap) / ( 1.0 / compressibility - 0.5 * pressureCap ) * 0.5 * pressureCap;
    } );
  } );

  // apply pressure boundary condition in the explicit solver
  FieldSpecificationManager & fsManager = FieldSpecificationManager::get();
  fsManager.Apply( time_n + dt, &domain, "ElementRegions", viewKeyStruct::pressureString,
                    [&]( FieldSpecificationBase const * const fs,
                         string const &,
                         SortedArrayView<localIndex const> const & lset,
                         Group * subRegion,
                         string const & ) -> void
  {
    fs->ApplyFieldValue<FieldSpecificationEqual>( lset,
                                                  time_n + dt,
                                                  subRegion,
                                                  viewKeyStruct::pressureString );
  });

  // update state based on pressure
  forTargetSubRegions( mesh, [&]( localIndex const targetIndex,
                                  ElementSubRegionBase & subRegion )
  {
    updateState( domain );
  } );
  */
}

void SinglePhaseBase::updateFluidModel( ObjectManagerBase & dataGroup, localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  arrayView1d< real64 const > const pres = dataGroup.getExtrinsicData< extrinsicMeshData::flow::pressure >();
  arrayView1d< real64 const > const dPres = dataGroup.getExtrinsicData< extrinsicMeshData::flow::deltaPressure >();

  SingleFluidBase & fluid = getConstitutiveModel< SingleFluidBase >( dataGroup, m_fluidModelNames[targetIndex] );

  switch( m_timeIntegrationOption )
  {
    case TimeIntegrationOption::ExplicitTransient:
    {
      constitutiveUpdatePassThru( fluid, [&]( auto & castedFluid )
      {
        typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();
        FluidUpdateKernel::launch( fluidWrapper, pres );
      } );
      break;
    }
    default:
    {
    	  constitutiveUpdatePassThru( fluid, [&]( auto & castedFluid )
    	  {
    	    typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();
    	    FluidUpdateKernel::launch( fluidWrapper, pres, dPres );
    	  } );
    }
  }
}

void SinglePhaseBase::updateMobility( ObjectManagerBase & dataGroup, localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  // output

  arrayView1d< real64 > const mob =
    dataGroup.getExtrinsicData< extrinsicMeshData::flow::mobility >();

  arrayView1d< real64 > const dMob_dPres =
    dataGroup.getExtrinsicData< extrinsicMeshData::flow::dMobility_dPressure >();

  // input

  ConstitutiveBase & fluid = getConstitutiveModel( dataGroup, m_fluidModelNames[targetIndex] );
  FluidPropViews fluidProps = getFluidProperties( fluid );

  switch( m_timeIntegrationOption )
  {
    case TimeIntegrationOption::ExplicitTransient:
    {
      SinglePhaseBaseKernels::MobilityKernel::launch< parallelDevicePolicy<> >( dataGroup.size(),
                                                                                fluidProps.dens,
                                                                                fluidProps.visc,
                                                                                mob );
      break;
    }
    default:
    {
    	SinglePhaseBaseKernels::MobilityKernel::launch< parallelDevicePolicy<> >( dataGroup.size(),
                                                                            fluidProps.dens,
                                                                            fluidProps.dDens_dPres,
                                                                            fluidProps.visc,
                                                                            fluidProps.dVisc_dPres,
                                                                            mob,
                                                                            dMob_dPres );
    }
  }
}

void SinglePhaseBase::initializePostInitialConditionsPreSubGroups()
{
  GEOSX_MARK_FUNCTION;

  FlowSolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( string( extrinsicMeshData::flow::pressure::key() ) );

  CommunicationTools::getInstance().synchronizeFields( fieldNames, mesh, domain.getNeighbors(), false );

  // Compute hydrostatic equilibrium in the regions for which corresponding field specification tag has been specified
  computeHydrostaticEquilibrium();

  // Moved the following part from ImplicitStepSetup to here since it only needs to be initialized once
  // They will be updated in applySystemSolution and ImplicitStepComplete, respectively
  forTargetSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( mesh, [&]( localIndex const targetIndex,
                                                                                   auto & subRegion )
  {
    ConstitutiveBase const & fluid = getConstitutiveModel( subRegion, m_fluidModelNames[targetIndex] );

    real64 const defaultDensity = getFluidProperties( fluid ).defaultDensity;
    subRegion.template getWrapper< array1d< real64 > >( extrinsicMeshData::flow::densityOld::key() ).setDefaultValue( defaultDensity );

    // 1. update porosity, permeability, and density/viscosity

    updatePorosityAndPermeability( subRegion, targetIndex );
    updateFluidState( subRegion, targetIndex );

    // 2. save the initial density (for use in the single-phase poromechanics solver to compute the deltaBodyForce)

    if( dynamicCast< SingleFluidBase const * >( &fluid ) )
    {
      SingleFluidBase const & singleFluid = dynamicCast< SingleFluidBase const & >( fluid );
      singleFluid.initializeState();
    }

    // 3. save the initial/old porosity

    CoupledSolidBase const & porousSolid = getConstitutiveModel< CoupledSolidBase >( subRegion, m_solidModelNames[targetIndex] );

    porousSolid.initializeState();
  } );

  mesh.getElemManager().forElementRegions< SurfaceElementRegion >( targetRegionNames(),
                                                                   [&]( localIndex const targetIndex,
                                                                        SurfaceElementRegion & region )
  {
    region.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
    {
      ConstitutiveBase & fluid = getConstitutiveModel( subRegion, m_fluidModelNames[targetIndex] );
      real64 const defaultDensity = getFluidProperties( fluid ).defaultDensity;

      subRegion.getWrapper< real64_array >( extrinsicMeshData::flow::hydraulicAperture::key() ).
        setApplyDefaultValue( region.getDefaultAperture() );

      subRegion.getWrapper< real64_array >( FaceElementSubRegion::viewKeyStruct::creationMassString() ).
        setApplyDefaultValue( defaultDensity * region.getDefaultAperture() );
    } );
  } );

  // Save initial pressure field (needed by the poromechanics solvers to compute the deltaPressure needed by the total stress)
  forTargetSubRegions( mesh, [&]( localIndex const,
                                  ElementSubRegionBase & subRegion )
  {
    arrayView1d< real64 const > const pres = subRegion.getExtrinsicData< extrinsicMeshData::flow::pressure >();
    arrayView1d< real64 > const initPres = subRegion.getExtrinsicData< extrinsicMeshData::flow::initialPressure >();

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      initPres[ei] = pres[ei];
    } );
  } );

  backupFields( mesh );
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
  for( string const & regionName : targetRegionNames() )
  {
    regionFilter.insert( regionName );
  }

  // then start the actual table construction
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

    string const & fluidName = m_fluidModelNames[ targetRegionIndex( region.getName() ) ];

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

  real64 dt_return = dt;

  if( m_timeIntegrationOption == TimeIntegrationOption::ExplicitTransient )
  {
	  explicitStepSetup( time_n, dt, domain );
	  dt_return = explicitStep( time_n, dt, cycleNumber, domain );
  }
  else if( m_timeIntegrationOption == TimeIntegrationOption::ImplicitTransient)
  {
	  // setup dof numbers and linear system
	  setupSystem( domain, m_dofManager, m_localMatrix, m_rhs, m_solution );

	  implicitStepSetup( time_n, dt, domain );

	  // currently the only method is implicit time integration
	  dt_return = nonlinearImplicitStep( time_n, dt, cycleNumber, domain );

	  // final step for completion of timestep. typically secondary variable updates and cleanup.
	  implicitStepComplete( time_n, dt_return, domain );
  }
  return dt_return;
}


void SinglePhaseBase::explicitStepSetup( real64 const & time_n ,
                                         real64 const & dt,
                                         DomainPartition & domain )
{
	MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

	// This initialization should be done after running SurfaceGenerator
	static int setFlowSolverTimeStep = 0;
	if( setFlowSolverTimeStep == 0 )
	{
		/*
	  // initialize the fluidMass by defaultDensity (by ron)
		forTargetSubRegions( mesh, [&]( localIndex const targetIndex,
	                                        ElementSubRegionBase & subRegion )
		{
	          arrayView1d< real64 const > const vol = subRegion.getReference< array1d< real64 > >( CellBlock::viewKeyStruct::elementVolumeString );
	          arrayView1d< real64 const > const poro = subRegion.getReference< array1d< real64 > >( viewKeyStruct::porosityString );
	          arrayView1d< real64 > const mass = subRegion.getReference< array1d< real64 > >( viewKeyStruct::fluidMassString );
	          CompressibleSinglePhaseFluid & fluid = GetConstitutiveModel< CompressibleSinglePhaseFluid >( subRegion, m_fluidModelNames[targetIndex] );
	          real64 referenceDensity = fluid.referenceDensity();

	          forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
	          {
	          	mass[ei] =  referenceDensity * vol[ei] * poro[ei] * 1.0008;
	          } );
	    } );
	    updateState(domain);
	    */

	    forTargetSubRegions( mesh, [&]( localIndex const targetIndex,
                                    ElementSubRegionBase & subRegion )
		{

      arrayView1d< real64 > const totalCompressibility = subRegion.getReference< array1d< real64 > >( viewKeyStruct::totalCompressibilityString() );
      arrayView1d< real64 > const referencePressure = subRegion.getReference< array1d< real64 > >( viewKeyStruct::referencePressureString() );

      CoupledSolidBase const & porousSolid = getConstitutiveModel< CoupledSolidBase >( subRegion, m_solidModelNames[targetIndex] );
      SolidBase const & solid = getConstitutiveModel< SolidBase >( subRegion, porousSolid.getsolidModelName());

      //SolidBase & solid = GetConstitutiveModel< SolidBase >( subRegion, m_solidModelNames[targetIndex] );
      CompressibleSinglePhaseFluid & fluid = getConstitutiveModel< CompressibleSinglePhaseFluid >( subRegion, m_fluidModelNames[targetIndex] );

      totalCompressibility.setValues< parallelDevicePolicy<> >( solid.getCompressibility() + fluid.compressibility() );
      referencePressure.setValues< parallelDevicePolicy<> >( fluid.referencePressure() );

      updateEOSExplicit( time_n, dt, domain );
        } );

  }
	updateState(domain);

  // get the maxStableDt for the first time step
  if( setFlowSolverTimeStep == 0 )
  {
    assembleFluxTermsExplicit( time_n, dt, domain );
    setFlowSolverTimeStep = 1;
  }
}

real64 SinglePhaseBase::explicitStep( real64 const & time_n,
                                      real64 const & dt,
                                      integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                                      DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  calculateAndApplyMassFlux( time_n, dt, domain );

  updateEOSExplicit( time_n, dt, domain );

  return dt;
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
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  forTargetSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( mesh, [&]( localIndex const targetIndex,
                                                                                   auto & subRegion )
  {
    arrayView1d< real64 > const & dPres = subRegion.template getExtrinsicData< extrinsicMeshData::flow::deltaPressure >();
    arrayView1d< real64 > const & dVol = subRegion.template getExtrinsicData< extrinsicMeshData::flow::deltaVolume >();

    dPres.zero();
    dVol.zero();

    // This should fix NaN density in newly created fracture elements
    updatePorosityAndPermeability( subRegion, targetIndex );
    updateFluidState( subRegion, targetIndex );
  } );

  forTargetSubRegions< FaceElementSubRegion >( mesh, [&]( localIndex const targetIndex,
                                                          FaceElementSubRegion & subRegion )
  {
    arrayView1d< real64 const > const aper = subRegion.getExtrinsicData< extrinsicMeshData::flow::hydraulicAperture >();
    arrayView1d< real64 > const aper0 = subRegion.getExtrinsicData< extrinsicMeshData::flow::aperture0 >();

    aper0.setValues< parallelDevicePolicy<> >( aper );

    // Needed coz faceElems don't exist when initializing.
    CoupledSolidBase const & porousSolid = getConstitutiveModel< CoupledSolidBase >( subRegion, m_solidModelNames[targetIndex] );
    porousSolid.saveConvergedState();

    updatePorosityAndPermeability( subRegion, targetIndex );
    updateFluidState( subRegion, targetIndex );
  } );

  backupFields( mesh );
}

void SinglePhaseBase::implicitStepComplete( real64 const & time,
                                            real64 const & dt,
                                            DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  // note: we have to save the aquifer state **before** updating the pressure,
  // otherwise the aquifer flux is saved with the wrong pressure time level
  saveAquiferConvergedState( time, dt, domain );

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex,
                                  ElementSubRegionBase & subRegion )
  {
    arrayView1d< real64 const > const dPres = subRegion.getExtrinsicData< extrinsicMeshData::flow::deltaPressure >();
    arrayView1d< real64 const > const dVol = subRegion.getExtrinsicData< extrinsicMeshData::flow::deltaVolume >();

    arrayView1d< real64 > const pres = subRegion.getExtrinsicData< extrinsicMeshData::flow::pressure >();
    arrayView1d< real64 > const vol = subRegion.getReference< array1d< real64 > >( CellElementSubRegion::viewKeyStruct::elementVolumeString() );

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      pres[ei] += dPres[ei];
      vol[ei] += dVol[ei];
    } );

    CoupledSolidBase const & porousSolid = getConstitutiveModel< CoupledSolidBase >( subRegion, m_solidModelNames[targetIndex] );

    porousSolid.saveConvergedState();
  } );



  forTargetSubRegions< FaceElementSubRegion >( mesh, [&]( localIndex const,
                                                          FaceElementSubRegion & subRegion )
  {
    arrayView1d< integer const > const elemGhostRank = subRegion.ghostRank();
    arrayView1d< real64 const > const volume = subRegion.getElementVolume();
    arrayView1d< real64 const > const densOld = subRegion.getExtrinsicData< extrinsicMeshData::flow::densityOld >();
    arrayView1d< real64 > const creationMass = subRegion.getReference< real64_array >( FaceElementSubRegion::viewKeyStruct::creationMassString() );

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      if( elemGhostRank[ei] < 0 )
      {
        if( volume[ei] * densOld[ei] > 1.1 * creationMass[ei] )
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
}


void SinglePhaseBase::assembleSystem( real64 const time_n,
                                      real64 const dt,
                                      DomainPartition & domain,
                                      DofManager const & dofManager,
                                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                      arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  assembleAccumulationTerms< parallelDevicePolicy<> >( domain,
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

template< typename POLICY >
void SinglePhaseBase::accumulationLaunch( localIndex const targetIndex,
                                          CellElementSubRegion const & subRegion,
                                          DofManager const & dofManager,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs )
{
  string const dofKey = dofManager.getKey( extrinsicMeshData::flow::pressure::key() );
  globalIndex const rankOffset = dofManager.rankOffset();
  arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
  arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

  arrayView1d< real64 const > const densityOld =
    subRegion.getExtrinsicData< extrinsicMeshData::flow::densityOld >();

  arrayView1d< real64 const > const volume = subRegion.getElementVolume();

  ConstitutiveBase const & fluid = getConstitutiveModel( subRegion, fluidModelNames()[targetIndex] );
  FluidPropViews const fluidProps = getFluidProperties( fluid );
  arrayView2d< real64 const > const density = fluidProps.dens;
  arrayView2d< real64 const > const dDens_dPres = fluidProps.dDens_dPres;

  CoupledSolidBase const & solidModel = getConstitutiveModel< CoupledSolidBase >( subRegion, m_solidModelNames[targetIndex] );

  arrayView2d< real64 const > const & porosity    = solidModel.getPorosity();
  arrayView2d< real64 const > const & porosityOld = solidModel.getOldPorosity();
  arrayView2d< real64 const > const & dPoro_dPres = solidModel.getDporosity_dPressure();

  AccumulationKernel::template launch< POLICY >( subRegion.size(),
                                                 rankOffset,
                                                 dofNumber,
                                                 ghostRank,
                                                 volume,
                                                 porosityOld,
                                                 porosity,
                                                 dPoro_dPres,
                                                 densityOld,
                                                 density,
                                                 dDens_dPres,
                                                 localMatrix,
                                                 localRhs );

}

template< typename POLICY >
void SinglePhaseBase::accumulationLaunch( localIndex const targetIndex,
                                          SurfaceElementSubRegion const & subRegion,
                                          DofManager const & dofManager,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs )
{
  string const dofKey = dofManager.getKey( extrinsicMeshData::flow::pressure::key() );
  globalIndex const rankOffset = dofManager.rankOffset();
  arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
  arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();

  arrayView1d< real64 const > const & densityOld = subRegion.getExtrinsicData< extrinsicMeshData::flow::densityOld >();
  arrayView1d< real64 const > const & volume = subRegion.getElementVolume();
  arrayView1d< real64 const > const & deltaVolume = subRegion.getExtrinsicData< extrinsicMeshData::flow::deltaVolume >();

  ConstitutiveBase const & fluid = getConstitutiveModel( subRegion, fluidModelNames()[targetIndex] );
  FluidPropViews const fluidProps = getFluidProperties( fluid );
  arrayView2d< real64 const > const & density = fluidProps.dens;
  arrayView2d< real64 const > const & dDens_dPres = fluidProps.dDens_dPres;


#if !defined(ALLOW_CREATION_MASS)
  static_assert( true, "must have ALLOW_CREATION_MASS defined" );
#endif

#if ALLOW_CREATION_MASS
  arrayView1d< real64 const > const &
  creationMass = subRegion.getReference< real64_array >( SurfaceElementSubRegion::viewKeyStruct::creationMassString() );
#endif

  CoupledSolidBase const & solidModel = getConstitutiveModel< CoupledSolidBase >( subRegion, m_solidModelNames[targetIndex] );


  arrayView2d< real64 const > const & porosity    = solidModel.getPorosity();
  arrayView2d< real64 const > const & porosityOld = solidModel.getOldPorosity();
  arrayView2d< real64 const > const & dPoro_dPres = solidModel.getDporosity_dPressure();

  AccumulationKernel::template launch< POLICY >( subRegion.size(),
                                                 rankOffset,
                                                 dofNumber,
                                                 ghostRank,
                                                 volume,
                                                 deltaVolume,
                                                 porosityOld,
                                                 porosity,
                                                 dPoro_dPres,
                                                 densityOld,
                                                 density,
                                                 dDens_dPres,
#if ALLOW_CREATION_MASS
                                                 creationMass,
#endif
                                                 localMatrix,
                                                 localRhs );

}

template< typename POLICY >
void SinglePhaseBase::assembleAccumulationTerms( DomainPartition & domain,
                                                 DofManager const & dofManager,
                                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                 arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  forTargetSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( mesh,
                                                                        [&]( localIndex const targetIndex,
                                                                             auto const & subRegion )
  {
    accumulationLaunch< POLICY >( targetIndex, subRegion, dofManager, localMatrix, localRhs );
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
  string const dofKey = dofManager.getKey( extrinsicMeshData::flow::pressure::key() );

  fsManager.apply( time_n + dt,
                   domain,
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


    arrayView1d< globalIndex const > const dofNumber =
      subRegion.getReference< array1d< globalIndex > >( dofKey );

    arrayView1d< real64 const > const pres =
      subRegion.getReference< array1d< real64 > >( extrinsicMeshData::flow::pressure::key() );

    arrayView1d< real64 const > const dPres =
      subRegion.getReference< array1d< real64 > >( extrinsicMeshData::flow::deltaPressure::key() );

    // call the application of the boundary condition to alter the matrix and rhs
    fs.applyBoundaryConditionToSystem< FieldSpecificationEqual,
                                       parallelDevicePolicy<> >( lset,
                                                                 time_n + dt,
                                                                 subRegion,
                                                                 dofNumber,
                                                                 dofManager.rankOffset(),
                                                                 localMatrix,
                                                                 localRhs,
                                                                 [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      return pres[a] + dPres[a];
    } );
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
  string const dofKey = dofManager.getKey( extrinsicMeshData::flow::pressure::key() );

  fsManager.apply( time_n + dt,
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
                                   getName(), time_n+dt, SourceFluxBoundaryCondition::catalogName(),
                                   fs.getName(), setName, subRegion.getName(), fs.getScale(), numTargetElems ) );
    }

    arrayView1d< globalIndex const > const
    dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );

    fs.applyBoundaryConditionToSystem< FieldSpecificationAdd,
                                       parallelDevicePolicy<> >( targetSet.toViewConst(),
                                                                 time_n + dt,
                                                                 dt,
                                                                 subRegion,
                                                                 dofNumber,
                                                                 dofManager.rankOffset(),
                                                                 localMatrix,
                                                                 localRhs,
                                                                 [] GEOSX_HOST_DEVICE ( localIndex const )
    {
      return 0.0;
    } );

  } );
}

void SinglePhaseBase::updateFluidState( ObjectManagerBase & subRegion, localIndex targetIndex ) const
{
  updateFluidModel( subRegion, targetIndex );
  updateMobility( subRegion, targetIndex );
}

void SinglePhaseBase::updateState( DomainPartition & domain )
{
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  this->template forTargetSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( mesh, [&] ( localIndex const targetIndex,
                                                                                                   auto & subRegion )
  {
    updatePorosityAndPermeability( subRegion, targetIndex );
    updateFluidState( subRegion, targetIndex );
  } );
}

void SinglePhaseBase::solveSystem( DofManager const & dofManager,
                                   ParallelMatrix & matrix,
                                   ParallelVector & rhs,
                                   ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  rhs.scale( -1.0 );
  solution.zero();

  SolverBase::solveSystem( dofManager, matrix, rhs, solution );
}

void SinglePhaseBase::resetStateToBeginningOfStep( DomainPartition & domain )
{
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  forTargetSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( mesh, [&]( localIndex const targetIndex,
                                                                                   auto & subRegion )
  {
    arrayView1d< real64 > const & dPres =
      subRegion.template getExtrinsicData< extrinsicMeshData::flow::deltaPressure >();

    dPres.zero();

    updatePorosityAndPermeability( subRegion, targetIndex );
    updateFluidState( subRegion, targetIndex );
  } );
}

void SinglePhaseBase::backupFields( MeshLevel & mesh ) const
{
  forTargetSubRegions( mesh, [&]( localIndex const targetIndex,
                                  ElementSubRegionBase & subRegion )
  {
    ConstitutiveBase const & fluid = getConstitutiveModel( subRegion, m_fluidModelNames[targetIndex] );
    arrayView2d< real64 const > const & dens = getFluidProperties( fluid ).dens;

    arrayView1d< real64 > const & densOld = subRegion.getExtrinsicData< extrinsicMeshData::flow::densityOld >();

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      densOld[ei] = dens[ei][0];
    } );
  } );
}

} /* namespace geosx */
