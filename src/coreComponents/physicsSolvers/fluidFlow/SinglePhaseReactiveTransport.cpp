/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SinglePhaseReactiveTransport.cpp
 */

#include "SinglePhaseReactiveTransport.hpp"

#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "constitutive/fluid/singlefluid/reactive/ReactiveSingleFluid.hpp"
#include "constitutive/fluid/multifluid/reactive/ReactiveFluidSelector.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mesh/DomainPartition.hpp"

/**
 * @namespace the geos namespace that encapsulates the majority of the code
 */
namespace geos
{

using namespace dataRepository;
using namespace constitutive;

SinglePhaseReactiveTransport::SinglePhaseReactiveTransport( const string & name,
                                                            Group * const parent ):
  SinglePhaseBase( name, parent ),
  m_numPrimarySpecies( 0 )
{
  // To add modeling parameters we want to add here
}

// TODO: we need to update the class of ReactiveSingleFluid to be consistent with the chemistry module!!!
void SinglePhaseReactiveTransport::registerDataOnMesh( Group & meshBodies )
{
  using namespace fields::flow;

  SinglePhaseBase::registerDataOnMesh( meshBodies );

  DomainPartition const & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  ConstitutiveManager const & cm = domain.getConstitutiveManager();

  // 0. Find a reactive fluid model name (at this point, models are already attached to subregions)
  forDiscretizationOnMeshTargets( meshBodies, [&]( string const &,
                                                   MeshLevel & mesh,
                                                   string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      if( m_reactiveFluidModelName.empty() )
      {
        m_reactiveFluidModelName = getConstitutiveName< ReactiveSingleFluid >( subRegion );
      }
    } );
  } );

  // 1. Set key dimensions of the problem
  // Check needed to avoid errors when running in schema generation mode.
  if( !m_reactiveFluidModelName.empty() )
  {
    ReactiveSingleFluid const & reactiveFluid = cm.getConstitutiveRelation< ReactiveSingleFluid >( m_reactiveFluidModelName );
    m_numPrimarySpecies = reactiveFluid.numPrimarySpecies();
    m_isThermal = reactiveFluid.isThermal();
  }

  // n_c components + one pressure ( + one temperature if needed )
  m_numDofPerCell = m_isThermal ? m_numPrimarySpecies + 2 : m_numPrimarySpecies + 1;

  // 2. Register and resize all fields as necessary (to finish)
  forDiscretizationOnMeshTargets( meshBodies, [&]( string const &,
                                                   MeshLevel & mesh,
                                                   string_array const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                              [&]( localIndex const,
                                                                   ElementSubRegionBase & subRegion )
    {
      subRegion.registerField< logPrimarySpeciesConcentration >( getName() ).
        reference().resizeDimension< 1 >( m_numPrimarySpecies );

      subRegion.registerField< logPrimarySpeciesConcentration_n >( getName() ).
        reference().resizeDimension< 1 >( m_numPrimarySpecies );

      subRegion.registerField< totalPrimarySpeciesAmount >( getName() ).
        reference().resizeDimension< 1 >( m_numPrimarySpecies );

      subRegion.registerField< totalPrimarySpeciesAmount_n >( getName() ).
        reference().resizeDimension< 1 >( m_numPrimarySpecies );

      subRegion.registerField< bcLogPrimarySpeciesConcentration >( getName() ).
        reference().resizeDimension< 1 >( m_numPrimarySpecies );
    } );
  } );
}

void SinglePhaseReactiveTransport::setupDofs( DomainPartition const & domain,
                                              DofManager & dofManager ) const
{
  // add a field for the cell-centered degrees of freedom
  dofManager.addField( viewKeyStruct::elemDofFieldString(),
                       FieldLocation::Elem,
                       m_numDofPerCell,
                       getMeshTargets() );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  dofManager.addCoupling( viewKeyStruct::elemDofFieldString(), fluxApprox );
}

void SinglePhaseReactiveTransport::resetStateToBeginningOfStep( DomainPartition & domain )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                                                   auto & subRegion )
    {
      arrayView1d< real64 > const pres = subRegion.template getField< fields::flow::pressure >();
      arrayView1d< real64 const > const pres_n = subRegion.template getField< fields::flow::pressure_n >();
      pres.setValues< parallelDevicePolicy<> >( pres_n );

      arrayView2d< real64, compflow::USD_COMP > const logPrimarySpeciesConc = subRegion.template getField< fields::flow::logPrimarySpeciesConcentration >();
      arrayView2d< real64 const, compflow::USD_COMP > const logPrimarySpeciesConc_n = subRegion.template getField< fields::flow::logPrimarySpeciesConcentration_n >();
      logPrimarySpeciesConc.setValues< parallelDevicePolicy<> >( logPrimarySpeciesConc_n );

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
        updateEnergy( subRegion );
      }
    } );
  } );
}

// void SinglePhaseReactiveTransport::implicitStepComplete( real64 const & time,
//                                                          real64 const & dt,
//                                                          DomainPartition & domain )
// {

// }

void SinglePhaseReactiveTransport::assembleSystem( real64 const GEOS_UNUSED_PARAM( time_n ),
                                                   real64 const dt,
                                                   DomainPartition & domain,
                                                   DofManager const & dofManager,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                   arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  assembleAccumulationTermsInMassBalanceAndSpeciesAmountEqs( dt,
                                                             domain,
                                                             dofManager,
                                                             localMatrix,
                                                             localRhs );
  assembleFluxTerms( dt,
                     domain,
                     dofManager,
                     localMatrix,
                     localRhs );
}

void SinglePhaseReactiveTransport::assembleAccumulationTermsInMassBalanceAndSpeciesAmountEqs( real64 const dt,
                                                                                              DomainPartition & domain,
                                                                                              DofManager const & dofManager,
                                                                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                                              arrayView1d< real64 > const & localRhs ) const
{
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase const & subRegion )
    {
      geos::constitutive::ReactiveSingleFluid const & fluid =
        getConstitutiveModel< geos::constitutive::ReactiveSingleFluid >( subRegion, subRegion.template getReference< string >( viewKeyStruct::fluidNamesString() ) );
      geos::constitutive::CoupledSolidBase const & solid =
        getConstitutiveModel< geos::constitutive::CoupledSolidBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::solidNamesString() ) );

      string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

      if( m_isThermal )
      {
        singlePhaseReactiveBaseKernels::
          AccumulationKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( m_numPrimarySpecies,
                                                     dt,
                                                     dofManager.rankOffset(),
                                                     dofKey,
                                                     subRegion,
                                                     fluid,
                                                     solid,
                                                     localMatrix,
                                                     localRhs );
      }
      else
      {
        singlePhaseReactiveBaseKernels::
          AccumulationKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( m_numPrimarySpecies,
                                                     dt,
                                                     dofManager.rankOffset(),
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

void SinglePhaseReactiveTransport::assembleFluxTerms( real64 const dt,
                                                      DomainPartition const & domain,
                                                      DofManager const & dofManager,
                                                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                      arrayView1d< real64 > const & localRhs )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               string_array const & )
  {
    NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

    string const & dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

    fluxApprox.forAllStencils( mesh, [&] ( auto & stencil )
    {
      typename TYPEOFREF( stencil ) ::KernelWrapper stencilWrapper = stencil.createKernelWrapper();

      if( m_isThermal ) // To implement the thermal case
      {
        singlePhaseReactiveFVMKernels::
          FluxComputeKernelFactory::createAndLaunch< parallelDevicePolicy<> >( m_numPrimarySpecies,
                                                                               dofManager.rankOffset(),
                                                                               dofKey,
                                                                               getName(),
                                                                               mesh.getElemManager(),
                                                                               stencilWrapper,
                                                                               dt,
                                                                               localMatrix.toViewConstSizes(),
                                                                               localRhs.toView() );
      }
      else
      {
        singlePhaseReactiveFVMKernels::
          FluxComputeKernelFactory::createAndLaunch< parallelDevicePolicy<> >( m_numPrimarySpecies,
                                                                               dofManager.rankOffset(),
                                                                               dofKey,
                                                                               getName(),
                                                                               mesh.getElemManager(),
                                                                               stencilWrapper,
                                                                               dt,
                                                                               localMatrix.toViewConstSizes(),
                                                                               localRhs.toView() );
      }

      // To add diffusion
    } );
  } );
}

SinglePhaseBase::FluidPropViews SinglePhaseReactiveTransport::getFluidProperties( constitutive::ConstitutiveBase const & fluid ) const
{
  ReactiveSingleFluid const & reactiveFluid = dynamicCast< ReactiveSingleFluid const & >( fluid );
  return { reactiveFluid.density(),
           reactiveFluid.dDensity(),
           reactiveFluid.viscosity(),
           reactiveFluid.dViscosity(),
           reactiveFluid.defaultDensity(),
           reactiveFluid.defaultViscosity() };
}

void SinglePhaseReactiveTransport::updateState( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  SinglePhaseBase::updateState( domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                                                   auto & subRegion )
    {
      updateSpeciesAmount( subRegion );
    } );
  } );
}

void SinglePhaseReactiveTransport::updateSpeciesAmount( ElementSubRegionBase & subRegion ) const
{
  GEOS_MARK_FUNCTION;

  arrayView2d< real64, compflow::USD_COMP > const totalPrimarySpeciesAmount = subRegion.getField< fields::flow::totalPrimarySpeciesAmount >();
  arrayView2d< real64, compflow::USD_COMP > const totalPrimarySpeciesAmount_n = subRegion.getField< fields::flow::totalPrimarySpeciesAmount_n >();

  CoupledSolidBase const & porousSolid =
    getConstitutiveModel< CoupledSolidBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::solidNamesString() ) );
  arrayView2d< real64 const > const porosity = porousSolid.getPorosity();
  arrayView2d< real64 const > const porosity_n = porousSolid.getPorosity_n();

  arrayView1d< real64 const > const volume = subRegion.getElementVolume();
  arrayView1d< real64 > const deltaVolume = subRegion.getField< fields::flow::deltaVolume >();

  ReactiveSingleFluid & fluid =
    getConstitutiveModel< ReactiveSingleFluid >( subRegion, subRegion.getReference< string >( viewKeyStruct::fluidNamesString() ) );
  arrayView2d< real64 const, compflow::USD_COMP > const primarySpeciesAggregateConcentration = fluid.primarySpeciesAggregateConcentration();
  arrayView2d< real64 const, compflow::USD_COMP > const primarySpeciesAggregateConcentration_n = fluid.primarySpeciesAggregateConcentration_n();

  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    for( integer is = 0; is < m_numPrimarySpecies; ++is )
    {
      totalPrimarySpeciesAmount[ei][is] = porosity[ei][0] * ( volume[ei] + deltaVolume[ei] ) * primarySpeciesAggregateConcentration[ei][is];

      if( isZero( totalPrimarySpeciesAmount_n[ei][is] ) )
        totalPrimarySpeciesAmount_n[ei][is] = porosity_n[ei][0] * volume[ei] * primarySpeciesAggregateConcentration_n[ei][is];
    }
  } );
}

void SinglePhaseReactiveTransport::updateFluidModel( ObjectManagerBase & dataGroup ) const
{
  GEOS_MARK_FUNCTION;

  arrayView1d< real64 const > const pres = dataGroup.getField< fields::flow::pressure >();
  arrayView1d< real64 const > const temp = dataGroup.getField< fields::flow::temperature >();
  arrayView2d< real64 const, compflow::USD_COMP > const logPrimaryConc = dataGroup.getField< fields::flow::logPrimarySpeciesConcentration >();

  ReactiveSingleFluid & fluid =
    getConstitutiveModel< ReactiveSingleFluid >( dataGroup, dataGroup.getReference< string >( viewKeyStruct::fluidNamesString() ) );

  constitutive::constitutiveUpdatePassThru( fluid, [&]( auto & castedFluid )
  {
    typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();
    singlePhaseReactiveBaseKernels::FluidUpdateKernel::launch( fluidWrapper, pres, temp, logPrimaryConc );
  } );
}

void SinglePhaseReactiveTransport::initializeFluidState( MeshLevel & mesh, string_array const & regionNames )
{
  mesh.getElemManager().forElementSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                                                 auto & subRegion )
  {
    ReactiveSingleFluid const & fluid =
      getConstitutiveModel< ReactiveSingleFluid >( subRegion, subRegion.template getReference< string >( viewKeyStruct::fluidNamesString()));
    updateFluidState( subRegion );

    // 2. save the initial density (for use in the single-phase poromechanics solver to compute the deltaBodyForce)
    fluid.initializeState();

    SinglePhaseBase::updateMass( subRegion );
    updateSpeciesAmount( subRegion );
  } );
}

void SinglePhaseReactiveTransport::initializePostInitialConditionsPreSubGroups()
{
  GEOS_MARK_FUNCTION;

  FlowSolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addElementFields( { fields::flow::logPrimarySpeciesConcentration::key() },
                                     regionNames );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), false );
  } );

  FlowSolverBase::initializeState( domain );
}

void SinglePhaseReactiveTransport::applyBoundaryConditions( real64 const time_n,
                                                            real64 const dt,
                                                            DomainPartition & domain,
                                                            DofManager const & dofManager,
                                                            CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                            arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  // if( m_keepVariablesConstantDuringInitStep )
  // {
  //   // this function is going to force the current flow state to be constant during the time step
  //   // this is used when the poromechanics solver is performing the stress initialization
  //   // TODO: in the future, a dedicated poromechanics kernel should eliminate the flow vars to construct a reduced system
  //   //       which will remove the need for this brittle passing aroung of flag
  //   keepVariablesConstantDuringInitStep( time_n, dt, dofManager, domain, localMatrix.toViewConstSizes(), localRhs.toView() );
  // }
  // else
  // {
  // apply pressure boundary conditions.
  applyPresSpeciesDirichletBC( time_n, dt, domain, dofManager, localMatrix.toViewConstSizes(), localRhs.toView() );

  // // apply flux boundary conditions (To finish)
  // applySourceFluxBC( time_n, dt, dofManager, domain, localMatrix.toViewConstSizes(), localRhs.toView() );

  // // apply aquifer boundary conditions (To finish)
  // applyAquiferBC( time_n, dt, dofManager, domain, localMatrix.toViewConstSizes(), localRhs.toView() );
  // }
}

// // To finish
// void SinglePhaseReactiveTransport::applySourceFluxBC( real64 const time,
//                                                       real64 const dt,
//                                                       DofManager const & dofManager,
//                                                       DomainPartition & domain,
//                                                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
//                                                       arrayView1d< real64 > const & localRhs ) const
// {
//   GEOS_MARK_FUNCTION;

//   FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

//   string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

//   // Step 1: count individual source flux boundary conditions

//   std::map< string, localIndex > bcNameToBcId;
//   localIndex bcCounter = 0;

//   fsManager.forSubGroups< SourceFluxBoundaryCondition >( [&] ( SourceFluxBoundaryCondition const & bc )
//   {
//     // collect all the bc names to idx
//     bcNameToBcId[bc.getName()] = bcCounter;
//     bcCounter++;
//   } );

//   if( bcCounter == 0 )
//   {
//     return;
//   }

//   // Step 2: count the set size for each source flux (each source flux may have multiple target sets)

//   array1d< globalIndex > bcAllSetsSize( bcNameToBcId.size() );

//   computeSourceFluxSizeScalingFactor( time_n,
//                                       dt,
//                                       domain,
//                                       bcNameToBcId,
//                                       bcAllSetsSize.toView() );

//   // Step 3: we are ready to impose the boundary condition, normalized by the set size

//   forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
//                                                                MeshLevel & mesh,
//                                                                arrayView1d< string const > const & )
//   {
//     integer const isThermal = m_isThermal;

//     fsManager.apply< ElementSubRegionBase,
//                      SourceFluxBoundaryCondition >( time + dt,
//                                                     mesh,
//                                                     SourceFluxBoundaryCondition::catalogName(),
//                                                     [&, isThermal]( SourceFluxBoundaryCondition const & fs,
//                                                                     string const & setName,
//                                                                     SortedArrayView< localIndex const > const & targetSet,
//                                                                     ElementSubRegionBase & subRegion,
//                                                                     string const & )
//     {
//       if( fs.getLogLevel() >= 1 && m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
//       {
//         globalIndex const numTargetElems = MpiWrapper::sum< globalIndex >( targetSet.size() );
//         GEOS_LOG_RANK_0( GEOS_FMT( bcLogMessage,
//                                    getName(), time+dt, fs.getCatalogName(), fs.getName(),
//                                    setName, subRegion.getName(), fs.getScale(), numTargetElems ) );
//       }

//       if( targetSet.size() == 0 )
//       {
//         return;
//       }
//       if( !subRegion.hasWrapper( dofKey ) )
//       {
//         if( fs.getLogLevel() >= 1 )
//         {
//           GEOS_LOG_RANK( GEOS_FMT( "{}: trying to apply SourceFlux, but its targetSet named '{}' intersects with non-simulated region
// named '{}'.",
//                                    getDataContext(), setName, subRegion.getName() ) );
//         }
//         return;
//       }

//       arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
//       arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

//       // Step 3.1: get the values of the source boundary condition that need to be added to the rhs

//       array1d< globalIndex > dofArray( targetSet.size() );
//       array1d< real64 > rhsContributionArray( targetSet.size() );
//       arrayView1d< real64 > rhsContributionArrayView = rhsContributionArray.toView();
//       localIndex const rankOffset = dofManager.rankOffset();

//       RAJA::ReduceSum< parallelDeviceReduce, real64 > massProd( 0.0 );

//       // note that the dofArray will not be used after this step (simpler to use dofNumber instead)
//       fs.computeRhsContribution< FieldSpecificationAdd,
//                                  parallelDevicePolicy<> >( targetSet.toViewConst(),
//                                                            time + dt,
//                                                            dt,
//                                                            subRegion,
//                                                            dofNumber,
//                                                            rankOffset,
//                                                            localMatrix,
//                                                            dofArray.toView(),
//                                                            rhsContributionArrayView,
//                                                            [] GEOS_HOST_DEVICE ( localIndex const )
//       {
//         return 0.0;
//       } );

//       // Step 3.2: we are ready to add the right-hand side contributions, taking into account our equation layout

//       // get the normalizer
//       real64 const sizeScalingFactor = bcAllSetsSize[bcNameToBcId.at( fs.getName())];

//       if( isThermal )
//       {

//       }
//       else
//       {
//         integer const fluidComponentId = fs.getComponent();
//         integer const numFluidSpecies = m_numPrimarySpecies;
//         forAll< parallelDevicePolicy<> >( targetSet.size(), [sizeScalingFactor,
//                                                              targetSet,
//                                                              rankOffset,
//                                                              ghostRank,
//                                                              fluidComponentId,
//                                                              numFluidSpecies,
//                                                              dofNumber,
//                                                              rhsContributionArrayView,
//                                                              localRhs,
//                                                              massProd] GEOS_HOST_DEVICE ( localIndex const a )
//         {
//           // we need to filter out ghosts here, because targetSet may contain them
//           localIndex const ei = targetSet[a];
//           if( ghostRank[ei] >= 0 )
//           {
//             return;
//           }

//           real64 const rhsValue = rhsContributionArrayView[a] / sizeScalingFactor; // scale the contribution by the sizeScalingFactor
// here!
//           massProd += rhsValue;

//           globalIndex const totalMassBalanceRow   = dofNumber[ei] - rankOffset;
//           globalIndex const speciesMassBalanceRow = dofNumber[ei] - rankOffset + fluidComponentId + 1;
//           localRhs[totalMassBalanceRow] += rhsValue;
//         } );
//       }
//     } );
//   } );
// }

namespace
{
char const bcLogMessage[] =
  "SinglePhaseReactiveTransport {}: at time {}s, "
  "the <{}> boundary condition '{}' is applied to the element set '{}' in subRegion '{}'. "
  "\nThe scale of this boundary condition is {} and multiplies the value of the provided function (if any). "
  "\nThe total number of target elements (including ghost elements) is {}. "
  "\nNote that if this number is equal to zero for all subRegions, the boundary condition will not be applied on this element set.";
}

void SinglePhaseReactiveTransport::applyPresSpeciesDirichletBC( real64 const time_n,
                                                                real64 const dt,
                                                                DomainPartition & domain,
                                                                DofManager const & dofManager,
                                                                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                arrayView1d< real64 > const & localRhs ) const
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & )
  {
    // 1. Apply pressure Dirichlet BCs, store in a separate field
    applyFieldValue< ElementSubRegionBase >( time_n, dt, mesh, bcLogMessage,
                                             fields::flow::pressure::key(), fields::flow::bcPressure::key() );
    // 2. Apply primary species BC (log promary species concentration) and store them for constitutive call
    applyFieldValue< ElementSubRegionBase >( time_n, dt, mesh, bcLogMessage,
                                             fields::flow::logPrimarySpeciesConcentration::key(), fields::flow::bcLogPrimarySpeciesConcentration::key() );
    // 3. Apply temperature Dirichlet BCs, store in a separate field
    if( m_isThermal )
    {
      applyFieldValue< ElementSubRegionBase >( time_n, dt, mesh, bcLogMessage,
                                               fields::flow::temperature::key(), fields::flow::bcTemperature::key() );
    }

    globalIndex const rankOffset = dofManager.rankOffset();
    string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

    // 4. Apply pressure and log primary species concentration to the system
    fsManager.apply< ElementSubRegionBase >( time_n + dt,
                                             mesh,
                                             fields::flow::pressure::key(),
                                             [&] ( FieldSpecificationBase const &,
                                                   string const &,
                                                   SortedArrayView< localIndex const > const & targetSet,
                                                   ElementSubRegionBase & subRegion,
                                                   string const & )
    {
      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();
      arrayView1d< globalIndex const > const dofNumber =
        subRegion.getReference< array1d< globalIndex > >( dofKey );

      // in the isothermal case, we use the reservoir temperature to enforce the boundary condition
      // in the thermal case, the validation function guarantees that temperature has been provided
      arrayView1d< real64 const > const bcPres =
        subRegion.getReference< array1d< real64 > >( fields::flow::bcPressure::key() );
      arrayView2d< real64 const, compflow::USD_COMP > const bcLogPrimaryConc =
        subRegion.getReference< array2d< real64, compflow::LAYOUT_COMP > >( fields::flow::bcLogPrimarySpeciesConcentration::key() );

      arrayView1d< real64 const > const pres =
        subRegion.getReference< array1d< real64 > >( fields::flow::pressure::key() );
      arrayView2d< real64 const, compflow::USD_COMP > const logPrimaryConc =
        subRegion.getReference< array2d< real64, compflow::LAYOUT_COMP > >( fields::flow::logPrimarySpeciesConcentration::key() );

      integer const numPrimarySpecies = m_numPrimarySpecies;
      forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
      {
        localIndex const ei = targetSet[a];
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
                                                    bcPres[ei],
                                                    pres[ei] );
        localRhs[localRow] = rhsValue;

        integer const speciesDofBeginIndex = m_isThermal? 2:1;

        // 4.2. For each component, apply target global density value
        for( integer is = 0; is < numPrimarySpecies; ++is )
        {
          FieldSpecificationEqual::SpecifyFieldValue( dofIndex + is + speciesDofBeginIndex,
                                                      rankOffset,
                                                      localMatrix,
                                                      rhsValue,
                                                      bcLogPrimaryConc[ei][is],
                                                      logPrimaryConc[ei][is] );
          localRhs[localRow + is + speciesDofBeginIndex] = rhsValue;
        }
      } );
    } );

    // 5. Apply temperature to the system
    if( m_isThermal )
    {
      fsManager.apply< ElementSubRegionBase >( time_n + dt,
                                               mesh,
                                               fields::flow::temperature::key(),
                                               [&] ( FieldSpecificationBase const &,
                                                     string const &,
                                                     SortedArrayView< localIndex const > const & targetSet,
                                                     ElementSubRegionBase & subRegion,
                                                     string const & )
      {
        arrayView1d< integer const > const ghostRank =
          subRegion.getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString() );
        arrayView1d< globalIndex const > const dofNumber =
          subRegion.getReference< array1d< globalIndex > >( dofKey );
        arrayView1d< real64 const > const bcTemp =
          subRegion.getReference< array1d< real64 > >( fields::flow::bcTemperature::key() );
        arrayView1d< real64 const > const temp =
          subRegion.getReference< array1d< real64 > >( fields::flow::temperature::key() );

        forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
        {
          localIndex const ei = targetSet[a];
          if( ghostRank[ei] >= 0 )
          {
            return;
          }

          globalIndex const dofIndex = dofNumber[ei];
          localIndex const localRow = dofIndex - rankOffset;
          real64 rhsValue;

          // 4.2. Apply temperature value to the matrix/rhs
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

real64 SinglePhaseReactiveTransport::calculateResidualNorm( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                                            real64 const & GEOS_UNUSED_PARAM( dt ),
                                                            DomainPartition const & domain,
                                                            DofManager const & dofManager,
                                                            arrayView1d< real64 const > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  integer constexpr numNorm = 3; // total mass balance, energy balance, and species mass balance
  array1d< real64 > localResidualNorm;
  array1d< real64 > localResidualNormalizer;
  localResidualNorm.resize( numNorm );
  localResidualNormalizer.resize( numNorm );

  physicsSolverBaseKernels::NormType const normType = SinglePhaseBase::getNonlinearSolverParameters().normType();

  globalIndex const rankOffset = dofManager.rankOffset();
  string const dofKey = dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );

  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel const & mesh,
                                                                      string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase const & subRegion )
    {
      real64 subRegionResidualNorm[numNorm]{};
      real64 subRegionResidualNormalizer[numNorm]{};

      // step 1: compute the norm in the subRegion

      if( m_isThermal )
      {
        singlePhaseReactiveBaseKernels::
          ResidualNormKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( normType,
                                                     m_numPrimarySpecies,
                                                     rankOffset,
                                                     dofKey,
                                                     localRhs,
                                                     subRegion,
                                                     m_nonlinearSolverParameters.m_minNormalizer,
                                                     subRegionResidualNorm,
                                                     subRegionResidualNormalizer );
      }
      else
      {
        real64 subRegionFlowResidualNorm[2]{};
        real64 subRegionFlowResidualNormalizer[2]{};

        singlePhaseReactiveBaseKernels::
          ResidualNormKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( normType,
                                                     m_numPrimarySpecies,
                                                     rankOffset,
                                                     dofKey,
                                                     localRhs,
                                                     subRegion,
                                                     m_nonlinearSolverParameters.m_minNormalizer,
                                                     subRegionFlowResidualNorm,
                                                     subRegionFlowResidualNormalizer );
        subRegionResidualNorm[0] = subRegionFlowResidualNorm[0];
        subRegionResidualNorm[1] = subRegionFlowResidualNorm[1];
        subRegionResidualNormalizer[0] = subRegionFlowResidualNormalizer[0];
        subRegionResidualNormalizer[1] = subRegionFlowResidualNormalizer[1];
      }

      // step 2: first reduction across meshBodies/regions/subRegions

      if( normType == physicsSolverBaseKernels::NormType::Linf )
      {
        physicsSolverBaseKernels::LinfResidualNormHelper::
          updateLocalNorm< numNorm >( subRegionResidualNorm, localResidualNorm );
      }
      else
      {
        physicsSolverBaseKernels::L2ResidualNormHelper::
          updateLocalNorm< numNorm >( subRegionResidualNorm, subRegionResidualNormalizer, localResidualNorm, localResidualNormalizer );
      }
    } );
  } );

  // step 3: second reduction across MPI ranks

  real64 residualNorm = 0.0;
  array1d< real64 > globalResidualNorm;
  globalResidualNorm.resize( numNorm );
  if( m_isThermal )
  {
    if( normType == physicsSolverBaseKernels::NormType::Linf )
    {
      physicsSolverBaseKernels::LinfResidualNormHelper::
        computeGlobalNorm( localResidualNorm, globalResidualNorm );
    }
    else
    {
      physicsSolverBaseKernels::L2ResidualNormHelper::
        computeGlobalNorm( localResidualNorm, localResidualNormalizer, globalResidualNorm );
    }
    residualNorm = sqrt( globalResidualNorm[0] * globalResidualNorm[0] + globalResidualNorm[1] * globalResidualNorm[1]  + globalResidualNorm[2] * globalResidualNorm[2] );

    GEOS_LOG_LEVEL_INFO_RANK_0_NLR( logInfo::Convergence, GEOS_FMT( "        ( RtotalMass RspeciesAmount ) = ( {:4.2e} {:4.2e} )        ( Renergy ) = ( {:4.2e} )",
                                                                    globalResidualNorm[0], globalResidualNorm[2], globalResidualNorm[1] ));
  }
  else
  {
    if( normType == physicsSolverBaseKernels::NormType::Linf )
    {
      physicsSolverBaseKernels::LinfResidualNormHelper::
        computeGlobalNorm( localResidualNorm, globalResidualNorm );
    }
    else
    {
      physicsSolverBaseKernels::L2ResidualNormHelper::
        computeGlobalNorm( localResidualNorm, localResidualNormalizer, globalResidualNorm );
    }
    residualNorm = sqrt( globalResidualNorm[0] * globalResidualNorm[0] + globalResidualNorm[1] * globalResidualNorm[1] );

    GEOS_LOG_LEVEL_INFO_RANK_0_NLR( logInfo::Convergence, GEOS_FMT( "        ( RtotalMass RspeciesAmount ) = ( {:4.2e} {:4.2e} )",
                                                                    globalResidualNorm[0], globalResidualNorm[1] ) );
  }
  return residualNorm;
}

void SinglePhaseReactiveTransport::applySystemSolution( DofManager const & dofManager,
                                                        arrayView1d< real64 const > const & localSolution,
                                                        real64 const scalingFactor,
                                                        real64 const dt,
                                                        DomainPartition & domain )
{
  GEOS_UNUSED_VAR( dt );

  if( m_isThermal )
  {
    DofManager::CompMask pressureMask( m_numDofPerCell, 0, 1 );
    DofManager::CompMask temperatureMask( m_numDofPerCell, 1, 2 );
    DofManager::CompMask speciesMask( m_numDofPerCell, 2, m_numPrimarySpecies+2 );

    dofManager.addVectorToField( localSolution,
                                 viewKeyStruct::elemDofFieldString(),
                                 fields::flow::pressure::key(),
                                 scalingFactor,
                                 pressureMask );

    dofManager.addVectorToField( localSolution,
                                 viewKeyStruct::elemDofFieldString(),
                                 fields::flow::temperature::key(),
                                 scalingFactor,
                                 temperatureMask );

    dofManager.addVectorToField( localSolution,
                                 viewKeyStruct::elemDofFieldString(),
                                 fields::flow::logPrimarySpeciesConcentration::key(),
                                 scalingFactor,
                                 speciesMask );
  }
  else
  {
    DofManager::CompMask pressureMask( m_numDofPerCell, 0, 1 );
    DofManager::CompMask speciesMask( m_numDofPerCell, 1, m_numPrimarySpecies+1 );

    dofManager.addVectorToField( localSolution,
                                 viewKeyStruct::elemDofFieldString(),
                                 fields::flow::pressure::key(),
                                 scalingFactor,
                                 pressureMask );

    dofManager.addVectorToField( localSolution,
                                 viewKeyStruct::elemDofFieldString(),
                                 fields::flow::logPrimarySpeciesConcentration::key(),
                                 scalingFactor,
                                 speciesMask );
  }

  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel & mesh,
                                                                      string_array const & regionNames )
  {
    std::vector< string > fields{ fields::flow::pressure::key() };

    if( m_isThermal )
    {
      fields.emplace_back( fields::flow::temperature::key() );
    }

    fields.emplace_back( fields::flow::logPrimarySpeciesConcentration::key() );

    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addElementFields( fields, regionNames );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), true );
  } );
}

void SinglePhaseReactiveTransport::saveConvergedState( ElementSubRegionBase & subRegion ) const
{
  FlowSolverBase::saveConvergedState( subRegion );

  arrayView2d< real64 const, compflow::USD_COMP > const totalPrimarySpeciesAmount = subRegion.template getField< fields::flow::totalPrimarySpeciesAmount >();
  arrayView2d< real64, compflow::USD_COMP > const totalPrimarySpeciesAmount_n = subRegion.template getField< fields::flow::totalPrimarySpeciesAmount_n >();
  totalPrimarySpeciesAmount_n.setValues< parallelDevicePolicy<> >( totalPrimarySpeciesAmount );
}

void SinglePhaseReactiveTransport::assembleEDFMFluxTerms( real64 const GEOS_UNUSED_PARAM( time_n ),
                                                          real64 const GEOS_UNUSED_PARAM( dt ),
                                                          DomainPartition const & GEOS_UNUSED_PARAM( domain ),
                                                          DofManager const & GEOS_UNUSED_PARAM( dofManager ),
                                                          CRSMatrixView< real64, globalIndex const > const & GEOS_UNUSED_PARAM( localMatrix ),
                                                          arrayView1d< real64 > const & GEOS_UNUSED_PARAM( localRhs ),
                                                          string const & GEOS_UNUSED_PARAM( jumpDofKey ) )
{}

void SinglePhaseReactiveTransport::applyAquiferBC( real64 const GEOS_UNUSED_PARAM( time ),
                                                   real64 const GEOS_UNUSED_PARAM( dt ),
                                                   DomainPartition & GEOS_UNUSED_PARAM( domain ),
                                                   DofManager const & GEOS_UNUSED_PARAM( dofManager ),
                                                   CRSMatrixView< real64, globalIndex const > const & GEOS_UNUSED_PARAM( localMatrix ),
                                                   arrayView1d< real64 > const & GEOS_UNUSED_PARAM( localRhs ) ) const
{}

void SinglePhaseReactiveTransport::assembleStabilizedFluxTerms( real64 const GEOS_UNUSED_PARAM( dt ),
                                                                DomainPartition const & GEOS_UNUSED_PARAM( domain ),
                                                                DofManager const & GEOS_UNUSED_PARAM( dofManager ),
                                                                CRSMatrixView< real64, globalIndex const > const & GEOS_UNUSED_PARAM( localMatrix ),
                                                                arrayView1d< real64 > const & GEOS_UNUSED_PARAM( localRhs ) )
{}

REGISTER_CATALOG_ENTRY( PhysicsSolverBase, SinglePhaseReactiveTransport, string const &, Group * const )

} /* namespace geos */
