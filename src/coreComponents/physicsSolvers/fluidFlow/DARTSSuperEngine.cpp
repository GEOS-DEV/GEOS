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
 * @file DARTSSuperEngine.cpp
 */

#include "DARTSSuperEngine.hpp"

#include "common/DataTypes.hpp"
#include "common/MpiWrapper.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/fluid/MultiFluidBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "dataRepository/Group.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fieldSpecification/SourceFluxBoundaryCondition.hpp"
#include "finiteVolume/BoundaryStencil.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/fluidFlow/DARTSSuperEngineExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/DARTSSuperEngineKernels.hpp"
//#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFVMKernels.hpp"


namespace geosx
{

namespace
{

MultivariableTableFunction const * makeOBLOperatorsTable( string const & OBLOperatorsTableFile,
                                                          FunctionManager & functionManager )
{
  string const tableName = "OBL_operators_table";
  if( functionManager.hasGroup< MultivariableTableFunction >( tableName ) )
  {
    return functionManager.getGroupPointer< MultivariableTableFunction >( tableName );
  }
  else
  {
    MultivariableTableFunction * const table = dynamicCast< MultivariableTableFunction * >( functionManager.createChild( "MultivariableTableFunction", tableName ) );
    table->initializeFunctionFromFile ( OBLOperatorsTableFile );
    return table;
  }
}
}

using namespace dataRepository;
using namespace constitutive;
//using namespace CompositionalMultiphaseFVMKernels;
using namespace DARTSSuperEngineKernels;

DARTSSuperEngine::DARTSSuperEngine( const string & name,
                                    Group * const parent )
  :
  FlowSolverBase( name, parent )
{
  this->getWrapper< array1d< string > >( viewKeyStruct::fluidNamesString() ).
    setInputFlag( InputFlags::FALSE );

  this->registerWrapper( viewKeyStruct::numComponentsString(), &m_numComponents ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Number of components" );

  this->registerWrapper( viewKeyStruct::numPhasesString(), &m_numPhases ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Number of phases" );

  this->registerWrapper( viewKeyStruct::enableEnergyBalanceString(), &m_enableEnergyBalance ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Enable energy balance calculation and temperature degree of freedom" );

  this->registerWrapper( viewKeyStruct::OBLOperatorsTableFileString(), &m_OBLOperatorsTableFile ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "File containing OBL operator values" );

  this->registerWrapper( viewKeyStruct::componentNamesString(), &m_componentNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "List of component names" );

  this->registerWrapper( viewKeyStruct::phaseNamesString(), &m_phaseNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "List of fluid phases" );



  m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::compositionalMultiphaseFVM;
}

void DARTSSuperEngine::initializePreSubGroups()
{
  FlowSolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  if( !fvManager.hasGroup< FluxApproximationBase >( m_discretizationName ) )
  {
    GEOSX_ERROR( "A discretization deriving from FluxApproximationBase must be selected with DARTSSuperEngine" );
  }

}

void DARTSSuperEngine::setupDofs( DomainPartition const & domain,
                                  DofManager & dofManager ) const
{
  dofManager.addField( viewKeyStruct::elemDofFieldString(),
                       DofManager::Location::Elem,
                       m_numDofPerCell,
                       targetRegionNames() );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );
  dofManager.addCoupling( viewKeyStruct::elemDofFieldString(), fluxApprox );
}


void DARTSSuperEngine::implicitStepComplete( real64 const & time,
                                             real64 const & dt,
                                             DomainPartition & domain )
{
  GEOSX_UNUSED_VAR( time );
  GEOSX_UNUSED_VAR( dt );
  integer const numComp = m_numComponents;

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

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

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      pres[ei] += dPres[ei];
      for( localIndex ic = 0; ic < numComp; ++ic )
      {
        compDens[ei][ic] += dCompDens[ei][ic];
      }
    } );

    CoupledSolidBase const & porousMaterial = getConstitutiveModel< CoupledSolidBase >( subRegion, m_solidModelNames[targetIndex] );
    porousMaterial.saveConvergedState();
  } );

  if( m_computeCFLNumbers )
  {
    //computeCFLNumbers( dt, domain );
  }
}

void DARTSSuperEngine::postProcessInput()
{
  // need to override to skip the check for fluidModel, which is enabled in FlowSolverBase
  SolverBase::postProcessInput();
  checkModelNames( m_solidModelNames, viewKeyStruct::solidNamesString() );
  checkModelNames( m_permeabilityModelNames, viewKeyStruct::permeabilityNamesString() );

  GEOSX_ERROR_IF_GT_MSG( m_maxCompFracChange, 1.0,
                         "The maximum absolute change in component fraction must smaller or equal to 1.0" );
  GEOSX_ERROR_IF_LT_MSG( m_maxCompFracChange, 0.0,
                         "The maximum absolute change in component fraction must larger or equal to 0.0" );

  m_OBLOperatorsTable = makeOBLOperatorsTable( m_OBLOperatorsTableFile, FunctionManager::getInstance());

  // Equations: [NC] Molar mass balance, ([1] energy balance if enabled)
  // Primary variables: [1] pressure, [NC-1] global component fractions, ([1] temperature)
  m_numDofPerCell = m_numComponents + m_enableEnergyBalance;

  m_numOBLOperators = COMPUTE_NUM_OPS( m_numPhases, m_numComponents, m_enableEnergyBalance );

  GEOSX_ERROR_IF_NE_MSG( m_numDofPerCell, m_OBLOperatorsTable->getNumDims(),
                         "The number of degrees of freedom per element used in solver and table should match" );

  GEOSX_ERROR_IF_NE_MSG( m_numOBLOperators, m_OBLOperatorsTable->getNumOps(),
                         "The number of operators per element used in solver and table should match" );

}

void DARTSSuperEngine::registerDataOnMesh( Group & meshBodies )
{
  using namespace extrinsicMeshData::flow;

  FlowSolverBase::registerDataOnMesh( meshBodies );

  DomainPartition const & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  ConstitutiveManager const & cm = domain.getConstitutiveManager();



  // 2. Register and resize all fields as necessary
  meshBodies.forSubGroups< MeshBody >( [&]( MeshBody & meshBody )
  {
    MeshLevel & mesh = meshBody.getMeshLevel( 0 );

    forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion )
    {

      subRegion.registerExtrinsicData< pressure >( getName() );
      subRegion.registerExtrinsicData< initialPressure >( getName() );
      subRegion.registerExtrinsicData< deltaPressure >( getName() );
      subRegion.registerExtrinsicData< bcPressure >( getName() );

      subRegion.registerExtrinsicData< temperature >( getName() );

      subRegion.registerExtrinsicData< OBLOperatorValues >( getName() ).
        reference().resizeDimension< 1 >( m_numOBLOperators );
      subRegion.registerExtrinsicData< OBLOperatorValuesOld >( getName() ).
        reference().resizeDimension< 1 >( m_numOBLOperators );
      subRegion.registerExtrinsicData< OBLOperatorDerivatives >( getName() ).
        reference().resizeDimension< 1, 2 >( m_numOBLOperators, m_numDofPerCell );

      // we need to register this fiels in any case (if energy balance is enabled or not)
      // to be able to pass the view to OBLOperatorsKernel
      subRegion.registerExtrinsicData< deltaTemperature >( getName() );

      // The resizing of the arrays needs to happen here, before the call to initializePreSubGroups,
      // to make sure that the dimensions are properly set before the timeHistoryOutput starts its initialization.
      subRegion.registerExtrinsicData< globalCompFraction >( getName() ).
        setDimLabels( 1, m_componentNames ).
        reference().resizeDimension< 1 >( m_numComponents );

      // we need to register this fiels in any case (if there is a single component or not)
      // to be able to pass the view to OBLOperatorsKernel
      subRegion.registerExtrinsicData< deltaGlobalCompFraction >( getName() ).
        reference().resizeDimension< 1 >( m_numComponents );

    } );

    FaceManager & faceManager = mesh.getFaceManager();
    {
      faceManager.registerExtrinsicData< facePressure >( getName() );
    }

  } );
}

real64 DARTSSuperEngine::calculateResidualNorm( DomainPartition const & domain,
                                                DofManager const & dofManager,
                                                arrayView1d< real64 const > const & localRhs )
{
  GEOSX_MARK_FUNCTION;
  GEOSX_UNUSED_VAR( localRhs );

  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  real64 localResidualNorm = 0.0;

  //globalIndex const rankOffset = dofManager.rankOffset();
  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase const & subRegion )
  {
    arrayView1d< globalIndex const > dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
    // arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();
    // arrayView1d< real64 const > const & volume = subRegion.getElementVolume();
    // arrayView1d< real64 const > const & totalDensOld = subRegion.getExtrinsicData< extrinsicMeshData::flow::totalDensityOld >();

    //CoupledSolidBase const & solidModel = getConstitutiveModel< CoupledSolidBase >( subRegion, m_solidModelNames[targetIndex] );

    // arrayView1d< real64 const > const & referencePorosity = solidModel.getReferencePorosity();

    real64 subRegionResidualNorm = 0.0;
    // ResidualNormKernel::launch< parallelDevicePolicy<>,
    //                             parallelDeviceReduce >( localRhs,
    //                                                     rankOffset,
    //                                                     numFluidComponents(),
    //                                                     dofNumber,
    //                                                     elemGhostRank,
    //                                                     referencePorosity,
    //                                                     volume,
    //                                                     totalDensOld,
    //                                                     subRegionResidualNorm );
    localResidualNorm += subRegionResidualNorm;

  } );

  // compute global residual norm
  real64 const residual = std::sqrt( MpiWrapper::sum( localResidualNorm ) );

  if( getLogLevel() >= 1 && logger::internal::rank == 0 )
  {
    std::cout << GEOSX_FMT( "    ( Rfluid ) = ( {:4.2e} ) ;", residual );
  }

  return residual;
}

real64 DARTSSuperEngine::scalingForSystemSolution( DomainPartition const & domain,
                                                   DofManager const & dofManager,
                                                   arrayView1d< real64 const > const & localSolution )
{
  GEOSX_MARK_FUNCTION;

  // check if we want to rescale the Newton update
  if( m_maxCompFracChange >= 1.0 )
  {
    // no rescaling wanted, we just return 1.0;
    return 1.0;
  }

  // real64 constexpr eps = CompositionalMultiphaseBaseKernels::minDensForDivision;
  // real64 const maxCompFracChange = m_maxCompFracChange;

  // localIndex const NC = m_numComponents;

  // MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  // globalIndex const rankOffset = dofManager.rankOffset();
  // string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
  // real64 scalingFactor = 1.0;

  // forTargetSubRegions( mesh, [&]( localIndex const, ElementSubRegionBase const & subRegion )
  // {
  //   arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
  //   arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();

  //   arrayView2d< real64 const, compflow::USD_COMP > const & compDens =
  //     subRegion.getExtrinsicData< extrinsicMeshData::flow::globalCompDensity >();
  //   arrayView2d< real64 const, compflow::USD_COMP > const & dCompDens =
  //     subRegion.getExtrinsicData< extrinsicMeshData::flow::deltaGlobalCompDensity >();

  //   RAJA::ReduceMin< parallelDeviceReduce, real64 > minVal( 1.0 );

  //   forAll< parallelDevicePolicy<> >( dofNumber.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
  //   {
  //     if( elemGhostRank[ei] < 0 )
  //     {
  //       real64 prevTotalDens = 0;
  //       for( localIndex ic = 0; ic < NC; ++ic )
  //       {
  //         prevTotalDens += compDens[ei][ic] + dCompDens[ei][ic];
  //       }

  //       // compute the change in component densities and component fractions
  //       for( localIndex ic = 0; ic < NC; ++ic )
  //       {
  //         localIndex const lid = dofNumber[ei] + ic + 1 - rankOffset;

  //         // compute scaling factor based on relative change in component densities
  //         real64 const absCompDensChange = LvArray::math::abs( localSolution[lid] );
  //         real64 const maxAbsCompDensChange = maxCompFracChange * prevTotalDens;

  //         // This actually checks the change in component fraction, using a lagged total density
  //         // Indeed we can rewrite the following check as:
  //         //    | prevCompDens / prevTotalDens - newCompDens / prevTotalDens | > maxCompFracChange
  //         // Note that the total density in the second term is lagged (i.e, we use prevTotalDens)
  //         // because I found it more robust than using directly newTotalDens (which can vary also
  //         // wildly when the compDens change is large)
  //         if( absCompDensChange > maxAbsCompDensChange && absCompDensChange > eps )
  //         {
  //           minVal.min( maxAbsCompDensChange / absCompDensChange );
  //         }
  //       }
  //     }
  //   } );

  //   if( minVal.get() < scalingFactor )
  //   {
  //     scalingFactor = minVal.get();
  //   }
  // } );

  //return LvArray::math::max( MpiWrapper::min( scalingFactor, MPI_COMM_GEOSX ), m_minScalingFactor );

  // dummy return
  return 1;
}

bool DARTSSuperEngine::checkSystemSolution( DomainPartition const & domain,
                                            DofManager const & dofManager,
                                            arrayView1d< real64 const > const & localSolution,
                                            real64 const scalingFactor )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
  localIndex localCheck = 1;

  forTargetSubRegions( mesh, [&]( localIndex const, ElementSubRegionBase const & subRegion )
  {
    // arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
    // arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();

    // arrayView1d< real64 const > const & pres = subRegion.getReference< array1d< real64 > >( extrinsicMeshData::flow::pressure::key() );
    // arrayView1d< real64 const > const & dPres = subRegion.getReference< array1d< real64 > >(
    // extrinsicMeshData::flow::deltaPressure::key() );
    // arrayView2d< real64 const, compflow::USD_COMP > const & compDens =
    //   subRegion.getExtrinsicData< extrinsicMeshData::flow::globalCompDensity >();
    // arrayView2d< real64 const, compflow::USD_COMP > const & dCompDens =
    //   subRegion.getExtrinsicData< extrinsicMeshData::flow::deltaGlobalCompDensity >();

    localIndex const subRegionSolutionCheck = 0.0;

    // subRegionSolutionCheck = SolutionCheckKernel::launch< parallelDevicePolicy<>,
    //                                parallelDeviceReduce >( localSolution,
    //                                                        dofManager.rankOffset(),
    //                                                        numFluidComponents(),
    //                                                        dofNumber,
    //                                                        elemGhostRank,
    //                                                        pres,
    //                                                        dPres,
    //                                                        compDens,
    //                                                        dCompDens,
    //                                                        m_allowCompDensChopping,
    //                                                        scalingFactor );

    localCheck = std::min( localCheck, subRegionSolutionCheck );
  } );

  return MpiWrapper::min( localCheck, MPI_COMM_GEOSX );
}

void DARTSSuperEngine::applySystemSolution( DofManager const & dofManager,
                                            arrayView1d< real64 const > const & localSolution,
                                            real64 const scalingFactor,
                                            DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  DofManager::CompMask pressureMask( m_numDofPerCell, 0, 1 );
  DofManager::CompMask compFracMask( m_numDofPerCell, 1, m_numComponents );

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::elemDofFieldString(),
                               extrinsicMeshData::flow::deltaPressure::key(),
                               scalingFactor,
                               pressureMask );

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::elemDofFieldString(),
                               extrinsicMeshData::flow::deltaGlobalCompFraction::key(),
                               scalingFactor,
                               compFracMask );
  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( extrinsicMeshData::flow::deltaPressure::key() );
  fieldNames["elems"].emplace_back( extrinsicMeshData::flow::deltaGlobalCompFraction::key() );

  if( m_enableEnergyBalance )
  {
    DofManager::CompMask temperatureMask( m_numDofPerCell, m_numDofPerCell - 1, m_numDofPerCell );
    dofManager.addVectorToField( localSolution,
                                 viewKeyStruct::elemDofFieldString(),
                                 extrinsicMeshData::flow::deltaTemperature::key(),
                                 scalingFactor,
                                 temperatureMask );
    fieldNames["elems"].emplace_back( extrinsicMeshData::flow::deltaTemperature::key() );
  }


  CommunicationTools::getInstance().synchronizeFields( fieldNames, mesh, domain.getNeighbors(), true );
}

void DARTSSuperEngine::initializePostInitialConditionsPreSubGroups()
{
  GEOSX_MARK_FUNCTION;

  FlowSolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( extrinsicMeshData::flow::pressure::key() );
  fieldNames["elems"].emplace_back( extrinsicMeshData::flow::globalCompFraction::key() );

  CommunicationTools::getInstance().synchronizeFields( fieldNames, mesh, domain.getNeighbors(), false );

  // Initialize primary variables from applied initial conditions
  //initializeFluidState( mesh );
}


real64 DARTSSuperEngine::solverStep( real64 const & time_n,
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

void DARTSSuperEngine::backupFields( MeshLevel & mesh ) const
{
  GEOSX_MARK_FUNCTION;

  integer const numComp = m_numComponents;
  integer const numPhase = m_numPhases;

  // backup some fields used in time derivative approximation
  // forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion )
  // {
  //   arrayView1d< integer const > const elemGhostRank = subRegion.ghostRank();

  //   arrayView2d< real64 const, compflow::USD_PHASE > const phaseVolFrac =
  //     subRegion.getExtrinsicData< extrinsicMeshData::flow::phaseVolumeFraction >();

  //   arrayView2d< real64, compflow::USD_PHASE > const phaseVolFracOld =
  //     subRegion.getExtrinsicData< extrinsicMeshData::flow::phaseVolumeFractionOld >();

  //   forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
  //   {
  //     if( elemGhostRank[ei] >= 0 )
  //       return;

  //     for( localIndex ip = 0; ip < numPhase; ++ip )
  //     {
  //       phaseVolFracOld[ei][ip] = phaseVolFrac[ei][ip];

  //       for( localIndex ic = 0; ic < numComp; ++ic )
  //       {
  //         phaseCompFracOld[ei][ip][ic] = phaseCompFrac[ei][0][ip][ic];
  //       }
  //     }
  //   } );
  // } );
}

void
DARTSSuperEngine::implicitStepSetup( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                     real64 const & GEOSX_UNUSED_PARAM( dt ),
                                     DomainPartition & domain )
{
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  // set deltas to zero and recompute dependent quantities
  resetStateToBeginningOfStep( domain );

  // backup fields used in time derivative approximation
  backupFields( mesh );
}

void DARTSSuperEngine::assembleSystem( real64 const GEOSX_UNUSED_PARAM( time_n ),
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

  // assembleFluxTerms( dt,
  //                    domain,
  //                    dofManager,
  //                    localMatrix,
  //                    localRhs );


}

void DARTSSuperEngine::assembleAccumulationAndVolumeBalanceTerms( DomainPartition & domain,
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
    CoupledSolidBase const & solid = getConstitutiveModel< CoupledSolidBase >( subRegion, m_solidModelNames[targetIndex] );

    ElementBasedAssemblyKernelFactory::
      createAndLaunch< parallelDevicePolicy<> >( m_numComponents,
                                                 m_numPhases,
                                                 dofManager.rankOffset(),
                                                 dofKey,
                                                 subRegion,
                                                 solid,
                                                 localMatrix,
                                                 localRhs );

  } );
}

void DARTSSuperEngine::applyBoundaryConditions( real64 const time_n,
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
}

namespace internal
{
string const bcLogMessage = string( "DARTSSUperEngine {}: at time {}s, " )
                            + string( "the <{}> boundary condition '{}' is applied to the element set '{}' in subRegion '{}'. " )
                            + string( "\nThe scale of this boundary condition is {} and multiplies the value of the provided function (if any). " )
                            + string( "\nThe total number of target elements (including ghost elements) is {}. " )
                            + string( "\nNote that if this number is equal to zero for all subRegions, the boundary condition will not be applied on this element set." );
}

void DARTSSuperEngine::applySourceFluxBC( real64 const time,
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


void DARTSSuperEngine::applyDirichletBC( real64 const time,
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
}

void DARTSSuperEngine::solveSystem( DofManager const & dofManager,
                                    ParallelMatrix & matrix,
                                    ParallelVector & rhs,
                                    ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  rhs.scale( -1.0 );
  solution.zero();

  SolverBase::solveSystem( dofManager, matrix, rhs, solution );
}

void DARTSSuperEngine::chopNegativeDensities( DomainPartition & domain )
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

void DARTSSuperEngine::resetStateToBeginningOfStep( DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  forTargetSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( mesh, [&]( localIndex const targetIndex, auto & subRegion )
  {
    arrayView1d< real64 > const & dPres =
      subRegion.template getReference< array1d< real64 > >( extrinsicMeshData::flow::deltaPressure::key() );

    dPres.zero();

    // for single component problems (e.g., geothermal), global component fraction is not a primary variable
    if( m_numComponents > 1 )
    {
      arrayView2d< real64, compflow::USD_COMP > const & dCompFrac =
        subRegion.template getExtrinsicData< extrinsicMeshData::flow::deltaGlobalCompFraction >();
      dCompFrac.zero();
    }

    if( m_enableEnergyBalance )
    {
      arrayView1d< real64 > const & dTemp =
        subRegion.template getReference< array1d< real64 > >( extrinsicMeshData::flow::deltaTemperature::key() );
      dTemp.zero();
    }



    // update porosity and permeability
    updatePorosityAndPermeability( subRegion, targetIndex );

    // update operator values
    updateOBLOperators( subRegion );

  } );
}

void DARTSSuperEngine::updateOBLOperators( ObjectManagerBase & dataGroup ) const
{
  GEOSX_MARK_FUNCTION;

  OBLOperatorsKernelFactory::
    createAndLaunch< parallelDevicePolicy<> >( m_numPhases,
                                               m_numComponents,
                                               m_enableEnergyBalance,
                                               dataGroup,
                                               *m_OBLOperatorsTable );

}


void DARTSSuperEngine::updateState( DomainPartition & domain )
{
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  forTargetSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( mesh, [&]( localIndex const targetIndex, auto & subRegion )
  {
    // update porosity and permeability
    updatePorosityAndPermeability( subRegion, targetIndex );
  } );
}



//START_SPHINX_INCLUDE_01
REGISTER_CATALOG_ENTRY( SolverBase, DARTSSuperEngine, string const &, Group * const )
//END_SPHINX_INCLUDE_01
}// namespace geosx
