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
 * @file ReactiveCompositionalMultiphaseOBL.cpp
 */

#include "ReactiveCompositionalMultiphaseOBL.hpp"

#include "constitutive/solid/CoupledSolidBase.hpp"
#include "dataRepository/Group.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fieldSpecification/SourceFluxBoundaryCondition.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/fluidFlow/ReactiveCompositionalMultiphaseOBLExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/ReactiveCompositionalMultiphaseOBLKernels.hpp"


namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace ReactiveCompositionalMultiphaseOBLKernels;

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

ReactiveCompositionalMultiphaseOBL::ReactiveCompositionalMultiphaseOBL( const string & name,
                                                                        Group * const parent )
  :
  FlowSolverBase( name, parent ),
  m_numPhases( 0 ),
  m_numComponents( 0 ),
  m_maxCompFracChange( 1.0 ),
  m_minScalingFactor( 0.01 ),
  m_allowOBLChopping( 1 ),
  m_systemSetupDone( false )
{
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
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "File containing OBL operator values" );

  this->registerWrapper( viewKeyStruct::maxCompFracChangeString(), &m_maxCompFracChange ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Maximum (absolute) change in a component fraction between two Newton iterations" );

  this->registerWrapper( viewKeyStruct::transMultExpString(), &m_transMultExp ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Exponent of dynamic transmissibility multiplier" );

  this->registerWrapper( viewKeyStruct::allowLocalOBLChoppingString(), &m_allowOBLChopping ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Allow keeping solution within OBL limits" );

  this->registerWrapper( viewKeyStruct::useDARTSL2NormString(), &m_useDARTSL2Norm ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Use L2 norm calculation similar to one used DARTS" );

  this->registerWrapper( viewKeyStruct::componentNamesString(), &m_componentNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "List of component names" );

  this->registerWrapper( viewKeyStruct::phaseNamesString(), &m_phaseNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "List of fluid phases" );

  m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::compositionalMultiphaseFVM;
}

void ReactiveCompositionalMultiphaseOBL::initializePreSubGroups()
{
  FlowSolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  if( !fvManager.hasGroup< FluxApproximationBase >( m_discretizationName ) )
  {
    GEOSX_ERROR( "A discretization deriving from FluxApproximationBase must be selected with ReactiveCompositionalMultiphaseOBL" );
  }

}

void ReactiveCompositionalMultiphaseOBL::setupDofs( DomainPartition const & domain,
                                                    DofManager & dofManager ) const
{
  dofManager.addField( viewKeyStruct::elemDofFieldString(),
                       FieldLocation::Elem,
                       m_numDofPerCell,
                       getMeshTargets() );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );
  dofManager.addCoupling( viewKeyStruct::elemDofFieldString(), fluxApprox );
}


void ReactiveCompositionalMultiphaseOBL::implicitStepComplete( real64 const & time,
                                                               real64 const & dt,
                                                               DomainPartition & domain )
{
  GEOSX_UNUSED_VAR( time );
  GEOSX_UNUSED_VAR( dt );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      // for output purposes (visualization, etc) we update the last component fraction
      integer const numComp = m_numComponents;
      arrayView2d< real64, compflow::USD_COMP > const compFrac =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::globalCompFraction >();

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
      {
        compFrac[ei][numComp-1] = 1.0;
        for( integer ic = 0; ic < numComp - 1; ++ic )
        {
          compFrac[ei][numComp-1] -= compFrac[ei][ic];
        }
      } );

      // save converged porosity state
      string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );
      CoupledSolidBase const & porousMaterial = getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );
      porousMaterial.saveConvergedState();
    } );
  } );
}

void ReactiveCompositionalMultiphaseOBL::postProcessInput()
{
  // need to override to skip the check for fluidModel, which is enabled in FlowSolverBase
  SolverBase::postProcessInput();

  GEOSX_THROW_IF_GT_MSG( m_maxCompFracChange, 1.0,
                         GEOSX_FMT( "The maximum absolute change in component fraction is set to {}, while it must not be greater than 1.0", m_maxCompFracChange ),
                         InputError );
  GEOSX_THROW_IF_LT_MSG( m_maxCompFracChange, 0.0,
                         GEOSX_FMT( "The maximum absolute change in component fraction is set to {}, while it must not be lesser than 0.0", m_maxCompFracChange ),
                         InputError );

  m_OBLOperatorsTable = makeOBLOperatorsTable( m_OBLOperatorsTableFile, FunctionManager::getInstance());

  // Equations: [NC] Molar mass balance, ([1] energy balance if enabled)
  // Primary variables: [1] pressure, [NC-1] global component fractions, ([1] temperature)
  m_numDofPerCell = m_numComponents + m_enableEnergyBalance;

  m_numOBLOperators = COMPUTE_NUM_OPS( m_numPhases, m_numComponents, m_enableEnergyBalance );

  GEOSX_THROW_IF_NE_MSG( m_numDofPerCell, m_OBLOperatorsTable->numDims(),
                         GEOSX_FMT( "The number of degrees of freedom per element used in solver - {} - and in operator table - {} - should match", m_numDofPerCell, m_OBLOperatorsTable->numDims()),
                         InputError );

  GEOSX_THROW_IF_NE_MSG( m_numOBLOperators, m_OBLOperatorsTable->numOps(),
                         GEOSX_FMT( "The number of operators per element used in solver - {} - and in operator table - {} - should match", m_numOBLOperators, m_OBLOperatorsTable->numOps()),
                         InputError );

}

void ReactiveCompositionalMultiphaseOBL::registerDataOnMesh( Group & meshBodies )
{
  using namespace extrinsicMeshData::flow;
  // 1. Call base class method
  FlowSolverBase::registerDataOnMesh( meshBodies );

  // 2. Register and resize all fields as necessary
  forDiscretizationOnMeshTargets( meshBodies, [&]( string const &,
                                                   MeshLevel & mesh,
                                                   arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      string const solverName = getName();

      subRegion.registerExtrinsicData< pressure >( solverName );
      subRegion.registerExtrinsicData< initialPressure >( solverName );
      subRegion.registerExtrinsicData< pressure_n >( solverName );
      subRegion.registerExtrinsicData< bcPressure >( solverName );

      subRegion.registerExtrinsicData< temperature >( solverName );
      subRegion.registerExtrinsicData< bcTemperature >( solverName );

      subRegion.registerExtrinsicData< OBLOperatorValues >( solverName ).
        reference().resizeDimension< 1 >( m_numOBLOperators );
      subRegion.registerExtrinsicData< OBLOperatorValues_n >( solverName ).
        reference().resizeDimension< 1 >( m_numOBLOperators );
      subRegion.registerExtrinsicData< OBLOperatorDerivatives >( solverName ).
        reference().resizeDimension< 1, 2 >( m_numOBLOperators, m_numDofPerCell );

      // we need to register this fiels in any case (if energy balance is enabled or not)
      // to be able to pass the view to OBLOperatorsKernel
      subRegion.registerExtrinsicData< temperature_n >( solverName );

      // The resizing of the arrays needs to happen here, before the call to initializePreSubGroups,
      // to make sure that the dimensions are properly set before the timeHistoryOutput starts its initialization.
      subRegion.registerExtrinsicData< globalCompFraction >( solverName ).
        setDimLabels( 1, m_componentNames ).
        reference().resizeDimension< 1 >( m_numComponents );

      subRegion.registerExtrinsicData< bcGlobalCompFraction >( solverName ).
        reference().resizeDimension< 1 >( m_numComponents );

      // we need to register this fiels in any case (if there is a single component or not)
      // to be able to pass the view to OBLOperatorsKernel
      subRegion.registerExtrinsicData< globalCompFraction_n >( solverName ).
        reference().resizeDimension< 1 >( m_numComponents );
      // in principle, referencePorosity could be used directly from solid model,
      // but was duplicated to remove dependency on solid
      subRegion.registerExtrinsicData< referencePorosity >( solverName );

      // referencePoreVolume and referenceRockVolume are introduced for the sake of performance:
      // this way the multiplication of constant arrays (e.g., referencePorosity and volume) every Newton step is avoided
      subRegion.registerExtrinsicData< referencePoreVolume >( solverName );
      subRegion.registerExtrinsicData< referenceRockVolume >( solverName );

      // thermal rock properties (again, register in any case)
      // it is not possible to use specificHeatCapacity from solid model here, because specificHeatCapacity includes several quantities,
      // which are split in OBL framework: constant rock volume and,
      // hidden inside operator - therefore variable - rock compressibility and rock energy
      subRegion.registerExtrinsicData< rockVolumetricHeatCapacity >( solverName );
      subRegion.registerExtrinsicData< rockThermalConductivity >( solverName );
      subRegion.registerExtrinsicData< rockKineticRateFactor >( solverName );

    } );

    FaceManager & faceManager = mesh.getFaceManager();
    {
      faceManager.registerExtrinsicData< facePressure >( getName() );
    }

  } );
}

real64 ReactiveCompositionalMultiphaseOBL::calculateResidualNorm( real64 const & GEOSX_UNUSED_PARAM( time ),
                                                                  real64 const & GEOSX_UNUSED_PARAM( dt ),
                                                                  DomainPartition const & domain,
                                                                  DofManager const & dofManager,
                                                                  arrayView1d< real64 const > const & localRhs )
{
  GEOSX_MARK_FUNCTION;
  GEOSX_UNUSED_VAR( localRhs );

  real64 localResidualNorm = 0.0;

  globalIndex const rankOffset = dofManager.rankOffset();
  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase const & subRegion )
    {

      arrayView1d< globalIndex const > dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();
      arrayView1d< real64 const > const & refPoreVolume = subRegion.getExtrinsicData< extrinsicMeshData::flow::referencePoreVolume >();

      real64 subRegionResidualNorm = 0.0;

      if( m_useDARTSL2Norm )
      {
        arrayView2d< real64 const, compflow::USD_OBL_VAL > const OBLVals = subRegion.getExtrinsicData< extrinsicMeshData::flow::OBLOperatorValues >();

        ResidualDARTSL2NormKernel::launch< parallelDevicePolicy<> >( localRhs,
                                                                     rankOffset,
                                                                     m_numDofPerCell,
                                                                     dofNumber,
                                                                     elemGhostRank,
                                                                     refPoreVolume,
                                                                     OBLVals,
                                                                     subRegionResidualNorm );
        if( localResidualNorm < subRegionResidualNorm )
        {
          localResidualNorm = subRegionResidualNorm;
        }
      }
      else
      {
        arrayView2d< real64 const, compflow::USD_OBL_VAL > const OBLVals_n = subRegion.getExtrinsicData< extrinsicMeshData::flow::OBLOperatorValues_n >();

        ResidualNormKernel::launch< parallelDevicePolicy<> >( localRhs,
                                                              rankOffset,
                                                              m_numDofPerCell,
                                                              dofNumber,
                                                              elemGhostRank,
                                                              refPoreVolume,
                                                              OBLVals_n,
                                                              subRegionResidualNorm );
        localResidualNorm += subRegionResidualNorm;
      }
    } );
  } );

  real64 const residual = m_useDARTSL2Norm ? MpiWrapper::max( localResidualNorm ) : std::sqrt( MpiWrapper::sum( localResidualNorm ) );

  GEOSX_LOG_LEVEL_RANK_0( 1, GEOSX_FMT( "    ( Rflow ) = ( {:4.2e} ) ;", residual ) );

  return residual;
}


real64 ReactiveCompositionalMultiphaseOBL::scalingForSystemSolution( DomainPartition const & domain,
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

  real64 constexpr eps = ReactiveCompositionalMultiphaseOBLKernels::minValueForDivision;
  real64 const maxCompFracChange = m_maxCompFracChange;

  localIndex const NC = m_numComponents;

  globalIndex const rankOffset = dofManager.rankOffset();
  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
  real64 scalingFactor = 1.0;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase const & subRegion )
    {
      arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();

      RAJA::ReduceMin< parallelDeviceReduce, real64 > minVal( 1.0 );

      forAll< parallelDevicePolicy<> >( dofNumber.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
      {
        if( elemGhostRank[ei] < 0 )
        {
          real64 lastCompFracChange = 0.0;
          // process NC-1 components first (they have explicit solution)
          for( integer ic = 0; ic < NC - 1; ++ic )
          {
            localIndex const lid = dofNumber[ei] + ic + 1 - rankOffset;

            // compute scaling factor based on relative change in component densities
            real64 const CompFracChange = LvArray::math::abs( localSolution[lid] );
            lastCompFracChange -= localSolution[lid];
            if( CompFracChange > maxCompFracChange && CompFracChange > eps )
            {
              minVal.min( maxCompFracChange / CompFracChange );
            }
          }
          lastCompFracChange = LvArray::math::abs( lastCompFracChange );
          // now deal with the last component
          if( lastCompFracChange > maxCompFracChange && lastCompFracChange > eps )
          {
            minVal.min( maxCompFracChange / lastCompFracChange );
          }

        }
      } );

      if( minVal.get() < scalingFactor )
      {
        scalingFactor = minVal.get();
      }
    } );
  } );

  return LvArray::math::max( MpiWrapper::min( scalingFactor, MPI_COMM_GEOSX ), m_minScalingFactor );
}

bool ReactiveCompositionalMultiphaseOBL::checkSystemSolution( DomainPartition const & domain,
                                                              DofManager const & dofManager,
                                                              arrayView1d< real64 const > const & localSolution,
                                                              real64 const scalingFactor )
{
  GEOSX_MARK_FUNCTION;

  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
  localIndex localCheck = 1;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase const & subRegion )
    {
      arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();

      arrayView1d< real64 const > const & pres =
        subRegion.getReference< array1d< real64 > >( extrinsicMeshData::flow::pressure::key() );
      arrayView2d< real64 const, compflow::USD_COMP > const & compFrac =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::globalCompFraction >();
      arrayView1d< real64 const > const & temp =
        subRegion.getReference< array1d< real64 > >( extrinsicMeshData::flow::temperature::key() );

      localIndex const subRegionSolutionCheck =
        SolutionCheckKernel::launch< parallelDevicePolicy<> >( localSolution,
                                                               dofManager.rankOffset(),
                                                               m_numComponents,
                                                               m_enableEnergyBalance,
                                                               dofNumber,
                                                               elemGhostRank,
                                                               pres,
                                                               compFrac,
                                                               temp,
                                                               m_allowOBLChopping,
                                                               scalingFactor );

      localCheck = std::min( localCheck, subRegionSolutionCheck );
    } );
  } );

  return MpiWrapper::min( localCheck, MPI_COMM_GEOSX );
}

void ReactiveCompositionalMultiphaseOBL::applySystemSolution( DofManager const & dofManager,
                                                              arrayView1d< real64 const > const & localSolution,
                                                              real64 const scalingFactor,
                                                              DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  DofManager::CompMask pressureMask( m_numDofPerCell, 0, 1 );

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::elemDofFieldString(),
                               extrinsicMeshData::flow::pressure::key(),
                               scalingFactor,
                               pressureMask );

  if( m_numComponents > 1 )
  {
    DofManager::CompMask compFracMask( m_numDofPerCell, 1, m_numComponents );

    dofManager.addVectorToField( localSolution,
                                 viewKeyStruct::elemDofFieldString(),
                                 extrinsicMeshData::flow::globalCompFraction::key(),
                                 scalingFactor,
                                 compFracMask );
  }

  if( m_enableEnergyBalance )
  {
    DofManager::CompMask temperatureMask( m_numDofPerCell, m_numDofPerCell - 1, m_numDofPerCell );
    dofManager.addVectorToField( localSolution,
                                 viewKeyStruct::elemDofFieldString(),
                                 extrinsicMeshData::flow::temperature::key(),
                                 scalingFactor,
                                 temperatureMask );
  }

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addElementFields( {  extrinsicMeshData::flow::pressure::key() }, regionNames );
    if( m_numComponents > 1 )
    {
      fieldsToBeSync.addElementFields( {  extrinsicMeshData::flow::globalCompFraction::key() }, regionNames );
    }
    if( m_enableEnergyBalance )
    {
      fieldsToBeSync.addElementFields( {  extrinsicMeshData::flow::temperature::key() }, regionNames );
    }
    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), true );
  } );
}

void ReactiveCompositionalMultiphaseOBL::initializePostInitialConditionsPreSubGroups()
{
  GEOSX_MARK_FUNCTION;

  FlowSolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  // set mass fraction flag on fluid models
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {

    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addElementFields( { extrinsicMeshData::flow::pressure::key() }, regionNames );
    fieldsToBeSync.addElementFields( { extrinsicMeshData::flow::globalCompFraction::key() }, regionNames );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), false );

    mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const,
                                                                  ElementSubRegionBase & subRegion )
    {
      string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );
      CoupledSolidBase const & solid = getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );

      arrayView1d< integer const > const elemGhostRank = subRegion.ghostRank();
      arrayView1d< real64 const > const volume =  subRegion.getElementVolume();
      arrayView1d< real64 const > const refPorosity =  solid.getReferencePorosity();
      arrayView1d< real64 > const referenceRockVolume =  subRegion.getExtrinsicData< extrinsicMeshData::flow::referenceRockVolume >();
      arrayView1d< real64 > const referencePoreVolume =  subRegion.getExtrinsicData< extrinsicMeshData::flow::referencePoreVolume >();
      arrayView1d< real64 > const referencePorosity =  subRegion.getExtrinsicData< extrinsicMeshData::flow::referencePorosity >();

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
      {
        if( elemGhostRank[ei] >= 0 )
        {
          return;
        }

        referencePoreVolume[ei] = volume[ei] * refPorosity[ei];
        referencePorosity[ei] = refPorosity[ei];
        referenceRockVolume[ei] = volume[ei] - referencePoreVolume[ei];

      } );
    } );
  } );
}


real64 ReactiveCompositionalMultiphaseOBL::solverStep( real64 const & time_n,
                                                       real64 const & dt,
                                                       integer const cycleNumber,
                                                       DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  // Only build the sparsity pattern once
  // TODO: this should be triggered by a topology change indicator
  if( !m_systemSetupDone )
  {
    setupSystem( domain, m_dofManager, m_localMatrix, m_rhs, m_solution );
    m_systemSetupDone = true;
  }

  implicitStepSetup( time_n, dt, domain );

  // currently the only method is implicit time integration
  real64 const dt_return = nonlinearImplicitStep( time_n, dt, cycleNumber, domain );

  // final step for completion of timestep. typically secondary variable updates and cleanup.
  implicitStepComplete( time_n, dt_return, domain );

  return dt_return;
}


void
ReactiveCompositionalMultiphaseOBL::implicitStepSetup( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                                       real64 const & GEOSX_UNUSED_PARAM( dt ),
                                                       DomainPartition & domain )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {

    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 > const & pres_n =
        subRegion.template getExtrinsicData< extrinsicMeshData::flow::pressure_n >();
      arrayView1d< real64 const > const & pres =
        subRegion.template getExtrinsicData< extrinsicMeshData::flow::pressure >();
      pres_n.setValues< parallelDevicePolicy<> >( pres );

      // for single component problems (e.g., geothermal), global component fraction is not a primary variable
      if( m_numComponents > 1 )
      {
        arrayView2d< real64, compflow::USD_COMP > const & compFrac_n =
          subRegion.template getExtrinsicData< extrinsicMeshData::flow::globalCompFraction_n >();
        arrayView2d< real64 const, compflow::USD_COMP > const & compFrac =
          subRegion.template getExtrinsicData< extrinsicMeshData::flow::globalCompFraction >();
        compFrac_n.setValues< parallelDevicePolicy<> >( compFrac );
      }

      if( m_enableEnergyBalance )
      {
        arrayView1d< real64 > const & temp_n =
          subRegion.template getExtrinsicData< extrinsicMeshData::flow::temperature_n >();
        arrayView1d< real64 const > const & temp =
          subRegion.template getExtrinsicData< extrinsicMeshData::flow::temperature >();
        temp_n.setValues< parallelDevicePolicy<> >( temp );
      }

      // update operator values
      updateOBLOperators( subRegion );

      arrayView2d< real64 const, compflow::USD_OBL_VAL > const OBLOperatorValues =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::OBLOperatorValues >();
      arrayView2d< real64, compflow::USD_OBL_VAL > const OBLOperatorValues_n =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::OBLOperatorValues_n >();
      OBLOperatorValues_n.setValues< parallelDevicePolicy<> >( OBLOperatorValues );

    } );

  } );
}

void ReactiveCompositionalMultiphaseOBL::assembleSystem( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                                         real64 const dt,
                                                         DomainPartition & domain,
                                                         DofManager const & dofManager,
                                                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                         arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  assembleAccumulationTerms( dt,
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

void ReactiveCompositionalMultiphaseOBL::assembleAccumulationTerms( real64 const dt,
                                                                    DomainPartition & domain,
                                                                    DofManager const & dofManager,
                                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                    arrayView1d< real64 > const & localRhs ) const
{
  GEOSX_MARK_FUNCTION;


  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase const & subRegion )
    {
      string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

      ElementBasedAssemblyKernelFactory::
        createAndLaunch< parallelDevicePolicy<> >( m_numPhases,
                                                   m_numComponents,
                                                   m_enableEnergyBalance,
                                                   dt,
                                                   dofManager.rankOffset(),
                                                   dofKey,
                                                   subRegion,
                                                   localMatrix,
                                                   localRhs );

    } );
  } );
}

void ReactiveCompositionalMultiphaseOBL::assembleFluxTerms( real64 const dt,
                                                            DomainPartition const & domain,
                                                            DofManager const & dofManager,
                                                            CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                            arrayView1d< real64 > const & localRhs ) const
{
  GEOSX_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               arrayView1d< string const > const & )
  {
    NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

    string const & elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

    fluxApprox.forAllStencils( mesh, [&] ( auto & stencil )
    {
      typename TYPEOFREF( stencil ) ::KernelWrapper stencilWrapper = stencil.createKernelWrapper();

      FaceBasedAssemblyKernelFactory::
        createAndLaunch< parallelDevicePolicy<> >( m_numPhases,
                                                   m_numComponents,
                                                   m_enableEnergyBalance,
                                                   m_transMultExp,
                                                   dofManager.rankOffset(),
                                                   elemDofKey,
                                                   getName(),
                                                   mesh.getElemManager(),
                                                   stencilWrapper,
                                                   dt,
                                                   localMatrix.toViewConstSizes(),
                                                   localRhs.toView() );
    } );
  } );
}


void ReactiveCompositionalMultiphaseOBL::applyBoundaryConditions( real64 const time_n,
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

namespace
{
char const bcLogMessage[] =
  "CompositionalMultiphaseBase {}: at time {}s, "
  "the <{}> boundary condition '{}' is applied to the element set '{}' in subRegion '{}'. "
  "\nThe scale of this boundary condition is {} and multiplies the value of the provided function (if any). "
  "\nThe total number of target elements (including ghost elements) is {}. "
  "\nNote that if this number is equal to zero for all subRegions, the boundary condition will not be applied on this element set.";
}

void ReactiveCompositionalMultiphaseOBL::applySourceFluxBC( real64 const time,
                                                            real64 const dt,
                                                            DofManager const & dofManager,
                                                            DomainPartition & domain,
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

  computeSourceFluxSizeScalingFactor( time,
                                      dt,
                                      domain,
                                      bcNameToBcId,
                                      bcAllSetsSize.toView() );

  // Step 3: we are ready to impose the boundary condition, normalized by the set size

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & )
  {
    fsManager.apply< ElementSubRegionBase >( time + dt,
                                             mesh,
                                             FieldSpecificationBase::viewKeyStruct::fluxBoundaryConditionString(),
                                             [&]( FieldSpecificationBase const & fs,
                                                  string const & setName,
                                                  SortedArrayView< localIndex const > const & targetSet,
                                                  ElementSubRegionBase & subRegion,
                                                  string const & )
    {
      if( fs.getLogLevel() >= 1 && m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
      {
        globalIndex const numTargetElems = MpiWrapper::sum< globalIndex >( targetSet.size() );
        GEOSX_LOG_RANK_0( GEOSX_FMT( bcLogMessage,
                                     getName(), time+dt, SourceFluxBoundaryCondition::catalogName(),
                                     fs.getName(), setName, subRegion.getName(), fs.getScale(), numTargetElems ) );
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

      // Step 3.2: we are ready to add the right-hand side contributions

      // get the normalizer
      real64 const sizeScalingFactor = bcAllSetsSize[bcNameToBcId.at( fs.getName())];

      integer const fluidComponentId = fs.getComponent();
      forAll< parallelDevicePolicy<> >( targetSet.size(), [sizeScalingFactor,
                                                           targetSet,
                                                           rankOffset,
                                                           ghostRank,
                                                           fluidComponentId,
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

        // for all "fluid components", we add the value to the component mass balance equation
        globalIndex const compMassBalanceRow = dofNumber[ei] - rankOffset + fluidComponentId;
        localRhs[compMassBalanceRow] += rhsContributionArrayView[a] / sizeScalingFactor;
      } );
    } );
  } );
}

bool ReactiveCompositionalMultiphaseOBL::validateDirichletBC( DomainPartition & domain,
                                                              real64 const time ) const
{
  constexpr integer MAX_NC = ReactiveCompositionalMultiphaseOBL::MAX_NUM_COMPONENTS + 1; // +1 is for energy component
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  integer const numComp = m_numComponents;
  integer const enableEnergyBalance = m_enableEnergyBalance;

  bool bcConsistent = true;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    // map to check consistent application of BC
    // this is used as bcStatusMap[regionName][subRegionName][setName] which returns the corresponding ComponentMask (if is has been set)
    map< string, map< string, map< string, ComponentMask< MAX_NC > > > > bcStatusMap;
    integer const numCompWithEnergy = numComp + enableEnergyBalance;

    // 1. Check pressure Dirichlet BCs
    fsManager.apply< ElementSubRegionBase >( time,
                                             mesh,
                                             extrinsicMeshData::flow::pressure::key(),
                                             [&]( FieldSpecificationBase const &,
                                                  string const & setName,
                                                  SortedArrayView< localIndex const > const &,
                                                  ElementSubRegionBase & subRegion,
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
      subRegionSetMap[setName].setNumComp( numCompWithEnergy );
    } );

    // 2. Check composition BC (global component fraction)
    fsManager.apply< ElementSubRegionBase >( time,
                                             mesh,
                                             extrinsicMeshData::flow::globalCompFraction::key(),
                                             [&] ( FieldSpecificationBase const & fs,
                                                   string const & setName,
                                                   SortedArrayView< localIndex const > const &,
                                                   ElementSubRegionBase & subRegion,
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

    if( enableEnergyBalance )
    {
      // 3. Check temperature Dirichlet BCs
      fsManager.apply< ElementSubRegionBase >( time,
                                               mesh,
                                               extrinsicMeshData::flow::temperature::key(),
                                               [&]( FieldSpecificationBase const &,
                                                    string const & setName,
                                                    SortedArrayView< localIndex const > const &,
                                                    ElementSubRegionBase & subRegion,
                                                    string const & )
      {
        // 1.0. Check whether pressure has already been applied to this set
        string const & subRegionName = subRegion.getName();
        string const & regionName = subRegion.getParent().getParent().getName();

        auto & subRegionSetMap = bcStatusMap[regionName][subRegionName];
        if( subRegionSetMap.count( setName ) == 0 )
        {
          bcConsistent = false;
          GEOSX_WARNING( GEOSX_FMT( "Pressure boundary condition not prescribed on set {}/{}/{}", regionName, subRegionName, setName ) );
        }

        ComponentMask< MAX_NC > & compMask = subRegionSetMap[setName];
        // if energy balance is enabled, energy is the last component in the mask
        if( compMask[numComp] )
        {
          bcConsistent = false;
          GEOSX_WARNING( GEOSX_FMT( "Conflicting temperature boundary conditions on set {}/{}/{}", regionName, subRegionName, setName ) );
        }
        compMask.set( numComp );
      } );
    }

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
  } );
  return bcConsistent;
}

void ReactiveCompositionalMultiphaseOBL::applyDirichletBC( real64 const time,
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
    bool const bcConsistent = validateDirichletBC( domain, time + dt );
    GEOSX_ERROR_IF( !bcConsistent, GEOSX_FMT( "CompositionalMultiphaseBase {}: inconsistent boundary conditions", getName() ) );
  }

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {

    // 1. Apply pressure Dirichlet BCs, store in a separate field
    fsManager.apply< ElementSubRegionBase >( time + dt,
                                             mesh,
                                             extrinsicMeshData::flow::pressure::key(),
                                             [&]( FieldSpecificationBase const & fs,
                                                  string const & setName,
                                                  SortedArrayView< localIndex const > const & targetSet,
                                                  ElementSubRegionBase & subRegion,
                                                  string const & )
    {
      if( fs.getLogLevel() >= 1 && m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
      {
        globalIndex const numTargetElems = MpiWrapper::sum< globalIndex >( targetSet.size() );
        GEOSX_LOG_RANK_0( GEOSX_FMT( bcLogMessage,
                                     getName(), time+dt, FieldSpecificationBase::catalogName(),
                                     fs.getName(), setName, subRegion.getName(), fs.getScale(), numTargetElems ) );
      }

      fs.applyFieldValue< FieldSpecificationEqual, parallelDevicePolicy<> >( targetSet,
                                                                             time + dt,
                                                                             subRegion,
                                                                             extrinsicMeshData::flow::bcPressure::key() );
    } );

    // 2. Apply composition BC (global component fraction), store in a separate field
    fsManager.apply< ElementSubRegionBase >( time + dt,
                                             mesh,
                                             extrinsicMeshData::flow::globalCompFraction::key(),
                                             [&] ( FieldSpecificationBase const & fs,
                                                   string const &,
                                                   SortedArrayView< localIndex const > const & targetSet,
                                                   ElementSubRegionBase & subRegion,
                                                   string const & )
    {
      fs.applyFieldValue< FieldSpecificationEqual, parallelDevicePolicy<> >( targetSet,
                                                                             time + dt,
                                                                             subRegion,
                                                                             extrinsicMeshData::flow::bcGlobalCompFraction::key() );
    } );

    // 3. Apply temperature Dirichlet BCs, store in a separate field
    fsManager.apply< ElementSubRegionBase >( time + dt,
                                             mesh,
                                             extrinsicMeshData::flow::temperature::key(),
                                             [&]( FieldSpecificationBase const & fs,
                                                  string const & setName,
                                                  SortedArrayView< localIndex const > const & targetSet,
                                                  ElementSubRegionBase & subRegion,
                                                  string const & )
    {
      if( fs.getLogLevel() >= 1 && m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
      {
        globalIndex const numTargetElems = MpiWrapper::sum< globalIndex >( targetSet.size() );
        GEOSX_LOG_RANK_0( GEOSX_FMT( bcLogMessage,
                                     getName(), time+dt, FieldSpecificationBase::catalogName(),
                                     fs.getName(), setName, subRegion.getName(), fs.getScale(), numTargetElems ) );
      }

      fs.applyFieldValue< FieldSpecificationEqual, parallelDevicePolicy<> >( targetSet,
                                                                             time + dt,
                                                                             subRegion,
                                                                             extrinsicMeshData::flow::bcTemperature::key() );
    } );

    globalIndex const rankOffset = dofManager.rankOffset();
    string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );


    // 3. Apply to the system
    fsManager.apply< ElementSubRegionBase >( time + dt,
                                             mesh,
                                             extrinsicMeshData::flow::pressure::key(),
                                             [&] ( FieldSpecificationBase const &,
                                                   string const &,
                                                   SortedArrayView< localIndex const > const & targetSet,
                                                   ElementSubRegionBase & subRegion,
                                                   string const & )
    {
      arrayView1d< real64 const > const bcPres = subRegion.getExtrinsicData< extrinsicMeshData::flow::bcPressure >();
      arrayView2d< real64 const, compflow::USD_COMP > const bcCompFrac = subRegion.getExtrinsicData< extrinsicMeshData::flow::bcGlobalCompFraction >();
      arrayView1d< real64 const > const bcTemp = subRegion.getExtrinsicData< extrinsicMeshData::flow::bcTemperature >();

      arrayView1d< real64 const > const pres = subRegion.getExtrinsicData< extrinsicMeshData::flow::pressure >();
      arrayView2d< real64 const, compflow::USD_COMP > const compFrac = subRegion.getExtrinsicData< extrinsicMeshData::flow::globalCompFraction >();
      arrayView1d< real64 const > const temp = subRegion.getExtrinsicData< extrinsicMeshData::flow::temperature >();

      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();
      arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );

      integer const numComp = m_numComponents;
      integer const enableEnergyBalance = m_enableEnergyBalance;

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
                                                    pres[ei] );
        localRhs[localRow] = rhsValue;

        // 3.2. For each component (but the last), apply target global component fraction value
        //      Note: the FieldSpecification for the last component is ignored
        for( integer ic = 0; ic < numComp - 1; ++ic )
        {
          FieldSpecificationEqual::SpecifyFieldValue( dofIndex + ic + 1,
                                                      rankOffset,
                                                      localMatrix,
                                                      rhsValue,
                                                      bcCompFrac[ei][ic],
                                                      compFrac[ei][ic] );
          localRhs[localRow + ic + 1] = rhsValue;
        }

        if( enableEnergyBalance )
        {
          // 3.3. If energy balance is enabled, apply BC temperature to the last equation
          FieldSpecificationEqual::SpecifyFieldValue( dofIndex + numComp,
                                                      rankOffset,
                                                      localMatrix,
                                                      rhsValue,
                                                      bcTemp[ei],
                                                      temp[ei] );
          localRhs[localRow + numComp] = rhsValue;
        }
      } );
    } );
  } );
}

// to be changed into enforceOBLLimits - to chop all primary variables to be within OBL discretization space
void ReactiveCompositionalMultiphaseOBL::chopPrimaryVariablesToOBLLimits( DomainPartition & domain )
{
  GEOSX_UNUSED_VAR( domain );

  // not implemented yet
}

void ReactiveCompositionalMultiphaseOBL::resetStateToBeginningOfStep( DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion,
                                                SurfaceElementSubRegion >( regionNames,
                                                                           [&]( localIndex const,
                                                                                auto & subRegion )
    {
      arrayView1d< real64 const > const & pres_n =
        subRegion.template getExtrinsicData< extrinsicMeshData::flow::pressure_n >();
      arrayView1d< real64 > const & pres =
        subRegion.template getExtrinsicData< extrinsicMeshData::flow::pressure >();
      pres.setValues< parallelDevicePolicy<> >( pres_n );

      // for single component problems (e.g., geothermal), global component fraction is not a primary variable
      if( m_numComponents > 1 )
      {
        arrayView2d< real64 const, compflow::USD_COMP > const & compFrac_n =
          subRegion.template getExtrinsicData< extrinsicMeshData::flow::globalCompFraction_n >();
        arrayView2d< real64, compflow::USD_COMP > const & compFrac =
          subRegion.template getExtrinsicData< extrinsicMeshData::flow::globalCompFraction >();
        compFrac.setValues< parallelDevicePolicy<> >( compFrac_n );
      }

      if( m_enableEnergyBalance )
      {
        arrayView1d< real64 const > const & temp_n =
          subRegion.template getExtrinsicData< extrinsicMeshData::flow::temperature_n >();
        arrayView1d< real64 > const & temp =
          subRegion.template getExtrinsicData< extrinsicMeshData::flow::temperature >();
        temp.setValues< parallelDevicePolicy<> >( temp_n );
      }

      // update operator values
      updateOBLOperators( subRegion );

    } );
  } );
}

void ReactiveCompositionalMultiphaseOBL::updateOBLOperators( ObjectManagerBase & dataGroup ) const
{
  GEOSX_MARK_FUNCTION;

  OBLOperatorsKernelFactory::
    createAndLaunch< parallelDevicePolicy<> >( m_numPhases,
                                               m_numComponents,
                                               m_enableEnergyBalance,
                                               dataGroup,
                                               *m_OBLOperatorsTable );

}


void ReactiveCompositionalMultiphaseOBL::updateState( DomainPartition & domain )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion,
                                                SurfaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                             auto & subRegion )
    {
      // update operator values
      updateOBLOperators( subRegion );
    } );
  } );
}



//START_SPHINX_INCLUDE_01
REGISTER_CATALOG_ENTRY( SolverBase, ReactiveCompositionalMultiphaseOBL, string const &, Group * const )
//END_SPHINX_INCLUDE_01
}// namespace geosx
