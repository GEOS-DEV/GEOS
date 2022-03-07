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
 * @file OBLSuperEngine.cpp
 */

#include "OBLSuperEngine.hpp"

#include "constitutive/solid/CoupledSolidBase.hpp"
#include "dataRepository/Group.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
 #include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/fluidFlow/OBLSuperEngineExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseExtrinsicData.hpp"
#include "physicsSolvers/fluidFlow/OBLSuperEngineKernels.hpp"


namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace OBLSuperEngineKernels;

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

OBLSuperEngine::OBLSuperEngine( const string & name,
                                Group * const parent )
  :
  FlowSolverBase( name, parent ),
  m_numPhases( 0 ),
  m_numComponents( 0 ),
  m_maxCompFracChange( 1.0 ),
  m_minScalingFactor( 0.01 ),
  m_allowOBLChopping( 1 )
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

void OBLSuperEngine::initializePreSubGroups()
{
  FlowSolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  if( !fvManager.hasGroup< FluxApproximationBase >( m_discretizationName ) )
  {
    GEOSX_ERROR( "A discretization deriving from FluxApproximationBase must be selected with OBLSuperEngine" );
  }

}

void OBLSuperEngine::setupDofs( DomainPartition const & domain,
                                DofManager & dofManager ) const
{
  dofManager.addField( viewKeyStruct::elemDofFieldString(),
                       DofManager::Location::Elem,
                       m_numDofPerCell,
                       m_meshTargets );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );
  dofManager.addCoupling( viewKeyStruct::elemDofFieldString(), fluxApprox );
}


void OBLSuperEngine::implicitStepComplete( real64 const & time,
                                           real64 const & dt,
                                           DomainPartition & domain )
{
  GEOSX_UNUSED_VAR( time );
  GEOSX_UNUSED_VAR( dt );
  integer const numComp = m_numComponents;
  integer const enableEnergyBalance = m_enableEnergyBalance;

  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                               MeshLevel & mesh,
                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 const > const dPres =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::deltaPressure >();
      arrayView2d< real64 const, compflow::USD_COMP > const dCompFrac =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::deltaGlobalCompFraction >();
      arrayView1d< real64 const > const dTemp =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::deltaTemperature >();

      arrayView1d< real64 > const pres =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::pressure >();
      arrayView2d< real64, compflow::USD_COMP > const compFrac =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::globalCompFraction >();
      arrayView1d< real64 > const temp =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::temperature >();

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
      {
        pres[ei] += dPres[ei];
        for( integer ic = 0; ic < numComp; ++ic )
        {
          compFrac[ei][ic] += dCompFrac[ei][ic];
        }
        if( enableEnergyBalance )
        {
          temp[ei] += dTemp[ei];
        }
      } );

      string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );
      CoupledSolidBase const & porousMaterial = getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );
      porousMaterial.saveConvergedState();
    } );
  } );
}

void OBLSuperEngine::postProcessInput()
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

  GEOSX_THROW_IF_NE_MSG( m_numDofPerCell, m_OBLOperatorsTable->getNumDims(),
                         GEOSX_FMT( "The number of degrees of freedom per element used in solver - {} - and in operator table - {} - should match", m_numDofPerCell, m_OBLOperatorsTable->getNumDims()),
                         InputError );

  GEOSX_THROW_IF_NE_MSG( m_numOBLOperators, m_OBLOperatorsTable->getNumOps(),
                         GEOSX_FMT( "The number of operators per element used in solver - {} - and in operator table - {} - should match", m_numOBLOperators, m_OBLOperatorsTable->getNumOps()),
                         InputError );

}

void OBLSuperEngine::registerDataOnMesh( Group & meshBodies )
{
  using namespace extrinsicMeshData::flow;
  // 1. Call base class method
  FlowSolverBase::registerDataOnMesh( meshBodies );

  // 2. Register and resize all fields as necessary
  forMeshTargets( meshBodies, [&]( string const &,
                                   MeshLevel & mesh,
                                   arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames,
                                                [&]( localIndex const,
                                                     ElementSubRegionBase & subRegion )
    {
      subRegion.registerExtrinsicData< pressure >( getName() );
      subRegion.registerExtrinsicData< initialPressure >( getName() );
      subRegion.registerExtrinsicData< deltaPressure >( getName() );
      subRegion.registerExtrinsicData< bcPressure >( getName() );

      subRegion.registerExtrinsicData< temperature >( getName() );
      subRegion.registerExtrinsicData< bcTemperature >( getName() );

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

      subRegion.registerExtrinsicData< bcGlobalCompFraction >( getName() ).
        reference().resizeDimension< 1 >( m_numComponents );

      // we need to register this fiels in any case (if there is a single component or not)
      // to be able to pass the view to OBLOperatorsKernel
      subRegion.registerExtrinsicData< deltaGlobalCompFraction >( getName() ).
        reference().resizeDimension< 1 >( m_numComponents );
      // in principle, referencePorosity could be used directly from solid model,
      // but was duplicated to remove dependency on solid
      subRegion.registerExtrinsicData< referencePorosity >( getName() );

      // referencePoreVolume and referenceRockVolume are introduced for the sake of performance:
      // this way the multiplication of constant arrays (e.g., referencePorosity and volume) every Newton step is avoided
      subRegion.registerExtrinsicData< referencePoreVolume >( getName() );
      subRegion.registerExtrinsicData< referenceRockVolume >( getName() );

      // thermal rock properties (again, register in any case)
      // it is not possible to use specificHeatCapacity from solid model here, because specificHeatCapacity includes several quantities,
      // which are split in OBL framework: constant rock volume and,
      // hidden inside operator - therefore variable - rock compressibility and rock energy
      subRegion.registerExtrinsicData< rockVolumetricHeatCapacity >( getName() );
      subRegion.registerExtrinsicData< rockThermalConductivity >( getName() );
      subRegion.registerExtrinsicData< rockKineticRateFactor >( getName() );

    } );

    FaceManager & faceManager = mesh.getFaceManager();
    {
      faceManager.registerExtrinsicData< facePressure >( getName() );
    }

  } );
}

real64 OBLSuperEngine::calculateResidualNorm( DomainPartition const & domain,
                                              DofManager const & dofManager,
                                              arrayView1d< real64 const > const & localRhs )
{
  GEOSX_MARK_FUNCTION;
  GEOSX_UNUSED_VAR( localRhs );

  real64 localResidualNorm = 0.0;

  globalIndex const rankOffset = dofManager.rankOffset();
  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
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

        ResidualDARTSL2NormKernel::launch< parallelDevicePolicy<>,
                                           parallelDeviceReduce >( localRhs,
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
        arrayView2d< real64 const, compflow::USD_OBL_VAL > const OBLValsOld = subRegion.getExtrinsicData< extrinsicMeshData::flow::OBLOperatorValuesOld >();

        ResidualNormKernel::launch< parallelDevicePolicy<>,
                                    parallelDeviceReduce >( localRhs,
                                                            rankOffset,
                                                            m_numDofPerCell,
                                                            dofNumber,
                                                            elemGhostRank,
                                                            refPoreVolume,
                                                            OBLValsOld,
                                                            subRegionResidualNorm );
        localResidualNorm += subRegionResidualNorm;
      }
    } );
  } );

  real64 const residual = m_useDARTSL2Norm ? MpiWrapper::max( localResidualNorm ) : std::sqrt( MpiWrapper::sum( localResidualNorm ) );

  if( getLogLevel() >= 1 && logger::internal::rank == 0 )
  {
    std::cout << GEOSX_FMT( "    ( Rfluid ) = ( {:4.2e} ) ;", residual );
  }

  return residual;
}


real64 OBLSuperEngine::scalingForSystemSolution( DomainPartition const & domain,
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

  real64 constexpr eps = OBLSuperEngineKernels::minValueForDivision;
  real64 const maxCompFracChange = m_maxCompFracChange;

  localIndex const NC = m_numComponents;

  globalIndex const rankOffset = dofManager.rankOffset();
  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
  real64 scalingFactor = 1.0;

  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
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

bool OBLSuperEngine::checkSystemSolution( DomainPartition const & domain,
                                          DofManager const & dofManager,
                                          arrayView1d< real64 const > const & localSolution,
                                          real64 const scalingFactor )
{
  GEOSX_MARK_FUNCTION;

  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
  localIndex localCheck = 1;

  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
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
      arrayView1d< real64 const > const & dPres =
        subRegion.getReference< array1d< real64 > >( extrinsicMeshData::flow::deltaPressure::key() );
      arrayView2d< real64 const, compflow::USD_COMP > const & compFrac =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::globalCompFraction >();
      arrayView2d< real64 const, compflow::USD_COMP > const & dCompFrac =
        subRegion.getExtrinsicData< extrinsicMeshData::flow::deltaGlobalCompFraction >();
      arrayView1d< real64 const > const & temp =
        subRegion.getReference< array1d< real64 > >( extrinsicMeshData::flow::temperature::key() );
      arrayView1d< real64 const > const & dTemp =
        subRegion.getReference< array1d< real64 > >( extrinsicMeshData::flow::deltaTemperature::key() );

      localIndex const subRegionSolutionCheck = SolutionCheckKernel::launch< parallelDevicePolicy<>,
                                                                             parallelDeviceReduce >( localSolution,
                                                                                                     dofManager.rankOffset(),
                                                                                                     m_numComponents,
                                                                                                     m_enableEnergyBalance,
                                                                                                     dofNumber,
                                                                                                     elemGhostRank,
                                                                                                     pres,
                                                                                                     dPres,
                                                                                                     compFrac,
                                                                                                     dCompFrac,
                                                                                                     temp,
                                                                                                     dTemp,
                                                                                                     m_allowOBLChopping,
                                                                                                     scalingFactor );

      localCheck = std::min( localCheck, subRegionSolutionCheck );
    } );
  } );

  return MpiWrapper::min( localCheck, MPI_COMM_GEOSX );
}

void OBLSuperEngine::applySystemSolution( DofManager const & dofManager,
                                          arrayView1d< real64 const > const & localSolution,
                                          real64 const scalingFactor,
                                          DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  DofManager::CompMask pressureMask( m_numDofPerCell, 0, 1 );
  std::map< string, string_array > fieldNames;

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::elemDofFieldString(),
                               extrinsicMeshData::flow::deltaPressure::key(),
                               scalingFactor,
                               pressureMask );
  fieldNames["elems"].emplace_back( extrinsicMeshData::flow::deltaPressure::key() );

  if( m_numComponents > 1 )
  {
    DofManager::CompMask compFracMask( m_numDofPerCell, 1, m_numComponents );

    dofManager.addVectorToField( localSolution,
                                 viewKeyStruct::elemDofFieldString(),
                                 extrinsicMeshData::flow::deltaGlobalCompFraction::key(),
                                 scalingFactor,
                                 compFracMask );
    fieldNames["elems"].emplace_back( extrinsicMeshData::flow::deltaGlobalCompFraction::key() );

  }

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

void OBLSuperEngine::initializePostInitialConditionsPreSubGroups()
{
  GEOSX_MARK_FUNCTION;

  FlowSolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( extrinsicMeshData::flow::pressure::key() );
  fieldNames["elems"].emplace_back( extrinsicMeshData::flow::globalCompFraction::key() );

  // set mass fraction flag on fluid models
  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                               MeshLevel & mesh,
                                               arrayView1d< string const > const & regionNames )
  {
    CommunicationTools::getInstance().synchronizeFields( fieldNames, mesh, domain.getNeighbors(), false );

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


real64 OBLSuperEngine::solverStep( real64 const & time_n,
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

void OBLSuperEngine::backupFields( MeshLevel & mesh,
                                   arrayView1d< string const > const & regionNames ) const
{
  GEOSX_MARK_FUNCTION;

  //backup some fields used in time derivative approximation
  mesh.getElemManager().forElementSubRegions( regionNames,
                                              [&]( localIndex const,
                                                   ElementSubRegionBase & subRegion )
  {
    arrayView2d< real64 const, compflow::USD_OBL_VAL > const OBLOperatorValues =
      subRegion.getExtrinsicData< extrinsicMeshData::flow::OBLOperatorValues >();

    arrayView2d< real64, compflow::USD_OBL_VAL > const OBLOperatorValuesOld =
      subRegion.getExtrinsicData< extrinsicMeshData::flow::OBLOperatorValuesOld >();

    OBLOperatorValuesOld.setValues< parallelDevicePolicy<> >( OBLOperatorValues );
  } );
}

void
OBLSuperEngine::implicitStepSetup( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                   real64 const & GEOSX_UNUSED_PARAM( dt ),
                                   DomainPartition & domain )
{
  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                               MeshLevel & mesh,
                                               arrayView1d< string const > const & regionNames )
  {

    // set deltas to zero and recompute dependent quantities
    resetStateToBeginningOfStep( domain );

    // backup fields used in time derivative approximation
    backupFields( mesh, regionNames );
  } );
}

void OBLSuperEngine::assembleSystem( real64 const GEOSX_UNUSED_PARAM( time_n ),
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

void OBLSuperEngine::assembleAccumulationTerms( real64 const dt,
                                                DomainPartition & domain,
                                                DofManager const & dofManager,
                                                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                arrayView1d< real64 > const & localRhs ) const
{
  GEOSX_MARK_FUNCTION;


  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
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

void OBLSuperEngine::assembleFluxTerms( real64 const dt,
                                        DomainPartition const & domain,
                                        DofManager const & dofManager,
                                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                        arrayView1d< real64 > const & localRhs ) const
{
  GEOSX_MARK_FUNCTION;

  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                               MeshLevel const & mesh,
                                               arrayView1d< string const > const & )
  {
    NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
    FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
    FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

    string const & elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > elemDofNumber =
      mesh.getElemManager().constructArrayViewAccessor< globalIndex, 1 >( elemDofKey );
    elemDofNumber.setName( getName() + "/accessors/" + elemDofKey );

    fluxApprox.forAllStencils( mesh, [&] ( auto & stencil )
    {
      typename TYPEOFREF( stencil ) ::StencilWrapper stencilWrapper = stencil.createStencilWrapper();

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


void OBLSuperEngine::applyBoundaryConditions( real64 const time_n,
                                              real64 const dt,
                                              DomainPartition & domain,
                                              DofManager const & dofManager,
                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                              arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  // apply pressure boundary conditions.
  applyDirichletBC( time_n, dt, dofManager, domain, localMatrix.toViewConstSizes(), localRhs.toView() );
}

namespace internal
{
string const bcLogMessage = string( "OBLSuperEngine {}: at time {}s, " )
                            + string( "the <{}> boundary condition '{}' is applied to the element set '{}' in subRegion '{}'. " )
                            + string( "\nThe scale of this boundary condition is {} and multiplies the value of the provided function (if any). " )
                            + string( "\nThe total number of target elements (including ghost elements) is {}. " )
                            + string( "\nNote that if this number is equal to zero for all subRegions, the boundary condition will not be applied on this element set." );
}

namespace
{

bool validateDirichletBC( DomainPartition & domain,
                          integer const numComp,
                          integer const enableEnergyBalance,
                          real64 const time )
{
  constexpr integer MAX_NC = OBLSuperEngine::MAX_NUM_COMPONENTS + 1; // +1 is for energy component
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  map< string, map< string, map< string, ComponentMask< MAX_NC > > > > bcStatusMap;   // map to check consistent application of BC
  bool bcConsistent = true;
  integer const numCompWithEnergy = numComp + enableEnergyBalance;

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
    subRegionSetMap[setName].setNumComp( numCompWithEnergy );
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
      return;     // can't check next part with invalid component id
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
    fsManager.apply( time,
                     domain,
                     "ElementRegions",
                     extrinsicMeshData::flow::temperature::key(),
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

  return bcConsistent;
}

}


void OBLSuperEngine::applyDirichletBC( real64 const time,
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
    bool const bcConsistent = validateDirichletBC( domain, m_numComponents, m_enableEnergyBalance, time + dt );
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

  // 2. Apply composition BC (global component fraction), store in a separate field
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
                                                                           extrinsicMeshData::flow::bcGlobalCompFraction::key() );
  } );

  // 3. Apply temperature Dirichlet BCs, store in a separate field
  fsManager.apply( time + dt,
                   domain,
                   "ElementRegions",
                   extrinsicMeshData::flow::temperature::key(),
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
                                                                           extrinsicMeshData::flow::bcTemperature::key() );
  } );

  globalIndex const rankOffset = dofManager.rankOffset();
  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );


  // 3. Apply to the system
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
    arrayView1d< real64 const > const bcPres =
      subRegion.getReference< array1d< real64 > >( extrinsicMeshData::flow::bcPressure::key() );
    arrayView2d< real64 const, compflow::USD_COMP > const bcCompFrac =
      subRegion.getReference< array2d< real64, compflow::LAYOUT_COMP > >( extrinsicMeshData::flow::bcGlobalCompFraction::key() );
    arrayView1d< real64 const > const bcTemp =
      subRegion.getReference< array1d< real64 > >( extrinsicMeshData::flow::bcTemperature::key() );

    arrayView1d< real64 const > const pres =
      subRegion.getReference< array1d< real64 > >( extrinsicMeshData::flow::pressure::key() );
    arrayView2d< real64 const, compflow::USD_COMP > const compFrac =
      subRegion.getReference< array2d< real64, compflow::LAYOUT_COMP > >( extrinsicMeshData::flow::globalCompFraction::key() );
    arrayView1d< real64 const > const temp =
      subRegion.getReference< array1d< real64 > >( extrinsicMeshData::flow::temperature::key() );



    arrayView1d< real64 const > const dPres =
      subRegion.getReference< array1d< real64 > >( extrinsicMeshData::flow::deltaPressure::key() );
    arrayView2d< real64 const, compflow::USD_COMP > const dCompFrac =
      subRegion.getReference< array2d< real64, compflow::LAYOUT_COMP > >( extrinsicMeshData::flow::deltaGlobalCompFraction::key() );
    arrayView1d< real64 const > const dTemp =
      subRegion.getReference< array1d< real64 > >( extrinsicMeshData::flow::deltaTemperature::key() );

    arrayView1d< integer const > const ghostRank =
      subRegion.getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString() );
    arrayView1d< globalIndex const > const dofNumber =
      subRegion.getReference< array1d< globalIndex > >( dofKey );

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
                                                  pres[ei] + dPres[ei] );
      localRhs[localRow] = rhsValue;

      // 3.2. For each component (but the last), apply target global component fraction value
      for( integer ic = 0; ic < numComp - 1; ++ic )
      {
        FieldSpecificationEqual::SpecifyFieldValue( dofIndex + ic + 1,
                                                    rankOffset,
                                                    localMatrix,
                                                    rhsValue,
                                                    bcCompFrac[ei][ic],
                                                    compFrac[ei][ic] + dCompFrac[ei][ic] );
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
                                                    temp[ei] + dTemp[ei] );
        localRhs[localRow + numComp] = rhsValue;
      }
    } );
  } );
}

void OBLSuperEngine::solveSystem( DofManager const & dofManager,
                                  ParallelMatrix & matrix,
                                  ParallelVector & rhs,
                                  ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  rhs.scale( -1.0 );
  solution.zero();

  SolverBase::solveSystem( dofManager, matrix, rhs, solution );
}

// to be changed into enforceOBLLimits - to chop all primary variables to be within OBL discretization space
void OBLSuperEngine::chopPrimaryVariablesToOBLLimits( DomainPartition & domain )
{
  GEOSX_UNUSED_VAR( domain );
  // GEOSX_MARK_FUNCTION;

  // MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  // integer const numComp = m_numComponents;
  // forTargetSubRegions( mesh, [&]( localIndex const, ElementSubRegionBase & subRegion )
  // {
  //   arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

  //   arrayView2d< real64 const, compflow::USD_COMP > const compDens =
  //     subRegion.getExtrinsicData< extrinsicMeshData::flow::globalCompDensity >();
  //   arrayView2d< real64, compflow::USD_COMP > const dCompDens =
  //     subRegion.getExtrinsicData< extrinsicMeshData::flow::deltaGlobalCompDensity >();

  //   forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
  //   {
  //     if( ghostRank[ei] < 0 )
  //     {
  //       for( integer ic = 0; ic < numComp; ++ic )
  //       {
  //         real64 const newDens = compDens[ei][ic] + dCompDens[ei][ic];
  //         if( newDens < minValueForDivision )
  //         {
  //           dCompDens[ei][ic] = -compDens[ei][ic] + minValueForDivision;
  //         }
  //       }
  //     }
  //   } );
  // } );
}

void OBLSuperEngine::resetStateToBeginningOfStep( DomainPartition & domain )
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

      // update operator values
      updateOBLOperators( subRegion );

    } );
  } );
}

void OBLSuperEngine::updateOBLOperators( ObjectManagerBase & dataGroup ) const
{
  GEOSX_MARK_FUNCTION;

  OBLOperatorsKernelFactory::
    createAndLaunch< parallelDevicePolicy<> >( m_numPhases,
                                               m_numComponents,
                                               m_enableEnergyBalance,
                                               dataGroup,
                                               *m_OBLOperatorsTable );

}


void OBLSuperEngine::updateState( DomainPartition & domain )
{
  forMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                               MeshLevel & mesh,
                                               arrayView1d< string const > const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion,
                                                SurfaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                             auto & subRegion )
    {
      // update porosity and permeability
      updatePorosityAndPermeability( subRegion );

      // update operator values
      updateOBLOperators( subRegion );
    } );
  } );
}



//START_SPHINX_INCLUDE_01
REGISTER_CATALOG_ENTRY( SolverBase, OBLSuperEngine, string const &, Group * const )
//END_SPHINX_INCLUDE_01
}// namespace geosx
