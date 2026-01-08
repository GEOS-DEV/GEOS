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
 * @file WellSolverBase.cpp
 */

#include "WellSolverBase.hpp"
#include "physicsSolvers/LogLevelsInfo.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/PerforationFields.hpp"
#include "mesh/WellElementRegion.hpp"
#include "mesh/WellElementSubRegion.hpp"
#include "physicsSolvers/fluidFlow/wells/LogLevelsInfo.hpp"
#include "physicsSolvers/fluidFlow/wells/WellControls.hpp"
#include "physicsSolvers/fluidFlow/wells/WellSolverBaseFields.hpp"
#include "fileIO/Outputs/OutputBase.hpp"
#include "functions/FunctionManager.hpp"
namespace geos
{

using namespace dataRepository;
using namespace fields;

WellSolverBase::WellSolverBase( string const & name,
                                Group * const parent )
  : PhysicsSolverBase( name, parent ),
  m_numPhases( 0 ),
  m_numComponents( 0 ),
  m_numDofPerWellElement( 0 ),
  m_numDofPerResElement( 0 ),
  m_isThermal( 0 ),
  m_ratesOutputDir( joinPath( OutputBase::getOutputDirectory(), name + "_rates" ) ),
  m_keepVariablesConstantDuringInitStep( false ),
  m_writeSegDebug( 0 ),
  m_globalNumTimeSteps( -1 ),
  m_currentDt( -1.0 ),
  my_ctime( 0 ),
  m_nextDt( -1 )
{
  registerWrapper( viewKeyStruct::isThermalString(), &m_isThermal ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag indicating whether the problem is thermal or not." );

  this->getWrapper< string >( viewKeyStruct::discretizationString() ).
    setInputFlag( InputFlags::FALSE );

  this->registerWrapper( viewKeyStruct::writeCSVFlagString(), &m_writeCSV ).
    setApplyDefaultValue( 0 ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "When set to 1, write the rates into a CSV file." );

  this->registerWrapper( viewKeyStruct::timeStepFromTablesFlagString(), &m_timeStepFromTables ).
    setApplyDefaultValue( 0 ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Choose time step to honor rates/bhp tables time intervals" );

  this->registerWrapper( viewKeyStruct::writeSegDebugFlagString(), &m_writeSegDebug ).
    setApplyDefaultValue( 0 ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setDescription( "Write well seg/perf debug into CSV files" );

  addLogLevel< logInfo::WellControl >();
}

Group * WellSolverBase::createChild( string const & childKey, string const & childName )
{
  static std::set< string > const childTypes = {
    keys::wellControls,
    PhysicsSolverBase::groupKeyStruct::linearSolverParametersString(),
    PhysicsSolverBase::groupKeyStruct::nonlinearSolverParametersString(),
  };
  GEOS_ERROR_IF( childTypes.count( childKey ) == 0,
                 CatalogInterface::unknownTypeError( childKey, getDataContext(), childTypes ),
                 getDataContext() );
  if( childKey == keys::wellControls )
  {
    return &registerGroup< WellControls >( childName );
  }
  else
  {
    PhysicsSolverBase::createChild( childKey, childName );
    return nullptr;
  }
}

void WellSolverBase::expandObjectCatalogs()
{
  createChild( keys::wellControls, keys::wellControls );
}

WellSolverBase::~WellSolverBase() = default;

void WellSolverBase::postInputInitialization()
{
  PhysicsSolverBase::postInputInitialization();

  // 1. Set key dimensions of the problem
  m_numDofPerWellElement = m_isThermal ?    m_numComponents + 2 : m_numComponents + 1; // 1 pressure  connectionRate + temp if thermal
  m_numDofPerResElement = m_isThermal ? m_numComponents  + 1: m_numComponents;   // 1 pressure   + temp if thermal
  m_writeSegDebug=2;
  m_writeCSV=1;
  if( m_writeSegDebug > 0 )
  {
    if( m_writeCSV == 0 )
    {
      m_writeCSV=1;
    }
  }
  // create dir for rates output
  if( m_writeCSV > 0 )
  {
    if( MpiWrapper::commRank() == 0 )
    {
      makeDirsForPath( m_ratesOutputDir );
    }
    // wait till the dir is created by rank 0
    MPI_Barrier( MPI_COMM_WORLD );
  }
}

void WellSolverBase::registerDataOnMesh( Group & meshBodies )
{
  PhysicsSolverBase::registerDataOnMesh( meshBodies );

  // loop over the wells
  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & meshLevel,
                                                    string_array const & regionNames )
  {

    ElementRegionManager & elementRegionManager = meshLevel.getElemManager();
    elementRegionManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                                       [&]( localIndex const,
                                                                            WellElementSubRegion & subRegion )
    {
      subRegion.registerField< well::pressure >( getName() );
      subRegion.registerField< well::pressure_n >( getName() );

      subRegion.registerField< well::temperature >( getName() );
      if( isThermal() )
      {
        subRegion.registerField< well::temperature_n >( getName() );
      }

      subRegion.registerField< well::gravityCoefficient >( getName() );

      PerforationData * const perforationData = subRegion.getPerforationData();
      perforationData->registerField< well::gravityCoefficient >( getName() );
    } );
  } );
}

void WellSolverBase::initializePostSubGroups()
{
  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  FunctionManager & functionManager = FunctionManager::getInstance();
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();
    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )
    {
      validateWellConstraints( 0, 0, subRegion );

      // validate perforation status table
      PerforationData & perforationData = *subRegion.getPerforationData();
      string_array const & perfStatusTableName = perforationData.getPerfStatusTableName();
      for( integer i=0; i<perforationData.size(); i++ )
      {
        TableFunction * tableFunction =  functionManager.getGroupPointer< TableFunction >( perfStatusTableName[i] );
        GEOS_THROW_IF( tableFunction->getInterpolationMethod() != TableFunction::InterpolationType::Lower,
                       "WellSolverBase " << getDataContext() << ": The interpolation method for the perforation status table "
                                         << tableFunction->getName() << " should be TableFunction::InterpolationType::Lower",
                       InputError );
      }
    } );
  } );
}

void WellSolverBase::setupDofs( DomainPartition const & domain,
                                DofManager & dofManager ) const
{
  map< std::pair< string, string >, string_array > meshTargets;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const & meshBodyName,
                                                               MeshLevel const & meshLevel,
                                                               string_array const & regionNames )
  {
    string_array regions;
    ElementRegionManager const & elementRegionManager = meshLevel.getElemManager();
    elementRegionManager.forElementRegions< WellElementRegion >( regionNames,
                                                                 [&]( localIndex const,
                                                                      WellElementRegion const & region )
    {
      regions.emplace_back( region.getName() );
    } );
    auto const key = std::make_pair( meshBodyName, meshLevel.getName());
    meshTargets[key] = std::move( regions );
  } );

  dofManager.addField( wellElementDofName(),
                       FieldLocation::Elem,
                       numDofPerWellElement(),
                       meshTargets );

  dofManager.addCoupling( wellElementDofName(),
                          wellElementDofName(),
                          DofManager::Connector::Node );
}

void WellSolverBase::setPerforationStatus( real64 const & time_n, DomainPartition & domain )
{
  FunctionManager & functionManager = FunctionManager::getInstance();

  // Set well element/perf status
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();
    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   WellElementSubRegion & subRegion )
    {
      WellControls & wellControls = getWellControls( subRegion );

      // Set perforation status

      PerforationData & perforationData = *subRegion.getPerforationData();
      string_array const & perfStatusTableName = perforationData.getPerfStatusTableName();
      arrayView1d< integer > perfStatus = perforationData.getLocalPerfStatus();
      // for now set to open
      for( integer i=0; i<perforationData.size(); i++ )
      {
        TableFunction * tableFunction =  functionManager.getGroupPointer< TableFunction >( perfStatusTableName[i] );
        perfStatus[i]=PerforationData::PerforationStatus::OPEN;
        if( tableFunction->evaluate( &time_n ) < LvArray::NumericLimits< real64 >::epsilon )
        {
          perfStatus[i]=PerforationData::PerforationStatus::CLOSED;
        }
      }

      array1d< localIndex > const perfWellElemIndex = perforationData.getField< fields::perforation::wellElementIndex >();
      // global index local elements (size == subregion.size)
      arrayView1d< globalIndex const > globalWellElementIndex = subRegion.getGlobalWellElementIndex();

      arrayView1d< integer const > const elemGhostRank  = subRegion.ghostRank();
      array1d< integer > & currentStatus = subRegion.getWellElementStatus();
      // Local elements
      array1d< integer > & localElemStatus = subRegion.getWellLocalElementStatus();

      integer numLocalElements = subRegion.getNumLocalElements();
      array1d< integer > segStatus( numLocalElements );

      // Local perforations
      for( integer j = 0; j < perforationData.size(); j++ )
      {
        localIndex const iwelem = perfWellElemIndex[j];
        if( elemGhostRank[iwelem] < 0 )
        {
          if( perfStatus[j] )
          {
            segStatus[iwelem] +=1;
          }
        }
      }
      // Broadcast segment status so all cores have same well status
      subRegion.setElementStatus( segStatus );
      integer numOpenElements = 0;
      array1d< integer > const & updatedStatus = subRegion.getWellElementStatus();
      for( integer i=0; i<currentStatus.size(); i++ )
      {
        numOpenElements += updatedStatus[i];
      }
      numOpenElements>0 ?  wellControls.setWellStatus( time_n, WellControls::Status::OPEN ) :  wellControls.setWellStatus( time_n, WellControls::Status::CLOSED );


      // Set local well element status array
      for( integer i=0; i<subRegion.size(); i++ )
      {
        integer gi = globalWellElementIndex[i];
        localElemStatus[i] = currentStatus[gi];
      }
    } );

  } );
}
void WellSolverBase::implicitStepSetup( real64 const & time_n,
                                        real64 const & GEOS_UNUSED_PARAM( dt ),
                                        DomainPartition & domain )
{

  // Open close perfs
  setPerforationStatus( time_n, domain );
  initializeWells( domain, time_n );

}
void WellSolverBase::setupWellDofs( DomainPartition & domain )
{
  if( m_estimatorDoFManager.empty() )
  {

    map< std::pair< string, string >, string_array > meshTargets;
    string_array regions;
    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const & meshBodyName,
                                                                 MeshLevel & meshLevel,
                                                                 string_array const & regionNames )
    {
      ElementRegionManager & elementRegionManager = meshLevel.getElemManager();
      elementRegionManager.forElementRegions< WellElementRegion >( regionNames,
                                                                   [&]( localIndex const,
                                                                        WellElementRegion & region )
      {
        meshTargets.clear();
        regions.clear();
        regions.emplace_back( region.getName() );
        auto const key = std::make_pair( meshBodyName, meshLevel.getName());
        meshTargets[key] = std::move( regions );

        DofManager regionDoFManager( region.getName());
        regionDoFManager.setDomain( domain );
        regionDoFManager.addField( wellElementDofName(),
                                   FieldLocation::Elem,
                                   numDofPerWellElement(),
                                   meshTargets );

        regionDoFManager.addCoupling( wellElementDofName(),
                                      wellElementDofName(),
                                      DofManager::Connector::Node );

        regionDoFManager.reorderByRank();
        m_estimatorDoFManager.emplace( region.getName(), std::move( regionDoFManager ));
      } );
    } );
  }
}


void WellSolverBase::selectWellConstraint( real64 const & time_n,
                                           real64 const & dt,
                                           const integer coupledIterationNumber,
                                           DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( dt );
  GEOS_UNUSED_VAR( coupledIterationNumber );

  setupWellDofs( domain );
  integer cycleNumber=0;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const & meshBodyName,
                                                               MeshLevel & meshLevel,
                                                               string_array const & regionNames )
  {
    GEOS_UNUSED_VAR( meshBodyName );
    ElementRegionManager & elementRegionManager = meshLevel.getElemManager();
    elementRegionManager.forElementRegions< WellElementRegion >( regionNames,
                                                                 [&]( localIndex const,
                                                                      WellElementRegion & region )
    {
      WellElementSubRegion & subRegion = region.getGroup( ElementRegionBase::viewKeyStruct::elementSubRegions() )
                                           .getGroup< WellElementSubRegion >( region.getSubRegionName() );
      WellControls & wellControls = getWellControls( subRegion );
      if( wellControls.isWellOpen() )
      {
        if( !wellControls.getWellState() )
        {
          wellControls.setWellState( 1 );

          initializeWell( domain, meshLevel, subRegion, time_n );
        }
      }
      else
      {
        wellControls.setWellState( 0 );
      }


      if( wellControls.getWellState())
      {
        //GEOS_LOG_RANK( "**** Estimate Well Solution - Start **** " << subRegion.getName() );
        auto it = m_estimatorDoFManager.find( region.getName());
        if( it == m_estimatorDoFManager.end())
        {
          throw std::runtime_error( "DofManager for region " + region.getName() + " not found." );
        }
        DofManager & dofManager = it->second;

// Only build the sparsity pattern if the mesh has changed
        Timestamp const meshModificationTimestamp = getMeshModificationTimestamp( domain );

        if( meshModificationTimestamp > getSystemSetupTimestamp() )
        {
          // These are esitmator matrices
          setupWellSystem( domain, dofManager, m_localMatrix, m_rhs, m_solution );
          //setSystemSetupTimestamp( meshModificationTimestamp );

          //std::ostringstream oss;
          //m_dofManager.printFieldInfo( oss );
          //GEOS_LOG_LEVEL( logInfo::Fields, oss.str())
        }

        evaluateConstraints( time_n,
                             dt,
                             cycleNumber,
                             coupledIterationNumber,
                             domain,
                             meshLevel,
                             elementRegionManager,
                             subRegion,
                             dofManager );

        // If a well is opened and then timestep is cut resulting in the well being shut, if the well is opened
        // the well initialization code requires control type to by synced
        integer owner = -1;
        // Only subregion owner evaluates well control and control changes need to be broadcast to all ranks
        if( subRegion.isLocallyOwned() )
        {
          owner = MpiWrapper::commRank( MPI_COMM_GEOS );
        }
        owner = MpiWrapper::max( owner );
        WellControls::Control wellControl = wellControls.getControl();
        MpiWrapper::broadcast( wellControl, owner );
        wellControls.setControl( wellControl );
      }
    } );
  } );

}

void WellSolverBase::setupWellSystem( DomainPartition & domain,
                                      DofManager & dofManager,
                                      CRSMatrix< real64, globalIndex > & localMatrix,
                                      ParallelVector & rhs,
                                      ParallelVector & solution,
                                      bool const setSparsity )
{
  GEOS_MARK_FUNCTION;

  setupWellDofs( domain );

  if( setSparsity )
  {
    SparsityPattern< globalIndex > pattern;
    dofManager.setSparsityPattern( pattern );
    localMatrix.assimilate< parallelDevicePolicy<> >( std::move( pattern ) );
  }
  localMatrix.setName( this->getName() + "/matrix" );

  rhs.setName( this->getName() + "/rhs" );
  rhs.create( dofManager.numLocalDofs(), MPI_COMM_GEOS );

  solution.setName( this->getName() + "/solution" );
  solution.create( dofManager.numLocalDofs(), MPI_COMM_GEOS );
}
void WellSolverBase::updateState( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();
    elemManager.forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                WellElementSubRegion & subRegion )
    { updateSubRegionState( elemManager, subRegion ); } );
  } );
}

void WellSolverBase::assembleWellSystem( real64 const time_n,
                                         real64 const dt,
                                         ElementRegionManager const & elementRegionManager,
                                         WellElementSubRegion & subRegion,
                                         DofManager const & dofManager,
                                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                         arrayView1d< real64 > const & localRhs )
{
  assembleWellAccumulationTerms( time_n, dt, subRegion, dofManager, localMatrix.toViewConstSizes(), localRhs );

  assembleWellConstraintTerms( time_n, dt, subRegion, dofManager, localMatrix.toViewConstSizes(), localRhs );

  assembleWellPressureRelations( time_n, dt, subRegion, dofManager, localMatrix.toViewConstSizes(), localRhs );
  computeWellPerforationRates( time_n, dt, elementRegionManager, subRegion );
  assembleWellFluxTerms( time_n, dt, subRegion, dofManager, localMatrix.toViewConstSizes(), localRhs );
  my_ctime=my_ctime+1;

  //  auto iterInfo = currentIter( time_n, dt );
  //  outputWellDebug( time_n, dt, std::get< 0 >( iterInfo ), std::get< 1 >( iterInfo ), std::get< 2 >( iterInfo ),
  //                 domain, dofManager, localMatrix, localRhs );


}

void WellSolverBase::assembleSystem( real64 const time,
                                     real64 const dt,
                                     DomainPartition & domain,
                                     DofManager const & dofManager,
                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                     arrayView1d< real64 > const & localRhs )
{


  // selects constraints one of 2 ways
  //  wellEstimator flag set to 0 => orginal logic rates are computed during update state and constraints are selected every newton
  // iteration
  //  wellEstimator flag > 0 =>   well esitmator solved for each constraint and then selects the constraint
  //                         =>   estimator solve only performed first "wellEstimator" iterations
  NonlinearSolverParameters const & nonlinearParams =  getNonlinearSolverParameters();
  selectWellConstraint( time, dt, nonlinearParams.m_numNewtonIterations, domain );



  // assemble the accumulation term in the mass balance equations
  assembleAccumulationTerms( time, dt, domain, dofManager, localMatrix, localRhs );

  // then assemble the pressure relations between well elements
  assemblePressureRelations( time, dt, domain, dofManager, localMatrix, localRhs );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    ElementRegionManager & elementRegionManager = mesh.getElemManager();
    elementRegionManager.forElementRegions< WellElementRegion >( regionNames,
                                                                 [&]( localIndex const,
                                                                      WellElementRegion & region )
    {
      WellElementSubRegion & subRegion = region.getGroup( ElementRegionBase::viewKeyStruct::elementSubRegions() )
                                           .getGroup< WellElementSubRegion >( region.getSubRegionName() );
        assembleWellConstraintTerms( time, dt, subRegion, dofManager, localMatrix.toViewConstSizes(), localRhs );
    } );
  } );

  // then compute the perforation rates (later assembled by the coupled solver)
  computePerforationRates( time, dt, domain );

  // then assemble the flux terms in the mass balance equations
  // get a reference to the degree-of-freedom numbers
  // then assemble the flux terms in the mass balance equations
  assembleFluxTerms( time, dt, domain, dofManager, localMatrix, localRhs );

  // sort out how this work with well estimator
  auto iterInfo = currentIter( time, dt );
  outputWellDebug( time, dt, std::get< 0 >( iterInfo ), std::get< 1 >( iterInfo ), std::get< 2 >( iterInfo ),
                   domain, dofManager, localMatrix, localRhs );


  my_ctime=my_ctime+1;

}

std::tuple< integer, integer, integer >
WellSolverBase::currentIter( real64 const time, real64 const dt )
{
  if( isEqual( m_currentDt, -1.0 ) )
  {
    m_globalNumTimeSteps=0;
    m_currentTime=time;
    m_prevTime=time;
    m_currentDt=dt;
    m_prevDt=dt;
    m_numTimeStepCuts=0;
    m_currentNewtonIteration=0;
  }
  else
  {
    if( !isEqual( time, m_currentTime ) )
    {
      m_globalNumTimeSteps++;
      m_prevTime=m_currentTime;
      m_prevDt=m_currentDt;
      m_currentTime=time;
      m_currentDt=dt;
      m_currentNewtonIteration=0;
      m_numTimeStepCuts=0;
    }
    else
    {
      if( dt < m_currentDt )
      {
        // timestep cut
        m_globalNumTimeSteps++;
        m_prevTime=m_currentTime;
        m_prevDt=m_currentDt;
        m_currentTime=time;
        m_currentDt=dt;
        m_currentNewtonIteration=0;
        m_numTimeStepCuts++;
        m_currentNewtonIteration=0;
      }
      /*
         else if ( isEqual(dt,m_currentDt ) )
         {
         // next timestep
         m_globalNumTimeSteps++;
         m_prevTime=m_currentTime;
         m_prevDt=m_currentDt;
         m_currentTime=time;
         m_currentDt=dt;
         m_currentNewtonIteration=0;
         m_numTimeStepCuts=0;
         m_currentNewtonIteration=0;
         }*/
      else
      {
        // continuation of current timestep
        m_currentNewtonIteration++;
      }
    }
  }

  return std::tuple< integer, integer, integer >( m_globalNumTimeSteps, m_numTimeStepCuts, m_currentNewtonIteration );

}
void WellSolverBase::initializePostInitialConditionsPreSubGroups()
{
  PhysicsSolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  // make sure that nextWellElementIndex is up-to-date (will be used in well initialization and assembly)
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          WellElementSubRegion & subRegion )
    { subRegion.reconstructLocalConnectivity(); } );
  } );

  // Precompute solver-specific constant data (e.g. gravity-coefficient)
  precomputeData( domain );
}

void WellSolverBase::precomputeData( DomainPartition & domain )
{
  R1Tensor const gravVector = gravityVector();
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          WellElementSubRegion & subRegion )
    {
      PerforationData & perforationData = *subRegion.getPerforationData();
      WellControls & wellControls = getWellControls( subRegion );
      real64 const refElev = wellControls.getReferenceElevation();

      arrayView2d< real64 const > const wellElemLocation = subRegion.getElementCenter();
      arrayView1d< real64 > const wellElemGravCoef = subRegion.getField< well::gravityCoefficient >();

      arrayView2d< real64 const > const perfLocation = perforationData.getField< perforation::location >();
      arrayView1d< real64 > const perfGravCoef = perforationData.getField< well::gravityCoefficient >();

      forAll< serialPolicy >( perforationData.size(), [=]( localIndex const iperf )
      {
        // precompute the depth of the perforations
        perfGravCoef[iperf] = LvArray::tensorOps::AiBi< 3 >( perfLocation[iperf], gravVector );
      } );

      forAll< serialPolicy >( subRegion.size(), [=]( localIndex const iwelem )
      {
        // precompute the depth of the well elements
        wellElemGravCoef[iwelem] = LvArray::tensorOps::AiBi< 3 >( wellElemLocation[iwelem], gravVector );
      } );

      wellControls.forSubGroups< BHPConstraint >( [&]( auto & constraint )
      {
        // set the reference well element where the BHP control is applied
        real64 const refElev1 = constraint.getReferenceElevation();
        constraint.setReferenceGravityCoef( refElev1 * gravVector[2] );
      } );

      // set the reference well element where the BHP control is applied
      wellControls.setReferenceGravityCoef( refElev * gravVector[2] );  // tjb remove
    } );
  } );
}

WellControls & WellSolverBase::getWellControls( WellElementSubRegion const & subRegion )
{
  return this->getGroup< WellControls >( subRegion.getWellControlsName());
}

WellControls const & WellSolverBase::getWellControls( WellElementSubRegion const & subRegion ) const
{
  return this->getGroup< WellControls >( subRegion.getWellControlsName());
}

real64 WellSolverBase::setNextDt( real64 const & currentTime, const real64 & currentDt, geos::DomainPartition & domain )
{
  FunctionManager & functionManager = FunctionManager::getInstance();
  real64 nextDt = PhysicsSolverBase::setNextDt( currentTime, currentDt, domain );

  if( m_nextDt > 0 )
  {
    nextDt = m_nextDt;
    m_nextDt=-1;
  }
  if( m_timeStepFromTables )
  {
    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                 MeshLevel & mesh,
                                                                 string_array const & regionNames )
    {
      mesh.getElemManager().forElementSubRegions< WellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                            WellElementSubRegion & subRegion )
      {
        real64 nextDt_perf=nextDt;
        WellControls & wellControls = getWellControls( subRegion );
        // Find min dt from perf status tables
        PerforationData & perforationData = *subRegion.getPerforationData();
        string_array const & perfStatusTableName = perforationData.getPerfStatusTableName();

        // Get dt for local perforations
        for( integer i=0; i<perforationData.size(); i++ )
        {
          TableFunction * tableFunction =  functionManager.getGroupPointer< TableFunction >( perfStatusTableName[i] );
          WellControls::setNextDtFromTable( tableFunction, currentTime, nextDt_perf );
        }
        nextDt = MpiWrapper::min< real64 >( nextDt_perf );
        // Find min dt including rate and status tables
        real64 const nextDt_orig = nextDt;
        wellControls.setNextDtFromTables( currentTime, nextDt );
        if( m_nonlinearSolverParameters.getLogLevel() > 0 && nextDt < nextDt_orig )
          GEOS_LOG_RANK_0( GEOS_FMT( "{}: next time step based on tables coordinates = {}", getName(), nextDt ));
      } );
    } );
  }

  return nextDt;
}


bool WellSolverBase::solveNonlinearSystem( real64 const & time_n,
                                           real64 const & stepDt,
                                           integer const cycleNumber,
                                           DomainPartition & domain,
                                           MeshLevel & mesh,
                                           ElementRegionManager & elemManager,
                                           WellElementSubRegion & subRegion,
                                           DofManager const & dofManager )
{
  integer const maxNewtonIter = m_nonlinearSolverParameters.m_maxIterNewton;
  integer dtAttempt = m_nonlinearSolverParameters.m_numTimeStepAttempts;
  integer configurationLoopIter = m_nonlinearSolverParameters.m_numConfigurationAttempts;
  integer const minNewtonIter = m_nonlinearSolverParameters.m_minIterNewton;
  real64 const newtonTol = m_nonlinearSolverParameters.m_newtonTol;

// keep residual from previous iteration in case we need to do a line search
  real64 lastResidual = 1e99;
  integer newtonIter = 0;
  real64 scaleFactor = 1.0;
  //m_writeLinearSystem = 2;
  bool isNewtonConverged = false;

  outputSingleWellDebug( time_n, stepDt, 0, 0, 0,
                         mesh, subRegion, dofManager, m_localMatrix.toViewConstSizes(), m_rhs.values()  );
  for( newtonIter = 0; newtonIter < maxNewtonIter; ++newtonIter )
  {
    if( m_nonlinearSolverParameters.getLogLevel() > 4 )
      GEOS_LOG_LEVEL_RANK_0( logInfo::NonlinearSolver,
                             GEOS_FMT( " Well: {}   Est Attempt: NewtonIter: {:2}", subRegion.getName(), stepDt, newtonIter ));

    {
      Timer timer( m_timers.get_inserted( "assemble" ) );

// We sync the nonlinear convergence history. The coupled solver parameters are the one being
// used. We want to propagate the info to subsolvers. It can be important for solvers that
// have special treatment for specific iterations.
      synchronizeNonlinearSolverParameters();

// zero out matrix/rhs before assembly
      m_localMatrix.zero();
      m_rhs.zero();

      arrayView1d< real64 > const localRhs = m_rhs.open();

// call assemble to fill the matrix and the rhs
      assembleWellSystem( time_n,
                          stepDt,
                          elemManager,
                          subRegion,
                          dofManager,
                          m_localMatrix.toViewConstSizes(),
                          localRhs );

// apply boundary conditions to system
      applyWellBoundaryConditions( time_n,
                                   stepDt,
                                   elemManager,
                                   subRegion,
                                   dofManager,
                                   localRhs,
                                   m_localMatrix.toViewConstSizes() );

      m_rhs.close();

      if( m_assemblyCallback )
      {
// Make a copy of LA objects and ship off to the callback
        array1d< real64 > localRhsCopy( m_rhs.localSize() );
        localRhsCopy.setValues< parallelDevicePolicy<> >( m_rhs.values() );
        m_assemblyCallback( m_localMatrix, std::move( localRhsCopy ) );
      }
    }

    ;
    real64 residualNorm = 0;
    {
      Timer timer( m_timers.get_inserted( "convergence check" ) );

// get residual norm
      residualNorm = calculateWellResidualNorm( time_n, stepDt, subRegion, dofManager, m_rhs.values() );
      if( m_nonlinearSolverParameters.getLogLevel() > 4 )
        GEOS_LOG_LEVEL_RANK_0( logInfo::Convergence,
                               GEOS_FMT( "        ( R ) = ( {:4.2e} )", residualNorm ) );
    }
    //auto iterInfo = currentIter( time_n, dt );
    //outputSingleWellDebug( time_n, stepDt, 0, newtonIter, 0,
    //                       mesh, subRegion, dofManager, m_localMatrix.toViewConstSizes(), m_rhs.values()  );
// if the residual norm is less than the Newton tolerance we denote that we have
// converged and break from the Newton loop immediately.
    std::cout << " Well: " << subRegion.getName() << "   Est Attempt: " << dtAttempt
              << ", ConfigurationIter: " << configurationLoopIter
              << ", NewtonIter: " << newtonIter
              << ", Residual Norm: " << residualNorm << std::endl;
    if( residualNorm < newtonTol && newtonIter >= minNewtonIter )
    {
      std::cout << "converged " << std::endl;
      isNewtonConverged = true;
      break;
    }

// if the residual norm is above the max allowed residual norm, we break from
// the Newton loop to avoid crashes due to Newton divergence
    if( residualNorm > m_nonlinearSolverParameters.m_maxAllowedResidualNorm )
    {
      string const maxAllowedResidualNormString = NonlinearSolverParameters::viewKeysStruct::maxAllowedResidualNormString();
      //if( m_nonlinearSolverParameters.getLogLevel() > 4 )
      GEOS_LOG_LEVEL_RANK_0( logInfo::Convergence,
                             GEOS_FMT( "    The residual norm is above the {} of {}. Newton loop terminated.",
                                       maxAllowedResidualNormString,
                                       m_nonlinearSolverParameters.m_maxAllowedResidualNorm )  );
      std::cout << "Residual norm " << residualNorm << " exceeded max allowed " << m_nonlinearSolverParameters.m_maxAllowedResidualNorm << ". Newton loop terminated." << std::endl;
      isNewtonConverged = false;
      break;
    }


    // do line search in case residual has increased

    if( m_nonlinearSolverParameters.m_lineSearchAction != NonlinearSolverParameters::LineSearchAction::None
        && residualNorm > lastResidual * m_nonlinearSolverParameters.m_lineSearchResidualFactor
        && newtonIter >= m_nonlinearSolverParameters.m_lineSearchStartingIteration )
    {
      bool lineSearchSuccess = false;
      if( m_nonlinearSolverParameters.m_lineSearchInterpType == NonlinearSolverParameters::LineSearchInterpolationType::Linear )
      {
        residualNorm = lastResidual;
        lineSearchSuccess = lineSearch1( time_n,
                                         stepDt,
                                         cycleNumber,
                                         domain,
                                         elemManager,
                                         subRegion,
                                         mesh,
                                         dofManager,
                                         m_localMatrix.toViewConstSizes(),
                                         m_rhs,
                                         m_solution,
                                         scaleFactor,
                                         residualNorm );
      }
      else
      {
        lineSearchSuccess = lineSearchWithParabolicInterpolation( time_n,
                                                                  stepDt,
                                                                  cycleNumber,
                                                                  newtonIter,
                                                                  domain,
                                                                  dofManager,
                                                                  m_localMatrix.toViewConstSizes(),
                                                                  m_rhs,
                                                                  m_solution,
                                                                  scaleFactor,
                                                                  lastResidual,
                                                                  residualNorm );
      }

      if( !lineSearchSuccess )
      {
        if( m_nonlinearSolverParameters.m_lineSearchAction == NonlinearSolverParameters::LineSearchAction::Attempt )
        {
          if( m_nonlinearSolverParameters.getLogLevel() > 4 )
            GEOS_LOG_LEVEL_RANK_0( logInfo::LineSearch,
                                   "        Line search failed to produce reduced residual. Accepting iteration." );
        }
        else if( m_nonlinearSolverParameters.m_lineSearchAction == NonlinearSolverParameters::LineSearchAction::Require )
        {
// if line search failed, then break out of the main Newton loop. Timestep will be cut.
          if( m_nonlinearSolverParameters.getLogLevel() > 4 )
            GEOS_LOG_LEVEL_RANK_0( logInfo::LineSearch,
                                   "        Line search failed to produce reduced residual. Exiting Newton Loop." );
          break;
        }
      }
    }

    {
      Timer timer( m_timers.get_inserted( "linear solver total" ) );

// if using adaptive Krylov tolerance scheme, update tolerance.
      LinearSolverParameters::Krylov & krylovParams = m_linearSolverParameters.get().krylov;
      if( krylovParams.useAdaptiveTol )
      {
        krylovParams.relTolerance = newtonIter > 0 ? eisenstatWalker( residualNorm, lastResidual, krylovParams ) : krylovParams.weakestTol;
      }

// TODO: Trilinos currently requires this, re-evaluate after moving to Tpetra-based solvers
      if( m_precond )
      {
        m_precond->clear();
      }

      {
        Timer timer_setup( m_timers.get_inserted( "linear solver create" ) );

// Compose parallel LA matrix/rhs out of local LA matrix/rhs
//
        m_matrix.create( m_localMatrix.toViewConst(), dofManager.numLocalDofs(), MPI_COMM_GEOS );
      }

// Output the linear system matrix/rhs for debugging purposes
      string tag = "_"+std::to_string( my_ctime );
      debugOutputSystem( time_n, cycleNumber, newtonIter, m_matrix, m_rhs, tag );

// Solve the linear system
      try
      {
        // m_writeLinearSystem=2;
        //debugOutputSystem( time_n, cycleNumber, 0, m_matrix, m_rhs );
        solveLinearSystem( dofManager, m_matrix, m_rhs, m_solution );
      } catch( ... )
      {
        // m_writeLinearSystem=2;
        // debugOutputSystem( time_n, cycleNumber, 0, m_matrix, m_rhs );

      }
// Increment the solver statistics for reporting purposes
      getIterationStats().updateNonlinearIteration( m_linearSolverResult.numIterations );

// Output the linear system solution for debugging purposes
      debugOutputSolution( time_n, cycleNumber, newtonIter, m_solution, tag );
    }

    {
      Timer timer( m_timers.get_inserted( "apply solution" ) );

// Compute the scaling factor for the Newton update
      scaleFactor = scalingForWellSystemSolution( subRegion, dofManager, m_solution.values() );
      if( m_nonlinearSolverParameters.getLogLevel() > 4 )
        GEOS_LOG_LEVEL_RANK_0( logInfo::Solution,
                               GEOS_FMT( "        {}: Global solution scaling factor = {}", getName(), scaleFactor ) );

      if( !checkWellSystemSolution( subRegion, dofManager, m_solution.values(), scaleFactor ) )
      {
// TODO try chopping (similar to line search)
        if( m_nonlinearSolverParameters.getLogLevel() > 4 )
          GEOS_LOG_RANK_0( GEOS_FMT( "    {}: Solution check failed. Newton loop terminated.", getName()) );
        break;
      }

// apply the system solution to the fields/variables
      applyWellSystemSolution( dofManager, m_solution.values(), scaleFactor, stepDt, domain, mesh, subRegion );
    }

    {
      Timer timer( m_timers.get_inserted( "update state" ) );

      // update derived variables (constitutive models)
      updateWellState( elemManager, subRegion );
      outputSingleWellDebug( time_n, stepDt, 0, newtonIter+1, 0,
                             mesh, subRegion, dofManager, m_localMatrix.toViewConstSizes(), m_rhs.values()  );
    }

    lastResidual = residualNorm;
  }
  std::cout << "WellSolverBase::solveNewtonSystem completed with isNewtonConverged = " << isNewtonConverged << std::endl;
  return isNewtonConverged;
}

bool WellSolverBase::lineSearch1( real64 const & time_n,
                                  real64 const & dt,
                                  integer const GEOS_UNUSED_PARAM( cycleNumber ),
                                  DomainPartition & domain,
                                  ElementRegionManager & elemManager,
                                  WellElementSubRegion & subRegion,
                                  MeshLevel & mesh,
                                  DofManager const & dofManager,
                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                  ParallelVector & rhs,
                                  ParallelVector & solution,
                                  real64 const scaleFactor,
                                  real64 & lastResidual )
{
  Timer timer( m_timers["line search"] );

  integer const maxNumberLineSearchCuts = m_nonlinearSolverParameters.m_lineSearchMaxCuts;
  real64 const lineSearchCutFactor = m_nonlinearSolverParameters.m_lineSearchCutFactor;

  // flag to determine if we should solve the system and apply the solution. If the line
  // search fails we just bail.
  bool lineSearchSuccess = false;

  real64 residualNorm = lastResidual;

  // scale factor is value applied to the previous solution. In this case we want to
  // subtract a portion of the previous solution.
  real64 localScaleFactor = -scaleFactor;
  real64 cumulativeScale = scaleFactor;

  // main loop for the line search.
  for( integer lineSearchIteration = 0; lineSearchIteration < maxNumberLineSearchCuts; ++lineSearchIteration )
  {
    // cut the scale factor by half. This means that the scale factors will
    // have values of -0.5, -0.25, -0.125, ...
    localScaleFactor *= lineSearchCutFactor;
    cumulativeScale += localScaleFactor;

    if( !checkWellSystemSolution( subRegion, dofManager, m_solution.values(), localScaleFactor ) )
    {
      GEOS_LOG_LEVEL_RANK_0( logInfo::LineSearch,
                             GEOS_FMT( "        Line search {}, solution check failed", lineSearchIteration ) );
      continue;
    }


    applyWellSystemSolution( dofManager, solution.values(), localScaleFactor, dt, domain, mesh, subRegion );
    // update non-primary variables (constitutive models)

    updateWellState( elemManager, subRegion );
    // re-assemble system
    localMatrix.zero();
    rhs.zero();

    arrayView1d< real64 > const localRhs = rhs.open();

    // call assemble to fill the matrix and the rhs
    assembleWellSystem( time_n,
                        dt,
                        elemManager,
                        subRegion,
                        dofManager,
                        localMatrix,
                        localRhs );

// apply boundary conditions to system
    applyWellBoundaryConditions( time_n,
                                 dt,
                                 elemManager,
                                 subRegion,
                                 dofManager,
                                 localRhs,
                                 localMatrix );

    rhs.close();

    GEOS_LOG_LEVEL_RANK_0( logInfo::LineSearch,
                           GEOS_FMT( "        Line search @ {:0.3f}:      ", cumulativeScale ));

    // get residual norm
    residualNorm = calculateWellResidualNorm( time_n, dt, subRegion, dofManager, rhs.values() );
    GEOS_LOG_LEVEL_RANK_0( logInfo::LineSearch,
                           GEOS_FMT( "        ( R ) = ( {:4.2e} )", residualNorm ) );

    // if the residual norm is less than the last residual, we can proceed to the
    // solution step
    if( residualNorm < lastResidual )
    {
      lineSearchSuccess = true;
      break;
    }
  }

  lastResidual = residualNorm;
  return lineSearchSuccess;
}

} // namespace geos
