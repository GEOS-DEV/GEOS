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
 * @file SinglePhaseMixedMFD.cpp
 */

#include "SinglePhaseMixedMFD.hpp"

#include <algorithm>

#include "common/logger/Logger.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "constitutive/permeability/PermeabilityFields.hpp"
#include "constitutive/thermalConductivity/SinglePhaseThermalConductivityBase.hpp"
#include "fieldSpecification/AquiferBoundaryCondition.hpp"
#include "fieldSpecification/FieldSpecificationImpl.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mixedMimetic/MixedMimeticDiscretization.hpp"
#include "mixedMimetic/MixedMimeticDiscretizationManager.hpp"
#include "mixedMimetic/MixedMimeticFields.hpp"
#include "mixedMimetic/MixedMimeticBoundaryConditions.hpp"
#include "mixedMimetic/consistency/ConsistencyAdaptation.hpp"
#include "physicsSolvers/LogLevelsInfo.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/SinglePhaseMixedMFDKernels.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/ThermalSinglePhaseMixedMFDKernels.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/ResidualNormKernel.hpp"

namespace geos
{

using namespace dataRepository;
using namespace constitutive;
using namespace fields;
using namespace mimeticInnerProduct;
using BoundaryType = mixedMimeticBoundary::BoundaryType;
using CellDof = thermalSinglePhaseMixedMFDKernels::CellDof;

SinglePhaseMixedMFD::SinglePhaseMixedMFD( const string & name,
                                          Group * const parent ):
  SinglePhaseBase( name, parent ),
  m_areaRelTol( 1e-8 )
{
  // one cell-centered dof per cell
  m_numDofPerCell = 1;
  m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::singlePhaseMixedMFD;
}

void SinglePhaseMixedMFD::registerDataOnMesh( Group & meshBodies )
{
  SinglePhaseBase::registerDataOnMesh( meshBodies );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    string_array const & regionNames )
  {
    // 1) Register the cell-centered adaptation data
    ElementRegionManager & elemManager = mesh.getElemManager();
    elemManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                              [&]( localIndex const,
                                                                   ElementSubRegionBase & subRegion )
    {
      subRegion.registerField< mixedMimetic::mfdFlag >( getName() );
      subRegion.registerField< mixedMimetic::prescribedMfdFlag >( getName() );
      subRegion.registerField< mixedMimetic::consistencyIndicator >( getName() );
      subRegion.registerField< mixedMimetic::degeneracyIndicator >( getName() );
      if( m_isThermal )
      {
        subRegion.registerField< mixedMimetic::enthalpy >( getName() );
        subRegion.registerField< mixedMimetic::enthalpy_n >( getName() );
      }
    } );

    // 2) Register the face data
    FaceManager & faceManager = mesh.getFaceManager();
    {
      // primary variables: face mass fluxes
      faceManager.registerField< mixedMimetic::faceMassFlux >( getName() );
      faceManager.registerField< mixedMimetic::faceMassFlux_n >( getName() );

      // boundary conditions of the flow operator: the user values and the per-face condition
      faceManager.registerField< flow::bcPressure >( getName() );
      faceManager.registerField< mixedMimetic::bcMassFlux >( getName() );
      faceManager.registerField< mixedMimetic::flowBoundaryType >( getName() );
      faceManager.registerField< mixedMimetic::flowBoundaryValue >( getName() );
      faceManager.registerField< mixedMimetic::flowBoundaryCoefficient >( getName() );

      // face residual of the consistency layer
      faceManager.registerField< mixedMimetic::faceResidual >( getName() );

      // face classification driving the TPFA-face condensation and the MGR labels
      faceManager.registerField< mixedMimetic::faceStencilLabel >( getName() );
      faceManager.registerField< mixedMimetic::faceOrientationCell >( getName() );
      faceManager.registerField< mixedMimetic::faceDofScale >( getName() );

      if( m_isThermal )
      {
        // energy module: face heat flux unknown and the thermal boundary data
        faceManager.registerField< mixedMimetic::faceHeatFlux >( getName() );
        faceManager.registerField< mixedMimetic::faceHeatFlux_n >( getName() );
        faceManager.registerField< flow::bcTemperature >( getName() );
        faceManager.registerField< mixedMimetic::bcEnthalpy >( getName() );
        faceManager.registerField< mixedMimetic::bcHeatFlux >( getName() );
        faceManager.registerField< mixedMimetic::bcHeatTransferCoefficient >( getName() );
        faceManager.registerField< mixedMimetic::bcAmbientTemperature >( getName() );
        faceManager.registerField< mixedMimetic::heatBoundaryType >( getName() );
        faceManager.registerField< mixedMimetic::heatBoundaryValue >( getName() );
        faceManager.registerField< mixedMimetic::heatBoundaryCoefficient >( getName() );
        faceManager.registerField< mixedMimetic::isEnthalpyBcFace >( getName() );
        faceManager.registerField< mixedMimetic::faceHeatDofScale >( getName() );
      }
    }
  } );
}

void SinglePhaseMixedMFD::initializePreSubGroups()
{
  SinglePhaseBase::initializePreSubGroups();

  if( m_isThermal )
  {
    m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::thermalSinglePhaseMixedMFD;
  }

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  MixedMimeticDiscretizationManager const & mmManager = numericalMethodManager.getMixedMimeticDiscretizationManager();

  GEOS_THROW_IF( !mmManager.hasGroup< MixedMimeticDiscretization >( m_discretizationName ),
                 "A MixedMimeticDiscretization must be selected with SinglePhaseMixedMFD",
                 InputError, getDataContext() );
}

void SinglePhaseMixedMFD::initializePostInitialConditionsPreSubGroups()
{
  GEOS_MARK_FUNCTION;

  SinglePhaseBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {
    GEOS_UNUSED_VAR( regionNames );

    ElementRegionManager const & elemManager = mesh.getElemManager();

    // in the kernels, we need to make sure that we act only on the target regions
    for( string const & regionName : regionNames )
    {
      m_regionFilter.insert( elemManager.getRegions().getIndex( regionName ) );
    }

    fsManager.forSubGroups< AquiferBoundaryCondition >( [&] ( AquiferBoundaryCondition const & bc )
    {
      GEOS_UNUSED_VAR( bc );
      GEOS_WARNING( "The aquifer boundary condition was requested in the XML file. \n"
                    "This type of boundary condition is not yet supported by SinglePhaseMixedMFD and will be ignored",
                    getDataContext(), bc.getDataContext() );
    } );
  } );

  markBoundaryFaces( domain );
  if( m_isThermal )
  {
    initializeEnthalpy( domain );
  }

  computeFaceOrientation( domain );

  classifyCells( domain );
}

void SinglePhaseMixedMFD::markBoundaryFaces( DomainPartition & domain )
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {
    FaceManager & faceManager = mesh.getFaceManager();
    arrayView2d< localIndex const > const elemRegionList = faceManager.elementRegionList();
    arrayView2d< localIndex const > const elemSubRegionList = faceManager.elementSubRegionList();
    arrayView2d< localIndex const > const elemList = faceManager.elementList();
    SortedArrayView< localIndex const > const regionFilter = m_regionFilter.toViewConst();

    // the field specifications define the type of the condition on a face; a boundary face without specification
    // has the homogeneous Neumann condition (impervious, adiabatic), an interior face has no condition
    auto const initialize = [&]( arrayView1d< integer > const & type )
    {
      forAll< parallelDevicePolicy<> >( faceManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const kf )
      {
        integer numTargetCells = 0;
        for( integer k = 0; k < elemRegionList.size( 1 ); ++k )
        {
          localIndex const er = elemRegionList[kf][k];
          numTargetCells += ( er >= 0 && elemSubRegionList[kf][k] >= 0 && elemList[kf][k] >= 0 && regionFilter.contains( er ) ) ? 1 : 0;
        }
        type[kf] = ( numTargetCells == 1 ) ? BoundaryType::neumann : BoundaryType::interior;
      } );
    };
    auto const mark = [&]( string const & key, arrayView1d< integer > const & flag, integer const value )
    {
      fsManager.apply< FaceManager >( 0.0, mesh, key,
                                      [&] ( FieldSpecification const &, string const &,
                                            SortedArrayView< localIndex const > const & targetSet, FaceManager &, string const & )
      {
        forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
        {
          flag[targetSet[a]] = value;
        } );
      } );
    };

    arrayView1d< integer > const flowType = faceManager.getField< mixedMimetic::flowBoundaryType >();
    initialize( flowType );
    mark( flow::bcPressure::key(), flowType, BoundaryType::dirichlet );
    mark( mixedMimetic::bcMassFlux::key(), flowType, BoundaryType::neumann );

    if( m_isThermal )
    {
      arrayView1d< integer > const heatType = faceManager.getField< mixedMimetic::heatBoundaryType >();
      initialize( heatType );
      mark( flow::bcTemperature::key(), heatType, BoundaryType::dirichlet );
      mark( mixedMimetic::bcHeatFlux::key(), heatType, BoundaryType::neumann );
      mark( mixedMimetic::bcHeatTransferCoefficient::key(), heatType, BoundaryType::robin );
      mark( mixedMimetic::bcEnthalpy::key(), faceManager.getField< mixedMimetic::isEnthalpyBcFace >(), 1 );
    }
  } );
}

void SinglePhaseMixedMFD::initializeEnthalpy( DomainPartition & domain )
{
  // h = h_eos( p, T ) at the initial state
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                                        [&]( localIndex const,
                                                                             ElementSubRegionBase & subRegion )
    {
      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      SingleFluidBase const & fluid = getConstitutiveModel< SingleFluidBase >( subRegion, fluidName );
      arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > const eosEnthalpy = fluid.enthalpy();
      arrayView1d< real64 > const enthalpy = subRegion.getField< mixedMimetic::enthalpy >();
      arrayView1d< real64 > const enthalpy_n = subRegion.getField< mixedMimetic::enthalpy_n >();
      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
      {
        enthalpy[ei] = eosEnthalpy[ei][0];
        enthalpy_n[ei] = eosEnthalpy[ei][0];
      } );
    } );
  } );
}

void SinglePhaseMixedMFD::computeFaceOrientation( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {
    FaceManager & faceManager = mesh.getFaceManager();
    ElementRegionManager const & elemManager = mesh.getElemManager();
    arrayView2d< localIndex const > const elemRegionList = faceManager.elementRegionList();
    arrayView2d< localIndex const > const elemSubRegionList = faceManager.elementSubRegionList();
    arrayView2d< localIndex const > const elemList = faceManager.elementList();
    ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > const elemLocalToGlobal =
      elemManager.constructArrayViewAccessor< globalIndex, 1 >( ObjectManagerBase::viewKeyStruct::localToGlobalMapString() );
    ElementRegionManager::ElementViewConst< arrayView1d< globalIndex const > > const l2g = elemLocalToGlobal.toNestedViewConst();
    arrayView1d< globalIndex > const orientationCell = faceManager.getField< mixedMimetic::faceOrientationCell >();

    // the owner of a face sees both of its cells; ghost faces take the owner's value through the sync
    forAll< parallelHostPolicy >( faceManager.size(), [=]( localIndex const kf )
    {
      globalIndex gMin = -1;
      for( localIndex k = 0; k < 2; ++k )
      {
        localIndex const er = elemRegionList[kf][k];
        localIndex const esr = elemSubRegionList[kf][k];
        localIndex const ei = elemList[kf][k];
        if( er >= 0 && esr >= 0 && ei >= 0 )
        {
          globalIndex const g = l2g[er][esr][ei];
          gMin = ( gMin < 0 || g < gMin ) ? g : gMin;
        }
      }
      orientationCell[kf] = gMin;
    } );

    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addFields( FieldLocation::Face, { mixedMimetic::faceOrientationCell::key() } );
    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), false );
  } );
}

void SinglePhaseMixedMFD::classifyCells( DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  MixedMimeticDiscretizationManager const & mmManager = numericalMethodManager.getMixedMimeticDiscretizationManager();
  MixedMimeticDiscretization const & discretization = mmManager.getMixedMimeticDiscretization( m_discretizationName );

  ConsistencyAdaptation::Parameters params;
  params.adaptiveConsistency = discretization.isAdaptiveConsistency();
  params.consistencyTolerance = discretization.getConsistencyTolerance();
  params.degeneracyTolerance = discretization.getDegeneracyTolerance();
  params.lengthTolerance = domain.getMeshBody( 0 ).getGlobalLengthScale() * m_areaRelTol;
  params.effectiveTpfa = discretization.isTpfaInnerProduct();
  R1Tensor const gradient = discretization.getNominalGradient();
  for( int d = 0; d < 3; ++d )
  {
    params.nominalGradient[d] = gradient[d];
  }

  ConsistencyAdaptation::Report report;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {
    using PermeabilityAccessors = StencilMaterialAccessors< PermeabilityBase, fields::permeability::permeability >;
    PermeabilityAccessors const permAccessors( mesh.getElemManager(), getName() );
    report.add( ConsistencyAdaptation::classify( mesh, regionNames, m_regionFilter.toViewConst(),
                                                 permAccessors.get( fields::permeability::permeability {} ),
                                                 params, domain.getNeighbors() ) );
  } );

  globalIndex const numCells = MpiWrapper::sum< globalIndex >( report.numCells );
  globalIndex const numConsistent = MpiWrapper::sum< globalIndex >( report.numConsistent );
  globalIndex const numPrescribed0 = MpiWrapper::sum< globalIndex >( report.numPrescribed0 );
  globalIndex const numPrescribed1 = MpiWrapper::sum< globalIndex >( report.numPrescribed1 );
  globalIndex const numDegenerate = MpiWrapper::sum< globalIndex >( report.numDegenerate );
  globalIndex const numPrescribedDegenerate = MpiWrapper::sum< globalIndex >( report.numPrescribedDegenerate );
  globalIndex const numConsistentFinal = MpiWrapper::sum< globalIndex >( report.numConsistentFinal );
  string const consistencyLayer = params.adaptiveConsistency
                                  ? GEOS_FMT( "consistency layer (tolerance = {}) selected {}", params.consistencyTolerance, numConsistent )
                                  : "consistency layer off";
  GEOS_LOG_RANK_0( GEOS_FMT( "mixedMFD Flow: eta = 1 on {} / {} cells: {}, degeneracy layer (tolerance = {} %) switched {} free cells to eta = 0",
                             numConsistentFinal, numCells, consistencyLayer, params.degeneracyTolerance, numDegenerate ) );
  if( numPrescribed0 + numPrescribed1 > 0 )
  {
    GEOS_LOG_RANK_0( GEOS_FMT( "mixedMFD Flow: prescribed eta = 0 on {} cells and eta = 1 on {} cells, kept unchanged ({} of the latter below the degeneracy tolerance)",
                               numPrescribed0, numPrescribed1, numPrescribedDegenerate ) );
    GEOS_WARNING_IF( numPrescribed0 + numPrescribed1 == numCells,
                     GEOS_FMT( "{}: every cell is prescribed (no prescribedMfdFlag = -1): the consistency and degeneracy layers are inert",
                               getDataContext() ) );
  }
}

void SinglePhaseMixedMFD::implicitStepSetup( real64 const & time_n,
                                             real64 const & dt,
                                             DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  // setup the cell-centered fields
  SinglePhaseBase::implicitStepSetup( time_n, dt, domain );

  // setup the face fields
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {
    FaceManager & faceManager = mesh.getFaceManager();

    arrayView1d< real64 const > const & faceFlux =
      faceManager.getField< mixedMimetic::faceMassFlux >();
    arrayView1d< real64 > const & faceFlux_n =
      faceManager.getField< mixedMimetic::faceMassFlux_n >();
    faceFlux_n.setValues< parallelDevicePolicy<> >( faceFlux );

    if( m_isThermal )
    {
      faceManager.getField< mixedMimetic::faceHeatFlux_n >().setValues< parallelDevicePolicy<> >(
        faceManager.getField< mixedMimetic::faceHeatFlux >().toViewConst() );
      mesh.getElemManager().forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                                          [&]( localIndex const,
                                                                               ElementSubRegionBase & subRegion )
      {
        subRegion.getField< mixedMimetic::enthalpy_n >().setValues< parallelDevicePolicy<> >(
          subRegion.getField< mixedMimetic::enthalpy >().toViewConst() );
      } );
    }
  } );

  // evaluate the boundary face values used in the constitutive rows
  applyFaceBoundaryValues( time_n + dt, domain );
}

namespace
{
char const faceBcLogMessage[] =
  "SinglePhaseMixedMFD {}: at time {}s, "
  "the <{}> boundary condition '{}' is applied to the face set '{}' in '{}'. "
  "\nThe total number of target faces (including ghost faces) is {}. "
  "\nNote that if this number is equal to zero, the boundary condition will not be applied on this face set.";
}

void SinglePhaseMixedMFD::applyFaceBoundaryValues( real64 const time,
                                                   DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  stdVector< string > keys = { flow::bcPressure::key(), mixedMimetic::bcMassFlux::key() };
  if( m_isThermal )
  {
    for( string const & key : { flow::bcTemperature::key(), mixedMimetic::bcHeatFlux::key(), mixedMimetic::bcHeatTransferCoefficient::key(),
                                mixedMimetic::bcAmbientTemperature::key(), mixedMimetic::bcEnthalpy::key() } )
    {
      keys.push_back( key );
    }
  }

  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel & mesh,
                                                                      string_array const & )
  {
    // 1) the user values of the face specifications
    for( string const & key : keys )
    {
      fsManager.apply< FaceManager >( time,
                                      mesh,
                                      key,
                                      [&] ( FieldSpecification const & fs,
                                            string const & setName,
                                            SortedArrayView< localIndex const > const & targetSet,
                                            FaceManager & targetGroup,
                                            string const & )
      {
        // provide some logging at the first nonlinear iteration
        if( m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
        {
          globalIndex const numTargetFaces = MpiWrapper::sum< globalIndex >( targetSet.size() );
          GEOS_LOG_LEVEL_RANK_0_ON_GROUP( logInfo::BoundaryConditions,
                                          GEOS_FMT_RUNTIME( faceBcLogMessage,
                                                            this->getName(), time, fs.getCatalogName(), fs.getName(),
                                                            setName, targetGroup.getName(), numTargetFaces ),
                                          fs );
        }
        FieldSpecificationImpl::applyFieldValue< FieldSpecificationEqual,
                                                 parallelDevicePolicy<> >( fs, targetSet, time, targetGroup, key );
      } );
    }

    // 2) the condition of each operator on each face: (g, alpha) from the user values according to the type
    FaceManager & faceManager = mesh.getFaceManager();
    auto const fill = [&]( arrayView1d< integer const > const & type, arrayView1d< real64 > const & value, arrayView1d< real64 > const & coefficient,
                           arrayView1d< real64 const > const & dirichletValue, arrayView1d< real64 const > const & neumannValue,
                           arrayView1d< real64 const > const & robinValue, arrayView1d< real64 const > const & robinCoefficient )
    {
      forAll< parallelDevicePolicy<> >( faceManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const kf )
      {
        value[kf] = ( type[kf] == BoundaryType::dirichlet ) ? dirichletValue[kf]
                    : ( type[kf] == BoundaryType::neumann ) ? neumannValue[kf]
                    : ( type[kf] == BoundaryType::robin ) ? robinValue[kf] : 0.0;
        coefficient[kf] = ( type[kf] == BoundaryType::robin ) ? robinCoefficient[kf] : 0.0;
      } );
    };
    // TODO: Robin data of the flow operator from the Carter-Tracy aquifer. Over a time step the influx rate is
    // q = a - b ( p_f - p_f^n ), with a, b functions of the dimensionless time and of the cumulative influx W^n.
    // With f_A the fraction of the aquifer area on the face, sigma m_f = - rho f_A q = alpha |f| ( p_f - g ),
    // alpha |f| = rho f_A b, g = p_f^n + a / b. ( g, alpha ) is evaluated once per step; at convergence
    // W^{n+1} = W^n + dt q( p_f^{n+1} ). Until then no face of the flow operator has the Robin type
    arrayView1d< real64 const > const zero = faceManager.getField< mixedMimetic::flowBoundaryCoefficient >();
    fill( faceManager.getField< mixedMimetic::flowBoundaryType >(),
          faceManager.getField< mixedMimetic::flowBoundaryValue >(),
          faceManager.getField< mixedMimetic::flowBoundaryCoefficient >(),
          faceManager.getField< flow::bcPressure >(),
          faceManager.getField< mixedMimetic::bcMassFlux >(),
          zero, zero );
    if( m_isThermal )
    {
      fill( faceManager.getField< mixedMimetic::heatBoundaryType >(),
            faceManager.getField< mixedMimetic::heatBoundaryValue >(),
            faceManager.getField< mixedMimetic::heatBoundaryCoefficient >(),
            faceManager.getField< flow::bcTemperature >(),
            faceManager.getField< mixedMimetic::bcHeatFlux >(),
            faceManager.getField< mixedMimetic::bcAmbientTemperature >(),
            faceManager.getField< mixedMimetic::bcHeatTransferCoefficient >() );
    }
  } );
}

void SinglePhaseMixedMFD::setupSystem( DomainPartition & domain,
                                       DofManager & dofManager,
                                       CRSMatrix< real64, globalIndex > & localMatrix,
                                       ParallelVector & rhs,
                                       ParallelVector & solution,
                                       bool const setSparsity )
{
  SinglePhaseBase::setupSystem( domain, dofManager, localMatrix, rhs, solution, setSparsity );

  // with the dof numbering finalized, build the per-dof labels of the MGR reduction
  computeMgrPointMarkers( domain, dofManager );

  // and the sorted list of the ghost dofs, the off-rank columns of the local rows
  computeGhostDofs( domain, dofManager );
}

void SinglePhaseMixedMFD::computeMgrPointMarkers( DomainPartition const & domain,
                                                  DofManager const & dofManager )
{
  GEOS_MARK_FUNCTION;

  string const faceDofKey = dofManager.getKey( mixedMimetic::faceMassFlux::key() );
  string const elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
  globalIndex const rankOffset = dofManager.rankOffset();

  array1d< integer > & pointMarkers = m_linearSolverParameters.get().mgr.customPointMarkers;
  pointMarkers.resize( dofManager.numLocalDofs() );
  arrayView1d< integer > const markers = pointMarkers.toView();

  // flux dofs kept in the saddle point (label 1) and condensed into the pressure system (label 0)
  localIndex numLiveFaces = 0;
  localIndex numCondensedFaces = 0;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               string_array const & )
  {
    FaceManager const & faceManager = mesh.getFaceManager();
    arrayView1d< globalIndex const > const faceDofNumber =
      faceManager.getReference< array1d< globalIndex > >( faceDofKey );
    arrayView1d< integer const > const faceGhostRank = faceManager.ghostRank();
    arrayView1d< integer const > const faceStencilLabel = faceManager.getField< mixedMimetic::faceStencilLabel >();

    RAJA::ReduceSum< parallelHostReduce, localIndex > numLive( 0 );
    RAJA::ReduceSum< parallelHostReduce, localIndex > numCondensed( 0 );
    forAll< parallelHostPolicy >( faceManager.size(), [=]( localIndex const kf )
    {
      bool const owned = faceGhostRank[kf] < 0 && faceDofNumber[kf] >= 0;
      numLive += ( owned && faceStencilLabel[kf] == 1 ) ? 1 : 0;
      numCondensed += ( owned && faceStencilLabel[kf] == 0 ) ? 1 : 0;
    } );
    numLiveFaces += numLive.get();
    numCondensedFaces += numCondensed.get();
  } );
  globalIndex const globalNumLiveFaces = MpiWrapper::sum< globalIndex >( numLiveFaces );
  globalIndex const globalNumCondensedFaces = MpiWrapper::sum< globalIndex >( numCondensedFaces );
  GEOS_LOG_RANK_0( GEOS_FMT( "mixedMFD Flow: flux dofs {} non-condensed (saddle point), {} condensed (two-point closure)",
                             globalNumLiveFaces, globalNumCondensedFaces ) );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               string_array const & regionNames )
  {
    FaceManager const & faceManager = mesh.getFaceManager();
    ElementRegionManager const & elemManager = mesh.getElemManager();

    arrayView1d< globalIndex const > const faceDofNumber =
      faceManager.getReference< array1d< globalIndex > >( faceDofKey );
    arrayView1d< integer const > const faceGhostRank = faceManager.ghostRank();
    arrayView1d< integer const > const faceStencilLabel = faceManager.getField< mixedMimetic::faceStencilLabel >();

    // face-flux dofs: the face classification is the MGR label (0 = exactly-diagonal row)
    forAll< parallelHostPolicy >( faceManager.size(), [=]( localIndex const kf )
    {
      if( faceGhostRank[kf] >= 0 || faceDofNumber[kf] < 0 )
      {
        return;
      }
      markers[faceDofNumber[kf] - rankOffset] = faceStencilLabel[kf];
    } );

    // heat-flux dofs: labels 5 (condensed) and 6 (saddle point)
    if( m_isThermal )
    {
      arrayView1d< globalIndex const > const faceHeatDofNumber =
        faceManager.getReference< array1d< globalIndex > >( dofManager.getKey( mixedMimetic::faceHeatFlux::key() ) );
      forAll< parallelHostPolicy >( faceManager.size(), [=]( localIndex const kf )
      {
        if( faceGhostRank[kf] >= 0 || faceHeatDofNumber[kf] < 0 )
        {
          return;
        }
        markers[faceHeatDofNumber[kf] - rankOffset] = 5 + faceStencilLabel[kf];
      } );
    }

    // cell dofs: pressure 2, temperature 3, enthalpy 4
    integer const numCellDofs = numCellDofComponents();
    elemManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                              [&]( localIndex const,
                                                                   ElementSubRegionBase const & subRegion )
    {
      arrayView1d< globalIndex const > const elemDofNumber =
        subRegion.getReference< array1d< globalIndex > >( elemDofKey );
      arrayView1d< integer const > const elemGhostRank = subRegion.ghostRank();
      forAll< parallelHostPolicy >( subRegion.size(), [=]( localIndex const ei )
      {
        if( elemGhostRank[ei] < 0 )
        {
          for( integer c = 0; c < numCellDofs; ++c )
          {
            markers[elemDofNumber[ei] - rankOffset + c] = 2 + c;
          }
        }
      } );
    } );
  } );
}

void SinglePhaseMixedMFD::setupDofs( DomainPartition const & GEOS_UNUSED_PARAM( domain ),
                                     DofManager & dofManager ) const
{
  // face mass-flux unknowns: the dense NF x NF local mass matrices couple
  // the fluxes of the faces of each cell
  dofManager.addField( mixedMimetic::faceMassFlux::key(),
                       FieldLocation::Face,
                       1,
                       getMeshTargets() );

  dofManager.addCoupling( mixedMimetic::faceMassFlux::key(),
                          mixedMimetic::faceMassFlux::key(),
                          DofManager::Connector::Elem );

  // face heat flux unknowns, with the same structure as the mass fluxes
  if( m_isThermal )
  {
    dofManager.addField( mixedMimetic::faceHeatFlux::key(),
                         FieldLocation::Face,
                         1,
                         getMeshTargets() );
    dofManager.addCoupling( mixedMimetic::faceHeatFlux::key(),
                            mixedMimetic::faceHeatFlux::key(),
                            DofManager::Connector::Elem );
  }

  // cell unknowns (pressure, and temperature and enthalpy when thermal): the TPFA-face
  // condensation writes two-point stencil entries directly into the conservation rows,
  // so cell-to-cell coupling through the faces is required
  dofManager.addField( viewKeyStruct::elemDofFieldString(),
                       FieldLocation::Elem,
                       numCellDofComponents(),
                       getMeshTargets() );

  dofManager.addCoupling( viewKeyStruct::elemDofFieldString(),
                          viewKeyStruct::elemDofFieldString(),
                          DofManager::Connector::Face );

  // coupling between the face fluxes and the cell unknowns
  dofManager.addCoupling( mixedMimetic::faceMassFlux::key(),
                          viewKeyStruct::elemDofFieldString(),
                          DofManager::Connector::Elem );
  if( m_isThermal )
  {
    dofManager.addCoupling( mixedMimetic::faceHeatFlux::key(),
                            viewKeyStruct::elemDofFieldString(),
                            DofManager::Connector::Elem );
  }
}

void SinglePhaseMixedMFD::assembleFluxTerms( real64 const dt,
                                             DomainPartition const & domain,
                                             DofManager const & dofManager,
                                             CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                             arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  MixedMimeticDiscretizationManager const & mmManager = numericalMethodManager.getMixedMimeticDiscretizationManager();
  MixedMimeticDiscretization const & discretization = mmManager.getMixedMimeticDiscretization( m_discretizationName );
  MimeticInnerProductBase const & mimeticInnerProductBase =
    discretization.getReference< MimeticInnerProductBase >( MixedMimeticDiscretization::viewKeyStruct::innerProductString() );

  string const faceDofKey = dofManager.getKey( mixedMimetic::faceMassFlux::key() );
  string const elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  // tolerance for the mass matrix computations
  real64 const lengthTolerance = domain.getMeshBody( 0 ).getGlobalLengthScale() * m_areaRelTol;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               string_array const & regionNames )
  {
    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames,
                                                                        [&]( localIndex const,
                                                                             CellElementSubRegion const & subRegion )
    {
      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      SingleFluidBase const & fluid = getConstitutiveModel< SingleFluidBase >( subRegion, fluidName );

      string const & permName = subRegion.getReference< string >( viewKeyStruct::permeabilityNamesString() );
      PermeabilityBase const & permeability = getConstitutiveModel< PermeabilityBase >( subRegion, permName );

      if( m_isThermal )
      {
        string const & condName = subRegion.getReference< string >( viewKeyStruct::thermalConductivityNamesString() );
        SinglePhaseThermalConductivityBase const & conductivity = getConstitutiveModel< SinglePhaseThermalConductivityBase >( subRegion, condName );
        string const faceHeatDofKey = dofManager.getKey( mixedMimetic::faceHeatFlux::key() );

        thermalSinglePhaseMixedMFDKernels::
          ElementBasedAssemblyKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( dofManager.rankOffset(),
                                                     lengthTolerance,
                                                     elemDofKey,
                                                     faceDofKey,
                                                     faceHeatDofKey,
                                                     nodeManager,
                                                     faceManager,
                                                     mesh.getElemManager(),
                                                     subRegion,
                                                     mimeticInnerProductBase,
                                                     fluid,
                                                     permeability,
                                                     conductivity,
                                                     m_regionFilter.toViewConst(),
                                                     dt,
                                                     localMatrix,
                                                     localRhs );

        // cell closure of the enthalpy unknown
        thermalSinglePhaseMixedMFDKernels::
          EnthalpyClosureKernel::launch< parallelDevicePolicy<> >( dofManager.rankOffset(), elemDofKey, subRegion, fluid, localMatrix, localRhs );
      }
      else
      {
        singlePhaseMixedMFDKernels::
          ElementBasedAssemblyKernelFactory::
          createAndLaunch< parallelDevicePolicy<> >( dofManager.rankOffset(),
                                                     lengthTolerance,
                                                     elemDofKey,
                                                     faceDofKey,
                                                     nodeManager,
                                                     faceManager,
                                                     subRegion,
                                                     mimeticInnerProductBase,
                                                     fluid,
                                                     permeability,
                                                     dt,
                                                     localMatrix,
                                                     localRhs );
      }
    } );

    // condensed (label-0) faces: two-point flux contributions to the conservation
    // rows and one-way closure rows, assembled in a single face-based sweep
    if( m_isThermal )
    {
      thermalSinglePhaseMixedMFDKernels::
        TpfaCondensedFluxKernelFactory::
        createAndLaunch< parallelDevicePolicy<> >( dofManager.rankOffset(),
                                                   lengthTolerance,
                                                   elemDofKey,
                                                   faceDofKey,
                                                   dofManager.getKey( mixedMimetic::faceHeatFlux::key() ),
                                                   getName(),
                                                   nodeManager,
                                                   faceManager,
                                                   mesh.getElemManager(),
                                                   m_regionFilter.toViewConst(),
                                                   dt,
                                                   localMatrix,
                                                   localRhs );
    }
    else
    {
      singlePhaseMixedMFDKernels::
        TpfaCondensedFluxKernelFactory::
        createAndLaunch< parallelDevicePolicy<> >( dofManager.rankOffset(),
                                                   lengthTolerance,
                                                   elemDofKey,
                                                   faceDofKey,
                                                   getName(),
                                                   nodeManager,
                                                   faceManager,
                                                   mesh.getElemManager(),
                                                   m_regionFilter.toViewConst(),
                                                   dt,
                                                   localMatrix,
                                                   localRhs );
    }
  } );
}

void SinglePhaseMixedMFD::assembleStabilizedFluxTerms( real64 const dt,
                                                       DomainPartition const & domain,
                                                       DofManager const & dofManager,
                                                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                       arrayView1d< real64 > const & localRhs )
{
  // pressure stabilization not implemented
  GEOS_UNUSED_VAR( dt, domain, dofManager, localMatrix, localRhs );
  GEOS_ERROR( "Stabilized flux not available for this flow solver" );
}

void SinglePhaseMixedMFD::assembleEDFMFluxTerms( real64 const GEOS_UNUSED_PARAM( time_n ),
                                                 real64 const dt,
                                                 DomainPartition const & domain,
                                                 DofManager const & dofManager,
                                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                 arrayView1d< real64 > const & localRhs,
                                                 string const & jumpDofKey )
{
  GEOS_UNUSED_VAR( jumpDofKey );

  assembleFluxTerms( dt,
                     domain,
                     dofManager,
                     localMatrix,
                     localRhs );
}

void SinglePhaseMixedMFD::applyBoundaryConditions( real64 const time_n,
                                                   real64 const dt,
                                                   DomainPartition & domain,
                                                   DofManager const & dofManager,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                   arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  // the face pressure boundary values are applied inside the constitutive rows during assembly;
  // only the cell-centered boundary conditions (Dirichlet cells, source fluxes) remain to be applied here
  SinglePhaseBase::applyBoundaryConditions( time_n, dt, domain, dofManager, localMatrix, localRhs );

  if( m_isThermal )
  {
    removeEnthalpyReference( domain, dofManager, localMatrix, localRhs );
  }

  // the matrix is final: the residual-norm weights are a function of it and are recomputed with it
  computeResidualWeights( domain, dofManager, localMatrix.toViewConst() );
}

void SinglePhaseMixedMFD::removeEnthalpyReference( DomainPartition & domain,
                                                   DofManager const & dofManager,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                   arrayView1d< real64 > const & localRhs ) const
{
  GEOS_MARK_FUNCTION;

  // energy row <- energy row - h_K * mass row. A shift h_0 of the enthalpy reference adds h_0 * (mass row) to the
  // energy row: the reduced system is invariant under that shift, and the Newton update is unchanged
  globalIndex const rankOffset = dofManager.rankOffset();
  string const elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                                        [&]( localIndex const,
                                                                             ElementSubRegionBase const & subRegion )
    {
      arrayView1d< globalIndex const > const elemDofNumber = subRegion.getReference< array1d< globalIndex > >( elemDofKey );
      arrayView1d< integer const > const elemGhostRank = subRegion.ghostRank();
      arrayView1d< real64 const > const enthalpy = subRegion.getField< mixedMimetic::enthalpy >();

      RAJA::ReduceSum< parallelDeviceReduce, localIndex > numMissing( 0 );
      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
      {
        if( elemGhostRank[ei] >= 0 )
        {
          return;
        }
        localIndex const massRow = elemDofNumber[ei] - rankOffset + CellDof::pressure;
        localIndex const energyRow = elemDofNumber[ei] - rankOffset + CellDof::temperature;
        real64 const h = enthalpy[ei];

        arraySlice1d< globalIndex const > const massCols = localMatrix.getColumns( massRow );
        arraySlice1d< real64 const > const massVals = localMatrix.getEntries( massRow );
        arraySlice1d< globalIndex const > const energyCols = localMatrix.getColumns( energyRow );
        arraySlice1d< real64 > const energyVals = localMatrix.getEntries( energyRow );
        localIndex const numEnergy = localMatrix.numNonZeros( energyRow );

        // both rows are sorted by column
        localIndex j = 0;
        for( localIndex k = 0; k < localMatrix.numNonZeros( massRow ); ++k )
        {
          while( j < numEnergy && energyCols[j] < massCols[k] )
          {
            ++j;
          }
          if( j < numEnergy && energyCols[j] == massCols[k] )
          {
            energyVals[j] -= h * massVals[k];
          }
          else if( LvArray::math::abs( massVals[k] ) > 0.0 )
          {
            numMissing += 1;
          }
        }
        localRhs[energyRow] -= h * localRhs[massRow];
      } );
      GEOS_ERROR_IF( numMissing.get() > 0, getDataContext() << ": the energy rows do not contain the pattern of the mass rows" );
    } );
  } );
}

void SinglePhaseMixedMFD::computeGhostDofs( DomainPartition const & domain,
                                            DofManager const & dofManager )
{
  GEOS_MARK_FUNCTION;

  string const elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
  integer const numCellDofs = numCellDofComponents();
  stdVector< string > faceDofKeys = { dofManager.getKey( mixedMimetic::faceMassFlux::key() ) };
  if( m_isThermal )
  {
    faceDofKeys.push_back( dofManager.getKey( mixedMimetic::faceHeatFlux::key() ) );
  }

  array1d< globalIndex > dofs;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               string_array const & regionNames )
  {
    FaceManager const & faceManager = mesh.getFaceManager();
    arrayView1d< integer const > const faceGhostRank = faceManager.ghostRank();
    for( string const & key : faceDofKeys )
    {
      arrayView1d< globalIndex const > const faceDofNumber = faceManager.getReference< array1d< globalIndex > >( key );
      for( localIndex kf = 0; kf < faceManager.size(); ++kf )
      {
        if( faceGhostRank[kf] >= 0 && faceDofNumber[kf] >= 0 )
        {
          dofs.emplace_back( faceDofNumber[kf] );
        }
      }
    }
    mesh.getElemManager().forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                                        [&]( localIndex const,
                                                                             ElementSubRegionBase const & subRegion )
    {
      arrayView1d< globalIndex const > const elemDofNumber = subRegion.getReference< array1d< globalIndex > >( elemDofKey );
      arrayView1d< integer const > const elemGhostRank = subRegion.ghostRank();
      for( localIndex ei = 0; ei < subRegion.size(); ++ei )
      {
        if( elemGhostRank[ei] >= 0 && elemDofNumber[ei] >= 0 )
        {
          for( integer c = 0; c < numCellDofs; ++c )
          {
            dofs.emplace_back( elemDofNumber[ei] + c );
          }
        }
      }
    } );
  } );

  std::sort( dofs.begin(), dofs.end() );
  m_ghostDofs.clear();
  m_ghostDofs.insert( dofs.begin(), std::unique( dofs.begin(), dofs.end() ) );
}

void SinglePhaseMixedMFD::computeResidualWeights( DomainPartition & domain,
                                                  DofManager const & dofManager,
                                                  CRSMatrixView< real64 const, globalIndex const > const & localMatrix )
{
  GEOS_MARK_FUNCTION;

  // residual-norm weights w_i = ||A_i S||_inf = max_j |A_ij| s_j, S = diag(s_j) the
  // characteristic scales of the unknowns: each term |A_ij| s_j has the units of r_i,
  // and w_i is invariant under row rescaling and under the TPFA-face condensation
  localIndex const numRows = localMatrix.numRows();
  m_residualWeight.resize( numRows );
  m_dofScale.resize( numRows );
  arrayView1d< real64 > const residualWeight = m_residualWeight.toView();
  arrayView1d< real64 > const dofScale = m_dofScale.toView();
  residualWeight.zero();
  dofScale.zero();

  globalIndex const rankOffset = dofManager.rankOffset();
  string const elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
  integer const numCellDofs = numCellDofComponents();

  // the cell scales |p_n|, and |T_n|, |h_n| when thermal; the face fields with the cell scale
  // their rows are measured with and the face field storing their scale
  stdVector< string > cellScaleKeys = { flow::pressure_n::key() };
  struct FaceField { string dofKey; string cellScaleKey; string scaleFieldKey; };
  stdVector< FaceField > faceFields =
  { { dofManager.getKey( mixedMimetic::faceMassFlux::key() ), flow::pressure_n::key(), mixedMimetic::faceDofScale::key() } };
  if( m_isThermal )
  {
    cellScaleKeys.push_back( flow::temperature_n::key() );
    cellScaleKeys.push_back( mixedMimetic::enthalpy_n::key() );
    faceFields.push_back( { dofManager.getKey( mixedMimetic::faceHeatFlux::key() ), flow::temperature_n::key(), mixedMimetic::faceHeatDofScale::key() } );
  }

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();
    FaceManager & faceManager = mesh.getFaceManager();

    // s of the cell dofs
    elemManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                              [&]( localIndex const,
                                                                   ElementSubRegionBase const & subRegion )
    {
      arrayView1d< globalIndex const > const elemDofNumber =
        subRegion.getReference< array1d< globalIndex > >( elemDofKey );
      arrayView1d< integer const > const elemGhostRank = subRegion.ghostRank();
      for( integer c = 0; c < numCellDofs; ++c )
      {
        arrayView1d< real64 const > const valueN = subRegion.getReference< array1d< real64 > >( cellScaleKeys[c] );
        forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
        {
          if( elemGhostRank[ei] < 0 )
          {
            dofScale[elemDofNumber[ei] - rankOffset + c] = LvArray::math::abs( valueN[ei] );
          }
        } );
      }
    } );

    // s_f = x_scale / |M_ff| for the face dofs, x_scale the mean adjacent cell scale
    arrayView1d< integer const > const faceGhostRank = faceManager.ghostRank();
    arrayView2d< localIndex const > const elemRegionList = faceManager.elementRegionList();
    arrayView2d< localIndex const > const elemSubRegionList = faceManager.elementSubRegionList();
    arrayView2d< localIndex const > const elemList = faceManager.elementList();
    SortedArrayView< localIndex const > const regionFilter = m_regionFilter.toViewConst();
    for( FaceField const & ff : faceFields )
    {
      ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const valueNAccessor =
        elemManager.constructArrayViewAccessor< real64, 1 >( ff.cellScaleKey );
      ElementRegionManager::ElementViewConst< arrayView1d< real64 const > > const valueN = valueNAccessor.toNestedViewConst();
      arrayView1d< globalIndex const > const faceDofNumber = faceManager.getReference< array1d< globalIndex > >( ff.dofKey );
      arrayView1d< real64 > const faceScale = faceManager.getReference< array1d< real64 > >( ff.scaleFieldKey );

      forAll< parallelDevicePolicy<> >( faceManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const kf )
      {
        if( faceGhostRank[kf] >= 0 || faceDofNumber[kf] < 0 )
        {
          return;
        }
        real64 sum = 0.0;
        integer elemCounter = 0;
        for( integer k = 0; k < elemRegionList.size( 1 ); ++k )
        {
          localIndex const er  = elemRegionList[kf][k];
          localIndex const esr = elemSubRegionList[kf][k];
          localIndex const ei  = elemList[kf][k];
          if( er >= 0 && esr >= 0 && ei >= 0 && regionFilter.contains( er ) )
          {
            sum += LvArray::math::abs( valueN[er][esr][ei] );
            elemCounter++;
          }
        }
        real64 const xScale = elemCounter > 0 ? sum / elemCounter : 0.0;
        localIndex const localRow = faceDofNumber[kf] - rankOffset;
        real64 diag = 0.0;
        arraySlice1d< globalIndex const > const columns = localMatrix.getColumns( localRow );
        arraySlice1d< real64 const > const entries = localMatrix.getEntries( localRow );
        for( localIndex k = 0; k < localMatrix.numNonZeros( localRow ); ++k )
        {
          if( columns[k] == faceDofNumber[kf] )
          {
            diag = LvArray::math::abs( entries[k] );
          }
        }
        dofScale[localRow] = diag > 0.0 ? xScale / diag : xScale;
        faceScale[kf] = dofScale[localRow];
      } );
    }
  } );

  // s_j of the ghost dofs: a cell may own none of its faces, so the face scales come from their owners;
  // the scale of a ghost cell dof is its synchronized previous value
  m_ghostDofScale.resize( m_ghostDofs.size() );
  arrayView1d< real64 > const ghostDofScale = m_ghostDofScale.toView();
  ghostDofScale.zero();
  SortedArrayView< globalIndex const > const ghostDofs = m_ghostDofs.toViewConst();

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & regionNames )
  {
    FieldIdentifiers fieldsToBeSync;
    for( FaceField const & ff : faceFields )
    {
      fieldsToBeSync.addFields( FieldLocation::Face, { ff.scaleFieldKey } );
    }
    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), false );

    FaceManager const & faceManager = mesh.getFaceManager();
    arrayView1d< integer const > const faceGhostRank = faceManager.ghostRank();
    for( FaceField const & ff : faceFields )
    {
      arrayView1d< globalIndex const > const faceDofNumber = faceManager.getReference< array1d< globalIndex > >( ff.dofKey );
      arrayView1d< real64 const > const faceScale = faceManager.getReference< array1d< real64 > >( ff.scaleFieldKey );
      forAll< parallelDevicePolicy<> >( faceManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const kf )
      {
        if( faceGhostRank[kf] >= 0 && faceDofNumber[kf] >= 0 )
        {
          localIndex const slot = LvArray::sortedArrayManipulation::find( ghostDofs.begin(), ghostDofs.size(), faceDofNumber[kf] );
          ghostDofScale[slot] = faceScale[kf];
        }
      } );
    }

    mesh.getElemManager().forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                                        [&]( localIndex const,
                                                                             ElementSubRegionBase const & subRegion )
    {
      arrayView1d< globalIndex const > const elemDofNumber = subRegion.getReference< array1d< globalIndex > >( elemDofKey );
      arrayView1d< integer const > const elemGhostRank = subRegion.ghostRank();
      for( integer c = 0; c < numCellDofs; ++c )
      {
        arrayView1d< real64 const > const valueN = subRegion.getReference< array1d< real64 > >( cellScaleKeys[c] );
        forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
        {
          if( elemGhostRank[ei] >= 0 && elemDofNumber[ei] >= 0 )
          {
            localIndex const slot = LvArray::sortedArrayManipulation::find( ghostDofs.begin(), ghostDofs.size(), elemDofNumber[ei] + c );
            ghostDofScale[slot] = LvArray::math::abs( valueN[ei] );
          }
        } );
      }
    } );
  } );

  // w_i = max_j |A_ij| s_j over every column of row i
  arrayView1d< real64 const > const ghostScale = m_ghostDofScale.toViewConst();
  arrayView1d< real64 const > const localScale = m_dofScale.toViewConst();
  forAll< parallelDevicePolicy<> >( numRows, [=] GEOS_HOST_DEVICE ( localIndex const i )
  {
    real64 w = 0.0;
    arraySlice1d< globalIndex const > const columns = localMatrix.getColumns( i );
    arraySlice1d< real64 const > const entries = localMatrix.getEntries( i );
    for( localIndex k = 0; k < localMatrix.numNonZeros( i ); ++k )
    {
      globalIndex const localCol = columns[k] - rankOffset;
      real64 scale = 0.0;
      if( localCol >= 0 && localCol < numRows )
      {
        scale = localScale[localCol];
      }
      else
      {
        localIndex const slot = LvArray::sortedArrayManipulation::find( ghostDofs.begin(), ghostDofs.size(), columns[k] );
        scale = ( slot < ghostDofs.size() && ghostDofs[slot] == columns[k] ) ? ghostScale[slot] : 0.0;
      }
      w = LvArray::math::max( w, LvArray::math::abs( entries[k] ) * scale );
    }
    residualWeight[i] = w;
  } );
}

void SinglePhaseMixedMFD::applyAquiferBC( real64 const time,
                                          real64 const dt,
                                          DomainPartition & domain,
                                          DofManager const & dofManager,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs ) const
{
  GEOS_MARK_FUNCTION;
  // no contribution: the aquifer is a Robin condition of the flow operator (see applyFaceBoundaryValues)
  GEOS_UNUSED_VAR( time, dt, dofManager, domain, localMatrix, localRhs );
}

void SinglePhaseMixedMFD::saveAquiferConvergedState( real64 const & time,
                                                     real64 const & dt,
                                                     DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;
  // TODO: W^{n+1} = W^n + dt q( p_f^{n+1} ) on the faces of each aquifer (see applyFaceBoundaryValues)
  GEOS_UNUSED_VAR( time, dt, domain );
}

real64 SinglePhaseMixedMFD::calculateResidualNorm( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                                   real64 const & GEOS_UNUSED_PARAM( dt ),
                                                   DomainPartition const & GEOS_UNUSED_PARAM( domain ),
                                                   DofManager const & GEOS_UNUSED_PARAM( dofManager ),
                                                   arrayView1d< real64 const > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  // row-equilibrated residual norm ||D^{-1} r||, D = diag(max(eps, w_i)), w_i = ||A_i S||_inf
  physicsSolverBaseKernels::NormType const normType = getNonlinearSolverParameters().normType();
  real64 const minNormalizer = m_nonlinearSolverParameters.m_minNormalizer;

  arrayView1d< real64 const > const residualWeight = m_residualWeight.toViewConst();
  localIndex const numRows = localRhs.size();

  real64 residualNorm = 0.0;
  if( normType == physicsSolverBaseKernels::NormType::Linf )
  {
    RAJA::ReduceMax< parallelDeviceReduce, real64 > maxVal( 0.0 );
    forAll< parallelDevicePolicy<> >( numRows, [=] GEOS_HOST_DEVICE ( localIndex const i )
    {
      maxVal.max( LvArray::math::abs( localRhs[i] ) /
                  LvArray::math::max( minNormalizer, residualWeight[i] ) );
    } );
    real64 const localNorm = maxVal.get();
    physicsSolverBaseKernels::LinfResidualNormHelper::computeGlobalNorm( localNorm, residualNorm );
  }
  else
  {
    RAJA::ReduceSum< parallelDeviceReduce, real64 > sumVal( 0.0 );
    forAll< parallelDevicePolicy<> >( numRows, [=] GEOS_HOST_DEVICE ( localIndex const i )
    {
      real64 const r = localRhs[i] / LvArray::math::max( minNormalizer, residualWeight[i] );
      sumVal += r * r;
    } );
    real64 const localSum = sumVal.get();
    real64 const localCount = static_cast< real64 >( numRows );
    physicsSolverBaseKernels::L2ResidualNormHelper::computeGlobalNorm( localSum, localCount, residualNorm );
  }

  GEOS_LOG_LEVEL_RANK_0_NLR( logInfo::ResidualNorm,
                             GEOS_FMT( "        ( R{} ) = ( {:4.2e} )", coupledSolverAttributePrefix(), residualNorm ));
  getConvergenceStats().setResidualValue( GEOS_FMT( "R{}", coupledSolverAttributePrefix()), residualNorm );

  return residualNorm;
}

void SinglePhaseMixedMFD::applySystemSolution( DofManager const & dofManager,
                                               arrayView1d< real64 const > const & localSolution,
                                               real64 const scalingFactor,
                                               real64 const dt,
                                               DomainPartition & domain )
{
  GEOS_UNUSED_VAR( dt );

  // 1. apply the cell-centered update

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::elemDofFieldString(),
                               flow::pressure::key(),
                               scalingFactor,
                               DofManager::CompMask( numCellDofComponents(), CellDof::pressure, CellDof::pressure + 1 ) );

  // 2. apply the face-based update

  dofManager.addVectorToField( localSolution,
                               mixedMimetic::faceMassFlux::key(),
                               mixedMimetic::faceMassFlux::key(),
                               scalingFactor );

  if( m_isThermal )
  {
    dofManager.addVectorToField( localSolution,
                                 viewKeyStruct::elemDofFieldString(),
                                 flow::temperature::key(),
                                 scalingFactor,
                                 DofManager::CompMask( CellDof::num, CellDof::temperature, CellDof::temperature + 1 ) );
    dofManager.addVectorToField( localSolution,
                                 viewKeyStruct::elemDofFieldString(),
                                 mixedMimetic::enthalpy::key(),
                                 scalingFactor,
                                 DofManager::CompMask( CellDof::num, CellDof::enthalpy, CellDof::enthalpy + 1 ) );
    dofManager.addVectorToField( localSolution,
                                 mixedMimetic::faceHeatFlux::key(),
                                 mixedMimetic::faceHeatFlux::key(),
                                 scalingFactor );
  }

  // 3. synchronize
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {
    FieldIdentifiers fieldsToBeSync;

    fieldsToBeSync.addElementFields( { flow::pressure::key() }, regionNames );
    fieldsToBeSync.addFields( FieldLocation::Face, { mixedMimetic::faceMassFlux::key() } );
    if( m_isThermal )
    {
      fieldsToBeSync.addElementFields( { flow::temperature::key(), mixedMimetic::enthalpy::key() }, regionNames );
      fieldsToBeSync.addFields( FieldLocation::Face, { mixedMimetic::faceHeatFlux::key() } );
    }

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), true );
  } );
}

void SinglePhaseMixedMFD::resetStateToBeginningOfStep( DomainPartition & domain )
{
  // Reset the cell-centered fields
  SinglePhaseBase::resetStateToBeginningOfStep( domain );

  // Reset the face-based fields
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {
    FaceManager & faceManager = mesh.getFaceManager();

    arrayView1d< real64 > const & faceFlux =
      faceManager.getField< mixedMimetic::faceMassFlux >();
    arrayView1d< real64 const > const & faceFlux_n =
      faceManager.getField< mixedMimetic::faceMassFlux_n >();
    faceFlux.setValues< parallelDevicePolicy<> >( faceFlux_n );

    if( m_isThermal )
    {
      faceManager.getField< mixedMimetic::faceHeatFlux >().setValues< parallelDevicePolicy<> >(
        faceManager.getField< mixedMimetic::faceHeatFlux_n >().toViewConst() );
      mesh.getElemManager().forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                                          [&]( localIndex const,
                                                                               ElementSubRegionBase & subRegion )
      {
        subRegion.getField< mixedMimetic::enthalpy >().setValues< parallelDevicePolicy<> >(
          subRegion.getField< mixedMimetic::enthalpy_n >().toViewConst() );
      } );
    }
  } );
}

REGISTER_CATALOG_ENTRY( PhysicsSolverBase, SinglePhaseMixedMFD, string const &, Group * const )
} /* namespace geos */
