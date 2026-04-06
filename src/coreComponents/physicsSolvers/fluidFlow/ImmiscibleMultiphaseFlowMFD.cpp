/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * (c) GEOS/GEOSX Contributors
 * ------------------------------------------------------------------------------------------------------------
 */

#include "ImmiscibleMultiphaseFlowMFD.hpp"

#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "physicsSolvers/fluidFlow/ImmiscibleMultiphaseFlowFields.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseUtilities.hpp"
#include "physicsSolvers/fluidFlow/kernels/immiscibleMultiphase/ImmiscibleMultiphaseMFDKernels.hpp"
#include "physicsSolvers/fluidFlow/kernels/immiscibleMultiphase/ImmiscibleMultiphaseKernels.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/RelativePermeabilityUpdateKernel.hpp"
#include "physicsSolvers/fluidFlow/kernels/compositional/CapillaryPressureUpdateKernel.hpp"
#include "physicsSolvers/PhysicsSolverBaseKernels.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/HybridMimeticDiscretization.hpp"
#include "mesh/DomainPartition.hpp"
// Add explicit include for FluxApproximationBase to use hasGroup with this type
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fieldSpecification/SourceFluxBoundaryCondition.hpp"
#include "fieldSpecification/EquilibriumInitialCondition.hpp"
#include "constitutive/fluid/twophaseimmisciblefluid/TwoPhaseImmiscibleFluid.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "constitutive/permeability/PermeabilityBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "constitutive/capillaryPressure/CapillaryPressureBase.hpp"
#include "constitutive/solid/porosity/PorosityFields.hpp"
#include "physicsSolvers/fluidFlow/SourceFluxStatistics.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "constitutive/relativePermeability/RelativePermeabilitySelector.hpp"
#include "constitutive/capillaryPressure/CapillaryPressureSelector.hpp"

namespace geos
{

using namespace dataRepository;
using namespace fields;
using namespace immiscibleMultiphaseMFDKernels;
using namespace immiscibleMultiphaseKernels;
using namespace constitutive;
using namespace isothermalCompositionalMultiphaseBaseKernels; // for relperm / capillary update kernels

// --- Small local helpers implementations to reduce duplication ---
void ImmiscibleMultiphaseFlowMFD::copyPressureToPrev( ElementSubRegionBase & subRegion ) const
{
  arrayView1d< real64 > const pres_n = subRegion.getField< flow::pressure_n >();
  arrayView1d< real64 const > const pres = subRegion.getField< flow::pressure >();
  pres_n.template setValues< parallelDevicePolicy<> >( pres );
}

void ImmiscibleMultiphaseFlowMFD::copyPhaseStateToPrev( ElementSubRegionBase & subRegion ) const
{
  arrayView2d< real64, immiscibleFlow::USD_PHASE > const phaseVol = subRegion.getField< immiscibleMultiphaseFlow::phaseVolumeFraction >();
  arrayView2d< real64, immiscibleFlow::USD_PHASE > phaseVol_n = subRegion.getField< immiscibleMultiphaseFlow::phaseVolumeFraction_n >();
  phaseVol_n.template setValues< parallelDevicePolicy<> >( phaseVol );
  arrayView2d< real64, immiscibleFlow::USD_PHASE > const phaseMass = subRegion.getField< immiscibleMultiphaseFlow::phaseMass >();
  arrayView2d< real64, immiscibleFlow::USD_PHASE > phaseMass_n = subRegion.getField< immiscibleMultiphaseFlow::phaseMass_n >();
  phaseMass_n.template setValues< parallelDevicePolicy<> >( phaseMass );
}

void ImmiscibleMultiphaseFlowMFD::restoreStateFromPrev( ElementSubRegionBase & subRegion ) const
{
  arrayView1d< real64 > pres = subRegion.getField< flow::pressure >();
  arrayView1d< real64 const > const pres_n = subRegion.getField< flow::pressure_n >();
  pres.template setValues< parallelDevicePolicy<> >( pres_n );
  arrayView2d< real64, immiscibleFlow::USD_PHASE > phaseVol = subRegion.getField< immiscibleMultiphaseFlow::phaseVolumeFraction >();
  arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const phaseVol_n = subRegion.getField< immiscibleMultiphaseFlow::phaseVolumeFraction_n >();
  phaseVol.template setValues< parallelDevicePolicy<> >( phaseVol_n );
  arrayView2d< real64, immiscibleFlow::USD_PHASE > phaseMass = subRegion.getField< immiscibleMultiphaseFlow::phaseMass >();
  arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const phaseMass_n = subRegion.getField< immiscibleMultiphaseFlow::phaseMass_n >();
  phaseMass.template setValues< parallelDevicePolicy<> >( phaseMass_n );
}

void ImmiscibleMultiphaseFlowMFD::copyFacePressureToPrev( MeshLevel & mesh ) const
{
  FaceManager & fm = mesh.getFaceManager();
  arrayView1d< real64 const > const faceP = fm.getField< flow::facePressure >();
  arrayView1d< real64 > faceP_n = fm.getField< flow::facePressure_n >();
  faceP_n.setValues< parallelDevicePolicy<> >( faceP );
}

void ImmiscibleMultiphaseFlowMFD::restoreFacePressureFromPrev( MeshLevel & mesh ) const
{
  FaceManager & fm = mesh.getFaceManager();
  arrayView1d< real64 > faceP = fm.getField< flow::facePressure >();
  arrayView1d< real64 const > const faceP_n = fm.getField< flow::facePressure_n >();
  faceP.setValues< parallelDevicePolicy<> >( faceP_n );
}

ImmiscibleMultiphaseFlowMFD::ImmiscibleMultiphaseFlowMFD( const string & name,
                                                          Group * const parent )
  : FlowSolverBase( name, parent ),
    m_numPhases( 2 ),
    m_hasCapPressure( false ),
    m_useTotalMassEquation( 1 ),
    m_targetRelativePresChange( 0.2 ),
    m_targetPhaseVolFracChange( 0.2 ),
    m_solutionChangeScalingFactor( 0.5 ),
    m_gravityDensityScheme( GravityDensityScheme::ArithmeticAverage ),
    m_areaRelTol( 1e-8 )
{
  // basic wrappers (subset of ImmiscibleMultiphaseFlow for now)
  this->registerWrapper( viewKeyStruct::useTotalMassEquationString(), &m_useTotalMassEquation )
      .setSizedFromParent( 0 ).setInputFlag( InputFlags::OPTIONAL ).setApplyDefaultValue( 1 )
      .setDescription( "Flag indicating whether total mass equation is used" );
  this->registerWrapper( "gravityDensityScheme", &m_gravityDensityScheme )
      .setSizedFromParent( 0 ).setInputFlag( InputFlags::OPTIONAL )
      .setApplyDefaultValue( GravityDensityScheme::ArithmeticAverage )
      .setDescription( "Scheme for density treatment in gravity" );
  this->registerWrapper( "solutionChangeScalingFactor", &m_solutionChangeScalingFactor )
      .setSizedFromParent( 0 ).setInputFlag( InputFlags::OPTIONAL ).setApplyDefaultValue( 0.5 );
  this->registerWrapper( "targetRelativePressureChangeInTimeStep", &m_targetRelativePresChange )
      .setSizedFromParent( 0 ).setInputFlag( InputFlags::OPTIONAL ).setApplyDefaultValue( 0.2 );
  this->registerWrapper( "targetPhaseVolFractionChangeInTimeStep", &m_targetPhaseVolFracChange )
      .setSizedFromParent( 0 ).setInputFlag( InputFlags::OPTIONAL ).setApplyDefaultValue( 0.2 );

  // expose dependent phase index (0-based) to XML; for two-phase, 0=phase 0, 1=phase 1
  this->registerWrapper( viewKeyStruct::dependentPhaseIndexString(), &m_dependentPhaseIndex )
      .setSizedFromParent( 0 )
      .setInputFlag( InputFlags::OPTIONAL )
      .setApplyDefaultValue( 1 )
      .setDescription( "Index of the dependent phase saturation (0-based). For two-phase, 0 or 1; s_dep = 1 - s_ind." );

  // number of dofs per cell = num phases (pressure + (numPhases-1) saturations OR total+components)
  m_numDofPerCell = m_numPhases;
  m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::immiscibleMultiphaseMFD;
}

void ImmiscibleMultiphaseFlowMFD::setConstitutiveNames( ElementSubRegionBase & subRegion ) const
{
  setConstitutiveNamesCallSuper( subRegion );
  setConstitutiveName< TwoPhaseImmiscibleFluid >( subRegion, viewKeyStruct::fluidNamesString(), "two phase immiscible fluid" );
  setConstitutiveName< RelativePermeabilityBase >( subRegion, viewKeyStruct::relPermNamesString(), "relative permeability" );
}

void ImmiscibleMultiphaseFlowMFD::registerDataOnMesh( Group & meshBodies )
{
  FlowSolverBase::registerDataOnMesh( meshBodies );

  
  forDiscretizationOnMeshTargets( meshBodies, [&]( string const &, MeshLevel & mesh, string_array const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();
    elemManager.forElementSubRegions( regionNames, [&]( localIndex const, ElementSubRegionBase & subRegion )
    {
      // primary & secondary fields (phase volume fraction, mass, mobility)
      subRegion.registerField< immiscibleMultiphaseFlow::phaseVolumeFraction >( getName() )
              .reference().resizeDimension< 1 >( m_numPhases );
      subRegion.registerField< immiscibleMultiphaseFlow::phaseVolumeFraction_n >( getName() )
              .reference().resizeDimension< 1 >( m_numPhases );
      subRegion.registerField< immiscibleMultiphaseFlow::bcPhaseVolumeFraction >( getName() )
              .reference().resizeDimension< 1 >( m_numPhases );
      subRegion.registerField< immiscibleMultiphaseFlow::phaseMass >( getName() )
              .reference().resizeDimension< 1 >( m_numPhases );
      subRegion.registerField< immiscibleMultiphaseFlow::phaseMass_n >( getName() )
              .reference().resizeDimension< 1 >( m_numPhases );
      subRegion.registerField< immiscibleMultiphaseFlow::phaseMobility >( getName() )
              .reference().resizeDimension< 1 >( m_numPhases );
      subRegion.registerField< immiscibleMultiphaseFlow::dPhaseMobility >( getName() )
              .reference().resizeDimension< 1, 2 >( m_numPhases, 2 );

      if( m_hasCapPressure )
      {
        subRegion.registerWrapper< string >( viewKeyStruct::capPressureNamesString() )
                 .setPlotLevel( PlotLevel::NOPLOT ).setRestartFlags( RestartFlags::NO_WRITE )
                 .setSizedFromParent( 0 );
        string & capName = subRegion.getReference< string >( viewKeyStruct::capPressureNamesString() );
        capName = getConstitutiveName< CapillaryPressureBase >( subRegion );
      }
    } );
    
    // face fields (previous time level face pressure)
    FaceManager & faceManager = mesh.getFaceManager();
    if( !faceManager.hasWrapper( flow::facePressure_n::key() ) )
    {
      faceManager.registerField< flow::facePressure_n >( getName() );
    }
    // register face bcPressure field for enforcing Dirichlet BCs (mirrors SinglePhaseHybridFVM pattern)
    if( !faceManager.hasWrapper( flow::bcPressure::key() ) )
    {
      faceManager.registerField< flow::bcPressure >( getName() );
    }
  } );
}

void ImmiscibleMultiphaseFlowMFD::initializePreSubGroups()
{
  FlowSolverBase::initializePreSubGroups();
  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  // initialize temperature field (isothermal assumption for now)
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &, MeshLevel & mesh, string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const, ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 > const temp = subRegion.getField< flow::temperature >();
      temp.setValues< parallelHostPolicy >( m_inputTemperature );
    } );
  } );
}

void ImmiscibleMultiphaseFlowMFD::initializePostInitialConditionsPreSubGroups()
{
  FlowSolverBase::initializePostInitialConditionsPreSubGroups();
  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &, MeshLevel & mesh, string_array const & regionNames )
  {
    FieldIdentifiers f; f.addElementFields( { flow::pressure::key(), immiscibleMultiphaseFlow::phaseVolumeFraction::key() }, regionNames );
    CommunicationTools::getInstance().synchronizeFields( f, mesh, domain.getNeighbors(), false );
  } );
  // After initial fields are synchronized, compute initial phase mass and set previous state (_n)
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &, MeshLevel & mesh, string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( regionNames, [&]( localIndex const, auto & subRegion )
    {
      // Ensure rock and fluid states are up-to-date before computing masses
      updatePorosityAndPermeability( subRegion );
      updateFluidState( subRegion );
      // Copy current phase state into _n state so accumulation uses consistent initial masses
      copyPhaseStateToPrev( subRegion );
    } );
  } );
  // build region filter
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &, MeshLevel & mesh, string_array const & regionNames )
  {
    ElementRegionManager const & erm = mesh.getElemManager();
    for( auto const & r : regionNames ) m_regionFilter.insert( erm.getRegions().getIndex( r ) );
  } );

  // Pre-compute and initialize face transGgradZ once (only for HybridMimetic discretization)
  {
    NumericalMethodsManager const & nm = domain.getNumericalMethodManager();
    FiniteVolumeManager const & fvManager = nm.getFiniteVolumeManager();
    if( !fvManager.hasGroup< HybridMimeticDiscretization >( m_discretizationName ) )
    {
      return; // Not using hybrid mimetic discretization; nothing to precompute
    }
    HybridMimeticDiscretization const & hm = fvManager.getHybridMimeticDiscretization( m_discretizationName );
    mimeticInnerProduct::MimeticInnerProductBase const & ip = hm.getReference< mimeticInnerProduct::MimeticInnerProductBase >( HybridMimeticDiscretization::viewKeyStruct::innerProductString() );
    real64 const lengthTolerance = domain.getMeshBody( 0 ).getGlobalLengthScale() * m_areaRelTol;

    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &, MeshLevel & mesh, string_array const & regionNames )
    {
      FaceManager & faceManager = mesh.getFaceManager();
      // Temporary accumulators
      array1d< real64 > invSum( faceManager.size() ); invSum.setValues< parallelDevicePolicy<> >( 0.0 );
      array1d< integer > count( faceManager.size() ); count.setValues< parallelDevicePolicy<> >( 0 );

      NodeManager const & nodeManager = mesh.getNodeManager();
      mesh.getElemManager().forElementSubRegionsComplete< CellElementSubRegion >( regionNames, [&]( localIndex const, localIndex const er, localIndex const esr, ElementRegionBase const &, CellElementSubRegion const & subRegion )
      {
        GEOS_UNUSED_VAR( er );
        GEOS_UNUSED_VAR( esr );
        string const & permName = subRegion.getReference< string >( viewKeyStruct::permeabilityNamesString() );
        PermeabilityBase const & permeability = getConstitutiveModel< PermeabilityBase >( subRegion, permName );
        TransGgradZKernelFactory::createAndLaunch< parallelDevicePolicy<> >( ip,
                                                                                nodeManager,
                                                                                faceManager,
                                                                                subRegion,
                                                                                permeability,
                                                                                lengthTolerance,
                                                                                invSum.toView(),
                                                                                count.toView() );
      } );

      // Reduce to effective value per face and write to field
      arrayView1d< real64 > transEff = faceManager.getField< flow::transGgradZ >();
      arrayView1d< integer const > ghost = faceManager.ghostRank();
      forAll< parallelDevicePolicy<> >( faceManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const kf )
      {
        if( ghost[kf] >= 0 ) {
          return;
        }
        real64 const s = invSum[kf];
        integer const c = count[kf];
        transEff[kf] = (c > 0 && s > 0.0) ? static_cast< real64 >( c ) / s : 0.0;
      } );
    } );
  }
}

void ImmiscibleMultiphaseFlowMFD::implicitStepSetup( real64 const & time_n,
                                                     real64 const & dt,
                                                     DomainPartition & domain )
{
  GEOS_UNUSED_VAR( time_n );
  GEOS_UNUSED_VAR( dt );
  GEOS_UNUSED_VAR( domain );

  // save converged state, update porosity/permeability, update fluid state
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &, MeshLevel & mesh, string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( regionNames, [&]( localIndex const, auto & subRegion )
    {
      // save old state
      copyPressureToPrev( subRegion );
      updatePorosityAndPermeability( subRegion );
      updateFluidState( subRegion );
      // Copy current state into _n state so accumulation uses consistent initial masses
      copyPhaseStateToPrev( subRegion );
    } );
  } );
  // face previous
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &, MeshLevel & mesh, string_array const & )
  {
    copyFacePressureToPrev( mesh );
  } );
}

void ImmiscibleMultiphaseFlowMFD::implicitStepComplete( real64 const & time,
                                                        real64 const & dt,
                                                        DomainPartition & domain )
{
  GEOS_UNUSED_VAR( time );
  GEOS_UNUSED_VAR( dt );
  // save converged solid + relperm + capillary states
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &, MeshLevel & mesh, string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const, ElementSubRegionBase & subRegion )
    {
      updateFluidState(subRegion);
    } );
  } );
}

void ImmiscibleMultiphaseFlowMFD::resetStateToBeginningOfStep( DomainPartition & domain )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &, MeshLevel & mesh, string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( regionNames, [&]( localIndex const, auto & subRegion )
    {
      restoreStateFromPrev( subRegion );
      updatePorosityAndPermeability( subRegion );
      updateFluidState( subRegion );
    } );
  } );
  // faces
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &, MeshLevel & mesh, string_array const & )
  {
    restoreFacePressureFromPrev( mesh );
  } );
}

void ImmiscibleMultiphaseFlowMFD::setupDofs( DomainPartition const & domain, DofManager & dofManager ) const
{
  GEOS_UNUSED_VAR( domain );
  // cell unknowns: pressure + (numPhases-1) independent saturations
  dofManager.addField( viewKeyStruct::elemDofFieldString(), FieldLocation::Elem, m_numDofPerCell, getMeshTargets() );
  
  // lagrange multiplier pressure
  dofManager.addField( flow::facePressure::key(), FieldLocation::Face, 1, getMeshTargets() );
  
  // couple element unknowns through flux approximation (face connectivity)
  dofManager.addCoupling( viewKeyStruct::elemDofFieldString(), viewKeyStruct::elemDofFieldString(), DofManager::Connector::Face );

  // couple face unknowns mimetic inner product supported on elements (Element connectivity)
  dofManager.addCoupling( flow::facePressure::key(), flow::facePressure::key(), DofManager::Connector::Elem );
  
  // couple face and cell unknowns through weak gradient  (Element connectivity)
  dofManager.addCoupling( flow::facePressure::key(), viewKeyStruct::elemDofFieldString(), DofManager::Connector::Elem );
}

void ImmiscibleMultiphaseFlowMFD::assembleAccumulationTerm( DomainPartition & domain,
                                                            DofManager const & dofManager,
                                                            CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                            arrayView1d< real64 > const & localRhs ) const
{
  // assemble MFD-specific accumulation: total mass and independent-phase mass
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &, MeshLevel const & mesh, string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const, ElementSubRegionBase const & subRegion )
    {
      string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
      string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );
      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      TwoPhaseImmiscibleFluid const & fluid = getConstitutiveModel< TwoPhaseImmiscibleFluid >( subRegion, fluidName );
      CoupledSolidBase const & solid = getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );
      // this assumes a two-phase system where phase 0 is independent and phase 1 is dependent
      integer const indep = indepPhaseIndex();
      AccumulationMFDKernelFactory::createAndLaunch< parallelDevicePolicy<> >( dofManager.rankOffset(),
                                                                               indep,
                                                                               dofKey,
                                                                               subRegion,
                                                                               fluid,
                                                                               solid,
                                                                               localMatrix,
                                                                               localRhs );
    } );
  } );
}

void ImmiscibleMultiphaseFlowMFD::assembleFluxTerms( real64 const dt,
                                                     DomainPartition const & domain,
                                                     DofManager const & dofManager,
                                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                     arrayView1d< real64 > const & localRhs ) const
{
  NumericalMethodsManager const & nm = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = nm.getFiniteVolumeManager();
  // If the discretization name corresponds to a HybridMimeticDiscretization (and not a FluxApproximationBase), skip this TPFA-style flux assembly.
  if( !fvManager.hasGroup< FluxApproximationBase >( m_discretizationName ) )
  {
    return; // hybrid-only case: flux terms will be handled by assembleFluxTermsHybrid
  }
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );
  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &, MeshLevel const & mesh, string_array const & )
  {
    fluxApprox.forAllStencils( mesh, [&]( auto & stencil )
    {
      auto stencilWrapper = stencil.createKernelWrapper();
      FluxComputeKernelFactory::createAndLaunch< parallelDevicePolicy<> >( m_numPhases,
                                                                           dofManager.rankOffset(),
                                                                           dofKey,
                                                                           m_hasCapPressure,
                                                                           m_useTotalMassEquation,
                                                                           m_gravityDensityScheme == GravityDensityScheme::PhasePresence,
                                                                           getName(),
                                                                           mesh.getElemManager(),
                                                                           stencilWrapper,
                                                                           dt,
                                                                           localMatrix.toViewConstSizes(),
                                                                           localRhs.toView() );
    } );
  } );
}

void ImmiscibleMultiphaseFlowMFD::assembleSystem( real64 const time_n,
                                                  real64 const dt,
                                                  DomainPartition & domain,
                                                  DofManager const & dofManager,
                                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                  arrayView1d< real64 > const & localRhs )
{
  GEOS_UNUSED_VAR( time_n );
  // Ensure fluid and transport properties reflect current state before assembling
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &, MeshLevel & mesh, string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const, ElementSubRegionBase & subRegion )
    {
      // Enforce saturation complementarity before any property updates
      updateFluidState(subRegion);
    } );
  } );

  assembleAccumulationTerm( domain, dofManager, localMatrix, localRhs );
  assembleFluxTermsHybrid( dt, domain, dofManager, localMatrix, localRhs );
}

void ImmiscibleMultiphaseFlowMFD::assembleFluxTermsHybrid( real64 const dt,
                                                           DomainPartition const & domain,
                                                           DofManager const & dofManager,
                                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                           arrayView1d< real64 > const & localRhs )
{
  NumericalMethodsManager const & nm = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = nm.getFiniteVolumeManager();
  HybridMimeticDiscretization const & hm = fvManager.getHybridMimeticDiscretization( m_discretizationName );
  mimeticInnerProduct::MimeticInnerProductBase const & ip = hm.getReference< mimeticInnerProduct::MimeticInnerProductBase >( HybridMimeticDiscretization::viewKeyStruct::innerProductString() );
  string const faceDofKey = dofManager.getKey( flow::facePressure::key() );
  string const elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
  real64 const lengthTolerance = domain.getMeshBody( 0 ).getGlobalLengthScale() * m_areaRelTol;
  // For two-phase: independent saturation index (0 or 1)
  integer const indep = indepPhaseIndex();

  // 2) Launch hybrid assembly using computed properties
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &, MeshLevel const & mesh, string_array const & regionNames )
  {
    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager const & faceManager = mesh.getFaceManager();
    mesh.getElemManager().forElementSubRegionsComplete< CellElementSubRegion >( regionNames, [&]( localIndex const, localIndex const er, localIndex const esr, ElementRegionBase const &, CellElementSubRegion const & subRegion )
    {
      GEOS_UNUSED_VAR( er );
      GEOS_UNUSED_VAR( esr );
      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      TwoPhaseImmiscibleFluid const & fluid = getConstitutiveModel< TwoPhaseImmiscibleFluid >( subRegion, fluidName );
      string const & permName = subRegion.getReference< string >( viewKeyStruct::permeabilityNamesString() );
      PermeabilityBase const & permeability = getConstitutiveModel< PermeabilityBase >( subRegion, permName );

      ElementBasedAssemblyKernelFactory::createAndLaunch< parallelDevicePolicy<> >( dofManager.rankOffset(),
                                                                                     er, esr,
                                                                                     lengthTolerance,
                                                                                     elemDofKey,
                                                                                     faceDofKey,
                                                                                     getName(),
                                                                                     nodeManager,
                                                                                     faceManager,
                                                                                     mesh.getElemManager(),
                                                                                     subRegion,
                                                                                     ip,
                                                                                     fluid,
                                                                                     permeability,
                                                                                     m_regionFilter.toViewConst(),
                                                                                     indep,
                                                                                     dt,
                                                                                     /*assembleCellEq=*/true,
                                                                                     localMatrix,
                                                                                     localRhs );
    } );
  } );
}


void ImmiscibleMultiphaseFlowMFD::applyDirichletBC( real64 const time_n,
                                                    real64 const dt,
                                                    DofManager const & dofManager,
                                                    DomainPartition & domain,
                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                    arrayView1d< real64 > const & localRhs ) const
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
//  if( m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
//  {
//    GEOS_ERROR_IF( !validateDirichletBC( domain, time_n + dt ), getName() + ": inconsistent Dirichlet BCs" );
//  }
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &, MeshLevel & mesh, string_array const & )
  {
//    // ---- Element (cell) pressure Dirichlet BCs ----
//    fsManager.apply< ElementSubRegionBase >( time_n + dt, mesh, flow::pressure::key(), [&]( FieldSpecificationBase const & fs, string const & setName, SortedArrayView< localIndex const > const & target, ElementSubRegionBase & sr, string const & )
//    {
//      GEOS_UNUSED_VAR( fs, setName );
//      // populate bcPressure from pressure specification
//      fs.applyFieldValue< FieldSpecificationEqual, parallelDevicePolicy<> >( target, time_n + dt, sr, flow::bcPressure::key() );
//      arrayView1d< real64 const > bcPres = sr.getField< flow::bcPressure >();
//      arrayView1d< real64 const > pres = sr.getField< flow::pressure >();
//      string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
//      arrayView1d< globalIndex const > dof = sr.getReference< array1d< globalIndex > >( dofKey );
//      arrayView1d< integer const > ghost = sr.ghostRank();
//      globalIndex rankOffset = dofManager.rankOffset();
//      auto lm = localMatrix; auto rhs = localRhs; // capture views
//      forAll< parallelDevicePolicy<> >( target.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
//      {
//        localIndex const ei = target[a];
//        if( ghost[ei] >= 0 ) return;
//        globalIndex const g = dof[ei];
//        localIndex const row = g - rankOffset;
//        real64 rhsVal;
//        FieldSpecificationEqual::SpecifyFieldValue( g, rankOffset, lm, rhsVal, bcPres[ei], pres[ei] );
//        rhs[row] = rhsVal;
//      } );
//    } );

    // ---- Face pressure Dirichlet BCs (hybrid face dofs) ----
    string const faceDofKey = dofManager.getKey( flow::facePressure::key() );
    FaceManager & faceManager = mesh.getFaceManager();
    arrayView1d< real64 const > presFace = faceManager.getField< flow::facePressure >();
    arrayView1d< real64 const > presFaceBC = faceManager.getField< flow::bcPressure >();
    arrayView1d< globalIndex const > faceDofNumber = faceManager.getReference< array1d< globalIndex > >( faceDofKey );
    arrayView1d< integer const > faceGhostRank = faceManager.ghostRank();
    globalIndex const rankOffset = dofManager.rankOffset();

    fsManager.apply< FaceManager >( time_n + dt, mesh, flow::bcPressure::key(), [&]( FieldSpecificationBase const & fs, string const & setName, SortedArrayView< localIndex const > const & targetSet, FaceManager & targetGroup, string const & )
    {
      GEOS_UNUSED_VAR( setName );
      // populate face bcPressure from pressure specification
      fs.applyFieldValue< FieldSpecificationEqual, parallelDevicePolicy<> >( targetSet,
                                                                            time_n + dt,
                                                                            targetGroup,
                                                                            flow::bcPressure::key() );
      // enforce on system rows
      forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
      {
        localIndex const kf = targetSet[a];
        if( faceGhostRank[kf] >= 0 ) return;
        globalIndex const dofIndex = faceDofNumber[kf];
        localIndex const localRow = dofIndex - rankOffset;
        real64 rhsValue;
        FieldSpecificationEqual::SpecifyFieldValue( dofIndex,
                                                    rankOffset,
                                                    localMatrix,
                                                    rhsValue,
                                                    presFaceBC[kf],
                                                    presFace[kf] );
        localRhs[localRow] = rhsValue;
      } );
    } );
  } );
}

void ImmiscibleMultiphaseFlowMFD::applySourceFluxBC( real64 const time,
                                                     real64 const dt,
                                                     DofManager const & dofManager,
                                                     DomainPartition & domain,
                                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                     arrayView1d< real64 > const & localRhs ) const
{
  GEOS_UNUSED_VAR( localMatrix );
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  // Collect all SourceFlux BC names -> id
  std::map< string, localIndex > nameToId; localIndex count = 0;
  fsManager.forSubGroups< SourceFluxBoundaryCondition >( [&]( SourceFluxBoundaryCondition const & bc )
  { nameToId[ bc.getName() ] = count++; } );
  if( count == 0 ) return;

  // Pre-compute scaling factors (set size normalizers) for each source BC
  array1d< globalIndex > setSizes( count );
  computeSourceFluxSizeScalingFactor( time, dt, domain, nameToId, setSizes.toView() );

  // independent saturation component index (0 or 1) for two-phase system
  integer const indep = indepPhaseIndex();
  globalIndex const rankOffset = dofManager.rankOffset();

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &, MeshLevel & mesh, string_array const & )
  {
    fsManager.apply< ElementSubRegionBase, SourceFluxBoundaryCondition >( time + dt, mesh, SourceFluxBoundaryCondition::catalogName(),
      [&]( SourceFluxBoundaryCondition const & fs, string const &, SortedArrayView< localIndex const > const & target, ElementSubRegionBase & sr, string const & )
    {
      if( target.size() == 0 || !sr.hasWrapper( dofKey ) ) return;

      // Accessors
      arrayView1d< globalIndex const > const dof = sr.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< integer const > const ghost = sr.ghostRank();

      // Temporary buffers for rhs contribution computation (per-target element)
      array1d< globalIndex > tmpDof( target.size() );
      array1d< real64 > tmpRhs( target.size() );

      // Compute raw rhs contributions (unscaled) for this BC on its target set
      fs.computeRhsContribution< FieldSpecificationAdd, parallelDevicePolicy<> >( target.toViewConst(),
                                                                                  time + dt,
                                                                                  dt,
                                                                                  sr,
                                                                                  dof,
                                                                                  rankOffset,
                                                                                  localMatrix,
                                                                                  tmpDof.toView(),
                                                                                  tmpRhs.toView(),
                                                                                  [] GEOS_HOST_DEVICE ( localIndex const ){ return 0.0; } );

      // Normalization factor for this BC (total number of target dofs across sets)
      real64 const scale = setSizes[ nameToId.at( fs.getName() ) ];
      if( scale <= 0 ) return; // safety

      integer const phaseId = fs.getComponent();
      auto rhs = localRhs; // capture by value for device lambda

      // Apply source contributions to pressure (total mass) equation
      forAll< parallelDevicePolicy<> >( target.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
      {
        localIndex const ei = target[a];
        if( ghost[ei] >= 0 ) return; // skip ghosts
        globalIndex const base = dof[ei] - rankOffset; // pressure row
        real64 const val = tmpRhs[a] / scale;          // scaled source term
        rhs[base] += val;                              // total mass / pressure equation
      } );

      // Apply source contributions to independent saturation equation (if phase matches)
      if( phaseId == indep )
      {
        forAll< parallelDevicePolicy<> >( target.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
        {
          localIndex const ei = target[a];
          if( ghost[ei] >= 0 ) return;
          globalIndex const base = dof[ei] - rankOffset; // pressure row base
          real64 const val = tmpRhs[a] / scale;
          rhs[base + 1] += val; // saturation row (only one independent saturation dof)
        } );
      }
    } );
  } );
}

void ImmiscibleMultiphaseFlowMFD::applyBoundaryConditions( real64 const time_n,
                                                           real64 const dt,
                                                           DomainPartition & domain,
                                                           DofManager const & dofManager,
                                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                           arrayView1d< real64 > const & localRhs )
{
  applyDirichletBC( time_n, dt, dofManager, domain, localMatrix, localRhs );
  applySourceFluxBC( time_n, dt, dofManager, domain, localMatrix, localRhs );
}

real64 ImmiscibleMultiphaseFlowMFD::calculateResidualNorm( real64 const & time_n,
                                                           real64 const & dt,
                                                           DomainPartition const & domain,
                                                           DofManager const & dofManager,
                                                           arrayView1d< real64 const > const & localRhs )
{
  // For now compute element-based norm only (reuse single equation group). Extend to face rows later.
  GEOS_UNUSED_VAR( time_n, dt );
  array1d< real64 > localNorm( 1 ), localNormalizer( 1 );
  localNorm[0]=0; localNormalizer[0]=0;
  physicsSolverBaseKernels::NormType normType = getNonlinearSolverParameters().normType();
  globalIndex rankOffset = dofManager.rankOffset();
  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &, MeshLevel const & mesh, string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const, ElementSubRegionBase const & sr )
    {
      arrayView1d< globalIndex const > dof = sr.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< integer const > ghost = sr.ghostRank();
      RAJA::ReduceSum< parallelDeviceReduce, real64 > sumSq( 0.0 );
      RAJA::ReduceMax< parallelDeviceReduce, real64 > maxVal( 0.0 );
      RAJA::ReduceSum< parallelDeviceReduce, real64 > normer( 0.0 );
      forAll< parallelDevicePolicy<> >( sr.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
      {
        if( ghost[ei] < 0 )
        {
          globalIndex row = dof[ei]-rankOffset; real64 v = localRhs[row];
          if( normType == physicsSolverBaseKernels::NormType::Linf ) maxVal.max( fabs( v ) );
          else { sumSq += v*v; normer += 1.0; }
        }
      } );
      if( normType == physicsSolverBaseKernels::NormType::Linf ) localNorm[0] = std::max( localNorm[0], maxVal.get() );
      else { localNorm[0] += sumSq.get(); localNormalizer[0] += normer.get(); }
    } );
  } );
  real64 globalNorm;
  if( normType == physicsSolverBaseKernels::NormType::Linf ) physicsSolverBaseKernels::LinfResidualNormHelper::computeGlobalNorm( localNorm[0], globalNorm );
  else physicsSolverBaseKernels::L2ResidualNormHelper::computeGlobalNorm( localNorm[0], localNormalizer[0], globalNorm );
  return globalNorm;
}

void ImmiscibleMultiphaseFlowMFD::applySystemSolution( DofManager const & dofManager,
                                                       arrayView1d< real64 const > const & localSolution,
                                                       real64 const scalingFactor,
                                                       real64 const dt,
                                                       DomainPartition & domain )
{
  GEOS_UNUSED_VAR( dt );
  DofManager::CompMask pressureMask( m_numDofPerCell, 0, 1 );
  dofManager.addVectorToField( localSolution, viewKeyStruct::elemDofFieldString(), flow::pressure::key(), scalingFactor, pressureMask );
  // Manually add saturation increment to the configured independent phase component
  integer const indep = indepPhaseIndex();
  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &, MeshLevel & mesh, string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const, ElementSubRegionBase & subRegion )
    {
      if( !subRegion.hasWrapper( dofKey ) ) return;
      arrayView1d< globalIndex const > dof = subRegion.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< integer const > ghost = subRegion.ghostRank();
      arrayView2d< real64, immiscibleFlow::USD_PHASE > sat = subRegion.getField< immiscibleMultiphaseFlow::phaseVolumeFraction >();
      globalIndex const rankOffset = dofManager.rankOffset();
      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
      {
        if( ghost[ei] >= 0 ) return;
        globalIndex const base = dof[ei] - rankOffset;
        // saturation dof is component 1 of elem dofs
        sat[ei][indep] += scalingFactor * localSolution[base + 1];
        sat[ei][m_dependentPhaseIndex] = 1.0 - sat[ei][indep];
      } );
    } );
  } );
  // update face pressures
  dofManager.addVectorToField( localSolution, flow::facePressure::key(), flow::facePressure::key(), scalingFactor );
  // synchronize
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &, MeshLevel & mesh, string_array const & regionNames )
  {
    FieldIdentifiers f; f.addElementFields( { flow::pressure::key(), immiscibleMultiphaseFlow::phaseVolumeFraction::key() }, regionNames );
    f.addFields( FieldLocation::Face, { flow::facePressure::key() } );
    CommunicationTools::getInstance().synchronizeFields( f, mesh, domain.getNeighbors(), true );
  } );
}

void ImmiscibleMultiphaseFlowMFD::updateFluidModel( ObjectManagerBase & group ) const
{
  // Update fluid constitutive state (densities/viscosities) from current pressure
  arrayView1d< real64 const > const pres = group.getField< flow::pressure >();
  TwoPhaseImmiscibleFluid & fluid = getConstitutiveModel< TwoPhaseImmiscibleFluid >( group, group.getReference< string >( viewKeyStruct::fluidNamesString() ) );
  constitutive::constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
  {
    using FluidType = TYPEOFREF( castedFluid );
    typename FluidType::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();
    FluidUpdateKernel::launch< parallelDevicePolicy<> >( group.size(), fluidWrapper, pres );
  } );
}

void ImmiscibleMultiphaseFlowMFD::updateRelPermModel( ObjectManagerBase & group ) const
{
  if( !group.hasWrapper( viewKeyStruct::relPermNamesString() ) ) return;
  string const & name = group.getReference< string >( viewKeyStruct::relPermNamesString() );
  RelativePermeabilityBase & relperm = getConstitutiveModel< RelativePermeabilityBase >( group, name );
  arrayView2d< real64, immiscibleFlow::USD_PHASE > const phaseVol = group.getField< immiscibleMultiphaseFlow::phaseVolumeFraction >();
  constitutive::constitutiveUpdatePassThru( relperm, [&]( auto & casted )
  {
    auto wrapper = casted.createKernelWrapper();
    RelativePermeabilityUpdateKernel::launch< parallelDevicePolicy<> >( group.size(), wrapper, phaseVol );
  } );
}

void ImmiscibleMultiphaseFlowMFD::updatePhaseMobility( ObjectManagerBase & group ) const
{
  // Compute phase mobilities and their derivatives using relperm and fluid models
  if( !group.hasWrapper( viewKeyStruct::relPermNamesString() ) ) return;
  string const & relpermName = group.getReference< string >( viewKeyStruct::relPermNamesString() );
  string const & fluidName = group.getReference< string >( viewKeyStruct::fluidNamesString() );
  RelativePermeabilityBase & relperm = getConstitutiveModel< RelativePermeabilityBase >( group, relpermName );
  TwoPhaseImmiscibleFluid & fluid = getConstitutiveModel< TwoPhaseImmiscibleFluid >( group, fluidName );
  constitutive::constitutiveUpdatePassThru( relperm, [&]( auto & castedRelperm )
  {
    GEOS_UNUSED_VAR( castedRelperm );
    // Use existing kernel to populate phaseMobility and dPhaseMobility on the group
    immiscibleMultiphaseKernels::PhaseMobilityKernelFactory::createAndLaunch< parallelDevicePolicy<> >( m_numPhases,
                                                                                                        group,
                                                                                                        fluid,
                                                                                                        relperm );
  } );
}

void ImmiscibleMultiphaseFlowMFD::updateCapPressureModel( ObjectManagerBase & group ) const
{
  if( !m_hasCapPressure ) return;
  if( !group.hasWrapper( viewKeyStruct::capPressureNamesString() ) ) return;
  string const & name = group.getReference< string >( viewKeyStruct::capPressureNamesString() );
  CapillaryPressureBase & cap = getConstitutiveModel< CapillaryPressureBase >( group, name );
  auto phaseVol = group.getField< immiscibleMultiphaseFlow::phaseVolumeFraction >();
  constitutive::constitutiveUpdatePassThru( cap, [&]( auto & casted )
  {
    auto wrapper = casted.createKernelWrapper();
    CapillaryPressureUpdateKernel::launch< parallelDevicePolicy<> >( group.size(), wrapper, phaseVol );
  } );
}

// --- Multiphase helper updates ---

void ImmiscibleMultiphaseFlowMFD::updatePhaseMass( ElementSubRegionBase & subRegion ) const
{
  GEOS_MARK_FUNCTION;
  string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );
  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  TwoPhaseImmiscibleFluid const & fluid = getConstitutiveModel< TwoPhaseImmiscibleFluid >( subRegion, fluidName );
  CoupledSolidBase const & solid = getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );
  arrayView1d< real64 const > const volume = subRegion.getElementVolume();
  arrayView2d< real64 const > const porosity = solid.getPorosity();
  arrayView2d< real64 const, immiscibleFlow::USD_PHASE > const phaseVolFrac = subRegion.getField< fields::immiscibleMultiphaseFlow::phaseVolumeFraction >();
  arrayView3d< real64 const, constitutive::multifluid::USD_PHASE > const phaseDens = fluid.phaseDensity();
  arrayView2d< real64, immiscibleFlow::USD_PHASE > const phaseMass = subRegion.getField< fields::immiscibleMultiphaseFlow::phaseMass >();
  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    real64 const poreVolume = volume[ei] * porosity[ei][0];
    for( integer ip=0; ip<2; ++ip )
    {
      phaseMass[ei][ip] = poreVolume * phaseVolFrac[ei][ip] * phaseDens[ei][0][ip];
    }
  } );
}

void ImmiscibleMultiphaseFlowMFD::updateVolumeConstraint( ElementSubRegionBase & subRegion ) const
{
  auto phaseVol = subRegion.getField< immiscibleMultiphaseFlow::phaseVolumeFraction >();
  integer const dep = m_dependentPhaseIndex; // 0 or 1
  integer const ind = indepPhaseIndex();
  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    phaseVol[ei][dep] = 1.0 - phaseVol[ei][ind];
  } );
}

void ImmiscibleMultiphaseFlowMFD::updateFluidState( ElementSubRegionBase & subRegion ) const
{
  updateFluidModel( subRegion );
  updateVolumeConstraint( subRegion );
  updatePhaseMass( subRegion );
  
  updateRelPermModel( subRegion );
  updatePhaseMobility( subRegion );
  updateCapPressureModel( subRegion );
}

void ImmiscibleMultiphaseFlowMFD::updateState( DomainPartition & domain )
{
  // Mirror ImmiscibleMultiphaseFlow::updateState: update rock, enforce saturation constraint, then update fluid-related state.
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &, MeshLevel & mesh, string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( regionNames, [&]( localIndex const, auto & subRegion )
    {
      updatePorosityAndPermeability( subRegion );
      updateFluidState( subRegion );
    } );
  } );
}

REGISTER_CATALOG_ENTRY( PhysicsSolverBase, ImmiscibleMultiphaseFlowMFD, string const &, Group * const )

} // namespace geos
