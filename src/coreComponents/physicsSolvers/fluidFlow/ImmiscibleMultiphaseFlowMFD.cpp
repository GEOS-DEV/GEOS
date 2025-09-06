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

namespace
{
char const bcLogMessage[] =
  "ImmiscibleMultiphaseFlowMFD {}: at time {}s, the <{}> boundary condition '{}' is applied to set '{}' in subRegion '{}'. Total target elements (incl. ghosts) = {}";
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
  m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::immiscibleMultiphaseFVM; // placeholder
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

  // detect capillary pressure presence & register fields
  forDiscretizationOnMeshTargets( meshBodies, [&]( string const &, MeshLevel & mesh, string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const, ElementSubRegionBase & subRegion )
    {
      string const capPresName = getConstitutiveName< CapillaryPressureBase >( subRegion );
      if( !capPresName.empty() ) m_hasCapPressure = true;
    } );
  } );

  forDiscretizationOnMeshTargets( meshBodies, [&]( string const &, MeshLevel & mesh, string_array const & regionNames )
  {
    // element fields
    mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const, ElementSubRegionBase & subRegion )
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
              .reference().resizeDimension< 1, 2 >( m_numPhases, m_numPhases );
      // pressure gradient diagnostic
      if( !subRegion.hasWrapper( flow::pressureGradient::key() ) )
      {
        subRegion.registerField< flow::pressureGradient >( getName() ).reference().resizeDimension< 1 >( 3 );
      }
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
  // build region filter
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &, MeshLevel & mesh, string_array const & regionNames )
  {
    ElementRegionManager const & erm = mesh.getElemManager();
    for( auto const & r : regionNames ) m_regionFilter.insert( erm.getRegions().getIndex( r ) );
  } );
}

void ImmiscibleMultiphaseFlowMFD::implicitStepSetup( real64 const & time_n,
                                                     real64 const & dt,
                                                     DomainPartition & domain )
{
  // save converged state, update porosity/permeability, update fluid state
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &, MeshLevel & mesh, string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( regionNames, [&]( localIndex const, auto & subRegion )
    {
      // save old state
      arrayView1d< real64 > const pres_n = subRegion.template getField< flow::pressure_n >();
      arrayView1d< real64 const > const pres = subRegion.template getField< flow::pressure >();
      pres_n.template setValues< parallelDevicePolicy<> >( pres );
      updatePorosityAndPermeability( subRegion );
      updateVolumeConstraint( subRegion );
      updateFluidState( subRegion );
      // copy phase fields to _n
      auto const phaseVol = subRegion.template getField< immiscibleMultiphaseFlow::phaseVolumeFraction >();
      auto phaseVol_n = subRegion.template getField< immiscibleMultiphaseFlow::phaseVolumeFraction_n >();
      phaseVol_n.template setValues< parallelDevicePolicy<> >( phaseVol );
      auto const phaseMass = subRegion.template getField< immiscibleMultiphaseFlow::phaseMass >();
      auto phaseMass_n = subRegion.template getField< immiscibleMultiphaseFlow::phaseMass_n >();
      phaseMass_n.template setValues< parallelDevicePolicy<> >( phaseMass );
    } );
  } );
  // face previous
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &, MeshLevel & mesh, string_array const & )
  {
    FaceManager & fm = mesh.getFaceManager();
    auto faceP = fm.getField< flow::facePressure >();
    auto faceP_n = fm.getField< flow::facePressure_n >();
    faceP_n.setValues< parallelDevicePolicy<> >( faceP );
  } );
}

void ImmiscibleMultiphaseFlowMFD::implicitStepComplete( real64 const & time,
                                                        real64 const & dt,
                                                        DomainPartition & domain )
{
  // save converged solid + relperm + capillary states
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &, MeshLevel & mesh, string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const, ElementSubRegionBase & subRegion )
    {
      // update converged saturation for relperm hysteresis if model supports it
      auto phaseVol = subRegion.getField< immiscibleMultiphaseFlow::phaseVolumeFraction >();
      GEOS_UNUSED_VAR( phaseVol );
      GEOS_UNUSED_VAR( time ); GEOS_UNUSED_VAR( dt );
    } );
  } );
}

void ImmiscibleMultiphaseFlowMFD::resetStateToBeginningOfStep( DomainPartition & domain )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &, MeshLevel & mesh, string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion, SurfaceElementSubRegion >( regionNames, [&]( localIndex const, auto & subRegion )
    {
      auto pres = subRegion.template getField< flow::pressure >();
      auto pres_n = subRegion.template getField< flow::pressure_n >();
      pres.template setValues< parallelDevicePolicy<> >( pres_n );
      auto phaseVol = subRegion.template getField< immiscibleMultiphaseFlow::phaseVolumeFraction >();
      auto phaseVol_n = subRegion.template getField< immiscibleMultiphaseFlow::phaseVolumeFraction_n >();
      phaseVol.template setValues< parallelDevicePolicy<> >( phaseVol_n );
      auto phaseMass = subRegion.template getField< immiscibleMultiphaseFlow::phaseMass >();
      auto phaseMass_n = subRegion.template getField< immiscibleMultiphaseFlow::phaseMass_n >();
      phaseMass.template setValues< parallelDevicePolicy<> >( phaseMass_n );
      updatePorosityAndPermeability( subRegion );
      updateFluidState( subRegion );
    } );
  } );
  // faces
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &, MeshLevel & mesh, string_array const & )
  {
    FaceManager & fm = mesh.getFaceManager();
    auto faceP = fm.getField< flow::facePressure >();
    auto faceP_n = fm.getField< flow::facePressure_n >();
    faceP.setValues< parallelDevicePolicy<> >( faceP_n );
  } );
}

void ImmiscibleMultiphaseFlowMFD::setupDofs( DomainPartition const & domain, DofManager & dofManager ) const
{
  dofManager.addField( viewKeyStruct::elemDofFieldString(), FieldLocation::Elem, m_numDofPerCell, getMeshTargets() );
  // couple element unknowns through flux approximation (face connectivity)
  dofManager.addCoupling( viewKeyStruct::elemDofFieldString(), viewKeyStruct::elemDofFieldString(), DofManager::Connector::Face );
  // face pressure
  dofManager.addField( flow::facePressure::key(), FieldLocation::Face, 1, getMeshTargets() );
  dofManager.addCoupling( flow::facePressure::key(), flow::facePressure::key(), DofManager::Connector::Elem );
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
      integer const indep = 1 - m_dependentPhaseIndex;
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
  // Ensure fluid properties reflect current pressure before assembling
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &, MeshLevel & mesh, string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions( regionNames, [&]( localIndex const, ElementSubRegionBase & subRegion )
    {
      updateFluidModel( subRegion );
    } );
  } );
  assembleAccumulationTerm( domain, dofManager, localMatrix, localRhs );
  assembleFluxTerms( dt, domain, dofManager, localMatrix, localRhs ); // TPFA transport + pressure (will later separate)
  assembleFluxTermsHybrid( dt, domain, dofManager, localMatrix, localRhs ); // face constraints only (cell eq disabled)
}

void ImmiscibleMultiphaseFlowMFD::assembleFluxTermsHybrid( real64 const dt,
                                                           DomainPartition const & domain,
                                                           DofManager const & dofManager,
                                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                           arrayView1d< real64 > const & localRhs )
{
  NumericalMethodsManager const & nm = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = nm.getFiniteVolumeManager();
  // Skip if the discretization name corresponds to a FluxApproximation (i.e. not hybrid)
  if( !fvManager.hasGroup< HybridMimeticDiscretization >( m_discretizationName ) )
  {
    return;
  }
  HybridMimeticDiscretization const & hm = fvManager.getHybridMimeticDiscretization( m_discretizationName );
  mimeticInnerProduct::MimeticInnerProductBase const & ip = hm.getReference< mimeticInnerProduct::MimeticInnerProductBase >( HybridMimeticDiscretization::viewKeyStruct::innerProductString() );
  string const faceDofKey = dofManager.getKey( flow::facePressure::key() );
  string const elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
  real64 const lengthTolerance = domain.getMeshBody( 0 ).getGlobalLengthScale() * m_areaRelTol;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &, MeshLevel const & mesh, string_array const & regionNames )
  {
    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager const & faceManager = mesh.getFaceManager();
    mesh.getElemManager().forElementSubRegionsComplete< CellElementSubRegion >( regionNames, [&]( localIndex const, localIndex const er, localIndex const esr, ElementRegionBase const &, CellElementSubRegion const & subRegion )
    {
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
                                                                                     dt,
                                                                                     /*assembleCellEq=*/false,
                                                                                     localMatrix,
                                                                                     localRhs );
    } );
  } );
}

bool ImmiscibleMultiphaseFlowMFD::validateDirichletBC( DomainPartition & domain, real64 const time ) const
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  bool ok = true;
  constexpr integer MAX_NP = 2;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &, MeshLevel & mesh, string_array const & )
  {
    // map: region -> subRegion -> set -> component mask
    map< string, map< string, map< string, ComponentMask< MAX_NP > > > > bcStatus;
    // pressure specs
    fsManager.apply< ElementSubRegionBase >( time, mesh, flow::pressure::key(), [&]( FieldSpecificationBase const &, string const & setName, SortedArrayView< localIndex const > const &, ElementSubRegionBase & sr, string const & )
    {
      string const regionName = sr.getParent().getParent().getName();
      string const subRegionName = sr.getName();
      auto & setMap = bcStatus[regionName][subRegionName];
      if( setMap.count( setName ) > 0 )
      {
        ok = false; // duplicate pressure BC on same set
      }
      setMap[setName].setNumComp( m_numPhases );
    } );
    // phase volume fraction specs
    fsManager.apply< ElementSubRegionBase >( time, mesh, immiscibleMultiphaseFlow::phaseVolumeFraction::key(), [&]( FieldSpecificationBase const & fs, string const & setName, SortedArrayView< localIndex const > const &, ElementSubRegionBase & sr, string const & )
    {
      string const regionName = sr.getParent().getParent().getName();
      string const subRegionName = sr.getName();
      integer const comp = fs.getComponent();
      auto & setMap = bcStatus[regionName][subRegionName];
      if( setMap.count( setName ) == 0 )
      {
        ok = false; // saturation BC without pressure BC on same set
      }
      if( comp < 0 || comp >= m_numPhases )
      {
        ok = false; return;
      }
      ComponentMask< MAX_NP > & mask = setMap[setName];
      if( mask[comp] )
      {
        ok = false; // duplicate component specification
      }
      mask.set( comp );
    } );
  } );
  return ok;
}

void ImmiscibleMultiphaseFlowMFD::applyDirichletBC( real64 const time_n,
                                                    real64 const dt,
                                                    DofManager const & dofManager,
                                                    DomainPartition & domain,
                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                    arrayView1d< real64 > const & localRhs ) const
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  if( m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
  {
    GEOS_ERROR_IF( !validateDirichletBC( domain, time_n + dt ), getName() + ": inconsistent Dirichlet BCs" );
  }
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &, MeshLevel & mesh, string_array const & )
  {
    fsManager.apply< ElementSubRegionBase >( time_n + dt, mesh, flow::pressure::key(), [&]( FieldSpecificationBase const & fs, string const & setName, SortedArrayView< localIndex const > const & target, ElementSubRegionBase & sr, string const & )
    {
      GEOS_UNUSED_VAR( fs, setName );
      // apply to bcPressure then enforce
      fs.applyFieldValue< FieldSpecificationEqual, parallelDevicePolicy<> >( target, time_n + dt, sr, flow::bcPressure::key() );
      arrayView1d< real64 const > bcPres = sr.getField< flow::bcPressure >();
      arrayView1d< real64 const > pres = sr.getField< flow::pressure >();
      string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
      arrayView1d< globalIndex const > dof = sr.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< integer const > ghost = sr.ghostRank();
      globalIndex rankOffset = dofManager.rankOffset();
      auto lm = localMatrix; auto rhs = localRhs; // capture copies
      forAll< parallelDevicePolicy<> >( target.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
      {
        localIndex const ei = target[a]; if( ghost[ei] >= 0 ) return; globalIndex const g = dof[ei]; localIndex const row = g - rankOffset; real64 rhsVal; FieldSpecificationEqual::SpecifyFieldValue( g, rankOffset, lm, rhsVal, bcPres[ei], pres[ei] ); rhs[row] = rhsVal; }
      );
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
  std::map< string, localIndex > nameToId; localIndex count=0; fsManager.forSubGroups< SourceFluxBoundaryCondition >( [&]( SourceFluxBoundaryCondition const & bc ){ nameToId[bc.getName()] = count++; } );
  if( count==0 ) return;
  array1d< globalIndex > setSizes( count );
  computeSourceFluxSizeScalingFactor( time, dt, domain, nameToId, setSizes.toView() );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &, MeshLevel & mesh, string_array const & )
  {
    fsManager.apply< ElementSubRegionBase, SourceFluxBoundaryCondition >( time+dt, mesh, SourceFluxBoundaryCondition::catalogName(), [&]( SourceFluxBoundaryCondition const & fs, string const &, SortedArrayView< localIndex const > const & target, ElementSubRegionBase & sr, string const & )
    {
      if( target.size()==0 || !sr.hasWrapper( dofKey ) ) return;
      auto dof = sr.getReference< array1d< globalIndex > >( dofKey );
      auto ghost = sr.ghostRank();
      array1d< globalIndex > tmpDof( target.size() ); array1d< real64 > tmpRhs( target.size() );
      fs.computeRhsContribution< FieldSpecificationAdd, parallelDevicePolicy<> >( target.toViewConst(), time+dt, dt, sr, dof, dofManager.rankOffset(), localMatrix, tmpDof.toView(), tmpRhs.toView(), [] GEOS_HOST_DEVICE ( localIndex const ){ return 0.0; } );
      real64 scale = setSizes[nameToId.at( fs.getName() )]; integer phaseId = fs.getComponent();
      integer const indep = 1 - m_dependentPhaseIndex;
      auto rhs = localRhs; globalIndex rankOffset = dofManager.rankOffset();
      forAll< parallelDevicePolicy<> >( target.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
      { localIndex const ei=target[a]; if( ghost[ei] >=0 ) return; globalIndex const base = dof[ei]-rankOffset; real64 val = tmpRhs[a]/scale; // total mass equation
        rhs[base] += val; // independent phase equation only if phase matches independent
        if( phaseId == indep ) rhs[base+1] += val; } );
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
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  string const faceDofKey = dofManager.getKey( flow::facePressure::key() );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &, MeshLevel & mesh, string_array const & )
  {
    FaceManager & fm = mesh.getFaceManager();
    auto faceP = fm.getField< flow::facePressure >();
    auto faceDof = fm.getReference< array1d< globalIndex > >( faceDofKey );
    auto ghost = fm.ghostRank();
    fsManager.apply< FaceManager >( time_n+dt, mesh, flow::pressure::key(), [&]( FieldSpecificationBase const & fs, string const &, SortedArrayView< localIndex const > const & target, FaceManager & group, string const & )
    {
      fs.applyFieldValue< FieldSpecificationEqual, parallelDevicePolicy<> >( target, time_n+dt, group, flow::facePressure::key() );
      globalIndex rankOffset = dofManager.rankOffset(); auto lm = localMatrix; auto rhs = localRhs;
      forAll< parallelDevicePolicy<> >( target.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
      { localIndex const fi = target[a]; if( ghost[fi] >=0 ) return; globalIndex const g = faceDof[fi]; localIndex const row = g - rankOffset; real64 rhsVal; FieldSpecificationEqual::SpecifyFieldValue( g, rankOffset, lm, rhsVal, faceP[fi], faceP[fi] ); rhs[row] = rhsVal; });
    } );
  } );
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
  integer const indep = 1 - m_dependentPhaseIndex;
  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString );
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
      } );
      // enforce dependent phase relation s_dep = 1 - s_ind
      updateVolumeConstraint( subRegion );
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

// --- Multiphase helper updates ---

void ImmiscibleMultiphaseFlowMFD::updatePhaseMass( ElementSubRegionBase & subRegion ) const
{
  string const & solidName = subRegion.getReference< string >( viewKeyStruct::solidNamesString() );
  string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
  TwoPhaseImmiscibleFluid const & fluid = getConstitutiveModel< TwoPhaseImmiscibleFluid >( subRegion, fluidName );
  CoupledSolidBase const & solid = getConstitutiveModel< CoupledSolidBase >( subRegion, solidName );
  auto vol = subRegion.getElementVolume();
  auto por = solid.getPorosity();
  auto phaseVol = subRegion.getField< immiscibleMultiphaseFlow::phaseVolumeFraction >();
  auto dens = fluid.phaseDensity();
  auto phaseMass = subRegion.getField< immiscibleMultiphaseFlow::phaseMass >();
  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
  {
    real64 pv = vol[ei]*por[ei][0];
    for( integer ip=0; ip<2; ++ip ) phaseMass[ei][ip] = pv * phaseVol[ei][ip] * dens[ei][0][ip];
  } );
}

void ImmiscibleMultiphaseFlowMFD::updateVolumeConstraint( ElementSubRegionBase & subRegion ) const
{
  auto phaseVol = subRegion.getField< immiscibleMultiphaseFlow::phaseVolumeFraction >();
  integer const dep = ( m_dependentPhaseIndex == 0 ? 0 : 1 );
  integer const ind = 1 - dep;
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

REGISTER_CATALOG_ENTRY( PhysicsSolverBase, ImmiscibleMultiphaseFlowMFD, string const &, Group * const )

} // namespace geos
