/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * (c) GEOS/GEOSX Contributors
 * ------------------------------------------------------------------------------------------------------------
 */

#include "ImmiscibleMultiphaseFlowMFD.hpp"

#include "ImmiscibleMultiphaseFlow.hpp"
#include "physicsSolvers/fluidFlow/ImmiscibleMultiphaseFlowFields.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/HybridMimeticDiscretization.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/fluidFlow/kernels/immiscibleMultiphase/ImmiscibleMultiphaseMFDKernels.hpp"
#include "physicsSolvers/fluidFlow/kernels/singlePhase/ResidualNormKernel.hpp" // reuse residual norm helpers
#include "physicsSolvers/PhysicsSolverBaseKernels.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fieldSpecification/FieldSpecificationBase.hpp"
#include "fieldSpecification/DirichletBoundaryCondition.hpp"
#include "constitutive/fluid/twophaseimmisciblefluid/TwoPhaseImmiscibleFluid.hpp"
#include "constitutive/permeability/PermeabilityBase.hpp"

namespace geos
{

using namespace dataRepository;
using namespace fields;
using namespace immiscibleMultiphaseMFDKernels;

ImmiscibleMultiphaseFlowMFD::ImmiscibleMultiphaseFlowMFD( const string & name,
                                                          Group * const parent )
  : ImmiscibleMultiphaseFlow( name, parent ),
    m_areaRelTol( 1e-8 )
{
  // switch linear solver strategy to custom (placeholder)
  m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::immiscibleMultiphaseFVM; // reuse existing for now
}

void ImmiscibleMultiphaseFlowMFD::registerDataOnMesh( Group & meshBodies )
{
  ImmiscibleMultiphaseFlow::registerDataOnMesh( meshBodies );

  // Add face pressure unknowns + gradient field (mirrors SinglePhaseHybridFVM)
  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    string_array const & regionNames )
  {
    // cell centered gradient (optional diagnostic)
    mesh.getElemManager().forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                                        [&]( localIndex const,
                                                                             ElementSubRegionBase & subRegion )
    {
      if( !subRegion.hasField( flow::pressureGradient::key() ) )
      {
        subRegion.registerField< flow::pressureGradient >( getName() ).reference().resizeDimension< 1 >( 3 );
      }
    } );

    FaceManager & faceManager = mesh.getFaceManager();
    if( !faceManager.hasField( flow::facePressure_n::key() ) )
    {
      faceManager.registerField< flow::facePressure_n >( getName() );
    }
  } );
}

void ImmiscibleMultiphaseFlowMFD::initializePreSubGroups()
{
  ImmiscibleMultiphaseFlow::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  NumericalMethodsManager const & nm = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = nm.getFiniteVolumeManager();

  GEOS_THROW_IF( !fvManager.hasGroup< HybridMimeticDiscretization >( m_discretizationName ),
                 getCatalogName() << " " << getDataContext() << ": HybridMimeticDiscretization required for ImmiscibleMultiphaseFlowMFD",
                 InputError );
}

void ImmiscibleMultiphaseFlowMFD::initializePostInitialConditionsPreSubGroups()
{
  ImmiscibleMultiphaseFlow::initializePostInitialConditionsPreSubGroups();

  // build region filter
  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {
    ElementRegionManager const & erm = mesh.getElemManager();
    for( string const & regionName : regionNames )
    {
      m_regionFilter.insert( erm.getRegions().getIndex( regionName ) );
    }
  } );
}

void ImmiscibleMultiphaseFlowMFD::implicitStepSetup( real64 const & time_n,
                                                     real64 const & dt,
                                                     DomainPartition & domain )
{
  ImmiscibleMultiphaseFlow::implicitStepSetup( time_n, dt, domain );
  // initialize face pressures previous state
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {
    FaceManager & faceManager = mesh.getFaceManager();
    arrayView1d< real64 const > const faceP = faceManager.getField< flow::facePressure >();
    arrayView1d< real64 > const faceP_n = faceManager.getField< flow::facePressure_n >();
    faceP_n.setValues< parallelDevicePolicy<> >( faceP );
  } );
}

void ImmiscibleMultiphaseFlowMFD::implicitStepComplete( real64 const & time,
                                                        real64 const & dt,
                                                        DomainPartition & domain )
{
  ImmiscibleMultiphaseFlow::implicitStepComplete( time, dt, domain );
}

void ImmiscibleMultiphaseFlowMFD::resetStateToBeginningOfStep( DomainPartition & domain )
{
  ImmiscibleMultiphaseFlow::resetStateToBeginningOfStep( domain );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {
    FaceManager & faceManager = mesh.getFaceManager();
    arrayView1d< real64 > const faceP = faceManager.getField< flow::facePressure >();
    arrayView1d< real64 const > const faceP_n = faceManager.getField< flow::facePressure_n >();
    faceP.setValues< parallelDevicePolicy<> >( faceP_n );
  } );
}

void ImmiscibleMultiphaseFlowMFD::setupDofs( DomainPartition const & domain,
                                             DofManager & dofManager ) const
{
  // First call base to add element cell-centered dofs
  ImmiscibleMultiphaseFlow::setupDofs( domain, dofManager );

  // Add face pressure dofs (1 per face) and couplings like single-phase
  dofManager.addField( flow::facePressure::key(), FieldLocation::Face, 1, getMeshTargets() );
  dofManager.addCoupling( flow::facePressure::key(), flow::facePressure::key(), DofManager::Connector::Elem );
  dofManager.addCoupling( flow::facePressure::key(), ImmiscibleMultiphaseFlow::viewKeyStruct::elemDofFieldString(), DofManager::Connector::Elem );
}

void ImmiscibleMultiphaseFlowMFD::assembleSystem( real64 const time_n,
                                                  real64 const dt,
                                                  DomainPartition & domain,
                                                  DofManager const & dofManager,
                                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                  arrayView1d< real64 > const & localRhs )
{
  // Accumulation & (temporary) standard multiphase flux terms
  assembleAccumulationTerm( domain, dofManager, localMatrix, localRhs );

  // Multiphase transport (reuse base TPFA assembly for now)
  ImmiscibleMultiphaseFlow::assembleFluxTerms( dt, domain, dofManager, localMatrix, localRhs );

  // Hybrid pressure contribution (adds pressure & face equations). This currently adds extra pressure residual; TODO: split pressure flux from transport kernel above.
  assembleFluxTermsHybrid( dt, domain, dofManager, localMatrix, localRhs );
}

void ImmiscibleMultiphaseFlowMFD::assembleFluxTermsHybrid( real64 const dt,
                                                           DomainPartition const & domain,
                                                           DofManager const & dofManager,
                                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                           arrayView1d< real64 > const & localRhs )
{
  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  HybridMimeticDiscretization const & hmDiscretization = fvManager.getHybridMimeticDiscretization( m_discretizationName );
  MimeticInnerProductBase const & ip = hmDiscretization.getReference< MimeticInnerProductBase >( HybridMimeticDiscretization::viewKeyStruct::innerProductString() );

  string const faceDofKey = dofManager.getKey( flow::facePressure::key() );
  string const elemDofKey = dofManager.getKey( ImmiscibleMultiphaseFlow::viewKeyStruct::elemDofFieldString() );

  real64 const lengthTolerance = domain.getMeshBody( 0 ).getGlobalLengthScale() * m_areaRelTol;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               string_array const & regionNames )
  {
    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    mesh.getElemManager().forElementSubRegionsComplete< CellElementSubRegion >( regionNames,
                                                                                [&]( localIndex const,
                                                                                     localIndex const er,
                                                                                     localIndex const esr,
                                                                                     ElementRegionBase const &,
                                                                                     CellElementSubRegion const & subRegion )
    {
      // fluid & permeability models
      string const & fluidName = subRegion.getReference< string >( ImmiscibleMultiphaseFlow::viewKeyStruct::fluidNamesString() );
      TwoPhaseImmiscibleFluid const & fluid = getConstitutiveModel< TwoPhaseImmiscibleFluid >( subRegion, fluidName );
      string const & permName = subRegion.getReference< string >( ImmiscibleMultiphaseFlow::viewKeyStruct::permeabilityNamesString() );
      PermeabilityBase const & permeability = getConstitutiveModel< PermeabilityBase >( subRegion, permName );

      immiscibleMultiphaseMFDKernels::ElementBasedAssemblyKernelFactory::createAndLaunch< parallelDevicePolicy<> >( dofManager.rankOffset(),
                                                                                                                     er,
                                                                                                                     esr,
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

void ImmiscibleMultiphaseFlowMFD::applyBoundaryConditions( real64 const time_n,
                                                           real64 const dt,
                                                           DomainPartition & domain,
                                                           DofManager const & dofManager,
                                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                           arrayView1d< real64 > const & localRhs )
{
  // Apply base element (pressure + saturation) BCs
  ImmiscibleMultiphaseFlow::applyBoundaryConditions( time_n, dt, domain, dofManager, localMatrix, localRhs );

  // Apply face pressure Dirichlet (reuse pattern from SinglePhaseHybridFVM::applyFaceDirichletBC)
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  string const faceDofKey = dofManager.getKey( flow::facePressure::key() );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {
    FaceManager & faceManager = mesh.getFaceManager();
    arrayView1d< real64 const > const facePres = faceManager.getField< flow::facePressure >();
    arrayView1d< globalIndex const > const faceDofNumber = faceManager.getReference< array1d< globalIndex > >( faceDofKey );
    arrayView1d< integer const > const faceGhostRank = faceManager.ghostRank();
    globalIndex const rankOffset = dofManager.rankOffset();

    fsManager.apply< FaceManager >( time_n + dt,
                                    mesh,
                                    flow::pressure::key(),
                                    [&] ( FieldSpecificationBase const & fs,
                                          string const & setName,
                                          SortedArrayView< localIndex const > const & targetSet,
                                          FaceManager & targetGroup,
                                          string const & )
    {
      // populate facePressure values from specified pressure BCs
      fs.applyFieldValue< FieldSpecificationEqual,
                          parallelDevicePolicy<> >( targetSet,
                                                    time_n + dt,
                                                    targetGroup,
                                                    flow::facePressure::key() );

      // enforce in system (overwrite rows)
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
                                                    facePres[kf],
                                                    facePres[kf] );
        localRhs[localRow] = rhsValue;
      } );
    } );
  } );
}

real64 ImmiscibleMultiphaseFlowMFD::calculateResidualNorm( real64 const & time_n,
                                                           real64 const & dt,
                                                           DomainPartition const & domain,
                                                           DofManager const & dofManager,
                                                           arrayView1d< real64 const > const & localRhs )
{
  // Reuse base (element) residual norm + skip face residual contribution for now (TODO add face contribution similar to single-phase)
  return ImmiscibleMultiphaseFlow::calculateResidualNorm( time_n, dt, domain, dofManager, localRhs );
}

void ImmiscibleMultiphaseFlowMFD::applySystemSolution( DofManager const & dofManager,
                                                       arrayView1d< real64 const > const & localSolution,
                                                       real64 const scalingFactor,
                                                       real64 const dt,
                                                       DomainPartition & domain )
{
  GEOS_UNUSED_VAR( dt );
  // Pressure + saturation update (base)
  ImmiscibleMultiphaseFlow::applySystemSolution( dofManager, localSolution, scalingFactor, dt, domain );
  // Face pressure update
  dofManager.addVectorToField( localSolution,
                               flow::facePressure::key(),
                               flow::facePressure::key(),
                               scalingFactor );
  // recompute pressure gradient diagnostic
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {
    FaceManager const & faceManager = mesh.getFaceManager();
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&] ( localIndex const,
                                                                                           CellElementSubRegion & subRegion )
    {
      immiscibleMultiphaseMFDKernels::AveragePressureGradientKernelFactory::createAndLaunch< parallelDevicePolicy<> >( subRegion, faceManager );
    } );
  } );
  // sync faces
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {
    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addFields( FieldLocation::Face, { flow::facePressure::key() } );
    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), true );
  } );
}

void ImmiscibleMultiphaseFlowMFD::updatePhaseMobility( ObjectManagerBase & dataGroup ) const
{
  // First compute phase mobilities (base implementation)
  ImmiscibleMultiphaseFlow::updatePhaseMobility( dataGroup );
  // Optionally compute aggregated mobility if needed by later enhancements
  computeTotalMobility( dataGroup );
}

void ImmiscibleMultiphaseFlowMFD::computeTotalMobility( ObjectManagerBase & GEOS_UNUSED_PARAM( dataGroup ) ) const
{
  // Placeholder: we could accumulate total mobility into flow::mobility field if registered
}

REGISTER_CATALOG_ENTRY( PhysicsSolverBase, ImmiscibleMultiphaseFlowMFD, string const &, Group * const )

} // namespace geos
