/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file TwoPhasePhaseFlow.cpp
 */

#include "TwoPhaseBase.hpp"

#include "mpiCommunications/CommunicationTools.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "constitutive/fluid/MultiFluidBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"

/**
 * @namespace the geosx namespace that encapsulates the majority of the code
 */
namespace geosx
{

using namespace dataRepository;
using namespace constitutive;


TwoPhaseBase::TwoPhaseBase( const std::string& name,
                            Group * const parent ):
  FlowSolverBase(name, parent)
{
  m_ipw     = -1;
  m_ipnw  = -1;
  m_numDofPerCell =  2;
 
  this->registerWrapper( viewKeyStruct::relPermNameString,  &m_relPermName,  false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of the relative permeability constitutive model to use");

  this->registerWrapper( viewKeyStruct::relPermIndexString, &m_relPermIndex, false );

}

void TwoPhaseBase::RegisterDataOnMesh(Group * const MeshBodies)
{
  FlowSolverBase::RegisterDataOnMesh(MeshBodies);

  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    MeshLevel * meshLevel = Group::group_cast<MeshBody *>(mesh.second)->getMeshLevel(0);

    ElementRegionManager * const elemManager = meshLevel->getElemManager();

    elemManager->forElementSubRegions<CellElementSubRegion>( [&]( CellElementSubRegion * const subRegion )
    {
      // non-wetting phase pressure
      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::pressureString )->setPlotLevel(PlotLevel::LEVEL_0);
      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::deltaPressureString );

      // wetting phase saturation
      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::wettingPhaseSatString )->setPlotLevel(PlotLevel::LEVEL_0);
      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::deltaWettingPhaseSatString );
      subRegion->registerWrapper< array2d<real64> >( viewKeyStruct::phaseSatString );
      
      // auxiliary variables
      subRegion->registerWrapper< array2d<real64> >( viewKeyStruct::phaseMobilityString );
      subRegion->registerWrapper< array2d<real64> >( viewKeyStruct::dPhaseMobility_dPressureString );
      subRegion->registerWrapper< array2d<real64> >( viewKeyStruct::dPhaseMobility_dSaturationString );

      // backup fields
      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::porosityOldString );
      subRegion->registerWrapper< array2d<real64> >( viewKeyStruct::phaseDensityOldString );

    });

  }
}

void TwoPhaseBase::UpdateFluidModel(Group * const dataGroup ) const
{
  GEOSX_MARK_FUNCTION;

  MultiFluidBase * const fluid = GetConstitutiveModel<MultiFluidBase>( dataGroup, m_fluidName );

  // in this function, I currently use hard-coded relations  
  // TODO: find a way to call a constitutive model here
  //       (maybe some variant of SingleFluid with pressure-dependent densities and viscosities
  //       or, better a MultiFluid two-phase Dead-Oil)

  arrayView1d<real64 const> const & pres  =
    dataGroup->getReference< array1d<real64> >( viewKeyStruct::pressureString );
  arrayView1d<real64 const> const & dPres =
    dataGroup->getReference< array1d<real64> >( viewKeyStruct::deltaPressureString );
  
  // phase densities
  arrayView3d<real64> const & phaseDens =
    fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::phaseDensityString );
  arrayView3d<real64> const & dPhaseDens_dPres =
    fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseDensity_dPressureString );

  // phase viscosities
  arrayView3d<real64> const & phaseVisc =
    fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::phaseViscosityString );
  arrayView3d<real64> const & dPhaseVisc_dPres =
    fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseViscosity_dPressureString );

  real64 const p_ref       = 1e6;
  real64 const c_w         = 1e-12;
  real64 const c_nw        = 1e-9;
  real64 const dens_ref_w  = 1000;
  real64 const dens_ref_nw = 800;
  real64 const visc_ref_w  = 0.001;
  real64 const visc_ref_nw = 0.0005;
  forall_in_range( 0, dataGroup->size(), GEOSX_LAMBDA ( localIndex const a )
  {
    // densities
    phaseDens[a][0][m_ipw]         = dens_ref_w * std::exp( c_w * ( pres[a]+dPres[a] - p_ref ) );
    dPhaseDens_dPres[a][0][m_ipw]  = c_w * phaseDens[a][0][m_ipw];
    phaseDens[a][0][m_ipnw]        = dens_ref_nw * std::exp( c_nw * ( pres[a]+dPres[a] - p_ref ) );
    dPhaseDens_dPres[a][0][m_ipnw] = c_nw * phaseDens[a][0][m_ipnw];
    // viscosities
    phaseVisc[a][0][m_ipw]         = visc_ref_w * std::exp( c_w * ( pres[a]+dPres[a] - p_ref ) );
    dPhaseVisc_dPres[a][0][m_ipw]  = c_w * phaseVisc[a][0][m_ipw];
    phaseVisc[a][0][m_ipnw]        = visc_ref_nw * std::exp( c_w * ( pres[a]+dPres[a] - p_ref ) );
    dPhaseVisc_dPres[a][0][m_ipnw] = c_nw * phaseVisc[a][0][m_ipnw];;        
  });
}

void TwoPhaseBase::UpdateSolidModel(Group * const dataGroup ) const
{
  GEOSX_MARK_FUNCTION;

  ConstitutiveBase * const solid = GetConstitutiveModel<ConstitutiveBase>( dataGroup, m_solidName );

  arrayView1d<real64 const> const & pres  =
    dataGroup->getReference< array1d<real64> >( viewKeyStruct::pressureString );
  arrayView1d<real64 const> const & dPres =
    dataGroup->getReference< array1d<real64> >( viewKeyStruct::deltaPressureString );

  forall_in_range( 0, dataGroup->size(), GEOSX_LAMBDA ( localIndex const a )
  {
    solid->StateUpdatePointPressure( pres[a] + dPres[a], a, 0 );
  });
}

void TwoPhaseBase::UpdateRelPermModel( Group * const dataGroup ) const 
{
  GEOSX_MARK_FUNCTION;

  RelativePermeabilityBase * const relPerm = GetConstitutiveModel<RelativePermeabilityBase>( dataGroup, m_relPermName );

  arrayView1d<real64 const> const & sat =
    dataGroup->getReference< array1d<real64> >( viewKeyStruct::wettingPhaseSatString );
  arrayView1d<real64 const> const & dSat =
    dataGroup->getReference< array1d<real64> >( viewKeyStruct::deltaWettingPhaseSatString );

  // in this function, I recompute the phase saturations 
  // TODO: find a better way to do that (maybe work with the array2d phaseSat everywhere)
  
  arrayView2d<real64> const & phaseSat =
    dataGroup->getReference< array2d<real64> >( viewKeyStruct::phaseSatString );

  forall_in_range( 0, dataGroup->size(), GEOSX_LAMBDA ( localIndex const a )
  {
    phaseSat[a][m_ipw]    = sat[a] + dSat[a];
    phaseSat[a][m_ipnw] = 1-phaseSat[a][m_ipw];     
  });
  
  relPerm->BatchUpdate( phaseSat );
}
  
void TwoPhaseBase::UpdatePhaseMobility( Group * const dataGroup ) const
{
  GEOSX_MARK_FUNCTION;

  MultiFluidBase * const fluid = GetConstitutiveModel<MultiFluidBase>( dataGroup, m_fluidName );

  RelativePermeabilityBase * const relPerm = GetConstitutiveModel<RelativePermeabilityBase>( dataGroup, m_relPermName );
  
  // phase relperms
  arrayView3d<real64 const> const & phaseRelPerm =
    relPerm->getReference< array3d<real64> >( RelativePermeabilityBase::viewKeyStruct::phaseRelPermString );
  arrayView4d<real64 const> const & dPhaseRelPerm_dSat =
    relPerm->getReference< array4d<real64> >( RelativePermeabilityBase::viewKeyStruct::dPhaseRelPerm_dPhaseVolFractionString );
  
  // phase viscosities
  arrayView3d<real64 const> const & phaseVisc =
    fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::phaseViscosityString );
  arrayView3d<real64 const> const & dPhaseVisc_dPres =
    fluid->getReference<array3d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseViscosity_dPressureString );

  // phase mobilities
  arrayView2d<real64> const & phaseMob =
    dataGroup->getReference<array2d<real64>>( viewKeyStruct::phaseMobilityString );
  arrayView2d<real64> const & dPhaseMob_dPres =
    dataGroup->getReference<array2d<real64>>( viewKeyStruct::dPhaseMobility_dPressureString );
  arrayView2d<real64> const & dPhaseMob_dSat =
    dataGroup->getReference<array2d<real64>>( viewKeyStruct::dPhaseMobility_dSaturationString );

  // note: unlike in the other flow solvers, the phase mobilities computed here do not include densities
  // we do not include densities here to be able to easily compute the mobility ratios in the flux terms
  
  localIndex constexpr numPhases  = NUM_PHASES;
  forall_in_range( 0, dataGroup->size(), GEOSX_LAMBDA ( localIndex const a )
  {
    for (localIndex ip = 0; ip < numPhases; ++ip)
    {
      phaseMob[a][ip]        = phaseRelPerm[a][0][ip] / phaseVisc[a][0][ip];
      dPhaseMob_dPres[a][ip] = - phaseRelPerm[a][0][ip] * dPhaseVisc_dPres[a][0][ip]
                               / (phaseVisc[a][0][ip] * phaseVisc[a][0][ip]);
      dPhaseMob_dSat[a][ip]  =   dPhaseRelPerm_dSat[a][0][ip][ip] / phaseVisc[a][0][ip];
    }
    
    // this is needed because we need the derivative wrt the wetting-phase saturation
    dPhaseMob_dSat[a][m_ipnw] *= -1;

  }); 
}

void TwoPhaseBase::UpdateState( Group * dataGroup ) const
{
  GEOSX_MARK_FUNCTION;

  UpdateFluidModel( dataGroup );
  UpdateSolidModel( dataGroup );
  UpdateRelPermModel( dataGroup );  
  UpdatePhaseMobility( dataGroup );
}

void TwoPhaseBase::PostProcessInput()
{
  FlowSolverBase::PostProcessInput();
}
  
void TwoPhaseBase::InitializePreSubGroups( Group * const rootGroup )
{
  FlowSolverBase::InitializePreSubGroups( rootGroup );

  DomainPartition * const domain = rootGroup->GetGroup<DomainPartition>( keys::domain );
  ConstitutiveManager const * const cm = domain->getConstitutiveManager();

  
  MultiFluidBase const * fluid = cm->GetConstitutiveRelation<MultiFluidBase>( m_fluidName );

  RelativePermeabilityBase const * relPerm = cm->GetConstitutiveRelation<RelativePermeabilityBase>( m_relPermName );
  GEOSX_ERROR_IF( relPerm == nullptr, "Relative permeability model " + m_relPermName + " not found" );
  m_relPermIndex = relPerm->getIndexInParent();

  // Consistency check between the models
  GEOSX_ERROR_IF( fluid->numFluidPhases() != relPerm->numFluidPhases(),
                 "Number of fluid phases differs between fluid model '" << m_fluidName
                 << "' and relperm model '" << m_relPermName << "'" );

  localIndex constexpr numPhases  = NUM_PHASES;
  GEOSX_ERROR_IF( fluid->numFluidPhases() != numPhases,
                 "Number of fluid phases differs between fluid model '" << m_fluidName
                 << "' and relperm model '" << m_relPermName << "'" );
 

  for (localIndex ip = 0; ip < numPhases; ++ip)
  {
    string const & phase_fl = fluid->phaseName( ip );
    string const & phase_rp = relPerm->phaseName( ip );
    GEOSX_ERROR_IF( phase_fl != phase_rp, "Phase '" << phase_fl << "' in fluid model '"   << m_fluidName
                    << "' does not match phase '"   << phase_rp << "' in relperm model '" << m_relPermName << "'" );
  }


  // determine the indices of the wetting and non-wetting phases
  if ( (fluid->phaseName( 0 ) == "oil" && fluid->phaseName( 1 ) == "gas") ||
       (fluid->phaseName( 1 ) == "oil" && fluid->phaseName( 0 ) == "water") )
  {
    m_ipw    = 0;
    m_ipnw = 1;
  }
  else if ( (fluid->phaseName( 1 ) == "oil" && fluid->phaseName( 0 ) == "gas") ||
            (fluid->phaseName( 0 ) == "oil" && fluid->phaseName( 1 ) == "water") )
  {
    m_ipw    = 1;
    m_ipnw = 0;
  }
  GEOSX_ERROR_IF( m_ipw == -1 || m_ipnw == -1,
                  "TwoPhaseBase: the accepted phase names are water, oil, and gas");

  // fill the array mapping the phase index to the row offset in the residual 
  m_phaseToRow.resize(NUM_PHASES); 
  m_phaseToRow[m_ipw]    = RowOffset::WETTING;
  m_phaseToRow[m_ipnw] = RowOffset::NONWETTING;

  for( auto & mesh : domain->getMeshBodies()->GetSubGroups() )
  {
    MeshLevel * meshLevel = Group::group_cast<MeshBody *>(mesh.second)->getMeshLevel(0);
    ResizeFields( meshLevel );
  }
}

void TwoPhaseBase::ResizeFields( MeshLevel * const meshLevel )
{
  localIndex constexpr numPhases  = NUM_PHASES;

  applyToSubRegions( meshLevel, [&] ( ElementSubRegionBase * const subRegion )
  {
    // this is currently only used to compute the phase relative permeabilities 
    subRegion->getReference< array2d<real64> >(viewKeyStruct::phaseSatString).resizeDimension<1>( numPhases );  
    
    subRegion->getReference< array2d<real64> >(viewKeyStruct::phaseMobilityString).resizeDimension<1>( numPhases );
    subRegion->getReference< array2d<real64> >(viewKeyStruct::dPhaseMobility_dPressureString).resizeDimension<1>( numPhases );
    subRegion->getReference< array2d<real64> >(viewKeyStruct::dPhaseMobility_dSaturationString).resizeDimension<1>( numPhases );

    subRegion->getReference< array2d<real64> >(viewKeyStruct::phaseDensityOldString).resizeDimension<1>( numPhases );

  });
}
  
void TwoPhaseBase::InitializePostInitialConditions_PreSubGroups( Group * const rootGroup )
{
  GEOSX_MARK_FUNCTION;

  FlowSolverBase::InitializePostInitialConditions_PreSubGroups( rootGroup );

  DomainPartition * const domain = rootGroup->GetGroup<DomainPartition>(keys::domain);

  MeshLevel * const mesh = domain->getMeshBody(0)->getMeshLevel(0);

  //TODO this is a hack until the sets are fixed to include ghosts!!
  std::map<string, string_array> fieldNames;
  fieldNames["elems"].push_back( viewKeyStruct::pressureString );
  fieldNames["elems"].push_back( viewKeyStruct::wettingPhaseSatString );

  array1d<NeighborCommunicator> & comms =
    domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors );

  CommunicationTools::SynchronizeFields( fieldNames, mesh, comms );

  ConstitutiveManager * const constitutiveManager = domain->getConstitutiveManager();

  // TODO find a way to set this before constitutive model is duplicated and attached to subregions?
  {
    MultiFluidBase * const fluid = constitutiveManager->GetConstitutiveRelation<MultiFluidBase>( m_fluidName );
    fluid->setMassFlag( static_cast<bool>(true) );
  }
  
  // set mass fraction flag on subregion models
  applyToSubRegions( mesh, [&] ( ElementSubRegionBase * const subRegion )
  {
    MultiFluidBase * const fluid = GetConstitutiveModel<MultiFluidBase>( subRegion, m_fluidName );
    fluid->setMassFlag( static_cast<bool>(true) );
  });

  // initialize primary variables from applied initial conditions
  ResetViews( domain );

  // initialize fluid state
  applyToSubRegions( mesh, [&] ( ElementSubRegionBase * const subRegion )
  {
    UpdateState( subRegion );
  });
}

real64 TwoPhaseBase::SolverStep( real64 const& time_n,
                                 real64 const& dt,
                                 const int cycleNumber,
                                 DomainPartition * domain )
{
  GEOSX_MARK_FUNCTION;

  real64 dt_return;

  ImplicitStepSetup( time_n, dt, domain, m_dofManager, m_matrix, m_rhs, m_solution );

  // currently the only method is implicit time integration
  dt_return = this->NonlinearImplicitStep( time_n, dt, cycleNumber, domain, m_dofManager, m_matrix, m_rhs, m_solution );

  // final step for completion of timestep. typically secondary variable updates and cleanup.
  ImplicitStepComplete( time_n, dt_return, domain );

  return dt_return;
}


void TwoPhaseBase::ImplicitStepSetup( real64 const & GEOSX_UNUSED_ARG( time_n ),
                                      real64 const & GEOSX_UNUSED_ARG( dt ),
                                      DomainPartition * const domain,
                                      DofManager & dofManager,
                                      ParallelMatrix & matrix,
                                      ParallelVector & rhs,
                                      ParallelVector & solution )
{
  // bind the stored views to the current domain
  ResetViews( domain );

  // set deltas to zero and recompute dependent quantities
  ResetStateToBeginningOfStep( domain );

  // backup fields used in time derivative approximation
  BackupFields( domain );

  if( !m_coupledWellsFlag )
  {
    SetupSystem( domain, dofManager, matrix, rhs, solution );
  }

}

void TwoPhaseBase::ImplicitStepComplete( real64 const & GEOSX_UNUSED_ARG( time_n ),
                                         real64 const & GEOSX_UNUSED_ARG( dt ),
                                         DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBody(0)->getMeshLevel(0);

  applyToSubRegions( mesh, [&] ( localIndex const er, localIndex const esr,
                                 ElementRegionBase * const GEOSX_UNUSED_ARG( region ),
                                 ElementSubRegionBase * const subRegion )
  {
    arrayView1d<real64 const> const & dPres = m_deltaPressure[er][esr];
    arrayView1d<real64 const> const & dSat  = m_deltaWettingPhaseSat[er][esr];     

    arrayView1d<real64> const & pres = m_pressure[er][esr];
    arrayView1d<real64> const & sat  = m_wettingPhaseSat[er][esr];

    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex const ei )
    {
      pres[ei] += dPres[ei];
      sat[ei]  += dSat[ei]; 
    });
  });
}


void TwoPhaseBase::AssembleSystem( real64 const time_n, 
                                   real64 const dt,
                                   DomainPartition * const domain,
                                   DofManager const & dofManager,
                                   ParallelMatrix & matrix,
                                   ParallelVector & rhs )
{
  GEOSX_MARK_FUNCTION;

  matrix.zero();
  rhs.zero();

  matrix.open();
  rhs.open();

  AssembleAccumulationTerms( domain, &dofManager, &matrix, &rhs );

  AssembleFluxTerms( time_n, dt, domain, &dofManager, &matrix, &rhs );

  matrix.close();
  rhs.close();
  
  // Debug for logLevel >= 2
  GEOSX_LOG_LEVEL_RANK_0( 2, "After TwoPhaseBase::AssembleSystem" );
  GEOSX_LOG_LEVEL_RANK_0( 2, "\nJacobian:\n" << matrix );
  GEOSX_LOG_LEVEL_RANK_0( 2, "\nResidual:\n" << rhs );

  if( getLogLevel() >= 3 )
  {
    integer newtonIter = m_nonlinearSolverParameters.m_numNewtonIterations;

    string filename_mat = "matrix_" + std::to_string( time_n ) + "_" + std::to_string( newtonIter ) + ".mtx";
    matrix.write( filename_mat, true );

    string filename_rhs = "rhs_" + std::to_string( time_n ) + "_" + std::to_string( newtonIter ) + ".mtx";
    rhs.write( filename_rhs, true );

    GEOSX_LOG_RANK_0( "After TwoPhaseBase::AssembleSystem" );
    GEOSX_LOG_RANK_0( "Jacobian: written to " << filename_mat );
    GEOSX_LOG_RANK_0( "Residual: written to " << filename_rhs );
  }

}

void TwoPhaseBase::AssembleAccumulationTerms( DomainPartition const * const domain,
                                              DofManager const * const dofManager,
                                              ParallelMatrix * const matrix,
                                              ParallelVector * const rhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );

  localIndex constexpr numPhases = NUM_PHASES;
  localIndex constexpr numDof    = NUM_DOF;

  string const dofKey = dofManager->getKey( viewKeyStruct::elemDofFieldString );

  applyToSubRegions( mesh, [&] ( localIndex const er, localIndex const esr,
                                 ElementRegionBase const * const GEOSX_UNUSED_ARG( region ),
                                 ElementSubRegionBase const * const subRegion )
  {
    arrayView1d<globalIndex const> const & dofNumber = subRegion->getReference< array1d<globalIndex> >( dofKey );

    arrayView1d<integer     const> const & elemGhostRank = m_elemGhostRank[er][esr];

    arrayView1d<real64 const> const & volume           = m_volume[er][esr];
    arrayView1d<real64 const> const & porosityRef      = m_porosityRef[er][esr];
    arrayView2d<real64 const> const & pvMult           = m_pvMult[er][esr][m_solidIndex];
    arrayView2d<real64 const> const & dPvMult_dPres    = m_dPvMult_dPres[er][esr][m_solidIndex];
    
    arrayView1d<real64 const> const & sat              = m_wettingPhaseSat[er][esr];
    arrayView1d<real64 const> const & dSat             = m_deltaWettingPhaseSat[er][esr];

    arrayView3d<real64 const> const & phaseDens        = m_phaseDens[er][esr][m_fluidIndex];
    arrayView3d<real64 const> const & dPhaseDens_dPres = m_dPhaseDens_dPres[er][esr][m_fluidIndex];

    arrayView1d<real64 const> const & porosityOld      = m_porosityOld[er][esr];
    arrayView2d<real64 const> const & phaseDensOld     = m_phaseDensOld[er][esr];

    localIndex const dp = ColOffset::DPRES;
    localIndex const dS = ColOffset::DSAT;

    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      if (elemGhostRank[ei] < 0)
      {

        stackArray1d<globalIndex, numPhases>     eqnRowIndices( numPhases );
        stackArray1d<globalIndex, numDof>        dofColIndices( numDof );
        stackArray1d<real64, numPhases>          localAccum( numPhases );
        stackArray2d<real64, numPhases * numDof> localAccumJacobian( numPhases, numDof );

        real64 const poreVolNew        = volume[ei] * porosityRef[ei] * pvMult[ei][0];
        real64 const dPoreVolNew_dPres = volume[ei] * porosityRef[ei] * dPvMult_dPres[ei][0];
        real64 const poreVolOld        = volume[ei] * porosityOld[ei];

        stackArray1d<real64, numPhases> satNew( numPhases );
        stackArray1d<real64, numPhases> satOld( numPhases );
        stackArray1d<real64, numPhases> dSatNew_dS( numPhases );
        satNew[m_ipw]      = sat[ei] + dSat[ei];
        satNew[m_ipnw]     = 1 - satNew[m_ipw];
        satOld[m_ipw]      = sat[ei];
        satOld[m_ipnw]     = 1 - satOld[m_ipw];
        dSatNew_dS[m_ipw]  = 1;
        dSatNew_dS[m_ipnw] = -1;

        // dof numbers
        dofColIndices[dp] = dofNumber[ei] + dp;
        dofColIndices[dS] = dofNumber[ei] + dS;
        
        for (localIndex ip = 0; ip < numPhases; ++ip)
        {
          localIndex const rowId = m_phaseToRow[ip];
          eqnRowIndices[rowId] = dofNumber[ei] + m_phaseToRow[ip];
          
          // residual
          localAccum[rowId]  = poreVolNew * phaseDens[ei][0][ip] * satNew[ip];
          localAccum[rowId] -= poreVolOld * phaseDensOld[ei][ip] * satOld[ip];
        
          // jacobian
          localAccumJacobian[rowId][dp] = ( dPoreVolNew_dPres * phaseDens[ei][0][ip]   
                                          + poreVolNew        * dPhaseDens_dPres[ei][0][ip] ) * satNew[ip];
          localAccumJacobian[rowId][dS] =   poreVolNew        * phaseDens[ei][0][ip]          * dSatNew_dS[ip];
        }
          
        // add contribution to global residual and jacobian
        
        rhs->add( eqnRowIndices.data(),
                  localAccum.data(),
                  numPhases );

        matrix->add( eqnRowIndices.data(),
                     dofColIndices.data(),
                     localAccumJacobian.data(),
                     numPhases, numDof );
        
      }
    });
  });
  
}


bool TwoPhaseBase::CheckSystemSolution( DomainPartition const * const domain,
                                        DofManager const & dofManager,
                                        ParallelVector const & solution,
                                        real64 const scalingFactor )
{
  MeshLevel const * const mesh = domain->getMeshBody(0)->getMeshLevel(0);
  real64 const * localSolution = solution.extractLocalVector();
  int localCheck = 1;

  string const elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString );

  applyToSubRegions( mesh, [&] ( localIndex const er, localIndex const esr,
                                 ElementRegionBase const * const GEOSX_UNUSED_ARG( region ),
                                 ElementSubRegionBase const * const subRegion )
  {
    arrayView1d<globalIndex const> const & elemDofNumber =
      subRegion->getReference< array1d<globalIndex> >( elemDofKey );

    arrayView1d<integer const> const & elemGhostRank = m_elemGhostRank[er][esr];

    arrayView1d<real64 const> const & pres  = m_pressure[er][esr];
    arrayView1d<real64 const> const & dPres = m_deltaPressure[er][esr];
    arrayView1d<real64 const> const & sat   = m_wettingPhaseSat[er][esr];
    arrayView1d<real64 const> const & dSat  = m_deltaWettingPhaseSat[er][esr];

    forall_in_range<serialPolicy>( 0, subRegion->size(), [&] ( localIndex ei )
    {
      if (elemGhostRank[ei] >= 0)
      {
        return;
      }

      // extract solution and apply to dP
      {
        localIndex const lid = solution.getLocalRowID( elemDofNumber[ei] + ColOffset::DPRES );
        real64 const newPres = pres[ei] + dPres[ei] + scalingFactor * localSolution[lid];
        if (newPres < 0.0)
        {
          localCheck = 0;
        }
      }

      // extract solution and apply to dS 
      {
        localIndex const lid = solution.getLocalRowID( elemDofNumber[ei] + ColOffset::DSAT );
        real64 const newSat  = sat[ei] + dSat[ei] + scalingFactor * localSolution[lid];
        if (newSat < 0.0 || newSat > 1.0 )
        {
          localCheck = 0;
        }
      }
    });
  });
  int globalCheck;

  MpiWrapper::allReduce( &localCheck,
                         &globalCheck,
                         1,
                         MPI_MIN,
                         MPI_COMM_GEOSX );

  bool result = true;
  if (globalCheck == 0)
  {
    result = false;
  }
  return result;
}


void TwoPhaseBase::SolveSystem( DofManager const & dofManager,
                                ParallelMatrix & matrix,
                                ParallelVector & rhs,
                                ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  rhs.scale( -1.0 );
  solution.zero();

  SolverBase::SolveSystem( dofManager, matrix, rhs, solution );
  
  // Debug for logLevel >= 2
  GEOSX_LOG_LEVEL_RANK_0( 2, "After TwoPhaseBase::SolveSystem" );
  GEOSX_LOG_LEVEL_RANK_0( 2, "\nSolution:\n" << solution );

}

void TwoPhaseBase::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  MeshLevel * mesh = domain->getMeshBody(0)->getMeshLevel(0);

  applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegionBase * const GEOSX_UNUSED_ARG( region ),
                                 ElementSubRegionBase * const subRegion )
  {
    arrayView1d<real64> const & dPres = m_deltaPressure[er][esr];
    arrayView1d<real64> const & dSat  = m_deltaWettingPhaseSat[er][esr];

    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      dPres[ei] = 0.0;
      dSat[ei]  = 0.0; 
    } );

    UpdateState( subRegion );
  } );
}
  
void TwoPhaseBase::ResetViews( DomainPartition * const domain )
{
  FlowSolverBase::ResetViews( domain );

  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = mesh->getElemManager();
  ConstitutiveManager * const constitutiveManager = domain->getConstitutiveManager();


  // primary variables
  m_pressure =
    elemManager->ConstructViewAccessor< array1d<real64>, arrayView1d<real64> >( viewKeyStruct::pressureString );
  m_deltaPressure =
    elemManager->ConstructViewAccessor< array1d<real64>, arrayView1d<real64> >( viewKeyStruct::deltaPressureString );
  m_wettingPhaseSat =
    elemManager->ConstructViewAccessor< array1d<real64>, arrayView1d<real64> >( viewKeyStruct::wettingPhaseSatString );
  m_deltaWettingPhaseSat =
    elemManager->ConstructViewAccessor< array1d<real64>, arrayView1d<real64> >( viewKeyStruct::deltaWettingPhaseSatString );

  // auxiliary data
  m_phaseMob =
    elemManager->ConstructViewAccessor< array2d<real64>, arrayView2d<real64> >( viewKeyStruct::phaseMobilityString );
  m_dPhaseMob_dPres =
    elemManager->ConstructViewAccessor< array2d<real64>, arrayView2d<real64> >( viewKeyStruct::dPhaseMobility_dPressureString );
  m_dPhaseMob_dSat =
    elemManager->ConstructViewAccessor< array2d<real64>, arrayView2d<real64> >( viewKeyStruct::dPhaseMobility_dSaturationString );

  m_pvMult =
    elemManager->ConstructFullMaterialViewAccessor< array2d<real64>, arrayView2d<real64> >( ConstitutiveBase::viewKeyStruct::poreVolumeMultiplierString,
                                                                                            constitutiveManager );
  m_dPvMult_dPres =
    elemManager->ConstructFullMaterialViewAccessor< array2d<real64>, arrayView2d<real64> >( ConstitutiveBase::viewKeyStruct::dPVMult_dPresString,
                                                                                            constitutiveManager );
  m_phaseDens =
    elemManager->ConstructFullMaterialViewAccessor< array3d<real64>, arrayView3d<real64> >( MultiFluidBase::viewKeyStruct::phaseDensityString,
                                                                                            constitutiveManager );
  m_dPhaseDens_dPres =
    elemManager->ConstructFullMaterialViewAccessor< array3d<real64>, arrayView3d<real64> >( MultiFluidBase::viewKeyStruct::dPhaseDensity_dPressureString,
                                                                                            constitutiveManager );

  // backup data
  m_porosityOld =
    elemManager->ConstructViewAccessor< array1d<real64>, arrayView1d<real64> >( viewKeyStruct::porosityOldString );
  m_phaseDensOld =
    elemManager->ConstructViewAccessor< array2d<real64>, arrayView2d<real64> >( viewKeyStruct::phaseDensityOldString );

}

void TwoPhaseBase::BackupFields( DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBody(0)->getMeshLevel(0);

  localIndex constexpr numPhases  = NUM_PHASES;
  
  // backup some fields used in time derivative approximation
  applyToSubRegions( mesh, [&] ( localIndex const er, localIndex const esr,
                                 ElementRegionBase * const GEOSX_UNUSED_ARG( region ),
                                 ElementSubRegionBase * const subRegion )
  {
    arrayView1d<integer const> const & elemGhostRank = m_elemGhostRank[er][esr];
    arrayView1d<real64  const> const & poroRef       = m_porosityRef[er][esr];
    arrayView2d<real64  const> const & pvMult        = m_pvMult[er][esr][m_solidIndex];
    
    arrayView3d<real64  const> const & phaseDens     = m_phaseDens[er][esr][m_fluidIndex];

    arrayView2d<real64> const & phaseDensOld         = m_phaseDensOld[er][esr];
    arrayView1d<real64> const & poroOld              = m_porosityOld[er][esr];
    
    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      if (elemGhostRank[ei] < 0)
      {
        
        poroOld[ei] = poroRef[ei] * pvMult[ei][0];      
        for (localIndex ip = 0; ip < numPhases; ++ip)
        {
          phaseDensOld[ei][ip] = phaseDens[ei][0][ip];
        }
        
      }
    });
  });
}

  

} /* namespace geosx */
