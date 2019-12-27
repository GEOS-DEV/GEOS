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
  m_numDofPerCell = 2;
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
      subRegion->registerWrapper< array2d<real64> >( viewKeyStruct::wettingPhaseSatString )->setPlotLevel(PlotLevel::LEVEL_0);
      subRegion->registerWrapper< array2d<real64> >( viewKeyStruct::deltaWettingPhaseSatString );

      // auxiliary variables
      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::porosityString );
      subRegion->registerWrapper< array2d<real64> >( viewKeyStruct::phaseMobilityString );
      subRegion->registerWrapper< array2d<real64> >( viewKeyStruct::dPhaseMobility_dPressureString );
      subRegion->registerWrapper< array2d<real64> >( viewKeyStruct::dPhaseMobility_dSaturationString );

      // backup fields
      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::porosityOldString );
      subRegion->registerWrapper< array2d<real64> >( viewKeyStruct::phaseDensityOldString );

    });

  }
}

void TwoPhaseBase::UpdateFluidModel(Group * const GEOSX_UNUSED_ARG( dataGroup ) ) const
{
  GEOSX_MARK_FUNCTION;
}

void TwoPhaseBase::UpdateSolidModel(Group * const GEOSX_UNUSED_ARG( dataGroup ) ) const
{
  GEOSX_MARK_FUNCTION;
}

void TwoPhaseBase::UpdateMobility( Group * const GEOSX_UNUSED_ARG( dataGroup ) ) const
{
  GEOSX_MARK_FUNCTION;
}


void TwoPhaseBase::UpdateState( Group * dataGroup ) const
{
  GEOSX_MARK_FUNCTION;

  UpdateFluidModel( dataGroup );
  UpdateSolidModel( dataGroup );
  UpdateMobility( dataGroup );
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

  // setup dof numbers and linear system
  SetupSystem( domain, m_dofManager, m_matrix, m_rhs, m_solution );

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
                                      DofManager & GEOSX_UNUSED_ARG( dofManager ),
                                      ParallelMatrix & GEOSX_UNUSED_ARG( matrix ),
                                      ParallelVector & GEOSX_UNUSED_ARG( rhs ),
                                      ParallelVector & GEOSX_UNUSED_ARG( solution ) )
{
  // bind the stored views to the current domain
  ResetViews( domain );

  // set deltas to zero and recompute dependent quantities
  ResetStateToBeginningOfStep( domain );

  // backup fields used in time derivative approximation
  BackupFields( domain );
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

  localIndex constexpr numPhases = 2;
  localIndex constexpr maxNumDof = numPhases + 1;

  localIndex const NP   = m_numPhases;
  localIndex const NDOF = m_numDofPerCell;

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

    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      if (elemGhostRank[ei] < 0)
      {
        
        stackArray1d<globalIndex, maxNumDof>        localAccumDOF( NDOF );
        stackArray1d<real64, numPhases>             localAccum( NP );
        stackArray2d<real64, numPhases * maxNumDof> localAccumJacobian( NP, NDOF );

        real64 const poreVolNew        = volume[ei] * porosityRef[ei] * pvMult[ei][0];
        real64 const dPoreVolNew_dPres = volume[ei] * porosityRef[ei] * dPvMult_dPres[ei][0];
        real64 const poreVolOld        = volume[ei] * porosityOld[ei];

        real64 const satNew            = sat[ei] + dSat[ei];
        real64 const satOld            = sat[ei];
        
        // dof numbers
        
        localAccumDOF[ColOffset::DPRES] = dofNumber[ei] + ColOffset::DPRES;
        localAccumDOF[ColOffset::DSAT]  = dofNumber[ei] + ColOffset::DSAT;
        
        // residual

        // wetting phase 
        localAccum[m_wettingPh]     = poreVolNew * phaseDens[ei][0][m_wettingPh] * satNew;
        localAccum[m_wettingPh]    -= poreVolOld * phaseDensOld[ei][m_wettingPh] * satOld;
        // non-wetting phase
        localAccum[m_nonWettingPh]  = poreVolNew * phaseDens[ei][0][m_nonWettingPh] * (1-satNew);
        localAccum[m_nonWettingPh] -= poreVolOld * phaseDensOld[ei][m_nonWettingPh] * (1-satOld);
        
        // jacobian

        // wetting phase
        localAccumJacobian[m_wettingPh][ColOffset::DPRES]    = ( dPoreVolNew_dPres * phaseDens[ei][0][m_wettingPh]   
                                                               + poreVolNew        * dPhaseDens_dPres[ei][0][m_wettingPh] ) * satNew;
        localAccumJacobian[m_wettingPh][ColOffset::DSAT]     =   poreVolNew        * phaseDens[ei][0][m_wettingPh];
        // non-wetting phase
        localAccumJacobian[m_nonWettingPh][ColOffset::DPRES] = ( dPoreVolNew_dPres * phaseDens[ei][0][m_nonWettingPh]   
                                                               + poreVolNew        * dPhaseDens_dPres[ei][0][m_nonWettingPh] ) * (1-satNew);
        localAccumJacobian[m_nonWettingPh][ColOffset::DSAT]  = - poreVolNew        * phaseDens[ei][0][m_wettingPh];

        // add contribution to global residual and jacobian
        
        rhs->add( localAccumDOF.data(),
                  localAccum.data(),
                  NP );

        matrix->add( localAccumDOF.data(),
                     localAccumDOF.data(),
                     localAccumJacobian.data(),
                     NP, NDOF );
        
      }
    });
  });
  
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

  m_phaseDens =
    elemManager->ConstructFullMaterialViewAccessor< array3d<real64>, arrayView3d<real64> >( MultiFluidBase::viewKeyStruct::phaseDensityString,
                                                                                            constitutiveManager );
  m_dPhaseDens_dPres =
    elemManager->ConstructFullMaterialViewAccessor< array3d<real64>, arrayView3d<real64> >( MultiFluidBase::viewKeyStruct::dPhaseDensity_dPressureString,
                                                                                            constitutiveManager );
  
}

void TwoPhaseBase::BackupFields( DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBody(0)->getMeshLevel(0);

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
        for (localIndex ip = 0; ip < m_numPhases; ++ip)
        {
          phaseDensOld[ei][ip] = phaseDens[ei][0][ip];
        }
        
      }
    });
  });
}

  

} /* namespace geosx */
