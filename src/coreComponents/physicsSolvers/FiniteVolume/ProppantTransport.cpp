/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file ProppantTransport.cpp
 */

#include "ProppantTransport.hpp"

#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/Fluid/SlurryFluidBase.hpp"
#include "constitutive/Fluid/ParticleFluidBase.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "MPI_Communications/CommunicationTools.hpp"
#include "systemSolverInterface/LinearSolverWrapper.hpp"
#include "systemSolverInterface/EpetraBlockSystem.hpp"
#include "MPI_Communications/NeighborCommunicator.hpp"
#include "mesh/FaceElementRegion.hpp"

/**
 * @namespace the geosx namespace that encapsulates the majority of the code
 */
namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace systemSolverInterface;

ProppantTransport::ProppantTransport( const std::string& name,
                                  Group * const parent ):
  FlowSolverBase(name, parent)
{
   m_numDofPerCell = 2;

  this->registerWrapper( viewKeyStruct::proppantNameString,  &m_proppantName,  false )->setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of proppant constitutive object to use for this solver.");

  this->registerWrapper( viewKeyStruct::proppantIndexString, &m_proppantIndex, false );

  registerWrapper( viewKeyStruct::updatePermeabilityString, &m_updatePermeability, false )->setApplyDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Flag that enables/disables real-time fracture permeability update");

  registerWrapper( viewKeyStruct::updateProppantMobilityString, &m_updateProppantMobility, false )->setApplyDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Flag that enables/disables real-time proppant mobility update");
  
}

void ProppantTransport::RegisterDataOnMesh(Group * const MeshBodies)
{
  FlowSolverBase::RegisterDataOnMesh(MeshBodies);

  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    MeshLevel * meshLevel = Group::group_cast<MeshBody *>(mesh.second)->getMeshLevel(0);
    ElementRegionManager * const elemManager = meshLevel->getElemManager();

    elemManager->forElementSubRegions<FaceElementSubRegion>( [&]( FaceElementSubRegion * const subRegion )
    {
      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::pressureString )->
        setPlotLevel(PlotLevel::LEVEL_0);
      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::deltaPressureString );
      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::proppantConcentrationString )->
        setPlotLevel(PlotLevel::LEVEL_0);
      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::deltaProppantConcentrationString );      
      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::deltaVolumeString );
      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::densityString );
      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::oldProppantConcentrationString );
      subRegion->registerWrapper< array1d<globalIndex> >( viewKeyStruct::blockLocalDofNumberString );
      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::porosityString )->
        setPlotLevel(PlotLevel::LEVEL_1);
    } );

  }
}

void ProppantTransport::InitializePreSubGroups(Group * const rootGroup)
{
  FlowSolverBase::InitializePreSubGroups(rootGroup);

  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);
  
  ConstitutiveManager * const cm = domain->getConstitutiveManager();

  ConstitutiveBase const * proppant  = cm->GetConstitutiveRelation<ConstitutiveBase>( m_proppantName );
  GEOS_ERROR_IF( proppant == nullptr, "Proppant model " + m_proppantName + " not found" );
  m_proppantIndex = proppant->getIndexInParent();
  
}

void ProppantTransport::UpdateFluidModel(Group * const dataGroup)
{
  GEOSX_MARK_FUNCTION;

  SlurryFluidBase * const fluid = GetConstitutiveModel<SlurryFluidBase>( dataGroup, m_fluidName );

  arrayView1d<real64 const> const & pres = dataGroup->getReference<array1d<real64>>( viewKeyStruct::pressureString );
  arrayView1d<real64 const> const & dPres = dataGroup->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

  arrayView1d<real64 const> const & conc = dataGroup->getReference<array1d<real64>>( viewKeyStruct::proppantConcentrationString );
  arrayView1d<real64 const> const & dConc = dataGroup->getReference<array1d<real64>>( viewKeyStruct::deltaProppantConcentrationString );  


  // TODO replace with batch update (need up-to-date pressure and temperature fields)
  forall_in_range<RAJA::seq_exec>( 0, dataGroup->size(), GEOSX_LAMBDA ( localIndex const a )
  {
    fluid->PointUpdate( pres[a] + dPres[a], conc[a] + dConc[a], a, 0 );
  });
  //fluid->BatchUpdate( pres, temp, compFrac );
}

void ProppantTransport::UpdateProppantModel(Group * const dataGroup)
{
  GEOSX_MARK_FUNCTION;

  ParticleFluidBase * const particle = GetConstitutiveModel<ParticleFluidBase>( dataGroup, m_proppantName );

  arrayView1d<real64 const> const & conc = dataGroup->getReference<array1d<real64>>( viewKeyStruct::proppantConcentrationString );
  arrayView1d<real64 const> const & dConc = dataGroup->getReference<array1d<real64>>( viewKeyStruct::deltaProppantConcentrationString );  


  // TODO replace with batch update (need up-to-date pressure and temperature fields)
  forall_in_range<RAJA::seq_exec>( 0, dataGroup->size(), GEOSX_LAMBDA ( localIndex const a )
  {
    particle->PointUpdate( conc[a] + dConc[a], a);
  });

}

void ProppantTransport::UpdateProppantModelStep(Group * const dataGroup)
{
  GEOSX_MARK_FUNCTION;

  ParticleFluidBase * const particle = GetConstitutiveModel<ParticleFluidBase>( dataGroup, m_proppantName );

  arrayView1d<real64 const> const & conc = dataGroup->getReference<array1d<real64>>( viewKeyStruct::proppantConcentrationString );

  arrayView1d<real64 const> const & aperture = dataGroup->getReference<array1d<real64>>( FaceElementSubRegion::viewKeyStruct::elementApertureString );  

  // TODO replace with batch update (need up-to-date pressure and temperature fields)
  forall_in_range<RAJA::seq_exec>( 0, dataGroup->size(), GEOSX_LAMBDA ( localIndex const a )
  {
    particle->PointUpdateMob( conc[a], aperture[a], a);
  });

}

void ProppantTransport::UpdateState( Group * dataGroup )
{
  GEOSX_MARK_FUNCTION;

  UpdateFluidModel( dataGroup );

  UpdateProppantModel( dataGroup );  

}

void ProppantTransport::InitializePostInitialConditions_PreSubGroups( Group * const rootGroup )
{
  GEOSX_MARK_FUNCTION;

  FlowSolverBase::InitializePostInitialConditions_PreSubGroups( rootGroup );

  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);
  MeshLevel * mesh = domain->getMeshBody(0)->getMeshLevel(0);

  //TODO this is a hack until the sets are fixed to include ghosts!!
  std::map<string, string_array > fieldNames;
  fieldNames["elems"].push_back( viewKeyStruct::pressureString );
  fieldNames["elems"].push_back( viewKeyStruct::proppantConcentrationString );  

  array1d<NeighborCommunicator> & comms =
    domain->getReference< array1d<NeighborCommunicator>>( domain->viewKeys.neighbors );

  CommunicationTools::SynchronizeFields( fieldNames, mesh, comms );

  ResetViews( domain );

  // We have to redo the below loop after fractures are generated
  
  applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegionBase * const region,
                                 ElementSubRegionBase * const subRegion )
  {
    UpdateState( subRegion );

    arrayView2d<real64 const> const & dens    = m_density[er][esr][m_fluidIndex];
    arrayView1d<real64> const & densOld = m_densityOld[er][esr];

    arrayView1d<real64> const & conc = m_proppantConcentration[er][esr];
    arrayView1d<real64> const & concOld = m_proppantConcentrationOld[er][esr];    

    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
        densOld[ei] = dens[ei][0];
        concOld[ei] = conc[ei];
    });

  } );
}

real64 ProppantTransport::SolverStep( real64 const& time_n,
                                    real64 const& dt,
                                    const int cycleNumber,
                                    DomainPartition * domain )
{
  GEOSX_MARK_FUNCTION;

  FlowSolverBase::PrecomputeData(domain);

  MeshLevel * mesh = domain->getMeshBody(0)->getMeshLevel(0);

  NodeManager const * const nodeManager = mesh->getNodeManager();
  FaceManager const * const faceManager = mesh->getFaceManager();
  
  real64 dt_return = dt;

  ImplicitStepSetup( time_n,
                     dt,
                     domain,
                     m_dofManager,
                     m_matrix,
                     m_rhs,
                     m_solution );

  if(cycleNumber == 1) {

    /*  assign intitial and boundary conditions */
    FieldSpecificationManager const * boundaryConditionManager = FieldSpecificationManager::get();

    boundaryConditionManager->ApplyInitialConditions( domain );

    /* Below must be called after ImplicitStepSetup */

    applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegionBase * const region,
                                 ElementSubRegionBase * const subRegion )
    {

      arrayView1d<real64> const & dVol    = m_deltaVolume[er][esr];
      arrayView1d<real64> const & dPres   = m_deltaPressure[er][esr];
      arrayView1d<real64> const & dConc   = m_deltaProppantConcentration[er][esr];        

      arrayView2d<real64 const> const & dens = m_density[er][esr][m_fluidIndex];
      arrayView1d<real64> const & densOld = m_densityOld[er][esr];

      arrayView1d<real64> const & conc = m_proppantConcentration[er][esr];
      arrayView1d<real64> const & concOld = m_proppantConcentrationOld[er][esr];

      subRegion->CalculateElementGeometricQuantities( *nodeManager,
                                                      *faceManager );
      
      forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
      {
        dPres[ei] = 0.0;
        dConc[ei] = 0.0;
        dVol[ei] = 0.0;
        concOld[ei] = conc[ei];
      } );

      UpdateState( subRegion );

      forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
      {
        densOld[ei] = dens[ei][0];
      } );
    });
  }


  if(m_updateProppantMobility)
    {
      applyToSubRegions( mesh, [&] ( ElementSubRegionBase * const subRegion )
     {
       UpdateProppantModelStep( subRegion );
     });

    }
      
  // currently the only method is implicit time integration
  dt_return= this->NonlinearImplicitStep( time_n,
                                          dt,
                                          cycleNumber,
                                          domain,
                                          m_dofManager,
                                          m_matrix,
                                          m_rhs,
                                          m_solution );

  // final step for completion of timestep. typically secondary variable updates and cleanup.
  ImplicitStepComplete( time_n, dt_return, domain );

  return dt_return;

}


void ProppantTransport::ImplicitStepSetup( real64 const & time_n,
                                           real64 const & dt,
                                           DomainPartition * const domain,
                                           DofManager & dofManager,
                                           ParallelMatrix & matrix,
                                           ParallelVector & rhs,
                                           ParallelVector & solution )
{
  ResetViews( domain );

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  /* The loop below could be moved to SolverStep after ImplicitStepSetup */
  
  applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegionBase * const region,
                                 ElementSubRegionBase * const subRegion )
  {

    arrayView1d<real64> const & dVol    = m_deltaVolume[er][esr];
    arrayView1d<real64> const & dPres   = m_deltaPressure[er][esr];
    arrayView1d<real64> const & dConc   = m_deltaProppantConcentration[er][esr];        
    arrayView2d<real64 const> const & dens = m_density[er][esr][m_fluidIndex];
    arrayView1d<real64> const & densOld = m_densityOld[er][esr];

    arrayView1d<real64> const & conc = m_proppantConcentration[er][esr];
    arrayView1d<real64> const & concOld = m_proppantConcentrationOld[er][esr];        

    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      dPres[ei] = 0.0;
      dConc[ei] = 0.0;      
      dVol[ei] = 0.0;

      densOld[ei] = dens[ei][0];
      concOld[ei] = conc[ei];      
    } );
  } );

  // setup dof numbers and linear system
  SetupSystem( domain,
               m_dofManager,
               m_matrix,
               m_rhs,
               m_solution  );

}

void ProppantTransport::ImplicitStepComplete( real64 const & time_n,
                                            real64 const & dt,
                                            DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegionBase * const region,
                                 ElementSubRegionBase * const subRegion )
  {
    arrayView1d<real64> const & pres = m_pressure[er][esr];
    arrayView1d<real64> const & conc = m_proppantConcentration[er][esr];
    
    arrayView1d<real64> const & vol  = m_volume[er][esr];

    arrayView1d<real64 const> const & dPres = m_deltaPressure[er][esr];
    arrayView1d<real64 const> const & dConc = m_deltaProppantConcentration[er][esr];    
    arrayView1d<real64 const> const & dVol  = m_deltaVolume[er][esr];

    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      pres[ei] += dPres[ei];
      conc[ei] += dConc[ei];      
      vol[ei] += dVol[ei];
    } );

  } );

}

void ProppantTransport::SetupDofs( DomainPartition const * const domain,
                                   DofManager & dofManager ) const
{
  dofManager.addField( viewKeyStruct::pressureString,
                       DofManager::Location::Elem,
                       DofManager::Connectivity::Face,
                       m_numDofPerCell,
                       m_targetRegions );
}


void ProppantTransport::AssembleSystem( real64 const time,
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

  if( verboseLevel() >= 2 )
  {
    GEOS_LOG_RANK_0("After ProppantTransport::AssembleAccumulationTerms");
    GEOS_LOG_RANK_0("\nJacobian:\n");
    matrix.print(std::cout);
    GEOS_LOG_RANK_0("\nResidual:\n");
    rhs.print(std::cout);
  }

  AssembleFluxTerms( time,
                     dt,
                     domain,
                     &dofManager,
                     &matrix,
                     &rhs);

  matrix.close();
  rhs.close();

  if( verboseLevel() >= 2 )
  {
    GEOS_LOG_RANK_0("After ProppantTransport::AssembleSystem");
    GEOS_LOG_RANK_0("\nJacobian:\n");
    matrix.print(std::cout);
    GEOS_LOG_RANK_0("\nResidual:\n");
    rhs.print(std::cout);
  }
}

void ProppantTransport::AssembleAccumulationTerms( DomainPartition const * const domain,
                                                   DofManager const * const dofManager,
                                                   ParallelMatrix * const matrix,
                                                   ParallelVector * const rhs)
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  NodeManager const * const nodeManager = mesh->getNodeManager();
  EdgeManager const * const edgeManager = mesh->getEdgeManager();
  FaceManager const * const faceManager = mesh->getFaceManager();

  
  applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegionBase const * const region,
                                 ElementSubRegionBase const * const subRegion )
  {
    string const dofKey = dofManager->getKey( viewKeyStruct::pressureString );
    arrayView1d<globalIndex const> const & dofNumber = subRegion->getReference< array1d<globalIndex> >( dofKey );

    arrayView1d<integer const>     const & elemGhostRank = m_elemGhostRank[er][esr];

    arrayView1d<real64 const> const & densOld       = m_densityOld[er][esr];

    arrayView1d<real64 const> const & concOld       = m_proppantConcentrationOld[er][esr];    

    arrayView1d<real64 const> const & volume        = m_volume[er][esr];
    arrayView1d<real64 const> const & dVol          = m_deltaVolume[er][esr];
    arrayView2d<real64 const> const & dens          = m_density[er][esr][m_fluidIndex];

    arrayView1d<real64 const> const & conc          = m_proppantConcentration[er][esr];

    arrayView1d<real64 const> const & dConc          = m_deltaProppantConcentration[er][esr];    
    
    arrayView2d<real64 const> const & dDens_dPres   = m_dDens_dPres[er][esr][m_fluidIndex];
    arrayView2d<real64 const> const & dDens_dConc   = m_dDens_dConc[er][esr][m_fluidIndex];    

    arrayView1d<real64 const> const & dPres         = m_deltaPressure[er][esr];

    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      if (elemGhostRank[ei] < 0)
      {
        stackArray1d<globalIndex, MAX_NUM_COMPONENTS>         localAccumDOF( m_numDofPerCell );
        stackArray1d<real64, MAX_NUM_COMPONENTS>             localAccum( m_numDofPerCell );
        stackArray2d<real64, MAX_NUM_COMPONENTS * MAX_NUM_COMPONENTS> localAccumJacobian( m_numDofPerCell, m_numDofPerCell );

        real64 const densNew = dens[ei][0];
        real64 const concNew = conc[ei] + dConc[ei];


        // fluid-mixture mass conservation
        localAccum[0] = (densNew  - densOld[ei]) * volume[ei];

        // proppant mass conservation
        localAccum[1] = (concNew - concOld[ei]) * volume[ei];

        // Derivative of residual wrt to pressure and concentration in the cell
        localAccumJacobian[0][0] = dDens_dPres[ei][0] * volume[ei];
        localAccumJacobian[0][1] = dDens_dConc[ei][0] * volume[ei];

        localAccumJacobian[1][0] = 0.0;
        localAccumJacobian[1][1] = volume[ei];

        globalIndex const elemDOF = dofNumber[ei];
        for (localIndex idof = 0; idof < m_numDofPerCell; ++idof)
        {
          localAccumDOF[idof] = elemDOF + idof;
        }
        // add contribution to global residual and dRdP
        rhs->add( localAccumDOF.data(),
                  localAccum.data(),
                  m_numDofPerCell );

        matrix->add( localAccumDOF.data(),
                     localAccumDOF.data(),
                     localAccumJacobian.data(),
                     m_numDofPerCell, m_numDofPerCell );
      }
    } );
  } );
}


void ProppantTransport::AssembleFluxTerms( real64 const time_n,
                                           real64 const dt,
                                           DomainPartition const * const domain,
                                           DofManager const * const dofManager,
                                           ParallelMatrix * const matrix,
                                           ParallelVector * const rhs )
{
  GEOSX_MARK_FUNCTION;

   //  R1Tensor unitGravityVector = m_gravityVector.UnitVector(); not working 
  R1Tensor unitGravityVector = m_gravityVector;  
  unitGravityVector.Normalize();
  
  MeshLevel const * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  NodeManager const * const nodeManager = mesh->getNodeManager();
  EdgeManager const * const edgeManager = mesh->getEdgeManager();
  FaceManager const * const faceManager = mesh->getFaceManager();
  ElementRegionManager const * const elemManager = mesh->getElemManager();


  NumericalMethodsManager const * numericalMethodManager =
    domain->getParent()->GetGroup<NumericalMethodsManager>( keys::numericalMethodsManager );

  FiniteVolumeManager const * fvManager =
    numericalMethodManager->GetGroup<FiniteVolumeManager>( keys::finiteVolumeManager );

  FluxApproximationBase const * fluxApprox = fvManager->getFluxApproximation( m_discretizationName );

  string const dofKey = dofManager->getKey( viewKeyStruct::pressureString );
  ElementRegionManager::ElementViewAccessor< arrayView1d<globalIndex> >
  dofNumberAccessor = elemManager->ConstructViewAccessor< array1d<globalIndex>, arrayView1d<globalIndex> >( dofKey );

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  const & pres        = m_pressure;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  const & dPres       = m_deltaPressure;

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  const & conc       = m_proppantConcentration;

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  const & concOld       = m_proppantConcentrationOld;
  
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  const & dConc       = m_deltaProppantConcentration;
  
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  const & gravDepth   = m_gravDepth;

  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & dens        = m_density;

  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & dDens_dPres = m_dDens_dPres;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & dDens_dConc = m_dDens_dConc;  

  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & visc        = m_viscosity;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & dVisc_dPres = m_dVisc_dPres;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & dVisc_dConc = m_dVisc_dConc;  

  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & fluidDensity = m_fluidDensity;

  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & dFluidDens_dPres = m_dFluidDens_dPres;  
  
  ElementRegionManager::MaterialViewAccessor<arrayView1d<real64>> const & settlingFactor = m_settlingFactor;
  ElementRegionManager::MaterialViewAccessor<arrayView1d<real64>> const & dSettlingFactor_dConc = m_dSettlingFactor_dConc;

  ElementRegionManager::MaterialViewAccessor<arrayView1d<real64>> const & collisionFactor = m_collisionFactor;

  ElementRegionManager::MaterialViewAccessor<arrayView1d<real64>> const & dCollisionFactor_dConc = m_dCollisionFactor_dConc;

  ElementRegionManager::MaterialViewAccessor<arrayView1d<bool>> const & isProppantMobile = m_isProppantMobile;      

  ElementRegionManager::MaterialViewAccessor<arrayView1d<real64>> const & proppantPackPermeability = m_proppantPackPermeability;      
  
  arrayView1d<R1Tensor const> const & X = nodeManager->referencePosition();
  
  constexpr localIndex numElems = FluxApproximationBase::CellStencil::NUM_POINT_IN_FLUX;
  constexpr localIndex maxStencilSize = FluxApproximationBase::CellStencil::MAX_STENCIL_SIZE;

  constexpr  localIndex DOF1 = numElems * MAX_NUM_COMPONENTS;
  constexpr  localIndex DOF2 = maxStencilSize * MAX_NUM_COMPONENTS;

  localIndex DOF = numElems * m_numDofPerCell;
  
  const FluxApproximationBase::CellStencil & fractureStencil = fluxApprox->getReference<FluxApproximationBase::CellStencil>(FluxApproximationBase::viewKeyStruct::fractureStencilString);

  ArrayOfArraysView<FluxApproximationBase::CellStencil::Entry const, true> const & connections = fractureStencil.getConnections();

  const localIndex_array & fconnectorIndices = fractureStencil.getConnectorMeshIndices();

  FaceElementRegionBase * const fractureRegion = elemManager->GetRegion<FaceElementRegion>(m_targetRegions[0]);
  
  array1d<localIndex> const & fractureConnectorsToEdges = fractureRegion->getReference< array1d<localIndex > >( FaceElementRegion::viewKeyStruct::fractureConnectorsToEdgesString );

  FaceElementSubRegion const * fractureSubRegion = static_cast<const FaceElementSubRegion *>(fractureRegion->GetSubRegion(0));

  arrayView2d< localIndex const > const & elemsToFaces = fractureSubRegion->faceList();

  arrayView1d< real64 const > const & aperture = fractureSubRegion->getElementAperture();
    
  arrayView1d<R1Tensor const> const & faceCenters = faceManager->faceCenter();

  
  forall_in_range<stencilPolicy>( 0, connections.size(), GEOSX_LAMBDA ( localIndex iconn )
  {
    localIndex const stencilSize = connections.sizeOfArray(iconn);
    localIndex const edgeIndex = fractureConnectorsToEdges[fconnectorIndices[iconn] ];

    R1Tensor edgeCenter, edgeToFaceVec;
    edgeManager->calculateCenter( edgeIndex, X, edgeCenter );
    edgeManager->calculateLength( edgeIndex, X, edgeToFaceVec );

    real64 edgeLength = edgeToFaceVec.L2_Norm();


    // working arrays
    stackArray1d<globalIndex, DOF1> eqnRowIndices(DOF);
    stackArray1d<globalIndex, DOF2> dofColIndices(DOF);

    stackArray1d<double, DOF1> localFlux(DOF);
    stackArray2d<double, DOF1*DOF2> localFluxJacobian(DOF, DOF);

    // Need to update weight
    stackArray1d<real64, numElems> weight(numElems);

    // mixture density and fluid density in each face
    stackArray1d<real64, numElems> mixDens(numElems);
    stackArray1d<real64, numElems> dMixDens_dP(numElems);
    stackArray1d<real64, numElems> dMixDens_dC(numElems);

    stackArray1d<real64, numElems> fluidDens(numElems);
    stackArray1d<real64, numElems> dFluidDens_dP(numElems);


    // realted to slip velocity calculation
    stackArray1d<real64, numElems> settlingFac(numElems);
    stackArray1d<real64, numElems> dSettlingFac_dC(numElems);

    stackArray1d<real64, numElems> collisionFac(numElems);
    stackArray1d<real64, numElems> dCollisionFac_dC(numElems);

    stackArray1d<real64, numElems> transT(numElems);
    stackArray1d<real64, numElems> coefs(numElems);
    stackArray1d<bool, numElems> isProppantMob(numElems);

    real64 edgeDensity, edgeViscosity;
    stackArray1d<real64, numElems> dEdgeDens_dP(numElems);
    stackArray1d<real64, numElems> dEdgeDens_dC(numElems);

    stackArray1d<real64, numElems> dEdgeVisc_dP(numElems);
    stackArray1d<real64, numElems> dEdgeVisc_dC(numElems);

    stackArray1d<real64, numElems> dPe_dP(numElems);
    stackArray1d<real64, numElems> dPe_dC(numElems);
    stackArray1d<R1Tensor, numElems> Up(numElems);

    stackArray1d<real64, numElems> dCe_dP(numElems);
    stackArray1d<real64, numElems> dCe_dC(numElems);

    stackArray1d<real64, numElems> edgeToFaceFlux(numElems);
    stackArray2d<real64, numElems * numElems> dEdgeToFaceFlux_dP(numElems,numElems);
    stackArray2d<real64, numElems * numElems> dEdgeToFaceFlux_dC(numElems,numElems);

    stackArray1d<real64, numElems> edgeToFaceProppantFlux(numElems);
    stackArray2d<real64, numElems * numElems> dEdgeToFaceProppantFlux_dP(numElems,numElems);
    stackArray2d<real64, numElems * numElems> dEdgeToFaceProppantFlux_dC(numElems,numElems);

    stackArray1d<real64, numElems> P(numElems);
    stackArray1d<real64, numElems> C(numElems);

    // clear working arrays
    eqnRowIndices = -1;

    edgeDensity = 0.0;
    edgeViscosity = 0.0;

    dEdgeDens_dP = 0.0;
    dEdgeDens_dC = 0.0;

    dEdgeVisc_dP = 0.0;
    dEdgeVisc_dC = 0.0;

    dEdgeToFaceProppantFlux_dP = 0.0;
    dEdgeToFaceProppantFlux_dC = 0.0;

    // TODO: Need to calculate weight from FractureStencil
    weight = 0.5;

    //numElems == stencilSize;

    //get averaged edgeDensity and edgeViscosity

    for (localIndex i = 0; i < numElems; ++i) {

      FluxApproximationBase::CellStencil::Entry const & entry = connections(iconn, i);
      localIndex const er  = entry.index.region;
      localIndex const esr = entry.index.subRegion;
      localIndex const ei  = entry.index.index;


      for (localIndex j = 0; j < m_numDofPerCell; ++j)
      {

        eqnRowIndices[i * m_numDofPerCell + j] = dofNumber[er][esr][ei] * m_numDofPerCell + j;
        dofColIndices[i * m_numDofPerCell + j] = dofNumber[er][esr][ei] * m_numDofPerCell + j;

      }

      edgeDensity += weight[i] * dens[er][esr][m_fluidIndex][ei][0];
      dEdgeDens_dP[i] = weight[i] * dDens_dPres[er][esr][m_fluidIndex][ei][0];
      dEdgeDens_dC[i] = weight[i] * dDens_dConc[er][esr][m_fluidIndex][ei][0];

      edgeViscosity += weight[i] * visc[er][esr][m_fluidIndex][ei][0];
      dEdgeVisc_dP[i] = weight[i] * dVisc_dPres[er][esr][m_fluidIndex][ei][0];
      dEdgeVisc_dC[i] = weight[i] * dVisc_dConc[er][esr][m_fluidIndex][ei][0];

      P[i] = pres[er][esr][ei] + dPres[er][esr][ei];
      C[i] = conc[er][esr][ei] + dConc[er][esr][ei];

      mixDens[i] = dens[er][esr][m_fluidIndex][ei][0];
      dMixDens_dP[i] = dDens_dPres[er][esr][m_fluidIndex][ei][0];
      dMixDens_dC[i] = dDens_dConc[er][esr][m_fluidIndex][ei][0];

      fluidDens[i] = fluidDensity[er][esr][m_fluidIndex][ei][0];
      dFluidDens_dP[i] = dFluidDens_dPres[er][esr][m_fluidIndex][ei][0];

      settlingFac[i] = settlingFactor[er][esr][m_proppantIndex][ei];
      dSettlingFac_dC[i] = dSettlingFactor_dConc[er][esr][m_proppantIndex][ei];

      collisionFac[i] = collisionFactor[er][esr][m_proppantIndex][ei];
      dCollisionFac_dC[i] = dCollisionFactor_dConc[er][esr][m_proppantIndex][ei];

      isProppantMob[i] = isProppantMobile[er][esr][m_proppantIndex][ei];

      transT[i] = fabs(entry.weight);

      localIndex const faceIndex = elemsToFaces[ei][0];

      edgeToFaceVec = faceCenters[faceIndex];
      edgeToFaceVec -= edgeCenter;

      if(proppantPackPermeability[er][esr][m_proppantIndex][ei] > 0.0 && m_updatePermeability)
      {
        transT[i] = transT[i] * (1.0 - concOld[er][esr][ei]) + proppantPackPermeability[er][esr][m_proppantIndex][ei]  * edgeLength * aperture[ei] / edgeToFaceVec.L2_Norm() * concOld[er][esr][ei];

      }

      edgeToFaceVec.Normalize();

      coefs[i] = Dot(edgeToFaceVec, unitGravityVector) * edgeLength * aperture[ei];

    }


    real64 transTSum = 0.0;
    real64 Pe = 0.0;
    dPe_dP = 0.0;
    dPe_dC = 0.0;

    for (localIndex i = 0; i < numElems; ++i)
    {
      FluxApproximationBase::CellStencil::Entry const & entry = connections(iconn, i);

      localIndex const er  = entry.index.region;
      localIndex const esr = entry.index.subRegion;
      localIndex const ei  = entry.index.index;

      real64 const gravD    = gravDepth[er][esr][ei];
      real64 const gravTerm = edgeDensity * gravD;

      Pe += transT[i] * (P[i] - gravTerm);

      transTSum += transT[i];

      dPe_dP[i] += transT[i];

      for (localIndex j = 0; j < numElems; ++j) {
        dPe_dP[j] += -transT[i] * gravD * dEdgeDens_dP[j];
        dPe_dC[j] += -transT[i] * gravD * dEdgeDens_dC[j];
      }

    }

    for (localIndex i = 0; i < numElems; ++i)
    {

      dPe_dP[i] /= transTSum;
      dPe_dC[i] /= transTSum;

    }

    Pe /= transTSum;

    for (localIndex i = 0; i < numElems; ++i)
    {

      FluxApproximationBase::CellStencil::Entry const & entry = connections(iconn, i);

      localIndex const er  = entry.index.region;
      localIndex const esr = entry.index.subRegion;
      localIndex const ei  = entry.index.index;

      real64 const gravD    = gravDepth[er][esr][ei];
      real64 const gravTerm = edgeDensity * gravD;

      real64 const fluxTerm = Pe - (P[i] - gravTerm);

      edgeToFaceFlux[i] = transT[i] * fluxTerm / edgeViscosity;

      dEdgeToFaceFlux_dP[i][i] += -transT[i] / edgeViscosity;

      for (localIndex j = 0; j < numElems; ++j)
      {

        dEdgeToFaceFlux_dP[i][j] += -transT[i] * fluxTerm * dEdgeVisc_dP[j] / (edgeViscosity * edgeViscosity) + transT[i] * (dPe_dP[j] + dEdgeDens_dP[j] * gravD) / edgeViscosity;

        dEdgeToFaceFlux_dC[i][j] += -transT[i] * fluxTerm * dEdgeVisc_dC[j] / (edgeViscosity * edgeViscosity) + transT[i] * (dPe_dC[j] + dEdgeDens_dC[j] * gravD) / edgeViscosity;


      }

      //contribution from paricle slipping

      edgeToFaceProppantFlux[i] = (1.0 + fluidDens[i] / mixDens[i] * (1.0 - C[i]) * collisionFac[i]) *  edgeToFaceFlux[i] + fluidDens[i] / mixDens[i] * (1.0 - C[i]) * settlingFac[i] * coefs[i];

      dEdgeToFaceProppantFlux_dP[i][i] = (dFluidDens_dP[i] / mixDens[i] - fluidDens[i] * dMixDens_dP[i] / (mixDens[i] * mixDens[i])) * (1.0 - C[i]) * (collisionFac[i] * edgeToFaceFlux[i] + settlingFac[i] * coefs[i]);

      dEdgeToFaceProppantFlux_dC[i][i] = -fluidDens[i] / mixDens[i] * (1.0 + dMixDens_dC[i] / mixDens[i] * (1.0 - C[i])) * (collisionFac[i] *  edgeToFaceFlux[i] + settlingFac[i] * coefs[i]) + fluidDens[i] / mixDens[i] * (1.0 - C[i]) * (dCollisionFac_dC[i] *  edgeToFaceFlux[i] + dSettlingFac_dC[i] * coefs[i]);

      for (localIndex j = 0; j < numElems; ++j)
      {

        dEdgeToFaceProppantFlux_dP[i][j] += (1.0 + fluidDens[i] / mixDens[i] * (1.0 - C[i]) * collisionFac[i]) *  dEdgeToFaceFlux_dP[i][j];

        dEdgeToFaceProppantFlux_dC[i][j] += (1.0 + fluidDens[i] / mixDens[i] * (1.0 - C[i]) * collisionFac[i]) *  dEdgeToFaceFlux_dC[i][j];
      }

    }


    // get Ce

    real64 Ce = 0.0;
    dCe_dP = 0.0;
    dCe_dC = 0.0;

    real64 downStreamFlux = 0.0;
    stackArray1d<real64, numElems> dDownStreamFlux_dP(numElems);
    stackArray1d<real64, numElems> dDownStreamFlux_dC(numElems);

    dDownStreamFlux_dP = 0.0;
    dDownStreamFlux_dC = 0.0;

    for (localIndex i = 0; i < numElems; ++i)
    {

      if(!isProppantMob[i] && m_updateProppantMobility)
        continue;

      if(edgeToFaceProppantFlux[i] >= 0.0)
      {
        // downstream
        downStreamFlux += edgeToFaceProppantFlux[i];

        for(localIndex j = 0; j < numElems; ++j)
        {

          dDownStreamFlux_dP[j] += dEdgeToFaceProppantFlux_dP[i][j];
          dDownStreamFlux_dC[j] += dEdgeToFaceProppantFlux_dC[i][j];

        }

      }
      else
      {

        // upstream
        Ce += -edgeToFaceProppantFlux[i] * C[i];

        dCe_dC[i] += -edgeToFaceProppantFlux[i];

        for(localIndex j = 0; j < numElems; ++j)
        {

          dCe_dP[j] += -dEdgeToFaceProppantFlux_dP[i][j] * C[i];
          dCe_dC[j] += -dEdgeToFaceProppantFlux_dC[i][j] * C[i];

        }

      }

    }


    if(downStreamFlux > 0.0)
    {

      for (localIndex i = 0; i < numElems; ++i)
      {

        if(!isProppantMob[i] && m_updateProppantMobility)
          continue;

        dCe_dP[i] =  dCe_dP[i] / downStreamFlux - Ce * dDownStreamFlux_dP[i] / (downStreamFlux * downStreamFlux);
        dCe_dC[i] =  dCe_dC[i] / downStreamFlux - Ce * dDownStreamFlux_dC[i] / (downStreamFlux * downStreamFlux);;

      }

      Ce = Ce / downStreamFlux;

    }
    else
    {

      Ce = 0.0;
      for (localIndex i = 0; i < numElems; ++i)
      {

        if(!isProppantMob[i] && m_updateProppantMobility)
          continue;

        dCe_dP[i] =  0.0;
        dCe_dC[i] =  weight[i];
        Ce += C[i] * weight[i];

      }

    }

    for (localIndex i = 0; i < numElems; ++i)
    {

      localIndex idx1 = i * m_numDofPerCell;

      localFlux[idx1] = -edgeDensity * edgeToFaceFlux[i] * dt;

      for (localIndex j = 0; j < numElems; ++j)
      {

        localIndex idx2 = j * m_numDofPerCell;

        //dP

        localFluxJacobian[idx1][idx2]  = -(dEdgeDens_dP[j] * edgeToFaceFlux[i] + edgeDensity * dEdgeToFaceFlux_dP[i][j]) * dt;

        //dC

        localFluxJacobian[idx1][idx2 + 1]  = -(dEdgeDens_dC[j] * edgeToFaceFlux[i] + edgeDensity * dEdgeToFaceFlux_dC[i][j]) * dt;

      }


      if(isProppantMob[i] || !m_updateProppantMobility)
      {

        if(edgeToFaceProppantFlux[i] >= 0.0)
        {


          localFlux[idx1 + 1] = -Ce * edgeToFaceProppantFlux[i] * dt;
        }
        else
        {

          localFlux[idx1 + 1] = -C[i] * edgeToFaceProppantFlux[i] * dt;

        }

        for (localIndex j = 0; j < numElems; ++j)
        {

          localIndex idx2 = j * m_numDofPerCell;


          if(edgeToFaceProppantFlux[i] >= 0.0)
          {


            localFluxJacobian[idx1 + 1][idx2] = -(dCe_dP[j] * edgeToFaceProppantFlux[i] + Ce * dEdgeToFaceProppantFlux_dP[i][j]) * dt;

            localFluxJacobian[idx1 + 1][idx2 + 1] = -(dCe_dC[j] * edgeToFaceProppantFlux[i] + Ce * dEdgeToFaceProppantFlux_dC[i][j]) * dt;

          }
          else
          {


            localFluxJacobian[idx1 + 1][idx2] = -C[i] * dEdgeToFaceProppantFlux_dP[i][j] * dt;
            localFluxJacobian[idx1 + 1][idx2 + 1] = -C[i] * dEdgeToFaceProppantFlux_dC[i][j] * dt;

            if(i == j)
              localFluxJacobian[idx1 + 1][idx2 + 1] += -edgeToFaceProppantFlux[i] * dt;

          }

        }

      }
      else
      {

        localFlux[idx1 + 1] = 0.0;
        for (localIndex j = 0; j < numElems; ++j)
        {

          localIndex idx2 = j * m_numDofPerCell;
          localFluxJacobian[idx1 + 1][idx2] = 0.0;
          localFluxJacobian[idx1 + 1][idx2 + 1] = 0.0;
        }
      }
    }

    // Add to global residual/jacobian

    jacobian->SumIntoGlobalValues( integer_conversion<int>(DOF),
                                   eqnRowIndices.data(),
                                   integer_conversion<int>(DOF),
                                   eqnRowIndices.data(),
                                   localFluxJacobian.data(),
                                   Epetra_FECrsMatrix::ROW_MAJOR);


    residual->SumIntoGlobalValues( integer_conversion<int>(DOF), eqnRowIndices.data(), localFlux.data() );
    
  });
}

void ProppantTransport::ApplyBoundaryConditions( DomainPartition * const domain,
                                               EpetraBlockSystem * const blockSystem,
                                               real64 const time_n,
                                               real64 const dt )
{
  GEOSX_MARK_FUNCTION;


  FieldSpecificationManager * fsManager = FieldSpecificationManager::get();

  // call the BoundaryConditionManager::ApplyField function that will check to see
  // if the boundary condition should be applied to this subregion
  fsManager->Apply( time_n + dt, domain, "ElementRegions", "FLUX",
                    [&]( FieldSpecificationBase const * const fs,
                    string const &,
                    set<localIndex> const & lset,
                    Group * subRegion,
                    string const & ) -> void
  {
    arrayView1d<globalIndex const> const &
    dofNumber = subRegion->getReference< array1d<globalIndex> >( viewKeyStruct::blockLocalDofNumberString );

    fs->ApplyBoundaryConditionToSystem<FieldSpecificationAdd>( lset,
                                                               true,
                                                               time_n + dt,
                                                               dt,
                                                               subRegion,
                                                               dofNumber,
                                                               2,
                                                               blockSystem,
                                                               BlockIDs::proppantTransportBlock,
                                                               [&] (localIndex const a) -> real64
    {
      return 0;
    });
  });

  fsManager->Apply( time_n + dt, domain, "ElementRegions", viewKeyStruct::pressureString,
                    [&]( FieldSpecificationBase const * const fs,
                    string const &,
                    set<localIndex> const & lset,
                    Group * subRegion,
                    string const & ) -> void
  {
    arrayView1d<globalIndex const> const &
    dofNumber = subRegion->getReference< array1d<globalIndex> >( viewKeyStruct::blockLocalDofNumberString );

    //for now assume all the non-flux boundary conditions are Dirichlet type BC.

    arrayView1d<real64 const> const &
    pres = subRegion->getReference<array1d<real64> >( viewKeyStruct::pressureString );

    arrayView1d<real64 const> const &
    dPres = subRegion->getReference<array1d<real64> >( viewKeyStruct::deltaPressureString );

    // call the application of the boundary condition to alter the matrix and rhs
    fs->ApplyBoundaryConditionToSystem<FieldSpecificationEqual>( lset,
                                                                 false,
                                                                 time_n + dt,
                                                                 subRegion,
                                                                 dofNumber,
                                                                 2,
                                                                 blockSystem,
                                                                 BlockIDs::proppantTransportBlock,
                                                                 [&] (localIndex const a) -> real64
    {
      return pres[a] + dPres[a];
    });
  });


 fsManager->Apply( time_n + dt, domain, "ElementRegions", viewKeyStruct::proppantConcentrationString,
                    [&]( FieldSpecificationBase const * const fs,
                    string const &,
                    set<localIndex> const & lset,
                    Group * subRegion,
                    string const & ) -> void
  {
    arrayView1d<globalIndex const> const &
    dofNumber = subRegion->getReference< array1d<globalIndex> >( viewKeyStruct::blockLocalDofNumberString );

    //for now assume all the non-flux boundary conditions are Dirichlet type BC.

    arrayView1d<real64 const> const &
    conc = subRegion->getReference<array1d<real64> >( viewKeyStruct::proppantConcentrationString );

    arrayView1d<real64 const> const &
    dConc = subRegion->getReference<array1d<real64> >( viewKeyStruct::deltaProppantConcentrationString );

    // call the application of the boundary condition to alter the matrix and rhs
    fs->ApplyBoundaryConditionToSystem<FieldSpecificationEqual>( lset,
                                                                 false,
                                                                 time_n + dt,
                                                                 subRegion,
                                                                 dofNumber,
                                                                 2,
                                                                 blockSystem,
                                                                 BlockIDs::proppantTransportBlock,
                                                                 [&] (localIndex const a) -> real64
    {
      return conc[a] + dConc[a];
    });
  
  });
  

  if (verboseLevel() >= 2)
  {

    Epetra_FECrsMatrix * const jacobian = blockSystem->GetMatrix( BlockIDs::proppantTransportBlock, BlockIDs::proppantTransportBlock );
    Epetra_FEVector * const residual = blockSystem->GetResidualVector( BlockIDs::proppantTransportBlock );

    GEOS_LOG_RANK( "After ProppantTransport::ApplyBoundaryConditions" );
    GEOS_LOG_RANK( "\nJacobian\n" << *jacobian );
    GEOS_LOG_RANK( "\nResidual\n" << *residual );
  }

}

real64
ProppantTransport::
CalculateResidualNorm( DomainPartition const * const domain,
                       DofManager const & dofManager,
                       ParallelVector const & rhs )
{
  Epetra_FEVector const * const residual = blockSystem->GetResidualVector( BlockIDs::proppantTransportBlock );
  Epetra_Map      const * const rowMap   = blockSystem->GetRowMap( BlockIDs::proppantTransportBlock );

  MeshLevel const * const mesh = domain->getMeshBody(0)->getMeshLevel(0);

  // get a view into local residual vector
  int localSizeInt;
  double* localResidual = nullptr;
  residual->ExtractView(&localResidual, &localSizeInt);

  // compute the norm of local residual scaled by cell pore volume
  real64 localResidualNorm = 0.0;

  applyToSubRegions( mesh, [&] ( localIndex const er, localIndex const esr,
                                 ElementRegion const * const region,
                                 ElementSubRegionBase const * const subRegion )
  {
    arrayView1d<integer const> const & elemGhostRank = m_elemGhostRank[er][esr];
    arrayView1d<globalIndex const> const & dofNumber = m_dofNumber[er][esr];
    arrayView1d<real64 const> const & refPoro        = m_porosityRef[er][esr];
    arrayView1d<real64 const> const & volume         = m_volume[er][esr];

    localResidualNorm += sum_in_range( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex const a )
    {
      if (elemGhostRank[a] < 0)
      {
        real64 cell_norm = 0.0;
        globalIndex const offset = m_numDofPerCell * dofNumber[a];
        for (localIndex idof = 0; idof < m_numDofPerCell; ++idof)
        {
          int const lid = rowMap->LID(integer_conversion<int>(offset + idof));
          real64 const val = localResidual[lid] / volume[a];
          cell_norm += val * val;
        }
        return cell_norm;
      }
      return 0.0;
    } );
  } );

  // compute global residual norm
  real64 globalResidualNorm;
  MPI_Allreduce(&localResidualNorm, &globalResidualNorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_GEOSX);

  return sqrt(globalResidualNorm);
}

void ProppantTransport::ApplySystemSolution( DofManager const & dofManager,
                                             ParallelVector const & solution,
                                             real64 const scalingFactor,
                                             DomainPartition * const domain )
{
  
  Epetra_Map const * const rowMap        = blockSystem->GetRowMap( BlockIDs::proppantTransportBlock );
  Epetra_FEVector const * const solution = blockSystem->GetSolutionVector( BlockIDs::proppantTransportBlock );

  MeshLevel * mesh = domain->getMeshBody(0)->getMeshLevel(0);

  int dummy;
  double* local_solution = nullptr;
  solution->ExtractView( &local_solution, &dummy );

  applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegionBase * const region,
                                 ElementSubRegionBase * const subRegion )
  {
    arrayView1d<globalIndex const> const & dofNumber = m_dofNumber[er][esr];
    arrayView1d<integer const> const & elemGhostRank = m_elemGhostRank[er][esr];

    arrayView1d<real64> const & dPres = m_deltaPressure[er][esr];
    arrayView1d<real64> const & dConc = m_deltaProppantConcentration[er][esr];    
    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      if (elemGhostRank[ei] < 0)
      {
        int lid;
        lid = rowMap->LID( integer_conversion<int>( dofNumber[ei] * m_numDofPerCell) );
        dPres[ei] += scalingFactor * local_solution[lid];

        lid = rowMap->LID( integer_conversion<int>( dofNumber[ei] * m_numDofPerCell + 1) );
        dConc[ei] += scalingFactor * local_solution[lid];
      }
    } );
  } );

  std::map<string, string_array > fieldNames;
  fieldNames["elems"].push_back( viewKeyStruct::deltaPressureString );
  fieldNames["elems"].push_back( viewKeyStruct::deltaProppantConcentrationString );

  array1d<NeighborCommunicator> &
  comms = domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors );

  CommunicationTools::SynchronizeFields( fieldNames, mesh, comms );

  applyToSubRegions( mesh, [&] ( ElementSubRegionBase * subRegion )
  {
    UpdateState( subRegion );
  } );

}

void ProppantTransport::SolveSystem( DofManager const & dofManager,
                                     ParallelMatrix & matrix,
                                     ParallelVector & rhs,
                                     ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  Epetra_FEVector * const
  solution = blockSystem->GetSolutionVector( BlockIDs::proppantTransportBlock );

  Epetra_FEVector * const
  residual = blockSystem->GetResidualVector( BlockIDs::proppantTransportBlock );

  residual->Scale(-1.0);

  solution->Scale(0.0);
  
  m_linearSolverWrapper.SolveSingleBlockSystem( blockSystem,
                                                params,
                                                BlockIDs::proppantTransportBlock );

  if( verboseLevel() >= 2 )
  {
    GEOS_LOG_RANK("After ProppantTransport::SolveSystem");
    GEOS_LOG_RANK("\nsolution\n" << *solution);
  }

}

void ProppantTransport::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  MeshLevel * mesh = domain->getMeshBody(0)->getMeshLevel(0);

  applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegionBase * const region,
                                 ElementSubRegionBase * const subRegion )
  {
    arrayView1d<real64> const & dPres = m_deltaPressure[er][esr];
    arrayView1d<real64> const & dConc = m_deltaProppantConcentration[er][esr];    

    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      dPres[ei] = 0.0;
      dConc[ei] = 0.0;      
    } );

    UpdateState( subRegion );
  } );

}

void ProppantTransport::ResetViews(DomainPartition * const domain)
{
  FlowSolverBase::ResetViews(domain);

  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = mesh->getElemManager();
  ConstitutiveManager * const constitutiveManager = domain->getConstitutiveManager();

  m_dofNumber =
    elemManager->ConstructViewAccessor<array1d<globalIndex>, arrayView1d<globalIndex>>( viewKeyStruct::blockLocalDofNumberString );

  m_pressure =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( viewKeyStruct::pressureString );
  m_deltaPressure =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( viewKeyStruct::deltaPressureString );

  m_proppantConcentration =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( viewKeyStruct::proppantConcentrationString );
  m_deltaProppantConcentration =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( viewKeyStruct::deltaProppantConcentrationString );
  
  m_deltaVolume =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( viewKeyStruct::deltaVolumeString );

  /*
  m_slipVelocity =
    elemManager->ConstructViewAccessor<array1d<R1Tensor>, arrayView1d<R1Tensor>>( viewKeyStruct::slipVelocityString );
  */
  
  m_densityOld =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( viewKeyStruct::densityString );

  m_proppantConcentrationOld =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( viewKeyStruct::oldProppantConcentrationString );  

  m_density =
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64> >( SlurryFluidBase::viewKeyStruct::densityString, constitutiveManager );

  m_fluidDensity = 
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64> >( SlurryFluidBase::viewKeyStruct::fluidDensityString, constitutiveManager );

  m_dDens_dPres =
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64> >( SlurryFluidBase::viewKeyStruct::dDens_dPresString, constitutiveManager );

  m_dDens_dConc =
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64> >( SlurryFluidBase::viewKeyStruct::dDens_dConcString, constitutiveManager );

  m_dFluidDens_dPres =
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64> >( SlurryFluidBase::viewKeyStruct::dFluidDens_dPresString, constitutiveManager );

  m_viscosity =
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64> >( SlurryFluidBase::viewKeyStruct::viscosityString, constitutiveManager );
  
  m_dVisc_dPres =
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64> >( SlurryFluidBase::viewKeyStruct::dVisc_dPresString, constitutiveManager );

  m_dVisc_dConc =
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64> >( SlurryFluidBase::viewKeyStruct::dVisc_dConcString, constitutiveManager );

  m_pvMult =
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64> >( ConstitutiveBase::viewKeyStruct::poreVolumeMultiplierString,
                                                                                           constitutiveManager );
  m_dPvMult_dPres =
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64> >( ConstitutiveBase::viewKeyStruct::dPVMult_dPresString,
                                                                                           constitutiveManager );
  m_porosity =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( viewKeyStruct::porosityString );


  m_settlingFactor = 
    elemManager->ConstructFullMaterialViewAccessor<array1d<real64>, arrayView1d<real64> >( ParticleFluidBase::viewKeyStruct::settlingFactorString, constitutiveManager );

  m_dSettlingFactor_dConc = 
    elemManager->ConstructFullMaterialViewAccessor<array1d<real64>, arrayView1d<real64> >( ParticleFluidBase::viewKeyStruct::dSettlingFactor_dConcString, constitutiveManager );

  m_collisionFactor = 
    elemManager->ConstructFullMaterialViewAccessor<array1d<real64>, arrayView1d<real64> >( ParticleFluidBase::viewKeyStruct::collisionFactorString, constitutiveManager );

  m_dCollisionFactor_dConc = 
    elemManager->ConstructFullMaterialViewAccessor<array1d<real64>, arrayView1d<real64> >( ParticleFluidBase::viewKeyStruct::dCollisionFactor_dConcString, constitutiveManager );

  m_isProppantMobile = 
    elemManager->ConstructFullMaterialViewAccessor<array1d<bool>, arrayView1d<bool> >( ParticleFluidBase::viewKeyStruct::isProppantMobileString, constitutiveManager );
  
  m_proppantPackPermeability = 
    elemManager->ConstructFullMaterialViewAccessor<array1d<real64>, arrayView1d<real64> >( ParticleFluidBase::viewKeyStruct::proppantPackPermeabilityString, constitutiveManager );

}

REGISTER_CATALOG_ENTRY( SolverBase, ProppantTransport, std::string const &, Group * const )
} /* namespace geosx */
