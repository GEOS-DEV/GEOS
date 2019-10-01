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
#include "mpiCommunications/CommunicationTools.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"
#include "mesh/FaceElementRegion.hpp"

#include "physicsSolvers/FiniteVolume/ProppantTransportKernels.hpp"

/**
 * @namespace the geosx namespace that encapsulates the majority of the code
 */
namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace ProppantTransportKernels;
  
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
      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::aperture0String )->
        setDefaultValue(1e-5);
      
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
                                 ElementRegionBase * const GEOSX_UNUSED_ARG( region ),
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

  if(cycleNumber == 0) {

    /*  assign intitial and boundary conditions */
    FieldSpecificationManager const * boundaryConditionManager = FieldSpecificationManager::get();

    boundaryConditionManager->ApplyInitialConditions( domain );

    /* Below must be called after ImplicitStepSetup */

    applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegionBase * const GEOSX_UNUSED_ARG( region ),
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


void ProppantTransport::ImplicitStepSetup( real64 const & GEOSX_UNUSED_ARG( time_n ),
					   real64 const & GEOSX_UNUSED_ARG( dt ),
                                           DomainPartition * const domain,
                                           DofManager & GEOSX_UNUSED_ARG(dofManager),
                                           ParallelMatrix & GEOSX_UNUSED_ARG(matrix),
                                           ParallelVector & GEOSX_UNUSED_ARG(rhs),
					   ParallelVector & GEOSX_UNUSED_ARG(solution) )
{
  ResetViews( domain );

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  /* The loop below could be moved to SolverStep after ImplicitStepSetup */
  
  applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegionBase * const GEOSX_UNUSED_ARG( region ),
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

void ProppantTransport::ImplicitStepComplete( real64 const & GEOSX_UNUSED_ARG(time_n),
					      real64 const & GEOSX_UNUSED_ARG(dt),
                                            DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegionBase * const GEOSX_UNUSED_ARG( region ),
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

void ProppantTransport::SetupDofs( DomainPartition const * const GEOSX_UNUSED_ARG(domain),
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

  applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegionBase const * const GEOSX_UNUSED_ARG( region ),
                                 ElementSubRegionBase const * const subRegion )
  {
    string const dofKey = dofManager->getKey( viewKeyStruct::pressureString );
    arrayView1d<globalIndex const> const & dofNumber = subRegion->getReference< array1d<globalIndex> >( dofKey );

    arrayView1d<integer const>     const & elemGhostRank = m_elemGhostRank[er][esr];

    arrayView1d<real64 const> const & densOld       = m_densityOld[er][esr];

    arrayView1d<real64 const> const & concOld       = m_proppantConcentrationOld[er][esr];    

    arrayView1d<real64 const> const & volume        = m_volume[er][esr];

    arrayView2d<real64 const> const & dens          = m_density[er][esr][m_fluidIndex];

    arrayView1d<real64 const> const & conc          = m_proppantConcentration[er][esr];

    arrayView1d<real64 const> const & dConc          = m_deltaProppantConcentration[er][esr];    
    
    arrayView2d<real64 const> const & dDens_dPres   = m_dDens_dPres[er][esr][m_fluidIndex];
    arrayView2d<real64 const> const & dDens_dConc   = m_dDens_dConc[er][esr][m_fluidIndex];    

    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      if (elemGhostRank[ei] < 0)
      {

        stackArray1d<globalIndex, MAX_NUM_COMPONENTS>         localAccumDOF( m_numDofPerCell );
        stackArray1d<real64, MAX_NUM_COMPONENTS>             localAccum( m_numDofPerCell );
        stackArray2d<real64, MAX_NUM_COMPONENTS * MAX_NUM_COMPONENTS> localAccumJacobian( m_numDofPerCell, m_numDofPerCell );

	AccumulationKernel::Compute( dens[ei][0],
                                     densOld[ei],
				     dDens_dPres[ei][0],
				     dDens_dConc[ei][0],
				     conc[ei] + dConc[ei],
				     concOld[ei],
				     volume[ei],
				     localAccum,
				     localAccumJacobian );

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


void ProppantTransport::AssembleFluxTerms( real64 const GEOSX_UNUSED_ARG(time_n),
                                           real64 const dt,
                                           DomainPartition const * const domain,
                                           DofManager const * const dofManager,
                                           ParallelMatrix * const matrix,
                                           ParallelVector * const rhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager const * const elemManager = mesh->getElemManager();


  NumericalMethodsManager const * numericalMethodManager =
    domain->getParent()->GetGroup<NumericalMethodsManager>( keys::numericalMethodsManager );

  FiniteVolumeManager const * fvManager =
    numericalMethodManager->GetGroup<FiniteVolumeManager>( keys::finiteVolumeManager );

  FluxApproximationBase const * fluxApprox = fvManager->getFluxApproximation( m_discretizationName );

  string const dofKey = dofManager->getKey( viewKeyStruct::pressureString );

  ElementRegionManager::ElementViewAccessor< arrayView1d<globalIndex> >
  dofNumberAccessor = elemManager->ConstructViewAccessor< array1d<globalIndex>, arrayView1d<globalIndex> >( dofKey );

  FluxKernel::ElementView< arrayView1d<globalIndex const> > const & dofNumber = dofNumberAccessor.toViewConst();

  FluxKernel::ElementView < arrayView1d<real64 const> > const & pres        = m_pressure.toViewConst();
  
  FluxKernel::ElementView < arrayView1d<real64 const> > const & dPres       = m_deltaPressure.toViewConst();

  FluxKernel::ElementView < arrayView1d<real64 const> > const & conc       = m_proppantConcentration.toViewConst();

  FluxKernel::ElementView < arrayView1d<real64 const> > const & concOld       = m_proppantConcentrationOld.toViewConst();
  
  FluxKernel::ElementView < arrayView1d<real64 const> > const & dConc       = m_deltaProppantConcentration.toViewConst();
  
  FluxKernel::ElementView < arrayView1d<real64 const> > const & gravDepth   = m_gravDepth.toViewConst();

  FluxKernel::MaterialView< arrayView2d<real64 const> > const & dens        = m_density.toViewConst();

  FluxKernel::MaterialView< arrayView2d<real64 const> > const & dDens_dPres = m_dDens_dPres.toViewConst();
  FluxKernel::MaterialView< arrayView2d<real64 const> > const & dDens_dConc = m_dDens_dConc.toViewConst();  

  FluxKernel::MaterialView< arrayView2d<real64 const> > const & visc        = m_viscosity.toViewConst();
  FluxKernel::MaterialView< arrayView2d<real64 const> > const & dVisc_dPres = m_dVisc_dPres.toViewConst();
  FluxKernel::MaterialView< arrayView2d<real64 const> > const & dVisc_dConc = m_dVisc_dConc.toViewConst();  

  FluxKernel::MaterialView< arrayView2d<real64 const> > const & fluidDensity = m_fluidDensity.toViewConst();

  FluxKernel::MaterialView< arrayView2d<real64 const> > const & dFluidDens_dPres = m_dFluidDens_dPres.toViewConst();  

  
  FluxKernel::MaterialView< arrayView1d<real64 const> > const & settlingFactor = m_settlingFactor.toViewConst();

  FluxKernel::MaterialView< arrayView1d<real64 const> > const & dSettlingFactor_dConc = m_dSettlingFactor_dConc.toViewConst();

  FluxKernel::MaterialView< arrayView1d<real64 const> > const & collisionFactor = m_collisionFactor.toViewConst();

  FluxKernel::MaterialView< arrayView1d<real64 const> > const & dCollisionFactor_dConc = m_dCollisionFactor_dConc.toViewConst();

  FluxKernel::MaterialView< arrayView1d<bool const> > const & isProppantMobile = m_isProppantMobile.toViewConst();
  
  FluxKernel::MaterialView< arrayView1d<real64 const> > const & proppantPackPermeability = m_proppantPackPermeability.toViewConst();      

  FluxKernel::ElementView < arrayView1d<real64 const> > const & aperture0  = m_elementAperture0.toViewConst();

  FluxKernel::ElementView < arrayView1d<real64 const> > const & aperture  = m_elementAperture.toViewConst();
  
  localIndex const fluidIndex = m_fluidIndex;
  localIndex const proppantIndex = m_proppantIndex;  
  bool updateProppantMobilityFlag = m_updateProppantMobility;
  bool updatePermeabilityFlag = m_updatePermeability;  

  fluxApprox->forCellStencils( [&]( auto const & stencil )
  {

    FluxKernel::Launch( stencil,
			m_numDofPerCell,
                        dt,
                        fluidIndex,
			proppantIndex,
			updateProppantMobilityFlag,
			updatePermeabilityFlag,
                        dofNumber,
                        pres,
                        dPres,
			conc,
			concOld,
			dConc,
                        gravDepth,
                        dens,
                        dDens_dPres,
                        dDens_dConc,			
                        visc,
                        dVisc_dPres,
                        dVisc_dConc,
			fluidDensity,
			dFluidDens_dPres,
			settlingFactor,
			dSettlingFactor_dConc,
			collisionFactor,
			dCollisionFactor_dConc,
			isProppantMobile,
			proppantPackPermeability,
                        aperture0,
                        aperture,
                        matrix,
                        rhs );
  });

}

void ProppantTransport::ApplyBoundaryConditions(real64 const time_n,
						real64 const dt,
						DomainPartition * const domain,
						DofManager const & dofManager,
						ParallelMatrix & matrix,
						ParallelVector & rhs )
{
  GEOSX_MARK_FUNCTION;


  FieldSpecificationManager * fsManager = FieldSpecificationManager::get();
  string const dofKey = dofManager.getKey( viewKeyStruct::pressureString );
  
  fsManager->Apply( time_n + dt, domain, "ElementRegions", "FLUX",
                    [&]( FieldSpecificationBase const * const fs,
                    string const &,
                    set<localIndex> const & lset,
                    Group * subRegion,
                    string const & ) -> void
  {
    arrayView1d<globalIndex const> const &
    dofNumber = subRegion->getReference< array1d<globalIndex> >( dofKey );

    fs->ApplyBoundaryConditionToSystem<FieldSpecificationAdd, LAInterface>( lset,
                                                                            true,
                                                                            time_n + dt,
                                                                            dt,
                                                                            subRegion,
                                                                            dofNumber,
                                                                            2,
                                                                            matrix,
                                                                            rhs,
                                                                            [&]( localIndex const GEOSX_UNUSED_ARG(a) ) -> real64
    {
      return 0;
    } );

  } );

  fsManager->Apply( time_n + dt, domain, "ElementRegions", viewKeyStruct::pressureString,
                    [&]( FieldSpecificationBase const * const fs,
                    string const &,
                    set<localIndex> const & lset,
                    Group * subRegion,
                    string const & ) -> void
  {
    arrayView1d<globalIndex const> const &
    dofNumber = subRegion->getReference< array1d<globalIndex> >( dofKey );

    arrayView1d<real64 const> const &
    pres = subRegion->getReference<array1d<real64> >( viewKeyStruct::pressureString );

    arrayView1d<real64 const> const &
    dPres = subRegion->getReference<array1d<real64> >( viewKeyStruct::deltaPressureString );

    fs->ApplyBoundaryConditionToSystem<FieldSpecificationEqual, LAInterface>( lset,
                                                                              false,
                                                                              time_n + dt,
                                                                              subRegion,
                                                                              dofNumber,
                                                                              2,
                                                                              matrix,
                                                                              rhs,
                                                                              [&]( localIndex const a ) -> real64
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
    dofNumber = subRegion->getReference< array1d<globalIndex> >( dofKey );

    arrayView1d<real64 const> const &
    conc = subRegion->getReference<array1d<real64> >( viewKeyStruct::proppantConcentrationString );

    arrayView1d<real64 const> const &
    dConc = subRegion->getReference<array1d<real64> >( viewKeyStruct::deltaProppantConcentrationString );

    fs->ApplyBoundaryConditionToSystem<FieldSpecificationEqual, LAInterface>( lset,
                                                                              false,
                                                                              time_n + dt,
                                                                              subRegion,
                                                                              dofNumber,
                                                                              2,
                                                                              matrix,
                                                                              rhs,
                                                                              [&]( localIndex const a ) -> real64
    {
      return conc[a] + dConc[a];
    });
  
  });
  

  if( verboseLevel() >= 3 )
  {
    SystemSolverParameters * const solverParams = getSystemSolverParameters();
    integer newtonIter = solverParams->numNewtonIterations();

    string filename_mat = "matrix_bc_" + std::to_string( time_n ) + "_" + std::to_string( newtonIter ) + ".mtx";
    matrix.write( filename_mat, true );

    string filename_rhs = "rhs_bc_" + std::to_string( time_n ) + "_" + std::to_string( newtonIter ) + ".mtx";
    rhs.write( filename_rhs, true );

    GEOS_LOG_RANK_0( "After ProppantTransport::ApplyBoundaryConditions" );
    GEOS_LOG_RANK_0( "Jacobian: written to " << filename_mat );
    GEOS_LOG_RANK_0( "Residual: written to " << filename_rhs );
  } 

}

real64
ProppantTransport::
CalculateResidualNorm( DomainPartition const * const domain,
                       DofManager const & dofManager,
                       ParallelVector const & rhs )
{

  MeshLevel const * const mesh = domain->getMeshBody(0)->getMeshLevel(0);

  // get a view into local residual vector
  real64 const * localResidual = rhs.extractLocalVector();

  string const dofKey = dofManager.getKey( viewKeyStruct::pressureString );
  
  // compute the norm of local residual scaled by cell pore volume
  real64 localResidualNorm = 0.0;

  applyToSubRegions( mesh, [&] ( localIndex const er, localIndex const esr,
                                 ElementRegionBase const * const GEOSX_UNUSED_ARG( region ),
                                 ElementSubRegionBase const * const subRegion )
  {

    arrayView1d<globalIndex const> const & dofNumber = subRegion->getReference< array1d<globalIndex> >( dofKey );
    
    arrayView1d<integer const> const & elemGhostRank = m_elemGhostRank[er][esr];

    arrayView1d<real64 const> const & volume         = m_volume[er][esr];

    localIndex const subRegionSize = subRegion->size();
    for ( localIndex a = 0; a < subRegionSize; ++a )
    {

      if (elemGhostRank[a] < 0)
      {
        localIndex const offset = dofNumber[a];
        for (localIndex idof = 0; idof < m_numDofPerCell; ++idof)
        {
          localIndex const lid = rhs.getLocalRowID( offset + idof );
          real64 const val = localResidual[lid] / volume[a];
          localResidualNorm += val * val;
        }
      }
    }
    
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
  
  MeshLevel * mesh = domain->getMeshBody(0)->getMeshLevel(0);

  applyToSubRegions( mesh, [&] ( localIndex GEOSX_UNUSED_ARG(er), localIndex GEOSX_UNUSED_ARG(esr),
                                 ElementRegionBase * const GEOSX_UNUSED_ARG( region ),
                                 ElementSubRegionBase * const subRegion )
  {

    dofManager.addVectorToField( solution,
				 viewKeyStruct::pressureString,
                                 scalingFactor,
                                 subRegion,
                                 viewKeyStruct::deltaPressureString,
                                 0, 1 );

    dofManager.addVectorToField( solution,
				 viewKeyStruct::pressureString,
                                 scalingFactor,
                                 subRegion,
                                 viewKeyStruct::deltaProppantConcentrationString,
                                 1, m_numDofPerCell );
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

  rhs.scale( -1.0 );
  solution.zero();

  SolverBase::SolveSystem( dofManager, matrix, rhs, solution );
  
  if( verboseLevel() >= 2 )
  {
    GEOS_LOG_RANK("After ProppantTransport::SolveSystem");
    GEOS_LOG_RANK("\nsolution:\n" << solution);
  }

}

void ProppantTransport::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  MeshLevel * mesh = domain->getMeshBody(0)->getMeshLevel(0);

  applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegionBase * const GEOSX_UNUSED_ARG( region ),
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
