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

    elemManager->forElementSubRegions<CellElementSubRegion>( [&]( CellElementSubRegion * const subRegion )
    {

      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::proppantConcentrationString )->setDefaultValue(0.0)->setPlotLevel(PlotLevel::LEVEL_0);
      
      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::deltaProppantConcentrationString )->setDefaultValue(0.0);      
      
      subRegion->registerWrapper< array2d<real64> >( viewKeyStruct::componentConcentrationString )->setDefaultValue(0.0)->
        setPlotLevel(PlotLevel::LEVEL_0);

      subRegion->registerWrapper< array2d<real64> >( viewKeyStruct::deltaComponentConcentrationString )->setDefaultValue(0.0);      

      subRegion->registerWrapper< array2d<real64> >( viewKeyStruct::updatedComponentConcentrationString )->setDefaultValue(0.0);

      subRegion->registerWrapper< array1d<R1Tensor> >( viewKeyStruct::shearRateString );

      /*
      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::oldProppantConcentrationString );

      subRegion->registerWrapper< array2d<real64> >( viewKeyStruct::oldComponentDensityString );
      */
      
    });

    
    elemManager->forElementSubRegions<FaceElementSubRegion>( [&]( FaceElementSubRegion * const subRegion )
    {

      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::proppantConcentrationString )->
        setPlotLevel(PlotLevel::LEVEL_0);

      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::deltaProppantConcentrationString );

      subRegion->registerWrapper< array2d<real64> >( viewKeyStruct::componentConcentrationString )->
        setPlotLevel(PlotLevel::LEVEL_0);

      subRegion->registerWrapper< array2d<real64> >( viewKeyStruct::deltaComponentConcentrationString );            

      subRegion->registerWrapper< array2d<real64> >( viewKeyStruct::updatedComponentConcentrationString );
      
      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::oldProppantConcentrationString );

      subRegion->registerWrapper< array2d<real64> >( viewKeyStruct::oldComponentDensityString );
      
      subRegion->registerWrapper< array1d<R1Tensor> >( viewKeyStruct::shearRateString );
      
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

  SlurryFluidBase const * fluid = cm->GetConstitutiveRelation<SlurryFluidBase>( m_fluidName );

  m_numComponents = fluid->numFluidComponents();
  m_numDofPerCell = m_numComponents + 1;

  MeshLevel * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = meshLevel->getElemManager();

  localIndex const NC = m_numComponents;
  
  elemManager->forElementSubRegions<CellElementSubRegion>([&]( CellElementSubRegion * const subRegion )
  {

    subRegion->template getReference< array2d<real64> >(viewKeyStruct::componentConcentrationString).resizeDimension<1>(NC);
    subRegion->template getReference< array2d<real64> >(viewKeyStruct::deltaComponentConcentrationString).resizeDimension<1>(NC);
    /*
    subRegion->template getReference< array2d<real64> >(viewKeyStruct::updatedComponentConcentrationString).resizeDimension<1>(NC);
    subRegion->template getReference< array2d<real64> >(viewKeyStruct::oldComponentDensityString).resizeDimension<1>(NC);        
    */
  });

}

void ProppantTransport::ResizeFractureFields( real64 const & GEOSX_UNUSED_ARG( time_n ),
                                              real64 const & GEOSX_UNUSED_ARG( dt ),
                                              DomainPartition * const domain)
{

  localIndex const NC = m_numComponents;  

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  
  ElementRegionManager * const elemManager = mesh->getElemManager();
  
  elemManager->forElementSubRegions<FaceElementSubRegion>([&]( FaceElementSubRegion * const subRegion )
  {

    subRegion->template getReference< array2d<real64> >(viewKeyStruct::componentConcentrationString).resizeDimension<1>(NC);
    subRegion->template getReference< array2d<real64> >(viewKeyStruct::deltaComponentConcentrationString).resizeDimension<1>(NC);
    subRegion->template getReference< array2d<real64> >(viewKeyStruct::updatedComponentConcentrationString).resizeDimension<1>(NC);
    subRegion->template getReference< array2d<real64> >(viewKeyStruct::oldComponentDensityString).resizeDimension<1>(NC);    
    
  });

} 

void ProppantTransport::UpdateFluidModel(Group * const dataGroup)
{
  GEOSX_MARK_FUNCTION;

  SlurryFluidBase * const fluid = GetConstitutiveModel<SlurryFluidBase>( dataGroup, m_fluidName );

  arrayView1d<real64 const> const & pres = dataGroup->getReference<array1d<real64>>( viewKeyStruct::pressureString );
  arrayView1d<real64 const> const & dPres = dataGroup->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

  arrayView2d<real64 const> const & componentConc = dataGroup->getReference<array2d<real64>>( viewKeyStruct::componentConcentrationString );
  arrayView2d<real64 const> const & dComponentConc = dataGroup->getReference<array2d<real64>>( viewKeyStruct::deltaComponentConcentrationString );

  arrayView2d<real64> const & updatedComponentConc = dataGroup->getReference<array2d<real64>>( viewKeyStruct::updatedComponentConcentrationString );    

  forall_in_range<RAJA::seq_exec>( 0, dataGroup->size(), GEOSX_LAMBDA ( localIndex const a )
  {

    for(localIndex c = 0; c < m_numComponents; ++c)
      {
        updatedComponentConc[a][c] = componentConc[a][c] + dComponentConc[a][c];
      }
    
    fluid->PointUpdateFluidProperty( pres[a] + dPres[a], updatedComponentConc[a], 0.0, a, 0 );

  });

}

void ProppantTransport::UpdateProppantModel(Group * const dataGroup)
{
  
  GEOSX_MARK_FUNCTION;

  
  SlurryFluidBase * const fluid = GetConstitutiveModel<SlurryFluidBase>( dataGroup, m_fluidName );

  /*
  array1d<real64> const & nIndices = fluid->nIndex();
  array1d<real64> const & KIndices = fluid->KIndex();  
  */
  
  localIndex const NC = m_numComponents;

  ParticleFluidBase * const particle = GetConstitutiveModel<ParticleFluidBase>( dataGroup, m_proppantName );

  arrayView1d<real64 const> const & proppantConc = dataGroup->getReference<array1d<real64>>( viewKeyStruct::proppantConcentrationString );

  arrayView1d<real64 const> const & dProppantConc = dataGroup->getReference<array1d<real64>>( viewKeyStruct::deltaProppantConcentrationString );
  
  arrayView2d<real64 const> const & fluidDens = fluid->getReference<array2d<real64>>( SlurryFluidBase::viewKeyStruct::fluidDensityString );

  arrayView2d<real64 const> const & dFluidDens_dPres = fluid->getReference<array2d<real64>>( SlurryFluidBase::viewKeyStruct::dFluidDens_dPresString );

  arrayView3d<real64 const> const & dFluidDens_dCompConc = fluid->getReference<array3d<real64>>( SlurryFluidBase::viewKeyStruct::dFluidDens_dCompConcString );        

  arrayView2d<real64 const> const & fluidVisc = fluid->getReference<array2d<real64>>( SlurryFluidBase::viewKeyStruct::fluidViscosityString );

  arrayView2d<real64 const> const & dFluidVisc_dPres = fluid->getReference<array2d<real64>>( SlurryFluidBase::viewKeyStruct::dFluidVisc_dPresString );

  arrayView3d<real64 const> const & dFluidVisc_dCompConc = fluid->getReference<array3d<real64>>( SlurryFluidBase::viewKeyStruct::dFluidVisc_dCompConcString );        

  forall_in_range<RAJA::seq_exec>( 0, dataGroup->size(), GEOSX_LAMBDA ( localIndex const a )
  {     
    particle->PointUpdate(NC, proppantConc[a] + dProppantConc[a], fluidDens[a][0], dFluidDens_dPres[a][0], dFluidDens_dCompConc[a][0], fluidVisc[a][0], dFluidVisc_dPres[a][0], dFluidVisc_dCompConc[a][0], a);

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

  fieldNames["elems"].push_back( viewKeyStruct::proppantConcentrationString );
  fieldNames["elems"].push_back( viewKeyStruct::componentConcentrationString );
  
  array1d<NeighborCommunicator> & comms =
    domain->getReference< array1d<NeighborCommunicator>>( domain->viewKeys.neighbors );

  CommunicationTools::SynchronizeFields( fieldNames, mesh, comms );

  ResetViews( domain );  

  // We have to redo the below loop after fractures are generated

  localIndex const NC = m_numComponents;
  
  applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegionBase * const GEOSX_UNUSED_ARG( region ),
                                 ElementSubRegionBase * const subRegion )
  {
    UpdateState( subRegion );

    arrayView1d<real64> const & proppantConc = m_proppantConcentration[er][esr];
    arrayView1d<real64> const & proppantConcOld = m_proppantConcentrationOld[er][esr];

    arrayView3d<real64> const & componentDens = m_componentDensity[er][esr][m_fluidIndex];
    arrayView2d<real64> const & componentDensOld = m_componentDensityOld[er][esr];    


    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
        proppantConcOld[ei] = proppantConc[ei];

        for(localIndex c = 0; c < NC; ++c)
          componentDensOld[ei][c] = componentDens[ei][0][c];        
    });

  } );

  m_downVector = getGravityVector();
  m_downVector.Normalize();
  
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

    localIndex const NC = m_numComponents;
    
    applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegionBase * const GEOSX_UNUSED_ARG( region ),
                                 ElementSubRegionBase * const subRegion )
    {


      subRegion->CalculateElementGeometricQuantities( *nodeManager,
                                                      *faceManager );
      

      
      UpdateState( subRegion );

      arrayView1d<real64> const & dProppantConc   = m_deltaProppantConcentration[er][esr];
      arrayView2d<real64> const & dComponentConc   = m_deltaComponentConcentration[er][esr];              


      arrayView1d<real64> const & proppantConcOld = m_proppantConcentrationOld[er][esr];
      arrayView1d<real64> const & proppantConc = m_proppantConcentration[er][esr];      

      arrayView2d<real64> const & componentDensOld = m_componentDensityOld[er][esr];
      arrayView3d<real64> const & componentDens = m_componentDensity[er][esr][m_fluidIndex];      

      arrayView1d<R1Tensor> const & shearRate   = m_shearRate[er][esr];
      

      forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
      {
        dProppantConc[ei] = 0.0;

        proppantConcOld[ei] = proppantConc[ei];

        for(localIndex c = 0; c < NC; ++c)
          {

            dComponentConc[ei][c] = 0.0;        
            componentDensOld[ei][c] = componentDens[ei][0][c];

          }

        shearRate[ei] = 0.0;

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


void ProppantTransport::PreStepUpdate( real64 const&,
                                       real64 const&,
                                       const int,
                                       DomainPartition * domain )
{
  GEOSX_MARK_FUNCTION;

  FlowSolverBase::PrecomputeData(domain);

  MeshLevel * mesh = domain->getMeshBody(0)->getMeshLevel(0);

  if(m_updateProppantMobility)
    {
      applyToSubRegions( mesh, [&] ( ElementSubRegionBase * const subRegion )
     {
       UpdateProppantModelStep( subRegion );
     });

    }
      
}


void ProppantTransport::ImplicitStepSetup( real64 const & GEOSX_UNUSED_ARG( time_n ),
					   real64 const & GEOSX_UNUSED_ARG( dt ),
                                           DomainPartition * const domain,
                                           DofManager & GEOSX_UNUSED_ARG(dofManager),
                                           ParallelMatrix & GEOSX_UNUSED_ARG(matrix),
                                           ParallelVector & GEOSX_UNUSED_ARG(rhs),
					   ParallelVector & GEOSX_UNUSED_ARG(solution) )
{

  localIndex const NC = m_numComponents;  

  ResetViews( domain );

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  /* The loop below could be moved to SolverStep after ImplicitStepSetup */
  
  applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegionBase * const GEOSX_UNUSED_ARG( region ),
                                 ElementSubRegionBase * const subRegion )
  {

    arrayView1d<real64> const & dProppantConc   = m_deltaProppantConcentration[er][esr];
    arrayView2d<real64> const & dComponentConc   = m_deltaComponentConcentration[er][esr];        
    
    arrayView1d<real64> const & proppantConc = m_proppantConcentration[er][esr];
    arrayView1d<real64> const & proppantConcOld = m_proppantConcentrationOld[er][esr];

    arrayView3d<real64> const & componentDens = m_componentDensity[er][esr][m_fluidIndex];
    arrayView2d<real64> const & componentDensOld = m_componentDensityOld[er][esr];    

    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {

      dProppantConc[ei] = 0.0;
      proppantConcOld[ei] = proppantConc[ei];

      for(localIndex c = 0; c < NC; ++c)
        {

          dComponentConc[ei][c] = 0.0;                      
          componentDensOld[ei][c] = componentDens[ei][0][c];            

        }
          
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

  localIndex const NC = m_numComponents;
  
  applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegionBase * const GEOSX_UNUSED_ARG( region ),
                                 ElementSubRegionBase * const subRegion )
  {
    arrayView1d<real64> const & proppantConc = m_proppantConcentration[er][esr];
    arrayView1d<real64 const> const & dProppantConc = m_deltaProppantConcentration[er][esr];    

    arrayView2d<real64> const & componentConc = m_componentConcentration[er][esr];
    arrayView2d<real64 const> const & dComponentConc = m_deltaComponentConcentration[er][esr];    
    
    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      proppantConc[ei] += dProppantConc[ei];      

      for(localIndex c = 0; c < NC; ++c)
        componentConc[ei][c] += dComponentConc[ei][c];
      
    } );

  } );

}

void ProppantTransport::SetupDofs( DomainPartition const * const GEOSX_UNUSED_ARG(domain),
                                   DofManager & dofManager ) const
{
  dofManager.addField( viewKeyStruct::proppantConcentrationString,
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

  localIndex const NC = m_numComponents;
  
  MeshLevel const * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegionBase const * const GEOSX_UNUSED_ARG( region ),
                                 ElementSubRegionBase const * const subRegion )
  {
    string const dofKey = dofManager->getKey( viewKeyStruct::proppantConcentrationString );
    arrayView1d<globalIndex const> const & dofNumber = subRegion->getReference< array1d<globalIndex> >( dofKey );

    arrayView1d<integer const>     const & elemGhostRank = m_elemGhostRank[er][esr];

    arrayView1d<real64 const> const & proppantConcOld       = m_proppantConcentrationOld[er][esr];
    arrayView2d<real64 const> const & componentDensOld       = m_componentDensityOld[er][esr];        

    arrayView1d<real64 const> const & volume        = m_volume[er][esr];

    arrayView1d<real64 const> const & proppantConc          = m_proppantConcentration[er][esr];

    arrayView1d<real64 const> const & dProppantConc          = m_deltaProppantConcentration[er][esr];

    arrayView3d<real64 const> const & componentDens          = m_componentDensity[er][esr][m_fluidIndex];

    arrayView3d<real64 const> const & dCompDens_dPres          = m_dComponentDensity_dPressure[er][esr][m_fluidIndex];

    arrayView4d<real64 const> const & dCompDens_dCompConc          = m_dComponentDensity_dComponentConcentration[er][esr][m_fluidIndex];            

    arrayView1d<R1Tensor> const & shearRate        = m_shearRate[er][esr];    
    
    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      if (elemGhostRank[ei] < 0)
      {

        stackArray1d<globalIndex, MAX_NUM_COMPONENTS>         localAccumDOF( m_numDofPerCell );
        stackArray1d<real64, MAX_NUM_COMPONENTS>             localAccum( m_numDofPerCell );
        stackArray2d<real64, MAX_NUM_COMPONENTS * MAX_NUM_COMPONENTS> localAccumJacobian( m_numDofPerCell, m_numDofPerCell );

	AccumulationKernel::Compute(NC,
                                    proppantConcOld[ei],
				    proppantConc[ei] + dProppantConc[ei],
                                    componentDensOld[ei],
                                    componentDens[ei][0],
                                    dCompDens_dPres[ei][0],
                                    dCompDens_dCompConc[ei][0],
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

      shearRate[ei] = 0.0;

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

  string const dofKey = dofManager->getKey( viewKeyStruct::proppantConcentrationString );

  ElementRegionManager::ElementViewAccessor< arrayView1d<globalIndex> >
  dofNumberAccessor = elemManager->ConstructViewAccessor< array1d<globalIndex>, arrayView1d<globalIndex> >( dofKey );

  FluxKernel::ElementViewConst< arrayView1d<globalIndex const> > const & dofNumber = dofNumberAccessor.toViewConst();

  FluxKernel::ElementViewConst < arrayView1d<real64 const> > const & pres        = m_pressure.toViewConst();
  
  FluxKernel::ElementViewConst < arrayView1d<real64 const> > const & dPres       = m_deltaPressure.toViewConst();

  FluxKernel::ElementViewConst < arrayView1d<real64 const> > const & proppantConc       = m_proppantConcentration.toViewConst();

  FluxKernel::ElementViewConst < arrayView1d<real64 const> > const & proppantConcOld       = m_proppantConcentrationOld.toViewConst();
  
  FluxKernel::ElementViewConst < arrayView1d<real64 const> > const & dProppantConc       = m_deltaProppantConcentration.toViewConst();

  FluxKernel::ElementViewConst < arrayView1d<real64 const> > const & gravDepth   = m_gravDepth.toViewConst();

  FluxKernel::MaterialView< arrayView2d<real64 const> > const & dens        = m_density.toViewConst();

  FluxKernel::MaterialView< arrayView2d<real64 const> > const & dDens_dPres = m_dDensity_dPressure.toViewConst();

  FluxKernel::MaterialView< arrayView2d<real64 const> > const & dDens_dProppantConc = m_dDensity_dProppantConcentration.toViewConst();

  FluxKernel::MaterialView< arrayView3d<real64 const> > const & dDens_dComponentConc = m_dDensity_dComponentConcentration.toViewConst();    

  FluxKernel::MaterialView< arrayView2d<real64 const> > const & visc        = m_viscosity.toViewConst();
  FluxKernel::MaterialView< arrayView2d<real64 const> > const & dVisc_dPres = m_dViscosity_dPressure.toViewConst();
  FluxKernel::MaterialView< arrayView2d<real64 const> > const & dVisc_dProppantConc = m_dViscosity_dProppantConcentration.toViewConst();
  FluxKernel::MaterialView< arrayView3d<real64 const> > const & dVisc_dComponentConc = m_dViscosity_dComponentConcentration.toViewConst();    

  FluxKernel::MaterialView< arrayView3d<real64 const> > const & componentDens = m_componentDensity.toViewConst();
  FluxKernel::MaterialView< arrayView3d<real64 const> > const & dComponentDens_dPres = m_dComponentDensity_dPressure.toViewConst();
  FluxKernel::MaterialView< arrayView4d<real64 const> > const & dComponentDens_dComponentConc = m_dComponentDensity_dComponentConcentration.toViewConst();    

  FluxKernel::MaterialView< arrayView2d<real64 const> > const & fluidDensity = m_fluidDensity.toViewConst();

  FluxKernel::MaterialView< arrayView2d<real64 const> > const & dFluidDens_dPres = m_dFluidDensity_dPressure.toViewConst();
  FluxKernel::MaterialView< arrayView3d<real64 const> > const & dFluidDens_dComponentConc = m_dFluidDensity_dComponentConcentration.toViewConst();    

  
  FluxKernel::MaterialView< arrayView1d<real64 const> > const & settlingFactor = m_settlingFactor.toViewConst();

  FluxKernel::MaterialView< arrayView1d<real64 const> > const & dSettlingFactor_dPres = m_dSettlingFactor_dPressure.toViewConst();
  FluxKernel::MaterialView< arrayView1d<real64 const> > const & dSettlingFactor_dProppantConc = m_dSettlingFactor_dProppantConcentration.toViewConst();
  FluxKernel::MaterialView< arrayView2d<real64 const> > const & dSettlingFactor_dComponentConc = m_dSettlingFactor_dComponentConcentration.toViewConst();  

  FluxKernel::MaterialView< arrayView1d<real64 const> > const & collisionFactor = m_collisionFactor.toViewConst();

  FluxKernel::MaterialView< arrayView1d<real64 const> > const & dCollisionFactor_dProppantConc = m_dCollisionFactor_dProppantConcentration.toViewConst();

  FluxKernel::MaterialView< arrayView1d<integer const> > const & isProppantMobile = m_isProppantMobile.toViewConst();
  
  FluxKernel::MaterialView< arrayView1d<real64 const> > const & proppantPackPermeability = m_proppantPackPermeability.toViewConst();      

  FluxKernel::ElementViewConst < arrayView1d<real64 const> > const & volume  = m_volume.toViewConst();

  FluxKernel::ElementViewConst < arrayView1d<real64 const> > const & aperture  = m_elementAperture.toViewConst();

  FluxKernel::ElementView  < arrayView1d<R1Tensor> >  & shearRate  = m_shearRate.toView();  
  
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
                        m_downVector,
                        dofNumber,
                        pres,
                        dPres,
			proppantConc,
			proppantConcOld,
			dProppantConc,
                        componentDens,
                        dComponentDens_dPres,
                        dComponentDens_dComponentConc,
                        gravDepth,
                        dens,
                        dDens_dPres,
                        dDens_dProppantConc,
                        dDens_dComponentConc,
                        visc,
                        dVisc_dPres,
                        dVisc_dProppantConc,
                        dVisc_dComponentConc,                        
			fluidDensity,
			dFluidDens_dPres,
			dFluidDens_dComponentConc,
			settlingFactor,
			dSettlingFactor_dPres,                        
			dSettlingFactor_dProppantConc,
			dSettlingFactor_dComponentConc,
			collisionFactor,
			dCollisionFactor_dProppantConc,
			isProppantMobile,
			proppantPackPermeability,
                        volume,
                        aperture,
                        shearRate,
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
  string const dofKey = dofManager.getKey( viewKeyStruct::proppantConcentrationString );

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
    proppantConc = subRegion->getReference<array1d<real64> >( viewKeyStruct::proppantConcentrationString );

    arrayView1d<real64 const> const &
    dProppantConc = subRegion->getReference<array1d<real64> >( viewKeyStruct::deltaProppantConcentrationString );

    if(fs->GetFluxFlag())
      {
        
        fs->ApplyBoundaryConditionToSystem<FieldSpecificationAdd, LAInterface>( lset,
                                                                                true,
                                                                                time_n + dt,
                                                                                dt,
                                                                                subRegion,
                                                                                dofNumber,
                                                                                m_numDofPerCell,
                                                                                matrix,
                                                                                rhs,
                                                                                [&]( localIndex const GEOSX_UNUSED_ARG(a), localIndex const GEOSX_UNUSED_ARG(c) ) -> real64
        {
          return 0;
        } );
      }
    else
      {
        
        fs->ApplyBoundaryConditionToSystem<FieldSpecificationEqual, LAInterface>( lset,
                                                                                  false,
                                                                                  time_n + dt,
                                                                                  subRegion,
                                                                                  dofNumber,
                                                                                  m_numDofPerCell,
                                                                                  matrix,
                                                                                  rhs,
                                                                                  [&]( localIndex const a, localIndex const GEOSX_UNUSED_ARG(c)) -> real64
        {
          return proppantConc[a] + dProppantConc[a];
        });
      }

  } );
      

 fsManager->Apply( time_n + dt, domain, "ElementRegions", viewKeyStruct::componentConcentrationString,
                    [&]( FieldSpecificationBase const * const fs,
                    string const &,
                    set<localIndex> const & lset,
                    Group * subRegion,
                    string const & ) -> void
  {


    arrayView1d<globalIndex const> const &
    dofNumber = subRegion->getReference< array1d<globalIndex> >( dofKey );

    arrayView2d<real64 const> const &
    componentConc = subRegion->getReference<array2d<real64> >( viewKeyStruct::componentConcentrationString );

    arrayView2d<real64 const> const &
    dComponentConc = subRegion->getReference<array2d<real64> >( viewKeyStruct::deltaComponentConcentrationString );

    if(fs->GetFluxFlag())    
      {
        
        fs->ApplyBoundaryConditionToSystem<FieldSpecificationAdd, LAInterface>( lset,
                                                                                true,
                                                                                time_n + dt,
                                                                                dt,
                                                                                subRegion,
                                                                                dofNumber,
                                                                                m_numDofPerCell,
                                                                                matrix,
                                                                                rhs,
                                                                                [&]( localIndex const GEOSX_UNUSED_ARG(a), localIndex const GEOSX_UNUSED_ARG(c) ) -> real64
        {
          return 0;
        } );
      }
    else 
      {
        
        fs->ApplyBoundaryConditionToSystem<FieldSpecificationEqual, LAInterface>( lset,
                                                                                  false,
                                                                                  time_n + dt,
                                                                                  subRegion,
                                                                                  dofNumber,
                                                                                  m_numDofPerCell,
                                                                                  matrix,
                                                                                  rhs,
                                                                                  [&]( localIndex const a, localIndex const c ) -> real64
        {
          return componentConc[a][c-1] + dComponentConc[a][c-1];
        });
      }
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

  string const dofKey = dofManager.getKey( viewKeyStruct::proppantConcentrationString );
  
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
				 viewKeyStruct::proppantConcentrationString,
                                 scalingFactor,
                                 subRegion,
                                 viewKeyStruct::deltaProppantConcentrationString,
                                 0, 1 );


    if(m_numDofPerCell > 1) 
    dofManager.addVectorToField( solution,
				 viewKeyStruct::proppantConcentrationString,
                                 scalingFactor,
                                 subRegion,
                                 viewKeyStruct::deltaComponentConcentrationString,
                                 1, m_numDofPerCell );
  } );

  std::map<string, string_array > fieldNames;
  fieldNames["elems"].push_back( viewKeyStruct::deltaProppantConcentrationString );
  fieldNames["elems"].push_back( viewKeyStruct::deltaComponentConcentrationString );  

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

  localIndex const NC = m_numComponents;
  
  applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegionBase * const GEOSX_UNUSED_ARG( region ),
                                 ElementSubRegionBase * const subRegion )
  {
    arrayView1d<real64> const & dProppantConc = m_deltaProppantConcentration[er][esr];
    arrayView2d<real64> const & dComponentConc = m_deltaComponentConcentration[er][esr];        

    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      dProppantConc[ei] = 0.0;

      for(localIndex c = 0; c < NC; ++c)
        dComponentConc[ei][c] = 0.0;
      
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

  m_componentConcentration =
    elemManager->ConstructViewAccessor<array2d<real64>, arrayView2d<real64>>( viewKeyStruct::componentConcentrationString );
  m_deltaComponentConcentration =
    elemManager->ConstructViewAccessor<array2d<real64>, arrayView2d<real64>>( viewKeyStruct::deltaComponentConcentrationString );
  
  m_proppantConcentrationOld =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( viewKeyStruct::oldProppantConcentrationString );

  m_componentDensityOld =
    elemManager->ConstructViewAccessor<array2d<real64>, arrayView2d<real64>>( viewKeyStruct::oldComponentDensityString );

  m_updatedComponentConcentration =
    elemManager->ConstructViewAccessor<array2d<real64>, arrayView2d<real64>>( viewKeyStruct::updatedComponentConcentrationString );

  m_shearRate =
    elemManager->ConstructViewAccessor<array1d<R1Tensor>, arrayView1d<R1Tensor>>( viewKeyStruct::shearRateString );        

  m_density =
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64> >( SlurryFluidBase::viewKeyStruct::densityString, constitutiveManager );

  m_dDensity_dPressure =
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64> >( SlurryFluidBase::viewKeyStruct::dDens_dPresString, constitutiveManager );

  m_dDensity_dProppantConcentration =
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64> >( SlurryFluidBase::viewKeyStruct::dDens_dProppantConcString, constitutiveManager );

  m_dDensity_dComponentConcentration =
    elemManager->ConstructFullMaterialViewAccessor<array3d<real64>, arrayView3d<real64> >( SlurryFluidBase::viewKeyStruct::dDens_dCompConcString, constitutiveManager );

  m_componentDensity = 
    elemManager->ConstructFullMaterialViewAccessor<array3d<real64>, arrayView3d<real64> >( SlurryFluidBase::viewKeyStruct::componentDensityString, constitutiveManager );

  m_dComponentDensity_dPressure =
    elemManager->ConstructFullMaterialViewAccessor<array3d<real64>, arrayView3d<real64> >( SlurryFluidBase::viewKeyStruct::dCompDens_dPresString, constitutiveManager );

  m_dComponentDensity_dComponentConcentration =
    elemManager->ConstructFullMaterialViewAccessor<array4d<real64>, arrayView4d<real64> >( SlurryFluidBase::viewKeyStruct::dCompDens_dCompConcString, constitutiveManager );  

  m_fluidDensity = 
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64> >( SlurryFluidBase::viewKeyStruct::fluidDensityString, constitutiveManager );

  m_dFluidDensity_dPressure =
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64> >( SlurryFluidBase::viewKeyStruct::dFluidDens_dPresString, constitutiveManager );

  m_dFluidDensity_dComponentConcentration =
    elemManager->ConstructFullMaterialViewAccessor<array3d<real64>, arrayView3d<real64> >( SlurryFluidBase::viewKeyStruct::dFluidDens_dCompConcString, constitutiveManager );  

  m_viscosity =
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64> >( SlurryFluidBase::viewKeyStruct::viscosityString, constitutiveManager );
  
  m_dViscosity_dPressure =
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64> >( SlurryFluidBase::viewKeyStruct::dVisc_dPresString, constitutiveManager );

  m_dViscosity_dProppantConcentration =
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64> >( SlurryFluidBase::viewKeyStruct::dVisc_dProppantConcString, constitutiveManager );

  m_dViscosity_dComponentConcentration =
    elemManager->ConstructFullMaterialViewAccessor<array3d<real64>, arrayView3d<real64> >( SlurryFluidBase::viewKeyStruct::dVisc_dCompConcString, constitutiveManager );  

  m_settlingFactor = 
    elemManager->ConstructFullMaterialViewAccessor<array1d<real64>, arrayView1d<real64> >( ParticleFluidBase::viewKeyStruct::settlingFactorString, constitutiveManager );

  m_dSettlingFactor_dPressure = 
    elemManager->ConstructFullMaterialViewAccessor<array1d<real64>, arrayView1d<real64> >( ParticleFluidBase::viewKeyStruct::dSettlingFactor_dPressureString, constitutiveManager );

  m_dSettlingFactor_dProppantConcentration = 
    elemManager->ConstructFullMaterialViewAccessor<array1d<real64>, arrayView1d<real64> >( ParticleFluidBase::viewKeyStruct::dSettlingFactor_dProppantConcentrationString, constitutiveManager );

  m_dSettlingFactor_dComponentConcentration = 
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64> >( ParticleFluidBase::viewKeyStruct::dSettlingFactor_dComponentConcentrationString, constitutiveManager );  

  m_collisionFactor = 
    elemManager->ConstructFullMaterialViewAccessor<array1d<real64>, arrayView1d<real64> >( ParticleFluidBase::viewKeyStruct::collisionFactorString, constitutiveManager );

  m_dCollisionFactor_dProppantConcentration = 
    elemManager->ConstructFullMaterialViewAccessor<array1d<real64>, arrayView1d<real64> >( ParticleFluidBase::viewKeyStruct::dCollisionFactor_dProppantConcentrationString, constitutiveManager );

  m_isProppantMobile = 
    elemManager->ConstructFullMaterialViewAccessor<array1d<integer>, arrayView1d<integer> >( ParticleFluidBase::viewKeyStruct::isProppantMobileString, constitutiveManager );
  
  m_proppantPackPermeability = 
    elemManager->ConstructFullMaterialViewAccessor<array1d<real64>, arrayView1d<real64> >( ParticleFluidBase::viewKeyStruct::proppantPackPermeabilityString, constitutiveManager );

}


REGISTER_CATALOG_ENTRY( SolverBase, ProppantTransport, std::string const &, Group * const )
} /* namespace geosx */
