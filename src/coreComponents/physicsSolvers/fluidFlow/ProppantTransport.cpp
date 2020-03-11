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
 * @file ProppantTransport.cpp
 */

#include "ProppantTransport.hpp"

#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/fluid/SlurryFluidBase.hpp"
#include "constitutive/fluid/ParticleFluidBase.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "mpiCommunications/CommunicationTools.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"
#include "mesh/FaceElementRegion.hpp"

#include "physicsSolvers/fluidFlow/ProppantTransportKernels.hpp"

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

  registerWrapper( viewKeyStruct::bridgingFactorString, &m_bridgingFactor, false )->setApplyDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Bridging factor used for bridging/screen-out calculation");

  registerWrapper( viewKeyStruct::maxProppantConcentrationString, &m_maxProppantConcentration, false )->setApplyDefaultValue(0.6)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Maximum proppant concentration");

    registerWrapper( viewKeyStruct::proppantDiameterString, &m_proppantDiameter, false )->setApplyDefaultValue(0.4e-3)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Proppant diameter");

    registerWrapper( viewKeyStruct::proppantDensityString, &m_proppantDensity, false )->setApplyDefaultValue(2500.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Proppant density");

    registerWrapper( viewKeyStruct::criticalShieldsNumberString, &m_criticalShieldsNumber, false )->setApplyDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Critical Shields number");

    registerWrapper( viewKeyStruct::frictionCoefficientString, &m_frictionCoefficient, false )->setApplyDefaultValue(0.03)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Friction coefficient");        

    registerWrapper( viewKeyStruct::updateProppantPackingString, &m_updateProppantPacking, false )->setApplyDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Flag that enables/disables proppant-packing update");

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

      subRegion->registerWrapper< array1d<R1Tensor> >( viewKeyStruct::cellBasedFluxString )->setDefaultValue({0.0, 0.0, 0.0});

      subRegion->registerWrapper< array1d<integer> >( viewKeyStruct::isProppantBoundaryString )->setDefaultValue(0);

      subRegion->registerWrapper< array2d<real64> >( viewKeyStruct::bcComponentConcentrationString )->setDefaultValue(0.0);      
      
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
      
      subRegion->registerWrapper< array1d<R1Tensor> >( viewKeyStruct::cellBasedFluxString );

      subRegion->registerWrapper< array1d<integer> >( viewKeyStruct::isInterfaceElementString );

      subRegion->registerWrapper< array1d<integer> >( viewKeyStruct::isProppantBoundaryString );            

      subRegion->registerWrapper< array1d<integer> >( viewKeyStruct::isProppantMobileString );

      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::proppantPackVolumeFractionString )->
        setPlotLevel(PlotLevel::LEVEL_0);

      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::proppantExcessPackVolumeString );

      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::proppantLiftFluxString );                  

      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::poroMultiplierString )->setDefaultValue(1.0);

      subRegion->registerWrapper< array1d<R1Tensor> >( viewKeyStruct::transTMultiplierString )->setDefaultValue({1.0, 1.0, 1.0});            

      subRegion->registerWrapper< array2d<real64> >( viewKeyStruct::bcComponentConcentrationString )->setDefaultValue(0.0);
      
      
    } );

  }
}

void ProppantTransport::InitializePreSubGroups(Group * const rootGroup)
{
  FlowSolverBase::InitializePreSubGroups(rootGroup);

  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);
  
  ConstitutiveManager * const cm = domain->getConstitutiveManager();

  ConstitutiveBase const * proppant  = cm->GetConstitutiveRelation<ConstitutiveBase>( m_proppantName );
  GEOSX_ERROR_IF( proppant == nullptr, "Proppant model " + m_proppantName + " not found" );
  m_proppantIndex = proppant->getIndexInParent();

  SlurryFluidBase const * fluid = cm->GetConstitutiveRelation<SlurryFluidBase>( m_fluidName );

  m_numComponents = fluid->numFluidComponents();
  m_numDofPerCell = m_numComponents + 1;

  MeshLevel * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = meshLevel->getElemManager();

  localIndex const NC = m_numComponents;

  if(NC > 0)
    {
      
      elemManager->forElementSubRegions<CellElementSubRegion>([&]( CellElementSubRegion * const subRegion )
       {

         subRegion->template getReference< array2d<real64> >(viewKeyStruct::componentConcentrationString).resizeDimension<1>(NC);
         subRegion->template getReference< array2d<real64> >(viewKeyStruct::deltaComponentConcentrationString).resizeDimension<1>(NC);

     });
    }


}

void ProppantTransport::ResizeFractureFields( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                              real64 const & GEOSX_UNUSED_PARAM( dt ),
                                              DomainPartition * const domain)
{

  localIndex const NC = m_numComponents;  

  if(NC > 0)
    {
  
      MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  
      ElementRegionManager * const elemManager = mesh->getElemManager();
  
      elemManager->forElementSubRegions<FaceElementSubRegion>([&]( FaceElementSubRegion * const subRegion )
       {

         subRegion->template getReference< array2d<real64> >(viewKeyStruct::componentConcentrationString).resizeDimension<1>(NC);
         subRegion->template getReference< array2d<real64> >(viewKeyStruct::deltaComponentConcentrationString).resizeDimension<1>(NC);
         subRegion->template getReference< array2d<real64> >(viewKeyStruct::updatedComponentConcentrationString).resizeDimension<1>(NC);
         subRegion->template getReference< array2d<real64> >(viewKeyStruct::oldComponentDensityString).resizeDimension<1>(NC);

         subRegion->template getReference< array2d<real64> >(viewKeyStruct::bcComponentConcentrationString).resizeDimension<1>(NC);         
    
       });

    }
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

  forall_in_range<RAJA::seq_exec>( 0, dataGroup->size(), [=] ( localIndex const a )
  {

    for(localIndex c = 0; c < m_numComponents; ++c)
      {
        updatedComponentConc[a][c] = componentConc[a][c] + dComponentConc[a][c];
      }
    
    fluid->PointUpdateFluidProperty( pres[a] + dPres[a], updatedComponentConc[a], 0.0, a, 0 );

  });

}

void ProppantTransport::UpdateComponentDensity(Group * const dataGroup)
{
  GEOSX_MARK_FUNCTION;

  SlurryFluidBase * const fluid = GetConstitutiveModel<SlurryFluidBase>( dataGroup, m_fluidName );

  arrayView1d<real64 const> const & pres = dataGroup->getReference<array1d<real64>>( viewKeyStruct::pressureString );
  arrayView1d<real64 const> const & dPres = dataGroup->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

  arrayView2d<real64 const> const & componentConc = dataGroup->getReference<array2d<real64>>( viewKeyStruct::componentConcentrationString );
  arrayView2d<real64 const> const & dComponentConc = dataGroup->getReference<array2d<real64>>( viewKeyStruct::deltaComponentConcentrationString );

  arrayView2d<real64> const & updatedComponentConc = dataGroup->getReference<array2d<real64>>( viewKeyStruct::updatedComponentConcentrationString );    

  forall_in_range<RAJA::seq_exec>( 0, dataGroup->size(), [=] ( localIndex const a )
  {

    for(localIndex c = 0; c < m_numComponents; ++c)
      {
        updatedComponentConc[a][c] = componentConc[a][c] + dComponentConc[a][c];
      }
    
    fluid->PointUpdateComponentDensity( pres[a] + dPres[a], updatedComponentConc[a], a, 0 );

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

  forall_in_range<RAJA::seq_exec>( 0, dataGroup->size(), [=] ( localIndex const a )
  {     
    particle->PointUpdate(NC, proppantConc[a] + dProppantConc[a], fluidDens[a][0], dFluidDens_dPres[a][0], dFluidDens_dCompConc[a][0], fluidVisc[a][0], dFluidVisc_dPres[a][0], dFluidVisc_dCompConc[a][0], a);

  });

}

void ProppantTransport::UpdateProppantMobility(Group * const dataGroup)
{
  GEOSX_MARK_FUNCTION;

  arrayView1d<real64 const> const & conc = dataGroup->getReference<array1d<real64>>( viewKeyStruct::proppantConcentrationString );

  arrayView1d<real64 const> const & aperture = dataGroup->getReference<array1d<real64>>( FaceElementSubRegion::viewKeyStruct::elementApertureString );  

  arrayView1d<integer> const & isProppantMobile = dataGroup->getReference<array1d<integer>>( viewKeyStruct::isProppantMobileString );
  
  // TODO replace with batch update
  
  forall_in_range<RAJA::seq_exec>( 0, dataGroup->size(), [=] ( localIndex const a )
  {

    if(aperture[a] > m_minAperture && conc[a] < m_maxProppantConcentration)
      isProppantMobile[a] = 1;
    else
      isProppantMobile[a] = 0;      

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
  
  CommunicationTools::SynchronizeFields( fieldNames, mesh, domain->getNeighbors() );

  ResetViews( domain );  

  // We have to redo the below loop after fractures are generated

  localIndex const NC = m_numComponents;
  
  applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegionBase * const GEOSX_UNUSED_PARAM( region ),
                                 ElementSubRegionBase * const subRegion )
  {
    UpdateState( subRegion );

    arrayView1d<real64> const & proppantConc = m_proppantConcentration[er][esr];
    arrayView1d<real64> const & proppantConcOld = m_proppantConcentrationOld[er][esr];

    arrayView3d<real64> const & componentDens = m_componentDensity[er][esr][m_fluidIndex];
    arrayView2d<real64> const & componentDensOld = m_componentDensityOld[er][esr];    


    forall_in_range<serialPolicy>( 0, subRegion->size(), [=] ( localIndex ei )
    {
        proppantConcOld[ei] = proppantConc[ei];

        for(localIndex c = 0; c < NC; ++c)
          componentDensOld[ei][c] = componentDens[ei][0][c];        
    });

  } );

  m_downVector = gravityVector();
  m_downVector.Normalize();

  m_minAperture = m_bridgingFactor * m_proppantDiameter;

  m_proppantPackPermeability = pow(m_proppantDiameter, 2.0) / 180.0 * ((1.0 - m_maxProppantConcentration) * (1.0 - m_maxProppantConcentration) * (1.0 - m_maxProppantConcentration))/(m_maxProppantConcentration * m_maxProppantConcentration);
  
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

    FieldSpecificationManager const & boundaryConditionManager = FieldSpecificationManager::get();

    boundaryConditionManager.ApplyInitialConditions( domain );

    localIndex const NC = m_numComponents;
    
    applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegionBase * const GEOSX_UNUSED_PARAM( region ),
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

      arrayView1d<R1Tensor> const & cellBasedFlux   = m_cellBasedFlux[er][esr];

      arrayView1d<real64> const & proppantLiftFlux   = m_proppantLiftFlux[er][esr];      

      arrayView1d<real64> const & packVf = m_proppantPackVolumeFraction[er][esr];

      arrayView1d<real64> const & poroMultiplier = m_poroMultiplier[er][esr];
      arrayView1d<R1Tensor> const & transTMultiplier = m_transTMultiplier[er][esr];            

      arrayView1d<real64> const & excessPackV = m_proppantExcessPackVolume[er][esr];            

      arrayView1d<integer> const & isInterfaceElement = m_isInterfaceElement[er][esr];

      arrayView1d<integer> const & isProppantMobile = m_isProppantMobile[er][esr];          

      
      forall_in_range<serialPolicy>( 0, subRegion->size(), [=] ( localIndex ei )
      {
        dProppantConc[ei] = 0.0;

        proppantConcOld[ei] = proppantConc[ei];


        
        for(localIndex c = 0; c < NC; ++c)
          {

            dComponentConc[ei][c] = 0.0;        
            componentDensOld[ei][c] = componentDens[ei][0][c];

          }

        cellBasedFlux[ei] = 0.0;
        proppantLiftFlux[ei] = 0.0;
        
        packVf[ei] = 0.0;
        excessPackV[ei] = 0.0;

        poroMultiplier[ei] = 1.0;        
        transTMultiplier[ei] = 1.0;
            
        isInterfaceElement[ei] = 0;
        isProppantMobile[ei] = 1;                

      } );

    });
  }

  applyToSubRegions( mesh, [&] ( ElementSubRegionBase * const subRegion )
    {
       UpdateProppantMobility( subRegion );
    });


  UpdateCellBasedFlux(time_n, domain);
  
      
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

  applyToSubRegions( mesh, [&] ( ElementSubRegionBase * const subRegion )
    {
       UpdateProppantMobility( subRegion );
    });

  if(m_updateProppantPacking == 1)  
    UpdateProppantPackVolume(time_n, dt_return, domain);

  
  return dt_return;

}


void ProppantTransport::PreStepUpdate( real64 const& time,
                                       real64 const& dt,
                                       const int cycleNumber,
                                       DomainPartition * domain )
{
  GEOSX_MARK_FUNCTION;

  FlowSolverBase::PrecomputeData(domain);

  MeshLevel * mesh = domain->getMeshBody(0)->getMeshLevel(0);

  NodeManager const * const nodeManager = mesh->getNodeManager();
  FaceManager const * const faceManager = mesh->getFaceManager();

  
  if(cycleNumber == 0) {

    /*  assign intitial and boundary conditions */
    /*
    FieldSpecificationManager const * boundaryConditionManager = FieldSpecificationManager::get();

    boundaryConditionManager->ApplyInitialConditions( domain );
    */
    
    /* Below must be called after ImplicitStepSetup */

    ResizeFractureFields( time, dt, domain);
    
    localIndex const NC = m_numComponents;
    
    applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegionBase * const GEOSX_UNUSED_PARAM( region ),
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

      arrayView1d<R1Tensor> const & cellBasedFlux   = m_cellBasedFlux[er][esr];

      arrayView1d<real64> const & proppantLiftFlux   = m_proppantLiftFlux[er][esr];      

      arrayView1d<real64> const & packVf = m_proppantPackVolumeFraction[er][esr];

      arrayView1d<real64> const & excessPackV = m_proppantExcessPackVolume[er][esr];            

      arrayView1d<real64> const & poroMultiplier = m_poroMultiplier[er][esr];
      arrayView1d<R1Tensor> const & transTMultiplier = m_transTMultiplier[er][esr];            

      
      arrayView1d<integer> const & isInterfaceElement = m_isInterfaceElement[er][esr];

      arrayView1d<integer> const & isProppantMobile = m_isProppantMobile[er][esr];          

       forall_in_range<serialPolicy>( 0, subRegion->size(), [=] ( localIndex ei )
      {
        dProppantConc[ei] = 0.0;

        proppantConcOld[ei] = proppantConc[ei];

        for(localIndex c = 0; c < NC; ++c)
          {

            dComponentConc[ei][c] = 0.0;        
            componentDensOld[ei][c] = componentDens[ei][0][c];

          }

        cellBasedFlux[ei] = 0.0;
        proppantLiftFlux[ei] = 0.0;
        
        packVf[ei] = 0.0;
        excessPackV[ei] = 0.0;

        poroMultiplier[ei] = 1.0;        
        transTMultiplier[ei] = 1.0;            

        isInterfaceElement[ei] = 0;
        isProppantMobile[ei] = 1;                

      } );

    });
  }

  applyToSubRegions( mesh, [&] ( ElementSubRegionBase * const subRegion )
    {
       UpdateProppantMobility( subRegion );
       UpdateState( subRegion );

    });


  UpdateCellBasedFlux(time, domain);

}

void ProppantTransport::PostStepUpdate( real64 const & time_n,
                                        real64 const & dt_return,
                                        const int,
                                        DomainPartition * domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel * mesh = domain->getMeshBody(0)->getMeshLevel(0);
  
  applyToSubRegions( mesh, [&] ( ElementSubRegionBase * const subRegion )
    {
       UpdateProppantMobility( subRegion );
    });

  if(m_updateProppantPacking == 1)
    UpdateProppantPackVolume(time_n, dt_return, domain);

}
  

void ProppantTransport::ImplicitStepSetup( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                           real64 const & GEOSX_UNUSED_PARAM( dt ),
                                           DomainPartition * const domain,
                                           DofManager & GEOSX_UNUSED_PARAM(dofManager),
                                           ParallelMatrix & GEOSX_UNUSED_PARAM(matrix),
                                           ParallelVector & GEOSX_UNUSED_PARAM(rhs),
                                           ParallelVector & GEOSX_UNUSED_PARAM(solution) )
{

  localIndex const NC = m_numComponents;  

  ResetViews( domain );

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  /* The loop below could be moved to SolverStep after ImplicitStepSetup */
  
  applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegionBase * const GEOSX_UNUSED_PARAM( region ),
                                 ElementSubRegionBase * const subRegion )
  {

    arrayView1d<real64> const & dProppantConc   = m_deltaProppantConcentration[er][esr];
    arrayView2d<real64> const & dComponentConc   = m_deltaComponentConcentration[er][esr];        
    
    arrayView1d<real64> const & proppantConc = m_proppantConcentration[er][esr];
    arrayView1d<real64> const & proppantConcOld = m_proppantConcentrationOld[er][esr];

    arrayView3d<real64> const & componentDens = m_componentDensity[er][esr][m_fluidIndex];
    arrayView2d<real64> const & componentDensOld = m_componentDensityOld[er][esr];    

    arrayView1d<real64> const & excessPackV = m_proppantExcessPackVolume[er][esr];

    arrayView1d<real64> const & proppantLiftFlux = m_proppantLiftFlux[er][esr];    

    arrayView1d<R1Tensor> const & cellBasedFlux = m_cellBasedFlux[er][esr];    
    
    forall_in_range<serialPolicy>( 0, subRegion->size(), [=] ( localIndex ei )
    {

      dProppantConc[ei] = 0.0;
      proppantConcOld[ei] = proppantConc[ei];

      for(localIndex c = 0; c < NC; ++c)
        {

          dComponentConc[ei][c] = 0.0;                      
          componentDensOld[ei][c] = componentDens[ei][0][c];            

        }

      excessPackV[ei] = 0.0;
      proppantLiftFlux[ei] = 0.0;      
      cellBasedFlux[ei] = 0.0;
      
    } );
  } );

}

void ProppantTransport::ImplicitStepComplete( real64 const & GEOSX_UNUSED_PARAM(time_n),
                                              real64 const & GEOSX_UNUSED_PARAM(dt),
                                            DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  localIndex const NC = m_numComponents;
  
  applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegionBase * const GEOSX_UNUSED_PARAM( region ),
                                 ElementSubRegionBase * const subRegion )
  {
    arrayView1d<real64> const & proppantConc = m_proppantConcentration[er][esr];
    arrayView1d<real64 const> const & dProppantConc = m_deltaProppantConcentration[er][esr];    

    arrayView2d<real64> const & componentConc = m_componentConcentration[er][esr];
    arrayView2d<real64 const> const & dComponentConc = m_deltaComponentConcentration[er][esr];

    arrayView1d<real64> const & proppantLiftFlux = m_proppantLiftFlux[er][esr];        
    
    forall_in_range<serialPolicy>( 0, subRegion->size(), [=] ( localIndex ei )
    {
      proppantConc[ei] += dProppantConc[ei];
      proppantLiftFlux[ei] = 0.0;

      for(localIndex c = 0; c < NC; ++c)
        componentConc[ei][c] += dComponentConc[ei][c];
      
    } );

  } );

}

void ProppantTransport::SetupDofs( DomainPartition const * const GEOSX_UNUSED_PARAM(domain),
                                   DofManager & dofManager ) const
{
  dofManager.addField( viewKeyStruct::proppantConcentrationString,
                       DofManager::Location::Elem,
                       m_numDofPerCell,
                       m_targetRegions );

  dofManager.addCoupling( viewKeyStruct::proppantConcentrationString,
                          viewKeyStruct::proppantConcentrationString,
                          DofManager::Connector::Face );
  
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

  if( getLogLevel() >= 3 )
  {
    GEOSX_LOG_RANK_0("After ProppantTransport::AssembleAccumulationTerms");
    GEOSX_LOG_RANK_0("\nJacobian:\n");
    matrix.print(std::cout);
    GEOSX_LOG_RANK_0("\nResidual:\n");
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

  if( getLogLevel() >= 3 )
  {
    GEOSX_LOG_RANK_0("After ProppantTransport::AssembleSystem");
    GEOSX_LOG_RANK_0("\nJacobian:\n");
    matrix.print(std::cout);
    GEOSX_LOG_RANK_0("\nResidual:\n");
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
                                 ElementRegionBase const * const GEOSX_UNUSED_PARAM( region ),
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

    arrayView1d<real64 const> const & proppantPackVf          = m_proppantPackVolumeFraction[er][esr];            

    forall_in_range<serialPolicy>( 0, subRegion->size(), [=] ( localIndex ei )
    {
      if (elemGhostRank[ei] < 0)
      {

        stackArray1d<globalIndex, MAX_NUM_COMPONENTS>         localAccumDOF( m_numDofPerCell );
        stackArray1d<real64, MAX_NUM_COMPONENTS>             localAccum( m_numDofPerCell );
        stackArray2d<real64, MAX_NUM_COMPONENTS * MAX_NUM_COMPONENTS> localAccumJacobian( m_numDofPerCell, m_numDofPerCell );

        real64 effectiveVolume = volume[ei];
        real64 packPoreVolume = 0.0;        

        if(proppantPackVf[ei] < 1.0)
          {
            effectiveVolume = volume[ei]  * (1.0 - proppantPackVf[ei]);
            packPoreVolume = volume[ei]  * proppantPackVf[ei] * (1.0 - m_maxProppantConcentration);

          }

        AccumulationKernel::Compute(NC,
                                    proppantConcOld[ei],
                                    proppantConc[ei] + dProppantConc[ei],
                                    componentDensOld[ei],
                                    componentDens[ei][0],
                                    dCompDens_dPres[ei][0],
                                    dCompDens_dCompConc[ei][0],
                                    effectiveVolume,
                                    packPoreVolume,
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


void ProppantTransport::AssembleFluxTerms( real64 const GEOSX_UNUSED_PARAM(time_n),
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

  ElementRegionManager::ElementViewAccessor< arrayView1d<globalIndex const> >
  dofNumberAccessor = elemManager->ConstructViewAccessor< array1d<globalIndex>, arrayView1d<globalIndex const> >( dofKey );

  FluxKernel::ElementViewConst< arrayView1d<globalIndex const> > const & dofNumber = dofNumberAccessor.toViewConst();

  FluxKernel::ElementViewConst < arrayView1d<real64 const> > const & pres        = m_pressure.toViewConst();
  
  FluxKernel::ElementViewConst < arrayView1d<real64 const> > const & dPres       = m_deltaPressure.toViewConst();

  FluxKernel::ElementViewConst < arrayView1d<real64 const> > const & proppantConc       = m_proppantConcentration.toViewConst();

  FluxKernel::ElementViewConst < arrayView1d<real64 const> > const & dProppantConc       = m_deltaProppantConcentration.toViewConst();

  FluxKernel::ElementViewConst < arrayView1d<real64 const> > const & gravCoef   = m_gravCoef.toViewConst();

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

  FluxKernel::ElementViewConst < arrayView1d<integer const> > const & isProppantMobile  = m_isProppantMobile.toViewConst();

  FluxKernel::ElementViewConst < arrayView1d<R1Tensor const> > const & transTMultiplier  = m_transTMultiplier.toViewConst();    
  
  FluxKernel::ElementViewConst < arrayView1d<real64 const> > const & aperture  = m_elementAperture.toViewConst();

  FluxKernel::ElementViewConst < arrayView1d<real64 const> > const & proppantPackVf  = m_proppantPackVolumeFraction.toViewConst();

  FluxKernel::ElementViewConst < arrayView1d<real64 const> > const & proppantLiftFlux  = m_proppantLiftFlux.toViewConst();  

  FluxKernel::ElementViewConst < arrayView1d<integer const> > const & isInterfaceElement  = m_isInterfaceElement.toViewConst();  

  //  FluxKernel::ElementViewConst < arrayView1d<integer const> > const & elemGhostRank = m_elemGhostRank.toViewConst();

  localIndex const fluidIndex = m_fluidIndex;
  localIndex const proppantIndex = m_proppantIndex;  

  fluxApprox->forCellStencils( [&]( auto const & stencil )
  {

    FluxKernel::Launch( stencil,
                        m_numDofPerCell,
                        dt,
                        fluidIndex,
                        proppantIndex,
                        transTMultiplier,
                        m_updateProppantPacking,                        
                        m_downVector,
                        dofNumber,
                        pres,
                        dPres,
                        proppantConc,
                        dProppantConc,
                        componentDens,
                        dComponentDens_dPres,
                        dComponentDens_dComponentConc,
                        gravCoef,
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
                        proppantPackVf,
                        aperture,
                        proppantLiftFlux,
                        isInterfaceElement,
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

  matrix.open();
  rhs.open();

  FieldSpecificationManager & fsManager = FieldSpecificationManager::get();
  string const dofKey = dofManager.getKey( viewKeyStruct::proppantConcentrationString );


  //  Apply Dirichlet BC for proppant concentration

 fsManager.Apply( time_n + dt, domain, "ElementRegions", viewKeyStruct::proppantConcentrationString,
                    [&]( FieldSpecificationBase const * const fs,
                    string const &,
                    SortedArrayView<localIndex const> const & lset,
                    Group * subRegion,
                    string const & ) -> void
  {


    arrayView1d<globalIndex const> const &
    dofNumber = subRegion->getReference< array1d<globalIndex> >( dofKey );

    arrayView1d<real64 const> const &
    proppantConc = subRegion->getReference<array1d<real64> >( viewKeyStruct::proppantConcentrationString );

    arrayView1d<real64 const> const &
    dProppantConc = subRegion->getReference<array1d<real64> >( viewKeyStruct::deltaProppantConcentrationString );

        
    fs->ApplyBoundaryConditionToSystem<FieldSpecificationEqual, LAInterface>( lset,
                                                                              time_n + dt,
                                                                              subRegion,
                                                                              dofNumber,
                                                                              m_numDofPerCell,
                                                                              matrix,
                                                                              rhs,
                                                                              [&]( localIndex const a) -> real64
        {
          return proppantConc[a] + dProppantConc[a];
        });
  
  });
 
 //  Apply Dirichlet BC for component concentration

  localIndex const NC = m_numComponents;

  if(NC > 0)
    {

      map< string, map< string, array1d<bool> > > bcStatusMap; // map to check consistent application of BC

      fsManager.Apply( time_n + dt,
                       domain,
                       "ElementRegions",
                       viewKeyStruct::proppantConcentrationString,
                       [&]( FieldSpecificationBase const * const GEOSX_UNUSED_PARAM( fs ),
                            string const & setName,
                            SortedArrayView<localIndex const> const & GEOSX_UNUSED_PARAM( targetSet ),
                        Group * subRegion,
                        string const & )
      {

        string const & subRegionName = subRegion->getName();
        GEOSX_ERROR_IF( bcStatusMap[subRegionName].count( setName ) > 0, "Conflicting proppant boundary conditions on set " << setName );
        bcStatusMap[subRegionName][setName].resize( NC );
        bcStatusMap[subRegionName][setName] = false;

      });

      fsManager.Apply( time_n + dt,
                       domain,
                       "ElementRegions",
                       viewKeyStruct::componentConcentrationString,
                       [&] ( FieldSpecificationBase const * const fs,
                             string const & setName,
                             SortedArrayView<localIndex const> const & targetSet,
                             Group * subRegion,
                             string const & )      
      {

         string const & subRegionName = subRegion->getName();
         localIndex const comp = fs->GetComponent();

         GEOSX_ERROR_IF( bcStatusMap[subRegionName].count( setName ) == 0, "Proppant boundary condition not prescribed on set '" << setName << "'" );
         GEOSX_ERROR_IF( bcStatusMap[subRegionName][setName][comp], "Conflicting composition[" << comp << "] boundary conditions on set '" << setName << "'" );
         bcStatusMap[subRegionName][setName][comp] = true;
         
         fs->ApplyFieldValue<FieldSpecificationEqual>( targetSet,
                                                       time_n + dt,
                                                       subRegion,
                                                       viewKeyStruct::bcComponentConcentrationString );

      });
    
      bool bcConsistent = true;
      for (auto const & bcStatusEntryOuter : bcStatusMap)
        {
          for( auto const & bcStatusEntryInner : bcStatusEntryOuter.second )
            {
              for( localIndex ic = 0 ; ic < NC ; ++ic )
                {
                  bcConsistent &= bcStatusEntryInner.second[ic];
                  GEOSX_WARNING_IF( !bcConsistent, "Composition boundary condition not applied to component " << ic
                                    << " on region '" << bcStatusEntryOuter.first << "',"
                                    << " set '" << bcStatusEntryInner.first << "'" );
                }
            }
        }

      GEOSX_ERROR_IF( !bcConsistent, "Inconsistent composition boundary conditions" );
    
      fsManager.Apply( time_n + dt,
                       domain,
                       "ElementRegions",
                       viewKeyStruct::proppantConcentrationString,
                       [&] ( FieldSpecificationBase const * const GEOSX_UNUSED_PARAM( bc ),
                         string const & GEOSX_UNUSED_PARAM( setName ),
                         SortedArrayView<localIndex const> const & targetSet,
                         Group * subRegion,
                         string const & )
      {      

        arrayView1d<globalIndex const> const & dofNumber = subRegion->getReference< array1d<globalIndex> >( dofKey );

        arrayView2d<real64 const> const & compConc  = subRegion->getReference< array2d<real64> >( viewKeyStruct::componentConcentrationString );

        arrayView2d<real64 const> const & deltaCompConc  = subRegion->getReference< array2d<real64> >( viewKeyStruct::deltaComponentConcentrationString );

        arrayView2d<real64 const> const & bcCompConc  = subRegion->getReference< array2d<real64> >( viewKeyStruct::bcComponentConcentrationString );        

        array1d<real64> rhsContribution( targetSet.size() * NC );
        array1d<globalIndex> dof( targetSet.size() * NC );

        integer counter = 0;

        for (localIndex a : targetSet)
          {

            for (localIndex ic = 0; ic < NC; ++ic)
              {
                
                dof[counter] = dofNumber[a] + ic + 1;

                FieldSpecificationEqual::SpecifyFieldValue<LAInterface>( dof[counter],
                                                                         matrix,
                                                                         rhsContribution[counter],
                                                                         bcCompConc[a][ic],
                                                                         compConc[a][ic] + deltaCompConc[a][ic] );

                ++counter;

              }

            FieldSpecificationEqual::PrescribeRhsValues<LAInterface>( rhs,
                                                                      counter,
                                                                      dof.data(),
                                                                      rhsContribution.data());


          }
            
      });

    }

  matrix.close();
  rhs.close();

  if( getLogLevel() == 2 )
  {
    GEOSX_LOG_RANK_0( "After ProppantTransport::AssembleSystem" );
    GEOSX_LOG_RANK_0("\nJacobian:\n");
    std::cout << matrix;
    GEOSX_LOG_RANK_0("\nResidual:\n");
    std::cout << rhs;
  }

  if( getLogLevel() >= 3 )
  {
    integer const newtonIter = m_nonlinearSolverParameters.m_numNewtonIterations;

    string filename_mat = "matrix_bc_" + std::to_string( time_n ) + "_" + std::to_string( newtonIter ) + ".mtx";
    matrix.write( filename_mat, LAIOutputFormat::MATRIX_MARKET );

    string filename_rhs = "rhs_bc_" + std::to_string( time_n ) + "_" + std::to_string( newtonIter ) + ".mtx";
    rhs.write( filename_rhs, LAIOutputFormat::MATRIX_MARKET );

    GEOSX_LOG_RANK_0( "After ProppantTransport::ApplyBoundaryConditions" );
    GEOSX_LOG_RANK_0( "Jacobian: written to " << filename_mat );
    GEOSX_LOG_RANK_0( "Residual: written to " << filename_rhs );
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
                                 ElementRegionBase const * const GEOSX_UNUSED_PARAM( region ),
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


  dofManager.addVectorToField( solution,
                               viewKeyStruct::proppantConcentrationString,
                               viewKeyStruct::deltaProppantConcentrationString,
                               scalingFactor,
                               0, 1 );


  if(m_numDofPerCell > 1)
    dofManager.addVectorToField( solution,
                                 viewKeyStruct::proppantConcentrationString,
                                 viewKeyStruct::deltaComponentConcentrationString,
                                 scalingFactor,
                                 1, m_numDofPerCell );

  std::map<string, string_array > fieldNames;
  fieldNames["elems"].push_back( viewKeyStruct::deltaProppantConcentrationString );
  fieldNames["elems"].push_back( viewKeyStruct::deltaComponentConcentrationString );  

  CommunicationTools::SynchronizeFields( fieldNames, mesh, domain->getNeighbors() );

  applyToSubRegions( mesh, [&] ( ElementSubRegionBase * subRegion )
  {
    //    UpdateState( subRegion );
    UpdateComponentDensity( subRegion );    
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
  
  if( getLogLevel() >= 2 )
  {
    GEOSX_LOG_RANK("After ProppantTransport::SolveSystem");
    GEOSX_LOG_RANK_0("\nsolution:\n");
    std::cout << solution;
  }

}

void ProppantTransport::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  MeshLevel * mesh = domain->getMeshBody(0)->getMeshLevel(0);

  localIndex const NC = m_numComponents;
  
  applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegionBase * const GEOSX_UNUSED_PARAM( region ),
                                 ElementSubRegionBase * const subRegion )
  {
    arrayView1d<real64> const & dProppantConc = m_deltaProppantConcentration[er][esr];
    arrayView2d<real64> const & dComponentConc = m_deltaComponentConcentration[er][esr];        

    forall_in_range<serialPolicy>( 0, subRegion->size(), [=] ( localIndex ei )
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

  m_cellBasedFlux =
    elemManager->ConstructViewAccessor<array1d<R1Tensor>, arrayView1d<R1Tensor>>( viewKeyStruct::cellBasedFluxString );

  m_proppantPackVolumeFraction =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( viewKeyStruct::proppantPackVolumeFractionString );

  m_proppantExcessPackVolume =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( viewKeyStruct::proppantExcessPackVolumeString );            

  m_proppantLiftFlux =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( viewKeyStruct::proppantLiftFluxString );

  m_isProppantMobile =
    elemManager->ConstructViewAccessor<array1d<integer>, arrayView1d<integer>>( viewKeyStruct::isProppantMobileString );              

  m_isInterfaceElement =
    elemManager->ConstructViewAccessor<array1d<integer>, arrayView1d<integer>>( viewKeyStruct::isInterfaceElementString );              

  m_isProppantBoundaryElement =
    elemManager->ConstructViewAccessor<array1d<integer>, arrayView1d<integer>>( viewKeyStruct::isProppantBoundaryString );              
  
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

  m_fluidViscosity = 
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64> >( SlurryFluidBase::viewKeyStruct::fluidViscosityString, constitutiveManager );
  
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

  m_poroMultiplier =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( viewKeyStruct::poroMultiplierString );

  m_transTMultiplier =
    elemManager->ConstructViewAccessor<array1d<R1Tensor>, arrayView1d<R1Tensor>>( viewKeyStruct::transTMultiplierString );  
  
}



void ProppantTransport::UpdateCellBasedFlux( real64 const GEOSX_UNUSED_PARAM(time_n),
                                             DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;
  
  MeshLevel * mesh = domain->getMeshBody(0)->getMeshLevel(0);  

  NumericalMethodsManager const * numericalMethodManager =
    domain->getParent()->GetGroup<NumericalMethodsManager>( keys::numericalMethodsManager );

  FiniteVolumeManager const * fvManager =
    numericalMethodManager->GetGroup<FiniteVolumeManager>( keys::finiteVolumeManager );

  FluxApproximationBase const * fluxApprox = fvManager->getFluxApproximation( m_discretizationName );

  FluxKernel::ElementViewConst < arrayView1d<real64 const> > const & pres        = m_pressure.toViewConst();
  
  //  FluxKernel::ElementViewConst < arrayView1d<real64 const> > const & dPres       = m_deltaPressure.toViewConst();

  FluxKernel::ElementViewConst < arrayView1d<real64 const> > const & gravCoef   = m_gravCoef.toViewConst();

  FluxKernel::MaterialView< arrayView2d<real64 const> > const & dens        = m_density.toViewConst();

  FluxKernel::MaterialView< arrayView2d<real64 const> > const & visc        = m_viscosity.toViewConst();

  FluxKernel::ElementViewConst < arrayView1d<real64 const> > const & aperture  = m_elementAperture.toViewConst();

  FluxKernel::ElementView  < arrayView1d<R1Tensor> >  & cellBasedFlux  = m_cellBasedFlux.toView();  

  FluxKernel::ElementView < arrayView1d<real64 const> > const & proppantPackVf  = m_proppantPackVolumeFraction.toViewConst();

  FluxKernel::ElementView < arrayView1d<R1Tensor const> > const & transTMultiplier  = m_transTMultiplier.toViewConst();  

  localIndex const fluidIndex = m_fluidIndex;

  fluxApprox->forCellStencils( [&]( auto const & stencil )
  {

    FluxKernel::LaunchCellBasedFluxCalculation( stencil,
                                                fluidIndex,
                                                transTMultiplier,
                                                m_downVector,
                                                pres,
                                                gravCoef,
                                                dens,
                                                visc,
                                                aperture,
                                                proppantPackVf,
                                                cellBasedFlux);
  });

  
  std::map<string, string_array > fieldNames;
  fieldNames["elems"].push_back( viewKeyStruct::cellBasedFluxString );  

  CommunicationTools::SynchronizeFields( fieldNames, mesh, domain->getNeighbors() );
  
}

void ProppantTransport::UpdateProppantPackVolume( real64 const GEOSX_UNUSED_PARAM(time_n),
                                                  real64 const dt,
                                                  DomainPartition * const domain)
{
  
  GEOSX_MARK_FUNCTION;

  MeshLevel * mesh = domain->getMeshBody(0)->getMeshLevel(0);  

  NumericalMethodsManager const * numericalMethodManager =
    domain->getParent()->GetGroup<NumericalMethodsManager>( keys::numericalMethodsManager );

  FiniteVolumeManager const * fvManager =
    numericalMethodManager->GetGroup<FiniteVolumeManager>( keys::finiteVolumeManager );

  FluxApproximationBase const * fluxApprox = fvManager->getFluxApproximation( m_discretizationName );

  ProppantPackVolumeKernel::ElementView < arrayView1d<real64> > const & conc       = m_proppantConcentration.toView();

  ProppantPackVolumeKernel::MaterialViewConst< arrayView1d<real64> > const & settlingFactor = m_settlingFactor.toViewConst();

  ProppantPackVolumeKernel::MaterialViewConst< arrayView2d<real64> > const & density = m_density.toViewConst();
  
  ProppantPackVolumeKernel::MaterialViewConst< arrayView2d<real64> > const & fluidDensity = m_fluidDensity.toViewConst();

  ProppantPackVolumeKernel::MaterialViewConst< arrayView2d<real64> > const & fluidViscosity = m_fluidViscosity.toViewConst();  
  
  ProppantPackVolumeKernel::ElementView < arrayView1d<integer> > const & isProppantMobile = m_isProppantMobile.toView();

  ProppantPackVolumeKernel::ElementViewConst < arrayView1d<integer> > const & isProppantBoundaryElement = m_isProppantBoundaryElement.toViewConst();  
  
  ProppantPackVolumeKernel::ElementView < arrayView1d<real64> > const & proppantPackVf  = m_proppantPackVolumeFraction.toView();

  ProppantPackVolumeKernel::ElementView < arrayView1d<real64> > const & proppantExcessPackV  = m_proppantExcessPackVolume.toView();

  ProppantPackVolumeKernel::ElementView < arrayView1d<integer> > const & isInterfaceElement = m_isInterfaceElement.toView();

  ProppantPackVolumeKernel::ElementView < arrayView1d<R1Tensor> > const & cellBasedFlux  = m_cellBasedFlux.toView();

  ProppantPackVolumeKernel::ElementView < arrayView1d<real64> > const & proppantLiftFlux  = m_proppantLiftFlux.toView();    

  ProppantPackVolumeKernel::ElementViewConst < arrayView1d<real64 const> > const & volume  = m_volume.toViewConst();

  ProppantPackVolumeKernel::ElementViewConst < arrayView1d<real64 const> > const & aperture  = m_elementAperture.toViewConst();

  ProppantPackVolumeKernel::ElementViewConst < arrayView1d<integer const> > const & elemGhostRank = m_elemGhostRank.toViewConst();

  
  localIndex const fluidIndex = m_fluidIndex;
  localIndex const proppantIndex = m_proppantIndex;    

  fluxApprox->forCellStencils( [&]( auto const & stencil )
  {

    ProppantPackVolumeKernel::LaunchProppantPackVolumeCalculation( stencil,
                                                                   dt,
                                                                   fluidIndex,
                                                                   proppantIndex,
                                                                   m_proppantDensity,
                                                                   m_proppantDiameter,
                                                                   m_maxProppantConcentration,
                                                                   m_downVector,
                                                                   m_criticalShieldsNumber,
                                                                   m_frictionCoefficient,
                                                                   conc,
                                                                   settlingFactor,
                                                                   density,
                                                                   fluidDensity,
                                                                   fluidViscosity,
                                                                   isProppantMobile,
                                                                   isProppantBoundaryElement,                                                                   
                                                                   proppantPackVf,
                                                                   proppantExcessPackV,
                                                                   aperture,
                                                                   volume,
                                                                   elemGhostRank,
                                                                   cellBasedFlux,
                                                                   proppantLiftFlux);

  });


  {
  
    std::map<string, string_array > fieldNames;
    fieldNames["elems"].push_back( viewKeyStruct::proppantConcentrationString );  
    fieldNames["elems"].push_back( viewKeyStruct::proppantPackVolumeFractionString );
    fieldNames["elems"].push_back( viewKeyStruct::proppantExcessPackVolumeString );
    fieldNames["elems"].push_back( viewKeyStruct::proppantLiftFluxString );      

    CommunicationTools::SynchronizeFields( fieldNames, mesh, domain->getNeighbors() );

  }

  applyToSubRegions( mesh, [&] ( ElementSubRegionBase * const subRegion )
    {
       UpdateProppantMobility( subRegion );
    });

    
  fluxApprox->forCellStencils( [&]( auto const & stencil )
  {

    ProppantPackVolumeKernel::LaunchProppantPackVolumeUpdate( stencil,
                                                              m_downVector,
                                                              m_maxProppantConcentration,
                                                              conc,
                                                              isProppantMobile,
                                                              proppantPackVf,
                                                              proppantExcessPackV);
  });


  {
  
    std::map<string, string_array > fieldNames;
    fieldNames["elems"].push_back( viewKeyStruct::proppantConcentrationString );  
    fieldNames["elems"].push_back( viewKeyStruct::proppantPackVolumeFractionString );

    CommunicationTools::SynchronizeFields( fieldNames, mesh, domain->getNeighbors() );

  }

  applyToSubRegions( mesh, [&] ( ElementSubRegionBase * const subRegion )
    {
      UpdateProppantMobility( subRegion );
    });


  
  fluxApprox->forCellStencils( [&]( auto const & stencil )
  {

    ProppantPackVolumeKernel::LaunchInterfaceElementUpdate( stencil,
                                                            m_downVector,
                                                            isProppantMobile,
                                                            isInterfaceElement);
  });

  {
  
    std::map<string, string_array > fieldNames;
    fieldNames["elems"].push_back( viewKeyStruct::isInterfaceElementString );  

    CommunicationTools::SynchronizeFields( fieldNames, mesh, domain->getNeighbors() );

  }

  // update poroMultiplier and transTMultiplier

  applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegionBase * const GEOSX_UNUSED_PARAM( region ),
                                 ElementSubRegionBase * const subRegion )
  {

    arrayView1d<real64 const> const & proppantPackVfNew = m_proppantPackVolumeFraction[er][esr];

    arrayView1d<real64 const> const & aperture0 = m_elementAperture[er][esr];        
    
    arrayView1d<real64> const & poroMultiplier = m_poroMultiplier[er][esr];
    arrayView1d<R1Tensor> const & transTMultiplier = m_transTMultiplier[er][esr];        

    forall_in_range<serialPolicy>( 0, subRegion->size(), [=] ( localIndex ei )
    {

      poroMultiplier[ei] = 1.0 - m_maxProppantConcentration * proppantPackVfNew[ei];

      //K0 horizontal

      transTMultiplier[ei][0] = (1.0 - proppantPackVfNew[ei]) + proppantPackVfNew[ei] * m_proppantPackPermeability * 12.0 / (aperture0[ei] * aperture0[ei]);

      //K1 vertical

      transTMultiplier[ei][1] = 1.0 / (proppantPackVfNew[ei] * aperture0[ei] * aperture0[ei] / 12.0 / m_proppantPackPermeability + (1.0 - proppantPackVfNew[ei]));

    } );
  
  });
  
}


REGISTER_CATALOG_ENTRY( SolverBase, ProppantTransport, std::string const &, Group * const )
} /* namespace geosx */
