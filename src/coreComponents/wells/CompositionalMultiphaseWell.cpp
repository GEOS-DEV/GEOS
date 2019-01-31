/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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

/*
 * @file CompositionalMultiphaseWell.cpp
 *
 */

#include "CompositionalMultiphaseWell.hpp"
#include "Segment.hpp"
#include "Connection.hpp"
#include "Perforation.hpp"

#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/Fluid/SingleFluidBase.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "math/TensorT/TensorBaseT.h"
#include "systemSolverInterface/LinearSolverWrapper.hpp"
#include "systemSolverInterface/EpetraBlockSystem.hpp"


namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace systemSolverInterface;
  
CompositionalMultiphaseWell::CompositionalMultiphaseWell(string const & name, dataRepository::ManagedGroup * const parent)
  : WellBase( name, parent ),
    m_bhp(),
    m_numPhases(),
    m_numComponents()
{
  RegisterViewWrapper( viewKeyStruct::bhpString, &m_bhp, false );

  // Most segment-based fields are registered by specific well models, rather than
  // SegmentManager, ConnectionManager, and PerforationManager that are is physics-agnostic

  // Segment-based primary variables
  SegmentManager * segManager = GetGroup<SegmentManager>( groupKeyStruct::segmentsString );
  segManager->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::pressureString );
  segManager->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::deltaPressureString );
  segManager->RegisterViewWrapper<array2d<real64>>( viewKeyStruct::globalCompDensityString );
  segManager->RegisterViewWrapper<array2d<real64>>( viewKeyStruct::deltaGlobalCompDensityString );
  
  // Connection-based primary variables
  ConnectionManager * connManager = GetGroup<ConnectionManager>( groupKeyStruct::connectionsString );
  connManager->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::mixtureVelocityString );
  connManager->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::deltaMixtureVelocityString );

  // Connection-based input variables
  connManager->RegisterViewWrapper<array2d<real64>>( viewKeyStruct::globalComponentFracString );
  connManager->RegisterViewWrapper<array3d<real64>>( viewKeyStruct::dGlobalComponentFrac_dPresString );
  connManager->RegisterViewWrapper<array3d<real64>>( viewKeyStruct::dGlobalComponentFrac_dCompString );

  // Perforation-based quantities 
  PerforationManager * perfManager = GetGroup<PerforationManager>( groupKeyStruct::perforationsString );
  perfManager->RegisterViewWrapper<array2d<real64>>( viewKeyStruct::phaseFlowRateString );

  // more variables added later: secondary variables (density, viscosity, etc) for segment, flow rates for connections, etc
}

CompositionalMultiphaseWell::~CompositionalMultiphaseWell()
{

}

void CompositionalMultiphaseWell::InitializePostSubGroups( ManagedGroup * const rootGroup )
{
  WellBase::InitializePostSubGroups( rootGroup );

  // generate the "all" set to enable application of BC
  set<localIndex> & setAll = this->sets()->RegisterViewWrapper<set<localIndex>>("all")->reference();
  setAll.insert(0);
}

// Initialize all variables in the well
void CompositionalMultiphaseWell::InitializeState( DomainPartition * const domain )
{

}
  
// Apply well solution after the linear solve
void CompositionalMultiphaseWell::ApplySolution( systemSolverInterface::EpetraBlockSystem const * const blockSystem,
                                                 real64 const scalingFactor,
                                                 DomainPartition * const domain )
{
  Epetra_Map const * const rowMap        = blockSystem->GetRowMap( BlockIDs::compositionalBlock );
  Epetra_FEVector const * const solution = blockSystem->GetSolutionVector( BlockIDs::compositionalBlock );

  int dummy;
  double * local_solution = nullptr;
  solution->ExtractView( &local_solution, &dummy );

  // the pressure and compositions live on segments
  // the velocity lives on connections
  arrayView1d<real64> const & deltaPressure =
    m_segmentManager.getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );
  arrayView2d<real64> const & deltaGlobalCompDensity =
    m_segmentManager.getReference<array2d<real64>>( viewKeyStruct::deltaGlobalCompDensityString );
  arrayView1d<real64> const & deltaMixtureVelocity =
    m_connectionManager.getReference<array1d<real64>>( viewKeyStruct::deltaMixtureVelocityString );

  real64 dummy_newton_update = 0;
  
  // loop over the segments
  for (localIndex iseg = 0; iseg < numSegmentsLocal(); ++iseg)
  {
    // update pressure on segment iseg
    deltaPressure[iseg] += scalingFactor * dummy_newton_update;
    
    // update compositions on segment iseg
    for (localIndex ic = 0; ic < m_numComponents; ++ic)
    {
      deltaGlobalCompDensity[iseg][ic] += scalingFactor * dummy_newton_update;
    }
  }

  // loop over the connection
  for (localIndex iconn = 0; iconn < numConnectionsLocal(); ++iconn)
  {
    // update velocity
    deltaMixtureVelocity[iconn] += scalingFactor * dummy_newton_update;
  }
  
  // synchronize fields

  // update the fluid models
  StateUpdate( domain );
}
  
// Compute the rate for all the perforations and form the control equation
void CompositionalMultiphaseWell::AssembleWellTerms( DomainPartition * const domain,
                                                     EpetraBlockSystem * const blockSystem,
                                                     real64 const time_n,
                                                     real64 const dt )
{
  // we assemble n_c mass conservation equations per segment
  // we also need 1 volume balance equation per segment

  // 1- Accumulation term on segments
  for (localIndex iseg = 0; iseg < numSegmentsLocal(); ++iseg)
  {
    // 1.1- Compute accumulation term of mass conservation equations for segment iseg 

    // 1.2- Assemble volume balance equation
  }


  // 2- Flux term through connection
  for (localIndex iconn = 0; iconn < numConnectionsLocal(); ++iconn)
  {
    // 2.1- Perform upwinding of fluid dependent quantities at connection iconn
    
    // 2.2- Compute flux term at connection iconn

    // 2.3- Add flux term to residual and jacobian corresponding to segments nextSegment and prevSegment
  }
    
  // 3- Source/Sink term through perforations
  for (localIndex iperf = 0; iperf < numPerforationsLocal(); ++iperf)
  {
    // 3.1 Upwinding fluid dependent quantities at perforation iperf
    // (choose between reservoir and segment vars)
    
    // 3.2 Compute the rate at perforation iperf using segment data and reservoir data
    
    // 3.3 Add to residual and jacobian corresponding to segment iseg
    
  }
  
  FormControlEquation( domain, blockSystem, time_n, dt );
}

  
// Check if the well control needs to be switched
void CompositionalMultiphaseWell::CheckControlSwitch()
{
  // to defined later
}

// Reset the well to its initial state at the beginning of the time step
void CompositionalMultiphaseWell::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  // the pressure and compositions live on segments
  // the velocity lives on connections
  arrayView1d<real64> const & deltaPressure =
    m_segmentManager.getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );
  arrayView2d<real64> const & deltaGlobalCompDensity =
    m_segmentManager.getReference<array2d<real64>>( viewKeyStruct::deltaGlobalCompDensityString );
  arrayView1d<real64> const & deltaMixtureVelocity =
    m_connectionManager.getReference<array1d<real64>>( viewKeyStruct::deltaMixtureVelocityString );

  
  // loop over the segments
  for (localIndex iseg = 0; iseg < numSegmentsLocal(); ++iseg)
  {
    // set deltaPressure = 0 on segment iseg
    deltaPressure[iseg] = 0;
    
    // set deltaGlobalCompDensity = 0 on segment iseg
    for (localIndex ic = 0; ic < m_numComponents; ++ic)
    {
      deltaGlobalCompDensity[iseg][ic] = 0;
    }
  }

  // loop over the connection
  for (localIndex iconn = 0; iconn < numConnectionsLocal(); ++iconn)
  {
    // set deltaMixtureVelocity = 0 on connection iconn
    deltaMixtureVelocity[iconn] = 0;
  }
  
  // synchronize fields

  // update the fluid models
  StateUpdate( domain );
}

real64 CompositionalMultiphaseWell::GetTotalFlowRate()
{
  return 0.0;
}

// form the well control equation based on the type of well
void CompositionalMultiphaseWell::FormControlEquation( DomainPartition const * const domain,
						       EpetraBlockSystem * const blockSystem,
                                                       real64 const time_n,
                                                       real64 const dt )
{
  // to be defined later
}

// update each connection pressure from bhp and hydrostatic head
void CompositionalMultiphaseWell::StateUpdate( DomainPartition const * domain )
{
  // 1- loop over the segments
  for (localIndex iseg = 0; iseg < numSegmentsLocal(); ++iseg)
  {
    // 1.1 Update phase density and viscosity

    /* it would be nice to do something like:
      
       MultiFluidBase const * fluid = constitutiveRelations[segmentRegion]->group_cast<MultiFluidBase const *>();
       fluid->PointUpdate( pressure[iseg] );
 
       to update dens[iseg] 

     */
    
    // 1.2 Update average mixture density

    // 1.3 Update saturation

    // 1.4 Update relperm and capillary pressure   
  }

  // loop over the connections
  for (localIndex iconn = 0; iconn < numConnectionsLocal(); ++iconn)
  {
    // if (injector and firstConnection)
    //    set globalComponentFrac to injection stream
    
    // not much to update if the upwinding is computed in AssembleWellTerms
  }
}


REGISTER_CATALOG_ENTRY( WellBase, CompositionalMultiphaseWell, string const &, ManagedGroup * const )

} //namespace geosx
