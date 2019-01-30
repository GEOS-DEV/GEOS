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
    m_numPhases(),
    m_numComponents(),
    m_bhp()
{
  RegisterViewWrapper( viewKeyStruct::bhpString, &m_bhp, false );

  // Most segment-based fields are registered by specific well models, rather than
  // SegmentManager itself, which is physics-agnostic
  SegmentManager * segManager = GetGroup<SegmentManager>( groupKeyStruct::segmentsString );
  segManager->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::pressureString );
  segManager->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::deltaPressureString );
  segManager->RegisterViewWrapper<array2d<real64>>( viewKeyStruct::globalCompDensityString );
  segManager->RegisterViewWrapper<array2d<real64>>( viewKeyStruct::deltaGlobalCompDensityString );
  segManager->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::mixtureVelocityString );
  segManager->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::deltaMixtureVelocityString );

  segManager->RegisterViewWrapper<array2d<real64>>( viewKeyStruct::globalComponentFracString );
  segManager->RegisterViewWrapper<array3d<real64>>( viewKeyStruct::dGlobalComponentFrac_dPresString );
  segManager->RegisterViewWrapper<array3d<real64>>( viewKeyStruct::dGlobalComponentFrac_dCompString );

  PerforationManager * perfManager = GetGroup<PerforationManager>( groupKeyStruct::perforationsString );
  perfManager->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::flowRateString );
}

CompositionalMultiphaseWell::~CompositionalMultiphaseWell()
{

}

void CompositionalMultiphaseWell::InitializePostSubGroups( ManagedGroup * const rootGroup )
{
  WellBase::InitializePostSubGroups( rootGroup );

  // vars owned by the well itself are scalar, but must be stored as arrays for BC, so resize to 1
  resize(1);

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

  // loop over the segments
  for (localIndex iseg = 0; iseg < numSegmentsLocal(); ++iseg)
  {
    // TODO: check if "ghost segment"

    globalIndex const offset = 0;
    // extract solution and apply to deltaPressure
    {
      int const lid = rowMap->LID( integer_conversion<int>( offset ) );
      //deltaPressure[iseg] += scalingFactor * local_solution[lid];
      int const lid = rowMap->LID( integer_conversion<int>( offset + 1 ) );
      //deltaPressure[iseg] += scalingFactor * local_solution[lid];
    }

    for (localIndex ic = 0; ic < m_numComponents; ++ic)
    {
      int const lid = rowMap->LID( integer_conversion<int>( offset + ic + 2 ) );
      //deltaGlobalCompDensity[iseg][ic] += scalingFactor * local_solution[lid];
    }
  }
    
  // TODO: synchronize fields

  // update the fluid models
  StateUpdate( domain, fluidIndex );
}
  
// Compute the rate for all the perforations and form the control equation
void CompositionalMultiphaseWell::AssembleWellTerms( DomainPartition * const domain,
                                                     EpetraBlockSystem * const blockSystem,
                                                     real64 const time_n,
                                                     real64 const dt )
{
  
}

  
// Check if the well control needs to be switched
void CompositionalMultiphaseWell::CheckControlSwitch()
{

}

// Reset the well to its initial state at the beginning of the time ste
void CompositionalMultiphaseWell::ResetStateToBeginningOfStep( DomainPartition * const domain )
{

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
  
}

// update each connection pressure from bhp and hydrostatic head
void CompositionalMultiphaseWell::StateUpdate( DomainPartition const * domain, localIndex fluidIndex )
{
 
}


REGISTER_CATALOG_ENTRY( WellBase, CompositionalMultiphaseWell, string const &, ManagedGroup * const )

} //namespace geosx
