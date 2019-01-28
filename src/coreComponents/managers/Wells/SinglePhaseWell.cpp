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
 * @file SimpleWell.cpp
 *
 */

#include "SinglePhaseWell.hpp"
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
  
SinglePhaseWell::SinglePhaseWell(string const & name, dataRepository::ManagedGroup * const parent)
  : WellBase( name, parent ),
    m_bhp()
{
  RegisterViewWrapper( viewKeyStruct::bhpString, &m_bhp, false );

  // Most perforations-based fields are registered by specific well models, rather than
  // PerforationManager itself, which is physics-agnostic
  PerforationManager * perfManager = GetGroup<PerforationManager>( groupKeyStruct::perforationsString );
  
  perfManager->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::pressureString );
  perfManager->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::deltaPressureString );
  perfManager->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::flowRateString );

  perfManager->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::avgDensityString );
  
  perfManager->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::densityString );
  perfManager->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::dDensity_dPresString );
  perfManager->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::viscosityString );
  perfManager->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::dViscosity_dPresString );
  perfManager->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::mobilityString );
  perfManager->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::dMobility_dPresString );
}

SinglePhaseWell::~SinglePhaseWell()
{

}

void SinglePhaseWell::InitializePostSubGroups( ManagedGroup * const rootGroup )
{
  WellBase::InitializePostSubGroups( rootGroup );

  // vars owned by the well itself are scalar, but must be stored as arrays for BC, so resize to 1
  resize(1);

  // generate the "all" set to enable application of BC
  set<localIndex> & setAll = this->sets()->RegisterViewWrapper<set<localIndex>>("all")->reference();
  setAll.insert(0);
}

// Initialize all variables in the well
void SinglePhaseWell::InitializeState( DomainPartition * const domain )
{

}
  
// Apply well solution after the linear solve
void SinglePhaseWell::ApplySolution( systemSolverInterface::EpetraBlockSystem const * const blockSystem,
                                     real64 const scalingFactor,
                                     DomainPartition * const domain )
{

}

// Compute the rate for all the perforations and form the control equation
void SinglePhaseWell::AssembleWellTerms( DomainPartition * const domain,
                                         EpetraBlockSystem * const blockSystem,
                                         real64 const time_n,
                                         real64 const dt )
{
  
}
  

// Check if the well control needs to be switched
void SinglePhaseWell::CheckControlSwitch()
{

}

// Reset the well to its initial state at the beginning of the time ste
void SinglePhaseWell::ResetStateToBeginningOfStep( DomainPartition * const domain )
{

}

real64 SinglePhaseWell::GetTotalFlowRate()
{
  arrayView1d<real64 const> const & flowRate =
    m_perfManager.getReference<array1d<real64>>( viewKeyStruct::flowRateString );

  real64 totalRate = 0.0;
  for (localIndex iconn = 0; iconn < numConnectionsLocal(); ++iconn)
  {
    totalRate += flowRate[iconn];
  }

  return totalRate;
}

// form the well control equation based on the type of well
void SinglePhaseWell::FormControlEquation(  DomainPartition const * const domain,
			                    Epetra_FECrsMatrix * const jacobian,
			                    Epetra_FEVector * const residual,
			                    real64 const time_n,
			                    real64 const dt )
{
  
}

// update each connection pressure from bhp and hydrostatic head
void SinglePhaseWell::StateUpdate( DomainPartition const * domain, localIndex fluidIndex )
{
 
}


REGISTER_CATALOG_ENTRY( WellBase, SinglePhaseWell, string const &, ManagedGroup * const )

} //namespace geosx
