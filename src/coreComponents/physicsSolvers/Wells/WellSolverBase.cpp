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
 * @file WellSolverBase.cpp
 */

#include "WellSolverBase.hpp"

#include "managers/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "physicsSolvers/FiniteVolume/FlowSolverBase.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
  
WellSolverBase::WellSolverBase( std::string const & name,
                                ManagedGroup * const parent )
  : SolverBase( name, parent ),
    m_gravityFlag(1),
    m_fluidName(),
    m_fluidIndex(),
    m_numDofPerElement(0),
    m_numDofPerConnection(0),
    m_numDofPerResElement(0),
    m_firstWellElemDofNumber(-1)
{
  RegisterViewWrapper( viewKeyStruct::gravityFlagString, &m_gravityFlag, false )->
    setApplyDefaultValue(1)->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Flag that enables/disables gravity");

  this->RegisterViewWrapper( viewKeyStruct::fluidNameString,  &m_fluidName,  false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of fluid constitutive object to use for this solver.");

  this->RegisterViewWrapper( viewKeyStruct::fluidIndexString, &m_fluidIndex, false );
}

void WellSolverBase::RegisterDataOnMesh( ManagedGroup * const meshBodies )
{
  SolverBase::RegisterDataOnMesh( meshBodies );

  WellManager * wellManager = meshBodies->getParent()->group_cast<DomainPartition *>()->getWellManager();
  wellManager->RegisterViewWrapper<array1d<real64>>( viewKeyStruct::bhpString );

}

void WellSolverBase::SetSparsityPattern( DomainPartition const * const domain,
                                         Epetra_FECrsGraph * const sparsity,
					 globalIndex firstWellElemDofNumber,
					 localIndex numDofPerResElement )
{
}

void WellSolverBase::SetNumRowsAndTrilinosIndices( DomainPartition const * const domain,
                                                   localIndex & numLocalRows,
                                                   globalIndex & numGlobalRows,
                                                   localIndex offset )
{
}
 
void WellSolverBase::InitializePreSubGroups(ManagedGroup * const rootGroup)
{
  SolverBase::InitializePreSubGroups(rootGroup);

  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);
  ConstitutiveManager * const cm = domain->getConstitutiveManager();

  ConstitutiveBase const * fluid  = cm->GetConstitituveRelation<ConstitutiveBase>( m_fluidName );
  GEOS_ERROR_IF( fluid == nullptr, "Fluid model " + m_fluidName + " not found" );
  m_fluidIndex = fluid->getIndexInParent();

}
  
void WellSolverBase::InitializePostInitialConditions_PreSubGroups(ManagedGroup * const rootGroup)
{
  SolverBase::InitializePostInitialConditions_PreSubGroups(rootGroup);

  std::cout << "WellSolverBase: InitializePostInitialConditions_PreSubGroups" << std::endl;
  
  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);

  // bind the stored reservoir views to the current domain
  ResetViews( domain );

  // Precompute solver-specific constant data (e.g. gravity-depth)
  PrecomputeData(domain);
}

void WellSolverBase::PrecomputeData(DomainPartition * const domain)
{
  // TODO
}

void WellSolverBase::ResetViews( DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = mesh->getElemManager();

  m_resGravDepth =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( FlowSolverBase::viewKeyStruct::gravityDepthString );
}

void WellSolverBase::FormControlEquation( DomainPartition * const domain,
                                          Epetra_FECrsMatrix * const jacobian,
                                          Epetra_FEVector * const residual )
{
  // TODO
}

  
WellSolverBase::~WellSolverBase() = default;

} // namespace geosx
