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
    m_resFluidIndex(),
    m_numDofPerElement(0),
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

  this->RegisterViewWrapper( viewKeyStruct::resFluidIndexString, &m_resFluidIndex, false );

}

void WellSolverBase::RegisterDataOnMesh( ManagedGroup * const meshBodies )
{
  SolverBase::RegisterDataOnMesh( meshBodies );
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

  m_resFluidIndex = fluid->getIndexInParent(); // WARNING: assume same index, not sure it is true
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
}

void WellSolverBase::ResetViews( DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = mesh->getElemManager();

  m_resGravDepth =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( FlowSolverBase::viewKeyStruct::gravityDepthString );
}

  
WellSolverBase::~WellSolverBase() = default;

globalIndex WellSolverBase::getElementOffset( globalIndex welemDofNumber ) const
{
  /*
   * firstWellElemDofNumber denotes the first DOF number of the well segments, for all the wells (i.e., first segment of first well)
   * currentElemDofNumber denotes the DOF number of the current segment
   *
   * The coordinates of this element's 2x2 block in J_WW can be accessed using:
   *
   * IndexRow = firstWellElemDofNumber * resNDOF ( = all the equations in J_RR)
   *          + (currentElemDofNumber - firstWellElemDofNumber ) * wellNDOF ( = offset of current segment in J_WW)
   *          + idof ( = local dofs for this segment, pressure and rate)
   *           
   * This is needed because resNDOF is not equal to wellNDOF
   */

  localIndex const resNDOF  = numDofPerResElement(); // dof is pressure
  localIndex const wellNDOF = numDofPerElement(); // dofs are pressure and rate
  
  globalIndex const firstElemDofNumber = getFirstWellElementDofNumber();
  globalIndex const currentElemOffset  = firstElemDofNumber * resNDOF // number of eqns in J_RR
                                       + (welemDofNumber - firstElemDofNumber) * wellNDOF; // number of eqns in J_WW, before this element's equations

  return currentElemOffset;
}
  
} // namespace geosx
