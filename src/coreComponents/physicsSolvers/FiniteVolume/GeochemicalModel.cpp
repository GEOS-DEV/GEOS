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
 * @file GeochemicalModel.cpp
 */

#include "GeochemicalModel.hpp"

#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/Fluid/ReactiveFluidBase.hpp"
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

/**
 * @namespace the geosx namespace that encapsulates the majority of the code
 */
namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace systemSolverInterface;

GeochemicalModel::GeochemicalModel( const std::string& name,
                                  ManagedGroup * const parent ):
  FlowSolverBase(name, parent)
{

  this->RegisterViewWrapper( viewKeyStruct::reactiveFluidNameString,  &m_reactiveFluidName,  false )->setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of chemical system constitutive object to use for this solver");

  this->RegisterViewWrapper( viewKeyStruct::reactiveFluidIndexString, &m_reactiveFluidIndex, false );

  this->RegisterViewWrapper( viewKeyStruct::outputSpeciesFileNameString, &m_outputSpeciesFileName, false )->
    setInputFlag(InputFlags::OPTIONAL)->    
    setDescription("Output species to file");

}

void GeochemicalModel::RegisterDataOnMesh(ManagedGroup * const MeshBodies)
{
  FlowSolverBase::RegisterDataOnMesh(MeshBodies);

  for( auto & mesh : MeshBodies->GetSubGroups() )
  {

    MeshLevel * const meshLevel = ManagedGroup::group_cast<MeshBody *>(mesh.second)->getMeshLevel(0);

    applyToSubRegions( meshLevel, [&] ( ElementSubRegionBase * const subRegion)
    {
      subRegion->RegisterViewWrapper< array1d<real64> >( viewKeyStruct::pressureString )->setPlotLevel(PlotLevel::LEVEL_0);
      subRegion->RegisterViewWrapper< array1d<real64> >( viewKeyStruct::deltaPressureString );

      subRegion->RegisterViewWrapper< array1d<real64> >( viewKeyStruct::temperatureString )->setPlotLevel(PlotLevel::LEVEL_0);
      subRegion->RegisterViewWrapper< array1d<real64> >( viewKeyStruct::deltaTemperatureString );      

      subRegion->RegisterViewWrapper< array2d<real64> >( viewKeyStruct::concentrationString )->setPlotLevel(PlotLevel::LEVEL_0);
      subRegion->RegisterViewWrapper< array2d<real64> >( viewKeyStruct::deltaConcentrationString );

      subRegion->RegisterViewWrapper< array2d<real64> >( viewKeyStruct::totalConcentrationString );
      subRegion->RegisterViewWrapper< array2d<real64> >( viewKeyStruct::concentrationNewString );                  

      subRegion->RegisterViewWrapper< array1d<globalIndex> >( viewKeyStruct::blockLocalDofNumberString );

    } );

  }
}

void GeochemicalModel::InitializePreSubGroups(ManagedGroup * const rootGroup)
{
  FlowSolverBase::InitializePreSubGroups(rootGroup);

  // set the blockID for the block system interface
  getLinearSystemRepository()->SetBlockID( BlockIDs::geochemicalModelBlock, this->getName() );

  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);
  
  ConstitutiveManager * const cm = domain->getConstitutiveManager();

  ReactiveFluidBase const * reactiveFluid  = cm->GetConstitutiveRelation<ReactiveFluidBase>( m_reactiveFluidName );
  
  GEOS_ERROR_IF( reactiveFluid == nullptr, "Geochemical model " + m_reactiveFluidName + " not found" );
  m_reactiveFluidIndex = reactiveFluid->getIndexInParent();

  m_numBasisSpecies     = reactiveFluid->numBasisSpecies();
  m_numDependentSpecies  = reactiveFluid->numDependentSpecies();  

  m_numDofPerCell = m_numBasisSpecies;

  for( auto & mesh : domain->getMeshBodies()->GetSubGroups() )
  {
    MeshLevel * meshLevel = ManagedGroup::group_cast<MeshBody *>(mesh.second)->getMeshLevel(0);
    ResizeFields( meshLevel );
  }
}

void GeochemicalModel::ResizeFields( MeshLevel * const meshLevel )
{
  localIndex const NC = m_numBasisSpecies;

  applyToSubRegions( meshLevel, [&] ( ElementSubRegionBase * const subRegion )
  {
    subRegion->getReference< array2d<real64> >(viewKeyStruct::concentrationString).resizeDimension<1>(NC);
    subRegion->getReference< array2d<real64> >(viewKeyStruct::deltaConcentrationString).resizeDimension<1>(NC);

    subRegion->getReference< array2d<real64> >(viewKeyStruct::totalConcentrationString).resizeDimension<1>(NC);

    subRegion->getReference< array2d<real64> >(viewKeyStruct::concentrationNewString).resizeDimension<1>(NC);    

  });

}
  
void GeochemicalModel::UpdateReactiveFluidModel(ManagedGroup * const dataGroup)
{
  GEOSX_MARK_FUNCTION;

  ReactiveFluidBase * const reactiveFluid = GetConstitutiveModel<ReactiveFluidBase>( dataGroup, m_reactiveFluidName );

  arrayView1d<real64 const> const & pres = dataGroup->getReference<array1d<real64>>( viewKeyStruct::pressureString );
  arrayView1d<real64 const> const & dPres = dataGroup->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

  arrayView1d<real64 const> const & temp = dataGroup->getReference<array1d<real64>>( viewKeyStruct::temperatureString );
  arrayView1d<real64 const> const & dTemp = dataGroup->getReference<array1d<real64>>( viewKeyStruct::deltaTemperatureString );  

  arrayView2d<real64 const> const & conc = dataGroup->getReference<array2d<real64>>( viewKeyStruct::concentrationString );
  arrayView2d<real64 const> const & dConc = dataGroup->getReference<array2d<real64>>( viewKeyStruct::deltaConcentrationString );  
  arrayView2d<real64> & concNew = dataGroup->getReference<array2d<real64>>( viewKeyStruct::concentrationNewString );
  
  forall_in_range<RAJA::seq_exec>( 0, dataGroup->size(), GEOSX_LAMBDA ( localIndex const a )
  {
    for(localIndex ic = 0; ic < m_numBasisSpecies; ++ic)
      concNew[a][ic] = conc[a][ic] + dConc[a][ic];
    
    reactiveFluid->PointUpdate( pres[a] + dPres[a], temp[a] + dTemp[a], concNew[a], a);
  });

}

void GeochemicalModel::UpdateState( ManagedGroup * dataGroup )
{
  GEOSX_MARK_FUNCTION;

  UpdateReactiveFluidModel( dataGroup );

}

void GeochemicalModel::InitializePostInitialConditions_PreSubGroups( ManagedGroup * const rootGroup )
{
  GEOSX_MARK_FUNCTION;

  FlowSolverBase::InitializePostInitialConditions_PreSubGroups( rootGroup );

  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);
  MeshLevel * mesh = domain->getMeshBody(0)->getMeshLevel(0);

  //TODO this is a hack until the sets are fixed to include ghosts!!
  std::map<string, string_array > fieldNames;
  fieldNames["elems"].push_back( viewKeyStruct::pressureString );
  fieldNames["elems"].push_back( viewKeyStruct::temperatureString );  
  fieldNames["elems"].push_back( viewKeyStruct::concentrationString );  

  array1d<NeighborCommunicator> & comms =
    domain->getReference< array1d<NeighborCommunicator>>( domain->viewKeys.neighbors );

  CommunicationTools::SynchronizeFields( fieldNames, mesh, comms );

  ResetViews( domain );

}

real64 GeochemicalModel::SolverStep( real64 const& time_n,
                                    real64 const& dt,
                                    const int cycleNumber,
                                    DomainPartition * domain )
{
  GEOSX_MARK_FUNCTION;

  //  FlowSolverBase::PrecomputeData(domain);

  MeshLevel * mesh = domain->getMeshBody(0)->getMeshLevel(0);

  NodeManager const * const nodeManager = mesh->getNodeManager();
  FaceManager const * const faceManager = mesh->getFaceManager();
  
  real64 dt_return = dt;

  ImplicitStepSetup( time_n, dt, domain, getLinearSystemRepository() );

  // currently the only method is implicit time integration
  dt_return= this->NonlinearImplicitStep( time_n,
                                          dt,
                                          cycleNumber,
                                          domain,
                                          getLinearSystemRepository() );

  // final step for completion of timestep. typically secondary variable updates and cleanup.
  ImplicitStepComplete( time_n, dt_return, domain );

  return dt_return;

}


void GeochemicalModel::ImplicitStepSetup( real64 const& time_n,
                                         real64 const& dt,
                                         DomainPartition * const domain,
                                         EpetraBlockSystem * const blockSystem )
{
  ResetViews( domain );

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  /* The loop below could be moved to SolverStep after ImplicitStepSetup */
  
  applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegion * const region,
                                 ElementSubRegionBase * const subRegion )
  {

    arrayView1d<real64> const & dPres   = m_deltaPressure[er][esr];
    arrayView1d<real64> const & dTemp   = m_deltaTemperature[er][esr];
    arrayView2d<real64> const & dConc   = m_deltaConcentration[er][esr];
    arrayView2d<real64> const & conc   = m_concentration[er][esr];
    arrayView2d<real64> const & totalConc   = m_totalConcentration[er][esr];
    

    forall_in_range<elemPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      dPres[ei] = 0.0;
      dTemp[ei] = 0.0;
      for (localIndex ic = 0; ic < m_numBasisSpecies; ++ic)
	{
	  dConc[ei][ic] = 0.0;
	  conc[ei][ic] = log10(totalConc[ei][ic]);
	}

    } );

    UpdateState( subRegion );
    
  } );

  // setup dof numbers and linear system
  SetupSystem( domain, blockSystem );

}

void GeochemicalModel::ImplicitStepComplete( real64 const & time_n,
                                            real64 const & dt,
                                            DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegion * const region,
                                 ElementSubRegionBase * const subRegion )
  {
    arrayView1d<real64> const & pres = m_pressure[er][esr];
    arrayView1d<real64> const & temp = m_temperature[er][esr];    

    arrayView2d<real64> const & conc = m_concentration[er][esr];
    
    arrayView1d<real64 const> const & dPres = m_deltaPressure[er][esr];
    arrayView1d<real64 const> const & dTemp = m_deltaTemperature[er][esr];    
    arrayView2d<real64 const> const & dConc = m_deltaConcentration[er][esr];    

    forall_in_range<elemPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {

      pres[ei] += dPres[ei];
      temp[ei] += dTemp[ei];      
      for (localIndex ic = 0; ic < m_numBasisSpecies; ++ic)
        conc[ei][ic] += dConc[ei][ic];

    } );

  } );

  WriteSpeciesToFile(domain);

}

void GeochemicalModel::SetNumRowsAndTrilinosIndices( MeshLevel * const meshLevel,
                                                    localIndex & numLocalRows,
                                                    globalIndex & numGlobalRows,
                                                    localIndex offset )
{
  int numMpiProcesses;
  MPI_Comm_size( MPI_COMM_GEOSX, &numMpiProcesses );

  int thisMpiProcess = 0;
  MPI_Comm_rank( MPI_COMM_GEOSX, &thisMpiProcess );

  localIndex numLocalRowsToSend = numLocalRows;
  array1d<localIndex> gather( numMpiProcesses );

  // communicate the number of local rows to each process
  m_linearSolverWrapper.m_epetraComm.GatherAll( &numLocalRowsToSend,
                                                &gather.front(),
                                                1 );

  GEOS_ERROR_IF( numLocalRows != numLocalRowsToSend, "number of local rows inconsistent" );

  // find the first local row on this partition, and find the number of total global rows.
  localIndex firstLocalRow = 0;
  numGlobalRows = 0;

  for( integer p=0 ; p<numMpiProcesses ; ++p )
  {
    numGlobalRows += gather[p];
    if( p < thisMpiProcess )
    {
      firstLocalRow += gather[p];
    }
  }

  // loop over all elements and set the dof number if the element is not a ghost
  localIndex localCount = 0;

  applyToSubRegions( meshLevel, [&] ( localIndex er, localIndex esr,
                                      ElementRegion * const region,
                                      ElementSubRegionBase * const subRegion )
  {
    arrayView1d<integer const> const & elemGhostRank = m_elemGhostRank[er][esr];
    arrayView1d<globalIndex> const & dofNumber = m_dofNumber[er][esr];

    forall_in_range<RAJA::seq_exec>( 0, subRegion->size(), [&] ( localIndex const a )
    {
      if( elemGhostRank[a] < 0 )
      {
        dofNumber[a] = firstLocalRow + localCount + offset;
        localCount += 1;
      }
      else
      {
        dofNumber[a] = -1;
      }
    } );
  } );

  GEOS_ERROR_IF(localCount != numLocalRows, "Number of DOF assigned does not match numLocalRows" );
}

void GeochemicalModel::SetupSystem ( DomainPartition * const domain,
                                    EpetraBlockSystem * const blockSystem )
{
  // assume that there is only a single MeshLevel for now
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  // for this solver, the dof are on the cell center, and the row corrosponds to an element

  localIndex numLocalRows  = 0;
  globalIndex numGlobalRows = 0;

  // get the number of local elements, and ghost elements...i.e. local rows and ghost rows
  applyToSubRegions( mesh, [&] ( ElementSubRegionBase * const subRegion )
  {
    localIndex subRegionGhosts = subRegion->GetNumberOfGhosts();

    numLocalRows += subRegion->size() - subRegionGhosts;

  } );

  SetNumRowsAndTrilinosIndices( mesh,
                                numLocalRows,
                                numGlobalRows,
                                0 );

  std::map<string, string_array > fieldNames;
  fieldNames["elems"].push_back( viewKeyStruct::blockLocalDofNumberString );

  array1d<NeighborCommunicator> & comms =
    domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors );

  CommunicationTools::SynchronizeFields( fieldNames, mesh, comms );

  // construct row map, and set a pointer to the row map
  Epetra_Map * const
  rowMap = blockSystem->
           SetRowMap( BlockIDs::geochemicalModelBlock,
                      std::make_unique<Epetra_Map>( numGlobalRows * m_numDofPerCell,
                                                    numLocalRows * m_numDofPerCell,
                                                    0,
                                                    m_linearSolverWrapper.m_epetraComm ) );

  // construct sparisty matrix, set a pointer to the sparsity pattern matrix
  Epetra_FECrsGraph * const
  sparsity = blockSystem->SetSparsity( BlockIDs::geochemicalModelBlock,
                                       BlockIDs::geochemicalModelBlock,
                                       std::make_unique<Epetra_FECrsGraph>(Copy,*rowMap,0) );


  // set the sparsity patter
  SetSparsityPattern( domain, sparsity );

  // assemble the global sparsity matrix
  sparsity->GlobalAssemble();
  sparsity->OptimizeStorage();

  // construct system matrix
  blockSystem->SetMatrix( BlockIDs::geochemicalModelBlock,
                          BlockIDs::geochemicalModelBlock,
                          std::make_unique<Epetra_FECrsMatrix>(Copy,*sparsity) );

  // construct solution vector
  blockSystem->SetSolutionVector( BlockIDs::geochemicalModelBlock,
                                  std::make_unique<Epetra_FEVector>(*rowMap) );

  // construct residual vector
  blockSystem->SetResidualVector( BlockIDs::geochemicalModelBlock,
                                  std::make_unique<Epetra_FEVector>(*rowMap) );

}

void GeochemicalModel::SetSparsityPattern( DomainPartition const * const domain,
                                          Epetra_FECrsGraph * const sparsity )
{

  MeshLevel const * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  
  applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegion const * const region,
                                 ElementSubRegionBase const * const subRegion )
  {
    arrayView1d<integer const>     const & elemGhostRank = m_elemGhostRank[er][esr];
    arrayView1d<globalIndex const> const & dofNumber     = m_dofNumber[er][esr];

    forall_in_range<elemPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      if (elemGhostRank[ei] < 0)
      {
	stackArray1d<globalIndex, ReactiveFluidBase::MAX_NUM_SPECIES> dofIndex(m_numDofPerCell);
	
        globalIndex const elemDOF = dofNumber[ei];
	globalIndex const offset = m_numDofPerCell * elemDOF;	

	for (localIndex idof = 0; idof < m_numDofPerCell; ++idof)
	  {
	
	    dofIndex[idof] = offset + idof;

	  }

	sparsity->InsertGlobalIndices(integer_conversion<int>(m_numDofPerCell),
				      dofIndex.data(),
				      integer_conversion<int>(m_numDofPerCell),
				      dofIndex.data());

      }
	
    } );
  } );


}

void GeochemicalModel::AssembleSystem( DomainPartition * const domain,
                                      EpetraBlockSystem * const blockSystem,
                                      real64 const time_n,
                                      real64 const dt )
{
  GEOSX_MARK_FUNCTION;

  Epetra_FECrsMatrix * const jacobian = blockSystem->GetMatrix( BlockIDs::geochemicalModelBlock, BlockIDs::geochemicalModelBlock );
  Epetra_FEVector * const residual = blockSystem->GetResidualVector( BlockIDs::geochemicalModelBlock );

  jacobian->Scale(0.0);
  residual->Scale(0.0);

  AssembleAccumulationTerms( domain, jacobian, residual, time_n, dt );

  if( verboseLevel() >= 2 )
  {
    GEOS_LOG_RANK("After GeochemicalModel::AssembleAccumulationTerms");
    GEOS_LOG_RANK("\nJacobian:\n" << *jacobian);
    GEOS_LOG_RANK("\nResidual:\n" << *residual);
  }

  AssembleFluxTerms( domain, jacobian, residual, time_n, dt );

  jacobian->GlobalAssemble(true);
  residual->GlobalAssemble();

  if( verboseLevel() >= 2 )
  {
    GEOS_LOG_RANK("After GeochemicalModel::AssembleSystem");
    GEOS_LOG_RANK("\nJacobian:\n" << *jacobian);
    GEOS_LOG_RANK("\nResidual:\n" << *residual);
  }

}

void GeochemicalModel::AssembleAccumulationTerms( DomainPartition * const domain,
                                                 Epetra_FECrsMatrix * const jacobian,
                                                 Epetra_FEVector * const residual,
                                                 real64 const time_n,
                                                 real64 const dt )
{
  GEOSX_MARK_FUNCTION;

  ConstitutiveManager * const cm = domain->getConstitutiveManager();

  ReactiveFluidBase const * reactiveFluid  = cm->GetConstitutiveRelation<ReactiveFluidBase>( m_reactiveFluidName );
  
  const array2d<real64> & stochMatrix = reactiveFluid->StochMatrix();
  const array1d<bool> & isHplus = reactiveFluid->IsHplus();   

  MeshLevel const * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  
  applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegion const * const region,
                                 ElementSubRegionBase const * const subRegion )
  {
    arrayView1d<integer const>     const & elemGhostRank = m_elemGhostRank[er][esr];
    arrayView1d<globalIndex const> const & dofNumber     = m_dofNumber[er][esr];

    arrayView2d<real64 const> const & conc          = m_concentration[er][esr];
    arrayView2d<real64 const> const & dConc         = m_deltaConcentration[er][esr];    
    arrayView2d<real64 const> const & totalConc     = m_totalConcentration[er][esr];

    arrayView2d<real64 const> const & dependentConc     = m_dependentConc[er][esr][m_reactiveFluidIndex];
    arrayView3d<real64 const> const & dDependentConc_dConc     = m_dDependentConc_dConc[er][esr][m_reactiveFluidIndex];         

    forall_in_range<elemPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      if (elemGhostRank[ei] < 0)
      {
	stackArray1d<globalIndex, ReactiveFluidBase::MAX_NUM_SPECIES> dofIndex(m_numDofPerCell);
        stackArray1d<real64, ReactiveFluidBase::MAX_NUM_SPECIES>  resid(m_numDofPerCell);
        stackArray2d<real64, ReactiveFluidBase::MAX_NUM_SPECIES * ReactiveFluidBase::MAX_NUM_SPECIES> matrix(m_numDofPerCell, m_numDofPerCell);
	
        globalIndex const elemDOF = dofNumber[ei];
	globalIndex const offset = m_numDofPerCell * elemDOF;	

	for (localIndex idof = 0; idof < m_numDofPerCell; ++idof)
	  {
	
	    dofIndex[idof] = offset + idof;

	  }

	matrix = 0.0;

	for (localIndex ic = 0; ic < m_numBasisSpecies; ++ic)
	  {

	    real64 concBasis = pow(10.0, conc[ei][ic]+dConc[ei][ic]);

	    resid[ic] = concBasis - totalConc[ei][ic];

	    matrix[ic][ic] = log(10.0) * concBasis;

	    if(isHplus[ic])
	      {
		resid[ic] = 0.0;
		matrix[ic][ic] = 1.0;
		continue;
	      }

	    for (localIndex id = 0; id < m_numDependentSpecies; ++id)
	      {
		real64 concDependent = pow(10.0, dependentConc[ei][id]);
		resid[ic] -= stochMatrix[ic][id] * concDependent;

		for (localIndex idc = 0; idc < m_numBasisSpecies; ++idc)
		  {

		    matrix[ic][idc] -= stochMatrix[ic][id] * log(10.0) * concDependent * dDependentConc_dConc[ei][id][idc];

		  }
	      }
	  }

        // add contribution to global residual and jacobian
        residual->SumIntoGlobalValues(integer_conversion<int>(m_numDofPerCell), dofIndex.data(), resid.data());

        jacobian->SumIntoGlobalValues(integer_conversion<int>(m_numDofPerCell), dofIndex.data(), integer_conversion<int>(m_numDofPerCell), dofIndex.data(), matrix.data(),Epetra_FECrsMatrix::ROW_MAJOR);

      }

    } );
  } );

}

void GeochemicalModel::AssembleFluxTerms(DomainPartition * const domain,
                                         Epetra_FECrsMatrix * const jacobian,
                                         Epetra_FEVector * const residual,
                                         real64 const time_n,
                                         real64 const dt )
{

}

void GeochemicalModel::ApplyBoundaryConditions( DomainPartition * const domain,
                                               EpetraBlockSystem * const blockSystem,
                                               real64 const time_n,
                                               real64 const dt )
{
}


real64
GeochemicalModel::
CalculateResidualNorm( EpetraBlockSystem const * const blockSystem,
                       DomainPartition * const domain )
{
  Epetra_FEVector const * const residual = blockSystem->GetResidualVector( BlockIDs::geochemicalModelBlock );
  Epetra_Map      const * const rowMap   = blockSystem->GetRowMap( BlockIDs::geochemicalModelBlock );

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

    localResidualNorm += sum_in_range( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex const a )
    {
      if (elemGhostRank[a] < 0)
      {
        real64 cell_norm = 0.0;
        globalIndex const offset = m_numDofPerCell * dofNumber[a];
        for (localIndex idof = 0; idof < m_numDofPerCell; ++idof)
        {
          int const lid = rowMap->LID(integer_conversion<int>(offset + idof));
	  real64 const val = localResidual[lid];	  
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

void GeochemicalModel::ApplySystemSolution( EpetraBlockSystem const * const blockSystem,
                                           real64 const scalingFactor,
                                           DomainPartition * const domain )
{
  
  Epetra_Map const * const rowMap        = blockSystem->GetRowMap( BlockIDs::geochemicalModelBlock );
  Epetra_FEVector const * const solution = blockSystem->GetSolutionVector( BlockIDs::geochemicalModelBlock );

  MeshLevel * mesh = domain->getMeshBody(0)->getMeshLevel(0);

  int dummy;
  double* local_solution = nullptr;
  solution->ExtractView( &local_solution, &dummy );

  applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegion * const region,
                                 ElementSubRegionBase * const subRegion )
  {
    arrayView1d<globalIndex const> const & dofNumber = m_dofNumber[er][esr];
    arrayView1d<integer const> const & elemGhostRank = m_elemGhostRank[er][esr];

    arrayView2d<real64> const & dConc = m_deltaConcentration[er][esr];

    arrayView2d<real64> const & conc = m_concentration[er][esr];        
    
    forall_in_range<elemPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      if (elemGhostRank[ei] < 0)
      {
	int lid;
	for(localIndex ic = 0; ic < m_numBasisSpecies; ++ic)
	  {
	
	    lid = rowMap->LID( integer_conversion<int>( dofNumber[ei] * m_numDofPerCell + ic) );
	    dConc[ei][ic] += scalingFactor * local_solution[lid];

	  }

      }
    } );
  } );

  std::map<string, string_array > fieldNames;
  fieldNames["elems"].push_back( viewKeyStruct::deltaPressureString );
  fieldNames["elems"].push_back( viewKeyStruct::deltaTemperatureString );  
  fieldNames["elems"].push_back( viewKeyStruct::deltaConcentrationString );

  array1d<NeighborCommunicator> & comms =
    domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors );

  CommunicationTools::SynchronizeFields( fieldNames, mesh, comms );

  applyToSubRegions( mesh, [&] ( ElementSubRegionBase * subRegion )
  {
    UpdateState( subRegion );
  } );

}

void GeochemicalModel::SolveSystem( EpetraBlockSystem * const blockSystem,
                                   SystemSolverParameters const * const params )
{
  GEOSX_MARK_FUNCTION;

  Epetra_FEVector * const
  solution = blockSystem->GetSolutionVector( BlockIDs::geochemicalModelBlock );

  Epetra_FEVector * const
  residual = blockSystem->GetResidualVector( BlockIDs::geochemicalModelBlock );

  residual->Scale(-1.0);

  solution->Scale(0.0);
  
  m_linearSolverWrapper.SolveSingleBlockSystem( blockSystem,
                                                params,
                                                BlockIDs::geochemicalModelBlock );

  if( verboseLevel() >= 2 )
  {
    GEOS_LOG_RANK("After GeochemicalModel::SolveSystem");
    GEOS_LOG_RANK("\nsolution\n" << *solution);
  }

}

void GeochemicalModel::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
}

void GeochemicalModel::ResetViews(DomainPartition * const domain)
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

  m_temperature =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( viewKeyStruct::temperatureString );
  m_deltaTemperature =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( viewKeyStruct::deltaTemperatureString );  

  m_concentration =
    elemManager->ConstructViewAccessor<array2d<real64>, arrayView2d<real64>>( viewKeyStruct::concentrationString );
  m_deltaConcentration =
    elemManager->ConstructViewAccessor<array2d<real64>, arrayView2d<real64>>( viewKeyStruct::deltaConcentrationString );
  m_totalConcentration =
    elemManager->ConstructViewAccessor<array2d<real64>, arrayView2d<real64>>( viewKeyStruct::totalConcentrationString );  

  m_concentrationNew =
    elemManager->ConstructViewAccessor<array2d<real64>, arrayView2d<real64>>( viewKeyStruct::concentrationNewString );

  
  m_dependentConc = 
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64> >( ReactiveFluidBase::viewKeyStruct::dependentConcString, constitutiveManager );

  m_dDependentConc_dConc = 
    elemManager->ConstructFullMaterialViewAccessor<array3d<real64>, arrayView3d<real64> >( ReactiveFluidBase::viewKeyStruct::dDependentConc_dConcString, constitutiveManager );

}

void GeochemicalModel::WriteSpeciesToFile(DomainPartition * const domain)
{

  GEOSX_MARK_FUNCTION;

  if(m_outputSpeciesFileName.empty())
    return;

  std::ofstream os(m_outputSpeciesFileName);
  GEOS_ERROR_IF(!os.is_open(), "Cannot open the species-output file");
  
  ConstitutiveManager * const cm = domain->getConstitutiveManager();

  ReactiveFluidBase const * reactiveFluid  = cm->GetConstitutiveRelation<ReactiveFluidBase>( m_reactiveFluidName );

  const string_array & basisSpeciesNames = reactiveFluid->basisiSpeciesNames();
  const string_array & dependentSpeciesNames = reactiveFluid->dependentSpeciesNames();  
  
  MeshLevel const * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  ElementRegionManager const * const elementManager = mesh->getElemManager();

  for( localIndex er = 0 ; er < elementManager->numRegions() ; ++er )
    {
      
      ElementRegion const * const elemRegion = elementManager->GetRegion(er);
      
      for( localIndex esr = 0 ; esr < elemRegion->numSubRegions() ; ++esr) 
	{

	  arrayView1d<integer const>     const & elemGhostRank = m_elemGhostRank[er][esr];

	  arrayView2d<real64 const> const & conc          = m_concentration[er][esr];

	  arrayView2d<real64 const> const & dependentConc     = m_dependentConc[er][esr][m_reactiveFluidIndex];

	  CellElementSubRegion const * const subRegion = elemRegion->GetSubRegion<CellElementSubRegion>(esr);
      
	  for( localIndex ei = 0 ; ei < subRegion->size() ; ++ei )
	    {

	      if (elemGhostRank[ei] < 0)
		{

		  array1d<localIndex> indices;
		  array1d<real64> speciesConc;
		  localIndex count  = 0;

		  for(localIndex ic = 0; ic < m_numBasisSpecies; ++ic)
		    {

		      indices.push_back(count++);

		      speciesConc.push_back(conc[ei][ic]);

		    }

		  for(localIndex ic = 0; ic < m_numDependentSpecies; ++ic)
		    {
	    
		      indices.push_back(count++);

		      speciesConc.push_back(dependentConc[ei][ic]);	    

		    }

		  std::sort( indices.begin(),indices.end(), [&](localIndex i,localIndex j){return speciesConc[i] > speciesConc[j];});

		  os << "   --- Distribution of Aqueous Solute Species ---" << std::endl;

		  os << std::endl;

		  os << "Species                   Molality            Log Molality" << std::endl;

		  os << std::endl;	
	
		  for(localIndex ic = 0; ic < indices.size(); ic++)
		    {

		      localIndex idx = indices[ic];
		      real64 spC, spLogC;
		      string spName;
	    
		      if(idx < m_numBasisSpecies)
			{
			  spName = basisSpeciesNames[idx];
			  spLogC = conc[ei][idx];
			}
		      else 
			{
			  idx -= m_numBasisSpecies;
			  spName = dependentSpeciesNames[idx];		
			  spLogC = dependentConc[ei][idx];
			}

		      auto found = spName.find("(g)");
		      if(found != std::string::npos)
			continue;

		      spC = pow(10.0, spLogC);
	    
		      if(fabs(spLogC) == 0 || spC < 1e-40)
			continue;
	    
		      os <<  std::left << std::setw(25) << spName << std::setw(10) << std::scientific << std::setprecision(4)<< std::right << spC << std::fixed << std::setw(20) << spLogC << std::endl;

		    }

		  os << std::endl;
		  os << std::endl;
	
		  os << "            --- Gas Fugacities ---" << std::endl;

		  os << std::endl;

		  os << "Gas                      Log Fugacity           Fugacity" << std::endl;

		  os << std::endl;

		  for(localIndex ic = 0; ic < indices.size(); ic++)
		    {

		      localIndex idx = indices[ic];
		      real64 spC, spLogC;
		      string spName;
	    
		      if(idx < m_numBasisSpecies)
			{
			  spName = basisSpeciesNames[idx];
			  spLogC = conc[ei][idx];
			}
		      else 
			{
			  idx -= m_numBasisSpecies;
			  spName = dependentSpeciesNames[idx];		
			  spLogC = dependentConc[ei][idx];
			}

		      auto found = spName.find("(g)");
		      if(found == std::string::npos)
			continue;

		      spC = pow(10.0, spLogC);

		      os <<  std::left << std::setw(20) << spName << std::setw(15) << std::fixed << std::setprecision(5)<< std::right << spLogC << std::scientific << std::setw(23) << spC << std::endl;    	    	    

		    }

		  os << std::endl;
		  os << std::endl;
		  os << std::endl;			

		  os << "           --- Ionic Strength ---" << std::endl;

		  os << std::endl;
		  
		  os <<  "Ionic Strength = " <<  std::fixed << std::setprecision(4) << dependentConc[ei][m_numDependentSpecies] << std::endl;    
	
		}
	    }
	}
    }
	
  os.close();

}  

REGISTER_CATALOG_ENTRY( SolverBase, GeochemicalModel, std::string const &, ManagedGroup * const )
} /* namespace geosx */
