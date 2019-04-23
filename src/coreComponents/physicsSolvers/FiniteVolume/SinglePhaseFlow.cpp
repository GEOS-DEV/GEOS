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
 * @file SinglePhaseFlow.cpp
 */

#include "SinglePhaseFlow.hpp"

#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/Fluid/SingleFluidBase.hpp"
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

SinglePhaseFlow::SinglePhaseFlow( const std::string& name,
                                  ManagedGroup * const parent ):
  FlowSolverBase(name, parent)
{
  m_numDofPerCell = 1;

  // set the blockID for the block system interface
  // To generate the schema, multiple solvers of that use this command are constructed
  // Doing this can cause an error in the block setup, so move it to InitializePreSubGroups
  // getLinearSystemRepository()->SetBlockID( BlockIDs::fluidPressureBlock, this->getName() );
}

void SinglePhaseFlow::RegisterDataOnMesh(ManagedGroup * const MeshBodies)
{
  FlowSolverBase::RegisterDataOnMesh(MeshBodies);

  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    MeshLevel * meshLevel = ManagedGroup::group_cast<MeshBody *>(mesh.second)->getMeshLevel(0);

    ElementRegionManager * const elemManager = meshLevel->getElemManager();
    elemManager->forElementSubRegions([&]( ManagedGroup * const cellBlock) -> void
    {
      cellBlock->RegisterViewWrapper< array1d<real64> >( viewKeyStruct::pressureString )->setPlotLevel(PlotLevel::LEVEL_0);
      cellBlock->RegisterViewWrapper< array1d<real64> >( viewKeyStruct::deltaPressureString );
      cellBlock->RegisterViewWrapper< array1d<real64> >( viewKeyStruct::deltaVolumeString );
      cellBlock->RegisterViewWrapper< array1d<real64> >( viewKeyStruct::densityString );
      cellBlock->RegisterViewWrapper< array1d<real64> >( viewKeyStruct::porosityString )->setPlotLevel(PlotLevel::LEVEL_1);
      cellBlock->RegisterViewWrapper< array1d<real64> >( viewKeyStruct::oldPorosityString );
      cellBlock->RegisterViewWrapper< array1d<globalIndex> >( viewKeyStruct::blockLocalDofNumberString );
    });

    FaceManager * const faceManager = meshLevel->getFaceManager();
    {
      faceManager->RegisterViewWrapper<array1d<real64> >( viewKeyStruct::facePressureString );
      faceManager->RegisterViewWrapper<array1d<real64> >( viewKeyStruct::densityString );
      faceManager->RegisterViewWrapper<array1d<real64> >( viewKeyStruct::viscosityString );
    }
  }
}

void SinglePhaseFlow::InitializePreSubGroups(ManagedGroup * const rootGroup)
{
  FlowSolverBase::InitializePreSubGroups(rootGroup);

  // set the blockID for the block system interface
  getLinearSystemRepository()->SetBlockID( BlockIDs::fluidPressureBlock, this->getName() );
}

void SinglePhaseFlow::UpdateFluidModel(ManagedGroup * const dataGroup)
{
  GEOSX_MARK_FUNCTION;

  SingleFluidBase * const fluid = GetConstitutiveModel<SingleFluidBase>( dataGroup, m_fluidName );

  arrayView1d<real64 const> const & pres = dataGroup->getReference<array1d<real64>>( viewKeyStruct::pressureString );
  arrayView1d<real64 const> const & dPres = dataGroup->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

  // TODO replace with batch update (need up-to-date pressure and temperature fields)
  forall_in_range<RAJA::seq_exec>( 0, dataGroup->size(), GEOSX_LAMBDA ( localIndex const a )
  {
    fluid->PointUpdate( pres[a] + dPres[a], a, 0 );
  });
  //fluid->BatchUpdate( pres, temp, compFrac );
}

void SinglePhaseFlow::UpdateSolidModel(ManagedGroup * const dataGroup)
{
  GEOSX_MARK_FUNCTION;

  ConstitutiveBase * const solid = GetConstitutiveModel<ConstitutiveBase>( dataGroup, m_solidName );

  arrayView1d<real64 const> const & pres  = dataGroup->getReference<array1d<real64>>( viewKeyStruct::pressureString );
  arrayView1d<real64 const> const & dPres = dataGroup->getReference<array1d<real64>>( viewKeyStruct::deltaPressureString );

  forall_in_range( 0, dataGroup->size(), GEOSX_LAMBDA ( localIndex const a )
  {
    solid->StateUpdatePointPressure( pres[a] + dPres[a], a, 0 );
  });
}


void SinglePhaseFlow::UpdateState( ManagedGroup * dataGroup )
{
  GEOSX_MARK_FUNCTION;

  UpdateFluidModel( dataGroup );
  UpdateSolidModel( dataGroup );
}

void SinglePhaseFlow::InitializePostInitialConditions_PreSubGroups( ManagedGroup * const rootGroup )
{
  GEOSX_MARK_FUNCTION;

  FlowSolverBase::InitializePostInitialConditions_PreSubGroups( rootGroup );

  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);

  //TODO this is a hack until the sets are fixed to include ghosts!!
  std::map<string, string_array > fieldNames;
  fieldNames["elems"].push_back( viewKeyStruct::pressureString );
  CommunicationTools::SynchronizeFields(fieldNames,
                              domain->getMeshBody(0)->getMeshLevel(0),
                              domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );

  ResetViews( domain );

  applyToSubRegions( domain, [&] ( ElementSubRegionBase * subRegion )
  {
    UpdateState( subRegion );
  } );

  // Moved the following part from ImplicitStepSetup to here since it only needs to be initialized once
  // They will be updated in ApplySystemSolution and ImplicitStepComplete, respectively
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  elemManager->forElementSubRegionsComplete( [&] ( localIndex er, localIndex esr,
                                            ElementRegion * const region,
                                            ElementSubRegionBase const * const subRegion )
  {
    arrayView1d<real64 const> const & poroRef = m_porosityRef[er][esr];
    arrayView2d<real64 const> const & dens    = m_density[er][esr][m_fluidIndex];
    arrayView2d<real64 const> const & pvmult  = m_pvMult[er][esr][m_solidIndex];

    arrayView1d<real64> const & poro = m_porosity[er][esr];
    arrayView1d<real64> const & densOld = m_densityOld[er][esr];
    arrayView1d<real64> const & poroOld = m_porosityOld[er][esr];

    if( pvmult.size() == poro.size() )
    {
      forall_in_range<elemPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
      {
        densOld[ei] = dens[ei][0];
        poro[ei] = poroRef[ei] * pvmult[ei][0];
        poroOld[ei] = poro[ei];
      } );
    }
    else
    {
      forall_in_range<elemPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
      {
        densOld[ei] = dens[ei][0];
        poro[ei] = poroRef[ei];
        poroOld[ei] = poro[ei];
      } );
    }
  } );
}

real64 SinglePhaseFlow::SolverStep( real64 const& time_n,
                                    real64 const& dt,
                                    const int cycleNumber,
                                    DomainPartition * domain )
{
  GEOSX_MARK_FUNCTION;

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


void SinglePhaseFlow::ImplicitStepSetup( real64 const& time_n,
                                         real64 const& dt,
                                         DomainPartition * const domain,
                                         EpetraBlockSystem * const blockSystem )
{
  ResetViews( domain );

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  elemManager->forElementSubRegionsComplete( [&] ( localIndex er, localIndex esr,
                                            ElementRegion * const region,
                                            ElementSubRegionBase * const subRegion )
  {
    arrayView2d<real64 const> const & dens = m_density[er][esr][m_fluidIndex];
    arrayView1d<real64 const> const & poro = m_porosity[er][esr];

    arrayView1d<real64> const & dPres   = m_deltaPressure[er][esr];
    arrayView1d<real64> const & dVol    = m_deltaVolume[er][esr];
    arrayView1d<real64> const & densOld = m_densityOld[er][esr];
    arrayView1d<real64> const & poroOld = m_porosityOld[er][esr];

    forall_in_range<elemPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      dPres[ei] = 0.0;
      dVol[ei] = 0.0;
      densOld[ei] = dens[ei][0];
      poroOld[ei] = poro[ei];
    } );
  } );

  // setup dof numbers and linear system
  if (!m_reservoirWellsSystemFlag)
    SetupSystem( domain, blockSystem );
}

void SinglePhaseFlow::ImplicitStepComplete( real64 const & time_n,
                                            real64 const & dt,
                                            DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  elemManager->forElementSubRegionsComplete( [&] ( localIndex er, localIndex esr,
                                            ElementRegion * const region,
                                            ElementSubRegionBase * const subRegion )
  {
    arrayView1d<real64> const & pres = m_pressure[er][esr];
    arrayView1d<real64> const & vol  = m_volume[er][esr];

    arrayView1d<real64 const> const & dPres = m_deltaPressure[er][esr];
    arrayView1d<real64 const> const & dVol  = m_deltaVolume[er][esr];

    forall_in_range<elemPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      pres[ei] += dPres[ei];
      vol[ei] += dVol[ei];
    } );
  } );
}



void SinglePhaseFlow::SetNumRowsAndTrilinosIndices( MeshLevel * const meshLevel,
                                                    localIndex & numLocalRows,
                                                    globalIndex & numGlobalRows,
                                                    localIndex offset )
{
  ElementRegionManager::ElementViewAccessor<arrayView1d<globalIndex>> const & dofNumber = m_dofNumber;
  ElementRegionManager::ElementViewAccessor<arrayView1d<integer>> const & elemGhostRank = m_elemGhostRank;

  int numMpiProcesses;
  MPI_Comm_size( MPI_COMM_GEOSX, &numMpiProcesses );

  int thisMpiProcess = 0;
  MPI_Comm_rank( MPI_COMM_GEOSX, &thisMpiProcess );

  localIndex numLocalRowsToSend = numLocalRows;
  array1d<localIndex> gather(numMpiProcesses);

  // communicate the number of local rows to each process
  m_linearSolverWrapper.m_epetraComm.GatherAll( &numLocalRowsToSend,
                                                &gather.front(),
                                                1 );

  GEOS_ERROR_IF( numLocalRows != numLocalRowsToSend, "number of local rows inconsistent" );

  // find the first local row on this partition, and find the number of total global rows.
  localIndex firstLocalRow = 0;
  numGlobalRows = 0;

  for( integer p=0 ; p<numMpiProcesses ; ++p)
  {
    numGlobalRows += gather[p];
    if(p<thisMpiProcess)
      firstLocalRow += gather[p];
  }

  // create trilinos dof indexing, setting initial values to -1 to indicate unset values.
  for( localIndex er=0 ; er<elemGhostRank.size() ; ++er )
  {
    for( localIndex esr=0 ; esr<elemGhostRank[er].size() ; ++esr )
    {
      dofNumber[er][esr] = -1;
    }
  }

  // loop over all elements and set the dof number if the element is not a ghost
  ReduceSum< reducePolicy, localIndex  > localCount(0);
  forAllElemsInMesh<RAJA::seq_exec>( meshLevel, GEOSX_LAMBDA ( localIndex const er,
                                                               localIndex const esr,
                                                               localIndex const k )
  {
    if( elemGhostRank[er][esr][k] < 0 )
    {
      dofNumber[er][esr][k] = firstLocalRow+localCount+offset;
      localCount += 1;
    }
    else
    {
      dofNumber[er][esr][k] = -1;
    }
  });

  GEOS_ERROR_IF(localCount != numLocalRows, "Number of DOF assigned does not match numLocalRows" );
}

void SinglePhaseFlow::SetupSystem ( DomainPartition * const domain,
                                    EpetraBlockSystem * const blockSystem )
{
  // assume that there is only a single MeshLevel for now
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elementRegionManager = mesh->getElemManager();

  // for this solver, the dof are on the cell center, and the row corrosponds to an element
  localIndex numGhostRows  = 0;
  localIndex numLocalRows  = 0;
  globalIndex numGlobalRows = 0;

  // get the number of local elements, and ghost elements...i.e. local rows and ghost rows
  elementRegionManager->forElementSubRegions( [&]( ObjectManagerBase * const subRegion )
    {
      localIndex subRegionGhosts = subRegion->GetNumberOfGhosts();
      numGhostRows += subRegionGhosts;
      numLocalRows += subRegion->size() - subRegionGhosts;
    });

  SetNumRowsAndTrilinosIndices( mesh,
                                numLocalRows,
                                numGlobalRows,
                                0 );

  //TODO element sync doesn't work yet
  std::map<string, string_array > fieldNames;
  fieldNames["elems"].push_back(viewKeysSinglePhaseFlow.blockLocalDofNumber.Key());
  CommunicationTools::
  SynchronizeFields(fieldNames,
                    mesh,
                    domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );


  // construct row map, and set a pointer to the row map
  Epetra_Map * const
  rowMap = blockSystem->
           SetRowMap( BlockIDs::fluidPressureBlock,
                      std::make_unique<Epetra_Map>( numGlobalRows,
                                                    numLocalRows,
                                                    0,
                                                    m_linearSolverWrapper.m_epetraComm ) );

  // construct sparisty matrix, set a pointer to the sparsity pattern matrix
  Epetra_FECrsGraph * const
  sparsity = blockSystem->SetSparsity( BlockIDs::fluidPressureBlock,
                                       BlockIDs::fluidPressureBlock,
                                       std::make_unique<Epetra_FECrsGraph>(Copy,*rowMap,0) );



  // set the sparsity patter
  SetSparsityPattern( domain, sparsity );

  // assemble the global sparsity matrix
  sparsity->GlobalAssemble();
  sparsity->OptimizeStorage();

  // construct system matrix
  blockSystem->SetMatrix( BlockIDs::fluidPressureBlock,
                          BlockIDs::fluidPressureBlock,
                          std::make_unique<Epetra_FECrsMatrix>(Copy,*sparsity) );

  // construct solution vector
  blockSystem->SetSolutionVector( BlockIDs::fluidPressureBlock,
                                  std::make_unique<Epetra_FEVector>(*rowMap) );

  // construct residual vector
  blockSystem->SetResidualVector( BlockIDs::fluidPressureBlock,
                                  std::make_unique<Epetra_FEVector>(*rowMap) );

}

void SinglePhaseFlow::SetSparsityPattern( DomainPartition const * const domain,
                                               Epetra_FECrsGraph * const sparsity )
{
  MeshLevel const * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  ElementRegionManager::ElementViewAccessor<arrayView1d<globalIndex>> const & dofNumber = m_dofNumber;
  ElementRegionManager::ElementViewAccessor<arrayView1d<integer>> const & elemGhostRank = m_elemGhostRank;

  NumericalMethodsManager const * numericalMethodManager = domain->
    getParent()->GetGroup<NumericalMethodsManager>(keys::numericalMethodsManager);

  FiniteVolumeManager const * fvManager = numericalMethodManager->
    GetGroup<FiniteVolumeManager>(keys::finiteVolumeManager);

  FluxApproximationBase const * fluxApprox = fvManager->getFluxApproximation(m_discretizationName);

  //**** loop over all faces. Fill in sparsity for all pairs of DOF/elem that are connected by face
  constexpr localIndex numElems = 2;
  globalIndex_array elementLocalDofIndexRow(numElems);
  globalIndex_array elementLocalDofIndexCol(numElems);
  fluxApprox->forCellStencils( [&]( FluxApproximationBase::CellStencil const & stencilCollection )
  {
    stencilCollection.forAll( [&] ( StencilCollection<CellDescriptor, real64>::Accessor stencil )
    {
      stencil.forConnected([&] (CellDescriptor const & cell, localIndex const i) -> void
      {
        elementLocalDofIndexRow[i] = dofNumber[cell.region][cell.subRegion][cell.index];
      });

      localIndex const stencilSize = stencil.size();
      elementLocalDofIndexCol.resize(stencilSize);
      stencil.forAll([&] (CellDescriptor const & cell, real64 w, localIndex const i) -> void
      {
        elementLocalDofIndexCol[i] = dofNumber[cell.region][cell.subRegion][cell.index];
      });

      sparsity->InsertGlobalIndices(integer_conversion<int>(numElems),
                                    elementLocalDofIndexRow.data(),
                                    integer_conversion<int>(stencilSize),
                                    elementLocalDofIndexCol.data());
    });
  });


  // loop over all elements and add all locals just in case the above connector loop missed some
  forAllElemsInMesh<RAJA::seq_exec>( meshLevel, GEOSX_LAMBDA ( localIndex const er,
                                                               localIndex const esr,
                                                               localIndex const k )
  {
    if (elemGhostRank[er][esr][k] < 0)
    {
      elementLocalDofIndexRow[0] = dofNumber[er][esr][k];

      sparsity->InsertGlobalIndices( 1,
                                     elementLocalDofIndexRow.data(),
                                     1,
                                     elementLocalDofIndexRow.data());
    }
  });

  // add additional connectivity resulting from boundary stencils
  fluxApprox->forBoundaryStencils( [&] (FluxApproximationBase::BoundaryStencil const & boundaryStencilCollection)
  {
    boundaryStencilCollection.forAll( GEOSX_LAMBDA ( StencilCollection<PointDescriptor, real64>::Accessor stencil )
    {
      elementLocalDofIndexRow.resize(1);
      stencil.forConnected( [&] ( PointDescriptor const & point, localIndex const i )
      {
        if (point.tag == PointDescriptor::Tag::CELL)
        {
          CellDescriptor const & c = point.cellIndex;
          elementLocalDofIndexRow[0] = dofNumber[c.region][c.subRegion][c.index];
        }
      });

      localIndex const stencilSize = stencil.size();
      elementLocalDofIndexCol.resize(stencilSize);
      integer counter = 0;
      stencil.forAll( [&] ( PointDescriptor const & point, real64 w, localIndex i )
      {
        if (point.tag == PointDescriptor::Tag::CELL)
        {
          CellDescriptor const & c = point.cellIndex;
          elementLocalDofIndexCol[counter++] = dofNumber[c.region][c.subRegion][c.index];
        }
      });

      sparsity->InsertGlobalIndices(1,
                                    elementLocalDofIndexRow.data(),
                                    integer_conversion<int>(counter),
                                    elementLocalDofIndexCol.data());
    });
  });
}

void SinglePhaseFlow::AssembleSystem( DomainPartition * const domain,
                                      EpetraBlockSystem * const blockSystem,
                                      real64 const time_n,
                                      real64 const dt )
{
  GEOSX_MARK_FUNCTION;

  Epetra_FECrsMatrix * const jacobian = blockSystem->GetMatrix( BlockIDs::fluidPressureBlock, BlockIDs::fluidPressureBlock );
  Epetra_FEVector * const residual = blockSystem->GetResidualVector( BlockIDs::fluidPressureBlock );

  if (!m_reservoirWellsSystemFlag)
  {
    jacobian->Scale(0.0);
    residual->Scale(0.0);
  }
  
  if (m_poroElasticFlag)
  {
    AssembleAccumulationTerms<true>( domain, jacobian, residual, time_n, dt );
  }
  else
  {
    AssembleAccumulationTerms<false>( domain, jacobian, residual, time_n, dt );
  }

  AssembleFluxTerms( domain, jacobian, residual, time_n, dt );

  if (!m_reservoirWellsSystemFlag)
  {
    jacobian->GlobalAssemble(true);
    residual->GlobalAssemble();
  }

  if( verboseLevel() >= 2 )
  {
    GEOS_LOG_RANK("After SinglePhaseFlow::AssembleSystem");
    GEOS_LOG_RANK("\nJacobian:\n" << *jacobian);
    GEOS_LOG_RANK("\nResidual:\n" << *residual);
  }
}

template< bool ISPORO >
struct AssembleAccumulationTermsHelper;

template<>
struct AssembleAccumulationTermsHelper<true>
{
  inline static constexpr void porosityUpdate( real64 & poro,
                                               real64 & dPoro_dPres,
                                               real64 const biotCoefficient,
                                               real64 const poroOld,
                                               real64 const bulkModulus,
                                               real64 const totalMeanStress,
                                               real64 const oldTotalMeanStress,
                                               real64 const dPres,
                                               real64 const poroRef,
                                               real64 const pvmult,
                                               real64 const dPVMult_dPres )
  {
    dPoro_dPres = (biotCoefficient - poroOld) / bulkModulus;
    poro = poroOld + dPoro_dPres * (totalMeanStress - oldTotalMeanStress + dPres);
  }
};

template<>
struct AssembleAccumulationTermsHelper<false>
{
  inline static constexpr void porosityUpdate( real64 & poro,
                                               real64 & dPoro_dPres,
                                               real64 const biotCoefficient,
                                               real64 const poroOld,
                                               real64 const bulkModulus,
                                               real64 const totalMeanStress,
                                               real64 const oldTotalMeanStress,
                                               real64 const dPres,
                                               real64 const poroRef,
                                               real64 const pvmult,
                                               real64 const dPVMult_dPres )
  {
    poro = poroRef * pvmult;
    dPoro_dPres = dPVMult_dPres * poroRef;
  }
};


template< bool ISPORO >
void SinglePhaseFlow::AssembleAccumulationTerms( DomainPartition const * const domain,
                                                 Epetra_FECrsMatrix * const jacobian,
                                                 Epetra_FEVector * const residual,
                                                 real64 const time_n,
                                                 real64 const dt )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager const * const elemManager = mesh->getElemManager();

  elemManager->forElementSubRegionsComplete( [&] ( localIndex er, localIndex esr,
                                                   ElementRegion const * const region,
                                                   ElementSubRegionBase const * const subRegion )
  {
    arrayView1d<integer const>     const & elemGhostRank = m_elemGhostRank[er][esr];
    arrayView1d<globalIndex const> const & dofNumber     = m_dofNumber[er][esr];

    arrayView1d<real64 const> const & densOld       = m_densityOld[er][esr];
    arrayView1d<real64>       const & poro          = m_porosity[er][esr];
    arrayView1d<real64 const> const & poroOld       = m_porosityOld[er][esr];
    arrayView1d<real64 const> const & poroRef       = m_porosityRef[er][esr];
    arrayView1d<real64 const> const & volume        = m_volume[er][esr];
    arrayView1d<real64 const> const & dVol          = m_deltaVolume[er][esr];
    arrayView2d<real64 const> const & dens          = m_density[er][esr][m_fluidIndex];
    arrayView2d<real64 const> const & dDens_dPres   = m_dDens_dPres[er][esr][m_fluidIndex];
    arrayView2d<real64 const> const & pvmult        = m_pvMult[er][esr][m_solidIndex];
    arrayView2d<real64 const> const & dPVMult_dPres = m_dPvMult_dPres[er][esr][m_solidIndex];

    arrayView1d<real64 const> const & dPres              = m_poroElasticFlag ? m_deltaPressure[er][esr]             : poroOld;
    arrayView1d<real64 const> const & oldTotalMeanStress = m_poroElasticFlag ? m_totalMeanStressOld[er][esr]        : poroOld;
    arrayView1d<real64 const> const & totalMeanStress    = m_poroElasticFlag ? m_totalMeanStress[er][esr]           : poroOld;
    arrayView2d<real64 const> const & bulkModulus        = m_poroElasticFlag ? m_bulkModulus[er][esr][m_solidIndex] : dPVMult_dPres;
    real64 const & biotCoefficient                       = m_poroElasticFlag ? m_biotCoefficient[er][esr][m_solidIndex] : 0;

    forall_in_range<elemPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      if (elemGhostRank[ei] < 0)
      {
        globalIndex const elemDOF = dofNumber[ei];
        real64 const densNew = dens[ei][0];
        real64 const volNew = volume[ei] + dVol[ei];

        real64 dPoro_dPres;
        AssembleAccumulationTermsHelper<ISPORO>::porosityUpdate( poro[ei], dPoro_dPres,
                                                                 biotCoefficient,
                                                                 poroOld[ei],
                                                                 bulkModulus[ei][0],
                                                                 totalMeanStress[ei],
                                                                 oldTotalMeanStress[ei],
                                                                 dPres[ei],
                                                                 poroRef[ei],
                                                                 pvmult[ei][0],
                                                                 dPVMult_dPres[ei][0] );


        // Residual contribution is mass conservation in the cell
        real64 const localAccum = poro[ei]    * densNew     * volNew
                                - poroOld[ei] * densOld[ei] * volume[ei];

        // Derivative of residual wrt to pressure in the cell
        real64 const localAccumJacobian = dPoro_dPres * densNew * volNew
                                        + dDens_dPres[ei][0] * poro[ei] * volNew;

        // add contribution to global residual and jacobian
        residual->SumIntoGlobalValues( 1, &elemDOF, &localAccum );
        jacobian->SumIntoGlobalValues( 1, &elemDOF, 1, &elemDOF, &localAccumJacobian );

      }
    } );
  } );
}


void SinglePhaseFlow::AssembleFluxTerms( DomainPartition const * const domain,
                                         Epetra_FECrsMatrix * const jacobian,
                                         Epetra_FEVector * const residual,
                                         real64 const time_n,
                                         real64 const dt )
{
  GEOSX_MARK_FUNCTION;

  NumericalMethodsManager const * numericalMethodManager =
    domain->getParent()->GetGroup<NumericalMethodsManager>( keys::numericalMethodsManager );

  FiniteVolumeManager const * fvManager =
    numericalMethodManager->GetGroup<FiniteVolumeManager>( keys::finiteVolumeManager );

  FluxApproximationBase const * fluxApprox = fvManager->getFluxApproximation( m_discretizationName );

  ElementRegionManager::ElementViewAccessor<arrayView1d<globalIndex>> const & dofNumber = m_dofNumber;

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  const & pres        = m_pressure;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  const & dPres       = m_deltaPressure;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  const & gravDepth   = m_gravDepth;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & dens        = m_density;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & dDens_dPres = m_dDens_dPres;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & visc        = m_viscosity;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & dVisc_dPres = m_dVisc_dPres;

  constexpr localIndex numElems = 2;
  constexpr localIndex maxStencilSize = StencilCollection<CellDescriptor, real64>::MAX_STENCIL_SIZE;

  fluxApprox->forCellStencils( [&]( FluxApproximationBase::CellStencil const & stencilCollection )
  {
    stencilCollection.forAll( [&] ( StencilCollection<CellDescriptor, real64>::Accessor stencil )
    {
      localIndex const stencilSize = stencil.size();

      // working arrays
      stackArray1d<globalIndex, numElems> eqnRowIndices(numElems);
      stackArray1d<globalIndex, maxStencilSize> dofColIndices(stencilSize);

      stackArray1d<real64, numElems> localFlux(numElems);
      stackArray2d<real64, numElems*maxStencilSize> localFluxJacobian(numElems, stencilSize);

      stackArray1d<real64, numElems> densWeight(numElems);
      stackArray1d<real64, numElems> mobility(numElems);
      stackArray1d<real64, numElems> dMobility_dP(numElems);
      stackArray1d<real64, maxStencilSize> dDensMean_dP(stencilSize);
      stackArray1d<real64, maxStencilSize> dFlux_dP(stencilSize);

      // clear working arrays
      eqnRowIndices = -1;
      dDensMean_dP = 0.0;

      // density averaging weights
      densWeight = 0.5;

      // calculate quantities on primary connected cells
      real64 densMean = 0.0;
      stencil.forConnected( [&] ( CellDescriptor const & cell, localIndex i )
      {
        localIndex const er  = cell.region;
        localIndex const esr = cell.subRegion;
        localIndex const ei  = cell.index;

        eqnRowIndices[i] = dofNumber[er][esr][ei];

        // density
        real64 const density   = dens[er][esr][m_fluidIndex][ei][0];
        real64 const dDens_dP  = dDens_dPres[er][esr][m_fluidIndex][ei][0];

        // viscosity
        real64 const viscosity = visc[er][esr][m_fluidIndex][ei][0];
        real64 const dVisc_dP  = dVisc_dPres[er][esr][m_fluidIndex][ei][0];

        // mobility
        mobility[i]  = density / viscosity;
        dMobility_dP[i]  = dDens_dP / viscosity - mobility[i] / viscosity * dVisc_dP;

        // average density
        densMean += densWeight[i] * density;
        dDensMean_dP[i] = densWeight[i] * dDens_dP;
      });

      //***** calculation of flux *****

      // compute potential difference MPFA-style
      real64 potDif = 0.0;
      stencil.forAll( [&] ( CellDescriptor const & cell, real64 w, localIndex i )
      {
        localIndex const er  = cell.region;
        localIndex const esr = cell.subRegion;
        localIndex const ei  = cell.index;

        dofColIndices[i] = dofNumber[er][esr][ei];

        real64 const gravD    = gravDepth[er][esr][ei];
        real64 const gravTerm = m_gravityFlag ? densMean * gravD : 0.0;
        real64 const dGrav_dP = m_gravityFlag ? dDensMean_dP[i] * gravD : 0.0;

        potDif += w * (pres[er][esr][ei] + dPres[er][esr][ei] - gravTerm);
        dFlux_dP[i] = w * (1.0 - dGrav_dP);
      });

      // upwinding of fluid properties (make this an option?)
      localIndex const k_up = (potDif >= 0) ? 0 : 1;

      // compute the final flux and derivatives
      real64 const flux = mobility[k_up] * potDif;
      for (localIndex ke = 0; ke < stencilSize; ++ke)
      {
        dFlux_dP[ke] *= mobility[k_up];
      }

      dFlux_dP[k_up] += dMobility_dP[k_up] * potDif;

      //***** end flux terms *****

      // populate local flux vector and derivatives
      localFlux[0] =  dt * flux;
      localFlux[1] = -localFlux[0];

      for (localIndex ke = 0; ke < stencilSize; ++ke)
      {
        localFluxJacobian[0][ke] =   dt * dFlux_dP[ke];
        localFluxJacobian[1][ke] = - dt * dFlux_dP[ke];
      }

      // Add to global residual/jacobian
      jacobian->SumIntoGlobalValues( integer_conversion<int>(numElems),
                                     eqnRowIndices.data(),
                                     integer_conversion<int>(stencilSize),
                                     dofColIndices.data(),
                                     localFluxJacobian.data() );

      residual->SumIntoGlobalValues( integer_conversion<int>(numElems), eqnRowIndices.data(), localFlux.data() );
    });
  });
}

void SinglePhaseFlow::ApplyBoundaryConditions( DomainPartition * const domain,
                                               EpetraBlockSystem * const blockSystem,
                                               real64 const time_n,
                                               real64 const dt )
{
  GEOSX_MARK_FUNCTION;

  // apply pressure boundary conditions.
  ApplyDirichletBC_implicit(domain, time_n, dt, blockSystem);
  ApplyFaceDirichletBC_implicit(domain, time_n, dt, blockSystem);

  if (verboseLevel() >= 2)
  {

    Epetra_FECrsMatrix * const jacobian = blockSystem->GetMatrix( BlockIDs::fluidPressureBlock,
                                                                  BlockIDs::fluidPressureBlock );
    Epetra_FEVector * const residual = blockSystem->GetResidualVector( BlockIDs::fluidPressureBlock );

    GEOS_LOG_RANK( "After SinglePhaseFlow::ApplyBoundaryConditions" );
    GEOS_LOG_RANK( "\nJacobian\n" << *jacobian );
    GEOS_LOG_RANK( "\nResidual\n" << *residual );
  }

}

void SinglePhaseFlow::ApplyDirichletBC_implicit( DomainPartition * domain,
                                                 real64 const time_n, real64 const dt,
                                                 EpetraBlockSystem * const blockSystem )
{
  FieldSpecificationManager * fsManager = FieldSpecificationManager::get();

  // call the BoundaryConditionManager::ApplyField function that will check to see
  // if the boundary condition should be applied to this subregion
  fsManager->Apply( time_n + dt, domain, "ElementRegions", viewKeyStruct::pressureString,
                    [&]( FieldSpecificationBase const * const fs,
                    string const &,
                    set<localIndex> const & lset,
                    ManagedGroup * subRegion,
                    string const & ) -> void
  {
    arrayView1d<globalIndex const> const &
    dofNumber = subRegion->getReference< array1d<globalIndex> >( viewKeyStruct::blockLocalDofNumberString );

    arrayView1d<real64 const> const &
    pres = subRegion->getReference<array1d<real64> >( viewKeyStruct::pressureString );

    arrayView1d<real64 const> const &
    dPres = subRegion->getReference<array1d<real64> >( viewKeyStruct::deltaPressureString );

    // call the application of the boundary condition to alter the matrix and rhs
    fs->ApplyBoundaryConditionToSystem<FieldSpecificationEqual>( lset,
                                                                 time_n + dt,
                                                                 subRegion,
                                                                 dofNumber,
                                                                 1,
                                                                 blockSystem,
                                                                 BlockIDs::fluidPressureBlock,
    [&] (localIndex const a) -> real64
    {
      return pres[a] + dPres[a];
    });
  } );
}

void SinglePhaseFlow::ApplyFaceDirichletBC_implicit(DomainPartition * domain,
                                                    real64 const time_n, real64 const dt,
                                                    EpetraBlockSystem * const blockSystem)
{
  FieldSpecificationManager * fsManager = FieldSpecificationManager::get();
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();
  FaceManager * const faceManager = mesh->getFaceManager();

  arrayView2d<localIndex> const & elemRegionList     = faceManager->elementRegionList();
  arrayView2d<localIndex> const & elemSubRegionList  = faceManager->elementSubRegionList();

  arrayView1d<integer> const & faceGhostRank =
    faceManager->getReference<array1d<integer>>(ObjectManagerBase::viewKeyStruct::ghostRankString);

  ConstitutiveManager * const constitutiveManager =
    domain->GetGroup<ConstitutiveManager>(keys::ConstitutiveManager);

  NumericalMethodsManager * const numericalMethodManager =
    domain->getParent()->GetGroup<NumericalMethodsManager>(keys::numericalMethodsManager);

  FiniteVolumeManager * const fvManager = numericalMethodManager->GetGroup<FiniteVolumeManager>(keys::finiteVolumeManager);

  FluxApproximationBase const * const fluxApprox = fvManager->getFluxApproximation( m_discretizationName );

  Epetra_FECrsMatrix * const jacobian = blockSystem->GetMatrix(BlockIDs::fluidPressureBlock, BlockIDs::fluidPressureBlock);
  Epetra_FEVector * const residual = blockSystem->GetResidualVector(BlockIDs::fluidPressureBlock);

  ElementRegionManager::ElementViewAccessor<arrayView1d<globalIndex>> const & dofNumber = m_dofNumber;

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  const & pres        = m_pressure;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  const & dPres       = m_deltaPressure;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>>  const & gravDepth   = m_gravDepth;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & dens        = m_density;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & dDens_dPres = m_dDens_dPres;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & visc        = m_viscosity;
  ElementRegionManager::MaterialViewAccessor<arrayView2d<real64>> const & dVisc_dPres = m_dVisc_dPres;

  ElementRegionManager::ConstitutiveRelationAccessor<ConstitutiveBase> constitutiveRelations =
    elemManager->ConstructFullConstitutiveAccessor<ConstitutiveBase>(constitutiveManager);

  // use ArrayView to make capture by value easy in lambdas
  arrayView1d<real64> & presFace      = faceManager->getReference<array1d<real64>>( viewKeyStruct::facePressureString );
  arrayView1d<real64> & densFace      = faceManager->getReference<array1d<real64>>( viewKeyStruct::densityString );
  arrayView1d<real64> & viscFace      = faceManager->getReference<array1d<real64>>( viewKeyStruct::viscosityString );
  arrayView1d<real64> & gravDepthFace = faceManager->getReference<array1d<real64>>( viewKeyStruct::gravityDepthString );

  dataRepository::ManagedGroup const * sets = faceManager->sets();

  // first, evaluate BC to get primary field values (pressure)
//  fsManager->ApplyField(faceManager, viewKeyStruct::facePressure, time + dt);
  fsManager->Apply( time_n + dt,
                    domain,
                    "faceManager",
                    viewKeyStruct::facePressureString,
                    [&] ( FieldSpecificationBase const * const fs,
                    string const &,
                    set<localIndex> const & targetSet,
                    ManagedGroup * const targetGroup,
                    string const fieldName )
  {
    fs->ApplyFieldValue<FieldSpecificationEqual>(targetSet,time_n + dt, targetGroup, fieldName);
  });


  // call constitutive models to get dependent quantities needed for flux (density, viscosity)
  fsManager->Apply( time_n + dt,
                    domain,
                    "faceManager",
                    viewKeyStruct::facePressureString,
                    [&] ( FieldSpecificationBase const * bc,
                    string const &,
                    set<localIndex> const & targetSet,
                    ManagedGroup * const,
                    string const & ) -> void
  {
    for (auto kf : targetSet)
    {
      if (faceGhostRank[kf] >= 0)
        continue;

      // since we don't have models on faces yet, we take them from an adjacent cell
      integer const ke = (elemRegionList[kf][0] >= 0) ? 0 : 1;
      localIndex const er  = elemRegionList[kf][ke];
      localIndex const esr = elemSubRegionList[kf][ke];

      real64 dummy; // don't need derivatives on faces

      SingleFluidBase * fluid = constitutiveRelations[er][esr][m_fluidIndex]->group_cast<SingleFluidBase *>();
      fluid->Compute( presFace[kf], densFace[kf], dummy, viscFace[kf], dummy );
    }
  });

  // *** assembly loop ***

  constexpr localIndex numElems = 2;
  constexpr localIndex maxStencilSize = StencilCollection<CellDescriptor, real64>::MAX_STENCIL_SIZE;

  real64 densWeight[numElems] = { 0.5, 0.5 };

  fsManager->Apply( time_n + dt,
                    domain,
                    "faceManager",
                    viewKeyStruct::facePressureString,
                    [&]( FieldSpecificationBase const * bc,
                    string const & setName,
                    set<localIndex> const &,
                    ManagedGroup * const,
                    string const & )
  {
    if (!sets->hasView(setName) || !fluxApprox->hasBoundaryStencil(setName))
      return;

    FluxApproximationBase::BoundaryStencil const & stencilCollection = fluxApprox->getBoundaryStencil(setName);

    stencilCollection.forAll( GEOSX_LAMBDA ( StencilCollection<PointDescriptor, real64>::Accessor stencil )
    {
      localIndex const stencilSize = stencil.size();

      stackArray1d<globalIndex, maxStencilSize> dofColIndices( stencilSize );

      stackArray1d<real64, numElems> mobility( numElems );
      stackArray1d<real64, numElems> dMobility_dP( numElems );
      stackArray1d<real64, maxStencilSize> dDensMean_dP( stencilSize );
      stackArray1d<real64, maxStencilSize> dFlux_dP( stencilSize );
      stackArray1d<real64, maxStencilSize> localFluxJacobian( stencilSize );

      // clear working arrays
      dDensMean_dP = 0.0;

      // calculate quantities on primary connected points
      real64 densMean = 0.0;
      globalIndex eqnRowIndex = -1;
      localIndex cell_order = -1;
      stencil.forConnected([&] (PointDescriptor const & point, localIndex i) -> void
      {
        real64 density = 0, dDens_dP = 0;
        real64 viscosity = 0, dVisc_dP = 0;
        switch (point.tag)
        {
          case PointDescriptor::Tag::CELL:
          {
            localIndex const er  = point.cellIndex.region;
            localIndex const esr = point.cellIndex.subRegion;
            localIndex const ei  = point.cellIndex.index;

            eqnRowIndex = dofNumber[er][esr][ei];

            density   = dens[er][esr][m_fluidIndex][ei][0];
            dDens_dP  = dDens_dPres[er][esr][m_fluidIndex][ei][0];

            viscosity = visc[er][esr][m_fluidIndex][ei][0];
            dVisc_dP  = dVisc_dPres[er][esr][m_fluidIndex][ei][0];

            cell_order = i; // mark position of the cell in connection for sign consistency later
            break;
          }
          case PointDescriptor::Tag::FACE:
          {
            localIndex const kf = point.faceIndex;

            density = densFace[kf];
            dDens_dP = 0.0;

            viscosity = viscFace[kf];
            dVisc_dP = 0.0;

            break;
          }
          default:
            GEOS_ERROR("Unsupported point type in stencil");
        }

        // mobility
        mobility[i]  = density / viscosity;
        dMobility_dP[i]  = dDens_dP / viscosity - mobility[i] / viscosity * dVisc_dP;

        // average density
        densMean += densWeight[i] * density;
        dDensMean_dP[i] = densWeight[i] * dDens_dP;
      });

      //***** calculation of flux *****

      // compute potential difference MPFA-style
      real64 potDif = 0.0;
      dofColIndices = -1;
      stencil.forAll([&] (PointDescriptor point, real64 w, localIndex i) -> void
      {
        real64 pressure = 0.0, gravD = 0.0;
        switch (point.tag)
        {
          case PointDescriptor::Tag::CELL:
          {
            localIndex const er = point.cellIndex.region;
            localIndex const esr = point.cellIndex.subRegion;
            localIndex const ei = point.cellIndex.index;

            dofColIndices[i] = dofNumber[er][esr][ei];
            pressure = pres[er][esr][ei] + dPres[er][esr][ei];
            gravD = gravDepth[er][esr][ei];

            break;
          }
          case PointDescriptor::Tag::FACE:
          {
            localIndex const kf = point.faceIndex;

            pressure = presFace[kf];
            gravD = gravDepthFace[kf];

            break;
          }
          default:
          GEOS_ERROR("Unsupported point type in stencil");
        }

        real64 const gravTerm = m_gravityFlag ? densMean * gravD : 0.0;
        real64 const dGrav_dP = m_gravityFlag ? dDensMean_dP[i] * gravD : 0.0;

        potDif += w * (pressure + gravTerm);
        dFlux_dP[i] = w * (1.0 + dGrav_dP);
      });

      // upwinding of fluid properties (make this an option?)
      localIndex const k_up = (potDif >= 0) ? 0 : 1;

      // compute the final flux and derivatives
      real64 const flux = mobility[k_up] * potDif;
      for (localIndex ke = 0; ke < stencilSize; ++ke)
        dFlux_dP[ke] *= mobility[k_up];
      dFlux_dP[k_up] += dMobility_dP[k_up] * potDif;

      //***** end flux terms *****

      // populate local flux vector and derivatives
      integer sign = (cell_order == 0 ? 1 : -1);
      real64 const localFlux =  dt * flux * sign;

      integer counter = 0;
      for (localIndex ke = 0; ke < stencilSize; ++ke)
      {
        // compress arrays, skipping face derivatives
        if (dofColIndices[ke] >= 0)
        {
          dofColIndices[counter] = dofColIndices[ke];
          localFluxJacobian[counter] = dt * dFlux_dP[ke] * sign;
          ++counter;
        }
      }

      // Add to global residual/jacobian
      jacobian->SumIntoGlobalValues(1, &eqnRowIndex,
                                    counter, dofColIndices.data(),
                                    localFluxJacobian.data());

      residual->SumIntoGlobalValues(1, &eqnRowIndex, &localFlux);
    });
  });
}

real64
SinglePhaseFlow::
CalculateResidualNorm( EpetraBlockSystem const * const blockSystem,
                       DomainPartition * const domain )
{
  Epetra_FEVector const * const residual = blockSystem->GetResidualVector( BlockIDs::fluidPressureBlock );
  Epetra_Map      const * const rowMap   = blockSystem->GetRowMap( BlockIDs::fluidPressureBlock );

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  ElementRegionManager::ElementViewAccessor<arrayView1d<globalIndex>> const & dofNumber = m_dofNumber;
  ElementRegionManager::ElementViewAccessor<arrayView1d<integer>> const & elemGhostRank = m_elemGhostRank;

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const refPoro = m_porosityRef;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const volume  = m_volume;

  // get a view into local residual vector
  int localSizeInt;
  double* localResidual = nullptr;
  residual->ExtractView(&localResidual, &localSizeInt);

  // compute the norm of local residual scaled by cell pore volume
  real64 localResidualNorm = sumOverElemsInMesh( mesh, [&] ( localIndex const er,
                                                             localIndex const esr,
                                                             localIndex const k ) -> real64
  {
    if (elemGhostRank[er][esr][k] < 0)
    {
      int const lid = rowMap->LID(dofNumber[er][esr][k]);
      real64 const val = localResidual[lid] / (refPoro[er][esr][k] * volume[er][esr][k]);
      return val * val;
    }
    return 0.0;
  });

  // compute global residual norm
  real64 globalResidualNorm;
  MPI_Allreduce(&localResidualNorm, &globalResidualNorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_GEOSX);

  return sqrt(globalResidualNorm);
}

void SinglePhaseFlow::ApplySystemSolution( EpetraBlockSystem const * const blockSystem,
                                           real64 const scalingFactor,
                                           DomainPartition * const domain )
{
  Epetra_Map const * const rowMap        = blockSystem->GetRowMap( BlockIDs::fluidPressureBlock );
  Epetra_FEVector const * const solution = blockSystem->GetSolutionVector( BlockIDs::fluidPressureBlock );

  int dummy;
  double* local_solution = nullptr;
  solution->ExtractView( &local_solution, &dummy );

  applyToSubRegions( domain, [&] ( localIndex er, localIndex esr,
                                   ElementRegion * const region,
                                   ElementSubRegionBase * const subRegion )
  {
    arrayView1d<globalIndex const> const & dofNumber = m_dofNumber[er][esr];
    arrayView1d<integer const> const & elemGhostRank = m_elemGhostRank[er][esr];

    arrayView1d<real64> const & dPres = m_deltaPressure[er][esr];

    forall_in_range<elemPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      if (elemGhostRank[ei] < 0)
      {
        // extract solution and apply to dP
        int const lid = rowMap->LID( integer_conversion<int>( dofNumber[ei] ) );
        dPres[ei] += scalingFactor * local_solution[lid];
      }
    } );
  } );

  // TODO Sync dP once element field syncing is reimplemented.
  std::map<string, string_array > fieldNames;
  fieldNames["elems"].push_back( viewKeyStruct::deltaPressureString );
  CommunicationTools::SynchronizeFields(fieldNames,
                              domain->getMeshBody(0)->getMeshLevel(0),
                              domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );

  applyToSubRegions( domain, [&] ( ElementSubRegionBase * subRegion )
  {
    UpdateState( subRegion );
  } );
}

void SinglePhaseFlow::SolveSystem( EpetraBlockSystem * const blockSystem,
                                   SystemSolverParameters const * const params )
{
  GEOSX_MARK_FUNCTION;

  Epetra_FEVector * const
  solution = blockSystem->GetSolutionVector( BlockIDs::fluidPressureBlock );

  Epetra_FEVector * const
  residual = blockSystem->GetResidualVector( BlockIDs::fluidPressureBlock );

  residual->Scale(-1.0);
  solution->Scale(0.0);

  m_linearSolverWrapper.SolveSingleBlockSystem( blockSystem,
                                                params,
                                                BlockIDs::fluidPressureBlock );

  if( verboseLevel() >= 2 )
  {
    GEOS_LOG_RANK("After SinglePhaseFlow::SolveSystem");
    GEOS_LOG_RANK("\nsolution\n" << *solution);
  }
}

void SinglePhaseFlow::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  applyToSubRegions( domain, [&] ( localIndex er, localIndex esr,
                                   ElementRegion * const region,
                                   ElementSubRegionBase * const subRegion )
  {
    arrayView1d<real64> const & dPres = m_deltaPressure[er][esr];

    forall_in_range<elemPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      dPres[ei] = 0.0;
    } );

    UpdateState( subRegion );
  } );
}

void SinglePhaseFlow::ResetViews(DomainPartition * const domain)
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
  m_deltaVolume =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( viewKeyStruct::deltaVolumeString );

  m_porosityOld =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( viewKeyStruct::oldPorosityString );
  m_densityOld =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( viewKeyStruct::densityString );

  m_pvMult =
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64> >( ConstitutiveBase::viewKeyStruct::poreVolumeMultiplierString,
                                                                                           constitutiveManager );
  m_dPvMult_dPres =
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64> >( ConstitutiveBase::viewKeyStruct::dPVMult_dPresString,
                                                                                           constitutiveManager );
  m_porosity =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( viewKeyStruct::porosityString );

  m_density =
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64> >( SingleFluidBase::viewKeyStruct::densityString,
                                                                                           constitutiveManager );
  m_dDens_dPres =
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64> >( SingleFluidBase::viewKeyStruct::dDens_dPresString,
                                                                                           constitutiveManager );
  m_viscosity =
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64> >( SingleFluidBase::viewKeyStruct::viscosityString,
                                                                                           constitutiveManager );
  m_dVisc_dPres =
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64> >( SingleFluidBase::viewKeyStruct::dVisc_dPresString,
                                                                                           constitutiveManager );

  if (m_poroElasticFlag)
  {
    // TODO where are these strings defined?
    m_totalMeanStressOld = elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( "oldTotalMeanStress" );
    m_totalMeanStress    = elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( "totalMeanStress" );

    m_bulkModulus = elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64> >( "BulkModulus",
                                                                                                       constitutiveManager );
    m_biotCoefficient = elemManager->ConstructFullMaterialViewAccessor<real64>( "BiotCoefficient",
                                                                            constitutiveManager );
  }
}


REGISTER_CATALOG_ENTRY( SolverBase, SinglePhaseFlow, std::string const &, ManagedGroup * const )
} /* namespace geosx */
