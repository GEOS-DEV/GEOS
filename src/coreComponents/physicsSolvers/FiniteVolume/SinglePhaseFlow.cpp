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
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/Fluid/SingleFluidBase.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "MPI_Communications/CommunicationTools.hpp"
#include "MPI_Communications/NeighborCommunicator.hpp"

#include "physicsSolvers/FiniteVolume/SinglePhaseFlowKernels.hpp"

/**
 * @namespace the geosx namespace that encapsulates the majority of the code
 */
namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace SinglePhaseFlowKernels;

SinglePhaseFlow::SinglePhaseFlow( const std::string& name,
                                  ManagedGroup * const parent ):
  FlowSolverBase(name, parent)
{
  m_numDofPerCell = 1;
}

void SinglePhaseFlow::RegisterDataOnMesh(ManagedGroup * const MeshBodies)
{
  FlowSolverBase::RegisterDataOnMesh(MeshBodies);

  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    MeshLevel * meshLevel = ManagedGroup::group_cast<MeshBody *>(mesh.second)->getMeshLevel(0);

    ElementRegionManager * const elemManager = meshLevel->getElemManager();

    elemManager->forElementSubRegions<CellElementSubRegion>( [&]( CellElementSubRegion * const subRegion )
    {
      subRegion->RegisterViewWrapper< array1d<real64> >( viewKeyStruct::pressureString )->setPlotLevel(PlotLevel::LEVEL_0);
      subRegion->RegisterViewWrapper< array1d<real64> >( viewKeyStruct::deltaPressureString );
      subRegion->RegisterViewWrapper< array1d<real64> >( viewKeyStruct::deltaVolumeString );
      subRegion->RegisterViewWrapper< array1d<real64> >( viewKeyStruct::mobilityString );
      subRegion->RegisterViewWrapper< array1d<real64> >( viewKeyStruct::dMobility_dPressureString );
      subRegion->RegisterViewWrapper< array1d<real64> >( viewKeyStruct::porosityString )->setPlotLevel(PlotLevel::LEVEL_1);
      subRegion->RegisterViewWrapper< array1d<real64> >( viewKeyStruct::porosityOldString );
      subRegion->RegisterViewWrapper< array1d<real64> >( viewKeyStruct::densityOldString );
      subRegion->RegisterViewWrapper< array1d<globalIndex> >( viewKeyStruct::blockLocalDofNumberString );
    });

    elemManager->forElementSubRegions<FaceElementSubRegion>( [&]( FaceElementSubRegion * const subRegion )
    {
      subRegion->RegisterViewWrapper< array1d<real64> >( viewKeyStruct::pressureString )->setPlotLevel(PlotLevel::LEVEL_0);
      subRegion->RegisterViewWrapper< array1d<real64> >( viewKeyStruct::deltaPressureString );
      subRegion->RegisterViewWrapper< array1d<real64> >( viewKeyStruct::deltaVolumeString );
      subRegion->RegisterViewWrapper< array1d<real64> >( viewKeyStruct::mobilityString );
      subRegion->RegisterViewWrapper< array1d<real64> >( viewKeyStruct::dMobility_dPressureString );
      subRegion->RegisterViewWrapper< array1d<real64> >( viewKeyStruct::porosityString )->
        setDefaultValue(1.0);
      subRegion->RegisterViewWrapper< array1d<real64> >( viewKeyStruct::porosityOldString )->
        setDefaultValue(1.0);
      subRegion->RegisterViewWrapper< array1d<real64> >( viewKeyStruct::densityOldString );
      subRegion->RegisterViewWrapper< array1d<globalIndex> >( viewKeyStruct::blockLocalDofNumberString );

      subRegion->RegisterViewWrapper< array1d<real64> >( viewKeyStruct::aperture0String )->
        setDefaultValue(0.0);

    } );

    // TODO restrict this to boundary sets
    FaceManager * const faceManager = meshLevel->getFaceManager();
    {
      faceManager->RegisterViewWrapper<array1d<real64> >( viewKeyStruct::facePressureString );
      faceManager->RegisterViewWrapper<array2d<real64> >( viewKeyStruct::faceDensityString )->reference().resizeDimension<1>(1);
      faceManager->RegisterViewWrapper<array2d<real64> >( viewKeyStruct::faceViscosityString )->reference().resizeDimension<1>(1);
      faceManager->RegisterViewWrapper<array1d<real64> >( viewKeyStruct::faceMobilityString );
    }
  }
}

void SinglePhaseFlow::UpdateFluidModel(ManagedGroup * const dataGroup) const
{
  GEOSX_MARK_FUNCTION;

  SingleFluidBase * const fluid = GetConstitutiveModel<SingleFluidBase>( dataGroup, m_fluidName );

  arrayView1d<real64 const> const & pres = dataGroup->getReference< array1d<real64> >( viewKeyStruct::pressureString );
  arrayView1d<real64 const> const & dPres = dataGroup->getReference< array1d<real64> >( viewKeyStruct::deltaPressureString );

  // TODO replace with batch update (need up-to-date pressure and temperature fields)
  forall_in_range<RAJA::seq_exec>( 0, dataGroup->size(), GEOSX_LAMBDA ( localIndex const a )
  {
    fluid->PointUpdate( pres[a] + dPres[a], a, 0 );
  });
  //fluid->BatchUpdate( pres, temp, compFrac );
}

void SinglePhaseFlow::UpdateSolidModel(ManagedGroup * const dataGroup) const
{
  GEOSX_MARK_FUNCTION;

  ConstitutiveBase * const solid = GetConstitutiveModel<ConstitutiveBase>( dataGroup, m_solidName );

  arrayView1d<real64 const> const & pres  = dataGroup->getReference< array1d<real64> >( viewKeyStruct::pressureString );
  arrayView1d<real64 const> const & dPres = dataGroup->getReference< array1d<real64> >( viewKeyStruct::deltaPressureString );

  forall_in_range( 0, dataGroup->size(), GEOSX_LAMBDA ( localIndex const a )
  {
    solid->StateUpdatePointPressure( pres[a] + dPres[a], a, 0 );
  });
}

void SinglePhaseFlow::UpdateMobility( ManagedGroup * const dataGroup ) const
{
  GEOSX_MARK_FUNCTION;

  // output

  arrayView1d<real64> const & mob =
    dataGroup->getReference< array1d<real64> >( viewKeyStruct::mobilityString );

  arrayView1d<real64> const & dMob_dPres =
    dataGroup->getReference< array1d<real64> >( viewKeyStruct::dMobility_dPressureString );

  // input

  SingleFluidBase * const fluid = GetConstitutiveModel<SingleFluidBase>( dataGroup, m_fluidName );

  arrayView2d<real64 const> const & dens =
    fluid->getReference< array2d<real64> >( SingleFluidBase::viewKeyStruct::densityString );

  arrayView2d<real64 const> const & dDens_dPres =
    fluid->getReference< array2d<real64> >( SingleFluidBase::viewKeyStruct::dDens_dPresString );

  arrayView2d<real64 const> const & visc =
    fluid->getReference< array2d<real64> >( SingleFluidBase::viewKeyStruct::viscosityString );

  arrayView2d<real64 const> const & dVisc_dPres =
    fluid->getReference< array2d<real64> >( SingleFluidBase::viewKeyStruct::dVisc_dPresString );

  MobilityKernel::Launch( 0, dataGroup->size(),
                          dens,
                          dDens_dPres,
                          visc,
                          dVisc_dPres,
                          mob,
                          dMob_dPres );
}


void SinglePhaseFlow::UpdateState( ManagedGroup * dataGroup ) const
{
  GEOSX_MARK_FUNCTION;

  UpdateFluidModel( dataGroup );
  UpdateSolidModel( dataGroup );
  UpdateMobility( dataGroup );
}

void SinglePhaseFlow::InitializePostInitialConditions_PreSubGroups( ManagedGroup * const rootGroup )
{
  GEOSX_MARK_FUNCTION;

  FlowSolverBase::InitializePostInitialConditions_PreSubGroups( rootGroup );

  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);
  MeshLevel * mesh = domain->getMeshBody(0)->getMeshLevel(0);

  //TODO this is a hack until the sets are fixed to include ghosts!!
  std::map<string, string_array > fieldNames;
  fieldNames["elems"].push_back( viewKeyStruct::pressureString );

  array1d<NeighborCommunicator> & comms =
    domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors );

  CommunicationTools::SynchronizeFields( fieldNames, mesh, comms );

  ResetViews( domain );

  // Moved the following part from ImplicitStepSetup to here since it only needs to be initialized once
  // They will be updated in ApplySystemSolution and ImplicitStepComplete, respectively

  applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegion * const region,
                                 ElementSubRegionBase * const subRegion )
  {
    UpdateState( subRegion );

    arrayView1d<real64 const> const & poroRef = m_porosityRef[er][esr];
    arrayView2d<real64 const> const & dens    = m_density[er][esr][m_fluidIndex];
    arrayView2d<real64 const> const & pvmult  = m_pvMult[er][esr][m_solidIndex];

    arrayView1d<real64> const & poro = m_porosity[er][esr];
    arrayView1d<real64> const & densOld = m_densityOld[er][esr];
    arrayView1d<real64> const & poroOld = m_porosityOld[er][esr];

    if( pvmult.size() == poro.size() )
    {
      forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
      {
        densOld[ei] = dens[ei][0];
        poro[ei] = poroRef[ei] * pvmult[ei][0];
        poroOld[ei] = poro[ei];
      } );
    }
    else
    {
      forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
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

  real64 dt_return;

  ImplicitStepSetup( time_n, dt, domain, m_dofManager, m_matrix, m_rhs, m_solution );

  // currently the only method is implicit time integration
  dt_return = this->NonlinearImplicitStep( time_n, dt, cycleNumber, domain, m_dofManager, m_matrix, m_rhs, m_solution );

  // final step for completion of timestep. typically secondary variable updates and cleanup.
  ImplicitStepComplete( time_n, dt_return, domain );

  return dt_return;
}


void SinglePhaseFlow::ImplicitStepSetup( real64 const & time_n,
                                         real64 const & dt,
                                         DomainPartition * const domain,
                                         DofManager & dofManager,
                                         ParallelMatrix & matrix,
                                         ParallelVector & rhs,
                                         ParallelVector & solution )
{
  ResetViews( domain );

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegion * const region,
                                 ElementSubRegionBase * const subRegion )
  {
    arrayView2d<real64 const> const & dens = m_density[er][esr][m_fluidIndex];
    arrayView1d<real64 const> const & poro = m_porosity[er][esr];

    arrayView1d<real64> const & dPres   = m_deltaPressure[er][esr];
    arrayView1d<real64> const & dVol    = m_deltaVolume[er][esr];
    arrayView1d<real64> const & densOld = m_densityOld[er][esr];
    arrayView1d<real64> const & poroOld = m_porosityOld[er][esr];

    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      dPres[ei] = 0.0;
      dVol[ei] = 0.0;
      densOld[ei] = dens[ei][0];
      poroOld[ei] = poro[ei];
    } );
  } );

  mesh->getElemManager()->
      forElementSubRegionsComplete<FaceElementSubRegion>( m_targetRegions,
                                                          [&] ( localIndex const er,
                                                                localIndex const esr,
                                                                ElementRegion *,
                                                                FaceElementSubRegion * subRegion )
  {
    arrayView1d<real64> const & aper0 = subRegion->getReference<array1d<real64>>( viewKeyStruct::aperture0String );
    arrayView1d<real64 const> const & aper = m_elementAperture[er][esr];

    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      aper0[ei] = aper[ei];
    } );


    UpdateMobility( subRegion );
  } );

  // setup dof numbers and linear system
  SetupSystem( domain, dofManager, matrix, rhs, solution );
}

void SinglePhaseFlow::ImplicitStepComplete( real64 const & time_n,
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
    arrayView1d<real64> const & vol  = m_volume[er][esr];

    arrayView1d<real64 const> const & dPres = m_deltaPressure[er][esr];
    arrayView1d<real64 const> const & dVol  = m_deltaVolume[er][esr];

    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      pres[ei] += dPres[ei];
      vol[ei] += dVol[ei];
    } );
  } );
}



void SinglePhaseFlow::SetNumRowsAndTrilinosIndices( MeshLevel * const meshLevel,
                                                    localIndex & numLocalRows,
                                                    globalIndex & numGlobalRows,
                                                    localIndex offset ) const
{
  int numMpiProcesses;
  MPI_Comm_size( MPI_COMM_GEOSX, &numMpiProcesses );

  int thisMpiProcess = 0;
  MPI_Comm_rank( MPI_COMM_GEOSX, &thisMpiProcess );

  localIndex numLocalRowsToSend = numLocalRows;
  array1d<localIndex> gather( numMpiProcesses );

  // communicate the number of local rows to each process
  CommunicationTools::allGather( numLocalRowsToSend, gather );

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

void SinglePhaseFlow::SetupSystem( DomainPartition * const domain,
                                   DofManager & dofManager,
                                   ParallelMatrix & matrix,
                                   ParallelVector & rhs,
                                   ParallelVector & solution )
{
  // assume that there is only a single MeshLevel for now
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  // for this solver, the dof are on the cell center, and the row corrosponds to an element
  localIndex numGhostRows  = 0;
  localIndex numLocalRows  = 0;
  globalIndex numGlobalRows = 0;

  // get the number of local elements, and ghost elements...i.e. local rows and ghost rows
  applyToSubRegions( mesh, [&] ( ElementSubRegionBase * const subRegion )
  {
    localIndex subRegionGhosts = subRegion->GetNumberOfGhosts();
    numGhostRows += subRegionGhosts;
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

  matrix.createWithLocalSize( numLocalRows, numLocalRows, 7, MPI_COMM_GEOSX );
  SetSparsityPattern( domain, &matrix );
  matrix.close();

  rhs.createWithLocalSize( numLocalRows, MPI_COMM_GEOSX );
  solution.createWithLocalSize( numLocalRows, MPI_COMM_GEOSX );
}

void SinglePhaseFlow::SetSparsityPattern( DomainPartition const * const domain,
                                          ParallelMatrix * const matrix ) const
{
  MeshLevel const * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  ElementRegionManager::ElementViewAccessor< arrayView1d<globalIndex> > const & dofNumber = m_dofNumber;

  NumericalMethodsManager const * numericalMethodManager =
    domain->getParent()->GetGroup<NumericalMethodsManager>( keys::numericalMethodsManager );

  FiniteVolumeManager const * fvManager =
    numericalMethodManager->GetGroup<FiniteVolumeManager>( keys::finiteVolumeManager );

  FluxApproximationBase const * fluxApprox = fvManager->getFluxApproximation( m_discretizationName );

  //**** loop over all faces. Fill in sparsity for all pairs of DOF/elem that are connected by face
  fluxApprox->forCellStencils( [&]( auto const & stencil )
  {
    typedef TYPEOFREF( stencil ) STENCIL_TYPE;
    typename STENCIL_TYPE::IndexContainerViewConstType const & eri = stencil.getElementRegionIndices();
    typename STENCIL_TYPE::IndexContainerViewConstType const & esri = stencil.getElementSubRegionIndices();
    typename STENCIL_TYPE::IndexContainerViewConstType const & ei = stencil.getElementIndices();

    forall_in_range<serialPolicy>( 0, stencil.size(), GEOSX_LAMBDA( localIndex iconn )
    {
      localIndex const stencilSize = stencil.stencilSize( iconn );
      localIndex const numFluxElems = stencilSize;
      stackArray1d<globalIndex, STENCIL_TYPE::NUM_POINT_IN_FLUX> dofIndexRow( numFluxElems );
      stackArray1d<globalIndex, STENCIL_TYPE::MAX_STENCIL_SIZE> dofIndexCol( stencilSize );

      stackArray2d<real64, STENCIL_TYPE::MAX_STENCIL_SIZE * STENCIL_TYPE::NUM_POINT_IN_FLUX>
      values( numFluxElems, stencilSize );

      values = 1.0;

      for( localIndex i = 0; i < numFluxElems; ++i )
      {
        dofIndexRow[i] = dofNumber[eri( iconn, i )][esri( iconn, i )][ei( iconn, i )];
      }

      for( localIndex i = 0; i < stencilSize; ++i )
      {
        dofIndexCol[i] = dofNumber[eri( iconn, i )][esri( iconn, i )][ei( iconn, i )];
      }

      matrix->insert( dofIndexRow.data(),
                      dofIndexCol.data(),
                      values.data(),
                      numFluxElems,
                      stencilSize );
    } );
  } );


  // add isolated elements not captured in the flux stencil
  applyToSubRegions( meshLevel, [&] ( localIndex er, localIndex esr,
                                      ElementRegion const * const region,
                                      ElementSubRegionBase const * const subRegion )
  {
    arrayView1d<integer const> const & elemGhostRank = m_elemGhostRank[er][esr];
    arrayView1d<globalIndex const> const & dofNumberSub = m_dofNumber[er][esr];

    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex const a )
    {
      if (elemGhostRank[a] < 0)
      {
        globalIndex const dofIndex = dofNumberSub[a];
        matrix->insert( dofIndex, dofIndex, 1.0 );
      }
    } );
  } );
}

void SinglePhaseFlow::AssembleSystem( real64 const time_n,
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

  if (m_poroElasticFlag)
  {
    AssembleAccumulationTerms<true>( domain, &dofManager, &matrix, &rhs );
  }
  else
  {
    AssembleAccumulationTerms<false>( domain, &dofManager, &matrix, &rhs );
  }

  AssembleFluxTerms( time_n, dt, domain, &dofManager, &matrix, &rhs );

  matrix.close();
  rhs.close();

  if( verboseLevel() == 2 )
  {
    GEOS_LOG_RANK_0( "After SinglePhaseFlow::AssembleSystem: " );
    GEOS_LOG_RANK_0( "\nJacobian:\n" );
    matrix.print(std::cout);
    rhs.print(std::cout);
  }

  if( verboseLevel() >= 3 )
  {
    SystemSolverParameters * const solverParams = getSystemSolverParameters();
    integer newtonIter = solverParams->numNewtonIterations();

    string filename_mat = "matrix_" + std::to_string( time_n ) + "_" + std::to_string( newtonIter ) + ".mtx";
    matrix.write( filename_mat, true );

    string filename_rhs = "rhs_" + std::to_string( time_n ) + "_" + std::to_string( newtonIter ) + ".mtx";
    rhs.write( filename_rhs, true );

    GEOS_LOG_RANK_0( "After SinglePhaseFlow::AssembleSystem" );
    GEOS_LOG_RANK_0( "Jacobian: written to " << filename_mat );
    GEOS_LOG_RANK_0( "Residual: written to " << filename_rhs );
  }
}

template< bool ISPORO >
void SinglePhaseFlow::AccumulationLaunch( localIndex const er,
                                          localIndex const esr,
                                          CellElementSubRegion const * const subRegion,
                                          ParallelMatrix * const matrix,
                                          ParallelVector * const rhs )
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
  arrayView1d<real64 const> const & bulkModulus        = m_poroElasticFlag ? m_bulkModulus[er][esr][m_solidIndex] : poroOld;
  real64 const & biotCoefficient                       = m_poroElasticFlag ? m_biotCoefficient[er][esr][m_solidIndex] : 0;



  forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
  {
    if (elemGhostRank[ei] < 0)
    {
      real64 localAccum, localAccumJacobian;
      globalIndex const elemDOF = dofNumber[ei];

      AccumulationKernel<CellElementSubRegion>::template Compute<ISPORO>( dPres[ei],
                                           dens[ei][0],
                                           densOld[ei],
                                           dDens_dPres[ei][0],
                                           volume[ei],
                                           dVol[ei],
                                           poroRef[ei],
                                           poroOld[ei],
                                           pvmult[ei][0],
                                           dPVMult_dPres[ei][0],
                                           biotCoefficient,
                                           bulkModulus[ei],
                                           totalMeanStress[ei],
                                           oldTotalMeanStress[ei],
                                           poro[ei],
                                           localAccum,
                                           localAccumJacobian );

        // add contribution to global residual and jacobian
      matrix->add( elemDOF, elemDOF, localAccumJacobian );
      rhs->add( elemDOF, localAccum );
    }
  } );

}

template< bool ISPORO >
void SinglePhaseFlow::AccumulationLaunch( localIndex const er,
                                          localIndex const esr,
                                          FaceElementSubRegion const * const subRegion,
                                          ParallelMatrix * const matrix,
                                          ParallelVector * const rhs )
{


  arrayView1d<integer const>     const & elemGhostRank = m_elemGhostRank[er][esr];
  arrayView1d<globalIndex const> const & dofNumber     = m_dofNumber[er][esr];

  arrayView1d<real64 const> const & densOld       = m_densityOld[er][esr];
  arrayView1d<real64 const> const & volume        = m_volume[er][esr];
  arrayView1d<real64 const> const & dVol          = m_deltaVolume[er][esr];
  arrayView2d<real64 const> const & dens          = m_density[er][esr][m_fluidIndex];
  arrayView2d<real64 const> const & dDens_dPres   = m_dDens_dPres[er][esr][m_fluidIndex];

  forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
  {
    if (elemGhostRank[ei] < 0)
    {
      real64 localAccum, localAccumJacobian;
      globalIndex const elemDOF = dofNumber[ei];

      AccumulationKernel<FaceElementSubRegion>::template Compute<ISPORO>( dens[ei][0],
                                                                          densOld[ei],
                                                                          dDens_dPres[ei][0],
                                                                          volume[ei],
                                                                          dVol[ei],
                                                                          localAccum,
                                                                          localAccumJacobian );

      // add contribution to global residual and jacobian
      matrix->add( elemDOF, elemDOF, localAccumJacobian );
      rhs->add( elemDOF, localAccum );
    }
  } );
}

template< bool ISPORO >
void SinglePhaseFlow::AssembleAccumulationTerms( DomainPartition const * const domain,
                                                 DofManager const * const dofManager,
                                                 ParallelMatrix * const matrix,
                                                 ParallelVector * const rhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);


  ElementRegionManager const * const elemManager = mesh->getElemManager();

  elemManager->forElementSubRegionsComplete<CellElementSubRegion,
                                            FaceElementSubRegion>( this->m_targetRegions,
                                                                   [&] ( localIndex er,
                                                                         localIndex esr,
                                                                         ElementRegion const * const region,
                                                                         auto const * const subRegion )
  {
    AccumulationLaunch<ISPORO>( er, esr, subRegion, matrix, rhs );
  } );
}


void SinglePhaseFlow::AssembleFluxTerms( real64 const time_n,
                                         real64 const dt,
                                         DomainPartition const * const domain,
                                         DofManager const * const dofManager,
                                         ParallelMatrix * const matrix,
                                         ParallelVector * const rhs )
{
  GEOSX_MARK_FUNCTION;

  NumericalMethodsManager const * numericalMethodManager =
    domain->getParent()->GetGroup<NumericalMethodsManager>( keys::numericalMethodsManager );

  FiniteVolumeManager const * fvManager =
    numericalMethodManager->GetGroup<FiniteVolumeManager>( keys::finiteVolumeManager );

  FluxApproximationBase const * fluxApprox = fvManager->getFluxApproximation( m_discretizationName );

  FluxKernel::ElementView < arrayView1d<globalIndex const> > const & dofNumber = m_dofNumber.toViewConst();

  FluxKernel::ElementView < arrayView1d<real64 const> > const & dPres       = m_deltaPressure.toViewConst();
  FluxKernel::ElementView < arrayView1d<real64 const> > const & pres        = m_pressure.toViewConst();
  FluxKernel::ElementView < arrayView1d<real64 const> > const & gravDepth   = m_gravDepth.toViewConst();
  FluxKernel::MaterialView< arrayView2d<real64 const> > const & dens        = m_density.toViewConst();
  FluxKernel::MaterialView< arrayView2d<real64 const> > const & dDens_dPres = m_dDens_dPres.toViewConst();
  FluxKernel::ElementView < arrayView1d<real64 const> > const & mob         = m_mobility.toViewConst();
  FluxKernel::ElementView < arrayView1d<real64 const> > const & dMob_dPres  = m_dMobility_dPres.toViewConst();

  FluxKernel::ElementView < arrayView1d<real64 const> > const & aperture0  = m_elementAperture0.toViewConst();
  FluxKernel::ElementView < arrayView1d<real64 const> > const & aperture  = m_elementAperture.toViewConst();

  integer const gravityFlag = m_gravityFlag;
  localIndex const fluidIndex = m_fluidIndex;


  fluxApprox->forCellStencils( [&]( auto const & stencil )
  {

//    typedef TYPEOFREF( stencil ) STENCIL_TYPE;

    FluxKernel::Launch( stencil,
                        dt,
                        fluidIndex,
                        gravityFlag,
                        dofNumber,
                        pres,
                        dPres,
                        gravDepth,
                        dens,
                        dDens_dPres,
                        mob,
                        dMob_dPres,
                        aperture0,
                        aperture,
                        matrix,
                        rhs );
  });

}

void
SinglePhaseFlow::ApplyBoundaryConditions( real64 const time_n,
                                          real64 const dt,
                                          DomainPartition * const domain,
                                          DofManager const & dofManager,
                                          ParallelMatrix & matrix,
                                          ParallelVector & rhs )
{
  GEOSX_MARK_FUNCTION;

  matrix.open();
  rhs.open();

  FieldSpecificationManager * fsManager = FieldSpecificationManager::get();

  // call the BoundaryConditionManager::ApplyField function that will check to see
  // if the boundary condition should be applied to this subregion
  fsManager->Apply( time_n + dt, domain, "ElementRegions", "FLUX",
                    [&]( FieldSpecificationBase const * const fs,
                         string const &,
                         set<localIndex> const & lset,
                         ManagedGroup * subRegion,
                         string const & ) -> void
  {
    arrayView1d<globalIndex const> const &
    dofNumber = subRegion->getReference< array1d<globalIndex> >( viewKeyStruct::blockLocalDofNumberString );

    fs->ApplyBoundaryConditionToSystem<FieldSpecificationAdd, LAInterface>( lset,
                                                                            true,
                                                                            time_n + dt,
                                                                            dt,
                                                                            subRegion,
                                                                            dofNumber,
                                                                            1,
                                                                            matrix,
                                                                            rhs,
                                                                            [&]( localIndex const a ) -> real64
    {
      return 0;
    } );

  } );


  fsManager->Apply( time_n + dt, domain, "ElementRegions", viewKeyStruct::pressureString,
                    [&]( FieldSpecificationBase const * const fs,
                         string const &,
                         set<localIndex> const & lset,
                         ManagedGroup * subRegion,
                         string const & ) -> void
  {
    arrayView1d<globalIndex const> const &
    dofNumber = subRegion->getReference< array1d<globalIndex> >( viewKeyStruct::blockLocalDofNumberString );

    //for now assume all the non-flux boundary conditions are Dirichlet type BC.

    arrayView1d<real64 const> const &
    pres = subRegion->getReference<array1d<real64> >( viewKeyStruct::pressureString );

    arrayView1d<real64 const> const &
    dPres = subRegion->getReference<array1d<real64> >( viewKeyStruct::deltaPressureString );

    // call the application of the boundary condition to alter the matrix and rhs
    fs->ApplyBoundaryConditionToSystem<FieldSpecificationEqual, LAInterface>( lset,
                                                                              false,
                                                                              time_n + dt,
                                                                              subRegion,
                                                                              dofNumber,
                                                                              1,
                                                                              matrix,
                                                                              rhs,
                                                                              [&]( localIndex const a ) -> real64
    {
      return pres[a] + dPres[a];
    });
  });


  ApplyFaceDirichletBC_implicit( time_n, dt, &dofManager, domain, &matrix, &rhs );

  matrix.close();
  rhs.close();

  if( verboseLevel() == 2 )
  {
    GEOS_LOG_RANK_0( "After SinglePhaseFlow::ApplyBoundaryConditions" );
    GEOS_LOG_RANK_0("\nJacobian:\n");
    matrix.print(std::cout);
    rhs.print(std::cout);
  }

  if( verboseLevel() >= 3 )
  {
    SystemSolverParameters * const solverParams = getSystemSolverParameters();
    integer newtonIter = solverParams->numNewtonIterations();

    string filename_mat = "matrix_bc_" + std::to_string( time_n ) + "_" + std::to_string( newtonIter ) + ".mtx";
    matrix.write( filename_mat, true );

    string filename_rhs = "rhs_bc_" + std::to_string( time_n ) + "_" + std::to_string( newtonIter ) + ".mtx";
    rhs.write( filename_rhs, true );

    GEOS_LOG_RANK_0( "After SinglePhaseFlow::ApplyBoundaryConditions" );
    GEOS_LOG_RANK_0( "Jacobian: written to " << filename_mat );
    GEOS_LOG_RANK_0( "Residual: written to " << filename_rhs );
  }
}

void SinglePhaseFlow::ApplyFaceDirichletBC_implicit( real64 const time_n,
                                                     real64 const dt,
                                                     DofManager const * const dofManager,
                                                     DomainPartition * const domain,
                                                     ParallelMatrix * const matrix,
                                                     ParallelVector * const rhs )
{
  FieldSpecificationManager * fsManager = FieldSpecificationManager::get();
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();
  FaceManager * const faceManager = mesh->getFaceManager();

  arrayView2d<localIndex> const & elemRegionList     = faceManager->elementRegionList();
  arrayView2d<localIndex> const & elemSubRegionList  = faceManager->elementSubRegionList();

  ConstitutiveManager * const constitutiveManager =
    domain->GetGroup<ConstitutiveManager>(keys::ConstitutiveManager);

  NumericalMethodsManager * const numericalMethodManager =
    domain->getParent()->GetGroup<NumericalMethodsManager>(keys::numericalMethodsManager);

  FiniteVolumeManager * const fvManager = numericalMethodManager->GetGroup<FiniteVolumeManager>(keys::finiteVolumeManager);

  FluxApproximationBase const * const fluxApprox = fvManager->getFluxApproximation( m_discretizationName );

  // make a list of region indices to be included
  set<localIndex> regionFilter;
  for (string const & regionName : m_targetRegions)
  {
    regionFilter.insert( elemManager->GetRegions().getIndex( regionName ) );
  }

  ElementRegionManager::ElementViewAccessor< arrayView1d<globalIndex> > const & dofNumber = m_dofNumber;

  ElementRegionManager::ElementViewAccessor< arrayView1d<real64> >  const & pres        = m_pressure;
  ElementRegionManager::ElementViewAccessor< arrayView1d<real64> >  const & dPres       = m_deltaPressure;
  ElementRegionManager::ElementViewAccessor< arrayView1d<real64> >  const & gravDepth   = m_gravDepth;
  ElementRegionManager::MaterialViewAccessor< arrayView2d<real64> > const & dens        = m_density;
  ElementRegionManager::MaterialViewAccessor< arrayView2d<real64> > const & dDens_dPres = m_dDens_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView1d<real64> >  const & mob         = m_mobility;
  ElementRegionManager::ElementViewAccessor< arrayView1d<real64> >  const & dMob_dPres  = m_dMobility_dPres;

  ElementRegionManager::ConstitutiveRelationAccessor<ConstitutiveBase> constitutiveRelations =
    elemManager->ConstructFullConstitutiveAccessor<ConstitutiveBase>(constitutiveManager);

  // use ArrayView to make capture by value easy in lambdas
  arrayView1d<real64 const> const & presFace      = faceManager->getReference< array1d<real64> >( viewKeyStruct::facePressureString );
  arrayView2d<real64>       const & densFace      = faceManager->getReference< array2d<real64> >( viewKeyStruct::faceDensityString );
  arrayView2d<real64>       const & viscFace      = faceManager->getReference< array2d<real64> >( viewKeyStruct::faceViscosityString );
  arrayView1d<real64>       const & mobFace       = faceManager->getReference< array1d<real64> >( viewKeyStruct::faceMobilityString );
  arrayView1d<real64 const> const & gravDepthFace = faceManager->getReference< array1d<real64> >( viewKeyStruct::gravityDepthString );

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
                          string const & )
  {
    for (auto kf : targetSet)
    {
      // since we don't have models on faces yet, we take them from an adjacent cell
      integer ke;
      for (ke = 0; ke < 2; ++ke)
      {
        if (elemRegionList[kf][ke] >= 0 && regionFilter.contains(elemRegionList[kf][ke]))
        {
          break;
        }
      }
      GEOS_ERROR_IF( ke > 1, "Face not adjacent to target regions: " << kf );
      localIndex const er  = elemRegionList[kf][ke];
      localIndex const esr = elemSubRegionList[kf][ke];

      real64 dummy; // don't need derivatives on faces

      SingleFluidBase * fluid = constitutiveRelations[er][esr][m_fluidIndex]->group_cast<SingleFluidBase *>();
      fluid->Compute( presFace[kf], densFace[kf][0], dummy, viscFace[kf][0], dummy );
    }

    MobilityKernel::Launch( targetSet, densFace, viscFace, mobFace );
  });

  // *** assembly loop ***

  constexpr localIndex numElems = CellElementStencilTPFA::NUM_POINT_IN_FLUX;
  constexpr localIndex maxStencilSize = CellElementStencilTPFA::MAX_STENCIL_SIZE;

  real64 densWeight[numElems] = { 0.5, 0.5 };

  fsManager->Apply( time_n + dt,
                    domain,
                    "faceManager",
                    viewKeyStruct::facePressureString,
                    [&] ( FieldSpecificationBase const * bc,
                          string const & setName,
                          set<localIndex> const &,
                          ManagedGroup * const,
                          string const & )
  {
    if (!sets->hasView(setName) || !fluxApprox->hasBoundaryStencil(setName))
      return;

    FluxApproximationBase::BoundaryStencil const & stencil = fluxApprox->getBoundaryStencil(setName);
    ArrayOfArraysView<FluxApproximationBase::BoundaryStencil::Entry const, true> const & connections = stencil.getConnections();

    forall_in_range<serialPolicy>( 0, connections.size(), GEOSX_LAMBDA ( localIndex iconn )
    {
      localIndex const stencilSize = connections.sizeOfArray(iconn);

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

      for (localIndex i = 0; i < numElems; ++i)
      {
        PointDescriptor const & point = connections(iconn, i).index;

        real64 density = 0, dDens_dP = 0;
        switch (point.tag)
        {
          case PointDescriptor::Tag::CELL:
          {
            localIndex const er  = point.cellIndex.region;
            localIndex const esr = point.cellIndex.subRegion;
            localIndex const ei  = point.cellIndex.index;

            eqnRowIndex = dofNumber[er][esr][ei];

            density  = dens[er][esr][m_fluidIndex][ei][0];
            dDens_dP = dDens_dPres[er][esr][m_fluidIndex][ei][0];

            mobility[i]     = mob[er][esr][ei];
            dMobility_dP[i] = dMob_dPres[er][esr][ei];

            cell_order = i; // mark position of the cell in connection for sign consistency later
            break;
          }
          case PointDescriptor::Tag::FACE:
          {
            density  = densFace[point.faceIndex][0];
            dDens_dP = 0.0;

            mobility[i]     = mobFace[point.faceIndex];
            dMobility_dP[i] = 0.0;
            break;
          }
          default:
            GEOS_ERROR("Unsupported point type in stencil");
        }

        // average density
        densMean += densWeight[i] * density;
        dDensMean_dP[i] = densWeight[i] * dDens_dP;
      }

      //***** calculation of flux *****

      // compute potential difference MPFA-style
      real64 potDif = 0.0;
      dofColIndices = -1;
      for (localIndex i = 0; i < stencilSize; ++i)
      {
        FluxApproximationBase::BoundaryStencil::Entry const & entry = connections(iconn, i);
        PointDescriptor const & point = entry.index;

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

        potDif += entry.weight * (pressure + gravTerm);
        dFlux_dP[i] = entry.weight * (1.0 + dGrav_dP);
      }

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
      matrix->add( eqnRowIndex, dofColIndices.data(), localFluxJacobian.data(), counter );
      rhs->add( eqnRowIndex, localFlux );
    } );
  } );
}

real64 SinglePhaseFlow::CalculateResidualNorm( DomainPartition const * const domain,
                                               DofManager const & dofManager,
                                               ParallelVector const & rhs )
{
  MeshLevel const * const mesh = domain->getMeshBody(0)->getMeshLevel(0);

  // get a view into local residual vector
  real64 * localResidual = nullptr;
  rhs.extractLocalVector( &localResidual );

  // compute the norm of local residual scaled by cell pore volume
  real64 localResidualNorm = 0.0;

  applyToSubRegions( mesh, [&] ( localIndex const er, localIndex const esr,
                                 ElementRegion const * const region,
                                 ElementSubRegionBase const * const subRegion )
  {
    arrayView1d<integer const> const & elemGhostRank = m_elemGhostRank[er][esr];
    arrayView1d<globalIndex const> const & dofNumber = m_dofNumber[er][esr];
    arrayView1d<real64 const> const & refPoro        = m_porosityRef[er][esr];
    arrayView1d<real64 const> const & volume         = m_volume[er][esr];
    arrayView1d<real64 const> const & dVol           = m_deltaVolume[er][esr];
//    arrayView1d<real64 const> const & dens           = m_density[er][esr][m_fluidIndex].dimReduce();
    arrayView2d<real64 const> const & dens           = m_density[er][esr][m_fluidIndex];

    localIndex const subRegionSize = subRegion->size();
    for ( localIndex a = 0; a < subRegionSize; ++a )
    {
      if (elemGhostRank[a] < 0)
      {
        localIndex const lid = rhs.getLocalRowID( dofNumber[a] );
        real64 const val = localResidual[lid] / (refPoro[a] * dens[a][0] * ( volume[a] + dVol[a]));
        localResidualNorm += val * val;
      }
    }
  });

  // compute global residual norm
  real64 globalResidualNorm;
  MPI_Allreduce(&localResidualNorm, &globalResidualNorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_GEOSX);

  return sqrt(globalResidualNorm);
}

void SinglePhaseFlow::ApplySystemSolution( DofManager const & dofManager,
                                           ParallelVector const & solution,
                                           real64 const scalingFactor,
                                           DomainPartition * const domain )
{
  MeshLevel * mesh = domain->getMeshBody(0)->getMeshLevel(0);

  real64 * localSolution = nullptr;
  solution.extractLocalVector( &localSolution );

  applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegion * const region,
                                 ElementSubRegionBase * const subRegion )
  {
    arrayView1d<globalIndex const> const & dofNumber = m_dofNumber[er][esr];
    arrayView1d<integer const> const & elemGhostRank = m_elemGhostRank[er][esr];

    arrayView1d<real64> const & dPres = m_deltaPressure[er][esr];

    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      if( elemGhostRank[ei] < 0 )
      {
        // extract solution and apply to dP
        localIndex const lid = solution.getLocalRowID( dofNumber[ei] );
        dPres[ei] += scalingFactor * localSolution[lid];
      }
    } );
  } );

  std::map<string, string_array> fieldNames;
  fieldNames["elems"].push_back( viewKeyStruct::deltaPressureString );

  array1d<NeighborCommunicator> & comms =
    domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors );

  CommunicationTools::SynchronizeFields( fieldNames, mesh, comms );

  applyToSubRegions( mesh, [&] ( ElementSubRegionBase * subRegion )
  {
    UpdateState( subRegion );
  } );
}

void SinglePhaseFlow::SolveSystem( DofManager const & dofManager,
                                   ParallelMatrix & matrix,
                                   ParallelVector & rhs,
                                   ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  rhs.scale( -1.0 );
  solution.zero();

  SolverBase::SolveSystem( dofManager, matrix, rhs, solution );

  if( verboseLevel() == 2 )
  {
    GEOS_LOG_RANK_0("After SinglePhaseFlow::SolveSystem\nSolution\n");
    GEOS_LOG(solution);
  }
}

void SinglePhaseFlow::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  MeshLevel * mesh = domain->getMeshBody(0)->getMeshLevel(0);

  applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegion * const region,
                                 ElementSubRegionBase * const subRegion )
  {
    arrayView1d<real64> const & dPres = m_deltaPressure[er][esr];

    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      dPres[ei] = 0.0;
    } );

    UpdateState( subRegion );
  } );
}

void SinglePhaseFlow::ResetViews( DomainPartition * const domain )
{
  FlowSolverBase::ResetViews( domain );

  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = mesh->getElemManager();
  ConstitutiveManager * const constitutiveManager = domain->getConstitutiveManager();

  m_dofNumber =
    elemManager->ConstructViewAccessor< array1d<globalIndex>, arrayView1d<globalIndex> >( viewKeyStruct::blockLocalDofNumberString );
  m_pressure =
    elemManager->ConstructViewAccessor< array1d<real64>, arrayView1d<real64> >( viewKeyStruct::pressureString );
  m_deltaPressure =
    elemManager->ConstructViewAccessor< array1d<real64>, arrayView1d<real64> >( viewKeyStruct::deltaPressureString );
  m_deltaVolume =
    elemManager->ConstructViewAccessor< array1d<real64>, arrayView1d<real64> >( viewKeyStruct::deltaVolumeString );

  m_mobility =
    elemManager->ConstructViewAccessor< array1d<real64>, arrayView1d<real64> >( viewKeyStruct::mobilityString );
  m_dMobility_dPres =
    elemManager->ConstructViewAccessor< array1d<real64>, arrayView1d<real64> >( viewKeyStruct::dMobility_dPressureString );

  m_porosityOld =
    elemManager->ConstructViewAccessor< array1d<real64>, arrayView1d<real64> >( viewKeyStruct::porosityOldString );
  m_densityOld =
    elemManager->ConstructViewAccessor< array1d<real64>, arrayView1d<real64> >( viewKeyStruct::densityOldString );

  m_pvMult =
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64> >( ConstitutiveBase::viewKeyStruct::poreVolumeMultiplierString,
                                                                                           constitutiveManager );
  m_dPvMult_dPres =
    elemManager->ConstructFullMaterialViewAccessor<array2d<real64>, arrayView2d<real64> >( ConstitutiveBase::viewKeyStruct::dPVMult_dPresString,
                                                                                           constitutiveManager );
  m_porosity =
    elemManager->ConstructViewAccessor< array1d<real64>, arrayView1d<real64> >( viewKeyStruct::porosityString );

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
    m_totalMeanStressOld = elemManager->ConstructViewAccessor< array1d<real64>, arrayView1d<real64> >( "oldTotalMeanStress" );
    m_totalMeanStress    = elemManager->ConstructViewAccessor< array1d<real64>, arrayView1d<real64> >( "totalMeanStress" );

    m_bulkModulus = elemManager->ConstructFullMaterialViewAccessor<array1d<real64>, arrayView1d<real64> >( "BulkModulus",
                                                                                                       constitutiveManager );
    m_biotCoefficient = elemManager->ConstructFullMaterialViewAccessor<real64>( "BiotCoefficient",
                                                                            constitutiveManager );
  }
}


REGISTER_CATALOG_ENTRY( SolverBase, SinglePhaseFlow, std::string const &, ManagedGroup * const )
} /* namespace geosx */
