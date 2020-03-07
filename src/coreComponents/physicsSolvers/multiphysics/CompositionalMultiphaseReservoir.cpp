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
 * @file CompositionalMultiphaseReservoir.cpp
 *
 */


#include "CompositionalMultiphaseReservoir.hpp"

#include "common/TimingMacros.hpp"
#include "constitutive/fluid/MultiFluidBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFlow.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseWell.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
  
CompositionalMultiphaseReservoir::CompositionalMultiphaseReservoir( const std::string& name,
                                                                    Group * const parent ):
  ReservoirSolverBase(name,parent)
{
}

CompositionalMultiphaseReservoir::~CompositionalMultiphaseReservoir()
{
}

void CompositionalMultiphaseReservoir::SetupDofs( DomainPartition const * const domain,
                                                  DofManager & dofManager ) const
{
  m_flowSolver->SetupDofs( domain, dofManager );
  m_wellSolver->SetupDofs( domain, dofManager );

  // TODO: add coupling when dofManager can support perforation connectors
}

void CompositionalMultiphaseReservoir::SetupSystem( DomainPartition * const domain,
                                                          DofManager & dofManager,
                                                          ParallelMatrix & matrix,
                                                          ParallelVector & rhs,
                                                          ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager const * const elemManager = meshLevel->getElemManager();

  dofManager.setMesh( domain, 0, 0 );
  SetupDofs( domain, dofManager );
  dofManager.reorderByRank();

  localIndex const numLocalDof = dofManager.numLocalDofs();

  matrix.createWithLocalSize( numLocalDof, numLocalDof, 8, MPI_COMM_GEOSX );
  rhs.createWithLocalSize( numLocalDof, MPI_COMM_GEOSX );
  solution.createWithLocalSize( numLocalDof, MPI_COMM_GEOSX );
  
  dofManager.setSparsityPattern( matrix, false ); // don't close the matrix

  // TODO: remove this and just call SolverBase::SetupSystem when DofManager can handle the coupling

  // Populate off-diagonal sparsity between well and reservoir

  string const resDofKey  = dofManager.getKey( m_wellSolver->ResElementDofName() );
  string const wellDofKey = dofManager.getKey( m_wellSolver->WellElementDofName() );

  localIndex const resNDOF = m_wellSolver->NumDofPerResElement();
  localIndex const wellNDOF = m_wellSolver->NumDofPerWellElement();

  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex constexpr maxNumDof  = maxNumComp + 1;

  ElementRegionManager::ElementViewAccessor<arrayView1d<globalIndex const>> const & resDofNumber =
    elemManager->ConstructViewAccessor<array1d<globalIndex>, arrayView1d<globalIndex const>>( resDofKey );

  elemManager->forElementSubRegions<WellElementSubRegion>( [&]( WellElementSubRegion const * const subRegion )
  {
    PerforationData const * const perforationData = subRegion->GetPerforationData();

    // get the well degrees of freedom and ghosting info
    arrayView1d< globalIndex const> const & wellElemDofNumber =
      subRegion->getReference< array1d<globalIndex> >( wellDofKey );

    // get the well element indices corresponding to each perforation
    arrayView1d< localIndex const > const & perfWellElemIndex =
      perforationData->getReference< array1d<localIndex> >( PerforationData::viewKeyStruct::wellElementIndexString );

    // get the element region, subregion, index
    arrayView1d< localIndex const > const & resElementRegion =
      perforationData->getReference< array1d<localIndex> >( PerforationData::viewKeyStruct::reservoirElementRegionString );

    arrayView1d< localIndex const > const & resElementSubRegion =
      perforationData->getReference< array1d<localIndex> >( PerforationData::viewKeyStruct::reservoirElementSubregionString );

    arrayView1d< localIndex const > const & resElementIndex =
      perforationData->getReference< array1d<localIndex> >( PerforationData::viewKeyStruct::reservoirElementIndexString );

    stackArray1d< globalIndex, maxNumDof > dofIndexRes( resNDOF );
    stackArray1d< globalIndex, maxNumDof > dofIndexWell( wellNDOF );
    stackArray2d< real64, maxNumDof * maxNumDof > values( resNDOF, wellNDOF );
    values = 1.0;

    // Insert the entries corresponding to reservoir-well perforations
    // This will fill J_WR, and J_RW
    for( localIndex iperf = 0; iperf < perforationData->size(); ++iperf )
    {
      // get the reservoir (sub)region and element indices
      localIndex const er = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei = resElementIndex[iperf];
      localIndex const iwelem = perfWellElemIndex[iperf];

      for( localIndex idof = 0; idof < resNDOF; ++idof )
      {
        dofIndexRes[idof] = resDofNumber[er][esr][ei] + idof;
      }

      for( localIndex idof = 0; idof < wellNDOF; ++idof )
      {
        dofIndexWell[idof] = wellElemDofNumber[iwelem] + idof;
      }

      // fill J_RW
      matrix.insert( dofIndexRes.data(),
                     dofIndexWell.data(),
                     values.data(),
                     resNDOF,
                     wellNDOF );

      // fill J_WR
      matrix.insert( dofIndexWell.data(),
                     dofIndexRes.data(),
                     values.data(),
                     wellNDOF,
                     resNDOF );
    }

   } );

  matrix.close();
}

void CompositionalMultiphaseReservoir::AssembleCouplingTerms( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                                              real64 const dt,
                                                              DomainPartition * const domain,
                                                              DofManager const * const dofManager,
                                                              ParallelMatrix * const matrix,
                                                              ParallelVector * const rhs )
{
  MeshLevel * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = meshLevel->getElemManager();

  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex constexpr maxNumDof  = maxNumComp + 1;
  
  localIndex const NC      = m_wellSolver->NumFluidComponents();
  localIndex const resNDOF = m_wellSolver->NumDofPerResElement();

  string const resDofKey  = dofManager->getKey( m_wellSolver->ResElementDofName() );
  string const wellDofKey = dofManager->getKey( m_wellSolver->WellElementDofName() );

  ElementRegionManager::ElementViewAccessor< arrayView1d<globalIndex const> > 
  resDofNumberAccessor = elemManager->ConstructViewAccessor< array1d<globalIndex>, arrayView1d<globalIndex const> >( resDofKey );

  ElementRegionManager::ElementViewAccessor< arrayView1d<globalIndex const> >::ViewTypeConst 
  resDofNumber = resDofNumberAccessor.toViewConst();

  elemManager->forElementSubRegions<WellElementSubRegion>( [&]( WellElementSubRegion * const subRegion )
  {

    PerforationData const * const perforationData = subRegion->GetPerforationData();

    // compute the local rates for this well
    ComputeAllPerforationRates( subRegion );
    
    // get the degrees of freedom 
    arrayView1d<globalIndex const> const & 
    wellElemDofNumber = subRegion->getReference<array1d<globalIndex>>( wellDofKey );
    
    // get well variables on perforations
    arrayView2d<real64 const> const & 
    compPerfRate = perforationData->getReference<array2d<real64>>( CompositionalMultiphaseWell::viewKeyStruct::compPerforationRateString );

    arrayView3d<real64 const> const & 
    dCompPerfRate_dPres = perforationData->getReference<array3d<real64>>( CompositionalMultiphaseWell::viewKeyStruct::dCompPerforationRate_dPresString );

    arrayView4d<real64 const> const & 
    dCompPerfRate_dComp = perforationData->getReference<array4d<real64>>( CompositionalMultiphaseWell::viewKeyStruct::dCompPerforationRate_dCompString );
    
    arrayView1d<localIndex const> const & 
    perfWellElemIndex = perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::wellElementIndexString );
    
    // get the element region, subregion, index
    arrayView1d<localIndex const> const & 
    resElementRegion = perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementRegionString );

    arrayView1d<localIndex const> const & 
    resElementSubRegion = perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementSubregionString );

    arrayView1d<localIndex const> const & 
    resElementIndex = perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementIndexString );

    // local working variables and arrays
    stackArray1d<long long, 2 * maxNumComp> eqnRowIndices( 2 * NC );
    stackArray1d<long long, 2 * maxNumDof>  dofColIndices( 2 * resNDOF );

    stackArray1d<double, 2 * maxNumComp>                 localPerf( 2 * NC );
    stackArray2d<double, 2 * maxNumComp * 2 * maxNumDof> localPerfJacobian( 2 * NC, 2 * resNDOF );

    // loop over the perforations and add the rates to the residual and jacobian
    for (localIndex iperf = 0; iperf < perforationData->size(); ++iperf)
    {

      // local working variables and arrays
      eqnRowIndices = -1;
      dofColIndices = -1;

      localPerf = 0;
      localPerfJacobian = 0;

      // get the reservoir (sub)region and element indices
      localIndex const er  = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei  = resElementIndex[iperf];

      // get the well element index for this perforation
      localIndex const iwelem = perfWellElemIndex[iperf];
      globalIndex const resOffset = resDofNumber[er][esr][ei];
      globalIndex const wellElemOffset = wellElemDofNumber[iwelem];

      for (localIndex ic = 0; ic < NC; ++ic)
      {
        eqnRowIndices[WellSolverBase::SubRegionTag::RES  * NC + ic] = resOffset + ic;
        eqnRowIndices[WellSolverBase::SubRegionTag::WELL * NC + ic] =
          wellElemOffset + CompositionalMultiphaseWell::RowOffset::MASSBAL + ic;
      }
      for (localIndex jdof = 0; jdof < resNDOF; ++jdof)
      {
        dofColIndices[WellSolverBase::SubRegionTag::RES  * resNDOF + jdof] = resOffset + jdof;
        dofColIndices[WellSolverBase::SubRegionTag::WELL * resNDOF + jdof] = 
          wellElemOffset + CompositionalMultiphaseWell::ColOffset::DPRES + jdof;
      }
      
      // populate local flux vector and derivatives
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        localPerf[WellSolverBase::SubRegionTag::RES  * NC + ic] =   dt * compPerfRate[iperf][ic];
        localPerf[WellSolverBase::SubRegionTag::WELL * NC + ic] = - dt * compPerfRate[iperf][ic];

        for (localIndex ke = 0; ke < 2; ++ke)
        {
          localIndex const localDofIndexPres = ke * resNDOF;
          localPerfJacobian[WellSolverBase::SubRegionTag::RES  * NC + ic][localDofIndexPres] =   
              dt * dCompPerfRate_dPres[iperf][ke][ic];
          localPerfJacobian[WellSolverBase::SubRegionTag::WELL * NC + ic][localDofIndexPres] = 
            - dt * dCompPerfRate_dPres[iperf][ke][ic];
        
          for (localIndex jc = 0; jc < NC; ++jc)
          {
            localIndex const localDofIndexComp = localDofIndexPres + jc + 1;
            localPerfJacobian[WellSolverBase::SubRegionTag::RES  * NC + ic][localDofIndexComp] =   
                dt * dCompPerfRate_dComp[iperf][ke][ic][jc];
            localPerfJacobian[WellSolverBase::SubRegionTag::WELL * NC + ic][localDofIndexComp] = 
              - dt * dCompPerfRate_dComp[iperf][ke][ic][jc];
          }
        }
      }

      rhs->add( eqnRowIndices.data(),
                localPerf.data(),
                2 * NC );

      matrix->add( eqnRowIndices.data(),
                   dofColIndices.data(),
                   localPerfJacobian.data(),
                   2 * NC, 2 * resNDOF );
    }
  });  
}

void CompositionalMultiphaseReservoir::ComputeAllPerforationRates( WellElementSubRegion * const subRegion ) 
{
  localIndex constexpr maxNumComp = MultiFluidBase::MAX_NUM_COMPONENTS; 
  
  localIndex const NC = m_wellSolver->NumFluidComponents();
  localIndex const NP = m_wellSolver->NumFluidPhases();
  
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const & resPressure  = m_resPressure;
  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const & dResPressure = m_deltaResPressure;
  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> const & resPhaseMob = m_resPhaseMob;
  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> const & dResPhaseMob_dPres = m_dResPhaseMob_dPres;
  ElementRegionManager::ElementViewAccessor<arrayView3d<real64>> const & dResPhaseMob_dComp = m_dResPhaseMob_dCompDens;
  ElementRegionManager::ElementViewAccessor<arrayView2d<real64>> const & dResPhaseVolFrac_dPres = m_dResPhaseVolFrac_dPres;
  ElementRegionManager::ElementViewAccessor<arrayView3d<real64>> const & dResPhaseVolFrac_dComp = m_dResPhaseVolFrac_dCompDens;
  ElementRegionManager::ElementViewAccessor<arrayView3d<real64>> const & dResCompFrac_dCompDens = m_dResCompFrac_dCompDens;

  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & resPhaseVisc        = m_resPhaseVisc;
  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & dResPhaseVisc_dPres = m_dResPhaseVisc_dPres;
  ElementRegionManager::MaterialViewAccessor<arrayView4d<real64>> const & dResPhaseVisc_dComp = m_dResPhaseVisc_dComp;
  ElementRegionManager::MaterialViewAccessor<arrayView4d<real64>> const & resPhaseCompFrac    = m_resPhaseCompFrac;
  ElementRegionManager::MaterialViewAccessor<arrayView4d<real64>> const & dResPhaseCompFrac_dPres = m_dResPhaseCompFrac_dPres;
  ElementRegionManager::MaterialViewAccessor<arrayView5d<real64>> const & dResPhaseCompFrac_dComp = m_dResPhaseCompFrac_dComp;
  ElementRegionManager::MaterialViewAccessor<arrayView3d<real64>> const & resPhaseRelPerm                = m_resPhaseRelPerm;
  ElementRegionManager::MaterialViewAccessor<arrayView4d<real64>> const & dResPhaseRelPerm_dPhaseVolFrac = m_dResPhaseRelPerm_dPhaseVolFrac;

  PerforationData * const perforationData = subRegion->GetPerforationData();

  // get depth
  arrayView1d<real64 const> const & 
  wellElemGravCoef = subRegion->getReference<array1d<real64>>( CompositionalMultiphaseWell::viewKeyStruct::gravityCoefString );
    
  // get well primary variables on well elements
  arrayView1d<real64 const> const & 
  wellElemPressure = subRegion->getReference<array1d<real64>>( CompositionalMultiphaseWell::viewKeyStruct::pressureString );

  arrayView1d<real64 const> const & 
  dWellElemPressure = subRegion->getReference<array1d<real64>>( CompositionalMultiphaseWell::viewKeyStruct::deltaPressureString );

  // get secondary well data on well elements
  arrayView1d<real64 const> const & 
  wellElemMixtureDensity = subRegion->getReference<array1d<real64>>( CompositionalMultiphaseWell::viewKeyStruct::mixtureDensityString );
  
  arrayView1d<real64 const> const & 
  dWellElemMixtureDensity_dPres = subRegion->getReference<array1d<real64>>( CompositionalMultiphaseWell::viewKeyStruct::dMixtureDensity_dPressureString );

  arrayView2d<real64 const> const & 
  dWellElemMixtureDensity_dComp = subRegion->getReference<array2d<real64>>( CompositionalMultiphaseWell::viewKeyStruct::dMixtureDensity_dGlobalCompDensityString );

  arrayView2d<real64 const> const & 
  wellElemCompFrac = subRegion->getReference<array2d<real64>>( CompositionalMultiphaseWell::viewKeyStruct::globalCompFractionString );

  arrayView3d<real64 const> const & 
  dWellElemCompFrac_dCompDens = subRegion->getReference<array3d<real64>>( CompositionalMultiphaseWell::viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

  // get well variables on perforations
  arrayView1d<real64 const> const & 
  perfGravCoef = perforationData->getReference<array1d<real64>>( CompositionalMultiphaseWell::viewKeyStruct::gravityCoefString );

  arrayView1d<localIndex const> const & 
  perfWellElemIndex = perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::wellElementIndexString );

  arrayView1d<real64 const> const & 
  perfTransmissibility = perforationData->getReference<array1d<real64>>( PerforationData::viewKeyStruct::transmissibilityString );

  arrayView2d<real64> const & 
  compPerfRate = perforationData->getReference<array2d<real64>>( CompositionalMultiphaseWell::viewKeyStruct::compPerforationRateString );

  arrayView3d<real64> const & 
  dCompPerfRate_dPres = perforationData->getReference<array3d<real64>>( CompositionalMultiphaseWell::viewKeyStruct::dCompPerforationRate_dPresString );

  arrayView4d<real64> const & 
  dCompPerfRate_dComp = perforationData->getReference<array4d<real64>>( CompositionalMultiphaseWell::viewKeyStruct::dCompPerforationRate_dCompString );

  // get the element region, subregion, index
  arrayView1d<localIndex const> const & 
  resElementRegion = perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementRegionString );

  arrayView1d<localIndex const> const & 
  resElementSubRegion = perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementSubregionString );

  arrayView1d<localIndex const> const & 
  resElementIndex = perforationData->getReference<array1d<localIndex>>( PerforationData::viewKeyStruct::reservoirElementIndexString );

  localIndex const & 
  resFluidIndex = m_wellSolver->getReference<localIndex>( WellSolverBase::viewKeyStruct::resFluidIndexString );


  // local working variables and arrays
  stackArray1d<real64, maxNumComp> dPhaseCompFrac_dCompDens( NC );

  stackArray1d<real64, 2>              pressure( 2 );
  stackArray1d<real64, 2>              dPressure_dP( 2 );
  stackArray2d<real64, 2 * maxNumComp> dPressure_dC( 2, NC );
      
  stackArray1d<real64, 2>              dFlux_dP( 2 );
  stackArray2d<real64, 2 * maxNumComp> dFlux_dC( 2, NC );

  stackArray1d<real64, 2>              dMult_dP( 2 );
  stackArray2d<real64, 2 * maxNumComp> dMult_dC( 2, NC );

  stackArray1d<real64, maxNumComp> dMixtureDensity_dC( NC );
  stackArray1d<real64, maxNumComp> dResTotalMobility_dC( NC );

  stackArray2d<real64, 2 * maxNumComp>              phaseCompFrac( 2, NC );
  stackArray2d<real64, 2 * maxNumComp>              dPhaseCompFrac_dP( 2, NC );
  stackArray3d<real64, 2 * maxNumComp * maxNumComp> dPhaseCompFrac_dC( 2, NC, NC );
    
  stackArray1d<real64, maxNumComp> dVisc_dC( NC );
  stackArray1d<real64, maxNumComp> dRelPerm_dC( NC );

  stackArray1d<real64, 2>              dPotDiff_dP( 2 );
  stackArray2d<real64, 2 * maxNumComp> dPotDiff_dC( 2, NC );

  stackArray1d<real64, 2> multiplier( 2 );


  // loop over the perforations to compute the perforation rates
  for (localIndex iperf = 0; iperf < perforationData->size(); ++iperf)
  {     

    // reset the perforation rates
    for (localIndex ic = 0; ic < NC; ++ic)
    {
      compPerfRate[iperf][ic] = 0.0;
      for (localIndex ke = 0; ke < 2; ++ke) 
      {
        dCompPerfRate_dPres[iperf][ke][ic] = 0.0;
        for (localIndex jc = 0; jc < NC; ++jc) 
        {
          dCompPerfRate_dComp[iperf][ke][ic][jc] = 0.0;
        }
      }
    }

    // clear working arrays
    pressure = 0.0;
    dPressure_dP = 0.0;
    dPressure_dC = 0.0;
        
    // 1) copy the variables from the reservoir and well element
     
    // a) get reservoir variables

    // get the reservoir (sub)region and element indices
    localIndex const er  = resElementRegion[iperf];
    localIndex const esr = resElementSubRegion[iperf];
    localIndex const ei  = resElementIndex[iperf];
    // get the index of the well elem
    localIndex const iwelem = perfWellElemIndex[iperf]; 

    pressure[CompositionalMultiphaseWell::SubRegionTag::RES] = 
      resPressure[er][esr][ei] + dResPressure[er][esr][ei];
    dPressure_dP[CompositionalMultiphaseWell::SubRegionTag::RES] = 1.0;

    // TODO: add a buoyancy term for the reservoir side here 

    multiplier[CompositionalMultiphaseWell::SubRegionTag::RES] = 1.0;

    // b) get well variables

    pressure[CompositionalMultiphaseWell::SubRegionTag::WELL] = 
      wellElemPressure[iwelem] + dWellElemPressure[iwelem];
    dPressure_dP[CompositionalMultiphaseWell::SubRegionTag::WELL] = 1.0;

    multiplier[CompositionalMultiphaseWell::SubRegionTag::WELL] = -1.0;

    real64 const gravD = ( perfGravCoef[iperf] - wellElemGravCoef[iwelem] );

    pressure[CompositionalMultiphaseWell::SubRegionTag::WELL] += 
      wellElemMixtureDensity[iwelem] * gravD;
    dPressure_dP[CompositionalMultiphaseWell::SubRegionTag::WELL] += 
      dWellElemMixtureDensity_dPres[iwelem] * gravD;
    for (localIndex ic = 0; ic < NC; ++ic)
    {
      dPressure_dC[CompositionalMultiphaseWell::SubRegionTag::WELL][ic] += 
        dWellElemMixtureDensity_dComp[iwelem][ic] * gravD;
    }

    // get transmissibility at the interface
    real64 const trans = perfTransmissibility[iperf]; 
        
    // 2) compute potential difference

    real64 potDiff = 0.0;
    dPotDiff_dP = 0.0;
    dPotDiff_dC = 0.0;

    for (localIndex i = 0; i < 2; ++i)
    {
      potDiff += multiplier[i] * trans * pressure[i]; // pressure = pres + dPres
      dPotDiff_dP[i] += multiplier[i] * trans * dPressure_dP[i];

      for (localIndex ic = 0; ic < NC; ++ic) 
      {
        dPotDiff_dC[i][ic] += multiplier[i] * trans * dPressure_dC[i][ic];
      }
    }

    real64 flux = 0.0;
    dFlux_dP = 0.0;
    dFlux_dC = 0.0;

    // 3) upwinding 

    if ( potDiff >= 0 ) // ** reservoir cell is upstream **
    {

      phaseCompFrac = 0.0;
      dPhaseCompFrac_dP = 0.0;
      dPhaseCompFrac_dC = 0.0;

      dPhaseCompFrac_dCompDens = 0;

      // loop over phases, compute and upwind phase flux
      // and sum contributions to each component's perforation rate
      for (localIndex ip = 0; ip < NP; ++ip)
      {
      
        // compute the phase flux and derivatives using upstream cell mobility
        flux = resPhaseMob[er][esr][ei][ip] * potDiff;

        dFlux_dP[CompositionalMultiphaseWell::SubRegionTag::RES] = 
            dResPhaseMob_dPres[er][esr][ei][ip] * potDiff 
          + resPhaseMob[er][esr][ei][ip]        * dPotDiff_dP[CompositionalMultiphaseWell::SubRegionTag::RES];

        dFlux_dP[CompositionalMultiphaseWell::SubRegionTag::WELL] = 
            resPhaseMob[er][esr][ei][ip]        * dPotDiff_dP[CompositionalMultiphaseWell::SubRegionTag::WELL];

        for (localIndex ic = 0; ic < NC; ++ic)
        {
          dFlux_dC[CompositionalMultiphaseWell::SubRegionTag::RES][ic] = 
              dResPhaseMob_dComp[er][esr][ei][ip][ic] * potDiff
            + resPhaseMob[er][esr][ei][ip]            * dPotDiff_dC[CompositionalMultiphaseWell::SubRegionTag::RES][ic];

          dFlux_dC[CompositionalMultiphaseWell::SubRegionTag::WELL][ic] = 
              resPhaseMob[er][esr][ei][ip]            * dPotDiff_dC[CompositionalMultiphaseWell::SubRegionTag::RES][ic];
        }
        
        // increment component fluxes
        for (localIndex ic = 0; ic < NC; ++ic)
        {
          compPerfRate[iperf][ic] += flux * resPhaseCompFrac[er][esr][resFluidIndex][ei][0][ip][ic];
        
          dCompPerfRate_dPres[iperf][CompositionalMultiphaseWell::SubRegionTag::RES][ic] += 
            resPhaseCompFrac[er][esr][resFluidIndex][ei][0][ip][ic]
            * dFlux_dP[CompositionalMultiphaseWell::SubRegionTag::RES];

          dCompPerfRate_dPres[iperf][CompositionalMultiphaseWell::SubRegionTag::RES][ic] += 
            dResPhaseCompFrac_dPres[er][esr][resFluidIndex][ei][0][ip][ic] * flux;

          dCompPerfRate_dPres[iperf][CompositionalMultiphaseWell::SubRegionTag::WELL][ic] += 
            resPhaseCompFrac[er][esr][resFluidIndex][ei][0][ip][ic]
            * dFlux_dP[CompositionalMultiphaseWell::SubRegionTag::WELL];

          applyChainRule( NC, 
                          dResCompFrac_dCompDens[er][esr][ei], 
                          dResPhaseCompFrac_dComp[er][esr][resFluidIndex][ei][0][ip][ic], 
                          dPhaseCompFrac_dCompDens );
           
          for (localIndex jc = 0; jc < NC; ++jc)
          {
            dCompPerfRate_dComp[iperf][CompositionalMultiphaseWell::SubRegionTag::RES][ic][jc] += 
              dFlux_dC[CompositionalMultiphaseWell::SubRegionTag::RES][jc] 
              * resPhaseCompFrac[er][esr][resFluidIndex][ei][0][ip][ic];

            dCompPerfRate_dComp[iperf][CompositionalMultiphaseWell::SubRegionTag::RES][ic][jc] += 
              flux * dPhaseCompFrac_dCompDens[jc];

            dCompPerfRate_dComp[iperf][CompositionalMultiphaseWell::SubRegionTag::WELL][ic][jc] += 
              dFlux_dC[CompositionalMultiphaseWell::SubRegionTag::WELL][jc] 
              * resPhaseCompFrac[er][esr][resFluidIndex][ei][0][ip][ic];
          }
        }
      }
    }
    else // ** well is upstream **
    {

      real64 resTotalMobility     = 0.0;
      real64 dResTotalMobility_dP = 0.0;
      dResTotalMobility_dC        = 0.0;

      // first, compute the reservoir total mobitity (excluding phase density)    
      for (localIndex ip = 0; ip < NP; ++ip)
      {
        // viscosity
        real64 const resViscosity = resPhaseVisc[er][esr][resFluidIndex][ei][0][ip];
        real64 const dResVisc_dP  = dResPhaseVisc_dPres[er][esr][resFluidIndex][ei][0][ip];
        applyChainRule( NC, dResCompFrac_dCompDens[er][esr][ei], 
                        dResPhaseVisc_dComp[er][esr][resFluidIndex][ei][0][ip], 
                        dVisc_dC );
      
        // relative permeability
        localIndex const & resRelPermIndex = m_wellSolver->getReference<localIndex>( CompositionalMultiphaseWell::viewKeyStruct::resRelPermIndexString );
        real64 const resRelPerm = resPhaseRelPerm[er][esr][resRelPermIndex][ei][0][ip];
        real64 dResRelPerm_dP = 0.0;
        dRelPerm_dC = 0.0;
        for (localIndex jp = 0; jp < NP; ++jp)
        {
          real64 const dResRelPerm_dS = dResPhaseRelPerm_dPhaseVolFrac[er][esr][resRelPermIndex][ei][0][ip][jp];
          dResRelPerm_dP += dResRelPerm_dS * dResPhaseVolFrac_dPres[er][esr][ei][jp];

          for (localIndex jc = 0; jc < NC; ++jc)
          {
            dRelPerm_dC[jc] += dResRelPerm_dS * dResPhaseVolFrac_dComp[er][esr][ei][jp][jc];
          }
        }

        // increment total mobility
        resTotalMobility     += resRelPerm / resViscosity;
        dResTotalMobility_dP += ( dResRelPerm_dP * resViscosity - resRelPerm * dResVisc_dP )
                              / ( resViscosity * resViscosity );
        for (localIndex ic = 0; ic < NC; ++ic)
        {
          dResTotalMobility_dC[ic] += ( dRelPerm_dC[ic] * resViscosity - resRelPerm * dVisc_dC[ic] )
            / ( resViscosity * resViscosity );
        }
      }

      // compute a potdiff multiplier = wellElemMixtureDensity * resTotalMobility 
      real64 const mult = wellElemMixtureDensity[iwelem] * resTotalMobility;
      dMult_dP[CompositionalMultiphaseWell::SubRegionTag::RES]  = 
        wellElemMixtureDensity[iwelem] * dResTotalMobility_dP;
      dMult_dP[CompositionalMultiphaseWell::SubRegionTag::WELL] = 
        dWellElemMixtureDensity_dPres[iwelem] * resTotalMobility;

      for (localIndex ic = 0; ic < NC; ++ic)
      {
        dMult_dC[CompositionalMultiphaseWell::SubRegionTag::RES][ic]  = 
          wellElemMixtureDensity[iwelem] * dResTotalMobility_dC[ic];

        dMult_dC[CompositionalMultiphaseWell::SubRegionTag::WELL][ic] = 
          dWellElemMixtureDensity_dComp[iwelem][ic] * resTotalMobility;
      }

      // compute the volumetric flux and derivatives using upstream cell mobility
      flux = mult * potDiff;

      dFlux_dP[CompositionalMultiphaseWell::SubRegionTag::RES] =  
        dMult_dP[CompositionalMultiphaseWell::SubRegionTag::RES] * potDiff 
        + mult * dPotDiff_dP[CompositionalMultiphaseWell::SubRegionTag::RES];

      dFlux_dP[CompositionalMultiphaseWell::SubRegionTag::WELL] = 
        dMult_dP[CompositionalMultiphaseWell::SubRegionTag::WELL] * potDiff 
        + mult * dPotDiff_dP[CompositionalMultiphaseWell::SubRegionTag::WELL];

      for (localIndex ic = 0; ic < NC; ++ic)
      {
        dFlux_dC[CompositionalMultiphaseWell::SubRegionTag::RES][ic] = 
          dMult_dC[CompositionalMultiphaseWell::SubRegionTag::RES][ic] * potDiff 
          + mult * dPotDiff_dC[CompositionalMultiphaseWell::SubRegionTag::RES][ic];

        dFlux_dC[CompositionalMultiphaseWell::SubRegionTag::WELL][ic] = 
          dMult_dC[CompositionalMultiphaseWell::SubRegionTag::WELL][ic] * potDiff 
          + mult * dPotDiff_dC[CompositionalMultiphaseWell::SubRegionTag::RES][ic];
      }
      
      // compute component fluxes
      for (localIndex ic = 0; ic < NC; ++ic)
      {
        compPerfRate[iperf][ic] += wellElemCompFrac[iwelem][ic] * flux;
        
        dCompPerfRate_dPres[iperf][CompositionalMultiphaseWell::SubRegionTag::RES][ic] = 
          wellElemCompFrac[iwelem][ic] * dFlux_dP[CompositionalMultiphaseWell::SubRegionTag::RES];

        dCompPerfRate_dPres[iperf][CompositionalMultiphaseWell::SubRegionTag::WELL][ic] = 
          wellElemCompFrac[iwelem][ic] * dFlux_dP[CompositionalMultiphaseWell::SubRegionTag::WELL];
           
        for (localIndex jc = 0; jc < NC; ++jc)
        {
          dCompPerfRate_dComp[iperf][CompositionalMultiphaseWell::SubRegionTag::RES][ic][jc]  += 
            wellElemCompFrac[iwelem][ic] * dFlux_dC[CompositionalMultiphaseWell::SubRegionTag::RES][jc]; 

          dCompPerfRate_dComp[iperf][CompositionalMultiphaseWell::SubRegionTag::WELL][ic][jc] += 
            wellElemCompFrac[iwelem][ic] * dFlux_dC[CompositionalMultiphaseWell::SubRegionTag::WELL][jc]; 

          dCompPerfRate_dComp[iperf][CompositionalMultiphaseWell::SubRegionTag::WELL][ic][jc] += 
            dWellElemCompFrac_dCompDens[iwelem][ic][jc] * flux;
        }
      }
    }      
  }
}


void CompositionalMultiphaseReservoir::ResetViews( DomainPartition * const domain )
{
  ReservoirSolverBase::ResetViews(domain);

  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = mesh->getElemManager();
  ConstitutiveManager * const constitutiveManager = domain->getConstitutiveManager();

  m_resPressure =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::pressureString );

  m_deltaResPressure =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::deltaPressureString );

  m_resPhaseMob =
    elemManager->ConstructViewAccessor<array2d<real64>, arrayView2d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::phaseMobilityString );

  m_dResPhaseMob_dPres =
    elemManager->ConstructViewAccessor<array2d<real64>, arrayView2d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::dPhaseMobility_dPressureString );

  m_dResPhaseMob_dCompDens =
    elemManager->ConstructViewAccessor<array3d<real64>, arrayView3d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::dPhaseMobility_dGlobalCompDensityString );

  m_dResPhaseVolFrac_dPres =
    elemManager->ConstructViewAccessor<array2d<real64>, arrayView2d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::dPhaseVolumeFraction_dPressureString );

  m_dResPhaseVolFrac_dCompDens =
    elemManager->ConstructViewAccessor<array3d<real64>, arrayView3d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );

  m_dResCompFrac_dCompDens =
    elemManager->ConstructViewAccessor<array3d<real64>, arrayView3d<real64>>( CompositionalMultiphaseFlow::viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );



  m_resPhaseVisc =
    elemManager->ConstructFullMaterialViewAccessor<array3d<real64>, arrayView3d<real64>>( MultiFluidBase::viewKeyStruct::phaseViscosityString,
                                                                                          constitutiveManager );
  m_dResPhaseVisc_dPres =
    elemManager->ConstructFullMaterialViewAccessor<array3d<real64>, arrayView3d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseViscosity_dPressureString,
                                                                                          constitutiveManager );
  m_dResPhaseVisc_dComp =
    elemManager->ConstructFullMaterialViewAccessor<array4d<real64>, arrayView4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseViscosity_dGlobalCompFractionString,
                                                                                          constitutiveManager );

  m_resPhaseCompFrac =
    elemManager->ConstructFullMaterialViewAccessor<array4d<real64>, arrayView4d<real64>>( MultiFluidBase::viewKeyStruct::phaseCompFractionString,
                                                                                          constitutiveManager );
  m_dResPhaseCompFrac_dPres =
    elemManager->ConstructFullMaterialViewAccessor<array4d<real64>, arrayView4d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseCompFraction_dPressureString,
                                                                                          constitutiveManager );
  m_dResPhaseCompFrac_dComp =
    elemManager->ConstructFullMaterialViewAccessor<array5d<real64>, arrayView5d<real64>>( MultiFluidBase::viewKeyStruct::dPhaseCompFraction_dGlobalCompFractionString,
                                                                                          constitutiveManager );

  m_resPhaseRelPerm =
    elemManager->ConstructFullMaterialViewAccessor<array3d<real64>, arrayView3d<real64>>( RelativePermeabilityBase::viewKeyStruct::phaseRelPermString,
                                                                                          constitutiveManager );

  m_dResPhaseRelPerm_dPhaseVolFrac =
    elemManager->ConstructFullMaterialViewAccessor<array4d<real64>, arrayView4d<real64>>( RelativePermeabilityBase::viewKeyStruct::dPhaseRelPerm_dPhaseVolFractionString,
                                                                                          constitutiveManager );

}

REGISTER_CATALOG_ENTRY( SolverBase, CompositionalMultiphaseReservoir, std::string const &, Group * const )

} /* namespace geosx */

