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
#include "mpiCommunications/CommunicationTools.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"

#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"


/**
 * @namespace the geosx namespace that encapsulates the majority of the code
 */
namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

GeochemicalModel::GeochemicalModel( const std::string& name,
                                  Group * const parent ):
  FlowSolverBase(name, parent)
{

  this->registerWrapper( viewKeyStruct::reactiveFluidNameString,  &m_reactiveFluidName,  false )->setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of chemical system constitutive object to use for this solver");

  this->registerWrapper( viewKeyStruct::reactiveFluidIndexString, &m_reactiveFluidIndex, false );

  this->registerWrapper( viewKeyStruct::outputSpeciesFileNameString, &m_outputSpeciesFileName, false )->
    setInputFlag(InputFlags::OPTIONAL)->    
    setDescription("Output species to file");

}

void GeochemicalModel::RegisterDataOnMesh(Group * const MeshBodies)
{
  FlowSolverBase::RegisterDataOnMesh(MeshBodies);

  for( auto & mesh : MeshBodies->GetSubGroups() )
  {

    MeshLevel * const meshLevel = Group::group_cast<MeshBody *>(mesh.second)->getMeshLevel(0);

    applyToSubRegions( meshLevel, [&] ( ElementSubRegionBase * const subRegion)
    {
      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::pressureString )->setPlotLevel(PlotLevel::LEVEL_0);
      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::deltaPressureString );

      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::temperatureString )->setPlotLevel(PlotLevel::LEVEL_0);
      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::deltaTemperatureString );      

      subRegion->registerWrapper< array2d<real64> >( viewKeyStruct::concentrationString )->setPlotLevel(PlotLevel::LEVEL_0);
      subRegion->registerWrapper< array2d<real64> >( viewKeyStruct::deltaConcentrationString );

      subRegion->registerWrapper< array2d<real64> >( viewKeyStruct::totalConcentrationString );
      subRegion->registerWrapper< array2d<real64> >( viewKeyStruct::concentrationNewString );                  

      subRegion->registerWrapper< array1d<globalIndex> >( viewKeyStruct::blockLocalDofNumberString );

    } );

  }
}

void GeochemicalModel::InitializePreSubGroups(Group * const rootGroup)
{
  FlowSolverBase::InitializePreSubGroups(rootGroup);

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
    MeshLevel * meshLevel = Group::group_cast<MeshBody *>(mesh.second)->getMeshLevel(0);
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
  
void GeochemicalModel::UpdateReactiveFluidModel(Group * const dataGroup)
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

void GeochemicalModel::UpdateState( Group * dataGroup )
{
  GEOSX_MARK_FUNCTION;

  UpdateReactiveFluidModel( dataGroup );

}

void GeochemicalModel::InitializePostInitialConditions_PreSubGroups( Group * const rootGroup )
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

  real64 dt_return = dt;

  ImplicitStepSetup( time_n,
                     dt,
                     domain,
                     m_dofManager,
                     m_matrix,
                     m_rhs,
                     m_solution );


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

  return dt_return;

}


void GeochemicalModel::ImplicitStepSetup( real64 const & GEOSX_UNUSED_ARG( time_n ),
                                          real64 const & GEOSX_UNUSED_ARG( dt ),
                                           DomainPartition * const domain,
                                           DofManager & GEOSX_UNUSED_ARG(dofManager),
                                           ParallelMatrix & GEOSX_UNUSED_ARG(matrix),
                                           ParallelVector & GEOSX_UNUSED_ARG(rhs),
                                          ParallelVector & GEOSX_UNUSED_ARG(solution) )
{
  ResetViews( domain );

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  /* The loop below could be moved to SolverStep after ImplicitStepSetup */
  
  applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegionBase * const GEOSX_UNUSED_ARG( region ),
                                 ElementSubRegionBase * const subRegion )
  {

    arrayView1d<real64> const & dPres   = m_deltaPressure[er][esr];
    arrayView1d<real64> const & dTemp   = m_deltaTemperature[er][esr];
    arrayView2d<real64> const & dConc   = m_deltaConcentration[er][esr];
    arrayView2d<real64> const & conc   = m_concentration[er][esr];
    arrayView2d<real64> const & totalConc   = m_totalConcentration[er][esr];
    

    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
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
  SetupSystem( domain,
               m_dofManager,
               m_matrix,
               m_rhs,
               m_solution  );

}

void GeochemicalModel::ImplicitStepComplete(real64 const & GEOSX_UNUSED_ARG( time_n ),
                                            real64 const & GEOSX_UNUSED_ARG( dt ),
                                            DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegionBase * const GEOSX_UNUSED_ARG( region ),
                                 ElementSubRegionBase * const subRegion )
  {
    arrayView1d<real64> const & pres = m_pressure[er][esr];
    arrayView1d<real64> const & temp = m_temperature[er][esr];    

    arrayView2d<real64> const & conc = m_concentration[er][esr];
    
    arrayView1d<real64 const> const & dPres = m_deltaPressure[er][esr];
    arrayView1d<real64 const> const & dTemp = m_deltaTemperature[er][esr];    
    arrayView2d<real64 const> const & dConc = m_deltaConcentration[er][esr];    

    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {

      pres[ei] += dPres[ei];
      temp[ei] += dTemp[ei];      
      for (localIndex ic = 0; ic < m_numBasisSpecies; ++ic)
        conc[ei][ic] += dConc[ei][ic];

    } );

  } );

  WriteSpeciesToFile(domain);

}

void GeochemicalModel::SetupDofs( DomainPartition const * const GEOSX_UNUSED_ARG(domain),
                                   DofManager & dofManager ) const
{
  dofManager.addField( viewKeyStruct::pressureString,
                       DofManager::Location::Elem,
                       DofManager::Connectivity::Face,
                       m_numDofPerCell,
                       m_targetRegions );
}


void GeochemicalModel::AssembleSystem( real64 const time,
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

  AssembleFluxTerms( time,
                     dt,
                     domain,
                     &dofManager,
                     &matrix,
                     &rhs);

    
  matrix.close();
  rhs.close();


  if( verboseLevel() >= 2 )
  {
    GEOS_LOG_RANK("After GeochemicalModel::AssembleSystem");
    GEOS_LOG_RANK_0("\nJacobian:\n");
    matrix.print(std::cout);
    GEOS_LOG_RANK_0("\nResidual:\n");
    rhs.print(std::cout);
  }

}

void GeochemicalModel::AssembleAccumulationTerms( DomainPartition * const domain,
                                                   DofManager const * const dofManager,
                                                   ParallelMatrix * const matrix,
                                                   ParallelVector * const rhs)
{
  GEOSX_MARK_FUNCTION;

  ConstitutiveManager * const cm = domain->getConstitutiveManager();

  ReactiveFluidBase const * reactiveFluid  = cm->GetConstitutiveRelation<ReactiveFluidBase>( m_reactiveFluidName );
  
  const array2d<real64> & stochMatrix = reactiveFluid->StochMatrix();
  const array1d<bool> & isHplus = reactiveFluid->IsHplus();   

  MeshLevel const * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  
  applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegionBase const * const GEOSX_UNUSED_ARG( region ),
                                 ElementSubRegionBase const * const subRegion )
  {

    string const dofKey = dofManager->getKey( viewKeyStruct::pressureString );
    arrayView1d<globalIndex const> const & dofNumber = subRegion->getReference< array1d<globalIndex> >( dofKey );    
    
    arrayView1d<integer const>     const & elemGhostRank = m_elemGhostRank[er][esr];

    arrayView2d<real64 const> const & conc          = m_concentration[er][esr];
    arrayView2d<real64 const> const & dConc         = m_deltaConcentration[er][esr];    
    arrayView2d<real64 const> const & totalConc     = m_totalConcentration[er][esr];

    arrayView2d<real64 const> const & dependentConc     = m_dependentConc[er][esr][m_reactiveFluidIndex];
    arrayView3d<real64 const> const & dDependentConc_dConc     = m_dDependentConc_dConc[er][esr][m_reactiveFluidIndex];         

    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      if (elemGhostRank[ei] < 0)
      {
        stackArray1d<globalIndex, ReactiveFluidBase::MAX_NUM_SPECIES> localAccumDOF( m_numDofPerCell );
        stackArray1d<real64, ReactiveFluidBase::MAX_NUM_SPECIES> localAccum( m_numDofPerCell );
        stackArray2d<real64, ReactiveFluidBase::MAX_NUM_SPECIES * ReactiveFluidBase::MAX_NUM_SPECIES> localAccumJacobian( m_numDofPerCell, m_numDofPerCell );
        
        globalIndex const elemDOF = dofNumber[ei];
        
        for (localIndex idof = 0; idof < m_numDofPerCell; ++idof)
        {
          localAccumDOF[idof] = elemDOF + idof;
        }
        
        localAccumJacobian = 0.0;

        for (localIndex ic = 0; ic < m_numBasisSpecies; ++ic)
          {

            real64 concBasis = pow(10.0, conc[ei][ic]+dConc[ei][ic]);

            localAccum[ic] = concBasis - totalConc[ei][ic];

            localAccumJacobian[ic][ic] = log(10.0) * concBasis;

            if(isHplus[ic])
              {
                localAccum[ic] = 0.0;
                localAccumJacobian[ic][ic] = 1.0;
                continue;
              }

            for (localIndex id = 0; id < m_numDependentSpecies; ++id)
              {
                real64 concDependent = pow(10.0, dependentConc[ei][id]);
                localAccum[ic] -= stochMatrix[ic][id] * concDependent;

                for (localIndex idc = 0; idc < m_numBasisSpecies; ++idc)
                  {

                    localAccumJacobian[ic][idc] -= stochMatrix[ic][id] * log(10.0) * concDependent * dDependentConc_dConc[ei][id][idc];

                  }
              }
          }



        // add contribution to global residual and jacobian
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

void GeochemicalModel::AssembleFluxTerms( real64 const GEOSX_UNUSED_ARG(time_n),
                                          real64 const GEOSX_UNUSED_ARG(dt),
                                          DomainPartition const * const GEOSX_UNUSED_ARG(domain),
                                          DofManager const * const GEOSX_UNUSED_ARG(dofManager),
                                          ParallelMatrix * const GEOSX_UNUSED_ARG(matrix),
                                          ParallelVector * const GEOSX_UNUSED_ARG(rhs) )
{

}

void GeochemicalModel::ApplyBoundaryConditions(real64 const GEOSX_UNUSED_ARG(time_n),
                                               real64 const GEOSX_UNUSED_ARG(dt),
                                               DomainPartition * const GEOSX_UNUSED_ARG(domain),
                                               DofManager const & GEOSX_UNUSED_ARG(dofManager),
                                               ParallelMatrix & GEOSX_UNUSED_ARG(matrix),
                                               ParallelVector & GEOSX_UNUSED_ARG(rhs) )
{

}


real64
GeochemicalModel::
CalculateResidualNorm( DomainPartition const * const domain,
                       DofManager const & dofManager,
                       ParallelVector const & rhs )
{

  MeshLevel const * const mesh = domain->getMeshBody(0)->getMeshLevel(0);

  // get a view into local residual vector
  real64 const * localResidual = rhs.extractLocalVector();

  string const dofKey = dofManager.getKey( viewKeyStruct::pressureString );  

  // compute the norm of local residual scaled by cell pore volume
  real64 localResidualNorm = 0.0;

  applyToSubRegions( mesh, [&] ( localIndex const er, localIndex const esr,
                                 ElementRegionBase const * const GEOSX_UNUSED_ARG( region ),
                                 ElementSubRegionBase const * const subRegion )
  {

    arrayView1d<globalIndex const> const & dofNumber = subRegion->getReference< array1d<globalIndex> >( dofKey );
    
    arrayView1d<integer const> const & elemGhostRank = m_elemGhostRank[er][esr];

    localIndex const subRegionSize = subRegion->size();
    for ( localIndex a = 0; a < subRegionSize; ++a )
    {

      if (elemGhostRank[a] < 0)
      {
        localIndex const offset = dofNumber[a];
        for (localIndex idof = 0; idof < m_numDofPerCell; ++idof)
        {
          localIndex const lid = rhs.getLocalRowID( offset + idof );
          real64 const val = localResidual[lid];
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


void GeochemicalModel::ApplySystemSolution( DofManager const & dofManager,
                                             ParallelVector const & solution,
                                             real64 const scalingFactor,
                                             DomainPartition * const domain )
{
  
  MeshLevel * mesh = domain->getMeshBody(0)->getMeshLevel(0);

  applyToSubRegions( mesh, [&] ( localIndex GEOSX_UNUSED_ARG(er), localIndex GEOSX_UNUSED_ARG(esr),
                                 ElementRegionBase * const GEOSX_UNUSED_ARG( region ),
                                 ElementSubRegionBase * const subRegion )
  {

    dofManager.addVectorToField( solution,
                                 viewKeyStruct::pressureString,
                                 scalingFactor,
                                 subRegion,
                                 viewKeyStruct::deltaConcentrationString,
                                 0, m_numDofPerCell );
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

void GeochemicalModel::SolveSystem( DofManager const & dofManager,
                                     ParallelMatrix & matrix,
                                     ParallelVector & rhs,
                                     ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  rhs.scale( -1.0 );
  solution.zero();

  SolverBase::SolveSystem( dofManager, matrix, rhs, solution );

  if( verboseLevel() >= 2 )
  {
    GEOS_LOG_RANK("After GeochemicalModel::SolveSystem");
    GEOS_LOG_RANK("\nsolution\n" << solution);
  }

}

void GeochemicalModel::ResetStateToBeginningOfStep( DomainPartition * const GEOSX_UNUSED_ARG(domain) )
{
}

void GeochemicalModel::ResetViews(DomainPartition * const domain)
{
  FlowSolverBase::ResetViews(domain);

  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = mesh->getElemManager();
  ConstitutiveManager * const constitutiveManager = domain->getConstitutiveManager();

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
      
      ElementRegionBase const * const elemRegion = elementManager->GetRegion(er);
      
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
            
                      if(fabs(spLogC) < 1e-64 || spC < 1e-40)
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

REGISTER_CATALOG_ENTRY( SolverBase, GeochemicalModel, std::string const &, Group * const )
} /* namespace geosx */
