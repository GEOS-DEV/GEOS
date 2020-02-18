/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SinglePhaseBase.cpp
 */

#include "SinglePhaseBase.hpp"

#include "mpiCommunications/CommunicationTools.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseKernels.hpp"

/**
 * @namespace the geosx namespace that encapsulates the majority of the code
 */
namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace SinglePhaseKernels;

SinglePhaseBase::SinglePhaseBase( const std::string& name,
                                  Group * const parent ):
  FlowSolverBase(name, parent)
{
  m_numDofPerCell = 1;
}


void SinglePhaseBase::RegisterDataOnMesh(Group * const MeshBodies)
{
  FlowSolverBase::RegisterDataOnMesh(MeshBodies);

  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    MeshLevel * meshLevel = Group::group_cast<MeshBody *>(mesh.second)->getMeshLevel(0);

    ElementRegionManager * const elemManager = meshLevel->getElemManager();

    elemManager->forElementSubRegions<CellElementSubRegion>( [&]( CellElementSubRegion * const subRegion )
    {
      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::pressureString )->setPlotLevel(PlotLevel::LEVEL_0);
      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::deltaPressureString );
      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::deltaVolumeString );
      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::mobilityString );
      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::dMobility_dPressureString );
      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::porosityString )->setPlotLevel(PlotLevel::LEVEL_1);
      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::porosityOldString );
      subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::densityOldString );
    });

    elemManager->forElementRegions<FaceElementRegion>( [&] ( FaceElementRegion * const region )
    {
      region->forElementSubRegions<FaceElementSubRegion>( [&]( FaceElementSubRegion * const subRegion )
      {
        subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::pressureString )->setPlotLevel(PlotLevel::LEVEL_0);
        subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::deltaPressureString );
        subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::deltaVolumeString );
        subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::mobilityString );
        subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::dMobility_dPressureString );
        subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::porosityString )->
          setDefaultValue(1.0);
        subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::porosityOldString )->
          setDefaultValue(1.0);
        subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::densityOldString );
      });
    });

    // TODO restrict this to boundary sets
    FaceManager * const faceManager = meshLevel->getFaceManager();
    {
      faceManager->registerWrapper<array1d<real64> >( viewKeyStruct::boundaryFacePressureString );
      faceManager->registerWrapper<array2d<real64> >( viewKeyStruct::boundaryFaceDensityString )->reference().resizeDimension<1>(1);
      faceManager->registerWrapper<array2d<real64> >( viewKeyStruct::boundaryFaceViscosityString )->reference().resizeDimension<1>(1);
      faceManager->registerWrapper<array1d<real64> >( viewKeyStruct::boundaryFaceMobilityString );
    }
  }
}

void SinglePhaseBase::UpdateFluidModel(Group * const dataGroup) const
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

void SinglePhaseBase::UpdateSolidModel(Group * const dataGroup) const
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

void SinglePhaseBase::UpdateMobility( Group * const dataGroup ) const
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


void SinglePhaseBase::UpdateState( Group * dataGroup ) const
{
  GEOSX_MARK_FUNCTION;

  UpdateFluidModel( dataGroup );
  UpdateSolidModel( dataGroup );
  UpdateMobility( dataGroup );
}

void SinglePhaseBase::InitializePostInitialConditions_PreSubGroups( Group * const rootGroup )
{
  GEOSX_MARK_FUNCTION;

  FlowSolverBase::InitializePostInitialConditions_PreSubGroups( rootGroup );

  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);
  MeshLevel * mesh = domain->getMeshBody(0)->getMeshLevel(0);

  ConstitutiveManager * const constitutiveManager = domain->getConstitutiveManager();

  //TODO this is a hack until the sets are fixed to include ghosts!!
  std::map<string, string_array > fieldNames;
  fieldNames["elems"].push_back( viewKeyStruct::pressureString );

  array1d<NeighborCommunicator> & comms =
    domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors );

  CommunicationTools::SynchronizeFields( fieldNames, mesh, comms );

  ResetViews( domain );

  // Moved the following part from ImplicitStepSetup to here since it only needs to be initialized once
  // They will be updated in ApplySystemSolution and ImplicitStepComplete, respectively

  real64 const defaultDensity = constitutiveManager->GetConstitutiveRelation( m_fluidIndex )->
                                getWrapper< array2d<real64> >( SingleFluidBase::viewKeyStruct::densityString )->
                                getDefaultValue();

  applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegionBase * const GEOSX_UNUSED_PARAM( region ),
                                 ElementSubRegionBase * const subRegion )
  {

    subRegion->getWrapper< array1d<real64> >( viewKeyStruct::densityOldString )->
      setDefaultValue( defaultDensity );


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

  ElementRegionManager * const elemManager = mesh->getElemManager();
   elemManager->forElementRegions<FaceElementRegion>( [&] ( FaceElementRegion * const region )
   {
     region->forElementSubRegions<FaceElementSubRegion>( [&]( FaceElementSubRegion * const subRegion )
     {
       subRegion->getWrapper<real64_array>(viewKeyStruct::effectiveApertureString)->
         setApplyDefaultValue(region->getDefaultAperture() );

       subRegion->getWrapper<real64_array>(FaceElementSubRegion::viewKeyStruct::creationMassString)->
         setApplyDefaultValue(defaultDensity * region->getDefaultAperture() );
     });
   });
}

real64 SinglePhaseBase::SolverStep( real64 const& time_n,
                                    real64 const& dt,
                                    const int cycleNumber,
                                    DomainPartition * domain )
{
  GEOSX_MARK_FUNCTION;

  real64 dt_return;

  // setup dof numbers and linear system
  if( !m_coupledWellsFlag )
  {
    SetupSystem( domain, m_dofManager, m_matrix, m_rhs, m_solution );
  }

  ImplicitStepSetup( time_n, dt, domain, m_dofManager, m_matrix, m_rhs, m_solution );

  // currently the only method is implicit time integration
  dt_return = this->NonlinearImplicitStep( time_n, dt, cycleNumber, domain, m_dofManager, m_matrix, m_rhs, m_solution );

  // final step for completion of timestep. typically secondary variable updates and cleanup.
  ImplicitStepComplete( time_n, dt_return, domain );

  return dt_return;
}

void SinglePhaseBase::SetupSystem( DomainPartition * const domain,
                                   DofManager & dofManager,
                                   ParallelMatrix & matrix,
                                   ParallelVector & rhs,
                                   ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;
  ResetViews( domain );


  SolverBase::SetupSystem( domain,
                           dofManager,
                           matrix,
                           rhs,
                           solution );


  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  std::unique_ptr<CRSMatrix<real64,localIndex,localIndex> > &
  derivativeFluxResidual_dAperture = getRefDerivativeFluxResidual_dAperture();
  {

    localIndex numRows = 0;
    localIndex numCols = 0;
    string_array const & flowRegions = getTargetRegions();
    elemManager->forElementSubRegions( flowRegions, [&]( ElementSubRegionBase const * const elementSubRegion )
    {
      numRows += elementSubRegion->size();
      numCols += elementSubRegion->size();
    });

    derivativeFluxResidual_dAperture = std::make_unique<CRSMatrix<real64,localIndex,localIndex>>( numRows, numCols );

    derivativeFluxResidual_dAperture->reserveNonZeros( matrix.localNonzeros() );
    localIndex maxRowSize = -1;
    for( localIndex row=0 ; row<matrix.localRows() ; ++row )
    {
      localIndex const rowSize = matrix.localRowLength( row );
      maxRowSize = maxRowSize > rowSize ? maxRowSize : rowSize;

      derivativeFluxResidual_dAperture->reserveNonZeros( row,
                                                         rowSize );
    }
    for( localIndex row=matrix.localRows() ; row<numRows ; ++row )
    {
      derivativeFluxResidual_dAperture->reserveNonZeros( row,
                                                         maxRowSize );
    }


  }


  string const presDofKey = dofManager.getKey( FlowSolverBase::viewKeyStruct::pressureString );

  NumericalMethodsManager const *
  numericalMethodManager = domain->getParent()->GetGroup<NumericalMethodsManager>( keys::numericalMethodsManager );

  FiniteVolumeManager const *
  fvManager = numericalMethodManager->GetGroup<FiniteVolumeManager>( keys::finiteVolumeManager );

  FluxApproximationBase const * fluxApprox = fvManager->getFluxApproximation( getDiscretization() );


  fluxApprox->forStencils<FaceElementStencil>( [&]( FaceElementStencil const & stencil )
  {
    for( localIndex iconn=0 ; iconn<stencil.size() ; ++iconn)
    {
      localIndex const numFluxElems = stencil.stencilSize(iconn);
//      typename FaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
//      typename FaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
      typename FaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

//      FaceElementSubRegion const * const
//      elementSubRegion = elemManager->GetRegion(seri[iconn][0])->GetSubRegion<FaceElementSubRegion>(sesri[iconn][0]);

      for( localIndex k0=0 ; k0<numFluxElems ; ++k0 )
      {
        for( localIndex k1=0 ; k1<numFluxElems ; ++k1 )
        {
          derivativeFluxResidual_dAperture->insertNonZero( sei[iconn][k0],sei[iconn][k1], 0.0 );
        }
      }
    }
  });
}

void SinglePhaseBase::ImplicitStepSetup( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                         real64 const & GEOSX_UNUSED_PARAM( dt ),
                                         DomainPartition * const domain,
                                         DofManager & GEOSX_UNUSED_PARAM( dofManager ),
                                         ParallelMatrix & GEOSX_UNUSED_PARAM( matrix ),
                                         ParallelVector & GEOSX_UNUSED_PARAM( rhs ),
                                         ParallelVector & GEOSX_UNUSED_PARAM( solution ) )
{
  ResetViews( domain );

  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );

  applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegionBase * const GEOSX_UNUSED_PARAM( region ),
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
                                                                ElementRegionBase *,
                                                                FaceElementSubRegion * subRegion )
  {
    arrayView1d<real64> const & aper0 = subRegion->getReference<array1d<real64>>( viewKeyStruct::aperture0String );
    arrayView1d<real64 const> const & aper = m_effectiveAperture[er][esr];

    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      aper0[ei] = aper[ei];
    } );


    UpdateMobility( subRegion );
  } );

}

void SinglePhaseBase::ImplicitStepComplete( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                            real64 const & GEOSX_UNUSED_PARAM( dt ),
                                            DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );

  applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegionBase * const GEOSX_UNUSED_PARAM( region ),
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


  ElementRegionManager * const elemManager = mesh->getElemManager();
  elemManager->forElementSubRegionsComplete<FaceElementSubRegion>( this->m_targetRegions,
                                                                   [&] ( localIndex er,
                                                                         localIndex esr,
                                                                         ElementRegionBase const * const GEOSX_UNUSED_PARAM( region ),
                                                                         FaceElementSubRegion * const subRegion )
  {

    arrayView1d<integer const> const & elemGhostRank = m_elemGhostRank[er][esr];
    arrayView1d<real64 const> const & densOld = m_densityOld[er][esr];
    arrayView1d<real64 const> const & volume = m_volume[er][esr];
    arrayView1d<real64> const &
    creationMass = subRegion->getReference<real64_array>( FaceElementSubRegion::viewKeyStruct::creationMassString);

    forall_in_range<serialPolicy>( 0, subRegion->size(), GEOSX_LAMBDA ( localIndex ei )
    {
      if (elemGhostRank[ei] < 0)
      {
        if( volume[ei] * densOld[ei] > 1.1 * creationMass[ei] )
        {
          creationMass[ei] *= 0.75;
          if( creationMass[ei]<1.0e-20 )
          {
            creationMass[ei] = 0.0;
          }
        }
      }});
  });
}


void SinglePhaseBase::AssembleSystem( real64 const time_n,
                                      real64 const dt,
                                      DomainPartition * const domain,
                                      DofManager const & dofManager,
                                      ParallelMatrix & matrix,
                                      ParallelVector & rhs )
{
  GEOSX_MARK_FUNCTION;

//  MeshLevel * mesh = domain->getMeshBody(0)->getMeshLevel(0);
//  applyToSubRegions( mesh, [&] ( localIndex , localIndex ,
//                                 ElementRegionBase * const GEOSX_UNUSED_PARAM( region ),
//                                 ElementSubRegionBase * const subRegion )
//  {
//    UpdateState( subRegion );
//  } );



  matrix.zero();
  rhs.zero();

  matrix.open();
  rhs.open();

  if( m_derivativeFluxResidual_dAperture==nullptr )
  {
    m_derivativeFluxResidual_dAperture = std::make_unique<CRSMatrix<real64,localIndex,localIndex>>( matrix.localRows(),
                                                                                                    matrix.localCols() );
  }
  m_derivativeFluxResidual_dAperture->setValues(0.0);

  if (m_poroElasticFlag)
  {
    AssembleAccumulationTerms<true>( domain, &dofManager, &matrix, &rhs );
  }
  else
  {
    AssembleAccumulationTerms<false>( domain, &dofManager, &matrix, &rhs );
  }

  AssembleFluxTerms( time_n, dt, domain, &dofManager, &matrix, &rhs );
 
  
  if (!m_coupledWellsFlag)
  {
    // these functions will be called by the ReservoirSolver
    // when coupled wells are present
    matrix.close();
    rhs.close();
  }
  

  if( getLogLevel() == 2 )
  {
    GEOSX_LOG_RANK_0( "After SinglePhaseFlow::AssembleSystem" );
    GEOSX_LOG_RANK_0("\nJacobian:\n");
    std::cout<< matrix;
    GEOSX_LOG_RANK_0("\nResidual:\n");
    std::cout<< rhs;
  }


  if( getLogLevel() >= 3 )
  {

    integer newtonIter = m_nonlinearSolverParameters.m_numNewtonIterations;

    string filename_mat = "matrix_" + std::to_string( time_n ) + "_" + std::to_string( newtonIter ) + ".mtx";
    matrix.write( filename_mat, true );

    string filename_rhs = "rhs_" + std::to_string( time_n ) + "_" + std::to_string( newtonIter ) + ".mtx";
    rhs.write( filename_rhs, true );

    GEOSX_LOG_RANK_0( "After SinglePhaseBase::AssembleSystem" );
    GEOSX_LOG_RANK_0( "Jacobian: written to " << filename_mat );
    GEOSX_LOG_RANK_0( "Residual: written to " << filename_rhs );
  }

}

template< bool ISPORO >
void SinglePhaseBase::AccumulationLaunch( localIndex const er,
                                          localIndex const esr,
                                          CellElementSubRegion const * const subRegion,
                                          DofManager const * const dofManager,
                                          ParallelMatrix * const matrix,
                                          ParallelVector * const rhs )
{
  string const dofKey = dofManager->getKey( viewKeyStruct::pressureString );
  arrayView1d<globalIndex const> const & dofNumber = subRegion->getReference< array1d<globalIndex> >( dofKey );

  arrayView1d<integer const>     const & elemGhostRank = m_elemGhostRank[er][esr];

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
void SinglePhaseBase::AccumulationLaunch( localIndex const er,
                                          localIndex const esr,
                                          FaceElementSubRegion const * const subRegion,
                                          DofManager const * const dofManager,
                                          ParallelMatrix * const matrix,
                                          ParallelVector * const rhs )
{
  string const dofKey = dofManager->getKey( viewKeyStruct::pressureString );
  arrayView1d<globalIndex const> const & dofNumber = subRegion->getReference< array1d<globalIndex> >( dofKey );

  arrayView1d<integer const>     const & elemGhostRank = m_elemGhostRank[er][esr];

  arrayView1d<real64 const> const & densOld       = m_densityOld[er][esr];
  arrayView1d<real64 const> const & volume        = m_volume[er][esr];
  arrayView1d<real64 const> const & dVol          = m_deltaVolume[er][esr];
  arrayView2d<real64 const> const & dens          = m_density[er][esr][m_fluidIndex];
  arrayView2d<real64 const> const & dDens_dPres   = m_dDens_dPres[er][esr][m_fluidIndex];

#if !defined(ALLOW_CREATION_MASS)
  static_assert(true,"must have ALLOW_CREATION_MASS defined");
#endif
#if ALLOW_CREATION_MASS>0
  arrayView1d<real64 const> const &
  creationMass = subRegion->getReference<real64_array>(FaceElementSubRegion::viewKeyStruct::creationMassString);
#endif
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
#if !defined(ALLOW_CREATION_MASS)
  static_assert(true,"must have ALLOW_CREATION_MASS defined");
#endif
#if ALLOW_CREATION_MASS>0
      if( volume[ei] * densOld[ei] > 1.1 * creationMass[ei] )
      {
        localAccum += creationMass[ei] * 0.25;
      }
#endif
      // add contribution to global residual and jacobian
      matrix->add( elemDOF, elemDOF, localAccumJacobian );
      rhs->add( elemDOF, localAccum );
    }
  } );
}

template< bool ISPORO >
void SinglePhaseBase::AssembleAccumulationTerms( DomainPartition const * const domain,
                                                 DofManager const * const dofManager,
                                                 ParallelMatrix * const matrix,
                                                 ParallelVector * const rhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );


  ElementRegionManager const * const elemManager = mesh->getElemManager();

  elemManager->forElementSubRegionsComplete<CellElementSubRegion,
                                            FaceElementSubRegion>( this->m_targetRegions,
                                                                   [&] ( localIndex er,
                                                                         localIndex esr,
                                                                         ElementRegionBase const * const GEOSX_UNUSED_PARAM( region ),
                                                                         auto const * const subRegion )
  {
    AccumulationLaunch<ISPORO>( er, esr, subRegion, dofManager, matrix, rhs );
  } );
}

void SinglePhaseBase::SolveSystem( DofManager const & dofManager,
                                   ParallelMatrix & matrix,
                                   ParallelVector & rhs,
                                   ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  rhs.scale( -1.0 );
  solution.zero();

  SolverBase::SolveSystem( dofManager, matrix, rhs, solution );

  // Debug for logLevel >= 2
  GEOSX_LOG_LEVEL_RANK_0( 2, "After SinglePhaseBase::SolveSystem" );
  GEOSX_LOG_LEVEL_RANK_0( 2, "\nSolution:\n" << solution );

}

void SinglePhaseBase::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  MeshLevel * mesh = domain->getMeshBody(0)->getMeshLevel(0);

  applyToSubRegions( mesh, [&] ( localIndex er, localIndex esr,
                                 ElementRegionBase * const GEOSX_UNUSED_PARAM( region ),
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

void SinglePhaseBase::ResetViews( DomainPartition * const domain )
{
  FlowSolverBase::ResetViews( domain );

  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = mesh->getElemManager();
  ConstitutiveManager * const constitutiveManager = domain->getConstitutiveManager();

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


} /* namespace geosx */
