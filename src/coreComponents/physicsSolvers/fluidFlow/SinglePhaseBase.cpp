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
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseKernels.hpp"

/**
 * @namespace the geosx namespace that encapsulates the majority of the code
 */
namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace SinglePhaseBaseKernels;

SinglePhaseBase::SinglePhaseBase( const std::string & name,
                                  Group * const parent ):
  FlowSolverBase( name, parent )
{
  m_numDofPerCell = 1;
}


void SinglePhaseBase::RegisterDataOnMesh( Group * const MeshBodies )
{
  FlowSolverBase::RegisterDataOnMesh( MeshBodies );

  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    MeshLevel * meshLevel = Group::group_cast< MeshBody * >( mesh.second )->getMeshLevel( 0 );

    ElementRegionManager * const elemManager = meshLevel->getElemManager();

    elemManager->forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
    {
      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::pressureString )->setPlotLevel( PlotLevel::LEVEL_0 );

      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::deltaPressureString )->
        setRestartFlags( RestartFlags::NO_WRITE );

      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::deltaVolumeString )->
        setRestartFlags( RestartFlags::NO_WRITE );

      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::mobilityString );

      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::dMobility_dPressureString )->
        setRestartFlags( RestartFlags::NO_WRITE );

      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::porosityString )->setPlotLevel( PlotLevel::LEVEL_1 );

      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::porosityOldString )->
        setRestartFlags( RestartFlags::NO_WRITE );

      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::densityOldString )->
        setRestartFlags( RestartFlags::NO_WRITE );
    } );

    elemManager->forElementSubRegions< FaceElementSubRegion >( [&] ( FaceElementSubRegion & subRegion )
    {
      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::pressureString )->setPlotLevel( PlotLevel::LEVEL_0 );

      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::deltaPressureString )->
        setRestartFlags( RestartFlags::NO_WRITE );

      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::deltaVolumeString )->
        setRestartFlags( RestartFlags::NO_WRITE );

      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::mobilityString );

      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::dMobility_dPressureString )->
        setRestartFlags( RestartFlags::NO_WRITE );

      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::porosityString )->
        setDefaultValue( 1.0 );

      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::porosityOldString )->
        setDefaultValue( 1.0 )->
        setRestartFlags( RestartFlags::NO_WRITE );

      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::densityOldString )->
        setRestartFlags( RestartFlags::NO_WRITE );

      subRegion.registerWrapper< array1d< R1Tensor > >( viewKeyStruct::transTMultString )->
        setDefaultValue( {1.0, 1.0, 1.0} );
      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::poroMultString )->
        setDefaultValue( 1.0 );
    } );

    // TODO restrict this to boundary sets
    FaceManager * const faceManager = meshLevel->getFaceManager();
    {
      faceManager->registerWrapper< array1d< real64 > >( viewKeyStruct::boundaryFacePressureString );
      faceManager->registerWrapper< array2d< real64 > >( viewKeyStruct::boundaryFaceDensityString )->reference().resizeDimension< 1 >( 1 );
      faceManager->registerWrapper< array2d< real64 > >( viewKeyStruct::boundaryFaceViscosityString )->reference().resizeDimension< 1 >( 1 );
      faceManager->registerWrapper< array1d< real64 > >( viewKeyStruct::boundaryFaceMobilityString );
    }
  }
}

void SinglePhaseBase::ValidateFluidModels( DomainPartition const & domain ) const
{
  for( auto & mesh : domain.getMeshBodies()->GetSubGroups() )
  {
    MeshLevel const & meshLevel = *Group::group_cast< MeshBody const * >( mesh.second )->getMeshLevel( 0 );
    ValidateModelMapping< SingleFluidBase >( *meshLevel.getElemManager(), m_fluidModelNames );
  }
}

void SinglePhaseBase::InitializePreSubGroups( Group * const rootGroup )
{
  FlowSolverBase::InitializePreSubGroups( rootGroup );

  DomainPartition & domain = *rootGroup->GetGroup< DomainPartition >( keys::domain );
  ValidateFluidModels( domain );
}

void SinglePhaseBase::UpdateFluidModel( Group & dataGroup, localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  arrayView1d< real64 const > const & pres = dataGroup.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
  arrayView1d< real64 const > const & dPres = dataGroup.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

  SingleFluidBase & fluid = GetConstitutiveModel< SingleFluidBase >( dataGroup, m_fluidModelNames[targetIndex] );

  forAll< serialPolicy >( dataGroup.size(), [&] ( localIndex const a )
  {
    fluid.PointUpdate( pres[a] + dPres[a], a, 0 );
  } );
}

void SinglePhaseBase::UpdateSolidModel( Group & dataGroup, localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  ConstitutiveBase & solid = GetConstitutiveModel< ConstitutiveBase >( dataGroup, m_solidModelNames[targetIndex] );

  arrayView1d< real64 const > const & pres  = dataGroup.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
  arrayView1d< real64 const > const & dPres = dataGroup.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

  forAll< serialPolicy >( dataGroup.size(), [&] ( localIndex const a )
  {
    solid.StateUpdatePointPressure( pres[a] + dPres[a], a, 0 );
  } );
}


void SinglePhaseBase::UpdateState( Group & dataGroup, localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  UpdateFluidModel( dataGroup, targetIndex );
  UpdateSolidModel( dataGroup, targetIndex );
  UpdateMobility< SingleFluidBase >( dataGroup, targetIndex );
}

void SinglePhaseBase::InitializePostInitialConditions_PreSubGroups( Group * const rootGroup )
{
  GEOSX_MARK_FUNCTION;

  FlowSolverBase::InitializePostInitialConditions_PreSubGroups( rootGroup );

  DomainPartition * domain = rootGroup->GetGroup< DomainPartition >( keys::domain );
  MeshLevel * mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );

  std::map< string, string_array > fieldNames;
  fieldNames["elems"].push_back( viewKeyStruct::pressureString );

  CommunicationTools::SynchronizeFields( fieldNames, mesh, domain->getNeighbors() );

  ResetViews( domain );

  // Moved the following part from ImplicitStepSetup to here since it only needs to be initialized once
  // They will be updated in ApplySystemSolution and ImplicitStepComplete, respectively
  forTargetSubRegionsComplete( *mesh,
                               [&]( localIndex const targetIndex,
                                    localIndex const er,
                                    localIndex const esr,
                                    ElementRegionBase &,
                                    ElementSubRegionBase & subRegion )
  {
    ConstitutiveBase & fluid = GetConstitutiveModel( subRegion, m_fluidModelNames[targetIndex] );
    ConstitutiveBase & solid = GetConstitutiveModel( subRegion, m_solidModelNames[targetIndex] );

    real64 const
    defaultDensity = fluid.getWrapper< array2d< real64 > >( SingleFluidBase::viewKeyStruct::densityString )->getDefaultValue();
    subRegion.getWrapper< array1d< real64 > >( viewKeyStruct::densityOldString )->setDefaultValue( defaultDensity );

    UpdateState( subRegion, targetIndex );

    arrayView1d< real64 const > const & poroRef = m_porosityRef[er][esr];
    arrayView2d< real64 const > const & dens = m_density[er][esr];
    arrayView2d< real64 const > const & pvmult =
      solid.getReference< array2d< real64 > >( ConstitutiveBase::viewKeyStruct::poreVolumeMultiplierString );

    arrayView1d< real64 > const & poro = m_porosity[er][esr];
    arrayView1d< real64 > const & densOld = m_densityOld[er][esr];
    arrayView1d< real64 > const & poroOld = m_porosityOld[er][esr];

    if( pvmult.size() == poro.size() )
    {
      forAll< serialPolicy >( subRegion.size(), [=]( localIndex ei )
      {
        densOld[ei] = dens[ei][0];
        poro[ei] = poroRef[ei] * pvmult[ei][0];
        poroOld[ei] = poro[ei];
      } );
    }
    else
    {
      forAll< serialPolicy >( subRegion.size(), [=]( localIndex ei )
      {
        densOld[ei] = dens[ei][0];
        poro[ei] = poroRef[ei];
        poroOld[ei] = poro[ei];
      } );
    }
  } );

  mesh->getElemManager()->forElementRegions< FaceElementRegion >( targetRegionNames(),
                                                                  [&]( localIndex const targetIndex,
                                                                       FaceElementRegion & region )
  {
    region.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
    {
      ConstitutiveBase & fluid = GetConstitutiveModel( subRegion, m_fluidModelNames[targetIndex] );
      real64 const defaultDensity = fluid.getWrapper< array2d< real64 > >( SingleFluidBase::viewKeyStruct::densityString )->getDefaultValue();

      subRegion.getWrapper< real64_array >( viewKeyStruct::effectiveApertureString )->
        setApplyDefaultValue( region.getDefaultAperture() );

      subRegion.getWrapper< real64_array >( FaceElementSubRegion::viewKeyStruct::creationMassString )->
        setApplyDefaultValue( defaultDensity * region.getDefaultAperture() );
    } );
  } );
}

real64 SinglePhaseBase::SolverStep( real64 const & time_n,
                                    real64 const & dt,
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

  MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  forTargetSubRegionsComplete( mesh,
                               [&]( localIndex const targetIndex,
                                    localIndex const er,
                                    localIndex const esr,
                                    ElementRegionBase &,
                                    ElementSubRegionBase & subRegion )
  {
    arrayView2d< real64 const > const & dens = m_density[er][esr];
    arrayView1d< real64 const > const & poro = m_porosity[er][esr];

    arrayView1d< real64 > const & dPres = m_deltaPressure[er][esr];
    arrayView1d< real64 > const & dVol = m_deltaVolume[er][esr];
    arrayView1d< real64 > const & densOld = m_densityOld[er][esr];
    arrayView1d< real64 > const & poroOld = m_porosityOld[er][esr];

    forAll< serialPolicy >( subRegion.size(), [=]( localIndex ei )
    {
      dPres[ei] = 0.0;
      dVol[ei] = 0.0;
    } );

    // This should fix NaN density in newly created fracture elements
    UpdateState( subRegion, targetIndex );

    forAll< serialPolicy >( subRegion.size(), [=]( localIndex const ei )
    {
      densOld[ei] = dens[ei][0];
      poroOld[ei] = poro[ei];
    } );
  } );

  forTargetSubRegionsComplete< FaceElementSubRegion >( mesh,
                                                       [&]( localIndex const targetIndex,
                                                            localIndex const er,
                                                            localIndex const esr,
                                                            ElementRegionBase &,
                                                            FaceElementSubRegion & subRegion )
  {
    arrayView1d< real64 > const
    & aper0 = subRegion.getReference< array1d< real64 > >( viewKeyStruct::aperture0String );
    arrayView1d< real64 const > const & aper = m_effectiveAperture[er][esr];

    forAll< serialPolicy >( subRegion.size(), [=]( localIndex const ei )
    {
      aper0[ei] = aper[ei];
    } );

    // UpdateMobility( &subRegion );
    UpdateState( subRegion, targetIndex );
  } );

}

void SinglePhaseBase::ImplicitStepComplete( real64 const &  time_n,
                                            real64 const &  dt ,
                                            DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  forTargetSubRegionsComplete( mesh,
                               [&]( localIndex const,
                                    localIndex const er,
                                    localIndex const esr,
                                    ElementRegionBase &,
                                    ElementSubRegionBase & subRegion )
  {
    arrayView1d< real64 > const & pres = m_pressure[er][esr];
    arrayView1d< real64 > const & vol = m_volume[er][esr];

    arrayView1d< real64 const > const & dPres = m_deltaPressure[er][esr];
    arrayView1d< real64 const > const & dVol = m_deltaVolume[er][esr];

    forAll< serialPolicy >( subRegion.size(), [=]( localIndex const ei )
    {
      pres[ei] += dPres[ei];
      vol[ei] += dVol[ei];
    } );

    //TJ: write pressure
    if (subRegion.size() > 0 )
    {
      real64 const viscosityFile = domain->GetGroup("Constitutive")
					 ->GetGroup("water")
					   ->getReference<real64>("defaultViscosity");
      Group * elementSubRegions = domain->GetGroup("MeshBodies")
					->GetGroup<MeshBody>("mesh1")
					->GetGroup<MeshLevel>("Level0")
					->GetGroup<ElementRegionManager>("ElementRegions")
					->GetRegion< FaceElementRegion >( "Fracture" )
					->GetGroup("elementSubRegions");
      FaceElementSubRegion * fracSubRegion = elementSubRegions->GetGroup< FaceElementSubRegion >( "default" );

      std::stringstream stream;
      stream << std::fixed << std::setprecision(1) << time_n+dt;

      string fileName;
      if (viscosityFile < 2.0e-3)
	fileName = "ZeroViscosity_" + stream.str() + "s.hist";
      else
	fileName = "ZeroToughness_" + stream.str() + "s.hist";

      std::ofstream myFile;
      myFile.open(fileName, std::ios::out | std::ios::app);

      for (localIndex ei = 0; ei < subRegion.size(); ei++)
      {
	myFile << subRegion.getElementCenter()[ei][0] << "\t"  //x
	       << subRegion.getElementCenter()[ei][1] << "\t"  //y
	       << subRegion.getElementCenter()[ei][2] << "\t"  //z
	       << fracSubRegion->getElementAperture()[ei]  << "\t"  //aperture
	       << pres[ei] << "\n";                            //pressure
      }

      myFile.close();
    }

  } );


  ElementRegionManager * const elemManager = mesh.getElemManager();
  elemManager->forElementSubRegionsComplete< FaceElementSubRegion >( targetRegionNames(),
                                                                     [&] ( localIndex const,
                                                                           localIndex const er,
                                                                           localIndex const esr,
                                                                           ElementRegionBase &,
                                                                           FaceElementSubRegion & subRegion )
  {
    arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();
    arrayView1d< real64 const > const & densOld = m_densityOld[er][esr];
    arrayView1d< real64 const > const & volume = subRegion.getElementVolume();
    arrayView1d< real64 > const & creationMass =
      subRegion.getReference< real64_array >( FaceElementSubRegion::viewKeyStruct::creationMassString );


    //TJ: what is creationMass?
    forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const ei )
    {
      if( elemGhostRank[ei] < 0 )
      {
        if( volume[ei] * densOld[ei] > 1.1 * creationMass[ei] )
        {
          creationMass[ei] *= 0.75;
          if( creationMass[ei]<1.0e-20 )
          {
            creationMass[ei] = 0.0;
          }
        }
      }
    } );
  } );
}


void SinglePhaseBase::AssembleSystem( real64 const time_n,
                                      real64 const dt,
                                      DomainPartition * const domain,
                                      DofManager const & dofManager,
                                      ParallelMatrix & matrix,
                                      ParallelVector & rhs )
{
  GEOSX_MARK_FUNCTION;

  matrix.open();
  rhs.open();

  if( m_poroElasticFlag )
  {
    AssembleAccumulationTerms< true >( domain, &dofManager, &matrix, &rhs );
  }
  else
  {
    AssembleAccumulationTerms< false >( domain, &dofManager, &matrix, &rhs );
  }

  AssembleFluxTerms( time_n, dt, domain, &dofManager, &matrix, &rhs );

  matrix.close();
  rhs.close();

  if( getLogLevel() == 2 )
  {
    GEOSX_LOG_RANK_0( "After SinglePhaseFlow::AssembleSystem" );
    GEOSX_LOG_RANK_0( "\nJacobian:\n" );
    std::cout << matrix;
    GEOSX_LOG_RANK_0( "\nResidual:\n" );
    std::cout << rhs;
  }


  if( getLogLevel() >= 3 )
  {

    integer newtonIter = m_nonlinearSolverParameters.m_numNewtonIterations;

    string filename_mat = "matrix_" + std::to_string( time_n ) + "_" + std::to_string( newtonIter ) + ".mtx";
    matrix.write( filename_mat, LAIOutputFormat::MATRIX_MARKET );

    string filename_rhs = "rhs_" + std::to_string( time_n ) + "_" + std::to_string( newtonIter ) + ".mtx";
    rhs.write( filename_rhs, LAIOutputFormat::MATRIX_MARKET );

    GEOSX_LOG_RANK_0( "After SinglePhaseBase::AssembleSystem" );
    GEOSX_LOG_RANK_0( "Jacobian: written to " << filename_mat );
    GEOSX_LOG_RANK_0( "Residual: written to " << filename_rhs );
  }

}

template< bool ISPORO >
void SinglePhaseBase::AccumulationLaunch( localIndex const er,
                                          localIndex const esr,
                                          CellElementSubRegion const & subRegion,
                                          DofManager const & dofManager,
                                          ParallelMatrix * const matrix,
                                          ParallelVector * const rhs )
{
  string const dofKey = dofManager.getKey( viewKeyStruct::pressureString );
  arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );

  arrayView1d< integer const > const & elemGhostRank = m_elemGhostRank[er][esr];

  arrayView1d< real64 const > const & dPres         = m_deltaPressure[er][esr];
  arrayView1d< real64 const > const & densOld       = m_densityOld[er][esr];
  arrayView1d< real64 >       const & poro          = m_porosity[er][esr];
  arrayView1d< real64 const > const & poroOld       = m_porosityOld[er][esr];
  arrayView1d< real64 const > const & poroRef       = m_porosityRef[er][esr];
  arrayView1d< real64 const > const & volume        = m_volume[er][esr];
  arrayView1d< real64 const > const & dVol          = m_deltaVolume[er][esr];
  arrayView2d< real64 const > const & dens          = m_density[er][esr];
  arrayView2d< real64 const > const & dDens_dPres   = m_dDens_dPres[er][esr];
  arrayView2d< real64 const > const & pvmult        = m_pvMult[er][esr];
  arrayView2d< real64 const > const & dPVMult_dPres = m_dPvMult_dPres[er][esr];

  arrayView1d< real64 const > const & oldTotalMeanStress = m_poroElasticFlag ? m_totalMeanStressOld[er][esr] : poroOld;
  arrayView1d< real64 const > const & totalMeanStress    = m_poroElasticFlag ? m_totalMeanStress[er][esr]    : poroOld;
  arrayView1d< real64 const > const & bulkModulus        = m_poroElasticFlag ? m_bulkModulus[er][esr]: poroOld;
  real64 const & biotCoefficient = m_poroElasticFlag ? m_biotCoefficient[er][esr] : 0;

  forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const ei )
  {
    if( elemGhostRank[ei] < 0 )
    {
      real64 localAccum, localAccumJacobian;
      globalIndex const elemDOF = dofNumber[ei];

      AccumulationKernel< CellElementSubRegion >::template Compute< ISPORO >( dPres[ei],
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
                                          FaceElementSubRegion const & subRegion,
                                          DofManager const & dofManager,
                                          ParallelMatrix * const matrix,
                                          ParallelVector * const rhs )
{
  string const dofKey = dofManager.getKey( viewKeyStruct::pressureString );
  arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );

  arrayView1d< integer const > const & elemGhostRank = m_elemGhostRank[er][esr];

  arrayView1d< real64 const > const & densOld        = m_densityOld[er][esr];
  arrayView1d< real64 const > const & volume         = m_volume[er][esr];
  arrayView1d< real64 const > const & dVol           = m_deltaVolume[er][esr];
  arrayView2d< real64 const > const & dens           = m_density[er][esr];
  arrayView2d< real64 const > const & dDens_dPres    = m_dDens_dPres[er][esr];
  arrayView1d< real64 const > const & poroMultiplier = m_poroMultiplier[er][esr];

  //TJ:
/*  std::cout << "dens(0) = " << dens[0][0] << std::endl;
  std::cout << "densOld(0) = " << densOld[0] << std::endl;
  std::cout << "dDens_dPres(0) = " << dDens_dPres[0][0] << std::endl;

  if (dens.size()>=2)
  {
    std::cout << "dens(1) = " << dens[1][0] << std::endl;
    std::cout << "densOld(1) = " << densOld[1] << std::endl;
    std::cout << "dDens_dPres(1) = " << dDens_dPres[1][0] << std::endl;
  }
*/

#if !defined(ALLOW_CREATION_MASS)
  static_assert( true, "must have ALLOW_CREATION_MASS defined" );
#endif
#if ALLOW_CREATION_MASS>0
  arrayView1d< real64 const > const &
  creationMass = subRegion.getReference< real64_array >( FaceElementSubRegion::viewKeyStruct::creationMassString );
#endif
  forAll< serialPolicy >( subRegion.size(), [=] ( localIndex ei )
  {
    if( elemGhostRank[ei] < 0 )
    {
      real64 localAccum, localAccumJacobian;
      globalIndex const elemDOF = dofNumber[ei];

      real64 effectiveVolume = volume[ei] * poroMultiplier[ei];

      AccumulationKernel< FaceElementSubRegion >::template Compute< ISPORO >( dens[ei][0],
                                                                              densOld[ei],
                                                                              dDens_dPres[ei][0],
                                                                              effectiveVolume,
                                                                              dVol[ei],
                                                                              localAccum,
                                                                              localAccumJacobian );
#if !defined(ALLOW_CREATION_MASS)
      static_assert( true, "must have ALLOW_CREATION_MASS defined" );
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

  MeshLevel const & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  forTargetSubRegionsComplete< CellElementSubRegion, FaceElementSubRegion >( mesh,
                                                                             [&]( localIndex const,
                                                                                  localIndex const er,
                                                                                  localIndex const esr,
                                                                                  ElementRegionBase const &,
                                                                                  auto const & subRegion )
  {
    AccumulationLaunch< ISPORO >( er, esr, subRegion, *dofManager, matrix, rhs );
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
}

void SinglePhaseBase::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex,
                                  ElementSubRegionBase & subRegion )
  {
    arrayView1d< real64 > const & dPres =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

    forAll< serialPolicy >( subRegion.size(), [=]( localIndex const ei )
    {
      dPres[ei] = 0.0;
    } );

    UpdateState( subRegion, targetIndex );
  } );
}

void SinglePhaseBase::ResetViews( DomainPartition * const domain )
{
  FlowSolverBase::ResetViews( domain );

  MeshLevel * const mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = mesh->getElemManager();

  m_pressure =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::pressureString );
  m_deltaPressure =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::deltaPressureString );
  m_deltaVolume =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::deltaVolumeString );

  m_mobility =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::mobilityString );
  m_dMobility_dPres =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::dMobility_dPressureString );

  m_porosityOld =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::porosityOldString );
  m_densityOld =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::densityOldString );

  m_pvMult =
    elemManager->ConstructMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 > >( ConstitutiveBase::viewKeyStruct::poreVolumeMultiplierString,
                                                                                            targetRegionNames(),
                                                                                            solidModelNames() );
  m_dPvMult_dPres =
    elemManager->ConstructMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 > >( ConstitutiveBase::viewKeyStruct::dPVMult_dPresString,
                                                                                            targetRegionNames(),
                                                                                            solidModelNames() );
  m_porosity =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::porosityString );

  ResetViewsPrivate( elemManager );

  if( m_poroElasticFlag )
  {
    // TODO where are these strings defined?
    m_totalMeanStressOld = elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( "oldTotalMeanStress" );
    m_totalMeanStress    = elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( "totalMeanStress" );

    m_bulkModulus = elemManager->ConstructMaterialViewAccessor< array1d< real64 >, arrayView1d< real64 > >( "BulkModulus",
                                                                                                            targetRegionNames(),
                                                                                                            solidModelNames() );
    m_biotCoefficient = elemManager->ConstructMaterialViewAccessor< real64 >( "BiotCoefficient",
                                                                              targetRegionNames(),
                                                                              solidModelNames() );
  }
}

void SinglePhaseBase::ResetViewsPrivate( ElementRegionManager * const elemManager )
{

  m_density =
    elemManager->ConstructMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 > >( SingleFluidBase::viewKeyStruct::densityString,
                                                                                            targetRegionNames(),
                                                                                            fluidModelNames() );
  m_dDens_dPres =
    elemManager->ConstructMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 > >( SingleFluidBase::viewKeyStruct::dDens_dPresString,
                                                                                            targetRegionNames(),
                                                                                            fluidModelNames() );
  m_viscosity =
    elemManager->ConstructMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 > >( SingleFluidBase::viewKeyStruct::viscosityString,
                                                                                            targetRegionNames(),
                                                                                            fluidModelNames() );
  m_dVisc_dPres =
    elemManager->ConstructMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 > >( SingleFluidBase::viewKeyStruct::dVisc_dPresString,
                                                                                            targetRegionNames(),
                                                                                            fluidModelNames() );

  m_poroMultiplier =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::poroMultString );


  m_transTMultiplier =
    elemManager->ConstructViewAccessor< array1d< R1Tensor >, arrayView1d< R1Tensor > >( viewKeyStruct::transTMultString );

}

} /* namespace geosx */
