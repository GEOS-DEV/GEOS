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
#include "constitutive/fluid/singleFluidSelector.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseKernels.hpp"

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

  bool const success =
  constitutiveUpdatePassThru( fluid, [&]( auto & castedFluid )
  {
    typename TYPEOFREF( castedFluid )::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();
    FluidUpdateKernel::Launch( fluidWrapper, pres, dPres );
  } );
  GEOSX_ERROR_IF( !success, "Kernel not launched due to unknown fluid type" );
}

void SinglePhaseBase::UpdateSolidModel( Group & dataGroup, localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  arrayView1d< real64 const > const & pres  = dataGroup.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
  arrayView1d< real64 const > const & dPres = dataGroup.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

  ConstitutiveBase & solid = GetConstitutiveModel< ConstitutiveBase >( dataGroup, m_solidModelNames[targetIndex] );
  solid.StateUpdateBatchPressure( pres, dPres );
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
  forTargetSubRegions( *mesh, [&]( localIndex const targetIndex,
                                   ElementSubRegionBase & subRegion )
  {
    SingleFluidBase & fluid = GetConstitutiveModel< SingleFluidBase >( subRegion, m_fluidModelNames[targetIndex] );
    ConstitutiveBase & solid = GetConstitutiveModel( subRegion, m_solidModelNames[targetIndex] );

    real64 const
    defaultDensity = fluid.getWrapper< array2d< real64 > >( SingleFluidBase::viewKeyStruct::densityString )->getDefaultValue();
    subRegion.getWrapper< array1d< real64 > >( viewKeyStruct::densityOldString )->setDefaultValue( defaultDensity );

    UpdateState( subRegion, targetIndex );

    arrayView1d< real64 const > const & poroRef = subRegion.getReference< array1d< real64 > >( viewKeyStruct::referencePorosityString );
    arrayView2d< real64 const > const & dens = fluid.density();
    arrayView2d< real64 const > const & pvmult =
      solid.getReference< array2d< real64 > >( ConstitutiveBase::viewKeyStruct::poreVolumeMultiplierString );

    arrayView1d< real64 > const & poro = subRegion.getReference< array1d< real64 > >( viewKeyStruct::porosityString );
    arrayView1d< real64 > const & densOld = subRegion.getReference< array1d< real64 > >( viewKeyStruct::densityOldString );
    arrayView1d< real64 > const & poroOld = subRegion.getReference< array1d< real64 > >( viewKeyStruct::porosityOldString );

    if( pvmult.size() == poro.size() )
    {
      forAll< parallelDevicePolicy< 128 > >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
      {
        densOld[ei] = dens[ei][0];
        poro[ei] = poroRef[ei] * pvmult[ei][0];
        poroOld[ei] = poro[ei];
      } );
    }
    else
    {
      forAll< parallelDevicePolicy< 128 > >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
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

  SolverBase::SetupLocalSystem( domain,
                                dofManager,
                                m_localMatrix,
                                m_localRhs,
                                m_localSolution );

  // TODO: can we make the matrix not lose its callback name after stealing?
  m_localMatrix.setName( this->getName() + "/localMatrix" );
  m_localRhs.setName( this->getName() + "/localRhs" );

  matrix.create( m_localMatrix.toViewConst(), MPI_COMM_GEOSX );
  rhs.create( m_localRhs.toViewConst(), MPI_COMM_GEOSX );
  solution.create( m_localSolution.toViewConst(), MPI_COMM_GEOSX );
}

void SinglePhaseBase::ImplicitStepSetup( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                         real64 const & GEOSX_UNUSED_PARAM( dt ),
                                         DomainPartition * const domain,
                                         DofManager & dofManager,
                                         ParallelMatrix & GEOSX_UNUSED_PARAM( matrix ),
                                         ParallelVector & GEOSX_UNUSED_PARAM( rhs ),
                                         ParallelVector & GEOSX_UNUSED_PARAM( solution ) )
{
  ResetViews( domain );

  MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  string const dofKey = dofManager.getKey( viewKeyStruct::pressureString );
  m_pressureDofIndex.clear();
  m_pressureDofIndex = mesh.getElemManager()->ConstructViewAccessor< array1d< globalIndex >, arrayView1d< globalIndex const > >( dofKey );
  m_pressureDofIndex.setName( getName() + "/accessors/pressureDofIndex" );

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex,
                                  ElementSubRegionBase & subRegion )
  {
    arrayView1d< real64 const > const & poro = subRegion.getReference< array1d< real64 > >( viewKeyStruct::porosityString );

    arrayView1d< real64 > const & dPres = subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );
    arrayView1d< real64 > const & dVol = subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaVolumeString );
    arrayView1d< real64 > const & densOld = subRegion.getReference< array1d< real64 > >( viewKeyStruct::densityOldString );
    arrayView1d< real64 > const & poroOld = subRegion.getReference< array1d< real64 > >( viewKeyStruct::porosityOldString );

    dPres.setValues< parallelDevicePolicy< 128 > >( 0.0 );
    dVol.setValues< parallelDevicePolicy< 128 > >( 0.0 );

    // This should fix NaN density in newly created fracture elements
    UpdateState( subRegion, targetIndex );

    SingleFluidBase const & fluid = GetConstitutiveModel< SingleFluidBase >( subRegion, m_fluidModelNames[targetIndex] );
    arrayView2d< real64 const > const & dens = fluid.density();

    forAll< parallelDevicePolicy< 128 > >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      densOld[ei] = dens[ei][0];
      poroOld[ei] = poro[ei];
    } );
  } );

  forTargetSubRegions< FaceElementSubRegion >( mesh, [&]( localIndex const targetIndex,
                                                          FaceElementSubRegion & subRegion )
  {
    arrayView1d< real64 const > const & aper  = subRegion.getReference< array1d< real64 > >( viewKeyStruct::effectiveApertureString );
    arrayView1d< real64 > const & aper0 = subRegion.getReference< array1d< real64 > >( viewKeyStruct::aperture0String );

    aper0.setValues< parallelDevicePolicy< 128 > >( aper );

    // UpdateMobility( &subRegion );
    UpdateState( subRegion, targetIndex );
  } );

}

void SinglePhaseBase::ImplicitStepComplete( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                            real64 const & GEOSX_UNUSED_PARAM( dt ),
                                            DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  forTargetSubRegions( mesh, [&]( localIndex const,
                                  ElementSubRegionBase & subRegion )
  {
    arrayView1d< real64 const > const & dPres = subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );
    arrayView1d< real64 const > const & dVol = subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaVolumeString );

    arrayView1d< real64 > const & pres = subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
    arrayView1d< real64 > const & vol = subRegion.getElementVolume();

    forAll< parallelDevicePolicy< 128 > >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      pres[ei] += dPres[ei];
      vol[ei] += dVol[ei];
    } );
  } );

  forTargetSubRegions< FaceElementSubRegion >( mesh, [&]( localIndex const,
                                                          FaceElementSubRegion & subRegion )
  {
    arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();
    arrayView1d< real64 const > const & volume = subRegion.getElementVolume();
    arrayView1d< real64 const > const & densOld = subRegion.getReference< array1d< real64 > >( viewKeyStruct::densityOldString );
    arrayView1d< real64 > const & creationMass = subRegion.getReference< real64_array >( FaceElementSubRegion::viewKeyStruct::creationMassString );

    forAll< parallelDevicePolicy< 128 > >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
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

  GEOSX_UNUSED_VAR( matrix )
  GEOSX_UNUSED_VAR( rhs )

  m_localMatrix.setValues< parallelDevicePolicy< 128 > >( 0.0 );
  m_localRhs.setValues< parallelDevicePolicy< 128 > >( 0.0 );

  if( m_poroElasticFlag )
  {
    AssembleAccumulationTerms< true, parallelDevicePolicy< 128 > >( *domain,
                                                                    dofManager,
                                                                    m_localMatrix.toViewConstSizes(),
                                                                    m_localRhs.toView() );
  }
  else
  {
    AssembleAccumulationTerms< false, parallelDevicePolicy< 128 > >( *domain,
                                                                     dofManager,
                                                                     m_localMatrix.toViewConstSizes(),
                                                                     m_localRhs.toView() );
  }

  AssembleFluxTerms( time_n,
                     dt,
                     *domain,
                     dofManager,
                     m_localMatrix.toViewConstSizes(),
                     m_localRhs.toView() );

#if 0 // TODO: debug output should be done in NonlinearImplicitStep()

  matrix.create( m_localMatrix.toViewConst(), MPI_COMM_GEOSX );
  rhs.create( m_localRhs, MPI_COMM_GEOSX );

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
#endif
}

template< bool ISPORO, typename POLICY >
void SinglePhaseBase::AccumulationLaunch( localIndex const targetIndex,
                                          CellElementSubRegion & subRegion,
                                          DofManager const & dofManager,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs )
{
  string const dofKey = dofManager.getKey( viewKeyStruct::pressureString );
  globalIndex const rankOffset = dofManager.rankOffset();
  arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
  arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();

  arrayView1d< real64 const > const & deltaPressure =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );
  arrayView1d< real64 const > const & densityOld =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::densityOldString );
  arrayView1d< real64 > const & porosity =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::porosityString );
  arrayView1d< real64 const > const & porosityOld =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::porosityOldString );
  arrayView1d< real64 const > const & porosityRef =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::referencePorosityString );
  arrayView1d< real64 const > const & volume = subRegion.getElementVolume();
  arrayView1d< real64 const > const & deltaVolume =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaVolumeString );

  arrayView1d< real64 const > const & totalMeanStressOld =
    ISPORO ? subRegion.getReference< array1d< real64 > >( "oldTotalMeanStress" ) : porosityOld;
  arrayView1d< real64 const > const & totalMeanStress =
    ISPORO ? subRegion.getReference< array1d< real64 > >( "totalMeanStress" ) : porosityOld;

  SingleFluidBase const & fluid = GetConstitutiveModel< SingleFluidBase >( subRegion, fluidModelNames()[targetIndex] );
  arrayView2d< real64 const > const & density = fluid.density();
  arrayView2d< real64 const > const & dDens_dPres = fluid.dDensity_dPressure();

  ConstitutiveBase const & solid = GetConstitutiveModel( subRegion, solidModelNames()[targetIndex] );
  arrayView2d< real64 const > const & pvMult =
    solid.getReference< array2d< real64 > >( ConstitutiveBase::viewKeyStruct::poreVolumeMultiplierString );
  arrayView2d< real64 const > const & dPvMult_dPres =
    solid.getReference< array2d< real64 > >( ConstitutiveBase::viewKeyStruct::dPVMult_dPresString );
  arrayView1d< real64 const > const & bulkModulus =
    ISPORO ? solid.getReference< array1d< real64 > >( "BulkModulus" ) : porosityOld;
  real64 const biotCoefficient = ISPORO ? solid.getReference< real64 >( "BiotCoefficient" ) : 0.0;

  using Kernel = AccumulationKernel< CellElementSubRegion >;

  Kernel::template Launch< ISPORO, POLICY >( subRegion.size(),
                                             rankOffset,
                                             dofNumber,
                                             ghostRank,
                                             deltaPressure,
                                             densityOld,
                                             porosity,
                                             porosityOld,
                                             porosityRef,
                                             volume,
                                             deltaVolume,
                                             density,
                                             dDens_dPres,
                                             pvMult,
                                             dPvMult_dPres,
                                             totalMeanStressOld,
                                             totalMeanStress,
                                             bulkModulus,
                                             biotCoefficient,
                                             localMatrix,
                                             localRhs );
}

template< bool ISPORO, typename POLICY >
void SinglePhaseBase::AccumulationLaunch( localIndex const targetIndex,
                                          FaceElementSubRegion const & subRegion,
                                          DofManager const & dofManager,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs )
{
  string const dofKey = dofManager.getKey( viewKeyStruct::pressureString );
  globalIndex const rankOffset = dofManager.rankOffset();
  arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
  arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();

  arrayView1d< real64 const > const & densityOld = subRegion.getReference< array1d< real64 > >( viewKeyStruct::densityOldString );
  arrayView1d< real64 const > const & volume = subRegion.getElementVolume();
  arrayView1d< real64 const > const & deltaVolume = subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaVolumeString );
  arrayView1d< real64 const > const & poroMult = subRegion.getReference< array1d< real64 > >( viewKeyStruct::poroMultString );


  SingleFluidBase const & fluid = GetConstitutiveModel< SingleFluidBase >( subRegion, m_fluidModelNames[targetIndex] );
  arrayView2d< real64 const > const & density = fluid.density();
  arrayView2d< real64 const > const & dDens_dPres = fluid.dDensity_dPressure();

#if !defined(ALLOW_CREATION_MASS)
  static_assert( true, "must have ALLOW_CREATION_MASS defined" );
#endif

#if ALLOW_CREATION_MASS
  arrayView1d< real64 const > const &
  creationMass = subRegion.getReference< real64_array >( FaceElementSubRegion::viewKeyStruct::creationMassString );
#endif

  using Kernel = AccumulationKernel< FaceElementSubRegion >;

  Kernel::template Launch< ISPORO, POLICY >( subRegion.size(),
                                             rankOffset,
                                             dofNumber,
                                             ghostRank,
                                             densityOld,
                                             volume,
                                             deltaVolume,
                                             density,
                                             dDens_dPres,
                                             poroMult,
#if ALLOW_CREATION_MASS
                                             creationMass,
#endif
                                             localMatrix,
                                             localRhs );
}

template< bool ISPORO, typename POLICY >
void SinglePhaseBase::AssembleAccumulationTerms( DomainPartition & domain,
                                                 DofManager const & dofManager,
                                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                 arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  forTargetSubRegions< CellElementSubRegion, FaceElementSubRegion >( mesh,
                                                                     [&]( localIndex const targetIndex,
                                                                          auto & subRegion )
  {
    AccumulationLaunch< ISPORO, POLICY >( targetIndex, subRegion, dofManager, localMatrix, localRhs );
  } );
}

void SinglePhaseBase::ApplyBoundaryConditions( real64 time_n,
                                               real64 dt,
                                               DomainPartition * const domain,
                                               DofManager const & dofManager,
                                               ParallelMatrix & matrix,
                                               ParallelVector & rhs )
{
  GEOSX_MARK_FUNCTION;

  GEOSX_UNUSED_VAR( matrix )
  GEOSX_UNUSED_VAR( rhs )

  ApplySourceFluxBC( time_n, dt, *domain, dofManager, m_localMatrix.toViewConstSizes(), m_localRhs.toView() );
  ApplyDiricletBC( time_n, dt, *domain, dofManager, m_localMatrix.toViewConstSizes(), m_localRhs.toView() );
}

void SinglePhaseBase::ApplyDiricletBC( real64 const time_n,
                                       real64 const dt,
                                       DomainPartition & domain,
                                       DofManager const & dofManager,
                                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                       arrayView1d< real64 > const & localRhs ) const
{
  GEOSX_MARK_FUNCTION;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::get();
  string const dofKey = dofManager.getKey( viewKeyStruct::pressureString );

  fsManager.Apply( time_n + dt,
                   &domain,
                   "ElementRegions",
                   viewKeyStruct::pressureString,
                   [&]( FieldSpecificationBase const * const fs,
                        string const &,
                        SortedArrayView< localIndex const > const & lset,
                        Group * subRegion,
                        string const & )
  {
    arrayView1d< globalIndex const > const & dofNumber =
      subRegion->getReference< array1d< globalIndex > >( dofKey );

    arrayView1d< real64 const > const & pres =
      subRegion->getReference< array1d< real64 > >( viewKeyStruct::pressureString );

    arrayView1d< real64 const > const & dPres =
      subRegion->getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

    // call the application of the boundary condition to alter the matrix and rhs
    fs->ApplyBoundaryConditionToSystem< FieldSpecificationEqual,
                                        parallelDevicePolicy< 128 > >( lset,
                                                                       time_n + dt,
                                                                       subRegion,
                                                                       dofNumber,
                                                                       dofManager.rankOffset(),
                                                                       localMatrix,
                                                                       localRhs,
                                                                       [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      return pres[a] + dPres[a];
    } );
  } );
}

void SinglePhaseBase::ApplySourceFluxBC( real64 const time_n,
                                         real64 const dt,
                                         DomainPartition & domain,
                                         DofManager const & dofManager,
                                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                         arrayView1d< real64 > const & localRhs ) const
{
  GEOSX_MARK_FUNCTION;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::get();
  string const dofKey = dofManager.getKey( viewKeyStruct::pressureString );

  fsManager.Apply( time_n + dt, &domain,
                   "ElementRegions",
                   FieldSpecificationBase::viewKeyStruct::fluxBoundaryConditionString,
                   [&]( FieldSpecificationBase const * const fs,
                        string const &,
                        SortedArrayView< localIndex const > const & lset,
                        Group * subRegion,
                        string const & ) -> void
  {
    arrayView1d< globalIndex const > const &
    dofNumber = subRegion->getReference< array1d< globalIndex > >( dofKey );

    arrayView1d< integer const > const &
    ghostRank = subRegion->getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString );

    SortedArray< localIndex > localSet;
    for( localIndex const a : lset )
    {
      if( ghostRank[a] < 0 )
      {
        localSet.insert( a );
      }
    }

    fs->ApplyBoundaryConditionToSystem< FieldSpecificationAdd,
                                        parallelDevicePolicy< 128 > >( localSet.toViewConst(),
                                                                       time_n + dt,
                                                                       dt,
                                                                       subRegion,
                                                                       dofNumber,
                                                                       dofManager.rankOffset(),
                                                                       localMatrix,
                                                                       localRhs,
                                                                       [] GEOSX_HOST_DEVICE ( localIndex const )
                                                                       {
                                                                         return 0.0;
                                                                       } );

  } );
}

void SinglePhaseBase::SolveSystem( DofManager const & dofManager,
                                   ParallelMatrix & matrix,
                                   ParallelVector & rhs,
                                   ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  matrix.create( m_localMatrix.toViewConst(), MPI_COMM_GEOSX );
  rhs.create( m_localRhs.toViewConst(), MPI_COMM_GEOSX );

  rhs.scale( -1.0 );
  solution.zero();

  SolverBase::SolveSystem( dofManager, matrix, rhs, solution );

  m_localSolution.resize( solution.localSize() );
  solution.extract( m_localSolution );
}

void SinglePhaseBase::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex,
                                  ElementSubRegionBase & subRegion )
  {
    arrayView1d< real64 > const & dPres =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

    dPres.setValues< parallelDevicePolicy< 128 > >( 0.0 );

    UpdateState( subRegion, targetIndex );
  } );
}

void SinglePhaseBase::ResetViews( DomainPartition * const domain )
{
  FlowSolverBase::ResetViews( domain );

  MeshLevel const & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager const & elemManager = *mesh.getElemManager();

  m_pressure.clear();
  m_pressure = elemManager.ConstructArrayViewAccessor< real64, 1 >( viewKeyStruct::pressureString );
  m_pressure.setName( getName() + "/accessors/" + viewKeyStruct::pressureString );

  m_deltaPressure.clear();
  m_deltaPressure = elemManager.ConstructArrayViewAccessor< real64, 1 >( viewKeyStruct::deltaPressureString );
  m_deltaPressure.setName( getName() + "/accessors/" + viewKeyStruct::deltaPressureString );

  m_volume.clear();
  m_volume = elemManager.ConstructArrayViewAccessor< real64, 1 >( ElementSubRegionBase::viewKeyStruct::elementVolumeString );
  m_volume.setName( string( "accessors/" ) + ElementSubRegionBase::viewKeyStruct::elementVolumeString );

  m_deltaVolume.clear();
  m_deltaVolume = elemManager.ConstructArrayViewAccessor< real64, 1 >( viewKeyStruct::deltaVolumeString );
  m_deltaVolume.setName( getName() + "/accessors/" + viewKeyStruct::deltaVolumeString );

  m_mobility.clear();
  m_mobility = elemManager.ConstructArrayViewAccessor< real64, 1 >( viewKeyStruct::mobilityString );
  m_mobility.setName( getName() + "/accessors/" + viewKeyStruct::mobilityString );

  m_dMobility_dPres.clear();
  m_dMobility_dPres = elemManager.ConstructArrayViewAccessor< real64, 1 >( viewKeyStruct::dMobility_dPressureString );
  m_dMobility_dPres.setName( getName() + "/accessors/" + viewKeyStruct::dMobility_dPressureString );

  ResetViewsPrivate( elemManager );
}

void SinglePhaseBase::ResetViewsPrivate( ElementRegionManager const & elemManager )
{
  m_density.clear();
  m_density = elemManager.ConstructMaterialArrayViewAccessor< real64, 2 >( SingleFluidBase::viewKeyStruct::densityString,
                                                                           targetRegionNames(),
                                                                           fluidModelNames() );
  m_density.setName( getName() + "/accessors/" + SingleFluidBase::viewKeyStruct::densityString );

  m_dDens_dPres.clear();
  m_dDens_dPres = elemManager.ConstructMaterialArrayViewAccessor< real64, 2 >( SingleFluidBase::viewKeyStruct::dDens_dPresString,
                                                                               targetRegionNames(),
                                                                               fluidModelNames() );
  m_dDens_dPres.setName( getName() + "/accessors/" + SingleFluidBase::viewKeyStruct::dDens_dPresString );

  m_viscosity.clear();
  m_viscosity = elemManager.ConstructMaterialArrayViewAccessor< real64, 2 >( SingleFluidBase::viewKeyStruct::viscosityString,
                                                                             targetRegionNames(),
                                                                             fluidModelNames() );
  m_viscosity.setName( getName() + "/accessors/" + SingleFluidBase::viewKeyStruct::viscosityString );

  m_dVisc_dPres.clear();
  m_dVisc_dPres = elemManager.ConstructMaterialArrayViewAccessor< real64, 2 >( SingleFluidBase::viewKeyStruct::dVisc_dPresString,
                                                                               targetRegionNames(),
                                                                               fluidModelNames() );
  m_dVisc_dPres.setName( getName() + "/accessors/" + SingleFluidBase::viewKeyStruct::dVisc_dPresString );

  m_poroMultiplier.clear();
  m_poroMultiplier = elemManager.ConstructArrayViewAccessor< real64, 1 >( viewKeyStruct::poroMultString );
  m_poroMultiplier.setName( getName() + "/accessors/" + viewKeyStruct::poroMultString );

  m_transTMultiplier.clear();
  m_transTMultiplier = elemManager.ConstructArrayViewAccessor< R1Tensor, 1 >( viewKeyStruct::transTMultString );
  m_transTMultiplier.setName( getName() + "/accessors/" + viewKeyStruct::transTMultString );
}

} /* namespace geosx */
