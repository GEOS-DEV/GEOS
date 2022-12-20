/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SinglePhaseHybridFVM.cpp
 */

#include "SinglePhaseHybridFVM.hpp"

#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutivePassThru.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "fieldSpecification/AquiferBoundaryCondition.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "finiteVolume/HybridMimeticDiscretization.hpp"
#include "finiteVolume/MimeticInnerProductDispatch.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseFields.hpp"


/**
 * @namespace the geosx namespace that encapsulates the majority of the code
 */
namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace singlePhaseHybridFVMKernels;
using namespace mimeticInnerProduct;

SinglePhaseHybridFVM::SinglePhaseHybridFVM( const string & name,
                                            Group * const parent ):
  SinglePhaseBase( name, parent ),
  m_areaRelTol( 1e-8 )
{

  // one cell-centered dof per cell
  m_numDofPerCell = 1;
  m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::singlePhaseHybridFVM;

}


void SinglePhaseHybridFVM::registerDataOnMesh( Group & meshBodies )
{

  // 1) Register the cell-centered data
  SinglePhaseBase::registerDataOnMesh( meshBodies );

  // 2) Register the face data
  meshBodies.forSubGroups< MeshBody >( [&] ( MeshBody & meshBody )
  {
    MeshLevel & meshLevel = meshBody.getBaseDiscretization();
    FaceManager & faceManager = meshLevel.getFaceManager();

    // primary variables: face pressures at the previous converged time step
    faceManager.registerField< fields::flow::facePressure_n >( getName() );
  } );
}

void SinglePhaseHybridFVM::initializePreSubGroups()
{
  SinglePhaseBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();

  GEOSX_THROW_IF( !fvManager.hasGroup< HybridMimeticDiscretization >( m_discretizationName ),
                  catalogName() << " " << getName() <<
                  ": the HybridMimeticDiscretization must be selected with SinglePhaseHybridFVM",
                  InputError );
}

void SinglePhaseHybridFVM::initializePostInitialConditionsPreSubGroups()
{
  GEOSX_MARK_FUNCTION;

  SinglePhaseBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager const & elemManager = mesh.getElemManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    // in the flux kernel, we need to make sure that we act only on the target regions
    // for that, we need the following region filter
    for( string const & regionName : regionNames )
    {
      m_regionFilter.insert( elemManager.getRegions().getIndex( regionName ) );
    }

    // check that multipliers are stricly larger than 0, which would work with SinglePhaseFVM, but not with SinglePhaseHybridFVM.
    // To deal with a 0 multiplier, we would just have to skip the corresponding face in the FluxKernel
    arrayView1d< real64 const > const transMultiplier = faceManager.getField< fields::flow::transMultiplier >();

    RAJA::ReduceMin< parallelDeviceReduce, real64 > minVal( 1.0 );
    forAll< parallelDevicePolicy<> >( faceManager.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iface )
    {
      minVal.min( transMultiplier[iface] );
    } );

    GEOSX_THROW_IF_LE_MSG( minVal.get(), 0.0,
                           catalogName() << " " << getName() <<
                           "The transmissibility multipliers used in SinglePhaseHybridFVM must strictly larger than 0.0",
                           std::runtime_error );

    FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
    fsManager.forSubGroups< AquiferBoundaryCondition >( [&] ( AquiferBoundaryCondition const & bc )
    {
      GEOSX_LOG_RANK_0( catalogName() << " " << getName() <<
                        "An aquifer boundary condition named " << bc.getName() << " was requested in the XML file. \n"
                                                                                  "This type of boundary condition is not yet supported by SinglePhaseHybridFVM and will be ignored" );
    } );
  } );
}

void SinglePhaseHybridFVM::implicitStepSetup( real64 const & time_n,
                                              real64 const & dt,
                                              DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  // setup the cell-centered fields
  SinglePhaseBase::implicitStepSetup( time_n, dt, domain );

  // setup the face fields
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    FaceManager & faceManager = mesh.getFaceManager();

    // get the face-based pressures
    arrayView1d< real64 const > const & facePres =
      faceManager.getField< fields::flow::facePressure >();
    arrayView1d< real64 > const & facePres_n =
      faceManager.getField< fields::flow::facePressure_n >();
    facePres_n.setValues< parallelDevicePolicy<> >( facePres );
  } );
}

void SinglePhaseHybridFVM::setupDofs( DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
                                      DofManager & dofManager ) const
{

  // setup the connectivity of elem fields
  // we need Connectivity::Face because of the two-point upwinding
  // in AssembleOneSidedMassFluxes
  dofManager.addField( viewKeyStruct::elemDofFieldString(),
                       FieldLocation::Elem,
                       1,
                       getMeshTargets() );

  dofManager.addCoupling( viewKeyStruct::elemDofFieldString(),
                          viewKeyStruct::elemDofFieldString(),
                          DofManager::Connector::Face );

  // setup the connectivity of face fields
  dofManager.addField( fields::flow::facePressure::key(),
                       FieldLocation::Face,
                       1,
                       getMeshTargets() );

  dofManager.addCoupling( fields::flow::facePressure::key(),
                          fields::flow::facePressure::key(),
                          DofManager::Connector::Elem );

  // setup coupling between pressure and face pressure
  dofManager.addCoupling( fields::flow::facePressure::key(),
                          viewKeyStruct::elemDofFieldString(),
                          DofManager::Connector::Elem );
}

void SinglePhaseHybridFVM::assembleFluxTerms( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                              real64 const dt,
                                              DomainPartition const & domain,
                                              DofManager const & dofManager,
                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                              arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  HybridMimeticDiscretization const & hmDiscretization = fvManager.getHybridMimeticDiscretization( m_discretizationName );
  MimeticInnerProductBase const & mimeticInnerProductBase =
    hmDiscretization.getReference< MimeticInnerProductBase >( HybridMimeticDiscretization::viewKeyStruct::innerProductString() );

  string const faceDofKey = dofManager.getKey( fields::flow::facePressure::key() );
  string const elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  // tolerance for transmissibility calculation
  real64 const lengthTolerance = domain.getMeshBody( 0 ).getGlobalLengthScale() * m_areaRelTol;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               arrayView1d< string const > const & regionNames )
  {
    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    mesh.getElemManager().forElementSubRegionsComplete< CellElementSubRegion >( regionNames,
                                                                                [&]( localIndex const,
                                                                                     localIndex const er,
                                                                                     localIndex const esr,
                                                                                     ElementRegionBase const &,
                                                                                     CellElementSubRegion const & subRegion )
    {
      string const & fluidName = subRegion.getReference< string >( viewKeyStruct::fluidNamesString() );
      SingleFluidBase const & fluid = getConstitutiveModel< SingleFluidBase >( subRegion, fluidName );


      string const & permName = subRegion.getReference< string >( viewKeyStruct::permeabilityNamesString() );
      PermeabilityBase const & permeability = getConstitutiveModel< PermeabilityBase >( subRegion, permName );

      singlePhaseHybridFVMKernels::
        ElementBasedAssemblyKernelFactory::
        createAndLaunch< parallelDevicePolicy<> >( dofManager.rankOffset(),
                                                   er,
                                                   esr,
                                                   lengthTolerance,
                                                   elemDofKey,
                                                   faceDofKey,
                                                   getName(),
                                                   nodeManager,
                                                   faceManager,
                                                   mesh.getElemManager(),
                                                   subRegion,
                                                   mimeticInnerProductBase,
                                                   fluid,
                                                   permeability,
                                                   m_regionFilter.toViewConst(),
                                                   dt,
                                                   localMatrix,
                                                   localRhs );
    } );
  } );

}

void SinglePhaseHybridFVM::assemblePoroelasticFluxTerms( real64 const time_n,
                                                         real64 const dt,
                                                         DomainPartition const & domain,
                                                         DofManager const & dofManager,
                                                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                         arrayView1d< real64 > const & localRhs,
                                                         string const & jumpDofKey )
{
  GEOSX_UNUSED_VAR ( jumpDofKey );

  assembleFluxTerms( time_n,
                     dt,
                     domain,
                     dofManager,
                     localMatrix,
                     localRhs );
}

void SinglePhaseHybridFVM::assembleHydrofracFluxTerms( real64 const time_n,
                                                       real64 const dt,
                                                       DomainPartition const & domain,
                                                       DofManager const & dofManager,
                                                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                       arrayView1d< real64 > const & localRhs,
                                                       CRSMatrixView< real64, localIndex const > const & dR_dAper )
{
  GEOSX_UNUSED_VAR ( time_n );
  GEOSX_UNUSED_VAR ( dt );
  GEOSX_UNUSED_VAR ( domain );
  GEOSX_UNUSED_VAR ( dofManager );
  GEOSX_UNUSED_VAR ( localMatrix );
  GEOSX_UNUSED_VAR ( localRhs );
  GEOSX_UNUSED_VAR ( dR_dAper );

  GEOSX_ERROR( "Poroelastic fluxes with conforming fractures not yet implemented." );
}

void SinglePhaseHybridFVM::applyBoundaryConditions( real64 const time_n,
                                                    real64 const dt,
                                                    DomainPartition & domain,
                                                    DofManager const & dofManager,
                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                    arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  SinglePhaseBase::applyBoundaryConditions( time_n, dt, domain, dofManager, localMatrix, localRhs );
  applyFaceDirichletBC( time_n, dt, dofManager, domain, localMatrix, localRhs );
}

namespace
{
char const faceBcLogMessage[] =
  "SinglePhaseHybridFVM {}: at time {}s, "
  "the <{}> boundary condition '{}' is applied to the face set '{}' in '{}'. "
  "\nThe total number of target faces (including ghost faces) is {}. "
  "\nNote that if this number is equal to zero, the boundary condition will not be applied on this face set.";
}

void SinglePhaseHybridFVM::applyFaceDirichletBC( real64 const time_n,
                                                 real64 const dt,
                                                 DofManager const & dofManager,
                                                 DomainPartition & domain,
                                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                 arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  string const faceDofKey = dofManager.getKey( fields::flow::facePressure::key() );

  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                      MeshLevel & mesh,
                                                                      arrayView1d< string const > const & )
  {
    FaceManager & faceManager = mesh.getFaceManager();

    arrayView1d< real64 const > const presFace =
      faceManager.getField< fields::flow::facePressure >();
    arrayView1d< globalIndex const > const faceDofNumber =
      faceManager.getReference< array1d< globalIndex > >( faceDofKey );
    arrayView1d< integer const > const faceGhostRank = faceManager.ghostRank();

    globalIndex const rankOffset = dofManager.rankOffset();

    // take BCs defined for "pressure" field and apply values to "facePressure"
    // this is done this way for consistency with the standard TPFA scheme, which works in the same fashion
    fsManager.apply< FaceManager >( time_n + dt,
                                    mesh,
                                    fields::flow::pressure::key(),
                                    [&] ( FieldSpecificationBase const & fs,
                                          string const & setName,
                                          SortedArrayView< localIndex const > const & targetSet,
                                          FaceManager & targetGroup,
                                          string const & )
    {

      // provide some logging at the first nonlinear iteration
      if( fs.getLogLevel() >= 1 && m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
      {
        globalIndex const numTargetFaces = MpiWrapper::sum< globalIndex >( targetSet.size() );
        GEOSX_LOG_RANK_0( GEOSX_FMT( faceBcLogMessage,
                                     this->getName(), time_n+dt, FieldSpecificationBase::catalogName(),
                                     fs.getName(), setName, targetGroup.getName(), numTargetFaces ) );
      }

      // next, we use the field specification functions to apply the boundary conditions to the system

      // 1. first, populate the face pressure vector at the boundaries of the domain
      fs.applyFieldValue< FieldSpecificationEqual,
                          parallelDevicePolicy<> >( targetSet,
                                                    time_n + dt,
                                                    targetGroup,
                                                    fields::flow::facePressure::key() );

      // 2. second, modify the residual/jacobian matrix as needed to impose the boundary conditions
      forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
      {

        localIndex const kf = targetSet[a];
        if( faceGhostRank[kf] >= 0 )
        {
          return;
        }

        // 2.1 get the dof number of this face
        globalIndex const dofIndex = faceDofNumber[kf];
        localIndex const localRow = dofIndex - rankOffset;
        real64 rhsValue;

        // 2.2 apply field value to the matrix/rhs
        FieldSpecificationEqual::SpecifyFieldValue( dofIndex,
                                                    rankOffset,
                                                    localMatrix,
                                                    rhsValue,
                                                    presFace[kf],
                                                    presFace[kf] );
        localRhs[localRow] = rhsValue;
      } );

    } );

  } );

}

void SinglePhaseHybridFVM::applyAquiferBC( real64 const time,
                                           real64 const dt,
                                           DomainPartition & domain,
                                           DofManager const & dofManager,
                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                           arrayView1d< real64 > const & localRhs ) const
{
  GEOSX_MARK_FUNCTION;

  GEOSX_UNUSED_VAR( time, dt, dofManager, domain, localMatrix, localRhs );
}

void SinglePhaseHybridFVM::saveAquiferConvergedState( real64 const & time,
                                                      real64 const & dt,
                                                      DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  GEOSX_UNUSED_VAR( time, dt, domain );
}


real64 SinglePhaseHybridFVM::calculateResidualNorm( DomainPartition const & domain,
                                                    DofManager const & dofManager,
                                                    arrayView1d< real64 const > const & localRhs )
{
  // here we compute the cell-centered residual norm in the derived class
  // to avoid duplicating a synchronization point

  // get a view into local residual vector

  string const elemDofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );
  string const faceDofKey = dofManager.getKey( fields::flow::facePressure::key() );

  globalIndex const rankOffset = dofManager.rankOffset();

  // local residual
  real64 localResidualNorm[4] = { 0.0, 0.0, 0.0, 0.0 };
  real64 globalResidualNorm[4] = { 0.0, 0.0, 0.0, 0.0 };

  // 1. Compute the residual for the mass conservation equations

  // compute the norm of local residual scaled by cell pore volume

  real64 defaultViscosity = 0; // for the normalization of the face residuals
  localIndex subRegionCounter = 0;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    StencilAccessors< fields::elementVolume > flowAccessors( mesh.getElemManager(), getName() );
    FaceManager const & faceManager = mesh.getFaceManager();

    mesh.getElemManager().forElementSubRegions< ElementSubRegionBase >( regionNames, [&]( localIndex const,
                                                                                          ElementSubRegionBase const & subRegion )
    {

      arrayView1d< globalIndex const > const & elemDofNumber = subRegion.template getReference< array1d< globalIndex > >( elemDofKey );
      arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();
      arrayView1d< real64 const > const & volume = subRegion.getElementVolume();

      SingleFluidBase const & fluidModel =
        getConstitutiveModel< SingleFluidBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::fluidNamesString() ) );
      arrayView2d< real64 const > const & density_n = fluidModel.density_n();

      CoupledSolidBase const & solidModel =
        SolverBase::getConstitutiveModel< CoupledSolidBase >( subRegion, subRegion.template getReference< string >( viewKeyStruct::solidNamesString() ) );
      arrayView2d< real64 const > const & porosity_n = solidModel.getPorosity_n();

      singlePhaseBaseKernels::
        ResidualNormKernel::launch< parallelDevicePolicy<> >( localRhs,
                                                              rankOffset,
                                                              elemDofNumber,
                                                              elemGhostRank,
                                                              volume,
                                                              density_n,
                                                              porosity_n,
                                                              localResidualNorm );

      defaultViscosity += fluidModel.defaultViscosity();
      subRegionCounter++;
    } );

    arrayView1d< integer const > const & faceGhostRank = faceManager.ghostRank();
    arrayView1d< globalIndex const > const & faceDofNumber =
      faceManager.getReference< array1d< globalIndex > >( faceDofKey );

    arrayView2d< localIndex const > const & elemRegionList    = faceManager.elementRegionList();
    arrayView2d< localIndex const > const & elemSubRegionList = faceManager.elementSubRegionList();
    arrayView2d< localIndex const > const & elemList          = faceManager.elementList();

    defaultViscosity /= subRegionCounter;

    // 2. Compute the residual for the face-based constraints
    singlePhaseHybridFVMKernels::
      ResidualNormKernel::launch< parallelDevicePolicy<> >( localRhs,
                                                            rankOffset,
                                                            faceDofNumber.toNestedViewConst(),
                                                            faceGhostRank.toNestedViewConst(),
                                                            elemRegionList.toNestedViewConst(),
                                                            elemSubRegionList.toNestedViewConst(),
                                                            elemList.toNestedViewConst(),
                                                            flowAccessors.get( fields::elementVolume{} ),
                                                            defaultViscosity,
                                                            &localResidualNorm[3] );


  } );
  // 3. Combine the two norms

  // compute global residual norm
  MpiWrapper::allReduce( localResidualNorm,
                         globalResidualNorm,
                         4,
                         MPI_SUM,
                         MPI_COMM_GEOSX );

  real64 const elemResidualNorm = sqrt( globalResidualNorm[0] )
                                  / ( ( globalResidualNorm[1] + m_fluxEstimate ) / (globalResidualNorm[2]+1) );
  real64 const faceResidualNorm = sqrt( globalResidualNorm[3] );

  real64 const residualNorm = ( elemResidualNorm > faceResidualNorm )
                            ? elemResidualNorm
                            : faceResidualNorm;

  if( getLogLevel() >= 1 && logger::internal::rank == 0 )
  {
    GEOSX_FMT( "    ( R{} ) = ( {:4.2e} ) ; ", coupledSolverAttributePrefix(), residualNorm );
  }

  return residualNorm;
}

void SinglePhaseHybridFVM::applySystemSolution( DofManager const & dofManager,
                                                arrayView1d< real64 const > const & localSolution,
                                                real64 const scalingFactor,
                                                DomainPartition & domain )
{
  // here we apply the cell-centered update in the derived class
  // to avoid duplicating a synchronization point

  // 1. apply the cell-centered update

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::elemDofFieldString(),
                               fields::flow::pressure::key(),
                               scalingFactor );

  // 2. apply the face-based update

  dofManager.addVectorToField( localSolution,
                               fields::flow::facePressure::key(),
                               fields::flow::facePressure::key(),
                               scalingFactor );

  // 3. synchronize
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    FieldIdentifiers fieldsToBeSync;

    fieldsToBeSync.addElementFields( { fields::flow::pressure::key() }, regionNames );
    fieldsToBeSync.addFields( FieldLocation::Face, { fields::flow::facePressure::key() } );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync, mesh, domain.getNeighbors(), true );
  } );
}


void SinglePhaseHybridFVM::resetStateToBeginningOfStep( DomainPartition & domain )
{
  // 1. Reset the cell-centered fields
  SinglePhaseBase::resetStateToBeginningOfStep( domain );

  // 2. Reset the face-based fields
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    FaceManager & faceManager = mesh.getFaceManager();

    // get the face pressure update
    arrayView1d< real64 > const & facePres =
      faceManager.getField< fields::flow::facePressure >();
    arrayView1d< real64 const > const & facePres_n =
      faceManager.getField< fields::flow::facePressure_n >();
    facePres.setValues< parallelDevicePolicy<> >( facePres_n );
  } );
}

REGISTER_CATALOG_ENTRY( SolverBase, SinglePhaseHybridFVM, string const &, Group * const )
} /* namespace geosx */
