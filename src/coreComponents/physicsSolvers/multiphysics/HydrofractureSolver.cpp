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
 * @file HydrofractureSolver.cpp
 */


#include "HydrofractureSolver.hpp"

#include "constitutive/contact/ContactSelector.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/multiphysics/HydrofractureSolverKernels.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/surfaceGeneration/SurfaceGenerator.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

HydrofractureSolver::HydrofractureSolver( const string & name,
                                          Group * const parent )
  : Base( name, parent ),
  m_contactRelationName(),
  m_surfaceGeneratorName(),
  m_surfaceGenerator( nullptr ),
  m_maxNumResolves( 10 )
{
  registerWrapper( viewKeyStruct::surfaceGeneratorNameString(), &m_surfaceGeneratorName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the surface generator to use in the hydrofracture solver" );

  registerWrapper( viewKeyStruct::contactRelationNameString(), &m_contactRelationName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of contact relation to enforce constraints on fracture boundary." );

  registerWrapper( viewKeyStruct::maxNumResolvesString(), &m_maxNumResolves ).
    setApplyDefaultValue( 10 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Value to indicate how many resolves may be executed to perform surface generation after the execution of flow and mechanics solver. " );

  registerWrapper( viewKeyStruct::couplingTypeOptionStringString(), &m_couplingTypeOption ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Coupling method. Valid options:\n* " + EnumStrings< CouplingTypeOption >::concat( "\n* " ) );

  m_numResolves[0] = 0;

  m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::hydrofracture;
  m_linearSolverParameters.get().mgr.separateComponents = false;
  m_linearSolverParameters.get().mgr.displacementFieldName = keys::TotalDisplacement;
  m_linearSolverParameters.get().dofsPerNode = 3;
}


void HydrofractureSolver::registerDataOnMesh( dataRepository::Group & meshBodies )
{
  CoupledSolver::registerDataOnMesh( meshBodies );

  forMeshTargets( meshBodies, [&] ( string const &,
                                    MeshLevel & mesh,
                                    arrayView1d< string const > const & regionNames )
  {

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                              [&]( localIndex const,
                                                                   ElementSubRegionBase & subRegion )
    {
      subRegion.registerWrapper< string >( viewKeyStruct::porousMaterialNamesString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRestartFlags( RestartFlags::NO_WRITE ).
        setSizedFromParent( 0 );
    } );
  } );

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
  meshBodies.forSubGroups< MeshBody >( [&] ( MeshBody & meshBody )
  {
    MeshLevel & meshLevel = *meshBody.getMeshLevel( 0 );

    ElementRegionManager & elemManager = meshLevel.getElemManager();
    elemManager.forElementRegions< SurfaceElementRegion >( [&] ( SurfaceElementRegion * const region )
    {
      region->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion * const subRegion )
      {
        subRegion->registerWrapper< array1d< real64 > >( viewKeyStruct::separationCoeff0String() ).
          setRestartFlags( RestartFlags::NO_WRITE );
        subRegion->registerWrapper< array1d< real64 > >( viewKeyStruct::apertureAtFailureString() ).
          setApplyDefaultValue( -1.0 ).
          setPlotLevel( PlotLevel::LEVEL_0 );

        subRegion->registerWrapper< array1d< real64 > >( FaceElementSubRegion::viewKeyStruct::dSeparationCoeffdAperString() ).
          setRestartFlags( RestartFlags::NO_WRITE );
      } );
    } );
  } );
#endif
}

void HydrofractureSolver::initializePreSubGroups()
{
  CoupledSolver::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel & mesh,
                                                arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager & elementRegionManager = mesh.getElemManager();
    elementRegionManager.forElementSubRegions< ElementSubRegionBase >( regionNames,
                                                                       [&]( localIndex const,
                                                                            ElementSubRegionBase & subRegion )
    {
      string & porousName = subRegion.getReference< string >( viewKeyStruct::porousMaterialNamesString() );
      porousName = getConstitutiveName< CoupledSolidBase >( subRegion );
      GEOSX_ERROR_IF( porousName.empty(), GEOSX_FMT( "Solid model not found on subregion {}", subRegion.getName() ) );
    } );
  } );
}

void HydrofractureSolver::implicitStepSetup( real64 const & time_n,
                                             real64 const & dt,
                                             DomainPartition & domain )
{
  updateDeformationForCoupling( domain );
  CoupledSolver::implicitStepSetup( time_n, dt, domain );

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  mesh.getElemManager().forElementRegions< SurfaceElementRegion >( [&]( SurfaceElementRegion & faceElemRegion )
  {
    faceElemRegion.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
    {
      arrayView1d< real64 > const &
      separationCoeff0 = subRegion.getReference< array1d< real64 > >( viewKeyStruct::separationCoeff0String() );
      arrayView1d< real64 const > const &
      separationCoeff = subRegion.getSeparationCoefficient();
      for( localIndex k=0; k<separationCoeff0.size(); ++k )
      {
        separationCoeff0[k] = separationCoeff[k];
      }
    } );
  } );
#endif

}

void HydrofractureSolver::postProcessInput()
{
  CoupledSolver::postProcessInput();
  m_surfaceGenerator = &this->getParent().getGroup< SurfaceGenerator >( m_surfaceGeneratorName );
}

real64 HydrofractureSolver::solverStep( real64 const & time_n,
                                        real64 const & dt,
                                        int const cycleNumber,
                                        DomainPartition & domain )
{
  real64 dtReturn = dt;

  if( m_couplingTypeOption == CouplingTypeOption::SIM_FixedStress )
  {
    dtReturn = splitOperatorStep( time_n, dt, cycleNumber, domain );
  }
  else if( m_couplingTypeOption == CouplingTypeOption::FIM )
  {
    implicitStepSetup( time_n, dt, domain );

    int const maxIter = m_maxNumResolves + 1;
    m_numResolves[1] = m_numResolves[0];
    int solveIter;
    for( solveIter=0; solveIter<maxIter; ++solveIter )
    {
      int locallyFractured = 0;
      int globallyFractured = 0;

      setupSystem( domain,
                   m_dofManager,
                   m_localMatrix,
                   m_rhs,
                   m_solution );

      // currently the only method is implicit time integration
      dtReturn = nonlinearImplicitStep( time_n, dt, cycleNumber, domain );


      if( m_surfaceGenerator->solverStep( time_n, dt, cycleNumber, domain ) > 0 )
      {
        locallyFractured = 1;
      }
      MpiWrapper::allReduce( &locallyFractured,
                             &globallyFractured,
                             1,
                             MPI_MAX,
                             MPI_COMM_GEOSX );

      if( globallyFractured == 0 )
      {
        break;
      }
      else
      {
        FieldIdentifiers fieldsToBeSync;

        fieldsToBeSync.addElementFields( { extrinsicMeshData::flow::pressure::key(),
                                           extrinsicMeshData::flow::pressure_n::key(),
                                           SurfaceElementSubRegion::viewKeyStruct::elementApertureString() },
                                         { m_surfaceGenerator->getFractureRegionName() } );

        fieldsToBeSync.addFields( FieldLocation::Node,
                                  { keys::IncrementalDisplacement,
                                    keys::TotalDisplacement } );

        CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                             domain.getMeshBody( 0 ).getMeshLevel( 0 ),
                                                             domain.getNeighbors(),
                                                             false );

        this->updateState( domain );

        if( getLogLevel() >= 1 )
        {
          GEOSX_LOG_RANK_0( "++ Fracture propagation. Re-entering Newton Solve." );
        }
      }
    }

    // final step for completion of timestep. typically secondary variable updates and cleanup.
    implicitStepComplete( time_n, dtReturn, domain );
    m_numResolves[1] = solveIter;
  }

  return dtReturn;
}

void HydrofractureSolver::updateDeformationForCoupling( DomainPartition & domain )
{
  MeshLevel & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  ElementRegionManager & elemManager = meshLevel.getElemManager();
  NodeManager const & nodeManager = meshLevel.getNodeManager();
  FaceManager & faceManager = meshLevel.getFaceManager();

  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const u = nodeManager.totalDisplacement();
  arrayView2d< real64 const > const faceNormal = faceManager.faceNormal();
  // arrayView1d<real64 const> const faceArea = faceManager.faceArea();
  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
  {

    ContactBase const & contact = getConstitutiveModel< ContactBase >( subRegion, m_contactRelationName );

    arrayView1d< real64 > const aperture = subRegion.getElementAperture();
    arrayView1d< real64 > const hydraulicAperture = subRegion.getExtrinsicData< extrinsicMeshData::flow::hydraulicAperture >();
    arrayView1d< real64 const > const volume = subRegion.getElementVolume();
    arrayView1d< real64 > const deltaVolume = subRegion.getExtrinsicData< extrinsicMeshData::flow::deltaVolume >();
    arrayView1d< real64 const > const area = subRegion.getElementArea();
    arrayView2d< localIndex const > const elemsToFaces = subRegion.faceList();

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
    arrayView1d< real64 const > const &
    apertureF = subRegion.getReference< array1d< real64 > >( viewKeyStruct::apertureAtFailureString() );

    arrayView1d< real64 > const &
    separationCoeff = subRegion.getSeparationCoefficient();

    arrayView1d< real64 > const &
    dSeparationCoeff_dAper = subRegion.getReference< array1d< real64 > >( FaceElementSubRegion::viewKeyStruct::dSeparationCoeffdAperString() );
    arrayView1d< real64 const > const &
    separationCoeff0 = subRegion.getReference< array1d< real64 > >( viewKeyStruct::separationCoeff0String() );
#endif

    constitutiveUpdatePassThru( contact, [&] ( auto & castedContact )
    {
      using ContactType = TYPEOFREF( castedContact );
      typename ContactType::KernelWrapper contactWrapper = castedContact.createKernelWrapper();

      hydrofractureSolverKernels::DeformationUpdateKernel
        ::launch< parallelDevicePolicy<> >( subRegion.size(),
                                            contactWrapper,
                                            u,
                                            faceNormal,
                                            faceToNodeMap,
                                            elemsToFaces,
                                            area,
                                            volume,
                                            deltaVolume,
                                            aperture,
                                            hydraulicAperture
#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
                                            ,
                                            apertureF,
                                            separationCoeff,
                                            dSeparationCoeff_dAper,
                                            separationCoeff0
#endif
                                            );

    } );

//#if defined(USE_CUDA)
//    deltaVolume.move( LvArray::MemorySpace::cuda );
//    aperture.move( LvArray::MemorySpace::cuda );
//    hydraulicAperture.move( LvArray::MemorySpace::cuda );
//#endif
  } );
}

real64 HydrofractureSolver::splitOperatorStep( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                               real64 const & dt,
                                               integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                                               DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{
  GEOSX_ERROR( "Not implemented" );
  return dt;
}

real64 HydrofractureSolver::explicitStep( real64 const & time_n,
                                          real64 const & dt,
                                          const int cycleNumber,
                                          DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  solidMechanicsSolver()->explicitStep( time_n, dt, cycleNumber, domain );
  flowSolver()->solverStep( time_n, dt, cycleNumber, domain );

  return dt;
}

void HydrofractureSolver::setupCoupling( DomainPartition const & domain,
                                         DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;

  // restrict coupling to fracture regions only (as done originally in setupSystem)
  map< string, array1d< string > > meshTargets;
  forMeshTargets( domain.getMeshBodies(), [&] ( string const & meshBodyName,
                                                MeshLevel const & meshLevel,
                                                arrayView1d< string const > const & regionNames )
  {
    array1d< string > regions;
    ElementRegionManager const & elementRegionManager = meshLevel.getElemManager();
    elementRegionManager.forElementRegions< SurfaceElementRegion >( regionNames,
                                                                    [&]( localIndex const,
                                                                         SurfaceElementRegion const & region )
    {
      regions.emplace_back( region.getName() );
    } );
    meshTargets[meshBodyName] = std::move( regions );
  } );

  dofManager.addCoupling( keys::TotalDisplacement,
                          SinglePhaseBase::viewKeyStruct::elemDofFieldString(),
                          DofManager::Connector::Elem,
                          meshTargets );
}


void HydrofractureSolver::setupSystem( DomainPartition & domain,
                                       DofManager & dofManager,
                                       CRSMatrix< real64, globalIndex > & localMatrix,
                                       ParallelVector & rhs,
                                       ParallelVector & solution,
                                       bool const setSparsity )
{
  GEOSX_MARK_FUNCTION;

  GEOSX_UNUSED_VAR( setSparsity );

  dofManager.setDomain( domain );

  setupDofs( domain, dofManager );
  dofManager.reorderByRank();

  localIndex const numLocalRows = dofManager.numLocalDofs();

  SparsityPattern< globalIndex > patternOriginal;
  dofManager.setSparsityPattern( patternOriginal );

  // Get the original row lengths (diagonal blocks only)
  array1d< localIndex > rowLengths( patternOriginal.numRows() );
  for( localIndex localRow = 0; localRow < patternOriginal.numRows(); ++localRow )
  {
    rowLengths[localRow] = patternOriginal.numNonZeros( localRow );
  }

  // Add the number of nonzeros induced by coupling
  addFluxApertureCouplingNNZ( domain, dofManager, rowLengths.toView() );

  // Create a new pattern with enough capacity for coupled matrix
  SparsityPattern< globalIndex > pattern;
  pattern.resizeFromRowCapacities< parallelHostPolicy >( patternOriginal.numRows(),
                                                         patternOriginal.numColumns(),
                                                         rowLengths.data() );

  // Copy the original nonzeros
  for( localIndex localRow = 0; localRow < patternOriginal.numRows(); ++localRow )
  {
    globalIndex const * cols = patternOriginal.getColumns( localRow ).dataIfContiguous();
    pattern.insertNonZeros( localRow, cols, cols + patternOriginal.numNonZeros( localRow ) );
  }

  // Add the nonzeros from coupling
  addFluxApertureCouplingSparsityPattern( domain, dofManager, pattern.toView() );

  localMatrix.assimilate< parallelDevicePolicy<> >( std::move( pattern ) );
  localMatrix.setName( this->getName() + "/matrix" );

  rhs.setName( this->getName() + "/rhs" );
  rhs.create( numLocalRows, MPI_COMM_GEOSX );

  solution.setName( this->getName() + "/solution" );
  solution.create( numLocalRows, MPI_COMM_GEOSX );

  setUpDflux_dApertureMatrix( domain, dofManager, localMatrix );
}

void HydrofractureSolver::addFluxApertureCouplingNNZ( DomainPartition & domain,
                                                      DofManager & dofManager,
                                                      arrayView1d< localIndex > const & rowLengths ) const
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  ElementRegionManager const & elemManager = mesh.getElemManager();

  string const presDofKey = dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );

  globalIndex const rankOffset = dofManager.rankOffset();

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( flowSolver()->getDiscretizationName() );

  fluxApprox.forStencils< SurfaceElementStencil >( mesh, [&]( SurfaceElementStencil const & stencil )
  {
    for( localIndex iconn=0; iconn<stencil.size(); ++iconn )
    {
      localIndex const numFluxElems = stencil.stencilSize( iconn );
      typename SurfaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
      typename SurfaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
      typename SurfaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

      FaceElementSubRegion const & elementSubRegion =
        elemManager.getRegion( seri[iconn][0] ).getSubRegion< FaceElementSubRegion >( sesri[iconn][0] );

      ArrayOfArraysView< localIndex const > const elemsToNodes = elementSubRegion.nodeList().toViewConst();

      arrayView1d< globalIndex const > const faceElementDofNumber =
        elementSubRegion.getReference< array1d< globalIndex > >( presDofKey );

      for( localIndex k0=0; k0<numFluxElems; ++k0 )
      {
        globalIndex const activeFlowDOF = faceElementDofNumber[sei[iconn][k0]];
        globalIndex const rowNumber = activeFlowDOF - rankOffset;

        if( rowNumber >= 0 && rowNumber < rowLengths.size() )
        {
          for( localIndex k1=0; k1<numFluxElems; ++k1 )
          {
            // The coupling with the nodal displacements of the cell itself has already been added by the dofManager
            // so we only add the coupling with the nodal displacements of the neighbours.
            if( k1 != k0 )
            {
              localIndex const numNodesPerElement = elemsToNodes[sei[iconn][k1]].size();
              rowLengths[rowNumber] += 3*numNodesPerElement;
            }
          }
        }
      }
    }
  } );

}

void HydrofractureSolver::addFluxApertureCouplingSparsityPattern( DomainPartition & domain,
                                                                  DofManager & dofManager,
                                                                  SparsityPatternView< globalIndex > const & pattern ) const
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  NodeManager const & nodeManager = mesh.getNodeManager();
  ElementRegionManager const & elemManager = mesh.getElemManager();

  string const presDofKey = dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );
  string const dispDofKey = dofManager.getKey( keys::TotalDisplacement );

  globalIndex const rankOffset = dofManager.rankOffset();

  arrayView1d< globalIndex const > const & dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( flowSolver()->getDiscretizationName() );

  fluxApprox.forStencils< SurfaceElementStencil >( mesh, [&]( SurfaceElementStencil const & stencil )
  {
    for( localIndex iconn=0; iconn<stencil.size(); ++iconn )
    {
      localIndex const numFluxElems = stencil.stencilSize( iconn );
      typename SurfaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
      typename SurfaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
      typename SurfaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

      FaceElementSubRegion const & elementSubRegion =
        elemManager.getRegion( seri[iconn][0] ).getSubRegion< FaceElementSubRegion >( sesri[iconn][0] );

      ArrayOfArraysView< localIndex const > const elemsToNodes = elementSubRegion.nodeList().toViewConst();

      arrayView1d< globalIndex const > const faceElementDofNumber =
        elementSubRegion.getReference< array1d< globalIndex > >( presDofKey );

      for( localIndex k0=0; k0<numFluxElems; ++k0 )
      {
        globalIndex const activeFlowDOF = faceElementDofNumber[sei[iconn][k0]];

        globalIndex const rowIndex = activeFlowDOF - rankOffset;

        if( rowIndex >= 0 && rowIndex < pattern.numRows() )
        {
          for( localIndex k1=0; k1<numFluxElems; ++k1 )
          {
            // The coupling with the nodal displacements of the cell itself has already been added by the dofManager
            // so we only add the coupling with the nodal displacements of the neighbours.
            if( k1 != k0 )
            {
              localIndex const numNodesPerElement = elemsToNodes[sei[iconn][k1]].size();

              for( localIndex a=0; a<numNodesPerElement; ++a )
              {
                for( int d=0; d<3; ++d )
                {
                  globalIndex const colIndex = dispDofNumber[elemsToNodes[sei[iconn][k1]][a]] + d;
                  pattern.insertNonZero( rowIndex, colIndex );
                }
              }
            }
          }
        }
      }
    }
  } );
}

void HydrofractureSolver::assembleSystem( real64 const time,
                                          real64 const dt,
                                          DomainPartition & domain,
                                          DofManager const & dofManager,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  solidMechanicsSolver()->assembleSystem( time,
                                          dt,
                                          domain,
                                          dofManager,
                                          localMatrix,
                                          localRhs );

  flowSolver()->assembleAccumulationTerms( domain,
                                           dofManager,
                                           localMatrix,
                                           localRhs );
  flowSolver()->assembleHydrofracFluxTerms( time,
                                            dt,
                                            domain,
                                            dofManager,
                                            localMatrix,
                                            localRhs,
                                            getDerivativeFluxResidual_dAperture() );

  assembleForceResidualDerivativeWrtPressure( domain, localMatrix, localRhs );
  assembleFluidMassResidualDerivativeWrtDisplacement( domain, localMatrix );

  this->getRefDerivativeFluxResidual_dAperture()->zero();
}

void
HydrofractureSolver::
  assembleForceResidualDerivativeWrtPressure( DomainPartition & domain,
                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                              arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  FaceManager const & faceManager = mesh.getFaceManager();
  NodeManager & nodeManager = mesh.getNodeManager();
  ElementRegionManager const & elemManager = mesh.getElemManager();

  arrayView2d< real64 const > const & faceNormal = faceManager.faceNormal();
  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager.nodeList().toViewConst();

  arrayView2d< real64 > const &
  fext = nodeManager.getReference< array2d< real64 > >( SolidMechanicsLagrangianFEM::viewKeyStruct::forceExternalString() );
  fext.zero();

  string const presDofKey = m_dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );
  string const dispDofKey = m_dofManager.getKey( keys::TotalDisplacement );

  globalIndex const rankOffset = m_dofManager.rankOffset();
  arrayView1d< globalIndex const > const & dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );

  elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion const & subRegion )
  {

    arrayView1d< globalIndex const > const &
    pressureDofNumber = subRegion.getReference< array1d< globalIndex > >( presDofKey );

    if( subRegion.hasWrapper( extrinsicMeshData::flow::pressure::key() ) )
    {
      arrayView1d< real64 const > const & fluidPressure = subRegion.getReference< array1d< real64 > >( extrinsicMeshData::flow::pressure::key() );
      arrayView1d< real64 const > const & area = subRegion.getElementArea();
      arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList();

      forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const kfe )
      {
        constexpr int kfSign[2] = { -1, 1 };

        real64 Nbar[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( faceNormal[elemsToFaces[kfe][0]] );
        LvArray::tensorOps::subtract< 3 >( Nbar, faceNormal[elemsToFaces[kfe][1]] );
        LvArray::tensorOps::normalize< 3 >( Nbar );

        localIndex const kf0 = elemsToFaces[kfe][0];
        localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kf0 );

        // TODO make if work for any element type.
        globalIndex rowDOF[24];   // Here it assumes 8 nodes?
        real64 nodeRHS[24];   // Here it assumes 8 nodes?
        stackArray2d< real64, 12*12 > dRdP( numNodesPerFace*3, 1 );
        globalIndex colDOF = pressureDofNumber[kfe];

        real64 const Ja = area[kfe] / numNodesPerFace;

        real64 nodalForceMag = fluidPressure[kfe] * Ja;
        real64 nodalForce[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( Nbar );
        LvArray::tensorOps::scale< 3 >( nodalForce, nodalForceMag );

        for( localIndex kf=0; kf<2; ++kf )
        {
          localIndex const faceIndex = elemsToFaces[kfe][kf];


          for( localIndex a=0; a<numNodesPerFace; ++a )
          {

            for( int i=0; i<3; ++i )
            {
              rowDOF[3*a+i] = dispDofNumber[faceToNodeMap( faceIndex, a )] + i;
              nodeRHS[3*a+i] = nodalForce[i] * kfSign[kf];
              fext[faceToNodeMap( faceIndex, a )][i] += nodalForce[i] * kfSign[kf];

              dRdP( 3*a+i, 0 ) = Ja * Nbar[i] * kfSign[kf];
            }
          }

          for( localIndex a=0; a<numNodesPerFace; ++a )
          {
            localIndex const localRow = LvArray::integerConversion< localIndex >( rowDOF[3*a] - rankOffset );
            if( localRow >= 0 && localRow < localMatrix.numRows() )
            {
              for( int i=0; i<3; ++i )
              {
                // TODO: use parallel atomic when loop is parallel
                RAJA::atomicAdd( serialAtomic{}, &localRhs[localRow + i], nodeRHS[3*a+i] );
                localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( localRow + i,
                                                                                  &colDOF,
                                                                                  &dRdP[3*a+i][0],
                                                                                  1 );
              }
            }
          }

        }
      } );
    }
  } );
}

void
HydrofractureSolver::
  assembleFluidMassResidualDerivativeWrtDisplacement( DomainPartition const & domain,
                                                      CRSMatrixView< real64, globalIndex const > const & localMatrix )
{
  GEOSX_MARK_FUNCTION;

  string const presDofKey = m_dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );
  string const dispDofKey = m_dofManager.getKey( keys::TotalDisplacement );

  globalIndex const rankOffset = m_dofManager.rankOffset();

  CRSMatrixView< real64 const, localIndex const > const
  dFluxResidual_dAperture = getDerivativeFluxResidual_dAperture().toViewConst();

  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel const & mesh,
                                                arrayView1d< string const > const & regionNames )
  {
    FaceManager const & faceManager = mesh.getFaceManager();
    NodeManager const & nodeManager = mesh.getNodeManager();
    ElementRegionManager const & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   FaceElementSubRegion const & subRegion )
    {
      ContactBase const & contact = getConstitutiveModel< ContactBase >( subRegion, m_contactRelationName );

      string const & fluidName = subRegion.getReference< string >( FlowSolverBase::viewKeyStruct::fluidNamesString() );
      SingleFluidBase const & fluid = getConstitutiveModel< SingleFluidBase >( subRegion, fluidName );

      arrayView1d< globalIndex const > const presDofNumber = subRegion.getReference< array1d< globalIndex > >( presDofKey );
      arrayView1d< globalIndex const > const dispDofNumber = nodeManager.getReference< array1d< globalIndex > >( dispDofKey );

      arrayView2d< real64 const > const dens = fluid.density();

      arrayView1d< real64 const > const aperture = subRegion.getElementAperture();
      arrayView1d< real64 const > const area = subRegion.getElementArea();

      arrayView2d< localIndex const > const elemsToFaces = subRegion.faceList();
      ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

      arrayView2d< real64 const > const faceNormal = faceManager.faceNormal();

      constitutiveUpdatePassThru( contact, [&] ( auto & castedContact )
      {
        using ContactType = TYPEOFREF( castedContact );
        typename ContactType::KernelWrapper contactWrapper = castedContact.createKernelWrapper();

        hydrofractureSolverKernels::FluidMassResidualDerivativeAssemblyKernel::
          launch< parallelDevicePolicy<> >( subRegion.size(),
                                            rankOffset,
                                            contactWrapper,
                                            elemsToFaces,
                                            faceToNodeMap,
                                            faceNormal,
                                            area,
                                            aperture,
                                            presDofNumber,
                                            dispDofNumber,
                                            dens,
                                            dFluxResidual_dAperture,
                                            localMatrix );

      } );
    } );
  } );
}

void HydrofractureSolver::solveLinearSystem( DofManager const & dofManager,
                                             ParallelMatrix & matrix,
                                             ParallelVector & rhs,
                                             ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  rhs.scale( -1.0 );
  solution.zero();

  SolverBase::solveLinearSystem( dofManager, matrix, rhs, solution );
}

void HydrofractureSolver::updateState( DomainPartition & domain )
{
  updateDeformationForCoupling( domain );
  flowSolver()->updateState( domain );
}

void HydrofractureSolver::setNextDt( real64 const & currentDt,
                                     real64 & nextDt )
{

  if( m_numResolves[0] == 0 && m_numResolves[1] == 0 )
  {
    this->setNextDtBasedOnNewtonIter( currentDt, nextDt );
  }
  else
  {
    SolverBase & surfaceGenerator = this->getParent().getGroup< SolverBase >( "SurfaceGen" );
    nextDt = surfaceGenerator.getTimestepRequest() < 1e99 ? surfaceGenerator.getTimestepRequest() : currentDt;
  }

  GEOSX_LOG_LEVEL_RANK_0( 3, this->getName() << ": nextDt request is "  << nextDt );
}

void HydrofractureSolver::setUpDflux_dApertureMatrix( DomainPartition & domain,
                                                      DofManager const & dofManager,
                                                      CRSMatrix< real64, globalIndex > & localMatrix )
{
  std::unique_ptr< CRSMatrix< real64, localIndex > > &
  derivativeFluxResidual_dAperture = this->getRefDerivativeFluxResidual_dAperture();

  {
    localIndex numRows = 0;
    forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                  MeshLevel & mesh,
                                                  arrayView1d< string const > const & regionNames )
    {
      mesh.getElemManager().forElementSubRegions< FaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                            FaceElementSubRegion const & elementSubRegion )
      {
        numRows += elementSubRegion.size();
      } );
    } );

    derivativeFluxResidual_dAperture = std::make_unique< CRSMatrix< real64, localIndex > >( numRows, numRows );
    derivativeFluxResidual_dAperture->setName( this->getName() + "/derivativeFluxResidual_dAperture" );

    derivativeFluxResidual_dAperture->reserveNonZeros( localMatrix.numNonZeros() );
    localIndex maxRowSize = -1;
    for( localIndex row = 0; row < localMatrix.numRows(); ++row )
    {
      localIndex const rowSize = localMatrix.numNonZeros( row );
      maxRowSize = maxRowSize > rowSize ? maxRowSize : rowSize;
    }
    // TODO This is way too much. The With the full system rowSize is not a good estimate for this.
    for( localIndex row = 0; row < numRows; ++row )
    {
      derivativeFluxResidual_dAperture->reserveNonZeros( row, maxRowSize );
    }
  }

  string const presDofKey = dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( flowSolver()->getDiscretizationName() );
  forMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                MeshLevel const & mesh,
                                                arrayView1d< string const > const & )
  {
    fluxApprox.forStencils< SurfaceElementStencil >( mesh, [&]( SurfaceElementStencil const & stencil )
    {
      for( localIndex iconn = 0; iconn < stencil.size(); ++iconn )
      {
        localIndex const numFluxElems = stencil.stencilSize( iconn );
        typename SurfaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

        for( localIndex k0 = 0; k0 < numFluxElems; ++k0 )
        {
          for( localIndex k1 = 0; k1 < numFluxElems; ++k1 )
          {
            derivativeFluxResidual_dAperture->insertNonZero( sei[iconn][k0], sei[iconn][k1], 0.0 );
          }
        }
      }
    } );
  } );
}

REGISTER_CATALOG_ENTRY( SolverBase, HydrofractureSolver, string const &, Group * const )
} /* namespace geosx */
