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
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "physicsSolvers/multiphysics/HydrofractureSolverKernels.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"
#include "physicsSolvers/multiphysics/SinglePhasePoromechanics.hpp"
#include "physicsSolvers/multiphysics/MultiphasePoromechanics.hpp"


namespace geos
{

using namespace constitutive;
using namespace dataRepository;
using namespace fields;


namespace
{

// This is meant to be specialized to work, see below
template< typename POROMECHANICS_SOLVER > class
  HydrofractureSolverCatalogNames {};

// Class specialization for a POROMECHANICS_SOLVER set to SinglePhasePoromechanics
template<> class HydrofractureSolverCatalogNames< SinglePhasePoromechanics >
{
public:
  static string name() { return "Hydrofracture"; }
};

// Class specialization for a POROMECHANICS_SOLVER set to MultiphasePoromechanics
template<> class HydrofractureSolverCatalogNames< MultiphasePoromechanics >
{
public:
  static string name() { return "MultiphaseHydrofracture"; }
};
}

// provide a definition for catalogName()
template< typename POROMECHANICS_SOLVER >
string
HydrofractureSolver< POROMECHANICS_SOLVER >::
catalogName()
{
  return HydrofractureSolverCatalogNames< POROMECHANICS_SOLVER >::name();
}


template< typename POROMECHANICS_SOLVER >
HydrofractureSolver< POROMECHANICS_SOLVER >::HydrofractureSolver( const string & name,
                                                                  Group * const parent )
  : Base( name, parent ),
  m_contactRelationName(),
  m_surfaceGeneratorName(),
  m_surfaceGenerator( nullptr ),
  m_maxNumResolves( 10 ),
  m_isMatrixPoroelastic()
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

  registerWrapper( viewKeyStruct::isMatrixPoroelasticString(), &m_isMatrixPoroelastic ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL );

  registerWrapper( viewKeyStruct::useQuasiNewtonString(), &m_useQuasiNewton ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL );

  m_numResolves[0] = 0;

  // This may need to be different depending on whether poroelasticity is on or not.
  m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::hydrofracture;
  m_linearSolverParameters.get().mgr.separateComponents = false;
  m_linearSolverParameters.get().mgr.displacementFieldName = solidMechanics::totalDisplacement::key();
  m_linearSolverParameters.get().dofsPerNode = 3;

}

template< typename POROMECHANICS_SOLVER >
void HydrofractureSolver< POROMECHANICS_SOLVER >::registerDataOnMesh( dataRepository::Group & meshBodies )
{
  Base::registerDataOnMesh( meshBodies );

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
  meshBodies.forSubGroups< MeshBody >( [&] ( MeshBody & meshBody )
  {
    MeshLevel & meshLevel = *meshBody.getBaseDiscretization();

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

template< typename POROMECHANICS_SOLVER >
void HydrofractureSolver< POROMECHANICS_SOLVER >::implicitStepSetup( real64 const & time_n,
                                                                     real64 const & dt,
                                                                     DomainPartition & domain )
{
  updateDeformationForCoupling( domain );
  Base::implicitStepSetup( time_n, dt, domain );

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
  MeshLevel & mesh = domain.getMeshBody( 0 ).getBaseDiscretization();

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

template< typename POROMECHANICS_SOLVER >
void HydrofractureSolver< POROMECHANICS_SOLVER >::postProcessInput()
{
  Base::postProcessInput();

  static const std::set< integer > binaryOptions = { 0, 1 };
  GEOS_ERROR_IF( binaryOptions.count( m_isMatrixPoroelastic ) == 0, viewKeyStruct::isMatrixPoroelasticString() << " option can be either 0 (false) or 1 (true)" );

  m_surfaceGenerator = &this->getParent().template getGroup< SurfaceGenerator >( m_surfaceGeneratorName );

  flowSolver()->allowNegativePressure();
}

template< typename POROMECHANICS_SOLVER >
real64 HydrofractureSolver< POROMECHANICS_SOLVER >::fullyCoupledSolverStep( real64 const & time_n,
                                                                            real64 const & dt,
                                                                            int const cycleNumber,
                                                                            DomainPartition & domain )
{
  real64 dtReturn = dt;

  implicitStepSetup( time_n, dt, domain );

  int const maxIter = m_maxNumResolves + 1;
  m_numResolves[1] = m_numResolves[0];
  int solveIter;
  for( solveIter=0; solveIter<maxIter; ++solveIter )
  {
    GEOS_LOG_RANK_0( GEOS_FMT( "  Fracture propagation iteration {}", solveIter ) );

    int locallyFractured = 0;
    int globallyFractured = 0;

    Timestamp const meshModificationTimestamp = getMeshModificationTimestamp( domain );

    // Only build the sparsity pattern if the mesh has changed
    if( meshModificationTimestamp > getSystemSetupTimestamp() )
    {
      setupSystem( domain,
                   m_dofManager,
                   m_localMatrix,
                   m_rhs,
                   m_solution );
      setSystemSetupTimestamp( meshModificationTimestamp );
    }

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

      fieldsToBeSync.addElementFields( { flow::pressure::key(),
                                         flow::pressure_n::key(),
                                         SurfaceElementSubRegion::viewKeyStruct::elementApertureString() },
                                       { m_surfaceGenerator->getFractureRegionName() } );

      fieldsToBeSync.addFields( FieldLocation::Node,
                                { solidMechanics::incrementalDisplacement::key(),
                                  solidMechanics::totalDisplacement::key() } );

      CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                           domain.getMeshBody( 0 ).getBaseDiscretization(),
                                                           domain.getNeighbors(),
                                                           false );

      this->updateState( domain );

      if( getLogLevel() >= 1 )
      {
        GEOS_LOG_RANK_0( "++ Fracture propagation. Re-entering Newton Solve." );
      }
    }
  }

  // final step for completion of timestep. typically secondary variable updates and cleanup.
  implicitStepComplete( time_n, dtReturn, domain );
  m_numResolves[1] = solveIter;

  return dtReturn;
}
template< typename POROMECHANICS_SOLVER >
void HydrofractureSolver< POROMECHANICS_SOLVER >::updateDeformationForCoupling( DomainPartition & domain )
{
  MeshLevel & meshLevel = domain.getMeshBody( 0 ).getBaseDiscretization();
  ElementRegionManager & elemManager = meshLevel.getElemManager();
  NodeManager const & nodeManager = meshLevel.getNodeManager();
  FaceManager & faceManager = meshLevel.getFaceManager();

  solidMechanics::arrayViewConst2dLayoutTotalDisplacement const u =
    nodeManager.getField< solidMechanics::totalDisplacement >();
  arrayView2d< real64 const > const faceNormal = faceManager.faceNormal();
  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  real64 maxApertureChange( 0.0 );
  real64 maxHydraulicApertureChange( 0.0 );
  real64 minAperture( 1e10 );
  real64 maxAperture( -1e10 );
  real64 minHydraulicAperture( 1e10 );
  real64 maxHydraulicAperture( -1e10 );

  elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
  {

    ContactBase const & contact = this->template getConstitutiveModel< ContactBase >( subRegion, m_contactRelationName );

    arrayView1d< real64 > const aperture = subRegion.getElementAperture();
    arrayView1d< real64 > const hydraulicAperture = subRegion.getField< flow::hydraulicAperture >();
    arrayView1d< real64 const > const volume = subRegion.getElementVolume();
    arrayView1d< real64 > const deltaVolume = subRegion.getField< flow::deltaVolume >();
    arrayView1d< real64 const > const area = subRegion.getElementArea();
    ArrayOfArraysView< localIndex const > const elemsToFaces = subRegion.faceList().toViewConst();

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

      auto const statistics = hydrofractureSolverKernels::DeformationUpdateKernel
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

      maxApertureChange = std::max( maxApertureChange, std::get< 0 >( statistics ));
      maxHydraulicApertureChange = std::max( maxHydraulicApertureChange, std::get< 1 >( statistics ));
      minAperture = std::min( minAperture, std::get< 2 >( statistics ));
      maxAperture = std::max( maxAperture, std::get< 3 >( statistics ));
      minHydraulicAperture = std::min( minHydraulicAperture, std::get< 4 >( statistics ));
      maxHydraulicAperture = std::max( maxHydraulicAperture, std::get< 5 >( statistics ));

    } );

//#if defined(USE_CUDA)
//    deltaVolume.move( parallelDeviceMemorySpace );
//    aperture.move( parallelDeviceMemorySpace );
//    hydraulicAperture.move( parallelDeviceMemorySpace );
//#endif
  } );

  maxApertureChange = MpiWrapper::max( maxApertureChange );
  maxHydraulicApertureChange = MpiWrapper::max( maxHydraulicApertureChange );
  minAperture  = MpiWrapper::min( minAperture );
  maxAperture  = MpiWrapper::max( maxAperture );
  minHydraulicAperture  = MpiWrapper::min( minHydraulicAperture );
  maxHydraulicAperture  = MpiWrapper::max( maxHydraulicAperture );

  GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "        {}: Max aperture change: {} m, max hydraulic aperture change: {} m",
                                      this->getName(), fmt::format( "{:.{}f}", maxApertureChange, 6 ), fmt::format( "{:.{}f}", maxHydraulicApertureChange, 6 ) ) );
  GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "        {}: Min aperture: {} m, max aperture: {} m",
                                      this->getName(), fmt::format( "{:.{}f}", minAperture, 6 ), fmt::format( "{:.{}f}", maxAperture, 6 ) ) );
  GEOS_LOG_LEVEL_RANK_0( 1, GEOS_FMT( "        {}: Min hydraulic aperture: {} m, max hydraulic aperture: {} m",
                                      this->getName(), fmt::format( "{:.{}f}", minHydraulicAperture, 6 ), fmt::format( "{:.{}f}", maxHydraulicAperture, 6 ) ) );

}
template< typename POROMECHANICS_SOLVER >
void HydrofractureSolver< POROMECHANICS_SOLVER >::setupCoupling( DomainPartition const & domain,
                                                                 DofManager & dofManager ) const
{
  GEOS_MARK_FUNCTION;

  if( m_isMatrixPoroelastic )
  {
    Base::setupCoupling( domain, dofManager );
  }
  else
  {
    string const solidDiscretizationName = solidMechanicsSolver()->getDiscretizationName();
    string const flowDiscretizationName = flowSolver()->getDiscretizationName();

    // restrict coupling to fracture regions only (as done originally in setupSystem)
    map< std::pair< string, string >, array1d< string > > dispMeshTargets;
    map< std::pair< string, string >, array1d< string > > presMeshTargets;

    forDiscretizationOnMeshTargets( domain.getMeshBodies(),
                                    [&] ( string const & meshBodyName,
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

      dispMeshTargets[std::make_pair( meshBodyName, solidDiscretizationName )] = std::move( regions );
      presMeshTargets[std::make_pair( meshBodyName, flowDiscretizationName )] = std::move( regions );
    } );

    dofManager.addCoupling( solidMechanics::totalDisplacement::key(),
                            SinglePhaseBase::viewKeyStruct::elemDofFieldString(),
                            DofManager::Connector::Elem,
                            dispMeshTargets );
  }
}

template< typename POROMECHANICS_SOLVER >
void HydrofractureSolver< POROMECHANICS_SOLVER >::setupSystem( DomainPartition & domain,
                                                               DofManager & dofManager,
                                                               CRSMatrix< real64, globalIndex > & localMatrix,
                                                               ParallelVector & rhs,
                                                               ParallelVector & solution,
                                                               bool const setSparsity )
{
  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( setSparsity );

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

template< typename POROMECHANICS_SOLVER >
void HydrofractureSolver< POROMECHANICS_SOLVER >::addFluxApertureCouplingNNZ( DomainPartition & domain,
                                                                              DofManager & dofManager,
                                                                              arrayView1d< localIndex > const & rowLengths ) const
{
  GEOS_MARK_FUNCTION;

  MeshLevel & mesh = domain.getMeshBody( 0 ).getBaseDiscretization();

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
template< typename POROMECHANICS_SOLVER >
void HydrofractureSolver< POROMECHANICS_SOLVER >::addFluxApertureCouplingSparsityPattern( DomainPartition & domain,
                                                                                          DofManager & dofManager,
                                                                                          SparsityPatternView< globalIndex > const & pattern ) const
{
  GEOS_MARK_FUNCTION;

  MeshLevel & mesh = domain.getMeshBody( 0 ).getBaseDiscretization();

  NodeManager const & nodeManager = mesh.getNodeManager();
  ElementRegionManager const & elemManager = mesh.getElemManager();

  string const presDofKey = dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );
  string const dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );

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
template< typename POROMECHANICS_SOLVER >
void HydrofractureSolver< POROMECHANICS_SOLVER >::assembleSystem( real64 const time,
                                                                  real64 const dt,
                                                                  DomainPartition & domain,
                                                                  DofManager const & dofManager,
                                                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                  arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  if( m_isMatrixPoroelastic )
  {
    assembleElementBasedTerms( time,
                               dt,
                               domain,
                               dofManager,
                               localMatrix,
                               localRhs );


    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                 MeshLevel & mesh,
                                                                 arrayView1d< string const > const & regionNames )
    {
      mesh.getElemManager().forElementSubRegions< SurfaceElementSubRegion >( regionNames,
                                                                             [&]( localIndex const,
                                                                                  SurfaceElementSubRegion & subRegion )
      {
        flowSolver()->accumulationAssemblyLaunch( dofManager,
                                                  subRegion,
                                                  localMatrix,
                                                  localRhs );
      } );
    } );
  }
  else
  {

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
  }

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

template< typename POROMECHANICS_SOLVER >
void HydrofractureSolver< POROMECHANICS_SOLVER >::
assembleForceResidualDerivativeWrtPressure( DomainPartition & domain,
                                            CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                            arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;
  MeshLevel & mesh = domain.getMeshBody( 0 ).getBaseDiscretization();

  FaceManager const & faceManager = mesh.getFaceManager();
  NodeManager & nodeManager = mesh.getNodeManager();
  ElementRegionManager const & elemManager = mesh.getElemManager();

  arrayView2d< real64 const > const & faceNormal = faceManager.faceNormal();
  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager.nodeList().toViewConst();

  arrayView2d< real64 > const &
  fext = nodeManager.getField< solidMechanics::externalForce >();
  fext.zero();

  string const presDofKey = m_dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );
  string const dispDofKey = m_dofManager.getKey( solidMechanics::totalDisplacement::key() );

  globalIndex const rankOffset = m_dofManager.rankOffset();
  arrayView1d< globalIndex const > const & dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );

  elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion const & subRegion )
  {

    arrayView1d< globalIndex const > const &
    pressureDofNumber = subRegion.getReference< array1d< globalIndex > >( presDofKey );

    if( subRegion.hasField< flow::pressure >() )
    {
      arrayView1d< real64 const > const & fluidPressure = subRegion.getField< flow::pressure >();
      arrayView1d< real64 const > const & area = subRegion.getElementArea();
      ArrayOfArraysView< localIndex const > const & elemsToFaces = subRegion.faceList().toViewConst();

      // if matching on lassen/crusher, move to device policy
      using execPolicy = serialPolicy;
      forAll< execPolicy >( subRegion.size(), [=] ( localIndex const kfe )
      {
        if( elemsToFaces.sizeOfArray( kfe ) != 2 )
        {
          return;
        }

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
              RAJA::atomicAdd( AtomicPolicy< execPolicy >{}, &(fext[faceToNodeMap( faceIndex, a )][i]), nodalForce[i] * kfSign[kf] );

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
                RAJA::atomicAdd( AtomicPolicy< execPolicy >{}, &localRhs[localRow + i], nodeRHS[3*a+i] );
                localMatrix.addToRowBinarySearchUnsorted< AtomicPolicy< execPolicy > >( localRow + i,
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

template< typename POROMECHANICS_SOLVER >
void HydrofractureSolver< POROMECHANICS_SOLVER >::
assembleFluidMassResidualDerivativeWrtDisplacement( DomainPartition const & domain,
                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix )
{
  GEOS_MARK_FUNCTION;

  string const presDofKey = m_dofManager.getKey( SinglePhaseBase::viewKeyStruct::elemDofFieldString() );
  string const dispDofKey = m_dofManager.getKey( solidMechanics::totalDisplacement::key() );

  globalIndex const rankOffset = m_dofManager.rankOffset();

  CRSMatrixView< real64 const, localIndex const > const
  dFluxResidual_dAperture = getDerivativeFluxResidual_dAperture().toViewConst();

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
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
      ContactBase const & contact = this->template getConstitutiveModel< ContactBase >( subRegion, m_contactRelationName );

      string const & fluidName = subRegion.getReference< string >( FlowSolverBase::viewKeyStruct::fluidNamesString() );
      SingleFluidBase const & fluid = this->template getConstitutiveModel< SingleFluidBase >( subRegion, fluidName );

      arrayView1d< globalIndex const > const presDofNumber = subRegion.getReference< array1d< globalIndex > >( presDofKey );
      arrayView1d< globalIndex const > const dispDofNumber = nodeManager.getReference< array1d< globalIndex > >( dispDofKey );

      arrayView2d< real64 const > const dens = fluid.density();

      arrayView1d< real64 const > const aperture = subRegion.getElementAperture();
      arrayView1d< real64 const > const area = subRegion.getElementArea();

      ArrayOfArraysView< localIndex const > const elemsToFaces = subRegion.faceList().toViewConst();
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
                                            m_useQuasiNewton,
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
template< typename POROMECHANICS_SOLVER >
void HydrofractureSolver< POROMECHANICS_SOLVER >::updateState( DomainPartition & domain )
{
  // remove the contribution of the hydraulic aperture from the stencil weights
  flowSolver()->prepareStencilWeights( domain );

  updateDeformationForCoupling( domain );

  // update the stencil weights using the updated hydraulic aperture
  flowSolver()->updateStencilWeights( domain );

  flowSolver()->updateState( domain );
}
template< typename POROMECHANICS_SOLVER >
real64 HydrofractureSolver< POROMECHANICS_SOLVER >::setNextDt( real64 const & currentDt,
                                                               DomainPartition & domain )
{
  GEOS_UNUSED_VAR( domain );
  real64 nextDt = 0.0;

  if( m_numResolves[0] == 0 && m_numResolves[1] == 0 )
  {
    nextDt = this->setNextDtBasedOnNewtonIter( currentDt );
  }
  else
  {
    nextDt = m_surfaceGenerator->getTimestepRequest() < 1e99 ? m_surfaceGenerator->getTimestepRequest() : currentDt;
  }

  GEOS_LOG_LEVEL_RANK_0( 3, this->getName() << ": nextDt request is "  << nextDt );
  return nextDt;
}
template< typename POROMECHANICS_SOLVER >
void HydrofractureSolver< POROMECHANICS_SOLVER >::setUpDflux_dApertureMatrix( DomainPartition & domain,
                                                                              DofManager const & dofManager,
                                                                              CRSMatrix< real64, globalIndex > & localMatrix )
{
  std::unique_ptr< CRSMatrix< real64, localIndex > > &
  derivativeFluxResidual_dAperture = this->getRefDerivativeFluxResidual_dAperture();

  {
    localIndex numRows = 0;
    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
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
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
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

namespace
{
typedef HydrofractureSolver< SinglePhasePoromechanics > SinglePhaseHydrofracture;
// typedef HydrofractureSolver< MultiphasePoromechanics > MultiphaseHydrofracture;
REGISTER_CATALOG_ENTRY( SolverBase, SinglePhaseHydrofracture, string const &, Group * const )
// REGISTER_CATALOG_ENTRY( SolverBase, MultiphaseHydrofracture, string const &, Group * const )
}

} /* namespace geos */
