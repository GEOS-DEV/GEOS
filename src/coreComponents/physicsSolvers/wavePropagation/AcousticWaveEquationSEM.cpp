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
 * @file AcousticWaveEquationSEM.cpp
 */

#include "AcousticWaveEquationSEM.hpp"
#include "AcousticWaveEquationSEMKernel.hpp"

#include "finiteElement/FiniteElementDiscretization.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fieldSpecification/PerfectlyMatchedLayer.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/ElementType.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "WaveSolverUtils.hpp"

namespace geos
{

using namespace dataRepository;

AcousticWaveEquationSEM::AcousticWaveEquationSEM( const std::string & name,
                                                  Group * const parent ):
  WaveSolverBase( name,
                  parent )
{

  registerWrapper( viewKeyStruct::pressureNp1AtReceiversString(), &m_pressureNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Pressure value at each receiver for each timestep" );

}

AcousticWaveEquationSEM::~AcousticWaveEquationSEM()
{
  // TODO Auto-generated destructor stub
}

void AcousticWaveEquationSEM::initializePreSubGroups()
{
  WaveSolverBase::initializePreSubGroups();
}


void AcousticWaveEquationSEM::registerDataOnMesh( Group & meshBodies )
{
  WaveSolverBase::registerDataOnMesh( meshBodies );
  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    nodeManager.registerField< fields::Pressure_nm1,
                               fields::Pressure_n,
                               fields::Pressure_np1,
                               fields::PressureDoubleDerivative,
                               fields::ForcingRHS,
                               fields::MassVector,
                               fields::DampingVector,
                               fields::StiffnessVector,
                               fields::FreeSurfaceNodeIndicator >( getName() );

    /// register  PML auxiliary variables only when a PML is specified in the xml
    if( m_usePML )
    {
      nodeManager.registerField< fields::AuxiliaryVar1PML,
                                 fields::AuxiliaryVar2PML,
                                 fields::AuxiliaryVar3PML,
                                 fields::AuxiliaryVar4PML >( getName() );

      nodeManager.getField< fields::AuxiliaryVar1PML >().resizeDimension< 1 >( 3 );
      nodeManager.getField< fields::AuxiliaryVar2PML >().resizeDimension< 1 >( 3 );
    }

    FaceManager & faceManager = mesh.getFaceManager();
    faceManager.registerField< fields::FreeSurfaceFaceIndicator >( getName() );

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
    {
      subRegion.registerField< fields::MediumVelocity >( getName() );
      subRegion.registerField< fields::MediumDensity >( getName() );
      subRegion.registerField< fields::PartialGradient >( getName() );
    } );

  } );
}


void AcousticWaveEquationSEM::postProcessInput()
{

  WaveSolverBase::postProcessInput();

  m_pressureNp1AtReceivers.resize( m_nsamplesSeismoTrace, m_receiverCoordinates.size( 0 ) + 1 );
}

void AcousticWaveEquationSEM::precomputeSourceAndReceiverTerm( MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
{
  GEOS_MARK_FUNCTION;
  NodeManager const & nodeManager = mesh.getNodeManager();
  FaceManager const & faceManager = mesh.getFaceManager();

  arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const
  X32 = nodeManager.getField< fields::referencePosition32 >().toViewConst();

  arrayView2d< real64 const > const faceNormal  = faceManager.faceNormal();
  arrayView2d< real64 const > const faceCenter  = faceManager.faceCenter();


  arrayView2d< real64 const > const sourceCoordinates = m_sourceCoordinates.toViewConst();
  arrayView2d< localIndex > const sourceNodeIds = m_sourceNodeIds.toView();
  arrayView2d< real64 > const sourceConstants = m_sourceConstants.toView();
  arrayView1d< localIndex > const sourceIsAccessible = m_sourceIsAccessible.toView();
  sourceNodeIds.setValues< EXEC_POLICY >( -1 );
  sourceConstants.setValues< EXEC_POLICY >( -1 );
  sourceIsAccessible.zero();

  arrayView2d< real64 const > const receiverCoordinates = m_receiverCoordinates.toViewConst();
  arrayView2d< localIndex > const receiverNodeIds = m_receiverNodeIds.toView();
  arrayView2d< real64 > const receiverConstants = m_receiverConstants.toView();
  arrayView1d< localIndex > const receiverIsLocal = m_receiverIsLocal.toView();
  receiverNodeIds.setValues< EXEC_POLICY >( -1 );
  receiverConstants.setValues< EXEC_POLICY >( -1 );
  receiverIsLocal.zero();

  arrayView2d< real32 > const sourceValue = m_sourceValue.toView();
  real64 dt = 0;
  EventManager const & event = getGroupByPath< EventManager >( "/Problem/Events" );
  for( localIndex numSubEvent = 0; numSubEvent < event.numSubGroups(); ++numSubEvent )
  {
    EventBase const * subEvent = static_cast< EventBase const * >( event.getSubGroups()[numSubEvent] );
    if( subEvent->getEventName() == "/Solvers/" + getName() )
    {
      dt = subEvent->getReference< real64 >( EventBase::viewKeyStruct::forceDtString() );
    }
  }

  mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                        CellElementSubRegion & elementSubRegion )
  {
    GEOS_THROW_IF( elementSubRegion.getElementType() != ElementType::Hexahedron,
                   getDataContext() << ": Invalid type of element, the acoustic solver is designed for hexahedral meshes only (C3D8), using the SEM formulation",
                   InputError );

    arrayView2d< localIndex const > const elemsToFaces = elementSubRegion.faceList();
    arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();
    arrayView2d< real64 const > const elemCenter = elementSubRegion.getElementCenter();
    arrayView1d< integer const > const elemGhostRank = elementSubRegion.ghostRank();

    finiteElement::FiniteElementBase const &
    fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
    finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
    {
      using FE_TYPE = TYPEOFREF( finiteElement );

      localIndex const numFacesPerElem = elementSubRegion.numFacesPerElement();

      {
        GEOS_MARK_SCOPE( acousticWaveEquationSEMKernels::PrecomputeSourceAndReceiverKernel );
        acousticWaveEquationSEMKernels::
          PrecomputeSourceAndReceiverKernel::
          launch< EXEC_POLICY, FE_TYPE >
          ( elementSubRegion.size(),
          numFacesPerElem,
          X32,
          elemGhostRank,
          elemsToNodes,
          elemsToFaces,
          elemCenter,
          faceNormal,
          faceCenter,
          sourceCoordinates,
          sourceIsAccessible,
          sourceNodeIds,
          sourceConstants,
          receiverCoordinates,
          receiverIsLocal,
          receiverNodeIds,
          receiverConstants,
          sourceValue,
          dt,
          m_timeSourceFrequency,
          m_timeSourceDelay,
          m_rickerOrder );
      }
    } );
  } );
}

void AcousticWaveEquationSEM::addSourceToRightHandSide( integer const & cycleNumber, arrayView1d< real32 > const rhs )
{
  arrayView2d< localIndex const > const sourceNodeIds = m_sourceNodeIds.toViewConst();
  arrayView2d< real64 const > const sourceConstants   = m_sourceConstants.toViewConst();
  arrayView1d< localIndex const > const sourceIsAccessible = m_sourceIsAccessible.toViewConst();
  arrayView2d< real32 const > const sourceValue   = m_sourceValue.toViewConst();

  GEOS_THROW_IF( cycleNumber > sourceValue.size( 0 ),
                 getDataContext() << ": Too many steps compared to array size",
                 std::runtime_error );
  forAll< EXEC_POLICY >( sourceConstants.size( 0 ), [=] GEOS_HOST_DEVICE ( localIndex const isrc )
  {
    if( sourceIsAccessible[isrc] == 1 )
    {
      for( localIndex inode = 0; inode < sourceConstants.size( 1 ); ++inode )
      {
        real32 const localIncrement = sourceConstants[isrc][inode] * sourceValue[cycleNumber][isrc];
        RAJA::atomicAdd< ATOMIC_POLICY >( &rhs[sourceNodeIds[isrc][inode]], localIncrement );
      }
    }
  } );
}

void AcousticWaveEquationSEM::initializePostInitialConditionsPreSubGroups()
{
  GEOS_MARK_FUNCTION;
  {
    GEOS_MARK_SCOPE( WaveSolverBase::initializePostInitialConditionsPreSubGroups );
    WaveSolverBase::initializePostInitialConditionsPreSubGroups();
  }
  if( m_usePML )
  {
    AcousticWaveEquationSEM::initializePML();
  }

  DomainPartition & domain = getGroupByPath< DomainPartition >( "/Problem/domain" );

  applyFreeSurfaceBC( 0.0, domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    precomputeSourceAndReceiverTerm( mesh, regionNames );

    NodeManager & nodeManager = mesh.getNodeManager();
    FaceManager & faceManager = mesh.getFaceManager();
    ElementRegionManager & elemManager = mesh.getElemManager();

    /// get the array of indicators: 1 if the face is on the boundary; 0 otherwise
    arrayView1d< integer const > const & facesDomainBoundaryIndicator = faceManager.getDomainBoundaryIndicator();
    arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords = nodeManager.getField< fields::referencePosition32 >().toViewConst();

    /// get face to node map
    ArrayOfArraysView< localIndex const > const facesToNodes = faceManager.nodeList().toViewConst();

    // mass matrix to be computed in this function
    arrayView1d< real32 > const mass = nodeManager.getField< fields::MassVector >();
    {
      GEOS_MARK_SCOPE( mass_zero );
      mass.zero();
    }
    /// damping matrix to be computed for each dof in the boundary of the mesh
    arrayView1d< real32 > const damping = nodeManager.getField< fields::DampingVector >();
    {
      GEOS_MARK_SCOPE( damping_zero );
      damping.zero();
    }
    /// get array of indicators: 1 if face is on the free surface; 0 otherwise
    arrayView1d< localIndex const > const freeSurfaceFaceIndicator = faceManager.getField< fields::FreeSurfaceFaceIndicator >();

    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                CellElementSubRegion & elementSubRegion )
    {
      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

      arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes = elementSubRegion.nodeList();
      arrayView2d< localIndex const > const elemsToFaces = elementSubRegion.faceList();

      arrayView1d< real32 const > const velocity = elementSubRegion.getField< fields::MediumVelocity >();
      arrayView1d< real32 const > const density = elementSubRegion.getField< fields::MediumDensity >();

      /// Partial gradient if gradient as to be computed
      arrayView1d< real32 > grad = elementSubRegion.getField< fields::PartialGradient >();

      grad.zero();

      real64 dtCompute=0.0;

      finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        acousticWaveEquationSEMKernels::MassMatrixKernel< FE_TYPE > kernelM( finiteElement );
        kernelM.template launch< EXEC_POLICY, ATOMIC_POLICY >( elementSubRegion.size(),
                                                               nodeCoords,
                                                               elemsToNodes,
                                                               velocity,
                                                               density,
                                                               mass );
        {
          GEOS_MARK_SCOPE( DampingMatrixKernel );
          acousticWaveEquationSEMKernels::DampingMatrixKernel< FE_TYPE > kernelD( finiteElement );
          kernelD.template launch< EXEC_POLICY, ATOMIC_POLICY >( elementSubRegion.size(),
                                                                 nodeCoords,
                                                                 elemsToFaces,
                                                                 facesToNodes,
                                                                 facesDomainBoundaryIndicator,
                                                                 freeSurfaceFaceIndicator,
                                                                 velocity,
                                                                 density,
                                                                 damping );
        }

        if( m_preComputeDt==1 )
        {
          acousticWaveEquationSEMKernels::ComputeTimeStep< FE_TYPE > kernelT( finiteElement );

          dtCompute = kernelT.template launch< EXEC_POLICY, ATOMIC_POLICY >( elementSubRegion.size(),
                                                                             nodeManager.size(),
                                                                             nodeCoords,
                                                                             elemsToNodes,
                                                                             mass );

          real64 globaldt = MpiWrapper::min( dtCompute );

          printf( "dt=%f\n", globaldt );

          exit( 2 );

        }

      } );
    } );
  } );

  WaveSolverUtils::initTrace( "seismoTraceReceiver", getName(), m_outputSeismoTrace, m_receiverConstants.size( 0 ), m_receiverIsLocal );
}

//This function is only to give an easy accesss to the computation of the time-step for Pygeosx interface and avoid to exit the code when
// using Pygeosx

real64 AcousticWaveEquationSEM::computeTimeStep( real64 & dtOut )
{

  DomainPartition & domain = getGroupByPath< DomainPartition >( "/Problem/domain" );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {

    NodeManager & nodeManager = mesh.getNodeManager();
    ElementRegionManager & elemManager = mesh.getElemManager();

    /// get the array of indicators: 1 if the face is on the boundary; 0 otherwise
    arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords = nodeManager.getField< fields::referencePosition32 >().toViewConst();

    // mass matrix to be computed in this function
    arrayView1d< real32 > const mass = nodeManager.getField< fields::MassVector >();

    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                CellElementSubRegion & elementSubRegion )
    {
      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

      arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes = elementSubRegion.nodeList();

      real64 dtCompute=0.0;

      finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        acousticWaveEquationSEMKernels::ComputeTimeStep< FE_TYPE > kernelT( finiteElement );

        dtCompute = kernelT.template launch< EXEC_POLICY, ATOMIC_POLICY >( elementSubRegion.size(),
                                                                           nodeManager.size(),
                                                                           nodeCoords,
                                                                           elemsToNodes,
                                                                           mass );

        dtOut = MpiWrapper::min( dtCompute );



      } );
    } );
  } );
  return dtOut;
}



void AcousticWaveEquationSEM::applyFreeSurfaceBC( real64 time, DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  FunctionManager const & functionManager = FunctionManager::getInstance();

  FaceManager & faceManager = domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ).getFaceManager();
  NodeManager & nodeManager = domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ).getNodeManager();

  arrayView1d< real32 > const p_nm1 = nodeManager.getField< fields::Pressure_nm1 >();
  arrayView1d< real32 > const p_n = nodeManager.getField< fields::Pressure_n >();
  arrayView1d< real32 > const p_np1 = nodeManager.getField< fields::Pressure_np1 >();

  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  /// array of indicators: 1 if a face is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceFaceIndicator = faceManager.getField< fields::FreeSurfaceFaceIndicator >();

  /// array of indicators: 1 if a node is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceNodeIndicator = nodeManager.getField< fields::FreeSurfaceNodeIndicator >();

  //freeSurfaceFaceIndicator.zero();
  //freeSurfaceNodeIndicator.zero();

  fsManager.apply< FaceManager >( time,
                                  domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ),
                                  string( "FreeSurface" ),
                                  [&]( FieldSpecificationBase const & bc,
                                       string const &,
                                       SortedArrayView< localIndex const > const & targetSet,
                                       FaceManager &,
                                       string const & )
  {
    string const & functionName = bc.getFunctionName();

    if( functionName.empty() || functionManager.getGroup< FunctionBase >( functionName ).isFunctionOfTime() == 2 )
    {
      real64 const value = bc.getScale();

      for( localIndex i = 0; i < targetSet.size(); ++i )
      {
        localIndex const kf = targetSet[ i ];
        freeSurfaceFaceIndicator[kf] = 1;

        localIndex const numNodes = faceToNodeMap.sizeOfArray( kf );
        for( localIndex a=0; a < numNodes; ++a )
        {
          localIndex const dof = faceToNodeMap( kf, a );
          freeSurfaceNodeIndicator[dof] = 1;

          p_np1[dof] = value;
          p_n[dof]   = value;
          p_nm1[dof] = value;
        }
      }
    }
    else
    {
      GEOS_ERROR( "This option is not supported yet" );
    }
  } );
}

void AcousticWaveEquationSEM::initializePML()
{
  GEOS_MARK_FUNCTION;

  registerWrapper< parametersPML >( viewKeyStruct::parametersPMLString() ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Parameters needed to compute damping in the PML region" );

  parametersPML & param = getReference< parametersPML >( viewKeyStruct::parametersPMLString() );

  /// Get the default thicknesses and wave speeds in the PML regions from the PerfectlyMatchedLayer
  /// field specification parameters (from the xml)
  real32 minThicknessPML=0;
  real32 smallestXMinPML=0;
  real32 largestXMaxPML=0;
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  fsManager.forSubGroups< PerfectlyMatchedLayer >( [&] ( PerfectlyMatchedLayer const & fs )
  {
    param.xMinPML=fs.getMin();
    param.xMaxPML=fs.getMax();
    param.thicknessMinXYZPML=fs.getThicknessMinXYZ();
    param.thicknessMaxXYZPML=fs.getThicknessMaxXYZ();
    param.reflectivityPML = fs.getReflectivity();
    param.waveSpeedMinXYZPML=fs.getWaveSpeedMinXYZ();
    param.waveSpeedMaxXYZPML=fs.getWaveSpeedMaxXYZ();
    minThicknessPML=fs.minThickness;
    smallestXMinPML=fs.smallestXMin;
    largestXMaxPML=fs.largestXMax;
  } );

  /// Now compute the PML parameters above internally
  DomainPartition & domain = getGroupByPath< DomainPartition >( "/Problem/domain" );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {

    NodeManager & nodeManager = mesh.getNodeManager();
    /// WARNING: the array below is one of the PML auxiliary variables
    arrayView1d< real32 > const indicatorPML = nodeManager.getField< fields::AuxiliaryVar4PML >();
    arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const X32 = nodeManager.getField< fields::referencePosition32 >().toViewConst();
    indicatorPML.zero();

    real32 xInteriorMin[3]{};
    real32 xInteriorMax[3]{};
    real32 xGlobalMin[3]{};
    real32 xGlobalMax[3]{};
    real32 cMin[3]{};
    real32 cMax[3]{};
    integer counterMin[3]{};
    integer counterMax[3]{};

    /// Set a node-based flag in the PML regions
    /// WARNING: the array used as a flag is one of the PML
    /// auxiliary variables to save memory
    fsManager.apply< ElementSubRegionBase,
                     PerfectlyMatchedLayer >( 0.0,
                                              mesh,
                                              PerfectlyMatchedLayer::catalogName(),
                                              [&]( PerfectlyMatchedLayer const &,
                                                   string const &,
                                                   SortedArrayView< localIndex const > const & targetSet,
                                                   ElementSubRegionBase & subRegion,
                                                   string const & )

    {
      CellElementSubRegion::NodeMapType const & elemToNodes =
        subRegion.getReference< CellElementSubRegion::NodeMapType >( CellElementSubRegion::viewKeyStruct::nodeListString() );
      traits::ViewTypeConst< CellElementSubRegion::NodeMapType > const elemToNodesViewConst = elemToNodes.toViewConst();

      forAll< EXEC_POLICY >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const l )
      {
        localIndex const k = targetSet[ l ];
        localIndex const numNodesPerElem = elemToNodesViewConst[k].size();

        for( localIndex i=0; i<numNodesPerElem; ++i )
        {
          indicatorPML[elemToNodesViewConst[k][i]]=1.0;
        }
      } );
    } );


    /// find the interior and global coordinates limits
    RAJA::ReduceMin< parallelDeviceReduce, real32 > xMinGlobal( LvArray::NumericLimits< real32 >::max );
    RAJA::ReduceMin< parallelDeviceReduce, real32 > yMinGlobal( LvArray::NumericLimits< real32 >::max );
    RAJA::ReduceMin< parallelDeviceReduce, real32 > zMinGlobal( LvArray::NumericLimits< real32 >::max );
    RAJA::ReduceMax< parallelDeviceReduce, real32 > xMaxGlobal( -LvArray::NumericLimits< real32 >::max );
    RAJA::ReduceMax< parallelDeviceReduce, real32 > yMaxGlobal( -LvArray::NumericLimits< real32 >::max );
    RAJA::ReduceMax< parallelDeviceReduce, real32 > zMaxGlobal( -LvArray::NumericLimits< real32 >::max );
    RAJA::ReduceMin< parallelDeviceReduce, real32 > xMinInterior( LvArray::NumericLimits< real32 >::max );
    RAJA::ReduceMin< parallelDeviceReduce, real32 > yMinInterior( LvArray::NumericLimits< real32 >::max );
    RAJA::ReduceMin< parallelDeviceReduce, real32 > zMinInterior( LvArray::NumericLimits< real32 >::max );
    RAJA::ReduceMax< parallelDeviceReduce, real32 > xMaxInterior( -LvArray::NumericLimits< real32 >::max );
    RAJA::ReduceMax< parallelDeviceReduce, real32 > yMaxInterior( -LvArray::NumericLimits< real32 >::max );
    RAJA::ReduceMax< parallelDeviceReduce, real32 > zMaxInterior( -LvArray::NumericLimits< real32 >::max );

    forAll< EXEC_POLICY >( nodeManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      xMinGlobal.min( X32[a][0] );
      yMinGlobal.min( X32[a][1] );
      zMinGlobal.min( X32[a][2] );
      xMaxGlobal.max( X32[a][0] );
      yMaxGlobal.max( X32[a][1] );
      zMaxGlobal.max( X32[a][2] );
      if( !isZero( indicatorPML[a] - 1.0 ))
      {
        xMinInterior.min( X32[a][0] );
        yMinInterior.min( X32[a][1] );
        zMinInterior.min( X32[a][2] );
        xMaxInterior.max( X32[a][0] );
        yMaxInterior.max( X32[a][1] );
        zMaxInterior.max( X32[a][2] );
      }
    } );

    xGlobalMin[0] = xMinGlobal.get();
    xGlobalMin[1] = yMinGlobal.get();
    xGlobalMin[2] = zMinGlobal.get();
    xGlobalMax[0] = xMaxGlobal.get();
    xGlobalMax[1] = yMaxGlobal.get();
    xGlobalMax[2] = zMaxGlobal.get();
    xInteriorMin[0] = xMinInterior.get();
    xInteriorMin[1] = yMinInterior.get();
    xInteriorMin[2] = zMinInterior.get();
    xInteriorMax[0] = xMaxInterior.get();
    xInteriorMax[1] = yMaxInterior.get();
    xInteriorMax[2] = zMaxInterior.get();

    for( integer i=0; i<3; ++i )
    {
      xGlobalMin[i] = MpiWrapper::min( xGlobalMin[i] );
      xGlobalMax[i] = MpiWrapper::max( xGlobalMax[i] );
      xInteriorMin[i] = MpiWrapper::min( xInteriorMin[i] );
      xInteriorMax[i] = MpiWrapper::max( xInteriorMax[i] );
    }


    /// if the coordinates limits and PML thicknesses are not provided
    /// from the xml, replace them with the above
    for( integer i=0; i<3; ++i )
    {
      if( param.xMinPML[i]<smallestXMinPML )
        param.xMinPML[i] = xInteriorMin[i];
      if( param.xMaxPML[i]>largestXMaxPML )
        param.xMaxPML[i] = xInteriorMax[i];
      if( param.thicknessMinXYZPML[i]<0 )
        param.thicknessMinXYZPML[i] = xInteriorMin[i]-xGlobalMin[i];
      if( param.thicknessMaxXYZPML[i]<0 )
        param.thicknessMaxXYZPML[i] = xGlobalMax[i]-xInteriorMax[i];
    }

    /// Compute the average wave speeds in the PML regions internally
    /// using the actual velocity field
    fsManager.apply< ElementSubRegionBase,
                     PerfectlyMatchedLayer >( 0.0,
                                              mesh,
                                              PerfectlyMatchedLayer::catalogName(),
                                              [&]( PerfectlyMatchedLayer const &,
                                                   string const &,
                                                   SortedArrayView< localIndex const > const & targetSet,
                                                   ElementSubRegionBase & subRegion,
                                                   string const & )

    {
      CellElementSubRegion::NodeMapType const & elemToNodes =
        subRegion.getReference< CellElementSubRegion::NodeMapType >( CellElementSubRegion::viewKeyStruct::nodeListString() );
      traits::ViewTypeConst< CellElementSubRegion::NodeMapType > const elemToNodesViewConst = elemToNodes.toViewConst();
      arrayView1d< real32 const > const vel = subRegion.getReference< array1d< real32 > >( fields::MediumVelocity::key());
      finiteElement::FiniteElementBase const &
      fe = subRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

      real32 const xMin[3]{param.xMinPML[0], param.xMinPML[1], param.xMinPML[2]};
      real32 const xMax[3]{param.xMaxPML[0], param.xMaxPML[1], param.xMaxPML[2]};

      finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        acousticWaveEquationSEMKernels::
          waveSpeedPMLKernel< FE_TYPE > kernel( finiteElement );
        kernel.template launch< EXEC_POLICY, ATOMIC_POLICY >
          ( targetSet,
          X32,
          elemToNodesViewConst,
          vel,
          xMin,
          xMax,
          cMin,
          cMax,
          counterMin,
          counterMax );
      } );
    } );

    for( integer i=0; i<3; ++i )
    {
      cMin[i] = MpiWrapper::sum( cMin[i] );
      cMax[i] = MpiWrapper::sum( cMax[i] );
      counterMin[i] = MpiWrapper::sum( counterMin[i] );
      counterMax[i] = MpiWrapper::sum( counterMax[i] );
    }
    for( integer i=0; i<3; ++i )
    {
      cMin[i] /= std::max( 1, counterMin[i] );
      cMax[i] /= std::max( 1, counterMax[i] );
    }

    /// if the PML wave speeds are not provided from the xml
    /// replace them with the above
    for( integer i=0; i<3; ++i )
    {
      if( param.waveSpeedMinXYZPML[i]<0 )
        param.waveSpeedMinXYZPML[i] = cMin[i];
      if( param.waveSpeedMaxXYZPML[i]<0 )
        param.waveSpeedMaxXYZPML[i] = cMax[i];
    }

    /// add safeguards when PML thickness is negative or too small
    for( integer i=0; i<3; ++i )
    {
      if( param.thicknessMinXYZPML[i]<=minThicknessPML )
      {
        param.thicknessMinXYZPML[i]=LvArray::NumericLimits< real32 >::max;
        param.waveSpeedMinXYZPML[i]=0;
      }
      if( param.thicknessMaxXYZPML[i]<=minThicknessPML )
      {
        param.thicknessMaxXYZPML[i]=LvArray::NumericLimits< real32 >::max;
        param.waveSpeedMaxXYZPML[i]=0;
      }
    }

    /// WARNING: don't forget to reset the indicator to zero
    /// so it can be used by the PML application
    indicatorPML.zero();

    GEOS_LOG_LEVEL_RANK_0( 1,
                           "PML parameters are: \n"
                           << "\t inner boundaries xMin = "<<param.xMinPML<<"\n"
                           << "\t inner boundaries xMax = "<<param.xMaxPML<<"\n"
                           << "\t left, front, top max PML thicknesses  = "<<param.thicknessMinXYZPML<<"\n"
                           << "\t right, back, bottom max PML thicknesses  = "<<param.thicknessMaxXYZPML<<"\n"
                           << "\t left, front, top average wave speed  = "<<param.waveSpeedMinXYZPML<<"\n"
                           << "\t right, back, bottom average wave speed  = "<<param.waveSpeedMaxXYZPML<<"\n"
                           << "\t theoretical reflectivity = "<< param.reflectivityPML );

  } );
}



void AcousticWaveEquationSEM::applyPML( real64 const time, DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  parametersPML const & param = getReference< parametersPML >( viewKeyStruct::parametersPMLString() );

  /// Loop over the different mesh bodies; for wave propagation, there is only one mesh body
  /// which is the whole mesh
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {

    NodeManager & nodeManager = mesh.getNodeManager();

    /// Array views of the pressure p, PML auxiliary variables, and node coordinates
    arrayView1d< real32 const > const p_n = nodeManager.getField< fields::Pressure_n >();
    arrayView2d< real32 const > const v_n = nodeManager.getField< fields::AuxiliaryVar1PML >();
    arrayView2d< real32 > const grad_n = nodeManager.getField< fields::AuxiliaryVar2PML >();
    arrayView1d< real32 > const divV_n = nodeManager.getField< fields::AuxiliaryVar3PML >();
    arrayView1d< real32 const > const u_n = nodeManager.getField< fields::AuxiliaryVar4PML >();
    arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const X32 = nodeManager.getField< fields::referencePosition32 >().toViewConst();

    /// Select the subregions concerned by the PML (specified in the xml by the Field Specification)
    /// 'targetSet' contains the indices of the elements in a given subregion
    fsManager.apply< ElementSubRegionBase,
                     PerfectlyMatchedLayer >( time,
                                              mesh,
                                              PerfectlyMatchedLayer::catalogName(),
                                              [&]( PerfectlyMatchedLayer const &,
                                                   string const &,
                                                   SortedArrayView< localIndex const > const & targetSet,
                                                   ElementSubRegionBase & subRegion,
                                                   string const & )

    {

      /// Get the element to nodes mapping in the subregion
      CellElementSubRegion::NodeMapType const & elemToNodes =
        subRegion.getReference< CellElementSubRegion::NodeMapType >( CellElementSubRegion::viewKeyStruct::nodeListString() );

      /// Get a const ArrayView of the mapping above
      traits::ViewTypeConst< CellElementSubRegion::NodeMapType > const elemToNodesViewConst = elemToNodes.toViewConst();

      /// Array view of the wave speed
      arrayView1d< real32 const > const vel = subRegion.getReference< array1d< real32 > >( fields::MediumVelocity::key());

      /// Get the object needed to determine the type of the element in the subregion
      finiteElement::FiniteElementBase const &
      fe = subRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

      real32 xMin[3];
      real32 xMax[3];
      real32 dMin[3];
      real32 dMax[3];
      real32 cMin[3];
      real32 cMax[3];
      for( integer i=0; i<3; ++i )
      {
        xMin[i] = param.xMinPML[i];
        xMax[i] = param.xMaxPML[i];
        dMin[i] = param.thicknessMinXYZPML[i];
        dMax[i] = param.thicknessMaxXYZPML[i];
        cMin[i] = param.waveSpeedMinXYZPML[i];
        cMax[i] = param.waveSpeedMaxXYZPML[i];
      }
      real32 const r = param.reflectivityPML;

      /// Get the type of the elements in the subregion
      finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        /// apply the PML kernel
        acousticWaveEquationSEMKernels::
          PMLKernel< FE_TYPE > kernel( finiteElement );
        kernel.template launch< EXEC_POLICY, ATOMIC_POLICY >
          ( targetSet,
          X32,
          elemToNodesViewConst,
          vel,
          p_n,
          v_n,
          u_n,
          xMin,
          xMax,
          dMin,
          dMax,
          cMin,
          cMax,
          r,
          grad_n,
          divV_n );
      } );
    } );
  } );

}

real64 AcousticWaveEquationSEM::explicitStepForward( real64 const & time_n,
                                                     real64 const & dt,
                                                     integer cycleNumber,
                                                     DomainPartition & domain,
                                                     bool computeGradient )
{
  real64 dtOut = explicitStepInternal( time_n, dt, cycleNumber, domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(),
                                  [&] ( string const &,
                                        MeshLevel & mesh,
                                        arrayView1d< string const > const & GEOS_UNUSED_PARAM ( regionNames ) )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    arrayView1d< real32 > const p_nm1 = nodeManager.getField< fields::Pressure_nm1 >();
    arrayView1d< real32 > const p_n = nodeManager.getField< fields::Pressure_n >();
    arrayView1d< real32 > const p_np1 = nodeManager.getField< fields::Pressure_np1 >();

    if( computeGradient && cycleNumber >= 0 )
    {

      arrayView1d< real32 > const p_dt2 = nodeManager.getField< fields::PressureDoubleDerivative >();

      if( m_enableLifo )
      {
        if( !m_lifo )
        {
          int const rank = MpiWrapper::commRank( MPI_COMM_GEOSX );
          std::string lifoPrefix = GEOS_FMT( "lifo/rank_{:05}/pdt2_shot{:06}", rank, m_shotIndex );
          m_lifo = std::make_unique< LifoStorage< real32, localIndex > >( lifoPrefix, p_dt2, m_lifoOnDevice, m_lifoOnHost, m_lifoSize );
        }

        m_lifo->pushWait();
      }
      forAll< EXEC_POLICY >( nodeManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const nodeIdx )
      {
        p_dt2[nodeIdx] = (p_np1[nodeIdx] - 2*p_n[nodeIdx] + p_nm1[nodeIdx])/(dt*dt);
      } );

      if( m_enableLifo )
      {
        // Need to tell LvArray data is on GPU to avoir HtoD copy
        p_dt2.move( MemorySpace::cuda, false );
        m_lifo->pushAsync( p_dt2 );
      }
      else
      {
        GEOS_MARK_SCOPE ( DirectWrite );
        p_dt2.move( MemorySpace::host, false );
        int const rank = MpiWrapper::commRank( MPI_COMM_GEOSX );
        std::string fileName = GEOS_FMT( "lifo/rank_{:05}/pressuredt2_{:06}_{:08}.dat", rank, m_shotIndex, cycleNumber );
        int lastDirSeparator = fileName.find_last_of( "/\\" );
        std::string dirName = fileName.substr( 0, lastDirSeparator );
        if( string::npos != (size_t)lastDirSeparator && !directoryExists( dirName ))
        {
          makeDirsForPath( dirName );
        }

        //std::string fileName = GEOS_FMT( "pressuredt2_{:06}_{:08}_{:04}.dat", m_shotIndex, cycleNumber, rank );
        //const int fileDesc = open( fileName.c_str(), O_CREAT | O_WRONLY | O_DIRECT | O_TRUNC, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP |
        // S_IROTH | S_IWOTH );

        std::ofstream wf( fileName, std::ios::out | std::ios::binary );
        GEOS_THROW_IF( !wf,
                       getDataContext() << ": Could not open file "<< fileName << " for writing",
                       InputError );
        wf.write( (char *)&p_dt2[0], p_dt2.size()*sizeof( real32 ) );
        wf.close( );
        GEOS_THROW_IF( !wf.good(),
                       getDataContext() << ": An error occured while writing "<< fileName,
                       InputError );
      }

    }

    forAll< EXEC_POLICY >( nodeManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      p_nm1[a] = p_n[a];
      p_n[a]   = p_np1[a];
    } );
  } );

  return dtOut;
}


real64 AcousticWaveEquationSEM::explicitStepBackward( real64 const & time_n,
                                                      real64 const & dt,
                                                      integer cycleNumber,
                                                      DomainPartition & domain,
                                                      bool computeGradient )
{
  real64 dtOut = explicitStepInternal( time_n, dt, cycleNumber, domain );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(),
                                  [&] ( string const &,
                                        MeshLevel & mesh,
                                        arrayView1d< string const > const & regionNames )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    arrayView1d< real32 const > const mass = nodeManager.getField< fields::MassVector >();

    arrayView1d< real32 > const p_nm1 = nodeManager.getField< fields::Pressure_nm1 >();
    arrayView1d< real32 > const p_n = nodeManager.getField< fields::Pressure_n >();
    arrayView1d< real32 > const p_np1 = nodeManager.getField< fields::Pressure_np1 >();

    EventManager const & event = getGroupByPath< EventManager >( "/Problem/Events" );
    real64 const & maxTime = event.getReference< real64 >( EventManager::viewKeyStruct::maxTimeString() );
    int const maxCycle = int(round( maxTime / dt ));

    if( computeGradient && cycleNumber < maxCycle )
    {
      ElementRegionManager & elemManager = mesh.getElemManager();

      arrayView1d< real32 > const p_dt2 = nodeManager.getField< fields::PressureDoubleDerivative >();

      if( m_enableLifo )
      {
        m_lifo->pop( p_dt2 );
        if( m_lifo->empty() )
          delete m_lifo.release();
      }
      else
      {
        GEOS_MARK_SCOPE ( DirectRead );

        int const rank = MpiWrapper::commRank( MPI_COMM_GEOSX );
        std::string fileName = GEOS_FMT( "lifo/rank_{:05}/pressuredt2_{:06}_{:08}.dat", rank, m_shotIndex, cycleNumber );
        std::ifstream wf( fileName, std::ios::in | std::ios::binary );
        GEOS_THROW_IF( !wf,
                       getDataContext() << ": Could not open file "<< fileName << " for reading",
                       InputError );
        //std::string fileName = GEOS_FMT( "pressuredt2_{:06}_{:08}_{:04}.dat", m_shotIndex, cycleNumber, rank );
        //const int fileDesc = open( fileName.c_str(), O_RDONLY | O_DIRECT );
        //GEOS_ERROR_IF( fileDesc == -1,
        //                "Could not open file "<< fileName << " for reading: " << strerror( errno ) );
        // maybe better with registerTouch()
        p_dt2.move( MemorySpace::host, true );
        wf.read( (char *)&p_dt2[0], p_dt2.size()*sizeof( real32 ) );
        wf.close( );
        remove( fileName.c_str() );
      }
      elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                  CellElementSubRegion & elementSubRegion )
      {
        arrayView1d< real32 const > const velocity = elementSubRegion.getField< fields::MediumVelocity >();
        arrayView1d< real32 > grad = elementSubRegion.getField< fields::PartialGradient >();
        arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();
        constexpr localIndex numNodesPerElem = 8;
        arrayView1d< integer const > const elemGhostRank = elementSubRegion.ghostRank();
        GEOS_MARK_SCOPE ( updatePartialGradient );
        forAll< EXEC_POLICY >( elementSubRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const eltIdx )
        {
          if( elemGhostRank[eltIdx]<0 )
          {
            for( localIndex i = 0; i < numNodesPerElem; ++i )
            {
              localIndex nodeIdx = elemsToNodes[eltIdx][i];
              grad[eltIdx] += (-2/velocity[eltIdx]) * mass[nodeIdx]/8.0 * (p_dt2[nodeIdx] * p_n[nodeIdx]);
            }
          }
        } );
      } );
    }

    forAll< EXEC_POLICY >( nodeManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      p_nm1[a] = p_n[a];
      p_n[a]   = p_np1[a];
    } );
  } );

  return dtOut;
}

real64 AcousticWaveEquationSEM::explicitStepInternal( real64 const & time_n,
                                                      real64 const & dt,
                                                      integer cycleNumber,
                                                      DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  GEOS_LOG_RANK_0_IF( dt < epsilonLoc, "Warning! Value for dt: " << dt << "s is smaller than local threshold: " << epsilonLoc );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(),
                                  [&] ( string const &,
                                        MeshLevel & mesh,
                                        arrayView1d< string const > const & regionNames )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    arrayView1d< real32 const > const mass = nodeManager.getField< fields::MassVector >();
    arrayView1d< real32 const > const damping = nodeManager.getField< fields::DampingVector >();

    arrayView1d< real32 > const p_nm1 = nodeManager.getField< fields::Pressure_nm1 >();
    arrayView1d< real32 > const p_n = nodeManager.getField< fields::Pressure_n >();
    arrayView1d< real32 > const p_np1 = nodeManager.getField< fields::Pressure_np1 >();

    arrayView1d< localIndex const > const freeSurfaceNodeIndicator = nodeManager.getField< fields::FreeSurfaceNodeIndicator >();
    arrayView1d< real32 > const stiffnessVector = nodeManager.getField< fields::StiffnessVector >();
    arrayView1d< real32 > const rhs = nodeManager.getField< fields::ForcingRHS >();

    bool const usePML = m_usePML;

    auto kernelFactory = acousticWaveEquationSEMKernels::ExplicitAcousticSEMFactory( dt );

    finiteElement::
      regionBasedKernelApplication< EXEC_POLICY,
                                    constitutive::NullModel,
                                    CellElementSubRegion >( mesh,
                                                            regionNames,
                                                            getDiscretizationName(),
                                                            "",
                                                            kernelFactory );

    EventManager const & event = getGroupByPath< EventManager >( "/Problem/Events" );
    real64 const & minTime = event.getReference< real64 >( EventManager::viewKeyStruct::minTimeString() );
    integer const cycleForSource = int(round( -minTime / dt + cycleNumber ));

    addSourceToRightHandSide( cycleForSource, rhs );

    /// calculate your time integrators
    real64 const dt2 = dt*dt;

    if( !usePML )
    {
      GEOS_MARK_SCOPE ( updateP );
      forAll< EXEC_POLICY >( nodeManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
      {
        if( freeSurfaceNodeIndicator[a] != 1 )
        {
          p_np1[a] = p_n[a];
          p_np1[a] *= 2.0*mass[a];
          p_np1[a] -= (mass[a]-0.5*dt*damping[a])*p_nm1[a];
          p_np1[a] += dt2*(rhs[a]-stiffnessVector[a]);
          p_np1[a] /= mass[a]+0.5*dt*damping[a];
        }
      } );
    }
    else
    {
      parametersPML const & param = getReference< parametersPML >( viewKeyStruct::parametersPMLString() );
      arrayView2d< real32 > const v_n = nodeManager.getField< fields::AuxiliaryVar1PML >();
      arrayView2d< real32 > const grad_n = nodeManager.getField< fields::AuxiliaryVar2PML >();
      arrayView1d< real32 > const divV_n = nodeManager.getField< fields::AuxiliaryVar3PML >();
      arrayView1d< real32 > const u_n = nodeManager.getField< fields::AuxiliaryVar4PML >();
      arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const X32 = nodeManager.getField< fields::referencePosition32 >().toViewConst();

      real32 const xMin[ 3 ] = {param.xMinPML[0], param.xMinPML[1], param.xMinPML[2]};
      real32 const xMax[ 3 ] = {param.xMaxPML[0], param.xMaxPML[1], param.xMaxPML[2]};
      real32 const dMin[ 3 ] = {param.thicknessMinXYZPML[0], param.thicknessMinXYZPML[1], param.thicknessMinXYZPML[2]};
      real32 const dMax[ 3 ] = {param.thicknessMaxXYZPML[0], param.thicknessMaxXYZPML[1], param.thicknessMaxXYZPML[2]};
      real32 const cMin[ 3 ] = {param.waveSpeedMinXYZPML[0], param.waveSpeedMinXYZPML[1], param.waveSpeedMinXYZPML[2]};
      real32 const cMax[ 3 ] = {param.waveSpeedMaxXYZPML[0], param.waveSpeedMaxXYZPML[1], param.waveSpeedMaxXYZPML[2]};
      real32 const r = param.reflectivityPML;

      /// apply the main function to update some of the PML auxiliary variables
      /// Compute (divV) and (B.pressureGrad - C.auxUGrad) vectors for the PML region
      applyPML( time_n, domain );

      GEOS_MARK_SCOPE ( updatePWithPML );
      forAll< EXEC_POLICY >( nodeManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
      {
        if( freeSurfaceNodeIndicator[a] != 1 )
        {
          real32 sigma[3];
          real32 xLocal[ 3 ];

          for( integer i=0; i<3; ++i )
          {
            xLocal[i] = X32[a][i];
          }

          acousticWaveEquationSEMKernels::PMLKernelHelper::computeDampingProfilePML(
            xLocal,
            xMin,
            xMax,
            dMin,
            dMax,
            cMin,
            cMax,
            r,
            sigma );

          real32 const alpha = sigma[0] + sigma[1] + sigma[2];

          p_np1[a] = dt2*( (rhs[a] - stiffnessVector[a])/mass[a] - divV_n[a])
                     - (1 - 0.5*alpha*dt)*p_nm1[a]
                     + 2*p_n[a];

          p_np1[a] = p_np1[a] / (1 + 0.5*alpha*dt);

          for( integer i=0; i<3; ++i )
          {
            v_n[a][i] = (1 - dt*sigma[i])*v_n[a][i] - dt*grad_n[a][i];
          }
          u_n[a] += dt*p_n[a];
        }
      } );
    }

    /// synchronize pressure fields
    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addFields( FieldLocation::Node, { fields::Pressure_np1::key() } );

    if( usePML )
    {
      fieldsToBeSync.addFields( FieldLocation::Node, {
          fields::AuxiliaryVar1PML::key(),
          fields::AuxiliaryVar4PML::key() } );
    }

    CommunicationTools & syncFields = CommunicationTools::getInstance();
    syncFields.synchronizeFields( fieldsToBeSync,
                                  mesh,
                                  domain.getNeighbors(),
                                  true );

    /// compute the seismic traces since last step.
    arrayView2d< real32 > const pReceivers = m_pressureNp1AtReceivers.toView();
    computeAllSeismoTraces( time_n, dt, p_np1, p_n, pReceivers );
    incrementIndexSeismoTrace( time_n );

    /// prepare next step
    forAll< EXEC_POLICY >( nodeManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      stiffnessVector[a] = 0.0;
      rhs[a] = 0.0;
    } );

    if( usePML )
    {
      arrayView2d< real32 > const grad_n = nodeManager.getField< fields::AuxiliaryVar2PML >();
      arrayView1d< real32 > const divV_n = nodeManager.getField< fields::AuxiliaryVar3PML >();
      grad_n.zero();
      divV_n.zero();
    }

  } );

  return dt;
}

void AcousticWaveEquationSEM::cleanup( real64 const time_n,
                                       integer const cycleNumber,
                                       integer const eventCounter,
                                       real64 const eventProgress,
                                       DomainPartition & domain )
{
  // call the base class cleanup (for reporting purposes)
  SolverBase::cleanup( time_n, cycleNumber, eventCounter, eventProgress, domain );

  // compute the remaining seismic traces, if needed
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    NodeManager & nodeManager = mesh.getNodeManager();
    arrayView1d< real32 const > const p_n = nodeManager.getField< fields::Pressure_n >();
    arrayView1d< real32 const > const p_np1 = nodeManager.getField< fields::Pressure_np1 >();
    arrayView2d< real32 > const pReceivers = m_pressureNp1AtReceivers.toView();
    computeAllSeismoTraces( time_n, 0.0, p_np1, p_n, pReceivers );

    WaveSolverUtils::writeSeismoTrace( "seismoTraceReceiver", getName(), m_outputSeismoTrace, m_receiverConstants.size( 0 ),
                                       m_receiverIsLocal, m_nsamplesSeismoTrace, pReceivers );
  } );
}

REGISTER_CATALOG_ENTRY( SolverBase, AcousticWaveEquationSEM, string const &, dataRepository::Group * const )

} /* namespace geos */
