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

#include "dataRepository/KeyNames.hpp"
#include "finiteElement/FiniteElementDiscretization.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fieldSpecification/PerfectlyMatchedLayer.hpp"
#include "mainInterface/ProblemManager.hpp"
//#include "mesh/CellBlock.hpp"
#include "mesh/ElementType.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"

namespace geosx
{

using namespace dataRepository;

AcousticWaveEquationSEM::AcousticWaveEquationSEM( const std::string & name,
                                                  Group * const parent ):
  WaveSolverBase( name,
                  parent )
{

  registerWrapper( viewKeyStruct::sourceNodeIdsString(), &m_sourceNodeIds ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Indices of the nodes (in the right order) for each source point" );

  registerWrapper( viewKeyStruct::sourceConstantsString(), &m_sourceConstants ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Constant part of the source for the nodes listed in m_sourceNodeIds" );

  registerWrapper( viewKeyStruct::sourceIsAccessibleString(), &m_sourceIsAccessible ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Flag that indicates whether the source is accessible to this MPI rank" );

  registerWrapper( viewKeyStruct::receiverNodeIdsString(), &m_receiverNodeIds ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Indices of the nodes (in the right order) for each receiver point" );

  registerWrapper( viewKeyStruct::receiverConstantsString(), &m_receiverConstants ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Constant part of the receiver for the nodes listed in m_receiverNodeIds" );

  registerWrapper( viewKeyStruct::receiverIsLocalString(), &m_receiverIsLocal ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Flag that indicates whether the receiver is local to this MPI rank" );

  registerWrapper( viewKeyStruct::pressureNp1AtReceiversString(), &m_pressureNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Pressure value at each receiver for each timestep" );

}

AcousticWaveEquationSEM::~AcousticWaveEquationSEM()
{
  // TODO Auto-generated destructor stub
}

localIndex AcousticWaveEquationSEM::getNumNodesPerElem()
{
  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();

  FiniteElementDiscretizationManager const &
  feDiscretizationManager = numericalMethodManager.getFiniteElementDiscretizationManager();

  FiniteElementDiscretization const * const
  feDiscretization = feDiscretizationManager.getGroupPointer< FiniteElementDiscretization >( m_discretizationName );
  GEOSX_THROW_IF( feDiscretization == nullptr,
                  getName() << ": FE discretization not found: " << m_discretizationName,
                  InputError );

  localIndex numNodesPerElem = 0;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(),
                                  [&]( string const &,
                                       MeshLevel const & mesh,
                                       arrayView1d< string const > const & regionNames )
  {
    ElementRegionManager const & elemManager = mesh.getElemManager();
    elemManager.forElementRegions( regionNames,
                                   [&] ( localIndex const,
                                         ElementRegionBase const & elemRegion )
    {
      elemRegion.forElementSubRegions( [&]( ElementSubRegionBase const & elementSubRegion )
      {
        finiteElement::FiniteElementBase const &
        fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
        localIndex const numSupportPoints = fe.getNumSupportPoints();
        if( numSupportPoints > numNodesPerElem )
        {
          numNodesPerElem = numSupportPoints;
        }
      } );
    } );


  } );
  return numNodesPerElem;
}

void AcousticWaveEquationSEM::initializePreSubGroups()
{
  WaveSolverBase::initializePreSubGroups();

  localIndex numNodesPerElem = getNumNodesPerElem();

  localIndex const numSourcesGlobal = m_sourceCoordinates.size( 0 );
  m_sourceNodeIds.resize( numSourcesGlobal, numNodesPerElem );
  m_sourceConstants.resize( numSourcesGlobal, numNodesPerElem );
  //m_sourceIsLocal.resize( numSourcesGlobal );

  localIndex const numReceiversGlobal = m_receiverCoordinates.size( 0 );
  m_receiverNodeIds.resize( numReceiversGlobal, numNodesPerElem );
  m_receiverConstants.resize( numReceiversGlobal, numNodesPerElem );
  m_receiverIsLocal.resize( numReceiversGlobal );

  m_pressureNp1AtReceivers.resizeDimension< 1 >( numReceiversGlobal );



}


void AcousticWaveEquationSEM::registerDataOnMesh( Group & meshBodies )
{

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    nodeManager.registerExtrinsicData< extrinsicMeshData::Pressure_nm1,
                                       extrinsicMeshData::Pressure_n,
                                       extrinsicMeshData::Pressure_np1,
                                       extrinsicMeshData::PressureDoubleDerivative,
                                       extrinsicMeshData::ForcingRHS,
                                       extrinsicMeshData::MassVector,
                                       extrinsicMeshData::DampingVector,
                                       extrinsicMeshData::StiffnessVector,
                                       extrinsicMeshData::FreeSurfaceNodeIndicator >( this->getName() );

    /// register  PML auxiliary variables only when a PML is specified in the xml
    if( m_usePML )
    {
      nodeManager.registerExtrinsicData< extrinsicMeshData::AuxiliaryVar1PML,
                                         extrinsicMeshData::AuxiliaryVar2PML,
                                         extrinsicMeshData::AuxiliaryVar3PML,
                                         extrinsicMeshData::AuxiliaryVar4PML >( this->getName() );

      nodeManager.getExtrinsicData< extrinsicMeshData::AuxiliaryVar1PML >().resizeDimension< 1 >( 3 );
      nodeManager.getExtrinsicData< extrinsicMeshData::AuxiliaryVar2PML >().resizeDimension< 1 >( 3 );
    }

    FaceManager & faceManager = mesh.getFaceManager();
    faceManager.registerExtrinsicData< extrinsicMeshData::FreeSurfaceFaceIndicator >( this->getName() );

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
    {
      subRegion.registerExtrinsicData< extrinsicMeshData::MediumVelocity >( this->getName() );
      subRegion.registerExtrinsicData< extrinsicMeshData::PartialGradient >( this->getName() );
    } );

  } );
}


void AcousticWaveEquationSEM::postProcessInput()
{

  WaveSolverBase::postProcessInput();
  GEOSX_THROW_IF( m_sourceCoordinates.size( 1 ) != 3,
                  "Invalid number of physical coordinates for the sources",
                  InputError );

  GEOSX_THROW_IF( m_receiverCoordinates.size( 1 ) != 3,
                  "Invalid number of physical coordinates for the receivers",
                  InputError );

  EventManager const & event = this->getGroupByPath< EventManager >( "/Problem/Events" );
  real64 const & maxTime = event.getReference< real64 >( EventManager::viewKeyStruct::maxTimeString() );
  real64 dt = 0;
  for( localIndex numSubEvent = 0; numSubEvent < event.numSubGroups(); ++numSubEvent )
  {
    EventBase const * subEvent = static_cast< EventBase const * >( event.getSubGroups()[numSubEvent] );
    if( subEvent->getEventName() == "/Solvers/" + this->getName() )
    {
      dt = subEvent->getReference< real64 >( EventBase::viewKeyStruct::forceDtString() );
    }
  }

  GEOSX_THROW_IF( dt < epsilonLoc*maxTime, "Value for dt: " << dt <<" is smaller than local threshold: " << epsilonLoc, std::runtime_error );

  if( m_dtSeismoTrace > 0 )
  {
    m_nsamplesSeismoTrace = int( maxTime / m_dtSeismoTrace) + 1;
  }
  else
  {
    m_nsamplesSeismoTrace = 0;
  }
  localIndex const nsamples = int(maxTime/dt) + 1;


  localIndex numNodesPerElem = getNumNodesPerElem();

  localIndex const numSourcesGlobal = m_sourceCoordinates.size( 0 );
  m_sourceNodeIds.resize( numSourcesGlobal, numNodesPerElem );
  m_sourceConstants.resize( numSourcesGlobal, numNodesPerElem );
  m_sourceIsAccessible.resize( numSourcesGlobal );

  localIndex const numReceiversGlobal = m_receiverCoordinates.size( 0 );
  m_receiverNodeIds.resize( numReceiversGlobal, numNodesPerElem );
  m_receiverConstants.resize( numReceiversGlobal, numNodesPerElem );
  m_receiverIsLocal.resize( numReceiversGlobal );

  m_pressureNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal );
  m_sourceValue.resize( nsamples, numSourcesGlobal );

}

void AcousticWaveEquationSEM::precomputeSourceAndReceiverTerm( MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
{
  NodeManager const & nodeManager = mesh.getNodeManager();
  FaceManager const & faceManager = mesh.getFaceManager();

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const
  X = nodeManager.referencePosition().toViewConst();
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

  real32 const timeSourceFrequency = this->m_timeSourceFrequency;
  localIndex const rickerOrder = this->m_rickerOrder;
  arrayView2d< real32 > const sourceValue = m_sourceValue.toView();
  real64 dt = 0;
  EventManager const & event = this->getGroupByPath< EventManager >( "/Problem/Events" );
  for( localIndex numSubEvent = 0; numSubEvent < event.numSubGroups(); ++numSubEvent )
  {
    EventBase const * subEvent = static_cast< EventBase const * >( event.getSubGroups()[numSubEvent] );
    if( subEvent->getEventName() == "/Solvers/" + this->getName() )
    {
      dt = subEvent->getReference< real64 >( EventBase::viewKeyStruct::forceDtString() );
    }
  }

  // arrayView2d< real64 >  center;

  // arrayView2d< localIndex >  tmp;

  // mesh0.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
  //                                                                                        CellElementSubRegion & elementSubRegion0 )
  // {
  //   center = elementSubRegion0.getElementCenter();
  //   tmp = elementSubRegion0.faceList();
  // } );

  mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                        CellElementSubRegion & elementSubRegion )
  {
    GEOSX_THROW_IF( elementSubRegion.getElementType() != ElementType::Hexahedron,
                    "Invalid type of element, the acoustic solver is designed for hexahedral meshes only (C3D8) ",
                    InputError );

    arrayView2d< localIndex const > const elemsToFaces = elementSubRegion.faceList();
    arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();
    arrayView2d< real64 const > const elemCenter = elementSubRegion.getElementCenter();
    arrayView1d< integer const > const elemGhostRank = elementSubRegion.ghostRank();

    finiteElement::FiniteElementBase const &
    fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
    finiteElement::dispatch3D( fe,
                               [&]
                                 ( auto const finiteElement )
    {
      using FE_TYPE = TYPEOFREF( finiteElement );

      constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
      localIndex const numFacesPerElem = elementSubRegion.numFacesPerElement();

      acousticWaveEquationSEMKernels::
        PrecomputeSourceAndReceiverKernel::
        launch< EXEC_POLICY, FE_TYPE >
        ( elementSubRegion.size(),
        numNodesPerElem,
        numFacesPerElem,
        X,
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
        timeSourceFrequency,
        rickerOrder );
    } );
  } );

}


void AcousticWaveEquationSEM::addSourceToRightHandSide( integer const & cycleNumber, arrayView1d< real32 > const rhs )
{
  arrayView2d< localIndex const > const sourceNodeIds = m_sourceNodeIds.toViewConst();
  arrayView2d< real64 const > const sourceConstants   = m_sourceConstants.toViewConst();
  arrayView1d< localIndex const > const sourceIsAccessible = m_sourceIsAccessible.toViewConst();
  arrayView2d< real32 const > const sourceValue   = m_sourceValue.toViewConst();

  GEOSX_THROW_IF( cycleNumber > sourceValue.size( 0 ), "Too many steps compared to array size", std::runtime_error );
  forAll< EXEC_POLICY >( sourceConstants.size( 0 ), [=] GEOSX_HOST_DEVICE ( localIndex const isrc )
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

void AcousticWaveEquationSEM::computeSeismoTrace( real64 const time_n,
                                                  real64 const dt,
                                                  real64 const timeSeismo,
                                                  localIndex iSeismo,
                                                  arrayView1d< real32 const > const var_np1,
                                                  arrayView1d< real32 const > const var_n,
                                                  arrayView2d< real32 > varAtReceivers )
{
  real64 const time_np1 = time_n + dt;
  arrayView2d< localIndex const > const receiverNodeIds = m_receiverNodeIds.toViewConst();
  arrayView2d< real64 const > const receiverConstants   = m_receiverConstants.toViewConst();
  arrayView1d< localIndex const > const receiverIsLocal = m_receiverIsLocal.toViewConst();


  real32 const a1 = (std::abs( dt ) < epsilonLoc) ? 1.0 : (time_np1 - timeSeismo)/dt;
  real32 const a2 = 1.0 - a1;

  if( m_nsamplesSeismoTrace > 0 )
  {
    forAll< EXEC_POLICY >( receiverConstants.size( 0 ), [=] GEOSX_HOST_DEVICE ( localIndex const ircv )
    {
      if( receiverIsLocal[ircv] == 1 )
      {
        varAtReceivers[iSeismo][ircv] = 0.0;
        real32 vtmp_np1 = 0.0;
        real32 vtmp_n = 0.0;
        for( localIndex inode = 0; inode < receiverConstants.size( 1 ); ++inode )
        {
          vtmp_np1 += var_np1[receiverNodeIds[ircv][inode]] * receiverConstants[ircv][inode];
          vtmp_n += var_n[receiverNodeIds[ircv][inode]] * receiverConstants[ircv][inode];
        }
        // linear interpolation between the pressure value at time_n and time_(n+1)
        varAtReceivers[iSeismo][ircv] = a1*vtmp_n + a2*vtmp_np1;
      }
    } );
  }

  // TODO DEBUG: the following output is only temporary until our wave propagation kernels are finalized.
  // Output will then only be done via the previous code.
  if( iSeismo == m_nsamplesSeismoTrace - 1 )
  {
    forAll< serialPolicy >( receiverConstants.size( 0 ), [=] ( localIndex const ircv )
    {
      if( this->m_outputSeismoTrace == 1 )
      {
        if( receiverIsLocal[ircv] == 1 )
        {
          // Note: this "manual" output to file is temporary
          //       It should be removed as soon as we can use TimeHistory to output data not registered on the mesh
          // TODO: remove saveSeismo and replace with TimeHistory
          std::ofstream f( GEOSX_FMT( "seismoTraceReceiver{:03}.txt", ircv ), std::ios::app );
          for( localIndex iSample = 0; iSample < m_nsamplesSeismoTrace; ++iSample )
          {
            f << iSample << " " << varAtReceivers[iSample][ircv] << std::endl;
          }
          f.close();
        }
      }
    } );
  }

}

/// Use for now until we get the same functionality in TimeHistory
/// TODO: move implementation into WaveSolverBase
void AcousticWaveEquationSEM::saveSeismo( localIndex const iSeismo, real32 const val, string const & filename )
{
  std::ofstream f( filename, std::ios::app );
  f<< iSeismo << " " << val << std::endl;
  f.close();
}

void AcousticWaveEquationSEM::initializePostInitialConditionsPreSubGroups()
{

  WaveSolverBase::initializePostInitialConditionsPreSubGroups();
  if( m_usePML )
  {
    AcousticWaveEquationSEM::initializePML();
  }

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  real64 const time = 0.0;
  applyFreeSurfaceBC( time, domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    precomputeSourceAndReceiverTerm( mesh, regionNames );

    NodeManager & nodeManager = mesh.getNodeManager();
    FaceManager & faceManager = mesh.getFaceManager();
    //TODO: Please make a shallow copy of face normal in the constructo for the new mesh level
    FaceManager & faceManager0 = this->getGroupByPath< FaceManager >( "/Problem/domain/MeshBodies/mesh/meshLevels/Level0/faceManager" );

    /// get the array of indicators: 1 if the face is on the boundary; 0 otherwise
    arrayView1d< integer > const & facesDomainBoundaryIndicator = faceManager.getDomainBoundaryIndicator();
    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition().toViewConst();

    /// get table containing all the face normals
    arrayView2d< real64 const > const faceNormal  = faceManager0.faceNormal();
    ArrayOfArraysView< localIndex const > const facesToNodes = faceManager.nodeList().toViewConst();

    // mass matrix to be computed in this function
    arrayView1d< real32 > const mass = nodeManager.getExtrinsicData< extrinsicMeshData::MassVector >();
    mass.zero();
    /// damping matrix to be computed for each dof in the boundary of the mesh
    arrayView1d< real32 > const damping = nodeManager.getExtrinsicData< extrinsicMeshData::DampingVector >();
    damping.zero();

    /// get array of indicators: 1 if face is on the free surface; 0 otherwise
    arrayView1d< localIndex const > const freeSurfaceFaceIndicator = faceManager.getExtrinsicData< extrinsicMeshData::FreeSurfaceFaceIndicator >();

    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          CellElementSubRegion & elementSubRegion )
    {

      arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes = elementSubRegion.nodeList();
      arrayView2d< localIndex const > const elemsToFaces = elementSubRegion.faceList();
      arrayView1d< real32 const > const velocity = elementSubRegion.getExtrinsicData< extrinsicMeshData::MediumVelocity >();

      /// Partial gradient if gradient as to be computed
      arrayView1d< real64 > grad = elementSubRegion.getExtrinsicData< extrinsicMeshData::PartialGradient >();
      grad.zero();

      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
      finiteElement::dispatch3D( fe,
                                 [&]
                                   ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        localIndex const numFacesPerElem = elementSubRegion.numFacesPerElement();
        localIndex const numNodesPerFace = facesToNodes.sizeOfArray( 0 );

        acousticWaveEquationSEMKernels::MassAndDampingMatrixKernel< FE_TYPE > kernel( finiteElement );

        kernel.template launch< EXEC_POLICY, ATOMIC_POLICY >( elementSubRegion.size(),
                                                              numFacesPerElem,
                                                              numNodesPerFace,
                                                              X,
                                                              elemsToNodes,
                                                              elemsToFaces,
                                                              facesToNodes,
                                                              facesDomainBoundaryIndicator,
                                                              freeSurfaceFaceIndicator,
                                                              faceNormal,
                                                              velocity,
                                                              mass,
                                                              damping );
      } );
    } );
  } );

}


void AcousticWaveEquationSEM::applyFreeSurfaceBC( real64 time, DomainPartition & domain )
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  FunctionManager const & functionManager = FunctionManager::getInstance();

  FaceManager & faceManager = domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ).getFaceManager();
  NodeManager & nodeManager = domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ).getNodeManager();

  arrayView1d< real32 > const p_nm1 = nodeManager.getExtrinsicData< extrinsicMeshData::Pressure_nm1 >();
  arrayView1d< real32 > const p_n = nodeManager.getExtrinsicData< extrinsicMeshData::Pressure_n >();
  arrayView1d< real32 > const p_np1 = nodeManager.getExtrinsicData< extrinsicMeshData::Pressure_np1 >();

  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  /// array of indicators: 1 if a face is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceFaceIndicator = faceManager.getExtrinsicData< extrinsicMeshData::FreeSurfaceFaceIndicator >();

  /// array of indicators: 1 if a node is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceNodeIndicator = nodeManager.getExtrinsicData< extrinsicMeshData::FreeSurfaceNodeIndicator >();

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
      GEOSX_ERROR( "This option is not supported yet" );
    }
  } );
}

void AcousticWaveEquationSEM::initializePML()
{
  GEOSX_MARK_FUNCTION;

  registerWrapper< parametersPML >( viewKeyStruct::parametersPMLString() ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Parameters needed to compute damping in the PML region" );

  parametersPML & param = getReference< parametersPML >( viewKeyStruct::parametersPMLString() );

  /// Get the default thicknesses and wave speeds in the PML regions from the PerfectlyMatchedLayer
  /// field specification parameters (from the xml)
  real64 minThicknessPML=0;
  real64 smallestXMinPML=0;
  real64 largestXMaxPML=0;
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
  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {

    NodeManager & nodeManager = mesh.getNodeManager();
    /// WARNING: the array below is one of the PML auxiliary variables
    arrayView1d< real32 > const indicatorPML = nodeManager.getExtrinsicData< extrinsicMeshData::AuxiliaryVar4PML >();
    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition().toViewConst();
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

      forAll< EXEC_POLICY >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const l )
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

    forAll< EXEC_POLICY >( nodeManager.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      xMinGlobal.min( X[a][0] );
      yMinGlobal.min( X[a][1] );
      zMinGlobal.min( X[a][2] );
      xMaxGlobal.max( X[a][0] );
      yMaxGlobal.max( X[a][1] );
      zMaxGlobal.max( X[a][2] );
      if( !isZero( indicatorPML[a] - 1.0 ))
      {
        xMinInterior.min( X[a][0] );
        yMinInterior.min( X[a][1] );
        zMinInterior.min( X[a][2] );
        xMaxInterior.max( X[a][0] );
        yMaxInterior.max( X[a][1] );
        zMaxInterior.max( X[a][2] );
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
      arrayView1d< real32 const > const vel = subRegion.getReference< array1d< real32 > >( extrinsicMeshData::MediumVelocity::key());
      finiteElement::FiniteElementBase const &
      fe = subRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

      real32 const xMin[3]{param.xMinPML[0], param.xMinPML[1], param.xMinPML[2]};
      real32 const xMax[3]{param.xMaxPML[0], param.xMaxPML[1], param.xMaxPML[2]};

      finiteElement::dispatch3D( fe,
                                 [&]
                                   ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        acousticWaveEquationSEMKernels::
          waveSpeedPMLKernel< FE_TYPE > kernel( finiteElement );
        kernel.template launch< EXEC_POLICY, ATOMIC_POLICY >
          ( targetSet,
          X,
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

    GEOSX_LOG_LEVEL_RANK_0( 1,
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
  GEOSX_MARK_FUNCTION;

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
    arrayView1d< real32 const > const p_n = nodeManager.getExtrinsicData< extrinsicMeshData::Pressure_n >();
    arrayView2d< real32 const > const v_n = nodeManager.getExtrinsicData< extrinsicMeshData::AuxiliaryVar1PML >();
    arrayView2d< real32 > const grad_n = nodeManager.getExtrinsicData< extrinsicMeshData::AuxiliaryVar2PML >();
    arrayView1d< real32 > const divV_n = nodeManager.getExtrinsicData< extrinsicMeshData::AuxiliaryVar3PML >();
    arrayView1d< real32 const > const u_n = nodeManager.getExtrinsicData< extrinsicMeshData::AuxiliaryVar4PML >();
    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition().toViewConst();

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
      arrayView1d< real32 const > const vel = subRegion.getReference< array1d< real32 > >( extrinsicMeshData::MediumVelocity::key());

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
      finiteElement::dispatch3D( fe,
                                 [&]
                                   ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        /// apply the PML kernel
        acousticWaveEquationSEMKernels::
          PMLKernel< FE_TYPE > kernel( finiteElement );
        kernel.template launch< EXEC_POLICY, ATOMIC_POLICY >
          ( targetSet,
          X,
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
                                        arrayView1d< string const > const & GEOSX_UNUSED_PARAM ( regionNames ) )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    arrayView1d< real64 > const p_nm1 = nodeManager.getExtrinsicData< extrinsicMeshData::Pressure_nm1 >();
    arrayView1d< real64 > const p_n = nodeManager.getExtrinsicData< extrinsicMeshData::Pressure_n >();
    arrayView1d< real64 > const p_np1 = nodeManager.getExtrinsicData< extrinsicMeshData::Pressure_np1 >();

    if( computeGradient )
    {

      arrayView1d< real64 > const p_dt2 = nodeManager.getExtrinsicData< extrinsicMeshData::PressureDoubleDerivative >();

      forAll< EXEC_POLICY >( nodeManager.size(), [=] GEOSX_HOST_DEVICE ( localIndex const nodeIdx )
      {
        p_dt2[nodeIdx] = (p_np1[nodeIdx] - 2*p_n[nodeIdx] + p_nm1[nodeIdx])/(dt*dt);
      } );

      p_dt2.move( MemorySpace::host, false );
      int const rank = MpiWrapper::commRank( MPI_COMM_GEOSX );
      std::string fileName = GEOSX_FMT( "pressuredt2_{:06}_{:08}_{:04}.dat", m_shotIndex, cycleNumber, rank );
      std::ofstream wf( fileName, std::ios::out | std::ios::binary );
      GEOSX_THROW_IF( !wf,
                      "Could not open file "<< fileName << " for writting",
                      InputError );
      wf.write( (char *)&p_dt2[0], p_dt2.size()*sizeof( real64 ) );
      wf.close();
      GEOSX_THROW_IF( !wf.good(),
                      "An error occured while writting "<< fileName,
                      InputError );

    }

    forAll< EXEC_POLICY >( nodeManager.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
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

    arrayView1d< real64 const > const mass = nodeManager.getExtrinsicData< extrinsicMeshData::MassVector >();

    arrayView1d< real64 > const p_nm1 = nodeManager.getExtrinsicData< extrinsicMeshData::Pressure_nm1 >();
    arrayView1d< real64 > const p_n = nodeManager.getExtrinsicData< extrinsicMeshData::Pressure_n >();
    arrayView1d< real64 > const p_np1 = nodeManager.getExtrinsicData< extrinsicMeshData::Pressure_np1 >();

    if( computeGradient )
    {
      ElementRegionManager & elemManager = mesh.getElemManager();

      arrayView1d< real64 > const p_dt2 = nodeManager.getExtrinsicData< extrinsicMeshData::PressureDoubleDerivative >();

      int const rank = MpiWrapper::commRank( MPI_COMM_GEOSX );
      std::string fileName = GEOSX_FMT( "pressuredt2_{:06}_{:08}_{:04}.dat", m_shotIndex, cycleNumber, rank );
      std::ifstream wf( fileName, std::ios::in | std::ios::binary );
      GEOSX_THROW_IF( !wf,
                      "Could not open file "<< fileName << " for reading",
                      InputError );
      // maybe better with registerTouch()
      p_dt2.move( MemorySpace::host, true );
      wf.read( (char *)&p_dt2[0], p_dt2.size()*sizeof( real64 ) );
      wf.close();
      remove( fileName.c_str() );

      elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                  CellElementSubRegion & elementSubRegion )
      {
        arrayView1d< real64 const > const velocity = elementSubRegion.getExtrinsicData< extrinsicMeshData::MediumVelocity >();
        arrayView1d< real64 > grad = elementSubRegion.getExtrinsicData< extrinsicMeshData::PartialGradient >();
        arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();
        constexpr localIndex numNodesPerElem = 8;

        GEOSX_MARK_SCOPE ( updatePartialGradient );
        forAll< EXEC_POLICY >( elementSubRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const eltIdx )
        {
          for( localIndex i = 0; i < numNodesPerElem; ++i )
          {
            localIndex nodeIdx = elemsToNodes[eltIdx][i];
            grad[eltIdx] += (-2/velocity[eltIdx]) * mass[nodeIdx]/8.0 * (p_dt2[nodeIdx] * p_n[nodeIdx]);
          }
        } );
      } );
    }

    forAll< EXEC_POLICY >( nodeManager.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
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
  GEOSX_MARK_FUNCTION;

  GEOSX_LOG_RANK_0_IF( dt < epsilonLoc, "Warning! Value for dt: " << dt << "s is smaller than local threshold: " << epsilonLoc );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(),
                                  [&] ( string const &,
                                        MeshLevel & mesh,
                                        arrayView1d< string const > const & regionNames )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    arrayView1d< real32 const > const mass = nodeManager.getExtrinsicData< extrinsicMeshData::MassVector >();
    arrayView1d< real32 const > const damping = nodeManager.getExtrinsicData< extrinsicMeshData::DampingVector >();

    arrayView1d< real32 > const p_nm1 = nodeManager.getExtrinsicData< extrinsicMeshData::Pressure_nm1 >();
    arrayView1d< real32 > const p_n = nodeManager.getExtrinsicData< extrinsicMeshData::Pressure_n >();
    arrayView1d< real32 > const p_np1 = nodeManager.getExtrinsicData< extrinsicMeshData::Pressure_np1 >();

    arrayView1d< localIndex const > const freeSurfaceNodeIndicator = nodeManager.getExtrinsicData< extrinsicMeshData::FreeSurfaceNodeIndicator >();
    arrayView1d< real32 > const stiffnessVector = nodeManager.getExtrinsicData< extrinsicMeshData::StiffnessVector >();
    arrayView1d< real32 > const rhs = nodeManager.getExtrinsicData< extrinsicMeshData::ForcingRHS >();

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

    addSourceToRightHandSide( cycleNumber, rhs );

    /// calculate your time integrators
    real64 const dt2 = dt*dt;

    if( !usePML )
    {
      GEOSX_MARK_SCOPE ( updateP );
      forAll< EXEC_POLICY >( nodeManager.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
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
      arrayView2d< real32 > const v_n = nodeManager.getExtrinsicData< extrinsicMeshData::AuxiliaryVar1PML >();
      arrayView2d< real32 > const grad_n = nodeManager.getExtrinsicData< extrinsicMeshData::AuxiliaryVar2PML >();
      arrayView1d< real32 > const divV_n = nodeManager.getExtrinsicData< extrinsicMeshData::AuxiliaryVar3PML >();
      arrayView1d< real32 > const u_n = nodeManager.getExtrinsicData< extrinsicMeshData::AuxiliaryVar4PML >();
      arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition().toViewConst();

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

      GEOSX_MARK_SCOPE ( updatePWithPML );
      forAll< EXEC_POLICY >( nodeManager.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
      {
        if( freeSurfaceNodeIndicator[a] != 1 )
        {
          real32 sigma[3];
          real64 xLocal[ 3 ];

          for( integer i=0; i<3; ++i )
          {
            xLocal[i] = X[a][i];
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
    fieldsToBeSync.addFields( FieldLocation::Node, { extrinsicMeshData::Pressure_np1::key() } );

    if( usePML )
    {
      fieldsToBeSync.addFields( FieldLocation::Node, {
          extrinsicMeshData::AuxiliaryVar1PML::key(),
          extrinsicMeshData::AuxiliaryVar4PML::key() } );
    }

    CommunicationTools & syncFields = CommunicationTools::getInstance();
    syncFields.synchronizeFields( fieldsToBeSync,
                                  mesh,
                                  domain.getNeighbors(),
                                  true );

    // compute the seismic traces since last step.
    arrayView2d< real32 > const pReceivers   = m_pressureNp1AtReceivers.toView();

    computeAllSeismoTraces( time_n, dt, p_np1, p_n, pReceivers );

    /// prepare next step
    forAll< EXEC_POLICY >( nodeManager.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      stiffnessVector[a] = 0.0;
      rhs[a] = 0.0;
    } );

    if( usePML )
    {
      arrayView2d< real32 > const grad_n = nodeManager.getExtrinsicData< extrinsicMeshData::AuxiliaryVar2PML >();
      arrayView1d< real32 > const divV_n = nodeManager.getExtrinsicData< extrinsicMeshData::AuxiliaryVar3PML >();
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
    arrayView1d< real32 const > const p_n = nodeManager.getExtrinsicData< extrinsicMeshData::Pressure_n >();
    arrayView1d< real32 const > const p_np1 = nodeManager.getExtrinsicData< extrinsicMeshData::Pressure_np1 >();
    arrayView2d< real32 > const pReceivers   = m_pressureNp1AtReceivers.toView();
    computeAllSeismoTraces( time_n, 0, p_np1, p_n, pReceivers );
  } );
}

void AcousticWaveEquationSEM::computeAllSeismoTraces( real64 const time_n,
                                                      real64 const dt,
                                                      arrayView1d< real32 const > const var_np1,
                                                      arrayView1d< real32 const > const var_n,
                                                      arrayView2d< real32 > varAtReceivers )
{

  /*
   * In forward case we compute seismo if time_n + dt  is the first time
   * step after the timeSeismo to write.
   *
   *  time_n        timeSeismo    time_n + dt
   *   ---|--------------|-------------|
   *
   * In backward (time_n goes decreasing) case we compute seismo if
   * time_n is the last time step before the timeSeismo to write.
   *
   *  time_n - dt    timeSeismo    time_n
   *   ---|--------------|-------------|
   */
  for( real64 timeSeismo;
       (m_forward)?((timeSeismo = m_dtSeismoTrace*m_indexSeismoTrace) <= (time_n + dt + epsilonLoc) && m_indexSeismoTrace < m_nsamplesSeismoTrace):
       ((timeSeismo = m_dtSeismoTrace*(m_nsamplesSeismoTrace-m_indexSeismoTrace-1)) >= (time_n - dt -  epsilonLoc) && m_indexSeismoTrace < m_nsamplesSeismoTrace);
       m_indexSeismoTrace++ )
  {
    computeSeismoTrace( time_n, (m_forward)?dt:-dt, timeSeismo, (m_forward)?m_indexSeismoTrace:(m_nsamplesSeismoTrace-m_indexSeismoTrace-1), var_np1, var_n, varAtReceivers );
  }
}

REGISTER_CATALOG_ENTRY( SolverBase, AcousticWaveEquationSEM, string const &, dataRepository::Group * const )

} /* namespace geosx */
