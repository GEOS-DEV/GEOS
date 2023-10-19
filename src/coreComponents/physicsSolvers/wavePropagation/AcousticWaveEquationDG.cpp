/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOS Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


/**
 * @file AcousticWaveEquationDG.cpp
 */

#include "AcousticWaveEquationDG.hpp"
//#include "AcousticWaveEquationDGKernel.hpp"


#include "finiteElement/FiniteElementDiscretization.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/ElementType.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "WaveSolverUtils.hpp"

namespace geos
{

using namespace dataRepository;
using namespace fields;

AcousticWaveEquationDG::AcousticWaveEquationDG( const std::string & name,
                                                Group * const parent ):
  WaveSolverBase( name,
                  parent )
{

  registerWrapper( viewKeyStruct::pressureNp1AtReceiversString(), &m_pressureNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Pressure value at each receiver for each timestep" );

  registerWrapper( viewKeyStruct::sourceElemString(), &m_sourceElem ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Element containing the sources" );

  registerWrapper( viewKeyStruct::receiverElemString(), &m_rcvElem ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Element containing the receivers" );

  registerWrapper( viewKeyStruct::receiverRegionString(), &m_receiverRegion ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Region containing the receivers" );

}

AcousticWaveEquationDG::~AcousticWaveEquationDG()
{
  // TODO Auto-generated destructor stub
}

void AcousticWaveEquationDG::initializePreSubGroups()
{
  WaveSolverBase::initializePreSubGroups();
}


void AcousticWaveEquationDG::registerDataOnMesh( Group & meshBodies )
{

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    nodeManager.registerField< wavesolverfields::FreeSurfaceNodeIndicator >( this->getName() );

    FaceManager & faceManager = mesh.getFaceManager();
    faceManager.registerField< wavesolverfields::FreeSurfaceFaceIndicator >( this->getName() );

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
    {
      subRegion.registerField< wavesolverfields::MediumVelocity >( this->getName() );

      subRegion.registerField< wavesolverfields::PressureDG_nm1 >( this->getName() );
      subRegion.registerField< wavesolverfields::PressureDG_n >( this->getName() );
      subRegion.registerField< wavesolverfields::PressureDG_np1 >( this->getName() );
      subRegion.registerField< wavesolverfields::StiffnessVectorDG >( this->getName() );

      finiteElement::FiniteElementBase const & fe = subRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

      finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
      {

        using FE_TYPE = TYPEOFREF( finiteElement );

        constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;

        subRegion.getField< wavesolverfields::PressureDG_nm1 >().resizeDimension< 1 >( numNodesPerElem );
        subRegion.getField< wavesolverfields::PressureDG_n >().resizeDimension< 1 >( numNodesPerElem );
        subRegion.getField< wavesolverfields::PressureDG_np1 >().resizeDimension< 1 >( numNodesPerElem );
        subRegion.getField< wavesolverfields::StiffnessVectorDG >().resizeDimension< 1 >( numNodesPerElem );

      } );

    } );
  } );
}


void AcousticWaveEquationDG::postProcessInput()
{
  WaveSolverBase::postProcessInput();

  localIndex const numSourcesGlobal = m_sourceCoordinates.size( 0 );
  localIndex const numReceiversGlobal = m_receiverCoordinates.size( 0 );
  m_pressureNp1AtReceivers.resize( m_nsamplesSeismoTrace, numReceiversGlobal );
  m_sourceElem.resize( numSourcesGlobal );
  m_rcvElem.resize( numReceiversGlobal );
  m_receiverRegion.resize( numReceiversGlobal );

}

void AcousticWaveEquationDG::precomputeSourceAndReceiverTerm( MeshLevel & mesh, arrayView1d< string const > const & regionNames )
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
  arrayView1d< localIndex > const sourceElem = m_sourceElem.toView();
  sourceNodeIds.setValues< EXEC_POLICY >( -1 );
  sourceConstants.setValues< EXEC_POLICY >( -1 );
  sourceIsAccessible.zero();

  arrayView2d< real64 const > const receiverCoordinates = m_receiverCoordinates.toViewConst();
  arrayView2d< localIndex > const receiverNodeIds = m_receiverNodeIds.toView();
  arrayView2d< real64 > const receiverConstants = m_receiverConstants.toView();
  arrayView1d< localIndex > const receiverIsLocal = m_receiverIsLocal.toView();
  arrayView1d< localIndex > const rcvElem = m_rcvElem.toView();
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


  mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                        CellElementSubRegion & elementSubRegion )
  {
    GEOS_THROW_IF( elementSubRegion.getElementType() != ElementType::Hexahedron,
                   "Invalid type of element, the acoustic solver is designed for hexahedral meshes only (C3D8) ",
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

      constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
      localIndex const numFacesPerElem = elementSubRegion.numFacesPerElement();

      // AcousticWaveEquationDGKernels::
      //   PrecomputeSourceAndReceiverKernel::
      //   launch< EXEC_POLICY, FE_TYPE >
      //   ( elementSubRegion.size(),
      //   numNodesPerElem,
      //   numFacesPerElem,
      //   X,
      //   elemGhostRank,
      //   elemsToNodes,
      //   elemsToFaces,
      //   elemCenter,
      //   faceNormal,
      //   faceCenter,
      //   sourceCoordinates,
      //   sourceIsAccessible,
      //   sourceElem,
      //   sourceNodeIds,
      //   sourceConstants,
      //   receiverCoordinates,
      //   receiverIsLocal,
      //   rcvElem,
      //   receiverNodeIds,
      //   receiverConstants,
      //   sourceValue,
      //   dt,
      //   timeSourceFrequency,
      //   rickerOrder );
    } );
  } );

}

//TODO: Modify to use on discontinuous variable
void AcousticWaveEquationDG::applyFreeSurfaceBC( real64 const time, DomainPartition & domain )
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  FunctionManager const & functionManager = FunctionManager::getInstance();

  FaceManager & faceManager = domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ).getFaceManager();
  NodeManager & nodeManager = domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ).getNodeManager();

  arrayView1d< real32 > const p_np1 = nodeManager.getField< wavesolverfields::Pressure_np1 >();

  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  /// array of indicators: 1 if a face is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceFaceIndicator = faceManager.getField< wavesolverfields::FreeSurfaceFaceIndicator >();

  /// array of indicators: 1 if a node is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceNodeIndicator = nodeManager.getField< wavesolverfields::FreeSurfaceNodeIndicator >();


  freeSurfaceFaceIndicator.zero();
  freeSurfaceNodeIndicator.zero();

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
        }
      }
    }
    else
    {
      GEOS_ERROR( "This option is not supported yet" );
    }
  } );
}

// Here for retrocompatibily
real64 AcousticWaveEquationDG::explicitStepForward( real64 const & time_n,
                                                    real64 const & dt,
                                                    integer cycleNumber,
                                                    DomainPartition & domain,
                                                    bool GEOS_UNUSED_PARAM( computeGradient ) )
{
  real64 dtOut = explicitStepInternal( time_n, dt, cycleNumber, domain );
  return dtOut;
}



real64 AcousticWaveEquationDG::explicitStepBackward( real64 const & time_n,
                                                     real64 const & dt,
                                                     integer cycleNumber,
                                                     DomainPartition & domain,
                                                     bool GEOS_UNUSED_PARAM( computeGradient ) )
{
  GEOS_ERROR( "Backward propagation for the first-order wave propagator not yet implemented" );
  real64 dtOut = explicitStepInternal( time_n, dt, cycleNumber, domain );
  return dtOut;
}


real64 AcousticWaveEquationDG::explicitStepInternal( real64 const & time_n,
                                                     real64 const & dt,
                                                     integer const cycleNumber,
                                                     DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( time_n, dt, cycleNumber );

  arrayView2d< real64 const > const sourceConstants = m_sourceConstants.toView();
  arrayView1d< localIndex const > const sourceIsAccessible = m_sourceIsAccessible.toView();
  arrayView1d< localIndex const > const sourceElem = m_sourceElem.toView();
  arrayView2d< real32 const > const sourceValue = m_sourceValue.toView();

  GEOS_LOG_RANK_0_IF( dt < epsilonLoc, "Warning! Value for dt: " << dt << "s is smaller than local threshold: " << epsilonLoc );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(),
                                  [&] ( string const &,
                                        MeshLevel & mesh,
                                        arrayView1d< string const > const & regionNames )
  {

    //arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition().toViewConst();

    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const regionIndex,
                                                                                          CellElementSubRegion & elementSubRegion )
    {
      arrayView1d< real32 const > const velocity = elementSubRegion.getField< wavesolverfields::MediumVelocity >();
      arrayView2d< real32 > const p_np1 = elementSubRegion.getField< wavesolverfields::PressureDG_np1 >();
      arrayView2d< real32 > const p_n = elementSubRegion.getField< wavesolverfields::PressureDG_n >();
      arrayView2d< real32 > const p_nm1 = elementSubRegion.getField< wavesolverfields::PressureDG_nm1 >();
      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
      finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );
        //TODO: Add the launch calling for flux + volumic computation
      } );

      // compute the seismic traces since last step.
      arrayView2d< real32 > const pReceivers   = m_pressureNp1AtReceivers.toView();
      compute2dVariableAllSeismoTraces( time_n, regionIndex, dt, p_np1, p_n, pReceivers );


    } );

    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addElementFields( {wavesolverfields::PressureDG_nm1::key(), wavesolverfields::PressureDG_n::key(), wavesolverfields::Pressure_np1::key()}, regionNames );

    CommunicationTools & syncFields = CommunicationTools::getInstance();
    syncFields.synchronizeFields( fieldsToBeSync,
                                  mesh,
                                  domain.getNeighbors(),
                                  true );


    // increment m_indexSeismoTrace
    while( (m_dtSeismoTrace*m_indexSeismoTrace) <= (time_n + epsilonLoc) && m_indexSeismoTrace < m_nsamplesSeismoTrace )
    {
      m_indexSeismoTrace++;
    }

  } );
  return dt;
}

void AcousticWaveEquationDG::cleanup( real64 const time_n, integer const, integer const, real64 const, DomainPartition & domain )
{
  // compute the remaining seismic traces, if needed
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    NodeManager & nodeManager = mesh.getNodeManager();
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const regionIndex,
                                                                                          CellElementSubRegion & elementSubRegion )
    {

      arrayView2d< real32 const > const p_np1 = elementSubRegion.getField< wavesolverfields::PressureDG_np1 >();
      arrayView2d< real32 const > const p_n = elementSubRegion.getField< wavesolverfields::PressureDG_n >();

      arrayView2d< real32 > const pReceivers   = m_pressureNp1AtReceivers.toView();
      compute2dVariableAllSeismoTraces( regionIndex, time_n, 0, p_np1, p_n, pReceivers );

    } );
  } );

  // increment m_indexSeismoTrace
  while( (m_dtSeismoTrace*m_indexSeismoTrace) <= (time_n + epsilonLoc) && m_indexSeismoTrace < m_nsamplesSeismoTrace )
  {
    m_indexSeismoTrace++;
  }
}

void AcousticWaveEquationDG::compute2dVariableAllSeismoTraces( localIndex const regionIndex,
                                                               real64 const time_n,
                                                               real64 const dt,
                                                               arrayView2d< real32 const > const var_np1,
                                                               arrayView2d< real32 const > const var_n,
                                                               arrayView2d< real32 > varAtReceivers )
{
  localIndex indexSeismoTrace = m_indexSeismoTrace;
  for( real64 timeSeismo;
       (timeSeismo = m_dtSeismoTrace*indexSeismoTrace) <= (time_n + epsilonLoc) && indexSeismoTrace < m_nsamplesSeismoTrace;
       indexSeismoTrace++ )
  {
    WaveSolverUtils::compute2dVariableSeismoTrace( time_n, dt, regionIndex, m_receiverRegion, timeSeismo, indexSeismoTrace, m_rcvElem, m_receiverConstants, m_receiverIsLocal, m_nsamplesSeismoTrace,
                                                   m_outputSeismoTrace,
                                                   var_np1, var_n, varAtReceivers );
  }
}


void AcousticWaveEquationDG::initializePML()
{
  GEOS_ERROR( "PML for the first order acoustic wave propagator not yet implemented" );
}

void AcousticWaveEquationDG::applyPML( real64 const, DomainPartition & )
{
  GEOS_ERROR( "PML for the first order acoustic wave propagator not yet implemented" );
}

REGISTER_CATALOG_ENTRY( SolverBase, AcousticWaveEquationDG, string const &, dataRepository::Group * const )

} /* namespace geos */
