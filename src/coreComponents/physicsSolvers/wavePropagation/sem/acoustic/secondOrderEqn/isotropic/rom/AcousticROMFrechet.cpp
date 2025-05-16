
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
 * @file AcousticROMFrechet.cpp
 */

#include "AcousticROMFrechet.hpp"
#include "AcousticROMFrechetKernel.hpp"

#include "finiteElement/FiniteElementDiscretization.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "fieldSpecification/PerfectlyMatchedLayer.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/ElementType.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/wavePropagation/shared/WaveSolverUtils.hpp"
#include "physicsSolvers/wavePropagation/sem/acoustic/shared/AcousticTimeSchemeSEMKernel.hpp"
#include "physicsSolvers/wavePropagation/sem/acoustic/shared/AcousticMatricesSEMKernel.hpp"
#include "events/EventManager.hpp"
#include "physicsSolvers/wavePropagation/shared/PrecomputeSourcesAndReceiversKernel.hpp"
#include "denseLinearAlgebra/interfaces/blaslapack/BlasLapackLA.hpp"

namespace geos
{

using namespace dataRepository;
using namespace fields;

AcousticROMFrechet::AcousticROMFrechet( const std::string & name,
                                                  Group * const parent ):
  WaveSolverBase( name,
                  parent )
{

  registerWrapper( viewKeyStruct::pressureNp1AtReceiversString(), &m_pressureNp1AtReceivers ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Pressure value at each receiver for each timestep" );

  registerWrapper( viewKeyStruct::orderFrechetString(), &m_orderFrechet ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Frechet derivative order computation" );

  registerWrapper( viewKeyStruct::orderGSString(), &m_orderGS ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( -1 ).
    setDescription( "Gram-Schmidt order computation" );

  registerWrapper( viewKeyStruct::epsilonGSString(), &m_epsilonGS ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1e-2 ).
    setDescription( "Precision threshold in gramSchmidtStiffness process" );

  registerWrapper( viewKeyStruct::sizePOD_fString(), &m_sizePOD_f ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Number of basis function for each order" );

  registerWrapper( viewKeyStruct::sizePODString(), &m_sizePOD ).
    setInputFlag( InputFlags::FALSE ).
    setApplyDefaultValue( 0 ).
    setDescription( "Total number of basis function" );

  registerWrapper( viewKeyStruct::selectionOrderString(), &m_selectionOrder ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Order of selection of each basis functions" );

  registerWrapper( viewKeyStruct::cycleOrderString(), &m_cycleOrder ).
    setInputFlag( InputFlags::FALSE ).
    setSizedFromParent( 0 ).
    setDescription( "Cycle at which the basis functions is selected" );

  registerWrapper( viewKeyStruct::solverROMString(), &m_solverROM ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "Flag to wheter use the ROM solver (1) or not (0)" );

  registerWrapper( viewKeyStruct::alphaString(), &m_alpha ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0 ).
    setDescription( "scalar for perturbation" );
}

AcousticROMFrechet::~AcousticROMFrechet()
{
  // TODO Auto-generated destructor stub
}

void AcousticROMFrechet::initializePreSubGroups()
{
  WaveSolverBase::initializePreSubGroups();
}


void AcousticROMFrechet::registerDataOnMesh( Group & meshBodies )
{
  WaveSolverBase::registerDataOnMesh( meshBodies );
  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    arrayView1d< string const > const & )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    nodeManager.registerField< acousticfields::Pressure_nm1,
			       acousticfields::Pressure_n,
			       acousticfields::Pressure_np1,
			       acousticfields::PressureForward,
			       acousticfields::PressureFrechet_nm1,
			       acousticfields::PressureFrechet_n,
			       acousticfields::PressureFrechet_np1,
                               acousticfields::ForcingRHS,
			       acousticfields::ForcingRHS_fp1,
			       acousticfields::AcousticMassVector,
			       acousticfields::AcousticMassPerturbationVector,
                               acousticfields::DampingVector,
			       acousticfields::DampingPerturbationVector,
                               acousticfields::StiffnessVector,
                               acousticfields::AcousticFreeSurfaceNodeIndicator >( getName() );


    FaceManager & faceManager = mesh.getFaceManager();
    faceManager.registerField< acousticfields::AcousticFreeSurfaceFaceIndicator >( getName() );

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
    {
      subRegion.registerField< acousticfields::AcousticVelocity >( getName() );
      subRegion.registerField< acousticfields::AcousticDensity >( getName() );
      subRegion.registerField< acousticfields::PartialGradient >( getName() );
      subRegion.registerField< acousticfields::PartialGradient2 >( getName() );
      subRegion.registerField< acousticfields::Perturbation >( getName() );
    } );

  } );
}


void AcousticROMFrechet::postInputInitialization()
{
  WaveSolverBase::postInputInitialization();

  m_pressureNp1AtReceivers.resize( m_nsamplesSeismoTrace, m_receiverCoordinates.size( 0 ) + 1 );
  if( m_solverROM == 0)
  {
    m_sizePOD_f.resize( m_orderFrechet + 1 );
    m_selectionOrder.resize( 2, 2000 );
    m_cycleOrder.resize( m_orderFrechet + 1, 2000 );
  }
}

void AcousticROMFrechet::precomputeSourceAndReceiverTerm( MeshLevel & baseMesh, MeshLevel & mesh,
                                                               arrayView1d< string const > const & regionNames )
{
  GEOS_MARK_FUNCTION;
  arrayView1d< globalIndex const > const nodeLocalToGlobalMap = baseMesh.getNodeManager().localToGlobalMap().toViewConst();
  ArrayOfArraysView< localIndex const > const nodesToElements = baseMesh.getNodeManager().elementList().toViewConst();
  ArrayOfArraysView< localIndex const > const facesToNodes = baseMesh.getFaceManager().nodeList().toViewConst();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const nodeCoords = baseMesh.getNodeManager().referencePosition();

  arrayView2d< real64 const > const sourceCoordinates = m_sourceCoordinates.toViewConst();
  arrayView2d< localIndex > const sourceNodeIds = m_sourceNodeIds.toView();
  arrayView2d< real64 > sourceConstants;

  arrayView2d< real64 const > const receiverCoordinates = m_receiverCoordinates.toViewConst();
  arrayView2d< localIndex > const receiverNodeIds = m_receiverNodeIds.toView();
  arrayView2d< real64 > receiverConstants;
  
  if(m_solverROM)
  {
    m_sourceConstantsPOD.resize( sourceCoordinates.size( 0 ), m_sizePOD );
    sourceConstants = m_sourceConstantsPOD.toView();

    m_receiverConstantsPOD.resize( receiverCoordinates.size( 0 ), m_sizePOD );
    receiverConstants = m_receiverConstantsPOD.toView();
  }
  else
  {
    sourceConstants = m_sourceConstants.toView();

    receiverConstants = m_receiverConstants.toView();
  }
  arrayView1d< localIndex > const sourceIsAccessible = m_sourceIsAccessible.toView();
  sourceNodeIds.setValues< EXEC_POLICY >( -1 );
  sourceConstants.setValues< EXEC_POLICY >( 0 );
  sourceIsAccessible.zero();

  arrayView1d< localIndex > const receiverIsLocal = m_receiverIsLocal.toView();
  receiverNodeIds.setValues< EXEC_POLICY >( -1 );
  receiverConstants.setValues< EXEC_POLICY >( 0 );
  receiverIsLocal.zero();

  //Correct size for sourceValue
  EventManager const & event = getGroupByPath< EventManager >( "/Problem/Events" );
  real64 const & maxTime = event.getReference< real64 >( EventManager::viewKeyStruct::maxTimeString() );
  real64 const & minTime = event.getReference< real64 >( EventManager::viewKeyStruct::minTimeString() );
  real64 dt = 0;
  for( localIndex numSubEvent = 0; numSubEvent < event.numSubGroups(); ++numSubEvent )
  {
    EventBase const * subEvent = static_cast< EventBase const * >( event.getSubGroups()[numSubEvent] );
    if( subEvent->getEventName() == "/Solvers/" + this->getName() )
    {
      dt = subEvent->getReference< real64 >( EventBase::viewKeyStruct::forceDtString() );
    }
  }

  real64 dtCompute;

  localIndex nSubSteps = (int) ceil( dt/m_timeStep );
  dtCompute = dt/nSubSteps;

  localIndex const nsamples = int( (maxTime - minTime) / dtCompute) + 1;

  localIndex const numSourcesGlobal = m_sourceCoordinates.size( 0 );
  m_sourceValue.resize( nsamples, numSourcesGlobal );

  arrayView2d< real64 > const sourceValue = m_sourceValue.toView();

  mesh.getElemManager().forElementSubRegionsComplete< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                                localIndex const er,
                                                                                                localIndex const esr,
                                                                                                ElementRegionBase &,
                                                                                                CellElementSubRegion & elementSubRegion )
  {
    GEOS_THROW_IF( elementSubRegion.getElementType() != ElementType::Hexahedron,
                   getDataContext() << ": Invalid type of element, the acoustic solver is designed for hexahedral meshes only (C3D8), using the SEM formulation",
                   InputError );

    arrayView2d< localIndex const > const elemsToFaces = elementSubRegion.faceList();
    arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();
    arrayView2d< localIndex const, cells::NODE_MAP_USD > const & baseElemsToNodes = baseMesh.getElemManager().getRegion( er ).getSubRegion< CellElementSubRegion >( esr ).nodeList();
    arrayView2d< real64 const > const elemCenter = elementSubRegion.getElementCenter();
    arrayView1d< integer const > const elemGhostRank = elementSubRegion.ghostRank();
    arrayView1d< globalIndex const > const elemLocalToGlobalMap = elementSubRegion.localToGlobalMap().toViewConst();

    if(m_solverROM)
    {
      
      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
      finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
      {
	using FE_TYPE = TYPEOFREF( finiteElement );
	acousticROMFrechetKernels::
          PrecomputeSourceAndReceiverWithoutConstantsKernel::
          launch< EXEC_POLICY, FE_TYPE >
          ( elementSubRegion.size(),
            facesToNodes,
            nodeCoords,
            nodeLocalToGlobalMap,
            elemLocalToGlobalMap,
            nodesToElements,
            baseElemsToNodes,
            elemGhostRank,
            elemsToNodes,
            elemsToFaces,
            elemCenter,
            sourceCoordinates,
            sourceIsAccessible,
            receiverCoordinates,
	    receiverIsLocal,
            sourceValue,
            dtCompute,
            m_timeSourceFrequency,
            m_timeSourceDelay,
            m_rickerOrder );
	
      } );

      
      int const rank = MpiWrapper::commRank( MPI_COMM_GEOS );
      std::vector<std::string> matrixNames = {"sourceConstants", "receiverConstants"};

      for( localIndex i=0; i<2; ++i )
      {
	std::string fileName = GEOS_FMT( "phi/shot_{:05}/finalBases/rank_{:05}/{:}.dat", m_shotIndex, rank, matrixNames[i]);
	std::ifstream wf( fileName, std::ios::in | std::ios::binary );
	GEOS_THROW_IF( !wf,
		       ": Could not open file "<< fileName << " for reading",
		       InputError );
	switch( i )
        {
	case 0:
        {
	  wf.read( (char *)&sourceConstants[0][0], sourceConstants.size(0)*m_sizePOD*sizeof( real64 ) );
	}
	case 1:
        {
	  wf.read( (char *)&receiverConstants[0][0], receiverConstants.size(0)*m_sizePOD*sizeof( real64 ) );
	}
	}
	wf.close( );
      }
      
      MpiWrapper::allReduce( sourceConstants.data(),
			     sourceConstants.data(),
			     sourceConstants.size( 0 )*sourceConstants.size( 1 ),
			     MpiWrapper::getMpiOp( MpiWrapper::Reduction::Sum ),
			     MPI_COMM_GEOS );

      sourceValue.move( LvArray::MemorySpace::host, true );
      MpiWrapper::allReduce( sourceValue.data(),
			     sourceValue.data(),
			     sourceValue.size( 0 )*sourceValue.size( 1 ),
			     MpiWrapper::getMpiOp( MpiWrapper::Reduction::Sum ),
			     MPI_COMM_GEOS );

      receiverIsLocal.move( LvArray::MemorySpace::host, false );
    }
    else
    {
       finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
      finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
      {
	using FE_TYPE = TYPEOFREF( finiteElement );
	GEOS_MARK_SCOPE( acousticWaveEquationSEMKernels::PrecomputeSourceAndReceiverKernel );
	PreComputeSourcesAndReceivers::
	  Compute1DSourceAndReceiverConstants
	  < EXEC_POLICY, FE_TYPE >
	  ( elementSubRegion.size(),
	    facesToNodes,
	    nodeCoords,
	    nodeLocalToGlobalMap,
	    elemLocalToGlobalMap,
	    nodesToElements,
	    baseElemsToNodes,
	    elemGhostRank,
	    elemsToNodes,
	    elemsToFaces,
	    elemCenter,
	    sourceCoordinates,
	    sourceIsAccessible,
	    sourceNodeIds,
	    sourceConstants,
	    receiverCoordinates,
	    receiverIsLocal,
	    receiverNodeIds,
	    receiverConstants,
	    sourceValue,
	    dtCompute,
	    m_timeSourceFrequency,
	    m_timeSourceDelay,
	    m_rickerOrder );
      } );
    }
    elementSubRegion.faceList().freeOnDevice();
    baseMesh.getElemManager().getRegion( er ).getSubRegion< CellElementSubRegion >( esr ).nodeList().freeOnDevice();
    elementSubRegion.getElementCenter().freeOnDevice();
    elementSubRegion.ghostRank().freeOnDevice();
    elementSubRegion.localToGlobalMap().freeOnDevice();
  } );
  baseMesh.getNodeManager().localToGlobalMap().freeOnDevice();
  baseMesh.getNodeManager().elementList().toView().freeOnDevice();
  baseMesh.getFaceManager().nodeList().toView().freeOnDevice();
  baseMesh.getNodeManager().referencePosition().freeOnDevice();
  facesToNodes.freeOnDevice();
  nodesToElements.freeOnDevice();
}


void AcousticROMFrechet::addSourceToRightHandSide( integer const & cycleNumber, arrayView1d< real32 > const rhs )
{
  arrayView2d< localIndex const > const sourceNodeIds = m_sourceNodeIds.toViewConst();
  arrayView2d< real64 > sourceConstants;
  if(m_solverROM)
  {
    sourceConstants = m_sourceConstantsPOD.toView();
  }
  else
  {
    sourceConstants = m_sourceConstants.toView();
  }
  arrayView1d< localIndex const > const sourceIsAccessible = m_sourceIsAccessible.toViewConst();
  arrayView2d< real64 const > const sourceValue   = m_sourceValue.toViewConst();

  GEOS_THROW_IF( cycleNumber > sourceValue.size( 0 ),
                 getDataContext() << ": Too many steps compared to array size",
                 std::runtime_error );
  if(m_solverROM)
  {
    for( localIndex isrc = 0; isrc < sourceConstants.size( 0 ); ++isrc )
    {
      for( localIndex inode = 0; inode < sourceConstants.size( 1 ); ++inode )
      {
        rhs[inode] = sourceConstants[isrc][inode] * sourceValue[cycleNumber][isrc];
      }
    }
  }
  else
  {
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
}

void AcousticROMFrechet::initializePostInitialConditionsPreSubGroups()
{
  GEOS_MARK_FUNCTION;
  {
    GEOS_MARK_SCOPE( WaveSolverBase::initializePostInitialConditionsPreSubGroups );
    WaveSolverBase::initializePostInitialConditionsPreSubGroups();
  }

  DomainPartition & domain = getGroupByPath< DomainPartition >( "/Problem/domain" );

  if( m_solverROM == 0 )
  {
    applyFreeSurfaceBC( 0.0, domain );
  }

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshBodyName,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    MeshLevel & baseMesh = domain.getMeshBodies().getGroup< MeshBody >( meshBodyName ).getBaseDiscretization();
    
    if( m_solverROM )
    {
      EventManager const & event = getGroupByPath< EventManager >( "/Problem/Events" );
      real64 dt = 0;
      for( localIndex numSubEvent = 0; numSubEvent < event.numSubGroups(); ++numSubEvent )
      {
	EventBase const * subEvent = static_cast< EventBase const * >( event.getSubGroups()[numSubEvent] );
	if( subEvent->getEventName() == "/Solvers/" + this->getName() )
	{
	  dt = subEvent->getReference< real64 >( EventBase::viewKeyStruct::forceDtString() );
	}
      }

      int const rank = MpiWrapper::commRank( MPI_COMM_GEOS );
      int count = 0;
      localIndex ifile = 1;
      while(ifile > 0)
      {
	std::string filename = GEOS_FMT( "phi/shot_{:05}/finalBases/rank_{:05}/vector_{:03}.dat", m_shotIndex, rank, ifile);
	std::ifstream wf(filename.c_str());
	if( wf.good() )
	{
	  ++count;
	  ++ifile;
	}
	else
	{
	  ifile = -1;
	}
      }
      m_sizePOD = count;

      m_massPOD.resize(m_sizePOD, m_sizePOD);
      m_dampingPOD.resize(m_sizePOD, m_sizePOD);
      m_massPerturbationPOD.resize(m_sizePOD, m_sizePOD);
      m_dampingPerturbationPOD.resize(m_sizePOD, m_sizePOD);
      
      m_rhsPOD.resize(m_sizePOD);
      m_a_np1.resize(m_sizePOD);
      m_a_n.resize(m_sizePOD);
      m_a_nm1.resize(m_sizePOD);
      
      m_massPOD.zero();
      m_dampingPOD.zero();
      m_massPerturbationPOD.zero();
      m_dampingPerturbationPOD.zero();
      
      m_rhsPOD.zero();
      m_a_np1.zero();
      m_a_n.zero();
      m_a_nm1.zero();
      
      arrayView2d < real32 > const massPOD = m_massPOD.toView();
      arrayView2d < real32 > const dampingPOD = m_dampingPOD.toView();
      arrayView2d < real32 > const massPerturbationPOD = m_massPerturbationPOD.toView();
      arrayView2d < real32 > const dampingPerturbationPOD = m_dampingPerturbationPOD.toView();
      
      std::vector<std::string> matrixNames = {"mass", "massPerturbation", "damping", "dampingPerturbation"};
      
      for( localIndex i=0; i<4; ++i )
      {
	std::string fileName = GEOS_FMT( "phi/shot_{:05}/finalBases/rank_{:05}/{:}.dat", m_shotIndex, rank, matrixNames[i]);
	std::ifstream wf( fileName, std::ios::in | std::ios::binary );
	GEOS_THROW_IF( !wf,
		       ": Could not open file "<< fileName << " for reading",
		       InputError );
	switch( i )
	{
	case 0:
	{
	  //massPOD.move( LvArray::MemorySpace::host, true);
	  wf.read( (char *)&massPOD[0][0], m_sizePOD*m_sizePOD*sizeof( real32 ) );
	  MpiWrapper::allReduce( massPOD.data(),
				 massPOD.data(),
				 massPOD.size( 0 )*massPOD.size( 1 ),
				 MpiWrapper::getMpiOp( MpiWrapper::Reduction::Sum ),
				 MPI_COMM_GEOS );
	}
	case 1:
	{
	  //massPerturbationPOD.move( LvArray::MemorySpace::host, true);
	  wf.read( (char *)&massPerturbationPOD[0][0], m_sizePOD*m_sizePOD*sizeof( real32 ) );
	  MpiWrapper::allReduce( massPerturbationPOD.data(),
				 massPerturbationPOD.data(),
				 massPerturbationPOD.size( 0 )*massPerturbationPOD.size( 1 ),
				 MpiWrapper::getMpiOp( MpiWrapper::Reduction::Sum ),
				 MPI_COMM_GEOS );
	}
	case 2:
	{
	  //dampingPOD.move( LvArray::MemorySpace::host, true);
	  wf.read( (char *)&dampingPOD[0][0], m_sizePOD*m_sizePOD*sizeof( real32 ) );
	  MpiWrapper::allReduce( dampingPOD.data(),
				 dampingPOD.data(),
				 dampingPOD.size( 0 )*dampingPOD.size( 1 ),
				 MpiWrapper::getMpiOp( MpiWrapper::Reduction::Sum ),
				 MPI_COMM_GEOS );
	}
	case 3:
	{
	  //dampingPerturbationPOD.move( LvArray::MemorySpace::host, true);
	  wf.read( (char *)&dampingPerturbationPOD[0][0], m_sizePOD*m_sizePOD*sizeof( real32 ) );
	  MpiWrapper::allReduce( dampingPerturbationPOD.data(),
				 dampingPerturbationPOD.data(),
				 dampingPerturbationPOD.size( 0 )*dampingPerturbationPOD.size( 1 ),
				 MpiWrapper::getMpiOp( MpiWrapper::Reduction::Sum ),
				 MPI_COMM_GEOS );
	}
	}
	wf.close( );
      }
	
      array2d< real64 > A(m_sizePOD, m_sizePOD);
      for(localIndex m=0; m<m_sizePOD; ++m)
      {
	for(localIndex n=0; n<=m; ++n)
	{
	  A[m][n] = massPOD[m][n] + m_alpha * massPerturbationPOD[m][n] + dt/2 * (dampingPOD[m][n] + m_alpha * dampingPerturbationPOD[m][n]);
	  if( m!=n )
	  {
	    A[n][m] = A[m][n];
	  }
	}
      }
      m_OpPOD.resize(m_sizePOD,m_sizePOD);
      m_OpPOD.zero();

      geos::BlasLapackLA::matrixInverse(A, m_OpPOD);
    }
    else
    {
      NodeManager & nodeManager = mesh.getNodeManager();
      FaceManager & faceManager = mesh.getFaceManager();
      ElementRegionManager & elemManager = mesh.getElemManager();

      m_sizePOD = 0;
      m_sizePOD_f.zero();
      if(m_orderFrechet > 0)
      {
	nodeManager.getField< acousticfields::PressureFrechet_np1 >().resizeDimension< 1 >(m_orderFrechet);
	nodeManager.getField< acousticfields::PressureFrechet_n >().resizeDimension< 1 >(m_orderFrechet);
	nodeManager.getField< acousticfields::PressureFrechet_nm1 >().resizeDimension< 1 >(m_orderFrechet);
	
	nodeManager.getField< acousticfields::PressureFrechet_np1 >().zero();
	nodeManager.getField< acousticfields::PressureFrechet_n >().zero();
	nodeManager.getField< acousticfields::PressureFrechet_nm1 >().zero();
      }
      /// get the array of indicators: 1 if the face is on the boundary; 0 otherwise
      arrayView1d< integer const > const & facesDomainBoundaryIndicator = faceManager.getDomainBoundaryIndicator();
      arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords = nodeManager.getField< fields::referencePosition32 >().toViewConst();
      
      /// get face to node map
      ArrayOfArraysView< localIndex const > const facesToNodes = faceManager.nodeList().toViewConst();
      
      // mass matrix to be computed in this function
      arrayView1d< real32 > const mass = nodeManager.getField< acousticfields::AcousticMassVector >();
      {
	GEOS_MARK_SCOPE( mass_zero );
	mass.zero();
      }
      /// damping matrix to be computed for each dof in the boundary of the mesh
      arrayView1d< real32 > const damping = nodeManager.getField< acousticfields::DampingVector >();
      {
	GEOS_MARK_SCOPE( damping_zero );
	damping.zero();
      }
      
      arrayView1d< real32 > massPerturbation;
      arrayView1d< real32 > dampingPerturbation;
      if(m_orderFrechet > 0)
      {
	// mass matrix with gradient to be computed in this function
	massPerturbation = nodeManager.getField< acousticfields::AcousticMassPerturbationVector >();
	{
	  GEOS_MARK_SCOPE( mass_zero );
	  massPerturbation.zero();
	}
	// damping matrix with gradient to be computed for each dof in the boundary of the mesh
	dampingPerturbation = nodeManager.getField< acousticfields::DampingPerturbationVector >();
	{
	  GEOS_MARK_SCOPE( damping_zero );
	  dampingPerturbation.zero();
	}
      }
   
      /// get array of indicators: 1 if face is on the free surface; 0 otherwise
      arrayView1d< localIndex const > const freeSurfaceFaceIndicator = faceManager.getField< acousticfields::AcousticFreeSurfaceFaceIndicator >();
      
      elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
										  CellElementSubRegion & elementSubRegion )
      {
	
	finiteElement::FiniteElementBase const &
	  fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
	
	arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes = elementSubRegion.nodeList();
	arrayView2d< localIndex const > const elemsToFaces = elementSubRegion.faceList();
	
	computeTargetNodeSet( elemsToNodes, elementSubRegion.size(), fe.getNumQuadraturePoints() );
	
	arrayView1d< real32 const > const velocity = elementSubRegion.getField< acousticfields::AcousticVelocity >();
	arrayView1d< real32 const > const density = elementSubRegion.getField< acousticfields::AcousticDensity >();
	
	/// Partial gradient if gradient has to be computed
	arrayView1d< real32 > grad = elementSubRegion.getField< acousticfields::PartialGradient >();
	grad.zero();
	
	finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
	{
	  using FE_TYPE = TYPEOFREF( finiteElement );
	  
	  acousticROMFrechetKernels::MassMatrixKernel< FE_TYPE > kernelM( finiteElement );
	  kernelM.template launch< EXEC_POLICY, ATOMIC_POLICY >( elementSubRegion.size(),
								 nodeCoords,
								 elemsToNodes,
								 velocity,
								 density,
								 mass );
	  
	  GEOS_MARK_SCOPE( DampingMatrixKernel );
	  acousticROMFrechetKernels::DampingMatrixKernel< FE_TYPE > kernelD( finiteElement );
	  kernelD.template launch< EXEC_POLICY, ATOMIC_POLICY >( elementSubRegion.size(),
								 nodeCoords,
								 elemsToFaces,
								 facesToNodes,
								 facesDomainBoundaryIndicator,
								 freeSurfaceFaceIndicator,
								 velocity,
								 density,
								 damping );
	  
	} );

	if( m_orderFrechet>0 )
        {
          arrayView1d< real32 const > const perturbation = elementSubRegion.getField< acousticfields::Perturbation >();

          finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
          {
          using FE_TYPE = TYPEOFREF( finiteElement );

          acousticROMFrechetKernels::MassPerturbationMatrixKernel< FE_TYPE > kernelM( finiteElement );
          kernelM.template launch< EXEC_POLICY, ATOMIC_POLICY >( elementSubRegion.size(),
                                                                 nodeCoords,
                                                                 elemsToNodes,
                                                                 perturbation,
                                                                 massPerturbation );

          GEOS_MARK_SCOPE( DampingMatrixKernel );
          acousticROMFrechetKernels::DampingPerturbationMatrixKernel< FE_TYPE > kernelD( finiteElement );
          kernelD.template launch< EXEC_POLICY, ATOMIC_POLICY >( elementSubRegion.size(),
                                                                 nodeCoords,
                                                                 elemsToFaces,
                                                                 facesToNodes,
                                                                 facesDomainBoundaryIndicator,
                                                                 freeSurfaceFaceIndicator,
                                                                 velocity,
                                                                 perturbation,
                                                                 dampingPerturbation );
          } );
        }
      } );
    }
    precomputeSourceAndReceiverTerm(baseMesh, mesh, regionNames );
  } );

  WaveSolverUtils::initTrace( "seismoTraceReceiver", getName(), m_outputSeismoTrace, m_receiverConstants.size( 0 ), m_receiverIsLocal );
}


void AcousticROMFrechet::applyFreeSurfaceBC( real64 time, DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  FunctionManager const & functionManager = FunctionManager::getInstance();

  FaceManager & faceManager = domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ).getFaceManager();
  NodeManager & nodeManager = domain.getMeshBody( 0 ).getMeshLevel( m_discretizationName ).getNodeManager();

  arrayView1d< real32 > const p_nm1 = nodeManager.getField< acousticfields::Pressure_nm1 >();
  arrayView1d< real32 > const p_n = nodeManager.getField< acousticfields::Pressure_n >();
  arrayView1d< real32 > const p_np1 = nodeManager.getField< acousticfields::Pressure_np1 >();

  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  /// array of indicators: 1 if a face is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceFaceIndicator = faceManager.getField< acousticfields::AcousticFreeSurfaceFaceIndicator >();

  /// array of indicators: 1 if a node is on on free surface; 0 otherwise
  arrayView1d< localIndex > const freeSurfaceNodeIndicator = nodeManager.getField< acousticfields::AcousticFreeSurfaceNodeIndicator >();

  // freeSurfaceFaceIndicator.zero();
  // freeSurfaceNodeIndicator.zero();

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


real64 AcousticROMFrechet::explicitStepForward( real64 const & time_n,
                                                real64 const & dt,
						integer cycleNumber,
						DomainPartition & domain,
						bool computeGradient )
{
  real64 dtCompute = explicitStepInternal( time_n, dt, domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(),
                                  [&] ( string const &,
                                        MeshLevel & mesh,
                                        arrayView1d< string const > const & GEOS_UNUSED_PARAM ( regionNames ) )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    arrayView1d< real32 > const p_n = nodeManager.getField< acousticfields::Pressure_n >();

    if( computeGradient && cycleNumber >= 0 )
    {

      if( m_enableLifo )
      {
        if( !m_lifo )
        {
          int const rank = MpiWrapper::commRank( MPI_COMM_GEOS );
          std::string lifoPrefix = GEOS_FMT( "lifo/rank_{:05}/p_forward_shot{:06}", rank, m_shotIndex );
          m_lifo = std::make_unique< LifoStorage< real32, localIndex > >( lifoPrefix, p_n, m_lifoOnDevice, m_lifoOnHost, m_lifoSize );
        }

        m_lifo->pushWait();
      }

      if( m_enableLifo )
      {
        // Need to tell LvArray data is on GPU to avoir HtoD copy
        p_n.move( LvArray::MemorySpace::cuda, false );
        m_lifo->pushAsync( p_n );
      }
      else
      {
        GEOS_MARK_SCOPE ( DirectWrite );
        p_n.move( LvArray::MemorySpace::host, false );
        int const rank = MpiWrapper::commRank( MPI_COMM_GEOS );
        std::string fileName = GEOS_FMT( "lifo/rank_{:05}/pressure_forward_{:06}_{:08}.dat", rank, m_shotIndex, cycleNumber );
        int lastDirSeparator = fileName.find_last_of( "/\\" );
        std::string dirName = fileName.substr( 0, lastDirSeparator );
        if( string::npos != (size_t)lastDirSeparator && !directoryExists( dirName ))
        {
          makeDirsForPath( dirName );
        }

        std::ofstream wf( fileName, std::ios::out | std::ios::binary );
        GEOS_THROW_IF( !wf,
                       getDataContext() << ": Could not open file "<< fileName << " for writing",
                       InputError );
        wf.write( (char *)&p_n[0], p_n.size()*sizeof( real32 ) );
        wf.close( );
        GEOS_THROW_IF( !wf.good(),
                       getDataContext() << ": An error occured while writing "<< fileName,
                       InputError );
      }

    }

    prepareNextTimestep( mesh );
  } );

  return dtCompute;
}


real64 AcousticROMFrechet::explicitStepBackward( real64 const & time_n,
                                                      real64 const & dt,
                                                      integer cycleNumber,
                                                      DomainPartition & domain,
                                                      bool computeGradient )
{
  real64 dtCompute = explicitStepInternal( time_n, dt, domain );
  forDiscretizationOnMeshTargets( domain.getMeshBodies(),
                                  [&] ( string const &,
                                        MeshLevel & mesh,
                                        arrayView1d< string const > const & regionNames )
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    arrayView1d< real32 > const p_nm1 = nodeManager.getField< acousticfields::Pressure_nm1 >();
    arrayView1d< real32 > const p_n = nodeManager.getField< acousticfields::Pressure_n >();
    arrayView1d< real32 > const p_np1 = nodeManager.getField< acousticfields::Pressure_np1 >();

    //// Compute q_dt2 and store it in p_nm1
    SortedArrayView< localIndex const > const solverTargetNodesSet = m_solverTargetNodesSet.toViewConst();
    forAll< EXEC_POLICY >( solverTargetNodesSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const n )
    {
      localIndex const a = solverTargetNodesSet[n];
      p_nm1[a] = (p_np1[a] - 2*p_n[a] + p_nm1[a]) / pow( dt, 2 );
    } );

    EventManager const & event = getGroupByPath< EventManager >( "/Problem/Events" );
    real64 const & maxTime = event.getReference< real64 >( EventManager::viewKeyStruct::maxTimeString() );
    int const maxCycle = int(round( maxTime / dt ));

    if( computeGradient && cycleNumber < maxCycle )
    {
      ElementRegionManager & elemManager = mesh.getElemManager();

      arrayView1d< real32 > const p_forward = nodeManager.getField< acousticfields::PressureForward >();

      if( m_enableLifo )
      {
        m_lifo->pop( p_forward );
        if( m_lifo->empty() )
          delete m_lifo.release();
      }
      else
      {
        GEOS_MARK_SCOPE ( DirectRead );

        int const rank = MpiWrapper::commRank( MPI_COMM_GEOS );
        std::string fileName = GEOS_FMT( "lifo/rank_{:05}/pressure_forward_{:06}_{:08}.dat", rank, m_shotIndex, cycleNumber );
        std::ifstream wf( fileName, std::ios::in | std::ios::binary );
        GEOS_THROW_IF( !wf,
                       getDataContext() << ": Could not open file "<< fileName << " for reading",
                       InputError );

        p_forward.move( LvArray::MemorySpace::host, true );
        wf.read( (char *)&p_forward[0], p_forward.size()*sizeof( real32 ) );
        wf.close( );
        remove( fileName.c_str() );
      }
      elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                  CellElementSubRegion & elementSubRegion )
      {
        arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords = nodeManager.getField< fields::referencePosition32 >().toViewConst();
        arrayView1d< real32 > grad = elementSubRegion.getField< acousticfields::PartialGradient >();
        arrayView1d< real32 > grad2 = elementSubRegion.getField< acousticfields::PartialGradient2 >();
        arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();
        arrayView1d< integer const > const elemGhostRank = elementSubRegion.ghostRank();
        GEOS_MARK_SCOPE ( updatePartialGradient );

        //COMPUTE GRADIENTS with respect to K=1/rho*c2 (grad) and b=1/rho (grad2)
        finiteElement::FiniteElementBase const &
        fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

        finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
        {
          using FE_TYPE = TYPEOFREF( finiteElement );

          AcousticMatricesSEM::GradientKappaBuoyancy< FE_TYPE > kernelG( finiteElement );
          kernelG.template computeGradient< EXEC_POLICY, ATOMIC_POLICY >( elementSubRegion.size(),
                                                                          nodeCoords,
                                                                          elemsToNodes,
                                                                          elemGhostRank,
                                                                          p_nm1,
                                                                          p_n,
                                                                          p_forward,
                                                                          grad,
                                                                          grad2 );


        } );

        // // Change of variables to return grad with respect to c and rho
        // arrayView1d< real32 const > const velocity = elementSubRegion.getField< acousticfields::AcousticVelocity >();
        // arrayView1d< real32 const > const density = elementSubRegion.getField< acousticfields::AcousticDensity >();
        // forAll< EXEC_POLICY >( elementSubRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const eltIdx )
        // {
        //   if( elemGhostRank[eltIdx]<0 )
        //   {
        //     grad2[eltIdx] = -1/(pow(density[eltIdx]*velocity[eltIdx],2)) * grad[eltIdx] - 1/pow(density[eltIdx],2) * grad2[eltIdx];
        //     grad[eltIdx]= -2/(density[eltIdx]*pow(velocity[eltIdx],3)) * grad[eltIdx];
        //   }
        // } );
      } );
    }

    prepareNextTimestep( mesh );
  } );

  return dtCompute;
}

void AcousticROMFrechet::prepareNextTimestep( MeshLevel & mesh )
{
  NodeManager & nodeManager = mesh.getNodeManager();

  arrayView1d< real32 > const p_nm1 = nodeManager.getField< acousticfields::Pressure_nm1 >();
  arrayView1d< real32 > const p_n   = nodeManager.getField< acousticfields::Pressure_n >();
  arrayView1d< real32 > const p_np1 = nodeManager.getField< acousticfields::Pressure_np1 >();

  arrayView2d< real32 > const pf_nm1 = nodeManager.getField< acousticfields::PressureFrechet_nm1 >();
  arrayView2d< real32 > const pf_n = nodeManager.getField< acousticfields::PressureFrechet_n >();
  arrayView2d< real32 > const pf_np1 = nodeManager.getField< acousticfields::PressureFrechet_np1 >();
  localIndex const ordF = m_orderFrechet;

  arrayView1d< real32 > const stiffnessVector = nodeManager.getField< acousticfields::StiffnessVector >();
  arrayView1d< real32 > const rhs = nodeManager.getField< acousticfields::ForcingRHS >();
  //arrayView1d< real32 > const rhs_fp1 = nodeManager.getField< acousticfields::ForcingRHS_fp1 >();

  SortedArrayView< localIndex const > const solverTargetNodesSet = m_solverTargetNodesSet.toViewConst();

  if( m_solverROM )
  {
    arrayView1d< real32 > const a_nm1 = m_a_nm1.toView();
    arrayView1d< real32 > const a_n = m_a_n.toView();
    arrayView1d< real32 > const a_np1 = m_a_np1.toView();

    
    for( localIndex m = 0; m < a_np1.size(0); ++m )
    {
      a_nm1[m] = a_n[m];
      a_n[m]   = a_np1[m];
    }
  }
  else
  {
    forAll< EXEC_POLICY >( solverTargetNodesSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const n )
    {
      localIndex const a = solverTargetNodesSet[n];

      p_nm1[a] = p_n[a];
      p_n[a]   = p_np1[a];

      stiffnessVector[a] = rhs[a] = 0.0;

      for(localIndex f=0; f<ordF; ++f)
      {
	pf_nm1[a][f] = pf_n[a][f];
	pf_n[a][f]   = pf_np1[a][f];
      }
    } );
  }
}

void AcousticROMFrechet::computeUnknowns( real64 const & time_n,
					  real64 const & dt,
					  DomainPartition & domain,
					  MeshLevel & mesh,
					  arrayView1d< string const > const & regionNames )
{
  EventManager const & event = getGroupByPath< EventManager >( "/Problem/Events" );
  real64 const & minTime = event.getReference< real64 >( EventManager::viewKeyStruct::minTimeString() );
  localIndex const cycleNumber = int(round(time_n/dt));
  integer const cycleForSource = int(round( -minTime / dt + cycleNumber ));

  /// calculate your time integrators
  real64 const dt2 = pow( dt, 2 );

  if( m_solverROM )
  {
    arrayView1d< real32 > const rhs = m_rhsPOD.toView();
    addSourceToRightHandSide( cycleForSource, rhs );
    
    arrayView1d< real32 > const a_nm1 = m_a_nm1.toView();
    arrayView1d< real32 > const a_n = m_a_n.toView();
    arrayView1d< real32 > const a_np1 = m_a_np1.toView();

    //GEOS_MARK_SCOPE ( updateP );

    real32 const alpha = m_alpha;
    arrayView2d< real32 const > const damping = m_dampingPOD.toViewConst();
    arrayView2d< real32 const > const mass = m_massPOD.toViewConst();
    arrayView2d< real32 const > const massPerturbation = m_massPerturbationPOD.toViewConst();
    arrayView2d< real32 const > const dampingPerturbation = m_dampingPerturbationPOD.toViewConst();

    arrayView2d< real64 const > const Op = m_OpPOD.toViewConst();

    array1d< real32 > const b(m_sizePOD);
    arrayView1d< real32 > const bV = b.toView();

    for( localIndex m = 0; m < bV.size(0); ++m )
    {
      bV[m] = dt2 * ( rhs[m] - a_n[m]);
      for( localIndex n = 0; n < bV.size(0); ++n )
      {
        bV[m] += 2 * (mass[m][n] + alpha * massPerturbation[m][n]) * a_n[n];
	bV[m] -= (mass[m][n] + alpha * massPerturbation[m][n] - dt * 0.5 * (damping[m][n] + alpha * dampingPerturbation[m][n])) * a_nm1[n];
      }
      rhs[m] = 0.0;
    }
    for( localIndex m = 0; m < a_np1.size(0); ++m )
    {
      a_np1[m] = 0;
      for( localIndex n = 0; n < a_np1.size(0); ++n )
      {
	a_np1[m] += Op[m][n] * bV[n];
      }
    }
  }

  else
  {
    NodeManager & nodeManager = mesh.getNodeManager();

    arrayView1d< real32 const > const mass = nodeManager.getField< acousticfields::AcousticMassVector >();
    arrayView1d< real32 const > const damping = nodeManager.getField< acousticfields::DampingVector >();
    arrayView1d< real32 const > const massPerturbation = nodeManager.getField< acousticfields::AcousticMassPerturbationVector >();
    arrayView1d< real32 const > const dampingPerturbation = nodeManager.getField< acousticfields::DampingPerturbationVector >();

    arrayView1d< real32 > const p_nm1 = nodeManager.getField< acousticfields::Pressure_nm1 >();
    arrayView1d< real32 > const p_n = nodeManager.getField< acousticfields::Pressure_n >();
    arrayView1d< real32 > const p_np1 = nodeManager.getField< acousticfields::Pressure_np1 >();

    arrayView1d< localIndex const > const freeSurfaceNodeIndicator = nodeManager.getField< acousticfields::AcousticFreeSurfaceNodeIndicator >();
    arrayView1d< real32 > const stiffnessVector = nodeManager.getField< acousticfields::StiffnessVector >();
    arrayView1d< real32 > const rhs = nodeManager.getField< acousticfields::ForcingRHS >();
    //arrayView1d< real32 > const rhs_fp1 = nodeManager.getField< acousticfields::ForcingRHS_fp1 >();
        
    localIndex const ordF = m_orderFrechet;
    localIndex const ordGS = m_orderGS;
    arrayView2d< wsCoordType const, nodes::REFERENCE_POSITION_USD > const
      nodeCoords32 = nodeManager.getField< fields::referencePosition32 >().toViewConst();

    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
											  CellElementSubRegion & elementSubRegion )
    {
      finiteElement::FiniteElementBase const &
	fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

      arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();
      //arrayView1d< real32 const > const velocity = elementSubRegion.getField< acousticfields::AcousticVelocity >();
      //arrayView1d< real32 const > const perturbation = elementSubRegion.getField< acousticfields::Perturbation >();
      
      finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
      {
	using FE_TYPE = TYPEOFREF( finiteElement );
	/*
	acousticROMFrechetKernels::computeStiffnessFrechetRhs::launch< EXEC_POLICY, ATOMIC_POLICY, FE_TYPE >( elementSubRegion.size(),
													      nodeCoords32,
													      elemsToNodes,
													      p_n,
													      rhs_fp1,
													      stiffnessVector,
													      perturbation,
													      velocity);
	*/
	acousticROMFrechetKernels::computeStiffnessFrechet::launch< EXEC_POLICY, ATOMIC_POLICY, FE_TYPE >( elementSubRegion.size(),
													   nodeCoords32,
													   elemsToNodes,
													   p_n,
													   stiffnessVector);
	
      } );
    } );

    addSourceToRightHandSide( cycleForSource, rhs );
    SortedArrayView< localIndex const > const solverTargetNodesSet = m_solverTargetNodesSet.toViewConst();
    
    GEOS_MARK_SCOPE ( updateP );
    forAll< EXEC_POLICY >( solverTargetNodesSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const n )
    {
      localIndex const a = solverTargetNodesSet[n];
      if( freeSurfaceNodeIndicator[a] != 1 )
      {
	p_np1[a] = p_n[a];
	p_np1[a] *= 2.0 * mass[a];
	p_np1[a] -= (mass[a] - 0.5 * dt * damping[a]) * p_nm1[a];
	p_np1[a] += dt2 * (rhs[a] - stiffnessVector[a]);
	p_np1[a] /= mass[a] + 0.5 * dt * damping[a];
      }
    } );

    if( ordGS >= 0)
    {
      mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
											    CellElementSubRegion & elementSubRegion )
      {
	finiteElement::FiniteElementBase const &
	  fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

	arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();
	arrayView1d< integer const > const nodeGhostRank = nodeManager.ghostRank();
	if( cycleForSource%20 == 0 || m_epsilonGS == 0 )
	{
	  bool success = gramSchmidtROMStiffness(fe,
						 stiffnessVector,
						 p_n,
						 nodeGhostRank,
						 elementSubRegion.size(),
						 elemsToNodes,
						 nodeCoords32,
						 0);
	  if( success )
    	  {
	    localIndex nq = m_sizePOD - 1;
	    m_selectionOrder[0][nq] = 0;
	    m_selectionOrder[1][nq] = m_sizePOD_f[0];
	    
	    if( m_sizePOD_f[0]%20 == 0 )
	    {
	      int const rank = MpiWrapper::commRank( MPI_COMM_GEOS );
	      std::string path = GEOS_FMT( "phi/shot_{:05}/rank_{:05}/order_{:02}/", m_shotIndex, rank, 0);

	      localIndex fail = reorthogonalization(fe,
						    nodeGhostRank,
						    elementSubRegion.size(),
						    elemsToNodes,
						    nodeCoords32,
						    m_sizePOD_f[0],
						    path);
	      
	      if(fail > 0)
	      {
		int count=0;
		for(int i=nq; i>=0; --i)
	        {
		  if(m_selectionOrder[0][i] == 0)
		  {
		    m_selectionOrder[1][i] = 0;
		    m_selectionOrder[0][i] = -1;
		    ++count;
		  }
		  if(count == fail)
		  {
		    break;
		  }
		}
		m_sizePOD_f[0] -= fail; 
	      }
	    }
	    m_cycleOrder[0][m_sizePOD_f[0]] = cycleForSource;
	  }
	}
      } );
    }
    for( localIndex f=0; f<ordF; ++f )
    {
      arrayView2d< real32 > const pf_nm1 = nodeManager.getField< acousticfields::PressureFrechet_nm1 >();
      arrayView2d< real32 > const pf_n = nodeManager.getField< acousticfields::PressureFrechet_n >();
      arrayView2d< real32 > const pf_np1 = nodeManager.getField< acousticfields::PressureFrechet_np1 >();

      array1d< real32 > pf;
      pf.resizeWithoutInitializationOrDestruction(LvArray::MemorySpace::cuda, pf_n.size( 0 ));
      arrayView1d< real32 > pfV = pf.toView();
      
      array1d< real32 > p_dt2;
      p_dt2.resizeWithoutInitializationOrDestruction(LvArray::MemorySpace::cuda, pf_n.size( 0 ));
      arrayView1d< real32 > p_dt2V = p_dt2.toView();
      
      forAll< EXEC_POLICY >( solverTargetNodesSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const n )
      {
	localIndex const a = solverTargetNodesSet[n];
	if( freeSurfaceNodeIndicator[a] != 1 )
	{
	  pfV[a] = pf_n[a][f];
	  if(f==0)
	  {
	    p_dt2V[a] = (p_np1[a] - 2*p_n[a] + p_nm1[a])/dt2;
	  }
	  else
	  {
	    p_dt2V[a] = (pf_np1[a][f-1] -2*pf_n[a][f-1] + pf_nm1[a][f-1])/dt2;
	  }
	}
	
	stiffnessVector[a] = 0.0;
	//rhs[a] = rhs_fp1[a];
        //rhs_fp1[a] *= -(f+1);
      } );
      mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
											    CellElementSubRegion & elementSubRegion )
      {
	finiteElement::FiniteElementBase const &
	  fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
	
	arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();
	arrayView1d< integer const > const nodeGhostRank = nodeManager.ghostRank();
	
	//arrayView1d< real32 const > const velocity = elementSubRegion.getField< acousticfields::AcousticVelocity >();
        //arrayView1d< real32 const > const perturbation = elementSubRegion.getField< acousticfields::Perturbation >();

        finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
	{
	  using FE_TYPE = TYPEOFREF( finiteElement );
	  /*
	  acousticROMFrechetKernels::computeStiffnessFrechetRhs::launch< EXEC_POLICY, ATOMIC_POLICY, FE_TYPE >( elementSubRegion.size(),
														nodeCoords32,
														elemsToNodes,
														pfV,
														rhs_fp1,
														stiffnessVector,
														perturbation,
														velocity);
	  */
	  acousticROMFrechetKernels::computeStiffnessFrechet::launch< EXEC_POLICY, ATOMIC_POLICY, FE_TYPE >( elementSubRegion.size(),
													     nodeCoords32,
													     elemsToNodes,
													     pfV,
                                                                                                             stiffnessVector);
	  
	} );

	forAll< EXEC_POLICY >( solverTargetNodesSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const n )
	{
	  localIndex const a = solverTargetNodesSet[n];
	  if( freeSurfaceNodeIndicator[a] != 1 )
	  {
	    rhs[a] = -(f+1) * massPerturbation[a] * p_dt2V[a];
	    pf_np1[a][f] = pf_n[a][f];
	    pf_np1[a][f] *= 2.0 * mass[a];
	    pf_np1[a][f] -= (mass[a] - 0.5 * dt * damping[a]) * pf_nm1[a][f];
	    pf_np1[a][f] += dt2 * (rhs[a] - stiffnessVector[a]);
	    pf_np1[a][f] /= mass[a] + 0.5 * dt * damping[a];
	  }
	} );
    
	if( ordGS >= f+1 )
	{
	  if( cycleForSource%20 == 0 || m_epsilonGS == 0 )
	  {
	    bool success = gramSchmidtROMStiffness(fe,
						   stiffnessVector,
						   pfV,
						   nodeGhostRank,
						   elementSubRegion.size(),
						   elemsToNodes,
						   nodeCoords32,
						   f+1);
	    if( success )
	    {
	      localIndex nq = m_sizePOD - 1;
	      m_selectionOrder[0][nq] = f+1;
	      m_selectionOrder[1][nq] = m_sizePOD_f[f+1];
	      
	      if( m_sizePOD_f[f+1]%20 == 0 )
              {
		int const rank = MpiWrapper::commRank( MPI_COMM_GEOS );
		std::string path = GEOS_FMT( "phi/shot_{:05}/rank_{:05}/order_{:02}/", m_shotIndex, rank, f+1);
		
		localIndex fail = reorthogonalization(fe,
						      nodeGhostRank,
						      elementSubRegion.size(),
						      elemsToNodes,
						      nodeCoords32,
						      m_sizePOD_f[f+1],
						      path);

		if(fail > 0)
		{
		  int count=0;
		  for(int i=nq; i>=0; --i)
		  {
		    if(m_selectionOrder[0][i] == f+1)
		    {
		      m_selectionOrder[1][i] = 0;
		      ++count;
		    }
		    if(count == fail)
		    {
		      break;
		    }
		  }	
		  m_sizePOD_f[f+1] -= fail;
		}
	      }
	      m_cycleOrder[0][m_sizePOD_f[f+1]] = cycleForSource;
	    }
	  }
	}
      } );
    }

    real64 const & maxTime = event.getReference< real64 >( EventManager::viewKeyStruct::maxTimeString() );
    if( cycleNumber == int(maxTime / dt) - 1 && m_orderGS >= 0)
    {
      
      arrayView1d< integer const > const nodeGhostRank = nodeManager.ghostRank();
      
      mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
											    CellElementSubRegion & elementSubRegion )
      {
	finiteElement::FiniteElementBase const &
	  fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

	arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();
	/*
	gramSchmidtROMStiffnessFinal(fe,
				     nodeGhostRank,
				     elementSubRegion.size(),
				     elemsToNodes,
				     nodeCoords32);
	*/
	modifiedGramSchmidtROMStiffness(fe,
					nodeGhostRank,
					elementSubRegion.size(),
					elemsToNodes,
					nodeCoords32);
	
      } );
      
      if( MpiWrapper::commRank( MPI_COMM_GEOS ) == 0)
      {
	std::cout<<"Size final basis = "<<m_sizePOD<<std::endl;
      }
      
      m_massPOD.resize( m_sizePOD, m_sizePOD );
      m_dampingPOD.resize( m_sizePOD, m_sizePOD );
      m_massPerturbationPOD.resize( m_sizePOD, m_sizePOD );
      m_dampingPerturbationPOD.resize( m_sizePOD, m_sizePOD );
      m_sourceConstantsPOD.resize(m_sourceConstants.size(0), m_sizePOD );
      m_receiverConstantsPOD.resize(m_receiverConstants.size(0), m_sizePOD );

      m_massPOD.zero();
      m_dampingPOD.zero();
      m_massPerturbationPOD.zero();
      m_dampingPerturbationPOD.zero();
      m_sourceConstantsPOD.zero();
      m_receiverConstantsPOD.zero();
      
      arrayView2d< real32 > const dampingPOD = m_dampingPOD.toView();
      arrayView2d< real32 > const massPOD = m_massPOD.toView();
      arrayView2d< real32 > const massPerturbationPOD = m_massPerturbationPOD.toView();
      arrayView2d< real32 > const dampingPerturbationPOD = m_dampingPerturbationPOD.toView();
      arrayView2d< real64 > const sourceConstantsPOD = m_sourceConstantsPOD.toView();
      arrayView2d< real64 > const receiverConstantsPOD = m_receiverConstantsPOD.toView();
      
      computeReducedMatrices( massPOD,
			      massPerturbationPOD,
			      dampingPOD,
			      dampingPerturbationPOD,
			      sourceConstantsPOD,
			      receiverConstantsPOD,
			      mass,
			      massPerturbation,
			      damping,
			      dampingPerturbation,
			      m_sourceConstants,
			      m_sourceNodeIds,
			      m_sourceIsAccessible,
			      m_receiverConstants,
			      m_receiverNodeIds,
			      m_receiverIsLocal,
			      nodeGhostRank );
      
    }
  }
}


void AcousticROMFrechet::synchronizeUnknowns( real64 const & time_n,
					      real64 const & dt,
					      DomainPartition & domain,
					      MeshLevel & mesh,
					      arrayView1d< string const > const & )
{
  NodeManager & nodeManager = mesh.getNodeManager();

  /// compute the seismic traces since last step.
  arrayView2d< real32 > const pReceivers = m_pressureNp1AtReceivers.toView();
  
  /// synchronize pressure fields
  FieldIdentifiers fieldsToBeSync;
  if(m_solverROM == 0)
  {
    
    arrayView1d< real32 > const p_n = nodeManager.getField< acousticfields::Pressure_n >();
    arrayView1d< real32 > const p_np1 = nodeManager.getField< acousticfields::Pressure_np1 >();
    
    fieldsToBeSync.addFields( FieldLocation::Node, { acousticfields::Pressure_np1::key() } );
    if(m_orderFrechet > 0)
    {
      fieldsToBeSync.addFields( FieldLocation::Node, { acousticfields::PressureFrechet_np1::key() } );
    }

    CommunicationTools & syncFields = CommunicationTools::getInstance();
    syncFields.synchronizeFields( fieldsToBeSync,
				  mesh,
				  domain.getNeighbors(),
				  true );
    
    computeAllSeismoTraces( time_n, dt, p_np1, p_n, pReceivers );
  }
  else
  {
    arrayView1d< real32 const > const a_n = m_a_n.toViewConst();
    arrayView1d< real32 const > const a_np1 = m_a_np1.toViewConst();
    
    arrayView2d< real64 const > const receiverConstants = m_receiverConstantsPOD.toViewConst();
    arrayView1d< localIndex const > const receiverIsLocal = m_receiverIsLocal.toViewConst();
    
    integer const dir = m_forward ? +1 : -1;
    integer const beginIndex = m_forward ? m_indexSeismoTrace : m_nsamplesSeismoTrace-m_indexSeismoTrace;
    for( localIndex iSeismo = beginIndex; iSeismo < m_nsamplesSeismoTrace; iSeismo++ )
    {
      localIndex seismoIndex = m_forward ? iSeismo : m_nsamplesSeismoTrace-iSeismo;
      real64 const timeSeismo = m_dtSeismoTrace * seismoIndex;
      if( dir * timeSeismo > dir * time_n + epsilonLoc )
	break;
      
      computeSeismoTracePOD( time_n,
			     dir * dt,
			     timeSeismo,
			     seismoIndex,
			     receiverConstants,
			     receiverIsLocal,
			     m_nsamplesSeismoTrace,
			     a_np1,
			     a_n,
			     pReceivers );
    }
  }
  incrementIndexSeismoTrace( time_n );
}



real64 AcousticROMFrechet::explicitStepInternal( real64 const & time_n,
						 real64 const & dt,
						 DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  GEOS_LOG_RANK_0_IF( dt < epsilonLoc, "Warning! Value for dt: " << dt << "s is smaller than local threshold: " << epsilonLoc );

  real64 dtCompute;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & regionNames )
  {
    localIndex nSubSteps = (int) ceil( dt/m_timeStep );
    dtCompute = dt/nSubSteps;
    computeUnknowns( time_n, dtCompute, domain, mesh, regionNames );
    synchronizeUnknowns( time_n, dtCompute, domain, mesh, regionNames );
  } );

  return dtCompute;
}

void AcousticROMFrechet::cleanup( real64 const time_n,
				  integer const cycleNumber,
				  integer const eventCounter,
				  real64 const eventProgress,
				  DomainPartition & domain )
{
  // call the base class cleanup (for reporting purposes)
  PhysicsSolverBase::cleanup( time_n, cycleNumber, eventCounter, eventProgress, domain );

  // compute the remaining seismic traces, if needed
  arrayView2d< real32 > const pReceivers = m_pressureNp1AtReceivers.toView();
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                arrayView1d< string const > const & )
  {
    if( m_solverROM == 0 )
    {
      NodeManager & nodeManager = mesh.getNodeManager();
      arrayView1d< real32 const > const p_n = nodeManager.getField< acousticfields::Pressure_n >();
      arrayView1d< real32 const > const p_np1 = nodeManager.getField< acousticfields::Pressure_np1 >();
      computeAllSeismoTraces( time_n, 0.0, p_np1, p_n, pReceivers );
    }
    else
    {
      arrayView1d< real32 const > const a_n = m_a_n.toViewConst();
      arrayView1d< real32 const > const a_np1 = m_a_np1.toViewConst();
      arrayView2d< real64 const > const receiverConstants = m_receiverConstantsPOD.toViewConst();
      arrayView1d< localIndex const > const receiverIsLocal = m_receiverIsLocal.toViewConst();
      integer const dir = m_forward ? +1 : -1;
      integer const beginIndex = m_forward ? m_indexSeismoTrace : m_nsamplesSeismoTrace-m_indexSeismoTrace;
      for( localIndex iSeismo = beginIndex; iSeismo < m_nsamplesSeismoTrace; iSeismo++ )
      {
	localIndex seismoIndex = m_forward ? iSeismo : m_nsamplesSeismoTrace-iSeismo;
	real64 const timeSeismo = m_dtSeismoTrace * seismoIndex;
	if( dir * timeSeismo > dir * time_n + epsilonLoc )
	  break;
	computeSeismoTracePOD( time_n,
			       0.0,
			       timeSeismo,
			       seismoIndex,
			       receiverConstants,
			       receiverIsLocal,
			       m_nsamplesSeismoTrace,
			       a_np1,
			       a_n,
			       pReceivers );
      }
    }
    WaveSolverUtils::writeSeismoTrace( "seismoTraceReceiver", getName(), m_outputSeismoTrace, m_receiverConstants.size( 0 ),
                                       m_receiverIsLocal, m_nsamplesSeismoTrace, pReceivers );
  } );
}

bool AcousticROMFrechet::gramSchmidtROMStiffness(finiteElement::FiniteElementBase const & fe,
						 arrayView1d< real32 const > const Ku,
						 arrayView1d< real32 const > const u,
						 arrayView1d< integer const > const nodeghostrank,
						 localIndex const elemRegionSize,
						 arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
						 arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const X,
						 localIndex const ordF)
{
  localIndex const size = Ku.size();
  RAJA::ReduceSum< parallelDeviceReduce, real64 > val( 0.0 );
  bool success = true;
  array1d< real32 > const q(size);
  array1d< real32 > const q_new(size);
  array1d< real32 > const stiffnessVector_q(size);
  stiffnessVector_q.zero();

  arrayView1d< real32 > const qV = q.toView();
  arrayView1d< real32 > const q_newV = q_new.toView();
  arrayView1d< real32 > const stiffnessVector_qV = stiffnessVector_q.toView();
  real32 const eps = m_epsilonGS;
  integer const shotIndex = m_shotIndex;

  forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
  {
    if (nodeghostrank[a]< 0)
    {
      val += u[a]*Ku[a];
    }
    q_newV[a] = u[a];
  } );

  real64 normK = MpiWrapper::sum(val.get());
  real64 val_all = normK;
  int nq = m_sizePOD_f[ordF];
  //int nq = m_sizePOD;
  for(int iq=nq; iq>0; --iq)
  {
    GEOS_MARK_SCOPE ( DirectRead );
    int const rank = MpiWrapper::commRank( MPI_COMM_GEOS );
    std::string fileName = GEOS_FMT( "phi/shot_{:05}/rank_{:05}/order_{:02}/vector_{:03}.dat", shotIndex, rank, ordF, iq);
    //std::string fileName = GEOS_FMT( "phi/shot_{:05}/finalBases/rank_{:05}/vector_{:03}.dat", shotIndex, rank, iq);
    std::ifstream wf( fileName, std::ios::in | std::ios::binary );
    GEOS_THROW_IF( !wf,
		   getDataContext() << ": Could not open file "<< fileName << " for reading",
		   InputError );
    qV.move( LvArray::MemorySpace::host, true );
    wf.read( (char *)&qV[0], size*sizeof( real32 ) );
    wf.close( );

    RAJA::ReduceSum< parallelDeviceReduce, real64 > valq( 0.0 );
    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      if( nodeghostrank[a]< 0 )
      {
	valq += qV[a]*Ku[a];
      }
    } );
    real64 valq_all = MpiWrapper::sum(valq.get());

    val_all -= pow(valq_all,2);
    if( val_all > eps*normK )
    {
      forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
      {
        q_newV[a] -= valq_all * qV[a];
      } );
    }
    else
    {
      success = false;
      return success;
    }
  }

  if( success==true && val_all > eps*normK && val_all > epsilonLoc )
  {
    finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
    {
      using FE_TYPE = TYPEOFREF( finiteElement );
      acousticROMFrechetKernels::computeStiffnessFrechet::launch< EXEC_POLICY, ATOMIC_POLICY, FE_TYPE >( elemRegionSize,
													 X,
													 elemsToNodes,
													 q_newV,
													 stiffnessVector_qV);
    } );
    RAJA::ReduceSum< parallelDeviceReduce, real64 > val3( 0.0 );
    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      if( nodeghostrank[a]< 0)
      {
	val3 += q_newV[a] * stiffnessVector_qV[a];
      }
    } );
    normK=MpiWrapper::sum(val3.get());
    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      q_newV[a] /= sqrt(normK);
    } );
    GEOS_MARK_SCOPE ( DirectWrite );
    int const rank = MpiWrapper::commRank( MPI_COMM_GEOS );
    q_newV.move( LvArray::MemorySpace::host, false );

    std::string fileName = GEOS_FMT( "phi/shot_{:05}/rank_{:05}/order_{:02}/vector_{:03}.dat", shotIndex, rank, ordF, nq+1);
    //std::string fileName = GEOS_FMT( "phi/shot_{:05}/finalBases/rank_{:05}/vector_{:03}.dat", shotIndex, rank, nq+1);
    int lastDirSeparator = fileName.find_last_of( "/\\" );
    std::string dirName = fileName.substr( 0, lastDirSeparator );
    if( string::npos != (size_t)lastDirSeparator && !directoryExists( dirName ))
    {
      makeDirsForPath( dirName );
    }
    std::ofstream wf( fileName, std::ios::out | std::ios::binary );
    GEOS_THROW_IF( !wf,
		   getDataContext() << ": Could not open file "<< fileName << " for writing",
		   InputError );
    wf.write( (char *)&q_newV[0], size*sizeof( real32 ) );
    wf.close( );
    GEOS_THROW_IF( !wf.good(),
		   getDataContext() << ": An error occured while writing "<< fileName,
		   InputError );
    m_sizePOD_f[ordF] += 1;
    m_sizePOD += 1;

    /*
    if( m_sizePOD_f[ordF]%20 == 0 )
      //if( m_sizePOD%20 == 0 )
    {
      std::string path = GEOS_FMT( "phi/shot_{:05}/rank_{:05}/order_{:02}/", shotIndex, rank, ordF);
      //std::string path = GEOS_FMT( "phi/shot_{:05}/finalBases/rank_{:05}/", shotIndex, rank);
      
      localIndex fail = reorthogonalization(fe,
					    nodeghostrank,
					    elemRegionSize,
					    elemsToNodes,
					    X,
					    m_sizePOD_f[ordF],
					    path);

      int count=0;
      for(int i=m_sizePOD-1; i>=0; --i)
      {
	if(m_selectionOrder[0][i] == ordF)
	{
	  std::cout<<i<<" "<<m_selectionOrder[1][i]
	  m_selectionOrder[1][i] = 0;
	  ++count;
	}
	if(count == fail)
	{
	  break;
	}
      }
      m_sizePOD_f[ordF] -= fail;
 
    }
    */
    return success;
  }
  else
  {
    success = false;
    return success;
  }
}


void AcousticROMFrechet::gramSchmidtROMStiffnessFinal(finiteElement::FiniteElementBase const & fe,
						      arrayView1d< integer const > const nodeghostrank,
						      localIndex const elemRegionSize,
						      arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
						      arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const X)
{
  int size = nodeghostrank.size();

  array1d< real32 > const q1(size);
  array1d< real32 > const q2(size);
  array1d< real32 > const stiffnessVector_q(size);
  stiffnessVector_q.zero();

  arrayView1d< real32 > const q1V = q1.toView();
  arrayView1d< real32 > const q2V = q2.toView();
  arrayView1d< real32 > const stiffnessVector_qV = stiffnessVector_q.toView();
  localIndex shotIndex = m_shotIndex;
  real32 eps = m_epsilonGS;
  localIndex totCount = m_sizePOD;
  localIndex jf = 0;
  localIndex fail = 0;
  
  for( int j=0; j<totCount; ++j )
  {
    bool success = true;
    localIndex jordF = m_selectionOrder[0][j];
    localIndex jq = m_selectionOrder[1][j];

    GEOS_MARK_SCOPE ( DirectReadWrite );
    int const rank = MpiWrapper::commRank( MPI_COMM_GEOS );
    std::string fileName1 = GEOS_FMT( "phi/shot_{:05}/rank_{:05}/order_{:02}/vector_{:03}.dat", shotIndex, rank, jordF, jq);
    std::ifstream wf1( fileName1, std::ios::in | std::ios::binary );
    GEOS_THROW_IF( !wf1,
                   getDataContext() << ": Could not open file "<< fileName1 << " for reading",
                   InputError );
    q1V.move( LvArray::MemorySpace::host, true );
    wf1.read( (char *)&q1V[0], size*sizeof( real32 ) );
    wf1.close( );

    finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
    {
      using FE_TYPE = TYPEOFREF( finiteElement );
      acousticROMFrechetKernels::computeStiffnessFrechet::launch< EXEC_POLICY, ATOMIC_POLICY, FE_TYPE >( elemRegionSize,
                                                                                                         X,
                                                                                                         elemsToNodes,
                                                                                                         q1V,
                                                                                                         stiffnessVector_qV);
    } );
    real64 val_all = 1.0;

    for(int iq=jf;iq>0;--iq)
    {
      std::string fileName2 = GEOS_FMT( "phi/shot_{:05}/finalBases/rank_{:05}/vector_{:03}.dat", shotIndex, rank, iq);
      std::ifstream wf2( fileName2, std::ios::in | std::ios::binary );
      GEOS_THROW_IF( !wf2,
                     getDataContext() << ": Could not open file "<< fileName2 << " for reading",
                     InputError );
      q2V.move( LvArray::MemorySpace::host, true );
      wf2.read( (char *)&q2V[0], size*sizeof( real32 ) );
      wf2.close( );

      RAJA::ReduceSum< parallelDeviceReduce, real64 > valq( 0.0 );
      forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
      {
        if (nodeghostrank[a]< 0)
        {
          valq += q2V[a]*stiffnessVector_qV[a];
        }
      } );
      real64 valq_all = MpiWrapper::sum(valq.get());
      val_all -= pow(valq_all,2);
      if (val_all <= eps)
      {
	success = false;
	break;
      }
      forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
      {
        q1V[a] -= valq_all * q2V[a];
      } );
    }

    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      stiffnessVector_qV[a] = 0.0;
    } );

    if (success==true and val_all > eps)
    {
      finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
      {
	using FE_TYPE = TYPEOFREF( finiteElement );
	acousticROMFrechetKernels::computeStiffnessFrechet::launch< EXEC_POLICY, ATOMIC_POLICY, FE_TYPE >( elemRegionSize,
													   X,
													   elemsToNodes,
													   q1V,
													   stiffnessVector_qV);
      } );

      RAJA::ReduceSum< parallelDeviceReduce, real64 > val3( 0.0 );
      forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
      {
	if( nodeghostrank[a]< 0)
	{
	  val3 += q1V[a] * stiffnessVector_qV[a];
	}
      } );
      real64 normK = MpiWrapper::sum(val3.get());
      forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
      {
	q1V[a] /= sqrt(normK);
	stiffnessVector_qV[a] = 0.0;
      } );

      q1V.move( LvArray::MemorySpace::host, false );

      jf += 1;
      std::string fileName = GEOS_FMT( "phi/shot_{:05}/finalBases/rank_{:05}/vector_{:03}.dat", shotIndex, rank, jf);
      int lastDirSeparator = fileName.find_last_of( "/\\" );
      std::string dirName = fileName.substr( 0, lastDirSeparator );
      if( string::npos != (size_t)lastDirSeparator && !directoryExists( dirName ))
      {
	makeDirsForPath( dirName );
      }
      std::ofstream wf( fileName, std::ios::out | std::ios::binary );
      GEOS_THROW_IF( !wf,
		     getDataContext() << ": Could not open file "<< fileName << " for writing",
		     InputError );
      wf.write( (char *)&q1V[0], size*sizeof( real32 ) );
      wf.close( );
      GEOS_THROW_IF( !wf.good(),
		     getDataContext() << ": An error occured while writing "<< fileName,
		     InputError );

      if( jf%20 == 0 )
      {
        std::string path = GEOS_FMT( "phi/shot_{:05}/finalBases/rank_{:05}/", shotIndex, rank );
	int fail2 = reorthogonalization(fe,
					nodeghostrank,
					elemRegionSize,
					elemsToNodes,
					X,
					jf,
					path);
	jf -= fail2;
	m_sizePOD -= fail2;
      }

    }
    else
    {
      fail += 1;
    }
    remove( fileName1.c_str() );
  }
  m_sizePOD -= fail;
  int const rank = MpiWrapper::commRank( MPI_COMM_GEOS );
  std::string directory = GEOS_FMT( "phi/shot_{:05}/rank_{:05}", shotIndex, rank);
  remove( directory.c_str() );
}

void AcousticROMFrechet::modifiedGramSchmidtROMStiffness(finiteElement::FiniteElementBase const & fe,
							 arrayView1d< integer const > const nodeghostrank,
							 localIndex const elemRegionSize,
							 arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
							 arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const X)
{
  int size = nodeghostrank.size();

  array1d< real32 > const q1(size);
  array1d< real32 > const q2(size);
  array1d< real32 > const stiffnessVector_q(size);
  stiffnessVector_q.zero();

  arrayView1d< real32 > const q1V = q1.toView();
  arrayView1d< real32 > const q2V = q2.toView();
  arrayView1d< real32 > const stiffnessVector_qV = stiffnessVector_q.toView();
  localIndex shotIndex = m_shotIndex;
  real32 eps = m_epsilonGS;
  localIndex totCount = m_sizePOD;
  localIndex jf = 0;
  
  for( int j=0; j<totCount; ++j )
  {
    localIndex jordF = m_selectionOrder[0][j];
    localIndex jq = m_selectionOrder[1][j];
    if(jq == 0)
    {
      continue;
    }
    
    GEOS_MARK_SCOPE ( DirectReadWrite );
    int const rank = MpiWrapper::commRank( MPI_COMM_GEOS );
    std::string fileName1 = GEOS_FMT( "phi/shot_{:05}/rank_{:05}/order_{:02}/vector_{:03}.dat", shotIndex, rank, jordF, jq);
    std::ifstream wf1( fileName1, std::ios::in | std::ios::binary );
    GEOS_THROW_IF( !wf1,
                   getDataContext() << ": Could not open file "<< fileName1 << " for reading",
                   InputError );
    q1V.move( LvArray::MemorySpace::host, true );
    wf1.read( (char *)&q1V[0], size*sizeof( real32 ) );
    wf1.close( );

    finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
    {
      using FE_TYPE = TYPEOFREF( finiteElement );
      acousticROMFrechetKernels::computeStiffnessFrechet::launch< EXEC_POLICY, ATOMIC_POLICY, FE_TYPE >( elemRegionSize,
                                                                                                         X,
                                                                                                         elemsToNodes,
                                                                                                         q1V,
                                                                                                         stiffnessVector_qV);
    } );

    RAJA::ReduceSum< parallelDeviceReduce, real64 > valq( 0.0 );
    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      if( nodeghostrank[a]< 0 )
      {
        valq += q1V[a] * stiffnessVector_qV[a];
      }
    } );
    real64 valq_all = MpiWrapper::sum(valq.get());

    if( valq_all < eps )
    {
      continue;
    }

    real64 root = sqrt(valq_all);
    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      q1V[a] /= root;
      stiffnessVector_qV[a] /= root;
    } );

    q1V.move( LvArray::MemorySpace::host, false );

    jf += 1;
    std::string fileName = GEOS_FMT( "phi/shot_{:05}/finalBases/rank_{:05}/vector_{:03}.dat", shotIndex, rank, jf);
    int lastDirSeparator = fileName.find_last_of( "/\\" );
    std::string dirName = fileName.substr( 0, lastDirSeparator );
    if( string::npos != (size_t)lastDirSeparator && !directoryExists( dirName ) )
    {
      makeDirsForPath( dirName );
    }
    std::ofstream wf( fileName, std::ios::out | std::ios::binary );
    GEOS_THROW_IF( !wf,
		   getDataContext() << ": Could not open file "<< fileName << "for writing",
		   InputError );
    wf.write( (char *)&q1V[0], size*sizeof( real32 ) );
    wf.close( );
    GEOS_THROW_IF( !wf.good(),
		   getDataContext() << ": An error occured while writing "<< fileName,
		   InputError );
    
    for(int i=j+1; i<totCount; ++i)
    {
      localIndex iordF = m_selectionOrder[0][i];
      localIndex iq = m_selectionOrder[1][i];
      if(iq == 0)
      {
	continue;
      }
      std::string fileName2 = GEOS_FMT( "phi/shot_{:05}/rank_{:05}/order_{:02}/vector_{:03}.dat", shotIndex, rank, iordF, iq);
      std::ifstream wf2( fileName2, std::ios::in | std::ios::binary );
      GEOS_THROW_IF( !wf2,
                     getDataContext() << ": Could not open file "<< fileName2 << " for reading",
                     InputError );
      q2V.move( LvArray::MemorySpace::host, true );
      wf2.read( (char *)&q2V[0], size*sizeof( real32 ) );
      wf2.close( );

      RAJA::ReduceSum< parallelDeviceReduce, real64 > val3( 0.0 );
      forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
      {
	if( nodeghostrank[a]< 0 )
        {
	  val3 += q2V[a] * stiffnessVector_qV[a];
	}
      } );
      real64 normK = MpiWrapper::sum(val3.get());

      forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
      {
	q2V[a] -= q1V[a] * normK;
      } );

      std::ofstream wf3( fileName2, std::ios::out | std::ios::binary );
      GEOS_THROW_IF( !wf3,
                     getDataContext() << ": Could not open file "<< fileName2 <<" for reading",
                     InputError );
      q2V.move( LvArray::MemorySpace::host, true );
      wf3.write( (char *)&q2V[0], size*sizeof( real32 ) );
      wf3.close( );
    }
    
    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      stiffnessVector_qV[a] = 0.0;
    } );
    remove( fileName1.c_str() );
  }
  m_sizePOD = jf;

  int const rank = MpiWrapper::commRank( MPI_COMM_GEOS );
  localIndex ifile = m_sizePOD+1;
  while(ifile > 0)
  {
    std::string filename = GEOS_FMT( "phi/shot_{:05}/finalBases/rank_{:05}/vector_{:03}.dat", m_shotIndex, rank, ifile);
    std::ifstream wf(filename.c_str());
    if( wf.good() )
    {
      remove( filename.c_str() );
      ifile+=1;
    }
    else
    {
      ifile = -1;
    }
  }
}
  
int AcousticROMFrechet::reorthogonalization(finiteElement::FiniteElementBase const & fe,
					    arrayView1d< integer const > const nodeghostrank,
					    localIndex const elemRegionSize,
					    arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
					    arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const X,
					    localIndex const nq,
					    std::string path)
{
  localIndex const size = nodeghostrank.size();
  array1d< real32 > const q1(size);
  array1d< real32 > const q2(size);
  array1d< real32 > const stiffnessVector_q(size);
  stiffnessVector_q.zero();
  
  arrayView1d< real32 > const q1V = q1.toView();
  arrayView1d< real32 > const q2V = q2.toView();
  arrayView1d< real32 > const stiffnessVector_qV = stiffnessVector_q.toView();
  localIndex jj = 1;
  localIndex fail = 0;
  real32 const eps = m_epsilonGS;
  
  for(int jq=1; jq<nq+1 ;++jq)
  {
    bool success = 1;
    GEOS_MARK_SCOPE ( DirectReadWrite );
    std::string fileName1 = GEOS_FMT( path + "/vector_{:03}.dat", jq);
    std::ifstream wf1( fileName1, std::ios::in | std::ios::binary );
    GEOS_THROW_IF( !wf1,
                   getDataContext() << ": Could not open file "<< fileName1 << " for reading",
                   InputError );
    q1V.move( LvArray::MemorySpace::host, true );
    wf1.read( (char *)&q1V[0], size*sizeof( real32 ) );
    wf1.close( );

    finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
    {
      using FE_TYPE = TYPEOFREF( finiteElement );
      acousticROMFrechetKernels::computeStiffnessFrechet::launch< EXEC_POLICY, ATOMIC_POLICY, FE_TYPE >( elemRegionSize,
													 X,
													 elemsToNodes,
													 q1V,
													 stiffnessVector_qV);
    } );
    real64 val_all = 1.0;

    for(int iq=1;iq<jj;++iq)
    {
      std::string fileName2 = GEOS_FMT( path + "/vector_{:03}.dat", iq);
      std::ifstream wf2( fileName2, std::ios::in | std::ios::binary );
      GEOS_THROW_IF( !wf2,
		     getDataContext() << ": Could not open file "<< fileName2 << " for reading",
		     InputError );
      q2V.move( LvArray::MemorySpace::host, true );
      wf2.read( (char *)&q2V[0], size*sizeof( real32 ) );
      wf2.close( );

      RAJA::ReduceSum< parallelDeviceReduce, real64 > valq( 0.0 );
      forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
      {
	if (nodeghostrank[a]< 0)
	{
	  valq += q2V[a]*stiffnessVector_qV[a];
	}
      } );
      real64 valq_all = MpiWrapper::sum(valq.get());
      val_all -= pow(valq_all,2);

      if(val_all > eps)
      {
	forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
	{
	  q1V[a] -= valq_all * q2V[a];
	} );
      }
      else
      {
	success = 0;
	break;
      }
    }
    
    forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
    {
      stiffnessVector_qV[a] = 0.0;
    } );

    if(success)
    {
      finiteElement::FiniteElementDispatchHandler< SEM_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
      {
	using FE_TYPE = TYPEOFREF( finiteElement );
	acousticROMFrechetKernels::computeStiffnessFrechet::launch< EXEC_POLICY, ATOMIC_POLICY, FE_TYPE >( elemRegionSize,
													   X,
													   elemsToNodes,
													   q1V,
													   stiffnessVector_qV);
      } );
      
      RAJA::ReduceSum< parallelDeviceReduce, real64 > val3( 0.0 );
      forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
      {
	if( nodeghostrank[a]< 0)
        {
	  val3 += q1V[a] * stiffnessVector_qV[a];
	}
      } );
      real64 normK = MpiWrapper::sum(val3.get());

      forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
      {
	q1V[a] /= sqrt(normK);
	stiffnessVector_qV[a] = 0.0;
      } );

      q1V.move( LvArray::MemorySpace::host, false );
      std::string fileName = GEOS_FMT( path + "/vector_{:03}.dat", jj);
      int lastDirSeparator = fileName.find_last_of( "/\\" );
      std::string dirName = fileName.substr( 0, lastDirSeparator );
      if( string::npos != (size_t)lastDirSeparator && !directoryExists( dirName ))
      {
	makeDirsForPath( dirName );
      }
      std::ofstream wf( fileName, std::ios::out | std::ios::binary );
      GEOS_THROW_IF( !wf,
		     getDataContext() << ": Could not open file "<< fileName << " for writing",
		     InputError );
      wf.write( (char *)&q1V[0], size*sizeof( real32 ) );
      wf.close( );
      GEOS_THROW_IF( !wf.good(),
		     getDataContext() << ": An error occured while writing "<< fileName,
		     InputError );
      jj+=1;
    }
    else
    {
      fail+=1;
    }
  }
  for(int jq=jj; jq<nq+1; ++jq)
  {
    std::string fileName = GEOS_FMT( path + "/vector_{:03}.dat", jq);
    remove( fileName.c_str() );
  }
  return fail;
}



void AcousticROMFrechet::computeReducedMatrices( arrayView2d< real32 > const massPOD,
						 arrayView2d< real32 > const massPerturbationPOD,
						 arrayView2d< real32 > const dampingPOD,
						 arrayView2d< real32 > const dampingPerturbationPOD,
						 arrayView2d< real64 > const sourceConstantsPOD,
						 arrayView2d< real64 > const receiverConstantsPOD,
						 arrayView1d< real32 const > const mass,
						 arrayView1d< real32 const > const massPerturbation,
						 arrayView1d< real32 const > const damping,
						 arrayView1d< real32 const > const dampingPerturbation,
						 arrayView2d< real64 const > const sourceConstants,
						 arrayView2d< localIndex const > const sourceNodeIds,
						 arrayView1d< localIndex const > const sourceIsAccessible,
						 arrayView2d< real64 const > const receiverConstants,
                                                 arrayView2d< localIndex const > const receiverNodeIds,
						 arrayView1d< localIndex const > const receiverIsLocal,
						 arrayView1d< localIndex const > const nodesGhostRank)
{
  GEOS_LOG_RANK_0("Computing mass and damping ROM...");
  int size = mass.size();
  array1d< real32 > const phim( size );
  array1d< real32 > const phin( size );
  arrayView1d< real32 > phimV = phim.toView();
  arrayView1d< real32 > phinV = phin.toView();

  GEOS_MARK_SCOPE ( DirectRead );
  int const rank = MpiWrapper::commRank( MPI_COMM_GEOS );
  for( localIndex m=0; m<m_sizePOD; ++m )
  {
    std::string fileName1 = GEOS_FMT( "phi/shot_{:05}/finalBases/rank_{:05}/vector_{:03}.dat", m_shotIndex, rank, m+1);
    std::ifstream wf1( fileName1, std::ios::in | std::ios::binary );
    GEOS_THROW_IF( !wf1,
		   ": Could not open file "<< fileName1 << " for reading",
		   InputError );
    phimV.move( LvArray::MemorySpace::host, true );
    wf1.read( (char *)&phimV[0], size*sizeof( real32 ) );
    wf1.close( );

    for( localIndex n=0; n<m+1; ++n )
    {
      std::string fileName2 = GEOS_FMT( "phi/shot_{:05}/finalBases/rank_{:05}/vector_{:03}.dat", m_shotIndex, rank, n+1);
      std::ifstream wf2( fileName2, std::ios::in | std::ios::binary );
      GEOS_THROW_IF( !wf2,
		     ": Could not open file "<< fileName2 << " for reading",
		     InputError );
      phinV.move( LvArray::MemorySpace::host, true );
      wf2.read( (char *)&phinV[0], size*sizeof( real32 ) );
      wf2.close( );

      forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const a )
      {
	if( nodesGhostRank[a] < 0 )
        {
	  if( abs(phimV[a]) > 1e-12 && abs(phinV[a]) > 1e-12 )
	  {
	    real32 localIncrement = phimV[a] * mass[a] * phinV[a];
	    RAJA::atomicAdd< ATOMIC_POLICY >( &massPOD[m][n], localIncrement );

	    if( massPerturbation[a] != 0 )
	    {
	      localIncrement = phimV[a] * massPerturbation[a] * phinV[a];
              RAJA::atomicAdd< ATOMIC_POLICY >( &massPerturbationPOD[m][n], localIncrement );
            }

	    if( damping[a] != 0 )
	    {
	      localIncrement = phimV[a] * damping[a] * phinV[a];
	      RAJA::atomicAdd< ATOMIC_POLICY >( &dampingPOD[m][n], localIncrement );
	    }

	    if( dampingPerturbation[a] != 0 )
	    {
	      localIncrement = phimV[a] * dampingPerturbation[a] * phinV[a];
	      RAJA::atomicAdd< ATOMIC_POLICY >( &dampingPerturbationPOD[m][n], localIncrement );
	    }
	  }
	}
      } ); // end loop over element
      if( m!=n )
      {
	forAll< EXEC_POLICY >( 1, [=] GEOS_HOST_DEVICE ( localIndex const a )
	{
	  massPOD[n][m] = massPOD[m][n];
	  dampingPOD[n][m] = dampingPOD[m][n];
	  massPerturbationPOD[n][m] = massPerturbationPOD[m][n];
	  dampingPerturbationPOD[n][m] = dampingPerturbationPOD[m][n];
	} );
      }
    }
    for( localIndex isrc = 0; isrc<sourceNodeIds.size(0); ++isrc)
    {
      forAll< EXEC_POLICY >( sourceNodeIds.size(1), [=] GEOS_HOST_DEVICE ( localIndex const a )
      {
	if( sourceIsAccessible[isrc] == 1 && nodesGhostRank[sourceNodeIds[isrc][a]] < 0 )
	{
	  //real64 localIncrement = phimV[sourceNodeIds[isrc][a]] * sourceConstants[isrc][a];
	  sourceConstantsPOD[isrc][m] += phimV[sourceNodeIds[isrc][a]] * sourceConstants[isrc][a];
	}
      } );
    }
    forAll< EXEC_POLICY >( receiverNodeIds.size(0), [=] GEOS_HOST_DEVICE ( localIndex const ircv )
    {
      for( localIndex a = 0; a<receiverNodeIds.size(1); ++a)
      {
	if( receiverIsLocal[ircv] == 1 )
	{
	  //real64 localIncrement = phimV[receiverNodeIds[ircv][a]] * receiverConstants[ircv][a];
	  receiverConstantsPOD[ircv][m] += phimV[receiverNodeIds[ircv][a]] * receiverConstants[ircv][a];
	}
      } 
    } );   
  }

  massPOD.move( LvArray::MemorySpace::host, false );
  massPerturbationPOD.move( LvArray::MemorySpace::host, false );
  dampingPOD.move( LvArray::MemorySpace::host, false );
  dampingPerturbationPOD.move( LvArray::MemorySpace::host, false );
  sourceConstantsPOD.move( LvArray::MemorySpace::host, false );
  receiverConstantsPOD.move( LvArray::MemorySpace::host, false );
  
  std::vector<std::string> matrixNames = {"mass", "massPerturbation", "damping", "dampingPerturbation", "sourceConstants", "receiverConstants"};

  for( localIndex i=0; i<6; ++i )
  {
    std::string fileName = GEOS_FMT( "phi/shot_{:05}/finalBases/rank_{:05}/{:}.dat", m_shotIndex, rank, matrixNames[i]);
    int lastDirSeparator = fileName.find_last_of( "/\\" );
    std::string dirName = fileName.substr( 0, lastDirSeparator );
    if( string::npos != (size_t)lastDirSeparator && !directoryExists( dirName ))
    {
      makeDirsForPath( dirName );
    }
    std::ofstream wf( fileName, std::ios::out | std::ios::binary );
    GEOS_THROW_IF( !wf,
		   getDataContext() << ": Could not open file "<< fileName << " for writing",
		   InputError );

    switch( i )
    {
    case 0:
    {
      wf.write( (char *)&massPOD[0][0], m_sizePOD*m_sizePOD*sizeof( real32 ) );
    }
    case 1:
    {
      wf.write( (char *)&massPerturbationPOD[0][0], m_sizePOD*m_sizePOD*sizeof( real32 ) );
    }
    case 2:
    {
      wf.write( (char *)&dampingPOD[0][0], m_sizePOD*m_sizePOD*sizeof( real32 ) );
    }
    case 3:
    {
      wf.write( (char *)&dampingPerturbationPOD[0][0], m_sizePOD*m_sizePOD*sizeof( real32 ) );
    }
    case 4:
    {
      wf.write( (char *)&sourceConstantsPOD[0][0], sourceConstantsPOD.size(0)*m_sizePOD*sizeof( real64 ) );
    }
    case 5:
    {
      wf.write( (char *)&receiverConstantsPOD[0][0], receiverConstantsPOD.size(0)*m_sizePOD*sizeof( real64 ) );
    }
    }
    wf.close( );
    GEOS_THROW_IF( !wf.good(),
		   getDataContext() << ": An error occured while writing "<< fileName,
		   InputError );
  }



}

void AcousticROMFrechet::computeSeismoTracePOD( real64 const time_n,
						real64 const dt,
						real64 const timeSeismo,
						localIndex iSeismo,
						arrayView2d< real64 const > const receiverConstants,
						arrayView1d< localIndex const > const receiverIsLocal,
						localIndex const nsamplesSeismoTrace,
						arrayView1d< real32 const > const var_np1,
						arrayView1d< real32 const > const var_n,
						arrayView2d< real32 > varAtReceivers )
{
  real64 const time_np1 = time_n + dt;

  real32 const a1 = (LvArray::math::abs( dt ) < WaveSolverBase::epsilonLoc ) ? 1.0 : (time_np1 - timeSeismo)/dt;
  real32 const a2 = 1.0 - a1;
  
  if( nsamplesSeismoTrace > 0 )
  {
    for( localIndex ircv = 0; ircv<receiverIsLocal.size( 0 ); ++ircv )
    {
      if( receiverIsLocal[ircv] == 1 )
      {
	varAtReceivers[iSeismo][ircv] = 0.0;
	real32 vtmp_np1 = 0.0;
	real32 vtmp_n = 0.0;

	for( localIndex m = 0; m<var_np1.size( 0 ); ++m )
	{
	  vtmp_np1 += var_np1[m] * receiverConstants[ircv][m];
	  vtmp_n += var_n[m] * receiverConstants[ircv][m];
	}
	// linear interpolation between the pressure value at time_n and time_(n+1)
	varAtReceivers[iSeismo][ircv] = a1*vtmp_n + a2*vtmp_np1;
      }
    }
  }
}


REGISTER_CATALOG_ENTRY( PhysicsSolverBase, AcousticROMFrechet, string const &, dataRepository::Group * const )

} /* namespace geos */
