/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file HydrofractureSolver.cpp
 *
 */


#include "HydrofractureSolver.hpp"

#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/contact/ContactRelationBase.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "finiteElement/Kinematics.h"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "mesh/SurfaceElementRegion.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/surfaceGeneration/SurfaceGenerator.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"
#include "mesh/ExtrinsicMeshData.hpp"
#include <unordered_set>
#include <array>
#include <algorithm>

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

HydrofractureSolver::HydrofractureSolver( const std::string & name,
                                          Group * const parent ):
  SolverBase( name, parent ),
  m_solidSolverName(),
  m_flowSolverName(),
  m_surfaceGeneratorName(),
  m_regimeTypeOption(),
  m_couplingTypeOption( CouplingTypeOption::FIM ),
  m_solidSolver( nullptr ),
  m_flowSolver( nullptr ),
  m_surfaceGeneratorSolver( nullptr ),
  m_maxNumResolves( 10 )
{
  registerWrapper( viewKeyStruct::solidSolverNameString, &m_solidSolverName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of the solid mechanics solver to use in the hydrofracture solver" );

  registerWrapper( viewKeyStruct::fluidSolverNameString, &m_flowSolverName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of the fluid mechanics solver to use in the hydrofracture solver" );

  //TJ: register the surfaceGenerator solver
  registerWrapper( viewKeyStruct::surfaceGeneratorSolverNameString, &m_surfaceGeneratorName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of the surface generator solver to use in the hydrofracture solver" );

  registerWrapper( viewKeyStruct::regimeTypeOptionString, &m_regimeTypeOption)->
    setInputFlag( InputFlags::REQUIRED)->
    setDescription( "Propagation regime. Valid options:\n* " + EnumStrings< RegimeTypeOption >::concat( "\n* " ));

  registerWrapper( viewKeyStruct::couplingTypeOptionStringString, &m_couplingTypeOption )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Coupling method. Valid options:\n* " + EnumStrings< CouplingTypeOption >::concat( "\n* " ) );

  registerWrapper( viewKeyStruct::contactRelationNameString, &m_contactRelationName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of contact relation to enforce constraints on fracture boundary." );

  registerWrapper( viewKeyStruct::maxNumResolvesString, &m_maxNumResolves )->
    setApplyDefaultValue( 10 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Value to indicate how many resolves may be executed to perform surface generation after the execution of flow and mechanics solver. " );

  m_numResolves[0] = 0;

  m_linearSolverParameters.get().mgr.strategy = "Hydrofracture";
  m_linearSolverParameters.get().mgr.separateComponents = true;
  m_linearSolverParameters.get().mgr.displacementFieldName = keys::TotalDisplacement;
  m_linearSolverParameters.get().dofsPerNode = 3;
}

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
void HydrofractureSolver::RegisterDataOnMesh( dataRepository::Group * const MeshBodies )
{
  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    MeshLevel * meshLevel = Group::group_cast< MeshBody * >( mesh.second )->getMeshLevel( 0 );

    ElementRegionManager * const elemManager = meshLevel->getElemManager();
    elemManager->forElementRegions< SurfaceElementRegion >( [&] ( SurfaceElementRegion * const region )
    {
      region->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion * const subRegion )
      {
        subRegion->registerWrapper< array1d< real64 > >( viewKeyStruct::separationCoeff0String )->
          setRestartFlags( RestartFlags::NO_WRITE );
        subRegion->registerWrapper< array1d< real64 > >( viewKeyStruct::apertureAtFailureString )->
          setApplyDefaultValue( -1.0 )->
          setPlotLevel( PlotLevel::LEVEL_0 );

        subRegion->registerWrapper< array1d< real64 > >( FaceElementSubRegion::viewKeyStruct::dSeparationCoeffdAperString )->
          setRestartFlags( RestartFlags::NO_WRITE );
      } );
    } );
  }
}
#endif

void HydrofractureSolver::ImplicitStepSetup( real64 const & time_n,
                                             real64 const & dt,
                                             DomainPartition & domain )
{
  UpdateDeformationForCoupling( domain );
  m_solidSolver->ImplicitStepSetup( time_n, dt, domain );
  m_flowSolver->ImplicitStepSetup( time_n, dt, domain );

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  mesh.getElemManager()->forElementRegions< SurfaceElementRegion >( [&]( SurfaceElementRegion & faceElemRegion )
  {
    faceElemRegion.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
    {
      arrayView1d< real64 > const &
      separationCoeff0 = subRegion.getReference< array1d< real64 > >( viewKeyStruct::separationCoeff0String );
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

void HydrofractureSolver::ImplicitStepComplete( real64 const & time_n,
                                                real64 const & dt,
                                                DomainPartition & domain )
{
  m_flowSolver->ImplicitStepComplete( time_n, dt, domain );
  m_solidSolver->ImplicitStepComplete( time_n, dt, domain );
}

void HydrofractureSolver::PostProcessInput()
{
  m_solidSolver = this->getParent()->GetGroup< SolidMechanicsLagrangianFEM >( m_solidSolverName );
  GEOSX_ERROR_IF( m_solidSolver == nullptr, this->getName() << ": invalid solid solver name: " << m_solidSolverName );

  m_flowSolver = this->getParent()->GetGroup< FlowSolverBase >( m_flowSolverName );
  GEOSX_ERROR_IF( m_flowSolver == nullptr, this->getName() << ": invalid flow solver name: " << m_flowSolverName );

  //TJ: add the surface generator to the hydrofracture solver
  m_surfaceGeneratorSolver = this->getParent()->GetGroup< SurfaceGenerator >( m_surfaceGeneratorName );
  GEOSX_ERROR_IF( m_surfaceGeneratorSolver == nullptr, this->getName() << ": invalid surface generator solver name: " << m_surfaceGeneratorName );
}

void HydrofractureSolver::InitializePostInitialConditions_PreSubGroups( Group * const GEOSX_UNUSED_PARAM( problemManager ) )
{}

HydrofractureSolver::~HydrofractureSolver()
{
  // TODO Auto-generated destructor stub
}

void HydrofractureSolver::ResetStateToBeginningOfStep( DomainPartition & domain )
{
  m_flowSolver->ResetStateToBeginningOfStep( domain );
  m_solidSolver->ResetStateToBeginningOfStep( domain );
}

real64 HydrofractureSolver::SolverStep( real64 const & time_n,
                                        real64 const & dt,
                                        int const cycleNumber,
                                        DomainPartition & domain )
{
  real64 dtReturn = dt;

  //SolverBase * const surfaceGenerator = this->getParent()->GetGroup< SolverBase >( "SurfaceGen" );
  SurfaceGenerator * const surfaceGenerator = m_surfaceGeneratorSolver;
  GEOSX_ERROR_IF( surfaceGenerator == nullptr, this->getName() << ": invalid surface generator solver name: " << std::endl );

  if( m_couplingTypeOption == CouplingTypeOption::SIM_FixedStress )
  {
    dtReturn = SplitOperatorStep( time_n, dt, cycleNumber, domain );
  }
  else if( m_couplingTypeOption == CouplingTypeOption::FIM )
  {

    ImplicitStepSetup( time_n, dt, domain );

    int const maxIter = m_maxNumResolves + 1;
    m_numResolves[1] = m_numResolves[0];
    int solveIter;
    for( solveIter=0; solveIter<maxIter; ++solveIter )
    {
      int locallyFractured = 0;
      int globallyFractured = 0;

      SetupSystem( domain,
                   m_dofManager,
                   m_localMatrix,
                   m_localRhs,
                   m_localSolution );

      if( solveIter > 0 )
      {
        m_solidSolver->ResetStressToBeginningOfStep( domain );
      }

      // currently the only method is implicit time integration
      dtReturn = NonlinearImplicitStep( time_n, dt, cycleNumber, domain );


//      m_solidSolver->updateStress( domain );

      if( surfaceGenerator!=nullptr )
      {
        if( surfaceGenerator->SolverStep( time_n, dt, cycleNumber, domain ) > 0 )
        {
          locallyFractured = 1;
        }
        MpiWrapper::allReduce( &locallyFractured,
                               &globallyFractured,
                               1,
                               MPI_MAX,
                               MPI_COMM_GEOSX );
      }
      if( globallyFractured == 0 )
      {
        break;
      }
      else
      {
        std::map< string, string_array > fieldNames;
        fieldNames["node"].emplace_back( keys::IncrementalDisplacement );
        fieldNames["node"].emplace_back( keys::TotalDisplacement );
        fieldNames["elems"].emplace_back( string( FlowSolverBase::viewKeyStruct::pressureString ) );
        fieldNames["elems"].emplace_back( "elementAperture" );

        CommunicationTools::SynchronizeFields( fieldNames,
                                               domain.getMeshBody( 0 )->getMeshLevel( 0 ),
                                               domain.getNeighbors() );

        this->UpdateDeformationForCoupling( domain );

        if( getLogLevel() >= 1 )
        {
          GEOSX_LOG_RANK_0( "++ Fracture propagation. Re-entering Newton Solve." );
        }
        m_flowSolver->ResetViews( *domain.getMeshBody( 0 )->getMeshLevel( 0 ) );
      }
    }

    // final step for completion of timestep. typically secondary variable updates and cleanup.
    ImplicitStepComplete( time_n, dtReturn, domain );
    m_numResolves[1] = solveIter;
  }

  SortedArray<localIndex> const & myTrailingFaces = m_surfaceGeneratorSolver->getTrailingFaces();
  for (auto const & trailingFace : myTrailingFaces )
  {
    std::cout << "Trailing face: " << trailingFace << std::endl;
  }

  MeshLevel * const meshLevel = domain.getMeshBody(0)->getMeshLevel(0);
  NodeManager * const nodeManager = meshLevel->getNodeManager();
  FaceManager * const faceManager = meshLevel->getFaceManager();
  FaceManager::NodeMapType const & faceToNodes = faceManager->nodeList();
  array2d<real64> const & faceNormal = faceManager->faceNormal();
  ElementRegionManager * const elementRegionManager = meshLevel->getElemManager();

  array1d<localIndex> const & childNodeIndices = nodeManager->getExtrinsicData< extrinsicMeshData::ChildIndex >();
  for (localIndex i=0; i<nodeManager->size(); i++)
  {
    if (childNodeIndices(i) >= 0)
    {
      std::cout << "Parent node: " << i << " -> child node: " << childNodeIndices(i) << std::endl;
    }
  }
  array1d<localIndex> const & childFaceIndices = faceManager->getExtrinsicData< extrinsicMeshData::ChildIndex >();
  for (localIndex i=0; i<faceManager->size(); i++)
  {
    if (childFaceIndices(i) >= 0)
    {
      std::cout << "Parent face: " << i << " -> child face: " << childFaceIndices(i) << std::endl;
    }
  }

  array1d<real64> & signedNodeDistance = nodeManager->getExtrinsicData< extrinsicMeshData::SignedNodeDistance >();
  std::cout << signedNodeDistance.size() << std::endl;
  array2d<real64, nodes::TOTAL_DISPLACEMENT_PERM> const & totalDisplacement = nodeManager->totalDisplacement();
  std::cout << totalDisplacement.size() << std::endl;
  array2d<real64, nodes::REFERENCE_POSITION_PERM> const & referencePosition = nodeManager->referencePosition();
  std::cout << referencePosition.size() << std::endl;

  elementRegionManager->forElementSubRegions<FaceElementSubRegion>( [&] (FaceElementSubRegion & subRegion)
  {
    // Get the shear modulus and bulk modulus from the constitutive model of
    // the solid material. These two material properties are needed for the
    // tip asymptotic relation
    Group const * const constitutiveModels = subRegion.GetConstitutiveModels();
    real64 shearModulus;
    real64 bulkModulus;
    constitutiveModels->forSubGroups<LinearElasticIsotropic>( [&shearModulus, &bulkModulus]
							      (LinearElasticIsotropic const & solidModel)
    {
      std::cout << solidModel.getName() << std::endl;
      shearModulus = solidModel.getReference<real64>
                                (LinearElasticIsotropic::viewKeyStruct::defaultShearModulusString);
      bulkModulus = solidModel.getReference<real64>
                                (LinearElasticIsotropic::viewKeyStruct::defaultBulkModulusString);
    } );
    std::cout << "Shear modulus = " << shearModulus << std::endl;
    std::cout << "Bulk modulus = " << bulkModulus << std::endl;

    real64 const toughness = m_surfaceGeneratorSolver->getReference<real64>("rockToughness");
    std::cout << "Rock toughness = " << toughness << std::endl;

    real64 const nu = ( 1.5 * bulkModulus - shearModulus ) / ( 3.0 * bulkModulus + shearModulus );
    real64 const E = ( 9.0 * bulkModulus * shearModulus )/ ( 3.0 * bulkModulus + shearModulus );
    real64 const Eprime = E/(1.0-nu*nu);
    real64 const PI = 2 * acos(0.0);
    real64 const Kprime = 4.0*sqrt(2.0/PI)*toughness;

    std::cout << subRegion.getName() << ": " << subRegion.size() << std::endl;
    FaceElementSubRegion::FaceMapType const & faceElmtToFacesRelation = subRegion.faceList();
    // We need to modify the face element partially open status, so no const
    array1d<integer> & isFaceElmtPartiallyOpen = subRegion.getExtrinsicData< extrinsicMeshData::IsFaceElmtPartiallyOpen >();

    // Define two unordered set, one for the fully opened face elements
    // the other for the partially opened face elements
    std::unordered_set<localIndex> fullyOpenFaceElmts, partiallyOpenFaceElmts;

    // Loop over all the faceElements
    for (localIndex iFaceElmt=0; iFaceElmt<subRegion.size(); iFaceElmt++)
    {
      // Find the parent face of the two faces in a faceElement
      localIndex const kf0 = faceElmtToFacesRelation[iFaceElmt][0];
      localIndex const kf1 = faceElmtToFacesRelation[iFaceElmt][1];
      localIndex const childFaceOfFace0 = childFaceIndices[kf0];
      localIndex const parentFace = childFaceOfFace0 >= 0 ? kf0 : kf1;

      // Count how many nodes on the parent face have child nodes
      int nodeCount = 0;
      for (localIndex const node : faceToNodes[parentFace])
      {
	if (childNodeIndices[node] >= 0 )
	  nodeCount++;
      }
      // If all the nodes on the parent face have children,
      // the faceElement is fully open. Otherwise, the faceElement
      // is partially open
      GEOSX_ERROR_IF(faceToNodes[parentFace].size() != 4,
		     "Face element " << parentFace << " needs to have FOUR nodes, "
		     << "it only has " << faceToNodes[parentFace].size() << " nodes.");

      if (nodeCount == faceToNodes[parentFace].size())
      {
	isFaceElmtPartiallyOpen[iFaceElmt] = 0;
	fullyOpenFaceElmts.insert(iFaceElmt);
      }
      else
      {
	isFaceElmtPartiallyOpen[iFaceElmt] = 1;
	partiallyOpenFaceElmts.insert(iFaceElmt);
      }

      std::cout << "Face elmt " << iFaceElmt << ": Faces " << faceElmtToFacesRelation[iFaceElmt][0] << ", "
	                                           << faceElmtToFacesRelation[iFaceElmt][1] << std::endl;
    } // for iFaceElmt < subRegion.size()
    std::cout << "Fully opened face element size: " << fullyOpenFaceElmts.size() << std::endl;
    std::cout << "Partially opened face element size: " << partiallyOpenFaceElmts.size() << std::endl;

    FaceElementSubRegion::NodeMapType const & faceElmtToNodesRelation = subRegion.nodeList();
    for (localIndex iFaceElmt=0; iFaceElmt<subRegion.size(); iFaceElmt++)
    {
      localIndex const kf0 = faceElmtToFacesRelation[iFaceElmt][0];
      localIndex const kf1 = faceElmtToFacesRelation[iFaceElmt][1];
      localIndex const childFaceOfFace0 = childFaceIndices[kf0];
      localIndex const parentFace = childFaceOfFace0 >= 0 ? kf0 : kf1;

      for (localIndex const parentNode : faceToNodes[parentFace])
      {
	localIndex const childNode = childNodeIndices[parentNode];
	if (childNode >= 0)
	{
	  real64 temp[3] = {0.0, 0.0, 0.0};
	  LvArray::tensorOps::add< 3 >( temp, totalDisplacement[parentNode] );
	  LvArray::tensorOps::subtract< 3 >( temp, totalDisplacement[childNode] );

	  // assume a linear relationship between node opening and signed distance
	  // negative distance for nodes that have been split
          if (m_regimeTypeOption == RegimeTypeOption::ToughnessDominated)
          {
            signedNodeDistance[parentNode] = -pow(-LvArray::tensorOps::AiBi< 3 >(temp, faceNormal[parentFace])
                                                  *Eprime/Kprime, 2);
          }
          else if (m_regimeTypeOption == RegimeTypeOption::ViscosityDominated)
          {
            GEOSX_ERROR("Viscosity dominated case not implemented yet.");
          }

	  signedNodeDistance[childNode] = signedNodeDistance[parentNode];
	  std::cout << "Parent node " << parentNode
		    << ": coords = " << referencePosition(parentNode, 0) << ", "
				     << referencePosition(parentNode, 1) << ", "
				     << referencePosition(parentNode, 2) << std::endl;
	  std::cout << "Parent node " << parentNode
		    << ": disps = "  << totalDisplacement(parentNode, 0) << ", "
				     << totalDisplacement(parentNode, 1) << ", "
				     << totalDisplacement(parentNode, 2) << std::endl;
	  std::cout << "Parent node " << parentNode
		    << ": signed dist = " << signedNodeDistance[parentNode] << std::endl;

	  std::cout << "Child node " << childNode
		    << ": coords = " << referencePosition(childNode, 0) << ", "
				     << referencePosition(childNode, 1) << ", "
				     << referencePosition(childNode, 2) << std::endl;
	  std::cout << "Child node " << childNode
		    << ": disps = "  << totalDisplacement(childNode, 0) << ", "
				     << totalDisplacement(childNode, 1) << ", "
				     << totalDisplacement(childNode, 2) << std::endl;
	  std::cout << "Child node " << childNode
		    << ": signed dist = " << signedNodeDistance[childNode] << std::endl;
	}
      }
    }

    for (localIndex i=0; i<subRegion.size(); i++)
    {
      std::cout << "Face elmt " << i << ": Nodes ";
      for (localIndex j : faceElmtToNodesRelation[i])
      {
	std::cout << j << ", ";
      }
      std::cout << std::endl;
    }

    // call the Eikonal Equation solver to solve the signed distance fields
    // at nodes of partially opened face elements
    EikonalEquationSolver(domain, subRegion, partiallyOpenFaceElmts);
    // calculate the geometric quantities of partially opened face elements,
    // include fluid volume, partially opened area and its center
    CalculatePartiallyOpenElmtQuantities(domain, subRegion, partiallyOpenFaceElmts, Eprime, Kprime);

  } );

  return dtReturn;
}

real64 HydrofractureSolver::CalculateSignedDistance1stOrder(localIndex const node,
                                                            std::array<std::unordered_set<localIndex>, 2> const & neighbors,
                                                            std::unordered_map<localIndex, int> & nodeStatus,
                                                            NodeManager * const nodeManager)
{
  array2d<real64, nodes::REFERENCE_POSITION_PERM> const & referencePosition = nodeManager->referencePosition();
  array1d<real64> & signedNodeDistance = nodeManager->getExtrinsicData< extrinsicMeshData::SignedNodeDistance >();

  // only nodes that are frozen (-1) will contribute to the quadratic equation
  std::array<std::unordered_set<localIndex>, 2> tempNeighbors;
  for (int i=0; i<2; i++)
  {
    for (localIndex const neighborNode : neighbors[i])
    {
      if (nodeStatus[neighborNode] == -1)
      {
        tempNeighbors[i].insert(neighborNode);
      }
    }
  }

  real64 a, b;
  localIndex nodeInDirection0, nodeInDirection1;
  // direction-0 has no neighbor nodes
  if (tempNeighbors[0].empty())
  {
    // direction-1 only has one node
    if (tempNeighbors[1].size() == 1 )
    {
      std::unordered_set<localIndex>::const_iterator itr1 = tempNeighbors[1].begin();
      a = signedNodeDistance[*itr1];
      nodeInDirection1 = *itr1;
    }
    else
    {
      std::unordered_set<localIndex>::const_iterator itr1 = tempNeighbors[1].begin();
      std::unordered_set<localIndex>::const_iterator itr2 = ++tempNeighbors[1].begin();
      a = std::min( signedNodeDistance[*itr1],
                    signedNodeDistance[*itr2] );
      nodeInDirection1 = signedNodeDistance[*itr1] < signedNodeDistance[*itr2] ? *itr1 : *itr2;
    }
    real64 temp[3] = {0.0, 0.0, 0.0};
    LvArray::tensorOps::add< 3 >( temp, referencePosition[node] );
    LvArray::tensorOps::subtract< 3 >( temp, referencePosition[nodeInDirection1] );
    return a + LvArray::tensorOps::l2Norm< 3 >( temp );
  }
  // direction-1 has no neighbor nodes
  else if (tempNeighbors[1].empty())
  {
    // direction-1 only has one node
    if (tempNeighbors[0].size() == 1 )
    {
      std::unordered_set<localIndex>::const_iterator itr1 = tempNeighbors[0].begin();
      a = signedNodeDistance[*itr1];
      nodeInDirection0 = *itr1;
    }
    else
    {
      std::unordered_set<localIndex>::const_iterator itr1 = tempNeighbors[0].begin();
      std::unordered_set<localIndex>::const_iterator itr2 = ++tempNeighbors[0].begin();
      a = std::min( signedNodeDistance[*itr1],
                    signedNodeDistance[*itr2] );
      nodeInDirection0 = signedNodeDistance[*itr1] < signedNodeDistance[*itr2] ? *itr1 : *itr2;
    }
    real64 temp[3] = {0.0, 0.0, 0.0};
    LvArray::tensorOps::add< 3 >( temp, referencePosition[node] );
    LvArray::tensorOps::subtract< 3 >( temp, referencePosition[nodeInDirection0] );
    return a + LvArray::tensorOps::l2Norm< 3 >( temp );
  }
  // both directions have neighbor nodes
  else
  {
    if (tempNeighbors[0].size() == 1)
    {
      std::unordered_set<localIndex>::const_iterator itr1 = tempNeighbors[0].begin();
      a = signedNodeDistance[*itr1];
      nodeInDirection0 = *itr1;
    }
    else
    {
      std::unordered_set<localIndex>::const_iterator itr1 = tempNeighbors[0].begin();
      std::unordered_set<localIndex>::const_iterator itr2 = ++tempNeighbors[0].begin();
      a = std::min( signedNodeDistance[*itr1],
                    signedNodeDistance[*itr2] );
      nodeInDirection0 = signedNodeDistance[*itr1] < signedNodeDistance[*itr2] ? *itr1 : *itr2;
    }
    real64 temp0[3] = {0.0, 0.0, 0.0};
    LvArray::tensorOps::add< 3 >( temp0, referencePosition[node] );
    LvArray::tensorOps::subtract< 3 >( temp0, referencePosition[nodeInDirection0] );
    real64 delta0 = LvArray::tensorOps::l2Norm< 3 >( temp0 );

    if (tempNeighbors[1].size() == 1)
    {
      std::unordered_set<localIndex>::const_iterator itr1 = tempNeighbors[1].begin();
      b = signedNodeDistance[*itr1];
      nodeInDirection1 = *itr1;
    }
    else
    {
      std::unordered_set<localIndex>::const_iterator itr1 = tempNeighbors[1].begin();
      std::unordered_set<localIndex>::const_iterator itr2 = ++tempNeighbors[1].begin();
      b = std::min( signedNodeDistance[*itr1],
                    signedNodeDistance[*itr2] );
      nodeInDirection1 = signedNodeDistance[*itr1] < signedNodeDistance[*itr2] ? *itr1 : *itr2;
    }
    real64 temp1[3] = {0.0, 0.0, 0.0};
    LvArray::tensorOps::add< 3 >( temp1, referencePosition[node] );
    LvArray::tensorOps::subtract< 3 >( temp1, referencePosition[nodeInDirection1] );
    real64 delta1 = LvArray::tensorOps::l2Norm< 3 >( temp1 );

    real64 beta = delta0/delta1;
    real64 determinant = delta0*delta0 * (1.0 + beta*beta) - beta*beta*(a-b)*(a-b);
    if (determinant >= 0.0)
    {
      return (a + beta*beta*b + std::sqrt(determinant))/(1.0+beta*beta);
    }
    else
    {
      if (a <= b)
      {
        return a + delta0;
      }
      else
      {
        return b + delta1;
      }
    }
  }
}

real64 CalculatePartiallyOpenElmtArea(real64 const cornerNodeSignedDistance,
                                      real64 const alpha,
                                      real64 const deltaX,
                                      real64 const deltaY)
{
  real64 const PI = 2 * acos(0.0);
  real64 const tol = 1.0e-9;

  real64 const maxLength = deltaX*cos(alpha) + deltaY*sin(alpha);
  real64 dist = std::abs(cornerNodeSignedDistance);
  real64 xi0;
  real64 m, area;

  if ( (alpha > tol) && (alpha < (PI/2.0 - tol)) )
  {
    if (tan(alpha) <= deltaX/deltaY)
    {
      xi0 = deltaY*sin(alpha);
    }
    else
    {
      xi0 = deltaX*cos(alpha);
    }

    m = 1.0/(sin(alpha)*cos(alpha));

    if ( dist <= xi0 ) // triangle
    {
      area = m * dist*dist/2.0;
    }
    else if ( (dist > xi0) && (dist <= (maxLength - xi0)) ) // quadrilateral
    {
      area = m * (dist*dist - (dist-xi0)*(dist-xi0))/2.0;
    }
    else if ( (dist > (maxLength - xi0)) && (dist <= maxLength) ) // pentagon
    {
      area = m * (   dist*dist
                  - (dist-xi0)*(dist-xi0)
                  - (dist+xi0-maxLength)*(dist+xi0-maxLength)
                 )/2.0;
    }
    else
    {
      area = m * (   dist*dist
                  - (dist-xi0)*(dist-xi0)
                  - (dist+xi0-maxLength)*(dist+xi0-maxLength)
                  + (dist-maxLength)*(dist-maxLength)
                 )/2.0;
    }
  }
  else
  {
    if (std::abs(alpha) <= tol)
    {
      if (dist <= deltaX)
      {
        area = deltaY*dist;
      }
      else
      {
        area = deltaY*dist - deltaY*(dist-deltaX);
      }
    }
    else if (std::abs(alpha - PI/2.0) <= tol)
    {
      if (dist <= deltaY)
      {
        area = deltaX*dist;
      }
      else
      {
        area = deltaX*dist - deltaX*(dist-deltaY);
      }
    }
    else
    {
      area = -1.0;
    }
  }

  GEOSX_ERROR_IF_LE(area, 0.0);
  GEOSX_ERROR_IF_GT(area, 1.000001*deltaX*deltaY);

  return area;
}

real64 CalculatePartiallyOpenElmtFuildVolToughnessDominated(real64 const cornerNodeSignedDistance,
                                                            real64 const alpha,
                                                            real64 const deltaX,
                                                            real64 const deltaY,
                                                            real64 const Eprime,
                                                            real64 const Kprime)
{
  real64 const PI = 2 * acos(0.0);
  real64 const tol = 1.0e-9;

  real64 const maxLength = deltaX*cos(alpha) + deltaY*sin(alpha);
  real64 const dist = std::abs(cornerNodeSignedDistance);

  real64 const coeff1 = 4.0/15.0 * Kprime/Eprime;
  real64 const coeff2 = 2.0/3.0  * Kprime/Eprime;

  real64 xi0;
  real64 m, vol;

  if ( (alpha > tol) && (alpha < (PI/2.0 - tol)) )
  {
    if (tan(alpha) <= deltaX/deltaY)
    {
      xi0 = deltaY*sin(alpha);
    }
    else
    {
      xi0 = deltaX*cos(alpha);
    }

    m = 1.0/(sin(alpha)*cos(alpha));

    if ( dist <= xi0 ) // triangle
    {
      vol = m * coeff1 * pow(dist, 5.0/2.0);
    }
    else if ( (dist > xi0) && (dist <= (maxLength - xi0)) ) // quadrilateral
    {
      vol = m * coeff1 * ( pow(dist, 5.0/2.0) - pow(dist-xi0, 5.0/2.0) );
    }
    else if ( (dist > (maxLength - xi0)) && (dist <= maxLength) ) // pentagon
    {
      vol = m * coeff1 * ( pow(dist, 5.0/2.0)
                         - pow(dist-xi0, 5.0/2.0)
                         - pow(dist+xi0-maxLength, 5.0/2.0)
                         );
    }
    else
    {
      vol = m * coeff1 * ( pow(dist, 5.0/2.0)
                         - pow(dist-xi0, 5.0/2.0)
                         - pow(dist+xi0-maxLength, 5.0/2.0)
                         + pow(dist-maxLength, 5.0/2.0)
                         );
    }
  }
  else
  {
    if (std::abs(alpha) <= tol)
    {
      if (dist <= deltaX)
      {
        vol = deltaY * coeff2 * pow(dist, 3.0/2.0);
      }
      else
      {
        vol = deltaY * coeff2 * ( pow(dist, 3.0/2.0) - pow(dist-deltaX, 3.0/2.0) );
      }
    }
    else if (std::abs(alpha - PI/2.0) <= tol)
    {
      if (dist <= deltaY)
      {
        vol = deltaX * coeff2 * pow(dist, 3.0/2.0);
      }
      else
      {
        vol = deltaX * coeff2 * ( pow(dist, 3.0/2.0) - pow(dist-deltaY, 3.0/2.0) );
      }
    }
    else
    {
      vol = -1.0;
    }
  }

  GEOSX_ERROR_IF_LE(vol, 0.0);

  return vol;
}

void HydrofractureSolver::CalculatePartiallyOpenElmtQuantities(DomainPartition & domain,
                                                               FaceElementSubRegion & subRegion,
                                                               std::unordered_set<localIndex> const & partiallyOpenFaceElmts,
                                                               real64 const Eprime,
                                                               real64 const Kprime)
{
  MeshLevel * const meshLevel = domain.getMeshBody(0)->getMeshLevel(0);
  NodeManager * const nodeManager = meshLevel->getNodeManager();
  FaceManager * const faceManager = meshLevel->getFaceManager();
  EdgeManager const * const edgeManager = meshLevel->getEdgeManager();
  FaceManager::NodeMapType const & faceToNodes = faceManager->nodeList();
  FaceManager::EdgeMapType const & faceToEdges = faceManager->edgeList();
  EdgeManager::NodeMapType const & edgeToNodes = edgeManager->nodeList();
  array1d<real64> & signedNodeDistance = nodeManager->getExtrinsicData< extrinsicMeshData::SignedNodeDistance >();
  array1d<localIndex> const & childFaceIndices = faceManager->getExtrinsicData< extrinsicMeshData::ChildIndex >();
  array1d<real64> & partiallyOpenFaceElmtFluidVol = subRegion.getExtrinsicData< extrinsicMeshData::PartiallyOpenFaceElmtFluidVol >();
  array1d<real64> & partiallyOpenFaceElmtArea = subRegion.getExtrinsicData< extrinsicMeshData::PartiallyOpenFaceElmtArea >();
  FaceElementSubRegion::FaceMapType const & faceElmtToFacesRelation = subRegion.faceList();

  array2d<real64, nodes::REFERENCE_POSITION_PERM> const & referencePosition = nodeManager->referencePosition();

  real64 const PI = 2 * acos(0.0);
  real64 const tol = 1.0e-9;

  // initialize the fluid volume in all the face element
  for (int i=0; i<subRegion.size(); i++)
  {
    partiallyOpenFaceElmtFluidVol[i] = 0.0;
    partiallyOpenFaceElmtArea[i] = 0.0;
  }

  // an unordered map to store the parent face of each face element
  // key: face element number
  // value: the face number that is the parent face
  std::unordered_map<localIndex, localIndex> parentFaceOfFaceElmt;

  // loop over all the partially open face elements to calculate
  // the fluid volume in these elements
  for (localIndex const faceElmt : partiallyOpenFaceElmts)
  {
    // the angle that determinates the fluid volume in the partially opened elmt
    real64 alpha;
    localIndex const kf0 = faceElmtToFacesRelation[faceElmt][0];
    localIndex const kf1 = faceElmtToFacesRelation[faceElmt][1];
    localIndex const childFaceOfFace0 = childFaceIndices[kf0];
    localIndex const parentFace = childFaceOfFace0 >= 0 ? kf0 : kf1;
    parentFaceOfFaceElmt[faceElmt] = parentFace;

    localIndex cornerNodeIndex = 0;
    real64 cornerNodeSignedDistance = 0.0;
    for (localIndex const parentNode : faceToNodes[parentFace] )
    {
      if (signedNodeDistance[parentNode] < cornerNodeSignedDistance)
      {
        cornerNodeSignedDistance = signedNodeDistance[parentNode];
        cornerNodeIndex = parentNode;
      }
    }
    std::cout << "corner node " << cornerNodeIndex << " => " << cornerNodeSignedDistance << std::endl;

    std::vector<localIndex> cornerNodeNeighbors;
    for (localIndex const edge : faceToEdges[parentFace])
    {
      if (edgeManager->hasNode(edge, cornerNodeIndex))
      {
        localIndex const node0 = edgeToNodes[edge][0];
        localIndex const node1 = edgeToNodes[edge][1];
        localIndex neighborNode = (cornerNodeIndex == node0) ? node1 : node0;
        cornerNodeNeighbors.push_back(neighborNode);
      }
    }
    std::cout << "corner node neighbor size = " << cornerNodeNeighbors.size()
              << " : " << cornerNodeNeighbors[0]
              << " , " << cornerNodeNeighbors[1] << std::endl;

    real64 const dist0 = signedNodeDistance[cornerNodeNeighbors[0]];
    real64 const dist1 = signedNodeDistance[cornerNodeNeighbors[1]];

    real64 temp0[3] = {0.0, 0.0, 0.0};
    LvArray::tensorOps::add< 3 >( temp0, referencePosition[cornerNodeIndex] );
    LvArray::tensorOps::subtract< 3 >( temp0, referencePosition[cornerNodeNeighbors[0]] );
    real64 delta0 = LvArray::tensorOps::l2Norm< 3 >( temp0 );

    real64 temp1[3] = {0.0, 0.0, 0.0};
    LvArray::tensorOps::add< 3 >( temp1, referencePosition[cornerNodeIndex] );
    LvArray::tensorOps::subtract< 3 >( temp1, referencePosition[cornerNodeNeighbors[1]] );
    real64 delta1 = LvArray::tensorOps::l2Norm< 3 >( temp1 );

    real64 const theta = atan(delta1/delta0);
    real64 const diag = sqrt( delta0*delta0 + delta1*delta1 );

    if (std::abs(cornerNodeSignedDistance - dist0) < tol)
    {
      alpha = PI/2.0;
    }
    else if (std::abs(cornerNodeSignedDistance - dist1) < tol)
    {
      alpha = 0.0;
    }
    else
    {
      alpha = PI/2.0 - theta - asin( (dist0 - dist1)/diag );
    }
    std::cout << alpha << std::endl;

    partiallyOpenFaceElmtArea[faceElmt] = CalculatePartiallyOpenElmtArea(cornerNodeSignedDistance,
                                                                         alpha,
                                                                         delta0,
                                                                         delta1);
    if (m_regimeTypeOption == RegimeTypeOption::ToughnessDominated)
    {
      partiallyOpenFaceElmtFluidVol[faceElmt] = CalculatePartiallyOpenElmtFuildVolToughnessDominated(cornerNodeSignedDistance,
                                                                                                     alpha,
                                                                                                     delta0,
                                                                                                     delta1,
                                                                                                     Eprime,
                                                                                                     Kprime);
    }
    else if (m_regimeTypeOption == RegimeTypeOption::ViscosityDominated)
    {
      GEOSX_ERROR("Viscosity-dominated case has not been implemented yet.");
    }

    std::cout << "Face elmt " << faceElmt << " partial fluid vol = "
              << partiallyOpenFaceElmtFluidVol[faceElmt] << ", "
              << partiallyOpenFaceElmtArea[faceElmt]<< std::endl;
  }

}

void HydrofractureSolver::EikonalEquationSolver(DomainPartition & domain,
						FaceElementSubRegion & subRegion,
						std::unordered_set<localIndex> const & partiallyOpenFaceElmts)
{
  MeshLevel * const meshLevel = domain.getMeshBody(0)->getMeshLevel(0);
  NodeManager * const nodeManager = meshLevel->getNodeManager();
  array2d<real64, nodes::REFERENCE_POSITION_PERM> const & referencePosition = nodeManager->referencePosition();
  array1d<real64> & signedNodeDistance = nodeManager->getExtrinsicData< extrinsicMeshData::SignedNodeDistance >();
  array1d<localIndex> const & childNodeIndices = nodeManager->getExtrinsicData< extrinsicMeshData::ChildIndex >();
  FaceManager * const faceManager = meshLevel->getFaceManager();
  FaceManager::NodeMapType const & faceToNodes = faceManager->nodeList();
  FaceManager::EdgeMapType const & faceToEdges = faceManager->edgeList();
  EdgeManager const * const edgeManager = meshLevel->getEdgeManager();
  EdgeManager::NodeMapType const & edgeToNodes = edgeManager->nodeList();
  array1d<localIndex> const & childFaceIndices = faceManager->getExtrinsicData< extrinsicMeshData::ChildIndex >();
  FaceElementSubRegion::FaceMapType const & faceElmtToFacesRelation = subRegion.faceList();
  FaceElementSubRegion::NodeMapType const & faceElmtToNodesRelation = subRegion.nodeList();

  // a set to store all the nodes located on the parent face of the partially opened
  // face element, some of these nodes are frozen with known signed distance as B.C. for
  // the Eikonal equation, other nodes are far
  std::unordered_set<localIndex> allNodeSet;
  std::unordered_set<localIndex> boundaryNodeSet;

  // an unordered map to store the parent face of each face element
  // key: face element number
  // value: the face number that is the parent face
  std::unordered_map<localIndex, localIndex> parentFaceOfFaceElmt;

  for (localIndex const faceElmt : partiallyOpenFaceElmts)
  {
    localIndex const kf0 = faceElmtToFacesRelation[faceElmt][0];
    localIndex const kf1 = faceElmtToFacesRelation[faceElmt][1];
    localIndex const childFaceOfFace0 = childFaceIndices[kf0];
    localIndex const parentFace = childFaceOfFace0 >= 0 ? kf0 : kf1;
    parentFaceOfFaceElmt[faceElmt] = parentFace;

    for (localIndex const parentNode : faceToNodes[parentFace] )
    {
      allNodeSet.insert(parentNode);
    }
    std::cout << "Face Elmt " << faceElmt << " :" << std::endl;
    for (localIndex const node : faceElmtToNodesRelation[faceElmt])
    {
      std::cout << "Node " << node << " : " << signedNodeDistance[node] << std::endl;
    }
  } // for faceElmt : partiallyOpenFaceElmts

  // an unordered map to store the node status:
  // key: node number;
  // value: -1->node is frozen; 0->node is in narrow band; 1->node is far.
  std::unordered_map<localIndex, int> nodeStatus;

  for (localIndex const & node : allNodeSet)
  {
    if (childNodeIndices(node) >=0)
    {
      // if node has child, it is frozen (-1).
      nodeStatus[node] = -1;
      boundaryNodeSet.insert(node);
    }
    else
    {
      // if node has no child, it is far (1).
      nodeStatus[node] = 1;
    }
  }

  // an unordered map to store the neighbors of each node
  // key: node number;
  // value: an array with two elements, each of which is an unordered set
  //        storing the neighbor elements in one direction. Since it is a
  //        two-dimensional problem, each unordered set stores neighbor nodes
  //        in one of the two directions.
  std::unordered_map< localIndex, std::array<std::unordered_set<localIndex>, 2> > nodeNeighbors;

  // find neighbor nodes for each node in the allNodeSet
  for (localIndex const faceElmt : partiallyOpenFaceElmts)
  {
    localIndex const parentFace = parentFaceOfFaceElmt[faceElmt];
    std::cout << "In face element " << faceElmt << std::endl;
    for (localIndex const node : faceToNodes[parentFace])
    {
      for (localIndex const edge : faceToEdges[parentFace])
      {
        if (edgeManager->hasNode(edge, node))
        {
          std::cout << "Node " << node << " is on Edge " << edge << std::endl;
          localIndex const node0 = edgeToNodes[edge][0];
          localIndex const node1 = edgeToNodes[edge][1];
          localIndex neighborNode = (node == node0) ? node1 : node0;
          std::cout << "Node " << node << " has neighbor node " << neighborNode << std::endl;

          if (nodeNeighbors[node][0].empty() && nodeNeighbors[node][1].empty())
          {
            nodeNeighbors[node][0].insert(neighborNode);
          }
          else if (nodeNeighbors[node][0].empty())
          {
            localIndex frontNode = *(nodeNeighbors[node][1].begin());
            real64 temp[3] = {0.0, 0.0, 0.0};
            LvArray::tensorOps::add< 3 >( temp, referencePosition[frontNode] );
            LvArray::tensorOps::subtract< 3 >( temp, referencePosition[node] );
            real64 temp1[3] = {0.0, 0.0, 0.0};
            LvArray::tensorOps::add< 3 >( temp1, referencePosition[neighborNode] );
            LvArray::tensorOps::subtract< 3 >( temp1, referencePosition[node] );
            real64 product = LvArray::tensorOps::AiBi< 3 >( temp, temp1 );
            if (std::abs(product) > 1.0e-6)
            {
              nodeNeighbors[node][1].insert(neighborNode);
            }
            else
            {
              nodeNeighbors[node][0].insert(neighborNode);
            }
          }
          else
          {
            localIndex frontNode = *(nodeNeighbors[node][0].begin());
            real64 temp[3] = {0.0, 0.0, 0.0};
            LvArray::tensorOps::add< 3 >( temp, referencePosition[frontNode] );
            LvArray::tensorOps::subtract< 3 >( temp, referencePosition[node] );
            real64 temp1[3] = {0.0, 0.0, 0.0};
            LvArray::tensorOps::add< 3 >( temp1, referencePosition[neighborNode] );
            LvArray::tensorOps::subtract< 3 >( temp1, referencePosition[node] );
            real64 product = LvArray::tensorOps::AiBi< 3 >( temp, temp1 );
            if (std::abs(product) > 1.0e-6)
            {
              nodeNeighbors[node][0].insert(neighborNode);
            }
            else
            {
              nodeNeighbors[node][1].insert(neighborNode);
            }
          }
        }
      }
    }
  }

  std::unordered_set<localIndex> narrowBandNodeSet;
  std::vector<std::pair<localIndex, real64>> narrowBandMinHeap;
  narrowBandMinHeap.reserve(allNodeSet.size());

  // step-1: initialize the nodes that are neighboring to the boundary
  //         nodes (split nodes)
  for (localIndex const bcNode : boundaryNodeSet)
  {
    // loop over two directions
    for (int i=0; i<2; i++)
    {
      for (localIndex const neighborNode : nodeNeighbors[bcNode][i])
      {
        // if neighbor node is far (1)
        if ( nodeStatus[neighborNode] == 1 )
        {
          // change the node status into narrow band (0)
          nodeStatus[neighborNode] = 0;

          // calculate the signed distance
          real64 signedDist = CalculateSignedDistance1stOrder(neighborNode,
                                                              nodeNeighbors[neighborNode],
                                                              nodeStatus,
                                                              nodeManager);
          signedNodeDistance[neighborNode] = signedDist;
          narrowBandNodeSet.insert(neighborNode);
          narrowBandMinHeap.emplace_back( std::pair<localIndex, real64>(neighborNode, signedDist) );
        }
      }
    }
  }
  // greater comparison to make a min heap
  std::make_heap(narrowBandMinHeap.begin(),
                 narrowBandMinHeap.end(),
                 [] (std::pair<localIndex, real64> const a, std::pair<localIndex, real64> const b)
                 {
                   return a.second > b.second;
                 });
  for (auto & item : narrowBandMinHeap)
  {
    std::cout << item.first << " => " << item.second << std::endl;
  }

  //step-2: calculate the signed distance at the nodes that are not split
  while(!narrowBandNodeSet.empty())
  {
    std::cout << "narrowBandNodeSet size = " << narrowBandNodeSet.size() << std::endl;

    for (auto & item : narrowBandMinHeap)
    {
      std::cout << item.first << " => " << item.second << std::endl;
    }
    std::pop_heap(narrowBandMinHeap.begin(),
                  narrowBandMinHeap.end(),
                  [] (std::pair<localIndex, real64> const a, std::pair<localIndex, real64> const b)
                  {
                    return a.second > b.second;
                  });
    for (auto & item : narrowBandMinHeap)
    {
      std::cout << item.first << " => " << item.second << std::endl;
    }
    std::pair<localIndex, real64> const minEntry = narrowBandMinHeap.back();
    narrowBandMinHeap.pop_back();
    std::cout << "minEntry: " << minEntry.first << " => " << minEntry.second << std::endl;
    for (auto & item : narrowBandMinHeap)
    {
      std::cout << item.first << " => " << item.second << std::endl;
    }
    // remove the node that has the smallest signed distance
    // from the narrow band node set, if the node has already
    // been removed, then do nothing
    if (narrowBandNodeSet.erase(minEntry.first))
    {
      // the node status must be narrow band
      GEOSX_ERROR_IF(nodeStatus[minEntry.first] != 0,
                           "The status of node " << minEntry.first << " is wrong!");
      // freeze this node
      nodeStatus[minEntry.first] = -1;
      // assign the signed distance
      signedNodeDistance[minEntry.first] = minEntry.second;
      // loop over two directions
      for (int i=0; i<2; i++)
      {
        for (localIndex const neighborNode : nodeNeighbors[minEntry.first][i])
        {
          // if the neighbor node is far
          if (nodeStatus[neighborNode] == 1)
          {
            // change its status to narrow band
            nodeStatus[neighborNode] = 0;
            // calculate the signed distance
            real64 signedDist = CalculateSignedDistance1stOrder(neighborNode,
                                                                nodeNeighbors[neighborNode],
                                                                nodeStatus,
                                                                nodeManager);
            signedNodeDistance[neighborNode] = signedDist;
            narrowBandNodeSet.insert(neighborNode);
            narrowBandMinHeap.emplace_back( std::pair<localIndex, real64>(neighborNode,
                                                                          signedDist) );
            std::push_heap(narrowBandMinHeap.begin(),
                           narrowBandMinHeap.end(),
                           [] (std::pair<localIndex, real64> const a, std::pair<localIndex, real64> const b)
                           {
                             return a.second > b.second;
                           });
          }
          // if the neighbor node is in narrow band
          else if (nodeStatus[neighborNode] == 0)
          {
            real64 newSignedDist = CalculateSignedDistance1stOrder(neighborNode,
                                                                   nodeNeighbors[neighborNode],
                                                                   nodeStatus,
                                                                   nodeManager);
            real64 oldSignedDist = signedNodeDistance[neighborNode];
            if (newSignedDist < oldSignedDist)
            {
              signedNodeDistance[neighborNode] = newSignedDist;
              narrowBandMinHeap.emplace_back( std::pair<localIndex, real64>(neighborNode,
                                                                            newSignedDist) );
              std::push_heap(narrowBandMinHeap.begin(),
                             narrowBandMinHeap.end(),
                             [] (std::pair<localIndex, real64> const a, std::pair<localIndex, real64> const b)
                             {
                               return a.second > b.second;
                             });

            }
          }
          // if the neighbor node is frozen, do nothing
          else
          {
            ;
          }
        }
      }
    } // if (narrowBandNodeSet.erase(minEntry.first))
  } // while(!narrowBandNodeSet.empty())




  std::unordered_map<localIndex, std::pair<real64, localIndex>> dict;
  std::vector<std::pair<real64, localIndex>*> vec;
  dict[100] = std::pair<real64, localIndex>(-0.3, 100);
  vec.push_back(&dict[100]);
  dict[80] = std::pair<real64, localIndex>(-0.1, 80);
  vec.push_back(&dict[80]);
  dict[60] = std::pair<real64, localIndex>(0.0, 60);
  vec.push_back(&dict[60]);

  for (auto item : vec)
  {
    std::cout << (*item).second << " => " << (*item).first << std::endl;
  }

  dict[80].first = 0.11;
  std::cout << "Change the value of 80" << std::endl;

  for (auto item : vec)
  {
    std::cout << (*item).second << " => " << (*item).first << std::endl;
  }

  std::sort(vec.begin(),
            vec.end(),
            [ ] (std::pair<real64, localIndex> const* a, std::pair<real64, localIndex> const* b)
            {
              return (*a).first < (*b).first;
            } );

  for (auto item : vec)
  {
    std::cout << (*item).second << " => " << (*item).first << std::endl;
  }

  std::sort(vec.begin(),
            vec.end(),
            [] (std::pair<real64, localIndex> const* a, std::pair<real64, localIndex> const* b)
            {
              return (*a).first > (*b).first;
            } );

  for (auto item : vec)
  {
    std::cout << (*item).second << " => " << (*item).first << std::endl;
  }

  std::cout << "make heap" << std::endl;
  // greater comparison to make a min heap
  std::make_heap(vec.begin(),
                 vec.end(),
                 [] (std::pair<real64, localIndex> const* a, std::pair<real64, localIndex> const* b)
                 {
                   return (*a).first > (*b).first;
                 });
  for (auto item : vec)
  {
    std::cout << (*item).second << " => " << (*item).first << std::endl;
  }

  dict[80].first = -0.5;
  std::cout << "make heap after change value of 80" << std::endl;
  // greater comparison to make a min heap
  std::make_heap(vec.begin(),
                 vec.end(),
                 [] (std::pair<real64, localIndex> const* a, std::pair<real64, localIndex> const* b)
                 {
                   return (*a).first > (*b).first;
                 });
  for (auto item : vec)
  {
    std::cout << (*item).second << " => " << (*item).first << std::endl;
  }


   vec.pop_back();
   std::cout << "Pop" << std::endl;
   for (auto item : vec)
   {
     std::cout << (*item).second << " => " << (*item).first << std::endl;
   }

   std::cout << "Dict" << std::endl;
   for (auto & item : dict)
   {
     std::cout << item.first << " -> " << (item.second).first
         << ", "<< (item.second).second << std::endl;
   }


/*
  dict[100] = -0.3;
  dict[80] = -0.1;
  std::unordered_map<localIndex, real64>::iterator itr1, itr2;
  itr1 = dict.begin();
  itr2 = ++dict.begin();

  std::cout << itr1->first << " => " << itr1->second << std::endl;
  std::cout << itr2->first << " => " << itr2->second << std::endl;

  dict[80] = -0.15;

  std::cout << itr1->first << " => " << itr1->second << std::endl;
  std::cout << itr2->first << " => " << itr2->second << std::endl;

  real64 * ptr = &(itr1->second);
  std::cout << *ptr << std::endl;

  dict[80] = -0.25;
  std::cout << *ptr << std::endl;
*/






  std::cout << "Size of all node set = " << allNodeSet.size() << std::endl;
  for (auto const & item : nodeStatus)
    std::cout << item.first << " -> " << item.second << std::endl;

  for (auto const & item : nodeNeighbors)
  {
    std::cout << "Node " << item.first << ": ";
    for (localIndex const i : (item.second)[0])
    {
      std::cout << i << ", ";
    }
    std::cout << std::endl;
    std::cout << "Node " << item.first << ": ";
    for (localIndex const i : (item.second)[1])
    {
      std::cout << i << ", ";
    }
    std::cout << std::endl;

  }

}

void HydrofractureSolver::UpdateDeformationForCoupling( DomainPartition & domain )
{
  MeshLevel * const meshLevel = domain.getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = meshLevel->getElemManager();
  NodeManager const * const nodeManager = meshLevel->getNodeManager();
  FaceManager * const faceManager = meshLevel->getFaceManager();

  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const u = nodeManager->totalDisplacement();
  arrayView2d< real64 const > const faceNormal = faceManager->faceNormal();
  // arrayView1d<real64 const> const faceArea = faceManager->faceArea();
  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager->nodeList().toViewConst();

  ConstitutiveManager const * const constitutiveManager = domain.getConstitutiveManager();

  ContactRelationBase const * const
  contactRelation = constitutiveManager->GetGroup< ContactRelationBase >( m_contactRelationName );

  elemManager->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
  {
    arrayView1d< real64 > const aperture = subRegion.getElementAperture();
    arrayView1d< real64 > const effectiveAperture = subRegion.getReference< array1d< real64 > >( FlowSolverBase::viewKeyStruct::effectiveApertureString );
    arrayView1d< real64 const > const volume = subRegion.getElementVolume();
    arrayView1d< real64 > const deltaVolume = subRegion.getReference< array1d< real64 > >( FlowSolverBase::viewKeyStruct::deltaVolumeString );
    arrayView1d< real64 const > const area = subRegion.getElementArea();
    arrayView2d< localIndex const > const elemsToFaces = subRegion.faceList();

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
    arrayView1d< real64 const > const &
    apertureF = subRegion.getReference< array1d< real64 > >( viewKeyStruct::apertureAtFailureString );

    arrayView1d< real64 > const &
    separationCoeff = subRegion.getSeparationCoefficient();

    arrayView1d< real64 > const &
    dSeparationCoeff_dAper = subRegion.getReference< array1d< real64 > >( FaceElementSubRegion::viewKeyStruct::dSeparationCoeffdAperString );
    arrayView1d< real64 const > const &
    separationCoeff0 = subRegion.getReference< array1d< real64 > >( viewKeyStruct::separationCoeff0String );
#endif

    forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const kfe )
    {
      localIndex const kf0 = elemsToFaces[kfe][0];
      localIndex const kf1 = elemsToFaces[kfe][1];
      localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kf0 );
      real64 temp[ 3 ] = { 0 };
      for( localIndex a=0; a<numNodesPerFace; ++a )
      {
        LvArray::tensorOps::add< 3 >( temp, u[ faceToNodeMap( kf0, a ) ] );
        LvArray::tensorOps::subtract< 3 >( temp, u[ faceToNodeMap( kf1, a ) ] );
      }

      // TODO this needs a proper contact based strategy for aperture
      aperture[kfe] = -LvArray::tensorOps::AiBi< 3 >( temp, faceNormal[ kf0 ] ) / numNodesPerFace;

      effectiveAperture[kfe] = contactRelation->effectiveAperture( aperture[kfe] );


#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
      real64 const s = aperture[kfe] / apertureF[kfe];
      if( separationCoeff0[kfe]<1.0 && s>separationCoeff0[kfe] )
      {
        if( s >= 1.0 )
        {
          separationCoeff[kfe] = 1.0;
          dSeparationCoeff_dAper[kfe] = 0.0;
        }
        else
        {
          separationCoeff[kfe] = s;
          dSeparationCoeff_dAper[kfe] = 1.0/apertureF[kfe];
        }
      }
#endif
      deltaVolume[kfe] = effectiveAperture[kfe] * area[kfe] - volume[kfe];
    } );

//#if defined(USE_CUDA)
//    deltaVolume.move( LvArray::MemorySpace::GPU );
//    aperture.move( LvArray::MemorySpace::GPU );
//    effectiveAperture.move( LvArray::MemorySpace::GPU );
//#endif
  } );
}

real64 HydrofractureSolver::SplitOperatorStep( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                               real64 const & dt,
                                               integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                                               DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{
  GEOSX_ERROR( "Not implemented" );
  real64 dtReturn = dt;
//  real64 dtReturnTemporary = dtReturn;
//
//  m_flowSolver->ImplicitStepSetup( time_n, dt, domain, getLinearSystemRepository() );
//  m_solidSolver->ImplicitStepSetup( time_n, dt, domain, getLinearSystemRepository() );
//  this->ImplicitStepSetup( time_n, dt, domain, getLinearSystemRepository() );
//
//
//
//  fluidSolver.ImplicitStepSetup( time_n, dt, domain,
//                                 fluidSolver.getDofManager(),
//                                 fluidSolver.getSystemMatrix(),
//                                 fluidSolver.getSystemRhs(),
//                                 fluidSolver.getSystemSolution() );
//
//  solidSolver.ImplicitStepSetup( time_n, dt, domain,
//                                 solidSolver.getDofManager(),
//                                 solidSolver.getSystemMatrix(),
//                                 solidSolver.getSystemRhs(),
//                                 solidSolver.getSystemSolution() );
//
//  this->UpdateDeformationForCoupling(domain);
//
//  int iter = 0;
//  while (iter < solverParams->maxIterNewton() )
//  {
//    if (iter == 0)
//    {
//      // reset the states of all child solvers if any of them has been reset
//      m_flowSolver->ResetStateToBeginningOfStep( domain );
//      m_solidSolver->ResetStateToBeginningOfStep( domain );
//      ResetStateToBeginningOfStep( domain );
//    }
//    LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1  << ", FlowSolver: " );
//
//    // call assemble to fill the matrix and the rhs
//    m_flowSolver->AssembleSystem( domain, getLinearSystemRepository(), time_n+dt, dt );
//
//    // apply boundary conditions to system
//    m_flowSolver->ApplyBoundaryConditions( domain, getLinearSystemRepository(), time_n, dt );
//
//    // call the default linear solver on the system
//    m_flowSolver->SolveSystem( getLinearSystemRepository(),
//                 getLinearSolverParameters() );
//
//    // apply the system solution to the fields/variables
//    m_flowSolver->ApplySystemSolution( getLinearSystemRepository(), 1.0, domain );
//
//    if (dtReturnTemporary < dtReturn)
//    {
//      iter = 0;
//      dtReturn = dtReturnTemporary;
//      continue;
//    }
//
////    if (m_fluidSolver->getLinearSolverParameters()->numNewtonIterations() == 0 && iter > 0 && getLogLevel() >= 1)
////    {
////      GEOSX_LOG_RANK_0( "***** The iterative coupling has converged in " << iter  << " iterations! *****\n" );
////      break;
////    }
//
//    if (getLogLevel() >= 1)
//    {
//      GEOSX_LOG_RANK_0( "\tIteration: " << iter+1  << ", MechanicsSolver: " );
//    }
//
//    // call assemble to fill the matrix and the rhs
//    m_solidSolver->AssembleSystem( domain, getLinearSystemRepository(), time_n+dt, dt );
//
//
//    ApplyFractureFluidCoupling( domain, *getLinearSystemRepository() );
//
//    // apply boundary conditions to system
//    m_solidSolver->ApplyBoundaryConditions( domain, getLinearSystemRepository(), time_n, dt );
//
//    // call the default linear solver on the system
//    m_solidSolver->SolveSystem( getLinearSystemRepository(),
//                 getLinearSolverParameters() );
//
//    // apply the system solution to the fields/variables
//    m_solidSolver->ApplySystemSolution( getLinearSystemRepository(), 1.0, domain );
//
//    if( m_flowSolver->CalculateResidualNorm( getLinearSystemRepository(), domain ) < solverParams->newtonTol() &&
//        m_solidSolver->CalculateResidualNorm( getLinearSystemRepository(), domain ) < solverParams->newtonTol() )
//    {
//      GEOSX_LOG_RANK_0( "***** The iterative coupling has converged in " << iter  << " iterations! *****\n" );
//      break;
//    }
//
//    if (dtReturnTemporary < dtReturn)
//    {
//      iter = 0;
//      dtReturn = dtReturnTemporary;
//      continue;
//    }
////    if (m_solidSolver->getLinearSolverParameters()->numNewtonIterations() > 0)
//    {
//      this->UpdateDeformationForCoupling(domain);
////      m_fluidSolver->UpdateState(domain);
//    }
//    ++iter;
//  }
//
//  this->ImplicitStepComplete( time_n, dt, domain );

  return dtReturn;
}

real64 HydrofractureSolver::ExplicitStep( real64 const & time_n,
                                          real64 const & dt,
                                          const int cycleNumber,
                                          DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->ExplicitStep( time_n, dt, cycleNumber, domain );
  m_flowSolver->SolverStep( time_n, dt, cycleNumber, domain );

  return dt;
}


void HydrofractureSolver::SetupDofs( DomainPartition const & domain,
                                     DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;

  m_solidSolver->SetupDofs( domain, dofManager );
  m_flowSolver->SetupDofs( domain, dofManager );

  // restrict coupling to fracture regions only (as done originally in SetupSystem)
  ElementRegionManager const & elemManager = *domain.getMeshBody( 0 )->getMeshLevel( 0 )->getElemManager();
  string_array fractureRegions;
  elemManager.forElementRegions< SurfaceElementRegion >( [&]( SurfaceElementRegion const & elementRegion )
  {
    fractureRegions.emplace_back( elementRegion.getName() );
  } );

  dofManager.addCoupling( keys::TotalDisplacement,
                          FlowSolverBase::viewKeyStruct::pressureString,
                          DofManager::Connector::Elem,
                          fractureRegions );

}

void HydrofractureSolver::SetupSystem( DomainPartition & domain,
                                       DofManager & dofManager,
                                       CRSMatrix< real64, globalIndex > & localMatrix,
                                       array1d< real64 > & localRhs,
                                       array1d< real64 > & localSolution,
                                       bool const setSparsity )
{
  GEOSX_MARK_FUNCTION;

  GEOSX_UNUSED_VAR( setSparsity );

  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  m_flowSolver->ResetViews( mesh );

  dofManager.setMesh( domain, 0, 0 );

  SetupDofs( domain, dofManager );
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

  localRhs.resize( numLocalRows );
  localSolution.resize( numLocalRows );

  localMatrix.setName( this->getName() + "/localMatrix" );
  localRhs.setName( this->getName() + "/localRhs" );
  localSolution.setName( this->getName() + "/localSolution" );

  m_flowSolver->setUpDflux_dApertureMatrix( domain, dofManager, localMatrix );

}

void HydrofractureSolver::addFluxApertureCouplingNNZ( DomainPartition & domain,
                                                      DofManager & dofManager,
                                                      arrayView1d< localIndex > const & rowLengths ) const
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  ElementRegionManager const & elemManager = *mesh.getElemManager();

  string const presDofKey = dofManager.getKey( FlowSolverBase::viewKeyStruct::pressureString );

  globalIndex const rankOffset = dofManager.rankOffset();

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_flowSolver->getDiscretization() );

  fluxApprox.forStencils< FaceElementStencil >( mesh, [&]( FaceElementStencil const & stencil )
  {
    for( localIndex iconn=0; iconn<stencil.size(); ++iconn )
    {
      localIndex const numFluxElems = stencil.stencilSize( iconn );
      typename FaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
      typename FaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
      typename FaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

      FaceElementSubRegion const & elementSubRegion =
        *elemManager.GetRegion( seri[iconn][0] )->GetSubRegion< FaceElementSubRegion >( sesri[iconn][0] );

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
    }//);
  } );

}

void HydrofractureSolver::addFluxApertureCouplingSparsityPattern( DomainPartition & domain,
                                                                  DofManager & dofManager,
                                                                  SparsityPatternView< globalIndex > const & pattern ) const
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  NodeManager const & nodeManager = *mesh.getNodeManager();
  ElementRegionManager const & elemManager = *mesh.getElemManager();

  string const presDofKey = dofManager.getKey( FlowSolverBase::viewKeyStruct::pressureString );
  string const dispDofKey = dofManager.getKey( keys::TotalDisplacement );

  globalIndex const rankOffset = dofManager.rankOffset();

  arrayView1d< globalIndex const > const & dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_flowSolver->getDiscretization() );

  fluxApprox.forStencils< FaceElementStencil >( mesh, [&]( FaceElementStencil const & stencil )
  {
    for( localIndex iconn=0; iconn<stencil.size(); ++iconn )
    {
      localIndex const numFluxElems = stencil.stencilSize( iconn );
      typename FaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
      typename FaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
      typename FaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

      FaceElementSubRegion const & elementSubRegion =
        *elemManager.GetRegion( seri[iconn][0] )->GetSubRegion< FaceElementSubRegion >( sesri[iconn][0] );

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
    }  //);
  } );
}

void HydrofractureSolver::AssembleSystem( real64 const time,
                                          real64 const dt,
                                          DomainPartition & domain,
                                          DofManager const & dofManager,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  m_solidSolver->AssembleSystem( time,
                                 dt,
                                 domain,
                                 dofManager,
                                 localMatrix,
                                 localRhs );

  m_flowSolver->ResetViews( *(domain.getMeshBody( 0 )->getMeshLevel( 0 ) ) );

  m_flowSolver->AssembleSystem( time,
                                dt,
                                domain,
                                dofManager,
                                localMatrix,
                                localRhs );


  AssembleForceResidualDerivativeWrtPressure( domain, localMatrix, localRhs );

  AssembleFluidMassResidualDerivativeWrtDisplacement( domain, localMatrix );
}

void HydrofractureSolver::ApplyBoundaryConditions( real64 const time,
                                                   real64 const dt,
                                                   DomainPartition & domain,
                                                   DofManager const & dofManager,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                   arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  m_solidSolver->ApplyBoundaryConditions( time,
                                          dt,
                                          domain,
                                          dofManager,
                                          localMatrix,
                                          localRhs );

  m_flowSolver->ApplyBoundaryConditions( time,
                                         dt,
                                         domain,
                                         dofManager,
                                         localMatrix,
                                         localRhs );
}

real64
HydrofractureSolver::
  CalculateResidualNorm( DomainPartition const & domain,
                         DofManager const & dofManager,
                         arrayView1d< real64 const > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  real64 const fluidResidual = m_flowSolver->CalculateResidualNorm( domain,
                                                                    dofManager,
                                                                    localRhs );

  real64 const solidResidual = m_solidSolver->CalculateResidualNorm( domain,
                                                                     dofManager,
                                                                     localRhs );

  if( getLogLevel() >= 1 && logger::internal::rank==0 )
  {
    char output[200] = {0};
    sprintf( output, "    ( Rfluid, Rsolid ) = ( %4.2e, %4.2e )", fluidResidual, solidResidual );
    std::cout << output << std::endl;
  }

  return std::sqrt( fluidResidual * fluidResidual + solidResidual * solidResidual );
}



void
HydrofractureSolver::
  AssembleForceResidualDerivativeWrtPressure( DomainPartition & domain,
                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                              arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;
  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  FaceManager const & faceManager = *mesh.getFaceManager();
  NodeManager & nodeManager = *mesh.getNodeManager();
  ElementRegionManager const & elemManager = *mesh.getElemManager();

  arrayView2d< real64 const > const & faceNormal = faceManager.faceNormal();
  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager.nodeList().toViewConst();

  arrayView2d< real64 > const &
  fext = nodeManager.getReference< array2d< real64 > >( SolidMechanicsLagrangianFEM::viewKeyStruct::forceExternal );
  fext.setValues< serialPolicy >( 0 );

  string const presDofKey = m_dofManager.getKey( FlowSolverBase::viewKeyStruct::pressureString );
  string const dispDofKey = m_dofManager.getKey( keys::TotalDisplacement );

  globalIndex const rankOffset = m_dofManager.rankOffset();
  arrayView1d< globalIndex const > const & dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );

  elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion const & subRegion )
  {

    arrayView1d< globalIndex const > const &
    pressureDofNumber = subRegion.getReference< array1d< globalIndex > >( presDofKey );

    if( subRegion.hasWrapper( "pressure" ) )
    {
      arrayView1d< real64 const > const & fluidPressure = subRegion.getReference< array1d< real64 > >( "pressure" );
      arrayView1d< real64 const > const & deltaFluidPressure = subRegion.getReference< array1d< real64 > >( "deltaPressure" );
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
        globalIndex rowDOF[24]; // Here it assumes 8 nodes?
        real64 nodeRHS[24];  // Here it assumes 8 nodes?
        stackArray2d< real64, 12*12 > dRdP( numNodesPerFace*3, 1 );
        globalIndex colDOF = pressureDofNumber[kfe];

        real64 const Ja = area[kfe] / numNodesPerFace;

        real64 nodalForceMag = ( fluidPressure[kfe]+deltaFluidPressure[kfe] ) * Ja;
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
  AssembleFluidMassResidualDerivativeWrtDisplacement( DomainPartition const & domain,
                                                      CRSMatrixView< real64, globalIndex const > const & localMatrix )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  FaceManager const & faceManager = *mesh.getFaceManager();
  NodeManager const & nodeManager = *mesh.getNodeManager();
  ConstitutiveManager const & constitutiveManager = *domain.getConstitutiveManager();

  string const presDofKey = m_dofManager.getKey( FlowSolverBase::viewKeyStruct::pressureString );
  string const dispDofKey = m_dofManager.getKey( keys::TotalDisplacement );

  globalIndex const rankOffset = m_dofManager.rankOffset();

  CRSMatrixView< real64 const, localIndex const > const
  dFluxResidual_dAperture = m_flowSolver->getDerivativeFluxResidual_dAperture().toViewConst();

  ContactRelationBase const * const
  contactRelation = constitutiveManager.GetGroup< ContactRelationBase >( m_contactRelationName );

  forTargetSubRegionsComplete< FaceElementSubRegion >( mesh,
                                                       [&]( localIndex const,
                                                            localIndex const,
                                                            localIndex const,
                                                            ElementRegionBase const & region,
                                                            FaceElementSubRegion const & subRegion )
  {
    string const & fluidName = m_flowSolver->fluidModelNames()[m_flowSolver->targetRegionIndex( region.getName() )];
    SingleFluidBase const & fluid = GetConstitutiveModel< SingleFluidBase >( subRegion, fluidName );

    arrayView1d< globalIndex const > const presDofNumber = subRegion.getReference< array1d< globalIndex > >( presDofKey );
    arrayView1d< globalIndex const > const dispDofNumber = nodeManager.getReference< array1d< globalIndex > >( dispDofKey );

    arrayView2d< real64 const > const dens = fluid.density();

    arrayView1d< real64 const > const aperture = subRegion.getElementAperture();
    arrayView1d< real64 const > const area = subRegion.getElementArea();

    arrayView2d< localIndex const > const elemsToFaces = subRegion.faceList();
    ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

    arrayView2d< real64 const > const faceNormal = faceManager.faceNormal();

    forAll< serialPolicy >( subRegion.size(), [=]( localIndex ei )
    {
      constexpr int kfSign[2] = { -1, 1 };

      globalIndex const elemDOF = presDofNumber[ei];
      // row number associated to the pressure dof
      globalIndex rowNumber = elemDOF - rankOffset;
      localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elemsToFaces[ei][0] );
      real64 const dAccumulationResidualdAperture = dens[ei][0] * area[ei];

      globalIndex nodeDOF[8 * 3];

      real64 Nbar[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( faceNormal[elemsToFaces[ei][0]] );
      LvArray::tensorOps::subtract< 3 >( Nbar, faceNormal[elemsToFaces[ei][1]] );
      LvArray::tensorOps::normalize< 3 >( Nbar );

      stackArray1d< real64, 24 > dRdU( 2 * numNodesPerFace * 3 );

      // Accumulation derivative
      for( localIndex kf = 0; kf < 2; ++kf )
      {
        for( localIndex a = 0; a < numNodesPerFace; ++a )
        {
          for( int i = 0; i < 3; ++i )
          {
            nodeDOF[kf * 3 * numNodesPerFace + 3 * a + i] = dispDofNumber[faceToNodeMap( elemsToFaces[ei][kf], a )] + i;
            real64 const dGap_dU = kfSign[kf] * Nbar[i] / numNodesPerFace;
            real64 const dAper_dU = contactRelation->dEffectiveAperture_dAperture( aperture[ei] ) * dGap_dU;
            dRdU( kf * 3 * numNodesPerFace + 3 * a + i ) = dAccumulationResidualdAperture * dAper_dU;
          }
        }
      }

      if( rowNumber >= 0  && rowNumber < localMatrix.numRows() )
      {
        localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( rowNumber,
                                                                          nodeDOF,
                                                                          dRdU.data(),
                                                                          2 * numNodesPerFace * 3 );
      }

      // flux derivative
      localIndex const numColumns = dFluxResidual_dAperture.numNonZeros( ei );
      arraySlice1d< localIndex const > const & columns = dFluxResidual_dAperture.getColumns( ei );
      arraySlice1d< real64 const > const & values = dFluxResidual_dAperture.getEntries( ei );

      for( localIndex kfe2 = 0; kfe2 < numColumns; ++kfe2 )
      {
        real64 dRdAper = values[kfe2];
        localIndex const ei2 = columns[kfe2];

        for( localIndex kf = 0; kf < 2; ++kf )
        {
          for( localIndex a = 0; a < numNodesPerFace; ++a )
          {
            for( int i = 0; i < 3; ++i )
            {
              nodeDOF[kf * 3 * numNodesPerFace + 3 * a + i] =
                dispDofNumber[faceToNodeMap( elemsToFaces[ei2][kf], a )] + i;
              real64 const dGap_dU = kfSign[kf] * Nbar[i] / numNodesPerFace;
              real64 const
              dAper_dU = contactRelation->dEffectiveAperture_dAperture( aperture[ei2] ) * dGap_dU;
              dRdU( kf * 3 * numNodesPerFace + 3 * a + i ) = dRdAper * dAper_dU;
            }
          }
        }

        if( rowNumber >= 0 && rowNumber < localMatrix.numRows() )
        {
          localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( rowNumber,
                                                                            nodeDOF,
                                                                            dRdU.data(),
                                                                            2 * numNodesPerFace * 3 );
        }
      }
    } );
  } );
}

void
HydrofractureSolver::
  ApplySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->ApplySystemSolution( dofManager,
                                      localSolution,
                                      scalingFactor,
                                      domain );
  m_flowSolver->ApplySystemSolution( dofManager,
                                     localSolution,
                                     -scalingFactor,
                                     domain );

  UpdateDeformationForCoupling( domain );
}


real64
HydrofractureSolver::ScalingForSystemSolution( DomainPartition const & domain,
                                               DofManager const & dofManager,
                                               arrayView1d< real64 const > const & localSolution )
{
  return m_solidSolver->ScalingForSystemSolution( domain,
                                                  dofManager,
                                                  localSolution );
}

void HydrofractureSolver::SetNextDt( real64 const & currentDt,
                                     real64 & nextDt )
{

  if( m_numResolves[0] == 0 && m_numResolves[1] == 0 )
  {
    this->SetNextDtBasedOnNewtonIter( currentDt, nextDt );
  }
  else
  {
    //SolverBase * const surfaceGenerator =  this->getParent()->GetGroup< SolverBase >( "SurfaceGen" );
    SurfaceGenerator * const surfaceGenerator = m_surfaceGeneratorSolver;

    GEOSX_ERROR_IF( surfaceGenerator == nullptr, this->getName() << ": invalid surface generator solver name: " << std::endl );

    nextDt = surfaceGenerator->GetTimestepRequest() < 1e99 ? surfaceGenerator->GetTimestepRequest() : currentDt;
  }
  GEOSX_LOG_LEVEL_RANK_0( 3, this->getName() << ": nextDt request is "  << nextDt );
}

void HydrofractureSolver::initializeNewFaceElements( DomainPartition const & )
{
//  m_flowSolver->
}

REGISTER_CATALOG_ENTRY( SolverBase, HydrofractureSolver, std::string const &, Group * const )
} /* namespace geosx */
