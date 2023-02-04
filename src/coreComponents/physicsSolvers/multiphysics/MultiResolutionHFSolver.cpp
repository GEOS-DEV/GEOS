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
 * @file MultiResolutionHFSolver.cpp
 *
 */


#include "MultiResolutionHFSolver.hpp"

#include "common/GEOS_RAJA_Interface.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/contact/ContactSelector.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "finiteElement/Kinematics.h"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/SurfaceElementRegion.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"
#include "mesh/mpiCommunications/NeighborCommunicator.hpp"
#include "physicsSolvers/contact/SolidMechanicsEFEMKernelsHelper.hpp"
#include "physicsSolvers/contact/SolidMechanicsEmbeddedFractures.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/surfaceGeneration/SurfaceGenerator.hpp"
#include "physicsSolvers/surfaceGeneration/EmbeddedSurfaceGenerator.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"



namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

MultiResolutionHFSolver::MultiResolutionHFSolver( const string & name,
                                                  Group * const parent ):
  SolverBase( name, parent ),
  m_baseSolverName(),
  m_patchSolverName(),
  m_nodeFixDamage(),
  m_nodeFixDisp(),
  m_fixedDispList(),
  m_maxNumResolves( 10 )
{
  registerWrapper( viewKeyStruct::baseSolverNameString(), &m_baseSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription(
    "Name of the EFEM (SolidMechanicsEmbeddedFractures) solver to be used as the base solver in the MultiResolution scheme" );

  registerWrapper( viewKeyStruct::patchSolverNameString(), &m_patchSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription(
    "Name of the Phase-Field (PhaseFieldFracture) solver to be used as patch solver in the MultiResolution scheme" );

  registerWrapper( viewKeyStruct::maxNumResolvesString(), &m_maxNumResolves ).
    setApplyDefaultValue( 10 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Value to indicate how many resolves may be executed to perform surface generation after the execution of base and patch scale solvers. " );

  registerWrapper( viewKeyStruct::initialBaseTipString(), &m_baseTip ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Position of the initial crack tip in the EFEM mesh." );

  registerWrapper( viewKeyStruct::initialTipElementIndexString(), &m_baseTipElementIndex ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Index of base element ahead of the crack tip." );

}

void MultiResolutionHFSolver::RegisterDataOnMesh( dataRepository::Group & MeshBodies )
{
  GEOSX_UNUSED_VAR( MeshBodies );
}

MultiResolutionHFSolver::~MultiResolutionHFSolver()
{
  // TODO Auto-generated destructor stub
}

void MultiResolutionHFSolver::implicitStepSetup( real64 const & time_n,
                                                 real64 const & dt,
                                                 DomainPartition & domain )
{
  m_baseSolver->implicitStepSetup( time_n, dt, domain );
  m_patchSolver->implicitStepSetup( time_n, dt, domain );
}

void MultiResolutionHFSolver::postProcessInput()
{
  SolverBase::postProcessInput();
  m_baseSolver = &this->getParent().getGroup< SolidMechanicsEmbeddedFractures >( m_baseSolverName );
  m_patchSolver = &this->getParent().getGroup< PhaseFieldFractureSolver >( m_patchSolverName );
}

void MultiResolutionHFSolver::initializePostInitialConditionsPreSubGroups()
{}

void MultiResolutionHFSolver::resetStateToBeginningOfStep( DomainPartition & domain )
{
  m_baseSolver->resetStateToBeginningOfStep( domain );
  m_patchSolver->resetStateToBeginningOfStep( domain );
}

real64 MultiResolutionHFSolver::solverStep( real64 const & time_n,
                                            real64 const & dt,
                                            int const cycleNumber,
                                            DomainPartition & domain )
{
  real64 dtReturn = dt;
  dtReturn = splitOperatorStep( time_n, dt, cycleNumber, domain );
  return dtReturn;
}

// this function will loop over all subdomain nodes and check their distance to the prescribed discrete crack. If the
// distance is smaller than 1 element size (subdomain), we set the damage in this node to be fixed at 1.
void MultiResolutionHFSolver::setInitialCrackDamageBCs( DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
                                                        CRSMatrixView< real64, globalIndex const > const & GEOSX_UNUSED_PARAM( localMatrix ),
                                                        MeshLevel const & GEOSX_UNUSED_PARAM( patch ),
                                                        MeshLevel const & base )
{

  ElementRegionManager const & baseElemManager = base.getElemManager();
  baseElemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & cellElementSubRegion )
  {
    SortedArrayView< localIndex const > const fracturedElements = cellElementSubRegion.fracturedElementsList();
    m_nodeFixDamage.resize( cellElementSubRegion.numNodesPerElement()*fracturedElements.size() );
    localIndex count = 0;
    for( localIndex a : fracturedElements )
    {
      //get all nodes of fracturedElements(a)
      for( localIndex b=0; b < cellElementSubRegion.numNodesPerElement(); b++ )
      {
        //TODO: if meshes arent identical, this should be nodeList(elem_mapped_to_patch(a),b) or check (2)
        localIndex c = cellElementSubRegion.nodeList( a, b );
        //append c = node b of element a
        if( std::find( m_nodeFixDamage.begin(), m_nodeFixDamage.end(), c ) == m_nodeFixDamage.end() ) //avoid repetition
        {
          //TODO: use node_mapped_to_patch(cellElementSubRegion.nodeList(a,b))
          m_nodeFixDamage( count ) = cellElementSubRegion.nodeList( a, b );
          ++count;
        }
      }
    }
    m_nodeFixDamage.resize( count );
  } );
}

// this function will read the patch solution and locate the crack tip to uptade the crack geometry in the base solver
void MultiResolutionHFSolver::findPhaseFieldTip( R1Tensor & tip,
                                                 MeshLevel const & patch )
{
  //reference point must be prescribed in base coordinate system (usually injection source)
  R1Tensor m_referencePoint = {-1.0, 0.0, 0.0}; //add this as input file parameter
  real64 threshold = 0.95;
  //get mpi communicator
  MPI_Comm const & comm = MPI_COMM_GEOSX;
  ElementRegionManager const & patchElemManager = patch.getElemManager();
  //this is risky, not sure FE_TYPE will come from patch
  //each rank has these arrays
  real64 rankMaxDist = 1e-20;
  R1Tensor rankFarthestCenter = {0.0, 0.0, 0.0};
  //loop over all elements in patch and compute elemental average damage
  //observe that patch coordinates are relative, so, reference point should be added
  //for elements in patch
  patchElemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & cellElementSubRegion )
  {
    localIndex numQuadraturePointsPerElem = 0;
    finiteElement::FiniteElementBase const & fe = cellElementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
    finiteElement::FiniteElementDispatchHandler< ALL_FE_TYPES >::dispatch3D( fe, [&] ( auto const finiteElement )
    {
      using FE_TYPE = TYPEOFREF( finiteElement );
      numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;
    } );
    string const & damageModelName = cellElementSubRegion.getReference< string >( PhaseFieldDamageFEM::viewKeyStruct::solidModelNamesString());
    //this call may need to constitutive pass-thru loop to be generalized to multiple damage types
    const constitutive::Damage< ElasticIsotropic > & damageUpdates = cellElementSubRegion.getConstitutiveModel< Damage< ElasticIsotropic > >( damageModelName );
    arrayView2d< const real64 > allElemCenters = cellElementSubRegion.getElementCenter();
    const arrayView2d< real64 const > qp_damage = damageUpdates.getDamage();
    forAll< serialPolicy >( cellElementSubRegion.size(), [&] ( localIndex const k )
    {
      //compute elemental averaged damage
      real64 average_d = 0;
      R1Tensor elemCenter;
      elemCenter[0] = allElemCenters[k][0];
      elemCenter[1] = allElemCenters[k][1];
      elemCenter[2] = allElemCenters[k][2]; //this is trying to create a R1Tensor from a array2d< real64 >
      for( localIndex q=0; q<numQuadraturePointsPerElem; q++ )
      {
        //get damage at quadrature point i
        average_d = average_d + qp_damage( k, q )/numQuadraturePointsPerElem;
      }
      //if elemental damage > 0.95 (or another threshold)
      if( average_d > threshold )
      {
        //check if this element is farther than current farthest
        R1Tensor elemVec = LVARRAY_TENSOROPS_INIT_LOCAL_3( m_referencePoint );
        LvArray::tensorOps::subtract< 3 >( elemVec, elemCenter );
        real64 dist = LvArray::tensorOps::l2Norm< 3 >( elemVec );
        if( dist > rankMaxDist )
        {
          rankMaxDist = dist;
          rankFarthestCenter = elemCenter;
        }
      }
    } );
  } );

  real64 globalMax = MpiWrapper::max< real64 >( rankMaxDist, comm );
  if( std::abs( rankMaxDist-globalMax )<1e-12 )
  {
    tip = rankFarthestCenter;
  }

}

//The function prepares a list of dofs and u values that will be used by the patch solver to set the boundary conditions
void MultiResolutionHFSolver::prepareSubProblemBCs( MeshLevel const & base,
                                                    MeshLevel & patch )
{

  // get list of nodes on the boundary of the patch
  FaceManager const & patchFaceManager = patch.getFaceManager();
  NodeManager & patchNodeManager = patch.getNodeManager();
  patchNodeManager.setIsExternal( patchFaceManager );
  SortedArray< localIndex > patchExternalSet = patchNodeManager.externalSet();
  arrayView1d< real64 const > const patchDamage = patchNodeManager.getReference< array1d< real64 > >( "Damage" );
  NodeManager const & baseNodeManager = base.getNodeManager();
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const baseDisp = baseNodeManager.getField< fields::solidMechanics::totalDisplacement >();
  real64 damage_threshold = 0.3;
  m_nodeFixDisp.resize( patchExternalSet.size() );
  m_fixedDispList.resize( patchExternalSet.size(), 3 );
  localIndex count=0;
  for( localIndex a : patchExternalSet )
  {
    if( patchDamage[a] < damage_threshold )
    {
      // NOTE: there needs to be a translation between patch and base mesh for the indices and values/weights.

      //append node a
      //USE REGISTERED ARRAY HERE
      m_nodeFixDisp( count ) = a;
      localIndex const numBaseNodes = 1;  //TODO: this wont be 1 in other cases, m_nodeMapIndices.sizeOfArray( a );
      for( localIndex b=0; b<numBaseNodes; ++b )
      {
        //write displacements from the base domain to become a boundary condition in the patch domain
        //TODO: a should be replaced by mapped(a) if meshes arent identical
        // mapped(a) is the node in base that has the same physical coordinates of a in patch
        m_fixedDispList( count, 0 ) = baseDisp( a, 0 );
        m_fixedDispList( count, 1 ) = baseDisp( a, 1 );
        m_fixedDispList( count, 2 ) = baseDisp( a, 2 );
      }

      ++count;
    }
  }
  //USE REGISTERED ARRAY HERE
  m_nodeFixDisp.resize( count );
  m_fixedDispList.resize( count, 3 );

}

real64 MultiResolutionHFSolver::splitOperatorStep( real64 const & time_n,
                                                   real64 const & dt,
                                                   integer const cycleNumber,
                                                   DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  real64 dtReturn = dt;
  real64 dtReturnTemporary;

  SolidMechanicsEmbeddedFractures &
  baseSolver = this->getParent().getGroup< SolidMechanicsEmbeddedFractures >( m_baseSolverName );

  PhaseFieldFractureSolver &
  patchSolver = this->getParent().getGroup< PhaseFieldFractureSolver >( m_patchSolverName );

  EmbeddedSurfaceGenerator &
  efemGenerator = this->getParent().getGroup< EmbeddedSurfaceGenerator >( "SurfaceGenerator" ); //this is hard coded

  PhaseFieldDamageFEM &
  patchDamageSolver = *patchSolver.damageSolver();

  SolidMechanicsLagrangianFEM &
  patchSolidSolver = *patchSolver.solidMechanicsSolver();


  baseSolver.setupSystem( domain,
                          baseSolver.getDofManager(),
                          baseSolver.getLocalMatrix(),
                          baseSolver.getSystemRhs(),
                          baseSolver.getSystemSolution(),
                          true );

  baseSolver.implicitStepSetup( time_n, dt, domain );

  this->implicitStepSetup( time_n, dt, domain );

  NonlinearSolverParameters & solverParams = getNonlinearSolverParameters();
  //although these iterations are not really Newton iterations, we will use this nomeclature to keep things consistent
  integer & iter = solverParams.m_numNewtonIterations;
  iter = 0;
  bool isConverged = false;
  while( iter < solverParams.m_maxIterNewton )
  {
    if( iter == 0 )
    {
      // reset the states of all slave solvers if any of them has been reset
      //TODO: this is potentially a code duplication since resetStateToBeginningOfStep(domain) already calls the slaves
      patchSolver.resetStateToBeginningOfStep( domain );
      baseSolver.resetStateToBeginningOfStep( domain );
      resetStateToBeginningOfStep( domain );
    }

    GEOSX_LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1 << ", BaseSolver: " );

    //we probably want to run a phase-field solve in the patch problem at timestep 0 to get a smooth initial crack. Also, re-run this
    // anytime the base crack changes
    map< std::pair< string, string >, array1d< string > > const & baseTargets = baseSolver.getReference< map< std::pair< string, string >, array1d< string > > >(
      SolverBase::viewKeyStruct::meshTargetsString());
    auto const baseTarget = baseTargets.begin()->first;
    map< std::pair< string, string >, array1d< string > > const & patchTargets = patchSolver.getReference< map< std::pair< string, string >, array1d< string > > >(
      SolverBase::viewKeyStruct::meshTargetsString());
    auto const patchTarget = patchTargets.begin()->first;
    CRSMatrix< real64, globalIndex > & patchDamageLocalMatrix = patchDamageSolver.getLocalMatrix();
    this->setInitialCrackDamageBCs( patchDamageSolver.getDofManager(), patchDamageLocalMatrix.toViewConstSizes(), domain.getMeshBody( patchTarget.first ).getBaseDiscretization(),
                                    domain.getMeshBody( baseTarget.first ).getBaseDiscretization() );
    patchDamageSolver.setInitialCrackNodes( m_nodeFixDamage );

    //now perform the subproblem run with no BCs on displacements, just to set the damage inital condition;
    if( iter == 0 )
    {
      //TODO: eventually, we will need to update the patch domain
      real64 dtUseless = patchSolver.solverStep( time_n,
                                                 dtReturn,
                                                 cycleNumber,
                                                 domain );
      GEOSX_UNUSED_VAR( dtUseless );
    }

    dtReturnTemporary = baseSolver.nonlinearImplicitStep( time_n,
                                                          dtReturn,
                                                          cycleNumber,
                                                          domain );

    if( dtReturnTemporary < dtReturn )
    {
      iter = 0;
      dtReturn = dtReturnTemporary;
      continue;
    }

    if( baseSolver.getNonlinearSolverParameters().m_numNewtonIterations >= 1 && iter > 0 )
    {
      GEOSX_LOG_LEVEL_RANK_0( 1, "***** The Global-Local iterative scheme has converged in " << iter << " iterations! *****\n" );
      isConverged = true;
      break;
    }

    //here, before calling the nonlinarImplicitStep of the patch solver, we must prescribe the displacement boundary conditions
    this->prepareSubProblemBCs( domain.getMeshBody( baseTarget.first ).getBaseDiscretization(), domain.getMeshBody( patchTarget.first ).getBaseDiscretization());

    //write disp BCs to local disp solver
    patchSolidSolver.setInternalBoundaryConditions( m_nodeFixDisp, m_fixedDispList );

    //finally, pass the shared boundary information to the patch solver
    //patchSolver.

    GEOSX_LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1 << ", PatchSolver: " );

    dtReturnTemporary = patchSolver.solverStep( time_n,
                                                dtReturn,
                                                cycleNumber,
                                                domain );
    this->findPhaseFieldTip( m_patchTip, domain.getMeshBody( patchTarget.first ).getBaseDiscretization());

    if( time_n > 0 )
    {
      efemGenerator.propagationStep( domain, m_baseTip, m_patchTip, m_baseTipElementIndex );
      baseSolver.setupSystem( domain,
                              baseSolver.getDofManager(),
                              baseSolver.getLocalMatrix(),
                              baseSolver.getSystemRhs(),
                              baseSolver.getSystemSolution(),
                              true );
      baseSolver.implicitStepSetup( time_n, dt, domain );
    }
    GEOSX_LOG_LEVEL_RANK_0( 2, "baseTipElement: "<<m_baseTipElementIndex );
    GEOSX_LOG_LEVEL_RANK_0( 2, "PFtipX: "<<m_patchTip[0] );
    GEOSX_LOG_LEVEL_RANK_0( 2, "PFtipY: "<<m_patchTip[1] );
    GEOSX_LOG_LEVEL_RANK_0( 2, "PFtipZ: "<<m_patchTip[2] );
    GEOSX_LOG_LEVEL_RANK_0( 2, "EFtipX: "<<m_baseTip[0] );
    GEOSX_LOG_LEVEL_RANK_0( 2, "EFtipY: "<<m_baseTip[1] );
    GEOSX_LOG_LEVEL_RANK_0( 2, "EFtipZ: "<<m_baseTip[2] );

    if( dtReturnTemporary < dtReturn )
    {
      iter = 0;
      dtReturn = dtReturnTemporary;
      continue;
    }

    ++iter;
  }

  GEOSX_ERROR_IF( !isConverged, "MultiResolutionHFSolver::SplitOperatorStep() did not converge" );

  baseSolver.implicitStepComplete( time_n, dt, domain );
  patchSolver.implicitStepComplete( time_n, dt, domain );

  return dtReturn;
}

REGISTER_CATALOG_ENTRY( SolverBase, MultiResolutionHFSolver, string const &, Group * const )
} /* namespace geosx */
