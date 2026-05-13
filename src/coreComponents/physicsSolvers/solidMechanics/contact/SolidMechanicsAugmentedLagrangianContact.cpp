/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/*
 * SolidMechanicsAugmentedLagrangianContact.cpp
 */

#include "SolidMechanicsAugmentedLagrangianContact.hpp"

#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBaseFields.hpp"

#include "physicsSolvers/solidMechanics/contact/kernels/SolidMechanicsConformingContactKernelsBase.hpp"
#include "physicsSolvers/solidMechanics/contact/kernels/SolidMechanicsALMKernels.hpp"
#include "physicsSolvers/solidMechanics/contact/kernels/SolidMechanicsALMKernelsBase.hpp"
#include "physicsSolvers/solidMechanics/contact/kernels/SolidMechanicsALMSimultaneousKernels.hpp"
#include "physicsSolvers/solidMechanics/contact/kernels/SolidMechanicsDisplacementJumpUpdateKernels.hpp"
#include "physicsSolvers/solidMechanics/contact/kernels/SolidMechanicsConformingPressureContributionKernels.hpp"
#include "physicsSolvers/solidMechanics/contact/kernels/SolidMechanicsContactFaceBubbleKernels.hpp"
#include "physicsSolvers/solidMechanics/contact/kernels/SolidMechanicsPressureFaceBubbleKernels.hpp"
#include "physicsSolvers/solidMechanics/contact/LogLevelsInfo.hpp"
#include "physicsSolvers/solidMechanics/contact/ContactFields.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsFields.hpp"
#include "physicsSolvers/LogLevelsInfo.hpp"
#include "constitutive/contact/FrictionSelector.hpp"
#include "constitutive/solid/PorousSolid.hpp"
#include "constitutive/solid/SolidFields.hpp"
#include "finiteElement/FiniteElementDiscretization.hpp"
#include "mesh/DomainPartition.hpp"

#include <stdio.h>

#if defined( GEOS_USE_CUDA )
#include <cuda_runtime.h>
#endif

namespace geos
{

using namespace constitutive;
using namespace dataRepository;
using namespace fields;

// Workaround for nvcc bug: forDiscretizationOnMeshTargets lambdas receive
// string_array const & (= std::vector<std::string>) as a parameter.  When the
// lambda body also contains device kernel launches, nvcc erroneously tries to
// generate a device-compatible destructor for std::vector<std::string> even
// though it is a reference parameter whose lifetime is not managed by the lambda.
GEOS_NV_HOST_DEVICE_DIAG_SUPPRESS

SolidMechanicsAugmentedLagrangianContact::SolidMechanicsAugmentedLagrangianContact( const string & name,
                                                                                    Group * const parent ):
  ContactSolverBase( name, parent )
{
  registerWrapper( viewKeyStruct::simultaneousString(), &m_simultaneous ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1 ).
    setDescription( "Flag to update the Lagrange multiplier at each Newton iteration (true), or only after the Newton loop has converged (false)" );

  registerWrapper( viewKeyStruct::symmetricString(), &m_symmetric ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1 ).
    setDescription( "Flag to neglect the non-symmetric contribution in the tangential matrix" );

  registerWrapper( viewKeyStruct::iterativePenaltyNFacString(), &m_iterPenaltyNFac ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 10.0 ).
    setDescription( "Factor for tuning the iterative penalty coefficient for normal traction" );

  registerWrapper( viewKeyStruct::iterativePenaltyTFacString(), &m_iterPenaltyTFac ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.1 ).
    setDescription( "Factor for tuning the iterative penalty coefficient for tangential traction" );

  registerWrapper( viewKeyStruct::tolJumpDispNFacString(), &m_tolJumpDispNFac ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1.e-07 ).
    setDescription( "Factor to adjust the tolerance for normal jump" );

  registerWrapper( viewKeyStruct::tolJumpDispTFacString(), &m_tolJumpDispTFac ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1.e-05 ).
    setDescription( "Factor to adjust the tolerance for tangential jump" );

  registerWrapper( viewKeyStruct::tolNormalTracFacString(), &m_tolNormalTracFac ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 0.5 ).
    setDescription( "Factor to adjust the tolerance for normal traction" );

  registerWrapper( viewKeyStruct::tolTauLimitString(), &m_slidingCheckTolerance ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 5.e-02 ).
    setDescription( "Tolerance for the sliding check" );

  // Set the default linear solver parameters
  LinearSolverParameters & linSolParams = m_linearSolverParameters.get();

  // Strategy: AMG with separate displacement components
  linSolParams.dofsPerNode = 3;
  linSolParams.isSymmetric = true;
  linSolParams.amg.separateComponents = true;

  // Strategy: static condensation of bubble dofs using MGR
  linSolParams.mgr.strategy = LinearSolverParameters::MGR::StrategyType::augmentedLagrangianContactMechanics;
  linSolParams.mgr.separateComponents = true;
}

SolidMechanicsAugmentedLagrangianContact::~SolidMechanicsAugmentedLagrangianContact()
{}

void SolidMechanicsAugmentedLagrangianContact::registerDataOnMesh( dataRepository::Group & meshBodies )
{
  ContactSolverBase::registerDataOnMesh( meshBodies );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & meshLevel,
                                                    string_array const & )
  {
    FaceManager & faceManager = meshLevel.getFaceManager();

    // Register the total bubble displacement
    faceManager.registerField< contact::totalBubbleDisplacement >( getName() ).
      reference().resizeDimension< 1 >( 3 );

    // Register the incremental bubble displacement
    faceManager.registerField< contact::incrementalBubbleDisplacement >( getName() ).
      reference().resizeDimension< 1 >( 3 );
  } );

  forFractureRegionOnMeshTargets( meshBodies, [&] ( SurfaceElementRegion & fractureRegion )
  {
    fractureRegion.forElementSubRegions< SurfaceElementSubRegion >( [&]( SurfaceElementSubRegion & subRegion )
    {
      subRegion.registerField< contact::deltaTraction >( getName() ).
        reference().resizeDimension< 1 >( 3 );

      // Register the rotation matrix
      subRegion.registerField< contact::rotationMatrix >( getName() ).
        reference().resizeDimension< 1, 2 >( 3, 3 );

      // Register the penalty coefficients for the iterative procedure
      subRegion.registerField< contact::iterativePenalty >( getName() ).
        reference().resizeDimension< 1 >( 5 );

      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::normalTractionToleranceString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRegisteringObjects( getName()).
        setDescription( "An array that holds the normal traction tolerance." );

      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::normalDisplacementToleranceString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRegisteringObjects( getName()).
        setDescription( "An array that holds the normal displacement tolerance." );

      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::slidingToleranceString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRegisteringObjects( getName()).
        setDescription( "An array that holds the sliding tolerance." );

      subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::dispJumpUpdPenaltyString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRegisteringObjects( getName()).
        setDescription( "An array that stores the displacement jumps used to update the penalty coefficients." ).
        reference().resizeDimension< 1 >( 3 );

      // Register pressure fields for sequential poromechanics coupling and stress initialization
      // In coupled poromechanics, the flow solver will overwrite these with actual values
      subRegion.registerField< flow::pressure >( getName() ).
        setApplyDefaultValue( 0.0 ).
        setPlotLevel( PlotLevel::NOPLOT );
      subRegion.registerField< flow::pressure_n >( getName() ).
        setApplyDefaultValue( 0.0 ).
        setPlotLevel( PlotLevel::NOPLOT );
    } );
  } );

}

void SolidMechanicsAugmentedLagrangianContact::initializePostInitialConditionsPreSubGroups()
{
  ContactSolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  validateTetrahedralQuadrature( domain.getMeshBodies() );
}

void SolidMechanicsAugmentedLagrangianContact::validateTetrahedralQuadrature( Group & meshBodies )
{
  string const discretizationName = getDiscretizationName();

  NumericalMethodsManager const & numericalMethodManager =
    this->getGroupByPath< DomainPartition >( "/Problem/domain" ).getNumericalMethodManager();
  FiniteElementDiscretizationManager const & feDiscretizationManager =
    numericalMethodManager.getFiniteElementDiscretizationManager();
  FiniteElementDiscretization const & feDiscretization =
    feDiscretizationManager.getGroup< FiniteElementDiscretization >( discretizationName );


  integer const useHighOrderQuadrature =
    feDiscretization.getReference< integer >( "useHighOrderQuadratureRule" );

  bool hasTetrahedra = false;
  forDiscretizationOnMeshTargets( meshBodies, [&]( string const &,
                                                   MeshLevel const & mesh,
                                                   string_array const & regionNames )
  {
    ElementRegionManager const & elemManager = mesh.getElemManager();
    elemManager.forElementRegions< CellElementRegion >( regionNames, [&]( localIndex const,
                                                                          CellElementRegion const & region )
    {
      region.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & subRegion )
      {
        if( subRegion.getElementType() == ElementType::Tetrahedron )
        {
          hasTetrahedra = true;
        }
      } );
    } );
  } );

  GEOS_ERROR_IF( hasTetrahedra && useHighOrderQuadrature != 1,
                 GEOS_FMT( "{}: Tetrahedral meshes require useHighOrderQuadratureRule=\"1\" for correct integration of bubble contributions. "
                           "Please add this attribute to your FiniteElements/{} XML block.",
                           getName(), discretizationName ) );
}

void SolidMechanicsAugmentedLagrangianContact::setupDofs( DomainPartition const & domain,
                                                          DofManager & dofManager ) const
{

  GEOS_MARK_FUNCTION;
  SolidMechanicsLagrangianFEM::setupDofs( domain, dofManager );

  map< std::pair< string, string >, string_array > meshTargets;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshBodyName,
                                                                MeshLevel const & meshLevel,
                                                                string_array const & )
  {
    string_array regions;
    regions.emplace_back( getUniqueFractureRegionName() );
    meshTargets[std::make_pair( meshBodyName, meshLevel.getName())] = std::move( regions );
  } );

  dofManager.addField( contact::totalBubbleDisplacement::key(),
                       FieldLocation::Face,
                       3,
                       meshTargets );

  // Add coupling between bubble
  // Useful to create connection between bubble dofs for Augmented Lagrangian formulation
  dofManager.addCoupling( contact::totalBubbleDisplacement::key(),
                          contact::totalBubbleDisplacement::key(),
                          DofManager::Connector::Elem );

}

void SolidMechanicsAugmentedLagrangianContact::setupSystem( DomainPartition & domain,
                                                            DofManager & dofManager,
                                                            CRSMatrix< real64, globalIndex > & localMatrix,
                                                            ParallelVector & rhs,
                                                            ParallelVector & solution,
                                                            bool const setSparsity )
{

  GEOS_MARK_FUNCTION;

  // Recompute geometric quantities (face normals, areas) after mesh topology changes.
  // This is critical for distorted/non-axis-aligned meshes and after fracture events.
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {
    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager & faceManager = mesh.getFaceManager();
    ElementRegionManager & elemManager = mesh.getElemManager();

    // Recompute face geometry (normals and areas)
    faceManager.computeGeometry( nodeManager );

    // Recompute element geometry for volume elements (centers and volumes)
    elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
    {
      subRegion.calculateElementGeometricQuantities( nodeManager, faceManager );
    } );

    // Recompute element geometry for face elements (uses updated face areas)
    elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
    {
      subRegion.calculateElementGeometricQuantities( nodeManager, faceManager );

      // Reorder kf1 nodes to match kf0 for conforming contact kernels.
      // flipFaceMap and fixNeighboringFacesNormals are already called by
      // ProblemManager::applyNumericalMethods after ghosting is complete.
      if( subRegion.size() > 0 )
      {
        subRegion.fixNeighboringFacesNormals( faceManager, elemManager );
        subRegion.orderKf1NodesConsistentlyWithKf0( faceManager, nodeManager );
      }
    } );
  } );

  // Create the lists of interface elements that have same type.
  createFaceTypeList( domain );

  // Create the lists of interface elements that have same type and same fracture state.
  updateStickSlipList( domain );

  // Create the list of cell elements that they are enriched with bubble functions.
  createBubbleCellList( domain );

  PhysicsSolverBase::setupSystem( domain, dofManager, localMatrix, rhs, solution, setSparsity );
}

void SolidMechanicsAugmentedLagrangianContact::postInputInitialization()
{
  ContactSolverBase::postInputInitialization();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteElementDiscretizationManager const & feDiscretizationManager =
    numericalMethodManager.getFiniteElementDiscretizationManager();
  FiniteElementDiscretization const & feDiscretization =
    feDiscretizationManager.getGroup< FiniteElementDiscretization >( getDiscretizationName() );

  m_faceTypeToFiniteElements.insert( {"Quadrilateral", feDiscretization.factory( ElementType::Quadrilateral )} );
  m_faceTypeToFiniteElements.insert( {"Triangle", feDiscretization.factory( ElementType::Triangle )} );

  GEOS_LOG_RANK_0( GEOS_FMT( "{} using finite element discretization {}:",
                             getName(), getDiscretizationName() ) );
  for( auto const & [name, fePtr] : m_faceTypeToFiniteElements )
  {
    GEOS_LOG_RANK_0( GEOS_FMT( "  {} face elements: {} quadrature points",
                               name, fePtr->getNumQuadraturePoints() ) );
  }
}

void SolidMechanicsAugmentedLagrangianContact::setSparsityPattern( DomainPartition & domain,
                                                                   DofManager & dofManager,
                                                                   CRSMatrix< real64, globalIndex > & GEOS_UNUSED_PARAM( localMatrix ),
                                                                   SparsityPattern< globalIndex > & pattern )
{
  // Set the sparsity pattern without the Abu and Aub blocks.
  SparsityPattern< globalIndex > patternDiag;
  dofManager.setSparsityPattern( patternDiag );

  // Get the original row lengths (diagonal blocks only)
  array1d< localIndex > rowLengths( patternDiag.numRows());
  for( localIndex localRow = 0; localRow < patternDiag.numRows(); ++localRow )
  {
    rowLengths[localRow] = patternDiag.numNonZeros( localRow );
  }

  // Add the number of nonzeros induced by coupling
  addCouplingNumNonzeros( domain, dofManager, rowLengths.toView());

  // Create a new pattern with enough capacity for coupled matrix
  pattern.resizeFromRowCapacities< parallelHostPolicy >( patternDiag.numRows(), patternDiag.numColumns(), rowLengths.data());

  // Copy the original nonzeros
  for( localIndex localRow = 0; localRow < patternDiag.numRows(); ++localRow )
  {
    globalIndex const * cols = patternDiag.getColumns( localRow ).dataIfContiguous();
    pattern.insertNonZeros( localRow, cols, cols + patternDiag.numNonZeros( localRow ));
  }

  // Add the nonzeros from coupling
  addCouplingSparsityPattern( domain, dofManager, pattern.toView());
}

void SolidMechanicsAugmentedLagrangianContact::implicitStepSetup( real64 const & time_n,
                                                                  real64 const & dt,
                                                                  DomainPartition & domain )
{

  SolidMechanicsLagrangianFEM::implicitStepSetup( time_n, dt, domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {

    FaceManager & faceManager = mesh.getFaceManager();
    ElementRegionManager & elemManager = mesh.getElemManager();

    SurfaceElementRegion & region = elemManager.getRegion< SurfaceElementRegion >( getUniqueFractureRegionName() );
    FaceElementSubRegion & subRegion = region.getUniqueSubRegion< FaceElementSubRegion >();

    arrayView2d< real64 const > const faceNormal = faceManager.faceNormal();
    arrayView2d< localIndex const > const elemsToFaces = subRegion.faceList().toViewConst();

    arrayView2d< real64 > const incrBubbleDisp =
      faceManager.getField< contact::incrementalBubbleDisplacement >();

    arrayView3d< real64 > const
    rotationMatrix = subRegion.getField< contact::rotationMatrix >().toView();

    arrayView2d< real64 > const unitNormal   = subRegion.getNormalVector();
    arrayView2d< real64 > const unitTangent1 = subRegion.getTangentVector1();
    arrayView2d< real64 > const unitTangent2 = subRegion.getTangentVector2();

    if( subRegion.size() > 0 )
    {
      solidMechanicsConformingContactKernels::ComputeRotationMatricesKernel::
        launch< parallelDevicePolicy<> >( subRegion.size(),
                                          faceNormal,
                                          elemsToFaces,
                                          rotationMatrix,
                                          unitNormal,
                                          unitTangent1,
                                          unitTangent2 );
    }

    // Set the tollerances
    computeTolerances( domain );

    // Initialize the traction from the stress in adjacent volume elements.
    // This is only done during stress initialization step to ensure the ALM solver
    // starts with a physically consistent traction rather than zero.
    // On subsequent time steps, the traction from the previous step is preserved.
    if( m_performStressInitialization )
    {
      initializeTractionFromAdjacentCellStress( domain );
    }

    // Set array to update penalty coefficients
    arrayView2d< real64 > const dispJumpUpdPenalty =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::dispJumpUpdPenaltyString() );

    arrayView2d< real64 > const
    iterativePenalty = subRegion.getField< contact::iterativePenalty >().toView();
    arrayView1d< integer const > const fractureState = subRegion.getField< contact::fractureState >();

    if( subRegion.size() > 0 )
    {
      if( m_simultaneous )
      {
        forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_DEVICE ( localIndex const k )
        {
          if( fractureState[k] == contact::FractureState::Stick )
          {
            iterativePenalty[k][2] = iterativePenalty[k][1];
            iterativePenalty[k][3] = iterativePenalty[k][1];
            iterativePenalty[k][4] = 0.0;
          }
          else
          {
            iterativePenalty[k][2] = 0.0;
            iterativePenalty[k][3] = 0.0;
            iterativePenalty[k][4] = 0.0;
          }
        } );
      }

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_DEVICE ( localIndex const k )
      {
        LvArray::tensorOps::fill< 3 >( dispJumpUpdPenalty[k], 0.0 );
        localIndex const kf0 = elemsToFaces[k][0];
        localIndex const kf1 = elemsToFaces[k][1];
        LvArray::tensorOps::fill< 3 >( incrBubbleDisp[kf0], 0.0 );
        LvArray::tensorOps::fill< 3 >( incrBubbleDisp[kf1], 0.0 );
      } );
    }
  } );

  // Sync iterativePenalty
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {
    FieldIdentifiers fieldsToBeSync;

    fieldsToBeSync.addElementFields( { contact::iterativePenalty::key() },
                                     { getUniqueFractureRegionName() } );

    // Synchronize bubble displacement fields to ensure ghost values are initialized
    // This prevents race conditions when computing residuals in parallel
    fieldsToBeSync.addFields( FieldLocation::Face,
                              { contact::incrementalBubbleDisplacement::key(),
                                contact::totalBubbleDisplacement::key() } );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                         mesh,
                                                         domain.getNeighbors(),
                                                         true );
  } );

}

void SolidMechanicsAugmentedLagrangianContact::assembleSystem( real64 const time,
                                                               real64 const dt,
                                                               DomainPartition & domain,
                                                               DofManager const & dofManager,
                                                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                               arrayView1d< real64 > const & localRhs )
{

  GEOS_MARK_FUNCTION;

  synchronizeFractureState( domain );

  SolidMechanicsLagrangianFEM::assembleSystem( time,
                                               dt,
                                               domain,
                                               dofManager,
                                               localMatrix,
                                               localRhs );

  assembleContact( time, dt, domain, dofManager, localMatrix, localRhs );

  // for sequential: add (fixed) pressure force contribution into residual (no derivatives)
  if( m_isFixedStressPoromechanicsUpdate || m_performStressInitialization )
  {
    assembleForceResidualPressureContribution( domain, dt, dofManager, localMatrix, localRhs );
  }

}

void SolidMechanicsAugmentedLagrangianContact::assembleContact( real64 const time,
                                                                real64 const dt,
                                                                DomainPartition & domain,
                                                                DofManager const & dofManager,
                                                                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                arrayView1d< real64 > const & localRhs )
{

  GEOS_MARK_FUNCTION;

  GEOS_UNUSED_VAR( time );

  // Loop for assembling contributes from interface elements (Aut*eps^-1*Atu and Aub*eps^-1*Abu)
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshName,
                                                                MeshLevel & mesh,
                                                                string_array const & )

  {

    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    string const & dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );
    string const & bubbleDofKey = dofManager.getKey( contact::totalBubbleDisplacement::key() );

    arrayView1d< globalIndex const > const dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
    arrayView1d< globalIndex const > const bubbleDofNumber = faceManager.getReference< globalIndex_array >( bubbleDofKey );

    string const & fractureRegionName = getUniqueFractureRegionName();

    forFiniteElementOnStickFractureSubRegions( meshName, [&] ( string const &,
                                                               finiteElement::FiniteElementBase const & subRegionFE,
                                                               arrayView1d< localIndex const > const & faceElementList,
                                                               bool const )
    {

      if( m_simultaneous )
      {
        solidMechanicsALMKernels::ALMSimultaneousFactory kernelFactory( dispDofNumber,
                                                                        bubbleDofNumber,
                                                                        dofManager.rankOffset(),
                                                                        localMatrix,
                                                                        localRhs,
                                                                        dt,
                                                                        faceElementList );

        real64 maxTraction = finiteElement::interfaceBasedKernelApplication< parallelDevicePolicy< >, CoulombFriction >( mesh,
                                                                                                                         fractureRegionName,
                                                                                                                         faceElementList,
                                                                                                                         subRegionFE,
                                                                                                                         viewKeyStruct::frictionLawNameString(),
                                                                                                                         kernelFactory );

        GEOS_UNUSED_VAR( maxTraction );

      }
      else
      {
        solidMechanicsALMKernels::ALMFactory kernelFactory( dispDofNumber,
                                                            bubbleDofNumber,
                                                            dofManager.rankOffset(),
                                                            localMatrix,
                                                            localRhs,
                                                            dt,
                                                            faceElementList,
                                                            m_symmetric );

        real64 maxTraction = finiteElement::
                               interfaceBasedKernelApplication
                             < parallelDevicePolicy< >,
                               CoulombFriction >( mesh,
                                                  fractureRegionName,
                                                  faceElementList,
                                                  subRegionFE,
                                                  viewKeyStruct::frictionLawNameString(),
                                                  kernelFactory );

        GEOS_UNUSED_VAR( maxTraction );
      }

    } );

    forFiniteElementOnSlipFractureSubRegions( meshName, [&] ( string const &,
                                                              finiteElement::FiniteElementBase const & subRegionFE,
                                                              arrayView1d< localIndex const > const & faceElementList,
                                                              bool const )
    {

      if( m_simultaneous )
      {
        solidMechanicsALMKernels::ALMSimultaneousFactory kernelFactory( dispDofNumber,
                                                                        bubbleDofNumber,
                                                                        dofManager.rankOffset(),
                                                                        localMatrix,
                                                                        localRhs,
                                                                        dt,
                                                                        faceElementList );

        real64 maxTraction = finiteElement::
                               interfaceBasedKernelApplication
                             < parallelDevicePolicy< >,
                               CoulombFriction >( mesh,
                                                  fractureRegionName,
                                                  faceElementList,
                                                  subRegionFE,
                                                  viewKeyStruct::frictionLawNameString(),
                                                  kernelFactory );

        GEOS_UNUSED_VAR( maxTraction );

      }
      else
      {
        solidMechanicsALMKernels::ALMFactory kernelFactory( dispDofNumber,
                                                            bubbleDofNumber,
                                                            dofManager.rankOffset(),
                                                            localMatrix,
                                                            localRhs,
                                                            dt,
                                                            faceElementList,
                                                            m_symmetric );

        real64 maxTraction = finiteElement::interfaceBasedKernelApplication< parallelDevicePolicy< >, CoulombFriction >( mesh,
                                                                                                                         fractureRegionName,
                                                                                                                         faceElementList,
                                                                                                                         subRegionFE,
                                                                                                                         viewKeyStruct::frictionLawNameString(),
                                                                                                                         kernelFactory );

        GEOS_UNUSED_VAR( maxTraction );
      }

    } );

  } );

  // Loop for assembling contributes of bubble elements (Abb, Abu, Aub)
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {
    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    string const & dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );
    string const & bubbleDofKey = dofManager.getKey( contact::totalBubbleDisplacement::key() );

    arrayView1d< globalIndex const > const dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
    arrayView1d< globalIndex const > const bubbleDofNumber = faceManager.getReference< globalIndex_array >( bubbleDofKey );

    real64 const gravityVectorData[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );


    solidMechanicsConformingContactKernels::FaceBubbleFactory kernelFactory( dispDofNumber,
                                                                             bubbleDofNumber,
                                                                             dofManager.rankOffset(),
                                                                             localMatrix,
                                                                             localRhs,
                                                                             dt,
                                                                             gravityVectorData );

    real64 maxTraction = finiteElement::regionBasedKernelApplication< parallelDevicePolicy< >, ElasticIsotropic, CellElementSubRegion >( mesh,
                                                                                                                                         regionNames,
                                                                                                                                         getDiscretizationName(),
                                                                                                                                         SolidMechanicsLagrangianFEM::viewKeyStruct::
                                                                                                                                           solidMaterialNamesString(),
                                                                                                                                         kernelFactory );

    GEOS_UNUSED_VAR( maxTraction );

  } );

}

void SolidMechanicsAugmentedLagrangianContact::assembleForceResidualPressureContribution( DomainPartition & domain,
                                                                                          real64 const & dt,
                                                                                          DofManager const & dofManager,
                                                                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                                          arrayView1d< real64 > const & localRhs )
{

  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshName,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {

    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    string const & dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );
    string const & bubbleDofKey = dofManager.getKey( contact::totalBubbleDisplacement::key() );

    arrayView1d< globalIndex const > const dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
    arrayView1d< globalIndex const > const bubbleDofNumber = faceManager.getReference< globalIndex_array >( bubbleDofKey );

    string const & fractureRegionName = this->getUniqueFractureRegionName();

    forFiniteElementOnFractureSubRegions( meshName, [&] ( string const &,
                                                          finiteElement::FiniteElementBase const & subRegionFE,
                                                          arrayView1d< localIndex const > const & faceElementList )
    {

      solidMechanicsConformingContactKernels::AssemblePressureContributionFactory
      kernelFactory( dispDofNumber,
                     bubbleDofNumber,
                     dofManager.rankOffset(),
                     localMatrix,
                     localRhs,
                     dt,
                     faceElementList );

      real64 maxTraction = finiteElement::
                             interfaceBasedKernelApplication
                           < parallelDevicePolicy< >,
                             constitutive::NullModel >( mesh,
                                                        fractureRegionName,
                                                        faceElementList,
                                                        subRegionFE,
                                                        "",
                                                        kernelFactory );

      GEOS_UNUSED_VAR( maxTraction );

    } );

    set< string > poromechanicsRegions;
    ElementRegionManager const & elementRegionManager = mesh.getElemManager();
    elementRegionManager.forElementSubRegions< CellElementSubRegion >( regionNames,
                                                                       [&]
                                                                         ( localIndex const regionIndex, auto & elementSubRegion )
    {
      if( elementSubRegion.template hasWrapper< string >( FlowSolverBase::viewKeyStruct::solidNamesString() ) )
      {
        poromechanicsRegions.insert( regionNames[regionIndex] );
      }
    } );

    string_array poromechanicsRegionNames;
    poromechanicsRegionNames.reserve( poromechanicsRegions.size());
    for( auto const & region : poromechanicsRegions )
    {
      poromechanicsRegionNames.emplace_back( region );
    }

    solidMechanicsConformingContactKernels::PressureFaceBubbleFactory kernelFactory( dispDofNumber,
                                                                                     bubbleDofNumber,
                                                                                     dofManager.rankOffset(),
                                                                                     localMatrix,
                                                                                     localRhs,
                                                                                     dt );


    real64 maxTraction = finiteElement::regionBasedKernelApplication
                         < parallelDevicePolicy< >,
                           PorousSolid< ElasticIsotropic, ConstantPermeability >,
                           CellElementSubRegion >( mesh,
                                                   poromechanicsRegionNames,
                                                   getDiscretizationName(),
                                                   FlowSolverBase::viewKeyStruct::solidNamesString(),
                                                   kernelFactory );

    GEOS_UNUSED_VAR( maxTraction );
  } );

}

void SolidMechanicsAugmentedLagrangianContact::implicitStepComplete( real64 const & time_n,
                                                                     real64 const & dt,
                                                                     DomainPartition & domain )
{

  SolidMechanicsLagrangianFEM::implicitStepComplete( time_n, dt, domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {
    //fractureRegion.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
    //{

    ElementRegionManager & elemManager = mesh.getElemManager();
    SurfaceElementRegion & region = elemManager.getRegion< SurfaceElementRegion >( getUniqueFractureRegionName() );
    FaceElementSubRegion & subRegion = region.getUniqueSubRegion< FaceElementSubRegion >();

    arrayView2d< real64 const > const dispJump = subRegion.getField< contact::dispJump >();
    arrayView2d< real64 > const oldDispJump = subRegion.getField< contact::oldDispJump >();
    arrayView2d< real64 > const deltaDispJump  = subRegion.getField< contact::deltaDispJump >();

    arrayView2d< real64 > const traction  = subRegion.getField< contact::traction >();

    arrayView1d< integer const > const fractureState = subRegion.getField< contact::fractureState >();
    arrayView1d< integer > const oldFractureState = subRegion.getField< contact::oldFractureState >();

    arrayView1d< real64 > const slip = subRegion.getField< contact::slip >();
    arrayView1d< real64 > const tangentialTraction  = subRegion.getField< contact::tangentialTraction >();

    forAll< parallelDevicePolicy<> >( subRegion.size(),
                                      [ = ]
                                      GEOS_DEVICE ( localIndex const kfe )
    {
      // Compute the slip
      real64 const shearDisp[2] = { dispJump[kfe][1],
                                    dispJump[kfe][2] };
      slip[kfe] = LvArray::tensorOps::l2Norm< 2 >( shearDisp );

      // Compute current Tau and limit Tau
      real64 const tau[2] = { traction[kfe][1],
                              traction[kfe][2] };
      tangentialTraction[kfe] = LvArray::tensorOps::l2Norm< 2 >( tau );

      LvArray::tensorOps::fill< 3 >( deltaDispJump[kfe], 0.0 );
      LvArray::tensorOps::copy< 3 >( oldDispJump[kfe], dispJump[kfe] );
      oldFractureState[kfe] = fractureState[kfe];

    } );
  } );
  // } );

}

real64 SolidMechanicsAugmentedLagrangianContact::calculateResidualNorm( real64 const & time,
                                                                        real64 const & dt,
                                                                        DomainPartition const & domain,
                                                                        DofManager const & dofManager,
                                                                        arrayView1d< real64 const > const & localRhs )
{

  GEOS_MARK_FUNCTION;

  real64 const solidResidualNorm = SolidMechanicsLagrangianFEM::calculateResidualNorm( time, dt, domain, dofManager, localRhs );

  string const bubbleDofKey = dofManager.getKey( contact::totalBubbleDisplacement::key() );

  globalIndex const rankOffset = dofManager.rankOffset();

  RAJA::ReduceSum< parallelDeviceReduce, real64 > localSum( 0.0 );

  // globalResidualNorm[0]: the sum of all the local sum(rhs^2).
  // globalResidualNorm[1]: max of max force of each rank. Basically max force globally
  real64 globalResidualNorm[2] = {0, 0};



  // Bubble residual
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                string_array const & )
  {

    FaceManager const & faceManager = mesh.getFaceManager();
    auto const & fGhost = faceManager.ghostRank();
    ElementRegionManager const & elemManager = mesh.getElemManager();

    SurfaceElementRegion const & region = elemManager.getRegion< SurfaceElementRegion >( getUniqueFractureRegionName() );
    FaceElementSubRegion const & subRegion = region.getUniqueSubRegion< FaceElementSubRegion >();

    arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

    arrayView2d< localIndex const > const elemsToFaces = subRegion.faceList().toViewConst();

    arrayView1d< globalIndex const > const bubbleDofNumber = faceManager.getReference< globalIndex_array >( bubbleDofKey );

    forAll< parallelDevicePolicy<> >( subRegion.size(),
                                      [ = ]
                                      GEOS_HOST_DEVICE ( localIndex const kfe )
    {

      if( ghostRank[kfe] < 0 )
      {
        for( int kk=0; kk<2; ++kk )
        {
          localIndex const k = elemsToFaces[kfe][kk];
          localIndex const localRow = LvArray::integerConversion< localIndex >( bubbleDofNumber[k] - rankOffset );
          if( fGhost[k] < 0 )
          {

            for( localIndex i = 0; i < 3; ++i )
            {
              localSum += localRhs[localRow + i] * localRhs[localRow + i];
            }
          }
        }
      }

    } );
    real64 const localResidualNorm[2] = { localSum.get(), SolidMechanicsLagrangianFEM::getMaxForce() };

    int const rank     = MpiWrapper::commRank( MPI_COMM_GEOS );
    int const numRanks = MpiWrapper::commSize( MPI_COMM_GEOS );
    array1d< real64 > globalValues( numRanks * 2 );

    // Everything is done on rank 0
    MpiWrapper::gather( localResidualNorm,
                        2,
                        globalValues.data(),
                        2,
                        0,
                        MPI_COMM_GEOS );

    if( rank==0 )
    {
      for( int r=0; r<numRanks; ++r )
      {
        // sum/max across all ranks
        globalResidualNorm[0] += globalValues[r*2];
        globalResidualNorm[1] = std::max( globalResidualNorm[1], globalValues[r*2+1] );
      }
    }

    MpiWrapper::bcast( globalResidualNorm, 2, 0, MPI_COMM_GEOS );

  } );

  real64 const bubbleResidualNorm = sqrt( globalResidualNorm[0] )/(globalResidualNorm[1]+1);  // the + 1 is for the first
  // time-step when maxForce = 0;

  GEOS_LOG_LEVEL_RANK_0_NLR( logInfo::ResidualNorm,
                             GEOS_FMT( "        ( RBubbleDisp ) = ( {:4.2e} )", bubbleResidualNorm ));


  real64 totalResidualNorm = sqrt( solidResidualNorm * solidResidualNorm + bubbleResidualNorm * bubbleResidualNorm );

  getConvergenceStats().setResidualValue( "RBubbleDisp", bubbleResidualNorm );

  return totalResidualNorm;
}

void SolidMechanicsAugmentedLagrangianContact::applySystemSolution( DofManager const & dofManager,
                                                                    arrayView1d< real64 const > const & localSolution,
                                                                    real64 const scalingFactor,
                                                                    real64 const dt,
                                                                    DomainPartition & domain )
{

  GEOS_MARK_FUNCTION;

  SolidMechanicsLagrangianFEM::applySystemSolution( dofManager,
                                                    localSolution,
                                                    scalingFactor,
                                                    dt,
                                                    domain );

  dofManager.addVectorToField( localSolution,
                               contact::totalBubbleDisplacement::key(),
                               contact::totalBubbleDisplacement::key(),
                               scalingFactor );

  dofManager.addVectorToField( localSolution,
                               contact::totalBubbleDisplacement::key(),
                               contact::incrementalBubbleDisplacement::key(),
                               scalingFactor );

  // Synchronize bubble displacements before computing displacement jump
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {
    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addFields( FieldLocation::Face,
                              { contact::incrementalBubbleDisplacement::key(),
                                contact::totalBubbleDisplacement::key() } );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                         mesh,
                                                         domain.getNeighbors(),
                                                         true );
  } );

  // Loop for updating the displacement jump
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshName,
                                                                MeshLevel & mesh,
                                                                string_array const & )

  {

    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    string const & dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );
    string const & bubbleDofKey = dofManager.getKey( contact::totalBubbleDisplacement::key() );

    arrayView1d< globalIndex const > const dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
    arrayView1d< globalIndex const > const bubbleDofNumber = faceManager.getReference< globalIndex_array >( bubbleDofKey );

    string const & fractureRegionName = getUniqueFractureRegionName();

    CRSMatrix< real64, globalIndex > const voidMatrix;
    array1d< real64 > const voidRhs;

    forFiniteElementOnFractureSubRegions( meshName, [&] ( string const &,
                                                          finiteElement::FiniteElementBase const & subRegionFE,
                                                          arrayView1d< localIndex const > const & faceElementList )
    {

      solidMechanicsConformingContactKernels::DispJumpUpdateFactory kernelFactory( dispDofNumber,
                                                                                   bubbleDofNumber,
                                                                                   dofManager.rankOffset(),
                                                                                   voidMatrix.toViewConstSizes(),
                                                                                   voidRhs.toView(),
                                                                                   dt,
                                                                                   faceElementList );

      real64 maxTraction = finiteElement::interfaceBasedKernelApplication< parallelDevicePolicy< >, NullModel >( mesh,
                                                                                                                 fractureRegionName,
                                                                                                                 faceElementList,
                                                                                                                 subRegionFE,
                                                                                                                 "",
                                                                                                                 kernelFactory );

      GEOS_UNUSED_VAR( maxTraction );

    } );
  } );

  // Synchronize displacement jump after computation
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {
    FieldIdentifiers fieldsToBeSync;

    fieldsToBeSync.addElementFields( { contact::dispJump::key(),
                                       contact::deltaDispJump::key() },
                                     { getUniqueFractureRegionName() } );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                         mesh,
                                                         domain.getNeighbors(),
                                                         true );
  } );

}

void SolidMechanicsAugmentedLagrangianContact::updateState( DomainPartition & domain )
{
  GEOS_UNUSED_VAR( domain );
}

bool SolidMechanicsAugmentedLagrangianContact::updateConfiguration( DomainPartition & domain,
                                                                    integer const GEOS_UNUSED_PARAM( configurationLoopIter ) )
{
  GEOS_MARK_FUNCTION;

  array1d< int > condConv;
  localIndex globalCondConv[5] = {0, 0, 0, 0, 0};

  array2d< real64 > traction_new;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                FaceElementSubRegion & subRegion )
    {

      string const & frictionLawName = subRegion.template getReference< string >( viewKeyStruct::frictionLawNameString() );
      FrictionBase const & frictionLaw = getConstitutiveModel< FrictionBase >( subRegion, frictionLawName );

      arrayView1d< integer const > const ghostRank = subRegion.ghostRank();
      arrayView2d< real64 > const traction = subRegion.getField< contact::traction >();
      arrayView2d< real64 const > const dispJump = subRegion.getField< contact::dispJump >();

      arrayView2d< real64 const > const deltaDispJump = subRegion.getField< contact::deltaDispJump >();
      arrayView2d< real64 const > const iterativePenalty = subRegion.getField< contact::iterativePenalty >();
      arrayView1d< integer const > const fractureState = subRegion.getField< contact::fractureState >();

      arrayView1d< real64 const > const normalDisplacementTolerance =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::normalDisplacementToleranceString() );
      arrayView1d< real64 > const & slidingTolerance =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::slidingToleranceString() );
      arrayView1d< real64 const > const & normalTractionTolerance =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::normalTractionToleranceString() );

      std::ptrdiff_t const sizes[ 2 ] = {subRegion.size(), 3};
      traction_new.resize( 2, sizes );
      arrayView2d< real64 > const traction_new_v = traction_new.toView();

      condConv.resize( subRegion.size());
      arrayView1d< int > const condConv_v = condConv.toView();

      // Update the traction field based on the displacement results from the nonlinear solve
      constitutiveUpdatePassThru( frictionLaw, [&] ( auto & castedFrictionLaw )
      {
        using FrictionType = TYPEOFREF( castedFrictionLaw );
        typename FrictionType::KernelWrapper frictionWrapper = castedFrictionLaw.createKernelUpdates();

        if( m_simultaneous )
        {
          solidMechanicsALMKernels::ComputeTractionSimultaneousKernel::
            launch< parallelDevicePolicy<> >( subRegion.size(),
                                              iterativePenalty,
                                              traction,
                                              dispJump,
                                              deltaDispJump,
                                              traction_new_v );
        }
        else
        {
          solidMechanicsALMKernels::ComputeTractionKernel::
            launch< parallelDevicePolicy<> >( subRegion.size(),
                                              frictionWrapper,
                                              iterativePenalty,
                                              traction,
                                              dispJump,
                                              deltaDispJump,
                                              traction_new_v );
        }
      } );

      real64 const slidingCheckTolerance = m_slidingCheckTolerance;

      constitutiveUpdatePassThru( frictionLaw, [&] ( auto & castedFrictionLaw )
      {
        using FrictionType = TYPEOFREF( castedFrictionLaw );
        typename FrictionType::KernelWrapper frictionWrapper = castedFrictionLaw.createKernelUpdates();

        solidMechanicsALMKernels::ConstraintCheckKernel::
          launch< parallelDevicePolicy<> >( subRegion.size(),
                                            frictionWrapper,
                                            ghostRank,
                                            traction,
                                            dispJump,
                                            deltaDispJump,
                                            normalTractionTolerance,
                                            normalDisplacementTolerance,
                                            slidingTolerance,
                                            slidingCheckTolerance,
                                            fractureState,
                                            condConv_v );
      } );

      RAJA::ReduceSum< parallelDeviceReduce, localIndex > localSum[5] =
      { RAJA::ReduceSum< parallelDeviceReduce, localIndex >( 0 ),
        RAJA::ReduceSum< parallelDeviceReduce, localIndex >( 0 ),
        RAJA::ReduceSum< parallelDeviceReduce, localIndex >( 0 ),
        RAJA::ReduceSum< parallelDeviceReduce, localIndex >( 0 ),
        RAJA::ReduceSum< parallelDeviceReduce, localIndex >( 0 ) };
      forAll< parallelDevicePolicy<> >( subRegion.size(), [ = ] GEOS_HOST_DEVICE ( localIndex const kfe )
      {
        if( ghostRank[kfe] < 0 )
        {
          localSum[condConv_v[kfe]] += 1;
        }
      } );

      localIndex const localConvCond[5] = { static_cast< localIndex >( localSum[0].get()),
                                            static_cast< localIndex >( localSum[1].get()),
                                            static_cast< localIndex >( localSum[2].get()),
                                            static_cast< localIndex >( localSum[3].get()),
                                            static_cast< localIndex >( localSum[4].get()) };

      int const rank     = MpiWrapper::commRank( MPI_COMM_GEOS );
      int const numRanks = MpiWrapper::commSize( MPI_COMM_GEOS );
      array1d< localIndex > globalValues( numRanks * 5 );

      // Everything is done on rank 0
      MpiWrapper::gather( localConvCond,
                          5,
                          globalValues.data(),
                          5,
                          0,
                          MPI_COMM_GEOS );

      if( rank==0 )
      {
        for( int r=0; r<numRanks; ++r )
        {
          // sum/max across all ranks
          for( int i=0; i<5; ++i )
          {
            globalCondConv[i] += globalValues[r*5+i];
          }
        }
      }

      MpiWrapper::bcast( globalCondConv, 5, 0, MPI_COMM_GEOS );

    } );
  } );

  localIndex totCondNotConv = 0;
  for( int i=0; i<4; ++i )
  {
    totCondNotConv+=globalCondConv[i+1];
  }

  int hasConfigurationConvergedGlobally = (totCondNotConv == 0) ? true : false;

  GEOS_LOG_LEVEL_RANK_0( logInfo::Convergence,
                         GEOS_FMT( "  ALM convergence summary:"
                                   " converged: {:6} | stick & gn>0: {:6} | interpenetration:  {:6} | stick & gt>lim:  {:6} | tau>tauLim:  {:6}\n",
                                   globalCondConv[0], globalCondConv[1], globalCondConv[2],
                                   globalCondConv[3], globalCondConv[4] ));

  if( hasConfigurationConvergedGlobally )
  {

    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                  MeshLevel & mesh,
                                                                  string_array const & regionNames )
    {
      ElementRegionManager & elemManager = mesh.getElemManager();

      elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                  FaceElementSubRegion & subRegion )
      {

        arrayView2d< real64 > const traction_new_v = traction_new.toView();
        arrayView2d< real64 > const traction = subRegion.getField< contact::traction >();

        forAll< parallelDevicePolicy<> >( subRegion.size(), [ = ] GEOS_HOST_DEVICE ( localIndex const kfe )
        {
          LvArray::tensorOps::copy< 3 >( traction[kfe], traction_new_v[kfe] );
        } );
      } );
    } );
  }
  else
  {
    forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                  MeshLevel & mesh,
                                                                  string_array const & regionNames )
    {
      ElementRegionManager & elemManager = mesh.getElemManager();

      elemManager.forElementSubRegions< FaceElementSubRegion >( regionNames, [m_symmetric=m_symmetric]( localIndex const,
                                                                                                        FaceElementSubRegion & subRegion )
      {

        string const & frictionLawName = subRegion.template getReference< string >( viewKeyStruct::frictionLawNameString() );
        FrictionBase const & frictionLaw = getConstitutiveModel< FrictionBase >( subRegion, frictionLawName );

        arrayView2d< real64 > const traction = subRegion.getField< contact::traction >();

        arrayView1d< real64 const > const normalTractionTolerance =
          subRegion.getReference< array1d< real64 > >( viewKeyStruct::normalTractionToleranceString() );

        arrayView2d< real64 > const iterativePenalty = subRegion.getField< contact::iterativePenalty >().toView();

        arrayView2d< real64 > const dispJumpUpdPenalty =
          subRegion.getReference< array2d< real64 > >( viewKeyStruct::dispJumpUpdPenaltyString() );

        arrayView1d< integer > const fractureState = subRegion.getField< contact::fractureState >();

        arrayView2d< real64 const > const dispJump = subRegion.getField< contact::dispJump >();

        arrayView2d< real64 const > const oldDispJump = subRegion.getField< contact::oldDispJump >();

        arrayView2d< real64 const > const deltaDispJump = subRegion.getField< contact::deltaDispJump >();

        constitutiveUpdatePassThru( frictionLaw, [&] ( auto & castedFrictionLaw )
        {
          using FrictionType = TYPEOFREF( castedFrictionLaw );
          typename FrictionType::KernelWrapper frictionWrapper = castedFrictionLaw.createKernelUpdates();

          solidMechanicsALMKernels::UpdateStateKernel::
            launch< parallelDevicePolicy<> >( subRegion.size(),
                                              frictionWrapper,
                                              oldDispJump,
                                              dispJump,
                                              iterativePenalty,
                                              m_symmetric,
                                              normalTractionTolerance,
                                              traction,
                                              fractureState );
        } );

        forAll< parallelDevicePolicy<> >( subRegion.size(), [ dispJumpUpdPenalty, dispJump, deltaDispJump ]
                                          GEOS_DEVICE ( localIndex const kfe )
        {
          dispJumpUpdPenalty[kfe][0] = dispJump[kfe][0];
          dispJumpUpdPenalty[kfe][1] = deltaDispJump[kfe][1];
          dispJumpUpdPenalty[kfe][2] = deltaDispJump[kfe][2];
        } );
      } );
    } );
  }

  // Need to synchronize the fracture state due to the use will be made of in AssemblyStabilization
  synchronizeFractureState( domain );

  // Update lists of stick and slip elements
  if( !hasConfigurationConvergedGlobally )
  {
    updateStickSlipList( domain );
  }

  return hasConfigurationConvergedGlobally;

}

void SolidMechanicsAugmentedLagrangianContact::updateStickSlipList( DomainPartition const & domain )
{

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshName,
                                                                MeshLevel const & mesh,
                                                                string_array const & )

  {

    ElementRegionManager const & elemManager = mesh.getElemManager();
    SurfaceElementRegion const & region = elemManager.getRegion< SurfaceElementRegion >( getUniqueFractureRegionName() );
    FaceElementSubRegion const & subRegion = region.getUniqueSubRegion< FaceElementSubRegion >();

    arrayView1d< integer const > const fractureState = subRegion.getField< contact::fractureState >();

    // Ensure the stick/slip maps contain an entry for this mesh even when
    // the fracture face type map is empty (no fracture elements).
    m_faceTypesToFaceElementsStick.get_inserted( meshName );
    m_faceTypesToFaceElementsSlip.get_inserted( meshName );

    forFiniteElementOnFractureSubRegions( meshName, [&] ( string const & finiteElementName,
                                                          finiteElement::FiniteElementBase const &,
                                                          arrayView1d< localIndex const > const & faceElementList )
    {

      array1d< localIndex > keys( subRegion.size());
      array1d< localIndex > vals( subRegion.size());
      array1d< localIndex > stickList;
      array1d< localIndex > slipList;
      RAJA::ReduceSum< ReducePolicy< parallelHostPolicy >, localIndex > nStick_r( 0 );
      RAJA::ReduceSum< ReducePolicy< parallelHostPolicy >, localIndex > nSlip_r( 0 );

      arrayView1d< localIndex > const keys_v = keys.toView();
      arrayView1d< localIndex > const vals_v = vals.toView();
      forAll< parallelHostPolicy >( faceElementList.size(),
                                    [ = ] GEOS_HOST ( localIndex const kfe )
      {
        localIndex const faceIndex = faceElementList[kfe];
        if( fractureState[faceIndex] == contact::FractureState::Stick )
        {
          keys_v[kfe]=0;
          vals_v[kfe]=faceIndex;
          nStick_r += 1;
        }
        else if(( fractureState[faceIndex] == contact::FractureState::Slip ) ||
                (fractureState[faceIndex] == contact::FractureState::NewSlip))
        {
          keys_v[kfe]=1;
          vals_v[kfe]=faceIndex;
          nSlip_r += 1;
        }
        else
        {
          keys_v[kfe] = 2;
        }

      } );

      localIndex nStick = static_cast< localIndex >(nStick_r.get());
      localIndex nSlip = static_cast< localIndex >(nSlip_r.get());

      // Sort vals according to keys to ensure that
      // elements of the same type are adjacent in the vals list.
      // This arrangement allows for efficient copying into the container
      // by leveraging parallelism.
      RAJA::sort_pairs< parallelHostPolicy >( keys_v, vals_v );

      stickList.resize( nStick );
      slipList.resize( nSlip );
      arrayView1d< localIndex > const stickList_v = stickList.toView();
      arrayView1d< localIndex > const slipList_v = slipList.toView();

      forAll< parallelHostPolicy >( nStick, [ = ] GEOS_HOST ( localIndex const kfe )
      {
        stickList_v[kfe] = vals_v[kfe];
      } );

      forAll< parallelHostPolicy >( nSlip, [ = ] GEOS_HOST ( localIndex const kfe )
      {
        slipList_v[kfe] = vals_v[nStick+kfe];
      } );

      // these operations are suspect for host/device issues.
      // the values of m_faceTypesToFaceElementsStick/Slip are host
      // but the values being inserted may not be.
      m_faceTypesToFaceElementsStick.get_inserted( meshName ).get_inserted( finiteElementName ) = stickList;
      m_faceTypesToFaceElementsSlip.get_inserted( meshName ).get_inserted( finiteElementName ) = slipList;

      GEOS_LOG_LEVEL_RANK_0( logInfo::ConfigurationStatistics, GEOS_FMT( "# stick elements: {}, # slip elements: {}", nStick, nSlip ))
    } );
  } );

}

void SolidMechanicsAugmentedLagrangianContact::createFaceTypeList( DomainPartition const & domain )
{

  // Generate lists containing elements of various face types
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshName,
                                                                MeshLevel const & mesh,
                                                                string_array const & )
  {
    FaceManager const & faceManager = mesh.getFaceManager();
    ElementRegionManager const & elemManager = mesh.getElemManager();
    ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

    SurfaceElementRegion const & region = elemManager.getRegion< SurfaceElementRegion >( getUniqueFractureRegionName() );
    FaceElementSubRegion const & subRegion = region.getUniqueSubRegion< FaceElementSubRegion >();

    // When there are no fracture face elements, insert empty lists so that
    // downstream .at( meshName ) lookups do not throw.
    if( subRegion.size() == 0 )
    {
      m_faceTypesToFaceElements.get_inserted( meshName );
      return;
    }

    arrayView2d< localIndex const > const elemsToFaces = subRegion.faceList().toViewConst();

    array1d< localIndex > keys( subRegion.size());
    array1d< localIndex > vals( subRegion.size());
    array1d< localIndex > quadList;
    array1d< localIndex > triList;
    RAJA::ReduceSum< ReducePolicy< parallelHostPolicy >, localIndex > nTri_r( 0 );
    RAJA::ReduceSum< ReducePolicy< parallelHostPolicy >, localIndex > nQuad_r( 0 );

    arrayView1d< localIndex > const keys_v = keys.toView();
    arrayView1d< localIndex > const vals_v = vals.toView();
    // Determine the size of the lists and generate the vector keys and vals for parallel indexing into lists.
    // (With RAJA, parallelizing this operation seems the most viable approach.)
    forAll< parallelHostPolicy >( subRegion.size(),
                                  [ = ] GEOS_HOST ( localIndex const kfe )
    {
      localIndex const kf0 = elemsToFaces[kfe][0];
      localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kf0 );
      if( numNodesPerFace == 3 )
      {
        keys_v[kfe]=0;
        vals_v[kfe]=kfe;
        nTri_r += 1;
      }
      else if( numNodesPerFace == 4 )
      {
        keys_v[kfe]=1;
        vals_v[kfe]=kfe;
        nQuad_r += 1;
      }
      else
      {
        GEOS_ERROR( "SolidMechanicsAugmentedLagrangianContact:: invalid face type", getDataContext() );
      }
    } );

    localIndex nQuad = static_cast< localIndex >(nQuad_r.get());
    localIndex nTri = static_cast< localIndex >(nTri_r.get());

    // Sort vals according to keys to ensure that
    // elements of the same type are adjacent in the vals list.
    // This arrangement allows for efficient copying into the container
    // by leveraging parallelism.
    RAJA::sort_pairs< parallelHostPolicy >( keys_v, vals_v );

    quadList.resize( nQuad );
    triList.resize( nTri );
    arrayView1d< localIndex > const quadList_v = quadList.toView();
    arrayView1d< localIndex > const triList_v = triList.toView();

    forAll< parallelHostPolicy >( nTri, [ = ] GEOS_HOST ( localIndex const kfe )
    {
      triList_v[kfe] = vals_v[kfe];
    } );

    forAll< parallelHostPolicy >( nQuad, [ = ] GEOS_HOST ( localIndex const kfe )
    {
      quadList_v[kfe] = vals_v[nTri+kfe];
    } );

    // these operations are suspect for host/device issues.
    m_faceTypesToFaceElements.get_inserted( meshName ).get_inserted( "Quadrilateral" ) = quadList;
    m_faceTypesToFaceElements.get_inserted( meshName ).get_inserted( "Triangle" ) = triList;
  } );

}

void SolidMechanicsAugmentedLagrangianContact::createBubbleCellList( DomainPartition & domain ) const
{

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    SurfaceElementRegion const & region = elemManager.getRegion< SurfaceElementRegion >( getUniqueFractureRegionName() );
    FaceElementSubRegion const & subRegion = region.getUniqueSubRegion< FaceElementSubRegion >();

    // Nothing to do when there are no fracture face elements.
    if( subRegion.size() == 0 )
    {
      return;
    }

    // Array to store face indexes
    array1d< localIndex > tmpSpace( 2*subRegion.size());
    SortedArray< localIndex > faceIdList;

    arrayView1d< localIndex > const tmpSpace_v = tmpSpace.toView();
    // Store indexes of faces in the temporany array.
    {
      arrayView2d< localIndex const > const elemsToFaces = subRegion.faceList().toViewConst();

      forAll< parallelHostPolicy >( subRegion.size(), [ = ] GEOS_HOST ( localIndex const kfe )
      {

        localIndex const kf0 = elemsToFaces[kfe][0], kf1 = elemsToFaces[kfe][1];
        tmpSpace_v[2*kfe] = kf0, tmpSpace_v[2*kfe+1] = kf1;

      } );
    }

    // Sort indexes to enable efficient searching using binary search.
    RAJA::stable_sort< parallelHostPolicy >( tmpSpace_v );
    // Move data back to host: after the device sort the buffer lives on the GPU,
    // but SortedArray::insert is a host-only operation that dereferences the iterators.
    faceIdList.insert( tmpSpace_v.begin(), tmpSpace_v.end());
    SortedArrayView< localIndex const > const faceIdList_v = faceIdList.toViewConst();

    // Search for bubble element on each CellElementSubRegion and
    // store element indexes, global and local face indexes.
    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames,
                                                              [&]( localIndex const,
                                                                   CellElementSubRegion & cellElementSubRegion )
    {

      arrayView2d< localIndex const > const elemsToFaces = cellElementSubRegion.faceList().toViewConst();

      RAJA::ReduceSum< ReducePolicy< parallelHostPolicy >, localIndex > nBubElems_r( 0 );

      localIndex const n_max = cellElementSubRegion.size() * elemsToFaces.size( 1 );
      array1d< localIndex > keys( n_max );
      array1d< localIndex > perms( n_max );
      array1d< localIndex > vals( n_max );
      array1d< localIndex > localFaceIds( n_max );

      arrayView1d< localIndex > const keys_v = keys.toView();
      arrayView1d< localIndex > const perms_v = perms.toView();
      arrayView1d< localIndex > const vals_v = vals.toView();
      arrayView1d< localIndex > const localFaceIds_v = localFaceIds.toView();

      forAll< parallelHostPolicy >( cellElementSubRegion.size(),
                                    [ = ] GEOS_HOST ( localIndex const kfe )
      {
        for( int i=0; i < elemsToFaces.size( 1 ); ++i )
        {
          perms_v[kfe*elemsToFaces.size( 1 )+i] = kfe*elemsToFaces.size( 1 )+i;
          if( faceIdList_v.contains( elemsToFaces[kfe][i] ))
          {
            keys_v[kfe*elemsToFaces.size( 1 )+i] = 0;
            vals_v[kfe*elemsToFaces.size( 1 )+i] = kfe;
            localFaceIds_v[kfe*elemsToFaces.size( 1 )+i] = i;
            nBubElems_r += 1;
          }
          else
          {
            keys_v[kfe*elemsToFaces.size( 1 )+i] = 1;
            vals_v[kfe*elemsToFaces.size( 1 )+i] = -1;
            localFaceIds_v[kfe*elemsToFaces.size( 1 )+i] = -1;
          }
        }
      } );

      // Sort perms according to keys to ensure that bubble elements are adjacent
      // and occupy the first positions of the list.
      // This arrangement allows for efficient copying into the container
      // by leveraging parallelism.
      localIndex nBubElems = static_cast< localIndex >(nBubElems_r.get());
      RAJA::sort_pairs< parallelHostPolicy >( keys_v, perms_v );

      array1d< localIndex > bubbleElemsList;
      bubbleElemsList.resize( nBubElems );

      arrayView1d< localIndex > const bubbleElemsList_v = bubbleElemsList.toView();

      forAll< parallelHostPolicy >( n_max, [ = ] GEOS_HOST ( localIndex const k )
      {
        keys_v[k] = vals_v[perms_v[k]];
      } );

      forAll< parallelHostPolicy >( nBubElems, [ = ] GEOS_HOST ( localIndex const k )
      {
        bubbleElemsList_v[k] = keys_v[k];
      } );

      // Get reference to the persistent storage and copy data to avoid dangling pointers
      array1d< localIndex > & bubbleCellsStorage =
        cellElementSubRegion.getReference< array1d< localIndex > >( CellElementSubRegion::viewKeyStruct::bubbleCellsString() );
      bubbleCellsStorage.resize( nBubElems );
      arrayView1d< localIndex > const bubbleCellsStorage_v = bubbleCellsStorage.toView();
      forAll< parallelHostPolicy >( nBubElems, [ = ] GEOS_HOST ( localIndex const k )
      {
        bubbleCellsStorage_v[k] = bubbleElemsList_v[k];
      } );

      forAll< parallelHostPolicy >( n_max, [ = ] GEOS_HOST ( localIndex const k )
      {
        keys_v[k] = localFaceIds_v[perms_v[k]];
      } );

      array2d< localIndex > faceElemsList;
      faceElemsList.resize( nBubElems, 2 );

      arrayView2d< localIndex > const faceElemsList_v = faceElemsList.toView();

      forAll< parallelHostPolicy >( nBubElems,
                                    [ = ] GEOS_HOST ( localIndex const k )
      {
        localIndex const kfe =  bubbleElemsList_v[k];
        faceElemsList_v[k][0] = elemsToFaces[kfe][keys_v[k]];
        faceElemsList_v[k][1] = keys_v[k];
      } );

      // Get reference to the persistent storage and copy data to avoid dangling pointers
      array2d< localIndex > & faceElemsStorage =
        cellElementSubRegion.getReference< array2d< localIndex > >( CellElementSubRegion::viewKeyStruct::toFaceElementsString() );
      faceElemsStorage.resize( nBubElems, 2 );
      arrayView2d< localIndex > const faceElemsStorage_v = faceElemsStorage.toView();
      forAll< parallelHostPolicy >( nBubElems, [ = ] GEOS_HOST ( localIndex const k )
      {
        faceElemsStorage_v[k][0] = faceElemsList_v[k][0];
        faceElemsStorage_v[k][1] = faceElemsList_v[k][1];
      } );

    } );

  } );

}

void SolidMechanicsAugmentedLagrangianContact::addCouplingNumNonzeros( DomainPartition & domain,
                                                                       DofManager & dofManager,
                                                                       arrayView1d< localIndex > const & rowLengths ) const
{

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                string_array const & regionNames )
  {

    ElementRegionManager const & elemManager = mesh.getElemManager();
    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

    globalIndex const rankOffset = dofManager.rankOffset();

    string const bubbleDofKey = dofManager.getKey( contact::totalBubbleDisplacement::key() );
    string const dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );

    arrayView1d< globalIndex const > const bubbleDofNumber = faceManager.getReference< globalIndex_array >( bubbleDofKey );
    arrayView1d< globalIndex const > const dispDofNumber =  nodeManager.getReference< globalIndex_array >( dispDofKey );

    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const, CellElementSubRegion const & cellElementSubRegion )
    {

      arrayView1d< localIndex const > const bubbleElemsList = cellElementSubRegion.bubbleElementsList();
      arrayView2d< localIndex const > const faceElemsList = cellElementSubRegion.faceElementsList();

      localIndex const numDispDof = 3*cellElementSubRegion.numNodesPerElement();

      for( localIndex bi=0; bi<bubbleElemsList.size(); ++bi )
      {
        localIndex const cellIndex = bubbleElemsList[bi];
        localIndex const k = faceElemsList[bi][0];

        localIndex const localRow = LvArray::integerConversion< localIndex >( bubbleDofNumber[k] - rankOffset );

        if( localRow >= 0 && localRow < rowLengths.size() )
        {
          for( localIndex i=0; i<3; ++i )
          {
            rowLengths[localRow + i] += numDispDof;
          }
        }

        for( localIndex a=0; a<cellElementSubRegion.numNodesPerElement(); ++a )
        {
          const localIndex & node = cellElementSubRegion.nodeList( cellIndex, a );
          localIndex const localDispRow = LvArray::integerConversion< localIndex >( dispDofNumber[node] - rankOffset );

          if( localDispRow >= 0 && localDispRow < rowLengths.size() )
          {
            for( int d=0; d<3; ++d )
            {
              rowLengths[localDispRow + d] += 3;
            }
          }
        }
      }

    } );

    SurfaceElementRegion const & region = elemManager.getRegion< SurfaceElementRegion >( getUniqueFractureRegionName() );
    FaceElementSubRegion const & subRegion = region.getUniqueSubRegion< FaceElementSubRegion >();
    arrayView2d< localIndex const > const elemsToFaces = subRegion.faceList().toViewConst();

    for( localIndex kfe=0; kfe<subRegion.size(); ++kfe )
    {
      localIndex const kf0 = elemsToFaces[kfe][0];
      localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kf0 );
      localIndex const numDispDof = 3*numNodesPerFace;

      for( int k=0; k<2; ++k )
      {
        localIndex const kf = elemsToFaces[kfe][k];

        localIndex const localRow = LvArray::integerConversion< localIndex >( bubbleDofNumber[kf] - rankOffset );

        if( localRow >= 0 && localRow < rowLengths.size() )
        {
          for( localIndex i=0; i<3; ++i )
          {
            rowLengths[localRow + i] += numDispDof;
          }
        }

        for( localIndex a=0; a<numNodesPerFace; ++a )
        {
          const localIndex & node = faceToNodeMap( kf, a );
          localIndex const localDispRow = LvArray::integerConversion< localIndex >( dispDofNumber[node] - rankOffset );

          if( localDispRow >= 0 && localDispRow < rowLengths.size() )
          {
            for( int d=0; d<3; ++d )
            {
              rowLengths[localDispRow + d] += 3;
            }
          }
        }
      }

    }

  } );
}

void SolidMechanicsAugmentedLagrangianContact::addCouplingSparsityPattern( DomainPartition const & domain,
                                                                           DofManager const & dofManager,
                                                                           SparsityPatternView< globalIndex > const & pattern ) const
{

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                string_array const & regionNames )
  {

    ElementRegionManager const & elemManager = mesh.getElemManager();
    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    globalIndex const rankOffset = dofManager.rankOffset();

    string const bubbleDofKey = dofManager.getKey( contact::totalBubbleDisplacement::key() );
    string const dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );

    arrayView1d< globalIndex const > const bubbleDofNumber = faceManager.getReference< globalIndex_array >( bubbleDofKey );
    arrayView1d< globalIndex const > const dispDofNumber =  nodeManager.getReference< globalIndex_array >( dispDofKey );

    static constexpr int maxNumDispDof = 3 * 8;

    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const, CellElementSubRegion const & cellElementSubRegion )
    {

      arrayView1d< localIndex const > const bubbleElemsList = cellElementSubRegion.bubbleElementsList();
      arrayView2d< localIndex const > const faceElemsList = cellElementSubRegion.faceElementsList();

      localIndex const numDispDof = 3*cellElementSubRegion.numNodesPerElement();

      for( localIndex bi=0; bi<bubbleElemsList.size(); ++bi )
      {
        localIndex const cellIndex = bubbleElemsList[bi];
        localIndex const k = faceElemsList[bi][0];

        // working arrays
        stackArray1d< globalIndex, maxNumDispDof > eqnRowIndicesDisp ( numDispDof );
        stackArray1d< globalIndex, 3 > eqnRowIndicesBubble( 3 );
        stackArray1d< globalIndex, maxNumDispDof > dofColIndicesDisp ( numDispDof );
        stackArray1d< globalIndex, 3 > dofColIndicesBubble( 3 );

        for( localIndex idof = 0; idof < 3; ++idof )
        {
          eqnRowIndicesBubble[idof] = bubbleDofNumber[k] + idof - rankOffset;
          dofColIndicesBubble[idof] = bubbleDofNumber[k] + idof;
        }

        for( localIndex a=0; a<cellElementSubRegion.numNodesPerElement(); ++a )
        {
          const localIndex & node = cellElementSubRegion.nodeList( cellIndex, a );
          for( localIndex idof = 0; idof < 3; ++idof )
          {
            eqnRowIndicesDisp[3*a + idof] = dispDofNumber[node] + idof - rankOffset;
            dofColIndicesDisp[3*a + idof] = dispDofNumber[node] + idof;
          }
        }

        for( localIndex i = 0; i < eqnRowIndicesDisp.size(); ++i )
        {
          if( eqnRowIndicesDisp[i] >= 0 && eqnRowIndicesDisp[i] < pattern.numRows() )
          {
            for( localIndex j = 0; j < dofColIndicesBubble.size(); ++j )
            {
              pattern.insertNonZero( eqnRowIndicesDisp[i], dofColIndicesBubble[j] );
            }
          }
        }

        for( localIndex i = 0; i < eqnRowIndicesBubble.size(); ++i )
        {
          if( eqnRowIndicesBubble[i] >= 0 && eqnRowIndicesBubble[i] < pattern.numRows() )
          {
            for( localIndex j=0; j < dofColIndicesDisp.size(); ++j )
            {
              pattern.insertNonZero( eqnRowIndicesBubble[i], dofColIndicesDisp[j] );
            }
          }
        }

      }

    } );

    SurfaceElementRegion const & region = elemManager.getRegion< SurfaceElementRegion >( getUniqueFractureRegionName() );
    FaceElementSubRegion const & subRegion = region.getUniqueSubRegion< FaceElementSubRegion >();
    arrayView2d< localIndex const > const elemsToFaces = subRegion.faceList().toViewConst();
    ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

    static constexpr int maxNumDispFaceDof = 3 * 4;

    for( localIndex kfe=0; kfe<subRegion.size(); ++kfe )
    {

      localIndex const kf0 = elemsToFaces[kfe][0];
      localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kf0 );
      localIndex const numDispDof = 3*numNodesPerFace;

      for( int k=0; k<2; ++k )
      {
        localIndex const kf = elemsToFaces[kfe][k];
        localIndex const kf_other = elemsToFaces[kfe][(1+k)%2];

        // working arrays
        stackArray1d< globalIndex, maxNumDispFaceDof > eqnRowIndicesDisp ( numDispDof );
        stackArray1d< globalIndex, 3 > eqnRowIndicesBubble( 3 );
        stackArray1d< globalIndex, maxNumDispFaceDof > dofColIndicesDisp ( numDispDof );
        stackArray1d< globalIndex, 3 > dofColIndicesBubble( 3 );

        for( localIndex idof = 0; idof < 3; ++idof )
        {
          eqnRowIndicesBubble[idof] = bubbleDofNumber[kf] + idof - rankOffset;
          dofColIndicesBubble[idof] = bubbleDofNumber[kf] + idof;
        }

        for( localIndex a=0; a<numNodesPerFace; ++a )
        {
          const localIndex & node = faceToNodeMap( kf_other, a );
          for( localIndex idof = 0; idof < 3; ++idof )
          {
            eqnRowIndicesDisp[3*a + idof] = dispDofNumber[node] + idof - rankOffset;
            dofColIndicesDisp[3*a + idof] = dispDofNumber[node] + idof;
          }
        }

        for( localIndex i = 0; i < eqnRowIndicesDisp.size(); ++i )
        {
          if( eqnRowIndicesDisp[i] >= 0 && eqnRowIndicesDisp[i] < pattern.numRows() )
          {
            for( localIndex j = 0; j < dofColIndicesBubble.size(); ++j )
            {
              pattern.insertNonZero( eqnRowIndicesDisp[i], dofColIndicesBubble[j] );
            }
          }
        }

        for( localIndex i = 0; i < eqnRowIndicesBubble.size(); ++i )
        {
          if( eqnRowIndicesBubble[i] >= 0 && eqnRowIndicesBubble[i] < pattern.numRows() )
          {
            for( localIndex j=0; j < dofColIndicesDisp.size(); ++j )
            {
              pattern.insertNonZero( eqnRowIndicesBubble[i], dofColIndicesDisp[j] );
            }
          }
        }

      }
    }
  } );

}

void SolidMechanicsAugmentedLagrangianContact::computeTolerances( DomainPartition & domain ) const
{
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {
    FaceManager const & faceManager = mesh.getFaceManager();
    ElementRegionManager & elemManager = mesh.getElemManager();

    // Get the "face to element" map (valid for the entire mesh)
    FaceManager::ElemMapType const & faceToElem = faceManager.toElementRelation();
    arrayView2d< localIndex const > const faceToElemRegion = faceToElem.m_toElementRegion;
    arrayView2d< localIndex const > const faceToElemSubRegion = faceToElem.m_toElementSubRegion;
    arrayView2d< localIndex const > const faceToElemIndex = faceToElem.m_toElementIndex;

    // Get the volume for all elements
    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const elemVolume =
      elemManager.constructViewAccessor< array1d< real64 >, arrayView1d< real64 const > >( ElementSubRegionBase::viewKeyStruct::elementVolumeString() );

    // Bulk modulus accessor
    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const bulkModulus =
      elemManager.constructMaterialViewAccessor< ElasticIsotropic, array1d< real64 >, arrayView1d< real64 const > >( fields::solid::bulkModulus::key() );
    // Shear modulus accessor
    ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const shearModulus =
      elemManager.constructMaterialViewAccessor< ElasticIsotropic, array1d< real64 >, arrayView1d< real64 const > >( fields::solid::shearModulus::key() );

    using NodeMapViewType = arrayView2d< localIndex const, cells::NODE_MAP_USD >;
    ElementRegionManager::ElementViewAccessor< NodeMapViewType > const elemToNode =
      elemManager.constructViewAccessor< CellElementSubRegion::NodeMapType, NodeMapViewType >( ElementSubRegionBase::viewKeyStruct::nodeListString() );

    elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
    {
      if( subRegion.hasField< contact::traction >() )
      {
        arrayView1d< real64 const > const faceArea = subRegion.getElementArea().toViewConst();
        arrayView3d< real64 const > const faceRotationMatrix = subRegion.getField< contact::rotationMatrix >().toView();
        arrayView2d< localIndex const > const elemsToFaces = subRegion.faceList().toViewConst();

        arrayView1d< real64 > const normalTractionTolerance =
          subRegion.getReference< array1d< real64 > >( viewKeyStruct::normalTractionToleranceString() );
        arrayView1d< real64 > const normalDisplacementTolerance =
          subRegion.getReference< array1d< real64 > >( viewKeyStruct::normalDisplacementToleranceString() );
        arrayView1d< real64 > const slidingTolerance =
          subRegion.getReference< array1d< real64 > >( viewKeyStruct::slidingToleranceString() );

        arrayView2d< real64 > const
        iterativePenalty = subRegion.getField< contact::iterativePenalty >().toView();

        arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

        forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
        {

          if( ghostRank[kfe] < 0 )
          {
            real64 const area = faceArea[kfe];
            // approximation of the stiffness along coordinate directions
            // ( first, second ) index -> ( element index, direction )
            // 1. T -> top (index 0), B -> bottom (index 1)
            // 2. the coordinate direction (x, y, z)
            real64 stiffDiagApprox[ 2 ][ 3 ];
            real64 averageYoungModulus = 0.0;
            real64 averageConstrainedModulus = 0.0;
            real64 averageCharLength = 0.0;

            for( localIndex i = 0; i < 2; ++i )
            {
              localIndex const faceIndex = elemsToFaces[kfe][i];
              localIndex const er = faceToElemRegion[faceIndex][0];
              localIndex const esr = faceToElemSubRegion[faceIndex][0];
              localIndex const ei = faceToElemIndex[faceIndex][0];

              real64 const volume = elemVolume[er][esr][ei];

              // Get linear elastic isotropic constitutive parameters for the element
              real64 const K = bulkModulus[er][esr][ei];
              real64 const G = shearModulus[er][esr][ei];
              real64 const E = 9.0 * K * G / ( 3.0 * K + G );
              real64 const nu = ( 3.0 * K - 2.0 * G ) / ( 2.0 * ( 3.0 * K + G ) );
              real64 const M = K + 4.0 / 3.0 * G;

              real64 const charLength = pow( volume, 1.0 / 3.0 );

              // Combine E and nu to obtain a stiffness approximation (like it was an hexahedron)
              for( localIndex j = 0; j < 3; ++j )
              {
                stiffDiagApprox[ i ][ j ] = E / ( ( 1.0 + nu )*( 1.0 - 2.0*nu ) ) * 4.0 / 9.0 * ( 2.0 - 3.0 * nu ) * charLength;
              }

              averageYoungModulus += 0.5*E;
              averageConstrainedModulus += 0.5*M;
              averageCharLength += 0.5*charLength;
            }

            // Average the stiffness and compute the inverse
            real64 invStiffApprox[ 3 ][ 3 ] = { { 0 } };
            for( localIndex j = 0; j < 3; ++j )
            {
              invStiffApprox[ j ][ j ] = ( stiffDiagApprox[ 0 ][ j ] + stiffDiagApprox[ 1 ][ j ] ) / ( stiffDiagApprox[ 0 ][ j ] * stiffDiagApprox[ 1 ][ j ] );
            }

            // Rotate in the local reference system, computing R^T * (invK) * R
            real64 temp[ 3 ][ 3 ];
            LvArray::tensorOps::Rij_eq_AkiBkj< 3, 3, 3 >( temp, faceRotationMatrix[ kfe ], invStiffApprox );
            real64 rotatedInvStiffApprox[ 3 ][ 3 ];
            LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( rotatedInvStiffApprox, temp, faceRotationMatrix[ kfe ] );
            LvArray::tensorOps::scale< 3, 3 >( rotatedInvStiffApprox, area );

            // Finally, compute tolerances for the given fracture element
            normalDisplacementTolerance[kfe] = rotatedInvStiffApprox[ 0 ][ 0 ] * averageYoungModulus * m_tolJumpDispNFac;
            slidingTolerance[kfe] = sqrt( pow( rotatedInvStiffApprox[ 1 ][ 1 ], 2 ) +
                                          pow( rotatedInvStiffApprox[ 2 ][ 2 ], 2 )) * averageYoungModulus * m_tolJumpDispTFac;
            normalTractionTolerance[kfe] = m_tolNormalTracFac * (averageConstrainedModulus / averageCharLength) *
                                           (normalDisplacementTolerance[kfe]);

            iterativePenalty[kfe][0] = m_iterPenaltyNFac*averageConstrainedModulus/(averageCharLength);
            iterativePenalty[kfe][1] = m_iterPenaltyTFac*averageConstrainedModulus/(averageCharLength);
          }
        } );
      }
    } );
  } );
}

void SolidMechanicsAugmentedLagrangianContact::initializeTractionFromAdjacentCellStress( DomainPartition & domain ) const
{
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {
    FaceManager & faceManager = mesh.getFaceManager();
    ElementRegionManager & elemManager = mesh.getElemManager();

    // Get the "face to element" map (valid for the entire mesh)
    FaceManager::ElemMapType const & faceToElem = faceManager.toElementRelation();
    arrayView2d< localIndex const > const faceToElemRegion = faceToElem.m_toElementRegion;
    arrayView2d< localIndex const > const faceToElemSubRegion = faceToElem.m_toElementSubRegion;
    arrayView2d< localIndex const > const faceToElemIndex = faceToElem.m_toElementIndex;

    arrayView2d< real64 > totalBubbleDisplacement = faceManager.getField< contact::totalBubbleDisplacement >();

    // Get stress accessor
    ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const, cells::RANK2_TENSOR_USD > > const
    avgElementStress = elemManager.constructViewAccessor< array2d< real64, cells::RANK2_TENSOR_PERM >,
                                                          arrayView2d< real64 const, cells::RANK2_TENSOR_USD > >( solidMechanics::averageStress::key() );

    ElementRegionManager::ElementViewConst< arrayView2d< real64 const, cells::RANK2_TENSOR_USD > > const
    avgElementStressView = avgElementStress.toNestedViewConst();

    elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
    {
      if( subRegion.hasField< contact::traction >() )
      {
        arrayView3d< real64 const > const faceRotationMatrix = subRegion.getField< contact::rotationMatrix >().toViewConst();
        arrayView2d< localIndex const > const elemsToFaces = subRegion.faceList().toViewConst();
        arrayView2d< real64 > const traction = subRegion.getField< contact::traction >();
        arrayView2d< real64 > const dispJump = subRegion.getField< contact::dispJump >();
        arrayView2d< real64 const > const iterativePenalty = subRegion.getField< contact::iterativePenalty >().toViewConst();
        arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

        // Get friction law parameters for Coulomb check
        string const & frictionLawName = subRegion.template getReference< string >( viewKeyStruct::frictionLawNameString() );
        FrictionBase const & frictionLaw = getConstitutiveModel< FrictionBase >( subRegion, frictionLawName );

        // Try to get Coulomb parameters if available
        bool const hasCoulombParams = frictionLaw.hasWrapper( CoulombFriction::viewKeyStruct::cohesionString() ) &&
                                      frictionLaw.hasWrapper( CoulombFriction::viewKeyStruct::frictionCoefficientString() );

        real64 const cohesion = hasCoulombParams ? frictionLaw.getReference< real64 >( CoulombFriction::viewKeyStruct::cohesionString() ) : 0.0;
        real64 const frictionCoefficient = hasCoulombParams ? frictionLaw.getReference< real64 >( CoulombFriction::viewKeyStruct::frictionCoefficientString() ) : 0.0;

        forAll< parallelHostPolicy >( subRegion.size(), [ ghostRank,
                                                          faceRotationMatrix,
                                                          elemsToFaces,
                                                          faceToElemRegion,
                                                          faceToElemSubRegion,
                                                          faceToElemIndex,
                                                          traction,
                                                          dispJump,
                                                          iterativePenalty,
                                                          totalBubbleDisplacement,
                                                          hasCoulombParams,
                                                          cohesion,
                                                          frictionCoefficient,
                                                          avgElementStressView] ( localIndex const kfe )
        {
          if( ghostRank[kfe] < 0 )
          {
            // Each fracture element has two adjacent faces (one on each side of the fracture)
            // We compute traction from each side and average them
            real64 tGlobalAvg[3]{};
            real64 avgSigma[6]{};
            integer numValidCells{};

            // Use the fracture element's normal from the rotation matrix (first column)
            real64 const n[3]{ faceRotationMatrix( kfe, 0, 0 ), faceRotationMatrix( kfe, 1, 0 ), faceRotationMatrix( kfe, 2, 0 ) };

            // Loop over both faces of the fracture element
            for( localIndex faceIdx = 0; faceIdx < 2; ++faceIdx )
            {
              localIndex const faceIndex = elemsToFaces( kfe, faceIdx );

              // Get the adjacent volume element for this face
              localIndex const er = faceToElemRegion( faceIndex, 0 );
              localIndex const esr = faceToElemSubRegion( faceIndex, 0 );
              localIndex const ei = faceToElemIndex( faceIndex, 0 );

              // Skip if no valid stress data available
              if( avgElementStressView[er].empty() || avgElementStressView[er][esr].empty() )
              {
                continue;
              }

              // Check principal stress components in 3D elements adjecent to the fracture element.
              // If the maximun principal stress is tensile, issue an error. Principal stresses are
              // sorted in ascending order
              real64 principalStresses[3];
              LvArray::tensorOps::symEigenvalues< 3 >( principalStresses, avgElementStressView[er][esr][ei] );
              GEOS_ERROR_IF( principalStresses[2] > 0.0,
                             GEOS_FMT(
                               "ERROR: Maximum principal stress is tensile in element adjacent to fracture element {} "
                               "connected to face {}. Principal stresses (sorted): ({:.6e}, {:.6e}, {:.6e})",
                               kfe, faceIdx, principalStresses[0], principalStresses[1], principalStresses[2] ) );

              // Accumulate stress for averaging (for warning message)
              LvArray::tensorOps::add< 6 >( avgSigma, avgElementStressView[er][esr][ei] );

              // Compute and accumulate traction in global coordinates: t = sigma * n
              LvArray::tensorOps::Ri_add_symAijBj< 3 >( tGlobalAvg, avgElementStressView[er][esr][ei], n );

              ++numValidCells;
            }

            // Average the tractions if we have valid data from both sides
            if( numValidCells > 0 )
            {
              real64 const invNumCells = 1.0 / static_cast< real64 >( numValidCells );
              LvArray::tensorOps::scale< 3 >( tGlobalAvg, invNumCells );
              LvArray::tensorOps::scale< 6 >( avgSigma, invNumCells );

              // Rotate averaged traction to local coordinate system: t_local = R^T * t_global
              real64 tLocal[3];
              LvArray::tensorOps::Ri_eq_AjiBj< 3, 3 >( tLocal, faceRotationMatrix[kfe], tGlobalAvg );

              // Store the traction
              LvArray::tensorOps::copy< 3 >( traction[kfe], tLocal );

              // Update displacement jump based on penalty stiffness
              real64 const kn = iterativePenalty[kfe][0];
              real64 const kt = iterativePenalty[kfe][1];
              real64 dJump[3] = { 0.0, 0.0, 0.0 };

              if( kn > LvArray::NumericLimits< real64 >::min )
              {
                dJump[0] = tLocal[0] / kn;
              }
              if( kt > LvArray::NumericLimits< real64 >::min )
              {
                dJump[1] = tLocal[1] / kt;
                dJump[2] = tLocal[2] / kt;
              }

              dispJump[kfe][0] = dJump[0];
              dispJump[kfe][1] = dJump[1];
              dispJump[kfe][2] = dJump[2];

              // Compute Global Bubble Displacement: b_global = R * [u]_local
              real64 bGlobal[3];
              LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( bGlobal, faceRotationMatrix[kfe], dJump );

              // Store to totalBubbleDisplacement for both faces
              for( localIndex faceIdx = 0; faceIdx < 2; ++faceIdx )
              {
                localIndex const faceIndex = elemsToFaces( kfe, faceIdx );
                totalBubbleDisplacement[faceIndex][0] = bGlobal[0];
                totalBubbleDisplacement[faceIndex][1] = bGlobal[1];
                totalBubbleDisplacement[faceIndex][2] = bGlobal[2];
              }

              // Check Coulomb friction consistency if parameters are available
              if( hasCoulombParams )
              {
                real64 const normalTraction = tLocal[0];  // Negative for compression
                real64 const tangentialTraction = LvArray::math::sqrt( tLocal[1] * tLocal[1] + tLocal[2] * tLocal[2] );

                // Coulomb criterion: |tau| <= cohesion - mu * sigma_n (sigma_n < 0 for compression)
                // For open fracture (sigma_n > 0): no traction should be applied
                real64 const tauLimit = cohesion - frictionCoefficient * normalTraction;

                bool isInvalid = false;
                string reason;

                if( normalTraction > 0.0 )
                {
                  // Tensile normal traction - fracture should be open, no traction
                  isInvalid = true;
                  reason = "tensile normal traction (fracture should be open)";
                }
                else if( tangentialTraction > tauLimit + 1.0e-10 * LvArray::math::abs( tauLimit ) )
                {
                  // Tangential traction exceeds Coulomb limit
                  isInvalid = true;
                  reason = "tangential traction exceeds Coulomb friction limit";
                }

                if( isInvalid )
                {
                  // Get tangent vectors from rotation matrix
                  real64 const t1[3]{ faceRotationMatrix( kfe, 0, 1 ), faceRotationMatrix( kfe, 1, 1 ), faceRotationMatrix( kfe, 2, 1 ) };
                  real64 const t2[3]{ faceRotationMatrix( kfe, 0, 2 ), faceRotationMatrix( kfe, 1, 2 ), faceRotationMatrix( kfe, 2, 2 ) };

                  GEOS_LOG_RANK( GEOS_FMT( "WARNING: Stress state inconsistent with Coulomb friction law at fracture element {}.\n"
                                           "  Reason: {}\n"
                                           "  Friction parameters: cohesion = {:.6e}, friction coefficient = {:.6f}\n"
                                           "  Normal vector (global):    n  = ({:.6f}, {:.6f}, {:.6f})\n"
                                           "  Tangent vector 1 (global): t1 = ({:.6f}, {:.6f}, {:.6f})\n"
                                           "  Tangent vector 2 (global): t2 = ({:.6f}, {:.6f}, {:.6f})\n"
                                           "  Average stress from {} adjacent cell(s) (Voigt: XX, YY, ZZ, YZ, XZ, XY):\n"
                                           "    sigma = ({:.6e}, {:.6e}, {:.6e}, {:.6e}, {:.6e}, {:.6e})\n"
                                           "  Traction (local: normal, tangent1, tangent2):\n"
                                           "    t = ({:.6e}, {:.6e}, {:.6e})\n"
                                           "  Coulomb check: |tau| = {:.6e}, tau_limit = {:.6e}",
                                           kfe, reason,
                                           cohesion, frictionCoefficient,
                                           n[0], n[1], n[2],
                                           t1[0], t1[1], t1[2],
                                           t2[0], t2[1], t2[2],
                                           numValidCells,
                                           avgSigma[0], avgSigma[1], avgSigma[2], avgSigma[3], avgSigma[4], avgSigma[5],
                                           tLocal[0], tLocal[1], tLocal[2],
                                           tangentialTraction, tauLimit ) );
                }
              }
            }
          }
        } );
      }
    } );
  } );

  // Synchronize the traction, displacement jump, and bubble displacement fields
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {
    FieldIdentifiers fieldsToBeSync;
    fieldsToBeSync.addElementFields( { contact::traction::key(),
                                       contact::dispJump::key() },
                                     { getUniqueFractureRegionName() } );

    fieldsToBeSync.addFields( FieldLocation::Face,
                              { contact::totalBubbleDisplacement::key() } );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                         mesh,
                                                         domain.getNeighbors(),
                                                         true );
  } );
}

REGISTER_CATALOG_ENTRY( PhysicsSolverBase, SolidMechanicsAugmentedLagrangianContact, string const &, Group * const )
} /* namespace geos */
