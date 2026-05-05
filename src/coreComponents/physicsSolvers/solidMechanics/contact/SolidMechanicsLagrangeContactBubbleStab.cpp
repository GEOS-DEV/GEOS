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

/**
 * @file SolidMechanicsLagrangeContactBubbleStab.cpp
 *
 */

#include "mesh/DomainPartition.hpp"
#include "SolidMechanicsLagrangeContactBubbleStab.hpp"

#include "physicsSolvers/solidMechanics/contact/kernels/SolidMechanicsConformingContactKernelsBase.hpp"
#include "physicsSolvers/solidMechanics/contact/kernels/SolidMechanicsLagrangeContactKernels.hpp"
#include "physicsSolvers/solidMechanics/contact/kernels/SolidMechanicsDisplacementJumpUpdateKernels.hpp"
#include "physicsSolvers/solidMechanics/contact/kernels/SolidMechanicsContactFaceBubbleKernels.hpp"
#include "physicsSolvers/solidMechanics/contact/LogLevelsInfo.hpp"
#include "physicsSolvers/LogLevelsInfo.hpp"
#include "finiteElement/FiniteElementDiscretization.hpp"
#include "constitutive/contact/FrictionSelector.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"


namespace geos
{

using namespace constitutive;
using namespace dataRepository;
using namespace fields;

SolidMechanicsLagrangeContactBubbleStab::SolidMechanicsLagrangeContactBubbleStab( const string & name,
                                                                                  Group * const parent ):
  ContactSolverBase( name, parent )
{
  LinearSolverParameters & linSolParams = m_linearSolverParameters.get();
  linSolParams.mgr.strategy = LinearSolverParameters::MGR::StrategyType::lagrangianContactMechanicsBubbleStab;
  linSolParams.mgr.separateComponents = true;
  linSolParams.dofsPerNode = 3;

  addLogLevel< logInfo::ResidualNorm >();
}

SolidMechanicsLagrangeContactBubbleStab::~SolidMechanicsLagrangeContactBubbleStab()
{
  // TODO Auto-generated destructor stub
}

real64 SolidMechanicsLagrangeContactBubbleStab::solverStep( real64 const & time_n,
                                                            real64 const & dt,
                                                            integer const cycleNumber,
                                                            DomainPartition & domain )
{
  if( cycleNumber == 0 )
  {
    /// Apply initial conditions to the Fault
    FieldSpecificationManager & fieldSpecificationManager = FieldSpecificationManager::getInstance();

    forDiscretizationOnMeshTargets ( domain.getMeshBodies(), [&]( string const &,
                                                                  MeshLevel & mesh,
                                                                  string_array const & )

    {
      fieldSpecificationManager.applyInitialConditions( mesh );
    } );
  }

  return ContactSolverBase::solverStep( time_n, dt, cycleNumber, domain );
}

void SolidMechanicsLagrangeContactBubbleStab::registerDataOnMesh( Group & meshBodies )
{
  ContactSolverBase::registerDataOnMesh( meshBodies );

  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & meshLevel,
                                                    string_array const & )
  {
    FaceManager & faceManager = meshLevel.getFaceManager();

    // Register the total bubble displacement
    faceManager.registerField< contact::totalBubbleDisplacement >( this->getName() ).
      reference().resizeDimension< 1 >( 3 );

    // Register the incremental bubble displacement
    faceManager.registerField< contact::incrementalBubbleDisplacement >( this->getName() ).
      reference().resizeDimension< 1 >( 3 );
  } );

  forFractureRegionOnMeshTargets( meshBodies, [&] ( SurfaceElementRegion & fractureRegion )
  {
    fractureRegion.forElementSubRegions< SurfaceElementSubRegion >( [&]( SurfaceElementSubRegion & subRegion )
    {
      // Register the rotation matrix is now handled by ContactSolverBase::registerDataOnMesh.

      subRegion.registerField< contact::deltaTraction >( getName() ).
        reference().resizeDimension< 1 >( 3 );

      subRegion.registerField< contact::targetIncrementalJump >( getName() ).
        reference().resizeDimension< 1 >( 3 );
    } );
  } );
}

void SolidMechanicsLagrangeContactBubbleStab::initializePostInitialConditionsPreSubGroups()
{
  ContactSolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  validateTetrahedralQuadrature( domain.getMeshBodies() );
}

void SolidMechanicsLagrangeContactBubbleStab::validateTetrahedralQuadrature( Group & meshBodies )
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

void SolidMechanicsLagrangeContactBubbleStab::setupDofs( DomainPartition const & domain,
                                                         DofManager & dofManager ) const
{
  GEOS_MARK_FUNCTION;

  SolidMechanicsLagrangianFEM::setupDofs( domain, dofManager );

  map< std::pair< string, string >, string_array > meshTargets;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshBodyName,
                                                                MeshLevel const & meshLevel,
                                                                string_array const & regionNames )
  {
    string_array regions;
    ElementRegionManager const & elementRegionManager = meshLevel.getElemManager();
    elementRegionManager.forElementRegions< SurfaceElementRegion >( regionNames,
                                                                    [&]( localIndex const,
                                                                         SurfaceElementRegion const & region )
    {
      regions.emplace_back( region.getName() );
    } );
    meshTargets[std::make_pair( meshBodyName, meshLevel.getName())] = std::move( regions );
  } );

  dofManager.addField( contact::totalBubbleDisplacement::key(),
                       FieldLocation::Face,
                       3,
                       meshTargets );

  dofManager.addField( contact::traction::key(),
                       FieldLocation::Elem,
                       3,
                       meshTargets );

  // Add coupling between bubble
  dofManager.addCoupling( contact::totalBubbleDisplacement::key(),
                          contact::totalBubbleDisplacement::key(),
                          DofManager::Connector::Elem );

  dofManager.addCoupling( contact::traction::key(),
                          contact::traction::key(),
                          DofManager::Connector::Elem );

  dofManager.addCoupling( solidMechanics::totalDisplacement::key(),
                          contact::traction::key(),
                          DofManager::Connector::Elem,
                          meshTargets );

  dofManager.addCoupling( contact::totalBubbleDisplacement::key(),
                          contact::traction::key(),
                          DofManager::Connector::Elem,
                          meshTargets );

  dofManager.addCoupling( solidMechanics::totalDisplacement::key(),
                          contact::totalBubbleDisplacement::key(),
                          DofManager::Connector::Elem,
                          meshTargets );
}

void SolidMechanicsLagrangeContactBubbleStab::postInputInitialization()
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

void SolidMechanicsLagrangeContactBubbleStab::setupSystem( DomainPartition & domain,
                                                           DofManager & dofManager,
                                                           CRSMatrix< real64, globalIndex > & localMatrix,
                                                           ParallelVector & rhs,
                                                           ParallelVector & solution,
                                                           bool const setSparsity )
{
  // Create the lists of interface elements that have same type.
  createFaceTypeList( domain );

  // Create the lists of interface elements that have same type and same fracture state.
  updateStickSlipList( domain );

  // Create the list of cell elements that they are enriched with bubble functions.
  createBubbleCellList( domain );

  // Enforce: f0 < f1 (global index) → n̂₀ = -n̂₁ → nodes(kf1) ~ nodes(kf0).
  // Then build the DOF structure via PhysicsSolverBase::setupSystem.
  ContactSolverBase::setupSystem( domain, dofManager, localMatrix, rhs, solution, setSparsity );

  computeRotationMatrices( domain );
}

void SolidMechanicsLagrangeContactBubbleStab::setSparsityPattern( DomainPartition & domain,
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
  this->addCouplingNumNonzeros( domain, dofManager, rowLengths.toView());

  // Create a new pattern with enough capacity for coupled matrix
  pattern.resizeFromRowCapacities< parallelHostPolicy >( patternDiag.numRows(), patternDiag.numColumns(), rowLengths.data());

  // Copy the original nonzeros
  for( localIndex localRow = 0; localRow < patternDiag.numRows(); ++localRow )
  {
    globalIndex const * cols = patternDiag.getColumns( localRow ).dataIfContiguous();
    pattern.insertNonZeros( localRow, cols, cols + patternDiag.numNonZeros( localRow ));
  }

  // Add the nonzeros from coupling
  this->addCouplingSparsityPattern( domain, dofManager, pattern.toView());
}

void SolidMechanicsLagrangeContactBubbleStab::implicitStepSetup( real64 const & time_n,
                                                                 real64 const & dt,
                                                                 DomainPartition & domain )
{
  SolidMechanicsLagrangianFEM::implicitStepSetup( time_n, dt, domain );
  computeRotationMatrices( domain );

  // Zero the incremental bubble displacement at the start of each implicit step.
  // This must be done here (not inside computeRotationMatrices) because the base-class
  // version of computeRotationMatrices knows nothing about bubble DOFs.
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & )
  {
    FaceManager & faceManager = mesh.getFaceManager();
    ElementRegionManager & elemManager = mesh.getElemManager();

    arrayView2d< real64 > const incrBubbleDisp =
      faceManager.getField< contact::incrementalBubbleDisplacement >();

    elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
    {
      arrayView2d< localIndex const > const elemsToFaces = subRegion.faceList().toViewConst();
      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
      {
        LvArray::tensorOps::fill< 3 >( incrBubbleDisp[elemsToFaces[k][0]], 0.0 );
        LvArray::tensorOps::fill< 3 >( incrBubbleDisp[elemsToFaces[k][1]], 0.0 );
      } );
    } );
  } );
}

void SolidMechanicsLagrangeContactBubbleStab::assembleSystem( real64 const time,
                                                              real64 const dt,
                                                              DomainPartition & domain,
                                                              DofManager const & dofManager,
                                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                              arrayView1d< real64 > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  SolidMechanicsLagrangianFEM::assembleSystem( time,
                                               dt,
                                               domain,
                                               dofManager,
                                               localMatrix,
                                               localRhs );

  assembleStabilization( dt, domain, dofManager, localMatrix, localRhs );

  assembleContact( dt, domain, dofManager, localMatrix, localRhs );
}

void SolidMechanicsLagrangeContactBubbleStab::assembleStabilization( real64 const dt,
                                                                     DomainPartition & domain,
                                                                     DofManager const & dofManager,
                                                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                     arrayView1d< real64 > const & localRhs )
{
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

    real64 maxTraction = finiteElement::
                           regionBasedKernelApplication
                         < parallelDevicePolicy< >,
                           ElasticIsotropic,
                           CellElementSubRegion >( mesh,
                                                   regionNames,
                                                   getDiscretizationName(),
                                                   SolidMechanicsLagrangianFEM::viewKeyStruct::solidMaterialNamesString(),
                                                   kernelFactory );

    GEOS_UNUSED_VAR( maxTraction );

  } );
}

void SolidMechanicsLagrangeContactBubbleStab::assembleContact( real64 const dt,
                                                               DomainPartition & domain,
                                                               DofManager const & dofManager,
                                                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                               arrayView1d< real64 > const & localRhs )
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const & meshName,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {
    NodeManager const & nodeManager = mesh.getNodeManager();
    FaceManager const & faceManager = mesh.getFaceManager();

    string const & dispDofKey = dofManager.getKey( solidMechanics::totalDisplacement::key() );
    string const & bubbleDofKey = dofManager.getKey( contact::totalBubbleDisplacement::key() );
    string const & tractionDofKey = dofManager.getKey( contact::traction::key() );

    arrayView1d< globalIndex const > const dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
    arrayView1d< globalIndex const > const bubbleDofNumber = faceManager.getReference< globalIndex_array >( bubbleDofKey );

    string const & fractureRegionName = this->getUniqueFractureRegionName();

    forFiniteElementOnStickFractureSubRegions( meshName, [&] ( string const &,
                                                               finiteElement::FiniteElementBase const & subRegionFE,
                                                               arrayView1d< localIndex const > const & faceElementList,
                                                               bool const )
    {
      solidMechanicsLagrangeContactKernels::LagrangeContactFactory kernelFactory( dispDofNumber,
                                                                                  bubbleDofNumber,
                                                                                  dofManager.rankOffset(),
                                                                                  localMatrix,
                                                                                  localRhs,
                                                                                  dt,
                                                                                  faceElementList,
                                                                                  tractionDofKey );

      real64 maxTraction = finiteElement::
                             interfaceBasedKernelApplication
                           < parallelDevicePolicy< >,
                             FrictionBase >( mesh,
                                             fractureRegionName,
                                             faceElementList,
                                             subRegionFE,
                                             viewKeyStruct::frictionLawNameString(),
                                             kernelFactory );

      GEOS_UNUSED_VAR( maxTraction );
    } );
  } );
}

void SolidMechanicsLagrangeContactBubbleStab::implicitStepComplete( real64 const & time_n,
                                                                    real64 const & dt,
                                                                    DomainPartition & domain )
{
  SolidMechanicsLagrangianFEM::implicitStepComplete( time_n, dt, domain );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {
    mesh.getElemManager().forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
    {
      arrayView2d< real64 > const deltaTraction  = subRegion.getField< contact::deltaTraction >();
      arrayView2d< real64 > const deltaDispJump    = subRegion.getField< contact::deltaDispJump >();
      arrayView2d< real64 const > const dispJump = subRegion.getField< contact::dispJump >();
      arrayView2d< real64 > const oldDispJump    = subRegion.getField< contact::oldDispJump >();

      forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
      {
        LvArray::tensorOps::fill< 3 >( deltaDispJump[kfe], 0.0 );
        LvArray::tensorOps::fill< 3 >( deltaTraction[kfe], 0.0 );
        LvArray::tensorOps::copy< 3 >( oldDispJump[kfe], dispJump[kfe] );
      } );
    } );
  } );
}

real64 SolidMechanicsLagrangeContactBubbleStab::calculateResidualNorm( real64 const & time,
                                                                       real64 const & dt,
                                                                       DomainPartition const & domain,
                                                                       DofManager const & dofManager,
                                                                       arrayView1d< real64 const > const & localRhs )
{
  GEOS_MARK_FUNCTION;

  real64 const solidResidual = SolidMechanicsLagrangianFEM::calculateResidualNorm( time, dt, domain, dofManager, localRhs );

  real64 const contactResidual = calculateContactResidualNorm( domain, dofManager, localRhs );

  real64 const totalResidual = sqrt( solidResidual * solidResidual + contactResidual * contactResidual );

  return totalResidual;
}

real64 SolidMechanicsLagrangeContactBubbleStab::calculateContactResidualNorm( DomainPartition const & domain,
                                                                              DofManager const & dofManager,
                                                                              arrayView1d< real64 const > const & localRhs )
{
  string const & dofKey = dofManager.getKey( contact::traction::key() );
  globalIndex const rankOffset = dofManager.rankOffset();

  real64 stickResidual = 0.0;
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel const & mesh,
                                                                string_array const & regionNames )
  {
    mesh.getElemManager().forElementSubRegions< FaceElementSubRegion >( regionNames,
                                                                        [&]( localIndex const, FaceElementSubRegion const & subRegion )
    {
      arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
      arrayView1d< real64 const > const & area = subRegion.getElementArea();

      RAJA::ReduceSum< parallelHostReduce, real64 > stickSum( 0.0 );
      forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const k )
      {
        if( ghostRank[k] < 0 )
        {
          localIndex const localRow = LvArray::integerConversion< localIndex >( dofNumber[k] - rankOffset );
          for( localIndex dim = 0; dim < 3; ++dim )
          {
            real64 const norm = localRhs[localRow + dim] / area[k];
            stickSum += norm * norm;
          }
        }
      } );

      stickResidual += stickSum.get();
    } );
  } );

  stickResidual = MpiWrapper::sum( stickResidual );
  stickResidual = sqrt( stickResidual );

  GEOS_LOG_LEVEL_RANK_0_NLR( logInfo::ResidualNorm,
                             GEOS_FMT( "        ( Rt  ) = ( {:15.6e}  )", stickResidual ));
  getConvergenceStats().setResidualValue( "Rt", stickResidual );

  return sqrt( stickResidual * stickResidual );
}


void SolidMechanicsLagrangeContactBubbleStab::applySystemSolution( DofManager const & dofManager,
                                                                   arrayView1d< real64 const > const & localSolution,
                                                                   real64 const scalingFactor,
                                                                   real64 const dt,
                                                                   DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;

  SolidMechanicsLagrangianFEM::applySystemSolution( dofManager, localSolution, scalingFactor, dt, domain );

  dofManager.addVectorToField( localSolution,
                               contact::traction::key(),
                               contact::deltaTraction::key(),
                               scalingFactor );

  dofManager.addVectorToField( localSolution,
                               contact::traction::key(),
                               contact::traction::key(),
                               scalingFactor );

  dofManager.addVectorToField( localSolution,
                               contact::totalBubbleDisplacement::key(),
                               contact::totalBubbleDisplacement::key(),
                               scalingFactor );

  dofManager.addVectorToField( localSolution,
                               contact::totalBubbleDisplacement::key(),
                               contact::incrementalBubbleDisplacement::key(),
                               scalingFactor );


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

    string const & fractureRegionName = this->getUniqueFractureRegionName();

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

      real64 maxTraction = finiteElement::
                             interfaceBasedKernelApplication
                           < parallelDevicePolicy< >,
                             NullModel >( mesh,
                                          fractureRegionName,
                                          faceElementList,
                                          subRegionFE,
                                          "",
                                          kernelFactory );

      GEOS_UNUSED_VAR( maxTraction );

    } );
  } );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {
    FieldIdentifiers fieldsToBeSync;

    fieldsToBeSync.addFields( FieldLocation::Face,
                              { contact::incrementalBubbleDisplacement::key(),
                                contact::totalBubbleDisplacement::key() } );

    fieldsToBeSync.addElementFields( { contact::traction::key(),
                                       contact::deltaTraction::key(),
                                       contact::dispJump::key() },
                                     { getUniqueFractureRegionName() } );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                         mesh,
                                                         domain.getNeighbors(),
                                                         true );
  } );
}

void SolidMechanicsLagrangeContactBubbleStab::addCouplingNumNonzeros( DomainPartition & domain,
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
      localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kfe );
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

void SolidMechanicsLagrangeContactBubbleStab::addCouplingSparsityPattern( DomainPartition const & domain,
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

      localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kfe );
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

void SolidMechanicsLagrangeContactBubbleStab::updateStickSlipList( DomainPartition const & domain )
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
    this->m_faceTypesToFaceElementsStick.get_inserted( meshName );
    this->m_faceTypesToFaceElementsSlip.get_inserted( meshName );

    forFiniteElementOnFractureSubRegions( meshName, [&] ( string const & finiteElementName,
                                                          finiteElement::FiniteElementBase const &,
                                                          arrayView1d< localIndex const > const & faceElementList )
    {

      array1d< localIndex > keys( subRegion.size());
      array1d< localIndex > vals( subRegion.size());
      array1d< localIndex > stickList;
      array1d< localIndex > slipList;
      RAJA::ReduceSum< ReducePolicy< parallelDevicePolicy<> >, localIndex > nStick_r( 0 );
      RAJA::ReduceSum< ReducePolicy< parallelDevicePolicy<> >, localIndex > nSlip_r( 0 );

      arrayView1d< localIndex > const keys_v = keys.toView();
      arrayView1d< localIndex > const vals_v = vals.toView();
      forAll< parallelDevicePolicy<> >( faceElementList.size(),
                                        [ = ]
                                        GEOS_HOST_DEVICE ( localIndex const kfe )
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
      RAJA::sort_pairs< parallelDevicePolicy<> >( keys_v, vals_v );

      stickList.resize( nStick );
      slipList.resize( nSlip );
      arrayView1d< localIndex > const stickList_v = stickList.toView();
      arrayView1d< localIndex > const slipList_v = slipList.toView();

      forAll< parallelDevicePolicy<> >( nStick, [ = ]
                                        GEOS_HOST_DEVICE ( localIndex const kfe )
      {
        stickList_v[kfe] = vals_v[kfe];
      } );

      forAll< parallelDevicePolicy<> >( nSlip, [ = ]
                                        GEOS_HOST_DEVICE ( localIndex const kfe )
      {
        slipList_v[kfe] = vals_v[nStick+kfe];
      } );

      this->m_faceTypesToFaceElementsStick.get_inserted( meshName ).get_inserted( finiteElementName ) = stickList;
      this->m_faceTypesToFaceElementsSlip.get_inserted( meshName ).get_inserted( finiteElementName ) = slipList;

      GEOS_LOG_LEVEL_RANK_0( logInfo::ConfigurationStatistics, GEOS_FMT( "# stick elements: {}, # slip elements: {}", nStick, nSlip ))
    } );
  } );

}

void SolidMechanicsLagrangeContactBubbleStab::createFaceTypeList( DomainPartition const & domain )
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
      this->m_faceTypesToFaceElements.get_inserted( meshName );
      return;
    }

    arrayView2d< localIndex const > const elemsToFaces = subRegion.faceList().toViewConst();

    array1d< localIndex > keys( subRegion.size());
    array1d< localIndex > vals( subRegion.size());
    array1d< localIndex > quadList;
    array1d< localIndex > triList;
    RAJA::ReduceSum< ReducePolicy< parallelDevicePolicy<> >, localIndex > nTri_r( 0 );
    RAJA::ReduceSum< ReducePolicy< parallelDevicePolicy<> >, localIndex > nQuad_r( 0 );

    arrayView1d< localIndex > const keys_v = keys.toView();
    arrayView1d< localIndex > const vals_v = vals.toView();
    // Determine the size of the lists and generate the vector keys and vals for parallel indexing into lists.
    // (With RAJA, parallelizing this operation seems the most viable approach.)
    forAll< parallelDevicePolicy<> >( subRegion.size(),
                                      [ = ] GEOS_HOST_DEVICE ( localIndex const kfe )
    {

      localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elemsToFaces[kfe][0] );
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
        GEOS_ERROR( "SolidMechanicsLagrangeContactBubbleStab:: invalid face type" );
      }
    } );

    localIndex nQuad = static_cast< localIndex >(nQuad_r.get());
    localIndex nTri = static_cast< localIndex >(nTri_r.get());

    // Sort vals according to keys to ensure that
    // elements of the same type are adjacent in the vals list.
    // This arrangement allows for efficient copying into the container
    // by leveraging parallelism.
    RAJA::sort_pairs< parallelDevicePolicy<> >( keys_v, vals_v );

    quadList.resize( nQuad );
    triList.resize( nTri );
    arrayView1d< localIndex > const quadList_v = quadList.toView();
    arrayView1d< localIndex > const triList_v = triList.toView();

    forAll< parallelDevicePolicy<> >( nTri, [ = ] GEOS_HOST_DEVICE ( localIndex const kfe )
    {
      triList_v[kfe] = vals_v[kfe];
    } );

    forAll< parallelDevicePolicy<> >( nQuad, [ = ] GEOS_HOST_DEVICE ( localIndex const kfe )
    {
      quadList_v[kfe] = vals_v[nTri+kfe];
    } );

    this->m_faceTypesToFaceElements.get_inserted( meshName ).get_inserted( "Quadrilateral" ) = quadList;
    this->m_faceTypesToFaceElements.get_inserted( meshName ).get_inserted( "Triangle" ) = triList;
  } );

}

void SolidMechanicsLagrangeContactBubbleStab::createBubbleCellList( DomainPartition & domain ) const
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

      forAll< parallelDevicePolicy<> >( subRegion.size(), [ = ] GEOS_HOST_DEVICE ( localIndex const kfe )
      {

        localIndex const kf0 = elemsToFaces[kfe][0], kf1 = elemsToFaces[kfe][1];
        tmpSpace_v[2*kfe] = kf0, tmpSpace_v[2*kfe+1] = kf1;

      } );
    }

    // Sort indexes to enable efficient searching using binary search.
    RAJA::stable_sort< parallelDevicePolicy<> >( tmpSpace_v );
    // Move data back to host: after the device sort the buffer lives on the GPU,
    // but SortedArray::insert is a host-only operation that dereferences the iterators.
    tmpSpace_v.move( LvArray::MemorySpace::host );
    faceIdList.insert( tmpSpace_v.begin(), tmpSpace_v.end());

    // Search for bubble element on each CellElementSubRegion and
    // store element indexes, global and local face indexes.
    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const, CellElementSubRegion & cellElementSubRegion )
    {

      arrayView2d< localIndex const > const elemsToFaces = cellElementSubRegion.faceList().toViewConst();

      RAJA::ReduceSum< ReducePolicy< parallelDevicePolicy<> >, localIndex > nBubElems_r( 0 );

      localIndex const n_max = cellElementSubRegion.size() * elemsToFaces.size( 1 );
      array1d< localIndex > keys( n_max );
      array1d< localIndex > perms( n_max );
      array1d< localIndex > vals( n_max );
      array1d< localIndex > localFaceIds( n_max );

      arrayView1d< localIndex > const keys_v = keys.toView();
      arrayView1d< localIndex > const perms_v = perms.toView();
      arrayView1d< localIndex > const vals_v = vals.toView();
      arrayView1d< localIndex > const localFaceIds_v = localFaceIds.toView();
      // Move faceIdList to device so that the SortedArrayView can be safely
      // accessed inside the GPU kernel (faceIdList was populated on the host).
      faceIdList.move( parallelDeviceMemorySpace );
      SortedArrayView< localIndex const > const faceIdList_v = faceIdList.toViewConst();

      forAll< parallelDevicePolicy<> >( cellElementSubRegion.size(),
                                        [ = ]
                                        GEOS_HOST_DEVICE ( localIndex const kfe )
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
      RAJA::sort_pairs< parallelDevicePolicy<> >( keys_v, perms_v );

      array1d< localIndex > bubbleElemsList;
      bubbleElemsList.resize( nBubElems );

      arrayView1d< localIndex > const bubbleElemsList_v = bubbleElemsList.toView();

      forAll< parallelDevicePolicy<> >( n_max, [ = ] GEOS_HOST_DEVICE ( localIndex const k )
      {
        keys_v[k] = vals_v[perms_v[k]];
      } );

      forAll< parallelDevicePolicy<> >( nBubElems, [ = ] GEOS_HOST_DEVICE ( localIndex const k )
      {
        bubbleElemsList_v[k] = keys_v[k];
      } );
      cellElementSubRegion.setBubbleElementsList( bubbleElemsList.toViewConst());

      forAll< parallelDevicePolicy<> >( n_max, [ = ] GEOS_HOST_DEVICE ( localIndex const k )
      {
        keys_v[k] = localFaceIds_v[perms_v[k]];
      } );

      array2d< localIndex > faceElemsList;
      faceElemsList.resize( nBubElems, 2 );

      arrayView2d< localIndex > const faceElemsList_v = faceElemsList.toView();

      forAll< parallelDevicePolicy<> >( nBubElems,
                                        [ = ]
                                        GEOS_HOST_DEVICE ( localIndex const k )
      {
        localIndex const kfe =  bubbleElemsList_v[k];
        faceElemsList_v[k][0] = elemsToFaces[kfe][keys_v[k]];
        faceElemsList_v[k][1] = keys_v[k];
      } );
      cellElementSubRegion.setFaceElementsList( faceElemsList.toViewConst());

    } );

  } );

}

REGISTER_CATALOG_ENTRY( PhysicsSolverBase, SolidMechanicsLagrangeContactBubbleStab, string const &, Group * const )

} /* namespace geos */
