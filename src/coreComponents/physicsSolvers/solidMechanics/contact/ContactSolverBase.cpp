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
 * ContactSolverBase.cpp
 */

#include "ContactSolverBase.hpp"

#include "common/TimingMacros.hpp"
#include "constitutive/contact/FrictionBase.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/SurfaceElementRegion.hpp"
#include "physicsSolvers/solidMechanics/contact/LogLevelsInfo.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "physicsSolvers/solidMechanics/contact/kernels/SolidMechanicsConformingContactKernelsBase.hpp"

namespace geos
{

using namespace dataRepository;
using namespace constitutive;
using namespace fields;

ContactSolverBase::ContactSolverBase( const string & name,
                                      Group * const parent ):
  SolidMechanicsLagrangianFEM( name, parent ),
  m_tangentRefDirection{ 1.0 / std::sqrt( 3.0 ), 1.0 / std::sqrt( 3.0 ), 1.0 / std::sqrt( 3.0 ) }
{
  this->getWrapper< string >( viewKeyStruct::contactRelationNameString() ).
    setInputFlag( dataRepository::InputFlags::FALSE );

  this->getWrapper< string >( viewKeyStruct::surfaceGeneratorNameString() ).
    setInputFlag( dataRepository::InputFlags::FALSE );

  addLogLevel< logInfo::ConfigurationStatistics >();
  addLogLevel< logInfo::ContactTolerance >();

  registerWrapper( viewKeyStruct::tangentRefDirectionString(), &m_tangentRefDirection ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_tangentRefDirection ).
    setDescription( "Reference direction used to build the canonical tangent frame on fracture elements "
                    "(t1 = normalize(refDir x normal), t2 = normal x t1). "
                    "Defaults to the isometric direction {1/sqrt(3), 1/sqrt(3), 1/sqrt(3)}. "
                    "Automatically normalised at run-time." );
}

void ContactSolverBase::postInputInitialization()
{
  SolidMechanicsLagrangianFEM::postInputInitialization();

  GEOS_THROW_IF( m_timeIntegrationOption != TimeIntegrationOption::QuasiStatic,
                 GEOS_FMT( "The attribute `{}` must be `{}`",
                           viewKeyStruct::timeIntegrationOptionString(),
                           EnumStrings< TimeIntegrationOption >::toString( TimeIntegrationOption::QuasiStatic ) ),
                 InputError, getDataContext() );

  // Normalise the user-supplied (or default) tangent reference direction.
  real64 const norm = LvArray::tensorOps::l2Norm< 3 >( m_tangentRefDirection.data );
  GEOS_THROW_IF( norm < 1.0e-10,
                 GEOS_FMT( "{}: `{}` must be a non-zero vector.",
                           getName(), viewKeyStruct::tangentRefDirectionString() ),
                 InputError, getDataContext() );
  LvArray::tensorOps::scale< 3 >( m_tangentRefDirection.data, 1.0 / norm );

  GEOS_LOG_RANK_0( GEOS_FMT( "{}: tangent reference direction = {{ {:.6f}, {:.6f}, {:.6f} }}",
                             getName(),
                             m_tangentRefDirection[0],
                             m_tangentRefDirection[1],
                             m_tangentRefDirection[2] ) );
}

void ContactSolverBase::registerDataOnMesh( dataRepository::Group & meshBodies )
{
  SolidMechanicsLagrangianFEM::registerDataOnMesh( meshBodies );

  setFractureRegions( meshBodies );

  string const labels[3] = { "normal", "tangent1", "tangent2" };
  string const labelsTangent[2] = { "tangent1", "tangent2" };

  forFractureRegionOnMeshTargets( meshBodies, [&] ( SurfaceElementRegion & fractureRegion )
  {
    fractureRegion.forElementSubRegions< SurfaceElementSubRegion >( [&]( SurfaceElementSubRegion & subRegion )
    {
      subRegion.registerField< contact::dispJump >( getName() ).
        setDimLabels( 1, labels ).
        reference().resizeDimension< 1 >( 3 );

      subRegion.registerField< contact::deltaDispJump >( getName() ).
        setDimLabels( 1, labels ).
        reference().resizeDimension< 1 >( 3 );

      subRegion.registerField< contact::oldDispJump >( getName() ).
        setDimLabels( 1, labels ).
        reference().resizeDimension< 1 >( 3 );

      subRegion.registerField< contact::dispJump_n >( getName() ).
        setDimLabels( 1, labels ).
        reference().resizeDimension< 1 >( 3 );

      subRegion.registerField< contact::traction >( getName() ).
        setDimLabels( 1, labels ).
        reference().resizeDimension< 1 >( 3 );

      subRegion.registerField< contact::fractureState >( getName() );

      subRegion.registerField< contact::oldFractureState >( getName() );

      subRegion.registerField< contact::slip >( getName() );

      subRegion.registerField< contact::tangentialTraction >( getName() );

      subRegion.registerField< contact::deltaSlip >( getName() ).
        setDimLabels( 1, labelsTangent ).reference().resizeDimension< 1 >( 2 );

      subRegion.registerField< contact::deltaSlip_n >( this->getName() ).
        setDimLabels( 1, labelsTangent ).reference().resizeDimension< 1 >( 2 );

      // Register the rotation matrix and unit frame vectors (shared by all contact solvers).
      subRegion.registerField< contact::rotationMatrix >( getName() ).
        reference().resizeDimension< 1, 2 >( 3, 3 );
    } );

  } );
}

void ContactSolverBase::setFractureRegions( dataRepository::Group const & meshBodies )
{
  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel const & mesh,
                                                    string_array const & regionNames )
  {
    mesh.getElemManager().forElementRegions< SurfaceElementRegion >( regionNames, [&] ( localIndex const, SurfaceElementRegion const & region )
    {
      m_fractureRegionNames.push_back( region.getName() );
    } );
  } );

  // TODO remove once multiple regions are fully supported
  GEOS_THROW_IF( m_fractureRegionNames.size() > 1,
                 "The number of fracture regions can not be more than one",
                 InputError, getDataContext() );
}

void ContactSolverBase::computeFractureStateStatistics( MeshLevel const & mesh,
                                                        globalIndex & numStick,
                                                        globalIndex & numNewSlip,
                                                        globalIndex & numSlip,
                                                        globalIndex & numOpen ) const
{
  using namespace contact;

  ElementRegionManager const & elemManager = mesh.getElemManager();

  array1d< globalIndex > localCounter( 4 );

  elemManager.forElementSubRegions< SurfaceElementSubRegion >( [&]( SurfaceElementSubRegion const & subRegion )
  {
    arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
    arrayView1d< integer const > const fractureState = subRegion.getField< contact::fractureState >();

    RAJA::ReduceSum< parallelHostReduce, localIndex > stickCount( 0 ), newSlipCount( 0 ), slipCount( 0 ), openCount( 0 );
    forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
    {
      if( ghostRank[kfe] < 0 )
      {
        switch( fractureState[kfe] )
        {
          case FractureState::Stick:
            {
              stickCount += 1;
              break;
            }
          case FractureState::NewSlip:
            {
              newSlipCount += 1;
              break;
            }
          case FractureState::Slip:
            {
              slipCount += 1;
              break;
            }
          case FractureState::Open:
            {
              openCount += 1;
              break;
            }
        }
      }
    } );

    localCounter[0] += stickCount.get();
    localCounter[1] += newSlipCount.get();
    localCounter[2] += slipCount.get();
    localCounter[3] += openCount.get();
  } );

  array1d< globalIndex > totalCounter( 4 );

  MpiWrapper::allReduce( localCounter,
                         totalCounter,
                         MpiWrapper::Reduction::Sum,
                         MPI_COMM_GEOS );

  numStick    = totalCounter[0];
  numNewSlip  = totalCounter[1];
  numSlip     = totalCounter[2];
  numOpen     = totalCounter[3];
}

void ContactSolverBase::outputConfigurationStatistics( DomainPartition const & domain ) const
{
  globalIndex numStick = 0;
  globalIndex numNewSlip  = 0;
  globalIndex numSlip  = 0;
  globalIndex numOpen  = 0;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel const & mesh,
                                                               string_array const & )
  {
    computeFractureStateStatistics( mesh, numStick, numNewSlip, numSlip, numOpen );

    GEOS_LOG_RANK_0( GEOS_FMT( "  Number of element for each fracture state:"
                               " stick = {:12} | new slip = {:12} | slip =  {:12} | open =  {:12}",
                               numStick, numNewSlip, numSlip, numOpen ) );
  } );
}

real64 ContactSolverBase::explicitStep( real64 const & GEOS_UNUSED_PARAM( time_n ),
                                        real64 const & dt,
                                        const int GEOS_UNUSED_PARAM( cycleNumber ),
                                        DomainPartition & GEOS_UNUSED_PARAM( domain ) )
{
  GEOS_MARK_FUNCTION;
  GEOS_ERROR( "ExplicitStep non available for contact solvers.", getDataContext() );
  return dt;
}

void ContactSolverBase::synchronizeFractureState( DomainPartition & domain ) const
{
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & )
  {
    FieldIdentifiers fieldsToBeSync;

    fieldsToBeSync.addElementFields( { contact::fractureState::key(),
                                       contact::traction::key() },
                                     { getUniqueFractureRegionName() } );

    CommunicationTools::getInstance().synchronizeFields( fieldsToBeSync,
                                                         mesh,
                                                         domain.getNeighbors(),
                                                         true );
  } );
}

void ContactSolverBase::setConstitutiveNamesCallSuper( ElementSubRegionBase & subRegion ) const
{
  if( dynamic_cast< CellElementSubRegion * >( &subRegion ) )
  {
    SolidMechanicsLagrangianFEM::setConstitutiveNamesCallSuper( subRegion );
  }
  else if( dynamic_cast< SurfaceElementSubRegion * >( &subRegion ) )
  {
    setConstitutiveName< FrictionBase >( subRegion, viewKeyStruct::frictionLawNameString(), "friction" );
  }
}

void ContactSolverBase::setupSystem( DomainPartition & domain,
                                     DofManager & dofManager,
                                     CRSMatrix< real64, globalIndex > & localMatrix,
                                     ParallelVector & rhs,
                                     ParallelVector & solution,
                                     bool const setSparsity )
{
  // Enforce: f0 < f1 (global index) → n̂₀ = -n̂₁ → nodes(kf1) ~ nodes(kf0).
  // This sequence must be preserved in this exact order and runs before
  // PhysicsSolverBase::setupSystem so the DOF structure is built on a consistent topology.
  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & )
  {
    FaceManager & faceManager = mesh.getFaceManager();
    ElementRegionManager & elemManager = mesh.getElemManager();

    NodeManager const & nodeManager = mesh.getNodeManager();

    elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
    {
      if( subRegion.size() > 0 )
      {
        // Step 1: ensure f0 is the face adjacent to the element with the smallest global index.
        subRegion.flipFaceMap( faceManager, elemManager );
        // Step 2: orient each face normal outward (may reverse node lists).
        subRegion.fixNeighboringFacesNormals( faceManager, elemManager );
        // Step 3: reorder kf1 nodes to match kf0 so conforming-contact kernels
        //         see a consistent node numbering on both sides of the fracture.
        subRegion.orderKf1NodesConsistentlyWithKf0( faceManager, nodeManager );
      }
    } );
  } );

  PhysicsSolverBase::setupSystem( domain, dofManager, localMatrix, rhs, solution, setSparsity );
}

void ContactSolverBase::computeRotationMatrices( DomainPartition & domain ) const
{
  GEOS_MARK_FUNCTION;

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                               MeshLevel & mesh,
                                                               string_array const & )
  {
    FaceManager const & faceManager = mesh.getFaceManager();
    ElementRegionManager & elemManager = mesh.getElemManager();

    arrayView2d< real64 const > const faceNormal = faceManager.faceNormal();

    // Iterate over ALL FaceElementSubRegions (not filtered by regionNames) so that
    // pre-split meshes — whose fracture region is not listed in the solver's
    // discretization target regions — are also processed.
    elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
    {
      arrayView2d< localIndex const > const elemsToFaces  = subRegion.faceList().toViewConst();
      arrayView3d< real64 >           const rotationMatrix = subRegion.getField< contact::rotationMatrix >().toView();
      arrayView2d< real64 >           const unitNormal    = subRegion.getNormalVector();
      arrayView2d< real64 >           const unitTangent1  = subRegion.getTangentVector1();
      arrayView2d< real64 >           const unitTangent2  = subRegion.getTangentVector2();

      solidMechanicsConformingContactKernels::ComputeRotationMatricesKernel::
        launch< parallelDevicePolicy<> >( subRegion.size(),
                                          faceNormal,
                                          elemsToFaces,
                                          m_tangentRefDirection,
                                          rotationMatrix,
                                          unitNormal,
                                          unitTangent1,
                                          unitTangent2 );
    } );
  } );
}

}   /* namespace geos */
