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
#include "finiteElement/elementFormulations/H1_TriangleFace_Lagrange1_Gauss.hpp"
#include "finiteElement/elementFormulations/H1_QuadrilateralFace_Lagrange1_GaussLegendre2.hpp"

namespace geos
{

using namespace dataRepository;
using namespace constitutive;
using namespace fields;
using namespace finiteElement;


ContactSolverBase::ContactSolverBase( const string & name,
                                      Group * const parent ):
  SolidMechanicsLagrangianFEM( name, parent )
{
  this->getWrapper< string >( viewKeyStruct::contactRelationNameString() ).
    setInputFlag( dataRepository::InputFlags::FALSE );

  this->getWrapper< string >( viewKeyStruct::surfaceGeneratorNameString() ).
    setInputFlag( dataRepository::InputFlags::FALSE );

  addLogLevel< logInfo::ConfigurationStatistics >();
  addLogLevel< logInfo::ContactTolerance >();
}

void ContactSolverBase::postInputInitialization()
{
  SolidMechanicsLagrangianFEM::postInputInitialization();

  GEOS_THROW_IF( m_timeIntegrationOption != TimeIntegrationOption::QuasiStatic,
                 GEOS_FMT( "The attribute `{}` must be `{}`",
                           viewKeyStruct::timeIntegrationOptionString(),
                           EnumStrings< TimeIntegrationOption >::toString( TimeIntegrationOption::QuasiStatic ) ),
                 InputError, getDataContext() );
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

void ContactSolverBase::computeFaceNodalArea( localIndex const kf0,
                                                          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
                                                          ArrayOfArraysView< localIndex const > const & faceToNodeMap,
                                                          ArrayOfArraysView< localIndex const > const & faceToEdgeMap,
                                                          arrayView2d< localIndex const > const & edgeToNodeMap,
                                                          arrayView2d< real64 const > const faceCenters,
                                                          arrayView2d< real64 const > const faceNormals,
                                                          arrayView1d< real64 const > const faceAreas,
                                                          stackArray1d< real64, FaceManager::maxFaceNodes() > & basisIntegrals ) const
{
  GEOS_MARK_FUNCTION;
  localIndex const TriangularPermutation[3] = { 0, 1, 2 };
  localIndex const QuadrilateralPermutation[4] = { 0, 1, 3, 2 };
  localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kf0 );

  basisIntegrals.resize( numNodesPerFace );
  for( localIndex a = 0; a < numNodesPerFace; ++a )
  {
    basisIntegrals[a] = 0.0;
  }
  localIndex const * const permutation = ( numNodesPerFace == 3 ) ? TriangularPermutation : QuadrilateralPermutation;
  if( numNodesPerFace == 3 )
  {
    real64 xLocal[3][3];
    for( localIndex a = 0; a < numNodesPerFace; ++a )
    {
      for( localIndex j = 0; j < 3; ++j )
      {
        xLocal[a][j] = nodePosition[faceToNodeMap( kf0, permutation[a] )][j];
      }
    }
    real64 N[3];
    for( localIndex q=0; q<H1_TriangleFace_Lagrange1_Gauss1::numQuadraturePoints; ++q )
    {
      real64 const detJ = H1_TriangleFace_Lagrange1_Gauss1::transformedQuadratureWeight( q, xLocal );
      H1_TriangleFace_Lagrange1_Gauss1::calcN( q, N );
      for( localIndex a = 0; a < numNodesPerFace; ++a )
      {
        basisIntegrals[a] += detJ * N[permutation[a]];
      }
    }
  }
  else if( numNodesPerFace == 4 )
  {
    real64 xLocal[4][3];
    for( localIndex a = 0; a < numNodesPerFace; ++a )
    {
      for( localIndex j = 0; j < 3; ++j )
      {
        xLocal[a][j] = nodePosition[faceToNodeMap( kf0, permutation[a] )][j];
      }
    }
    real64 N[4];
    for( localIndex q=0; q<H1_QuadrilateralFace_Lagrange1_GaussLegendre2::numQuadraturePoints; ++q )
    {
      real64 const detJ = H1_QuadrilateralFace_Lagrange1_GaussLegendre2::transformedQuadratureWeight( q, xLocal );
      H1_QuadrilateralFace_Lagrange1_GaussLegendre2::calcN( q, N );
      for( localIndex a = 0; a < numNodesPerFace; ++a )
      {
        basisIntegrals[a] += detJ * N[permutation[a]];
      }
    }
  }
  else if( numNodesPerFace > 4 && numNodesPerFace <= m_maxFaceNodes )
  {
    // we need to L2 projector based on VEM to approximate the quadrature weights
    // we need to use extra geometry information to computing L2 projector

    localIndex const MFN = m_maxFaceNodes; // Max number of face vertices.
    localIndex const faceIndex = kf0;
    localIndex const numFaceNodes = faceToNodeMap[ faceIndex ].size();

    // get the face center and normal.
    real64 const faceArea = faceAreas[ faceIndex ];
    localIndex faceToNodes[ MFN ];
    localIndex faceToEdges[ MFN ];
    for( localIndex i = 0; i < numFaceNodes; ++i )
    {
      faceToNodes[i] = faceToNodeMap[ faceIndex ][ i ];
      faceToEdges[i] = faceToEdgeMap[ faceIndex ][ i ];
    }
    // - get outward face normal and center
    real64 faceNormal[3] = { faceNormals[faceIndex][0],
                             faceNormals[faceIndex][1],
                             faceNormals[faceIndex][2] };
    real64 const faceCenter[3] { faceCenters[faceIndex][0],
                                 faceCenters[faceIndex][1],
                                 faceCenters[faceIndex][2] };
    // - compute integrals calling auxiliary method
    real64 threeDMonomialIntegrals[3] = { 0.0 };
    real64 const invCellDiameter = 0.0;
    real64 const cellCenter[3] { 0.0, 0.0, 0.0 };
    computeFaceIntegrals( nodePosition,
                          faceToNodes,
                          faceToEdges,
                          numFaceNodes,
                          faceArea,
                          faceCenter,
                          faceNormal,
                          edgeToNodeMap,
                          invCellDiameter,
                          cellCenter,
                          basisIntegrals,
                          threeDMonomialIntegrals );
  }
  else
  {
    GEOS_ERROR( GEOS_FMT( "Face with {} nodes. Only triangles and quadrilaterals and PEBI prisms up to 11 sides are supported.",
                          numNodesPerFace ),
                getDataContext() );
  }
}

void ContactSolverBase::computeFaceIntegrals( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodesCoords,
                                                          localIndex const (&faceToNodes)[11],
                                                          localIndex const (&faceToEdges)[11],
                                                          localIndex const & numFaceVertices,
                                                          real64 const & faceArea,
                                                          real64 const (&faceCenter)[3],
                                                          real64 const (&faceNormal)[3],
                                                          arrayView2d< localIndex const > const & edgeToNodes,
                                                          real64 const & invCellDiameter,
                                                          real64 const (&cellCenter)[3],
                                                          stackArray1d< real64, FaceManager::maxFaceNodes() > & basisIntegrals,
                                                          real64 (& threeDMonomialIntegrals)[3] ) const
{
  GEOS_MARK_FUNCTION;
  localIndex const MFN = m_maxFaceNodes; // Max number of face vertices.
  basisIntegrals.resize( numFaceVertices );
  // Rotate the face.
  //  - compute rotation matrix.
  real64 faceRotationMatrix[ 3 ][ 3 ];
  computationalGeometry::RotationMatrix_3D( faceNormal, faceRotationMatrix );
  //  - below we compute the diameter, the rotated vertices and the rotated center.
  real64 faceRotatedVertices[ MFN ][ 2 ];
  real64 faceDiameter = 0;

  for( localIndex numVertex = 0; numVertex < numFaceVertices; ++numVertex )
  {
    // apply the transpose (that is the inverse) of the rotation matrix to face vertices.
    // NOTE:
    // the second and third rows of the transpose of the rotation matrix rotate on the 2D face.
    faceRotatedVertices[numVertex][0] =
      faceRotationMatrix[ 0 ][ 1 ]*nodesCoords( faceToNodes[ numVertex ], 0 ) +
      faceRotationMatrix[ 1 ][ 1 ]*nodesCoords( faceToNodes[ numVertex ], 1 ) +
      faceRotationMatrix[ 2 ][ 1 ]*nodesCoords( faceToNodes[ numVertex ], 2 );
    faceRotatedVertices[numVertex][1] =
      faceRotationMatrix[ 0 ][ 2 ]*nodesCoords( faceToNodes[ numVertex ], 0 ) +
      faceRotationMatrix[ 1 ][ 2 ]*nodesCoords( faceToNodes[ numVertex ], 1 ) +
      faceRotationMatrix[ 2 ][ 2 ]*nodesCoords( faceToNodes[ numVertex ], 2 );
  }

  faceDiameter = computationalGeometry::computeDiameter< 2 >( faceRotatedVertices,
                                                              numFaceVertices );
  real64 const invFaceDiameter = 1.0/faceDiameter;
  // - rotate the face centroid as done for the vertices.
  real64 faceRotatedCentroid[2];
  faceRotatedCentroid[0] =
    faceRotationMatrix[ 0 ][ 1 ]*faceCenter[0] +
    faceRotationMatrix[ 1 ][ 1 ]*faceCenter[1] +
    faceRotationMatrix[ 2 ][ 1 ]*faceCenter[2];
  faceRotatedCentroid[1] =
    faceRotationMatrix[ 0 ][ 2 ]*faceCenter[0] +
    faceRotationMatrix[ 1 ][ 2 ]*faceCenter[1] +
    faceRotationMatrix[ 2 ][ 2 ]*faceCenter[2];
  // - compute edges' lengths, outward pointing normals and local edge-to-nodes map.
  real64 edgeOutwardNormals[ MFN ][ 2 ];
  real64 edgeLengths[ MFN ];
  localIndex localEdgeToNodes[ MFN ][ 2 ];

  for( localIndex numEdge = 0; numEdge < numFaceVertices; ++numEdge )
  {
    if( edgeToNodes( faceToEdges[numEdge], 0 ) == faceToNodes[ numEdge ] )
    {
      localEdgeToNodes[ numEdge ][ 0 ] = numEdge;
      localEdgeToNodes[ numEdge ][ 1 ] = (numEdge+1)%numFaceVertices;
    }
    else
    {
      localEdgeToNodes[ numEdge ][ 0 ] = (numEdge+1)%numFaceVertices;
      localEdgeToNodes[ numEdge ][ 1 ] = numEdge;
    }
    real64 edgeTangent[2];
    edgeTangent[0] = faceRotatedVertices[(numEdge+1)%numFaceVertices][0] -
                     faceRotatedVertices[numEdge][0];
    edgeTangent[1] = faceRotatedVertices[(numEdge+1)%numFaceVertices][1] -
                     faceRotatedVertices[numEdge][1];
    edgeOutwardNormals[numEdge][0] = edgeTangent[1];
    edgeOutwardNormals[numEdge][1] = -edgeTangent[0];
    real64 signTestVector[2];
    signTestVector[0] = faceRotatedVertices[numEdge][0] - faceRotatedCentroid[0];
    signTestVector[1] = faceRotatedVertices[numEdge][1] - faceRotatedCentroid[1];
    if( signTestVector[0]*edgeOutwardNormals[numEdge][0] +
        signTestVector[1]*edgeOutwardNormals[numEdge][1] < 0 )
    {
      edgeOutwardNormals[numEdge][0] = -edgeOutwardNormals[numEdge][0];
      edgeOutwardNormals[numEdge][1] = -edgeOutwardNormals[numEdge][1];
    }
    edgeLengths[numEdge] = LvArray::math::sqrt< real64 >( edgeTangent[0]*edgeTangent[0] +
                                                          edgeTangent[1]*edgeTangent[1] );
    edgeOutwardNormals[numEdge][0] /= edgeLengths[numEdge];
    edgeOutwardNormals[numEdge][1] /= edgeLengths[numEdge];
  }

  // Compute boundary quadrature weights (also equal to the integrals of basis functions on the
  // boundary).
  real64 boundaryQuadratureWeights[ MFN ];
  for( localIndex numWeight = 0; numWeight < numFaceVertices; ++numWeight )
    boundaryQuadratureWeights[numWeight] = 0.0;
  for( localIndex numEdge = 0; numEdge < numFaceVertices; ++numEdge )
  {
    boundaryQuadratureWeights[ localEdgeToNodes[ numEdge ][ 0 ] ] += 0.5*edgeLengths[numEdge];
    boundaryQuadratureWeights[ localEdgeToNodes[ numEdge ][ 1 ] ] += 0.5*edgeLengths[numEdge];
  }

  // Compute scaled monomials' integrals on edges.
  real64 monomBoundaryIntegrals[3] = { 0.0 };
  for( localIndex numVertex = 0; numVertex < numFaceVertices; ++numVertex )
  {
    monomBoundaryIntegrals[0] += boundaryQuadratureWeights[ numVertex ];
    monomBoundaryIntegrals[1] += (faceRotatedVertices[ numVertex ][ 0 ] - faceRotatedCentroid[0]) *
                                 invFaceDiameter*boundaryQuadratureWeights[ numVertex ];
    monomBoundaryIntegrals[2] += (faceRotatedVertices[ numVertex ][ 1 ] - faceRotatedCentroid[1]) *
                                 invFaceDiameter*boundaryQuadratureWeights[ numVertex ];
  }

  // Compute non constant 2D and 3D scaled monomials' integrals on the face.
  real64 monomInternalIntegrals[2] = { 0.0 };
  for( localIndex numSubTriangle = 0; numSubTriangle < numFaceVertices; ++numSubTriangle )
  {
    localIndex const nextVertex = (numSubTriangle+1)%numFaceVertices;
    // - compute value of 2D monomials at the quadrature point on the sub-triangle (the
    //   barycenter).
    //   The result is ((v(0)+v(1)+faceCenter)/3 - faceCenter) / faceDiameter =
    //   = (v(0) + v(1) - 2*faceCenter)/(3*faceDiameter).
    real64 monomialValues[2];
    for( localIndex i = 0; i < 2; ++i )
    {
      monomialValues[i] = (faceRotatedVertices[numSubTriangle][i] +
                           faceRotatedVertices[nextVertex][i] -
                           2.0*faceRotatedCentroid[i]) / (3.0*faceDiameter);
    }
    // compute value of 3D monomials at the quadrature point on the sub-triangle (the
    // barycenter).  The result is
    // ((v(0) + v(1) + faceCenter)/3 - cellCenter)/cellDiameter.
    real64 threeDMonomialValues[3];
    for( localIndex i = 0; i < 3; ++i )
    {
      threeDMonomialValues[i] = ( (faceCenter[i] +
                                   nodesCoords[faceToNodes[ numSubTriangle ]][i] +
                                   nodesCoords[faceToNodes[ nextVertex ]][i]) / 3.0 -
                                  cellCenter[i] ) * invCellDiameter;
    }
    // compute quadrature weight associated to the quadrature point (the area of the
    // sub-triangle).
    real64 edgesTangents[2][2];               // used to compute the area of the sub-triangle
    for( localIndex i = 0; i < 2; ++i )
    {
      edgesTangents[0][i] = faceRotatedVertices[numSubTriangle][i] - faceRotatedCentroid[i];
    }
    for( localIndex i = 0; i < 2; ++i )
    {
      edgesTangents[1][i] = faceRotatedVertices[nextVertex][i] - faceRotatedCentroid[i];
    }
    real64 subTriangleArea = 0.5*LvArray::math::abs
                               ( edgesTangents[0][0]*edgesTangents[1][1] -
                               edgesTangents[0][1]*edgesTangents[1][0] );
    // compute the integrals on the sub-triangle and add it to the global integrals
    for( localIndex i = 0; i < 2; ++i )
    {
      monomInternalIntegrals[ i ] += monomialValues[ i ]*subTriangleArea;
    }
    for( localIndex i = 0; i < 3; ++i )
    {
      // threeDMonomialIntegrals is assumed to be initialized to 0 by the caller
      threeDMonomialIntegrals[ i ] += threeDMonomialValues[ i ]*subTriangleArea;
    }
  }

  // Compute integral of basis functions times normal derivative of monomials on the boundary.
  real64 basisTimesMonomNormalDerBoundaryInt[ MFN ][ 2 ];
  for( localIndex numVertex = 0; numVertex < numFaceVertices; ++numVertex )
  {
    for( localIndex i = 0; i < 2; ++i )
    {
      basisTimesMonomNormalDerBoundaryInt[ numVertex ][ i ] = 0.0;
    }
  }
  for( localIndex numVertex = 0; numVertex < numFaceVertices; ++numVertex )
  {
    for( localIndex i = 0; i < 2; ++i )
    {
      real64 thisEdgeIntTimesNormal_i = edgeOutwardNormals[numVertex][i]*edgeLengths[numVertex];
      basisTimesMonomNormalDerBoundaryInt[ localEdgeToNodes[ numVertex ][ 0 ] ][i] += thisEdgeIntTimesNormal_i;
      basisTimesMonomNormalDerBoundaryInt[ localEdgeToNodes[ numVertex ][ 1 ] ][i] += thisEdgeIntTimesNormal_i;
    }
  }
  for( localIndex numVertex = 0; numVertex < numFaceVertices; ++numVertex )
  {
    for( localIndex i = 0; i < 2; ++i )
    {
      basisTimesMonomNormalDerBoundaryInt[ numVertex ][ i ] *= 0.5*invFaceDiameter;
    }
  }

  // Compute integral mean of basis functions on this face.
  real64 const invFaceArea = 1.0/faceArea;
  real64 const monomialDerivativeInverse = (faceDiameter*faceDiameter)*invFaceArea;
  for( localIndex numVertex = 0; numVertex < numFaceVertices; ++numVertex )
  {
    real64 piNablaDofs[ 3 ];
    piNablaDofs[ 1 ] = monomialDerivativeInverse *
                       basisTimesMonomNormalDerBoundaryInt[ numVertex ][ 0 ];
    piNablaDofs[ 2 ] = monomialDerivativeInverse *
                       basisTimesMonomNormalDerBoundaryInt[ numVertex ][ 1 ];
    piNablaDofs[ 0 ] = (boundaryQuadratureWeights[ numVertex ] -
                        piNablaDofs[ 1 ]*monomBoundaryIntegrals[ 1 ] -
                        piNablaDofs[ 2 ]*monomBoundaryIntegrals[ 2 ])/monomBoundaryIntegrals[ 0 ];
    basisIntegrals[ numVertex ] = piNablaDofs[ 0 ]*faceArea +
                                  (piNablaDofs[ 1 ]*monomInternalIntegrals[ 0 ] +
                                   piNablaDofs[ 2 ]*monomInternalIntegrals[ 1 ]);
  }
}

}   /* namespace geos */
