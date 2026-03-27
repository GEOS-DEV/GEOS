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
 * @file TractionBoundaryCondition.cpp
 */

#include "TractionBoundaryCondition.hpp"

#include "finiteElement/elementFormulations/H1_TriangleFace_Lagrange1_Gauss.hpp"
#include "finiteElement/elementFormulations/H1_QuadrilateralFace_Lagrange1_GaussLegendre2.hpp"
#include "functions/TableFunction.hpp"

namespace geos
{

using namespace dataRepository;

TractionBoundaryCondition::TractionBoundaryCondition( string const & name, Group * parent ):
  FieldSpecificationBase( name, parent ),
  m_tractionType( TractionType::vector ),
  m_inputStress{},
  m_scaleSet(),
  m_nodalScaleFlag( 0 )//,
//  m_stressFunctionNames(),
//  m_useStressFunctions(false),
//  m_stressFunctions{nullptr}
{
  registerWrapper( viewKeyStruct::tractionTypeString(), &m_tractionType ).
    setApplyDefaultValue( TractionType::vector ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 ).
    setDescription( "Type of traction boundary condition. Options are:\n" +
                    toString( TractionType::vector ) + " - traction is applied to the faces as specified from the scale and direction,\n" +
                    toString( TractionType::normal ) + " - traction is applied to the faces as a pressure specified from the product of scale and the outward face normal,\n" +
                    toString( TractionType::stress ) + " - traction is applied to the faces as specified by the inner product of input stress and face normal." );

  registerWrapper( viewKeyStruct::inputStressString(), &m_inputStress ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( string( "Input stress for " ) + viewKeyStruct::tractionTypeString() + " = " + toString( TractionType::stress ) );

//  registerWrapper( viewKeyStruct::stressFunctionString(), &m_stressFunctionNames ).
//    setInputFlag( InputFlags::OPTIONAL ).
//    setDescription( string("Function names for description of stress for ") + viewKeyStruct::tractionTypeString() +
//                    " = " + toString( TractionType::stress ) + ". Overrides " + viewKeyStruct::inputStressString() + "." );

  registerWrapper( viewKeyStruct::scaleSetString(), &m_scaleSet ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::nodalScaleFlagString(), &m_nodalScaleFlag ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "The flag to indicate whether to apply the nodal scale on the traction magnitude" );

  getWrapper< string >( FieldSpecificationBase::viewKeyStruct::fieldNameString() ).
    setInputFlag( InputFlags::FALSE );
  setFieldName( catalogName() );

  getWrapper< int >( FieldSpecificationBase::viewKeyStruct::componentString() ).
    setInputFlag( InputFlags::FALSE );
}


void TractionBoundaryCondition::postInputInitialization()
{
  if( m_tractionType == TractionType::vector )
  {
    GEOS_ERROR_IF( LvArray::tensorOps::l2Norm< 3 >( getDirection() ) < 1e-20,
                   GEOS_FMT( "{} is required for {} = {}, but appears to be unspecified",
                             viewKeyStruct::directionString(),
                             viewKeyStruct::tractionTypeString(),
                             EnumStrings< TractionType >::toString( TractionType::vector ) ),
                   getDataContext() );
  }
  else
  {
    GEOS_LOG_RANK_0_IF( LvArray::tensorOps::l2Norm< 3 >( getDirection() ) > 1e-20,
                        viewKeyStruct::directionString() << " is not required unless " <<
                        viewKeyStruct::tractionTypeString() << " = " <<
                        EnumStrings< TractionType >::toString( TractionType::vector ) <<
                        ", but appears to be specified" );
  }

  bool const inputStressRead = getWrapper< R2SymTensor >( viewKeyStruct::inputStressString() ).getSuccessfulReadFromInput();

  GEOS_LOG_RANK_0_IF( inputStressRead && m_tractionType != TractionType::stress,
                      viewKeyStruct::inputStressString() << " is specified, but " <<
                      viewKeyStruct::tractionTypeString() << " != " <<
                      EnumStrings< TractionType >::toString( TractionType::stress ) <<
                      ", so value of " << viewKeyStruct::inputStressString() << " is unused." );

  GEOS_ERROR_IF( !inputStressRead && m_tractionType == TractionType::stress,
                 GEOS_FMT( "{} = {}, but {} is not specified.",
                           viewKeyStruct::tractionTypeString(),
                           EnumStrings< TractionType >::toString( TractionType::stress ),
                           viewKeyStruct::inputStressString() ),
                 getDataContext() );


//  localIndex const numStressFunctionsNames = m_stressFunctionNames.size();
//  GEOS_ERROR_IF( numStressFunctionsNames > 0 && numStressFunctionsNames<6,
//                  "Either 0 or 6 stress functions must be specified using stressFunctions" );
//
//  if( numStressFunctionsNames==6 )
//  {
//    m_useStressFunctions = true;
//  }
}

void TractionBoundaryCondition::initializePreSubGroups()
{
//  FunctionManager const & functionManager = getGlobalState().getFunctionManager();
//
//  if( m_useStressFunctions )
//  {
//    for( localIndex i=0; i<6; ++i )
//    {
//      m_stressFunctions[i] = &( functionManager.getGroup< TableFunction >( m_stressFunctionNames[i] ) );
//    }
//  }
}


/**
 * @brief Integrate a uniform traction over a face using a standard FE face formulation.
 *
 * @tparam FE_TYPE The finite element face type (e.g. H1_TriangleFace_Lagrange1_Gauss1_impl,
 *                 H1_QuadrilateralFace_Lagrange1_GaussLegendre2_impl).
 *                 Must provide: numNodes, numQuadraturePoints, calcN(), transformedQuadratureWeight().
 * @param traction The traction vector [3] to integrate.
 * @param xFace    The face node coordinates in FE ordering [numNodes][3].
 * @param permutation Mapping from FE node index to mesh face-local node index.
 *                     Accounts for ordering differences (e.g. tensor-product vs CCW).
 * @param faceToNodeMap The face-to-node connectivity map.
 * @param kf       The face index.
 * @param blockLocalDofNumber Block-local DOF numbering.
 * @param dofRankOffset The rank offset for DOFs.
 * @param localRhs The local RHS vector to assemble into.
 */
template< typename FE_TYPE >
GEOS_HOST_DEVICE
inline
void integrateFaceTraction( real64 const ( &traction )[3],
                            real64 const ( &xFace )[FE_TYPE::numNodes][3],
                            int const ( &permutation )[FE_TYPE::numNodes],
                            ArrayOfArraysView< localIndex const > const & faceToNodeMap,
                            localIndex const kf,
                            arrayView1d< globalIndex const > const & blockLocalDofNumber,
                            globalIndex const dofRankOffset,
                            arrayView1d< real64 > const & localRhs )
{
  constexpr localIndex numNodes = FE_TYPE::numNodes;
  constexpr localIndex numQuadraturePoints = FE_TYPE::numQuadraturePoints;

  // Accumulate per-node integration weight: w_a = sum_q N_a(q) * |J(q)| * wq
  real64 nodalWeight[numNodes] {};

  for( localIndex q = 0; q < numQuadraturePoints; ++q )
  {
    real64 N[numNodes];
    FE_TYPE::calcN( q, N );

    real64 const detJxW = FE_TYPE::transformedQuadratureWeight( q, xFace );

    for( localIndex a = 0; a < numNodes; ++a )
    {
      nodalWeight[a] += N[a] * detJxW;
    }
  }

  // Assemble traction * nodalWeight into the RHS.
  // Use the permutation to map FE node 'a' back to the mesh face-local node index.
  for( localIndex a = 0; a < numNodes; ++a )
  {
    localIndex const meshNodeIdx = permutation[a];
    localIndex const dof = blockLocalDofNumber[ faceToNodeMap( kf, meshNodeIdx ) ] - dofRankOffset;
    if( dof < 0 || dof >= localRhs.size() )
      continue;
    RAJA::atomicAdd< parallelDeviceAtomic >( &localRhs[dof+0], traction[0] * nodalWeight[a] );
    RAJA::atomicAdd< parallelDeviceAtomic >( &localRhs[dof+1], traction[1] * nodalWeight[a] );
    RAJA::atomicAdd< parallelDeviceAtomic >( &localRhs[dof+2], traction[2] * nodalWeight[a] );
  }
}

/**
 * @brief Integrate a uniform traction over a general polygon face.
 *
 * This is a fallback for faces that are not triangles or quadrilaterals.
 *
 * @param traction The traction vector [3] to integrate.
 * @param faceArea The area of the face.
 * @param numNodes The number of nodes of the face.
 * @param faceToNodeMap The face-to-node connectivity map.
 * @param kf       The face index.
 * @param blockLocalDofNumber Block-local DOF numbering.
 * @param dofRankOffset The rank offset for DOFs.
 * @param localRhs The local RHS vector to assemble into.
 */
GEOS_HOST_DEVICE
inline
void integrateFaceTractionPolygon( real64 const ( &traction )[3],
                                   real64 const faceArea,
                                   localIndex const numNodes,
                                   ArrayOfArraysView< localIndex const > const & faceToNodeMap,
                                   localIndex const kf,
                                   arrayView1d< globalIndex const > const & blockLocalDofNumber,
                                   globalIndex const dofRankOffset,
                                   arrayView1d< real64 > const & localRhs )
{
  real64 const w = faceArea / numNodes;
  for( localIndex a = 0; a < numNodes; ++a )
  {
    localIndex const dof = blockLocalDofNumber[ faceToNodeMap( kf, a ) ] - dofRankOffset;
    if( dof < 0 || dof >= localRhs.size() )
      continue;
    RAJA::atomicAdd< parallelDeviceAtomic >( &localRhs[dof+0], traction[0] * w );
    RAJA::atomicAdd< parallelDeviceAtomic >( &localRhs[dof+1], traction[1] * w );
    RAJA::atomicAdd< parallelDeviceAtomic >( &localRhs[dof+2], traction[2] * w );
  }
}

/**
 * @brief Dispatch face traction integration to the appropriate FE formulation.
 *
 * Selects the FE face type based on numNodes:
 * - 3 nodes: H1_TriangleFace_Lagrange1_Gauss1_impl
 * - 4 nodes: H1_QuadrilateralFace_Lagrange1_GaussLegendre2_impl
 * - Other:   Equal distribution fallback
 *
 * @param traction The traction vector [3].
 * @param numNodes The number of nodes of the face.
 * @param faceArea The area of the face (only used for polygon fallback).
 * @param faceToNodeMap The face-to-node connectivity map.
 * @param kf       The face index.
 * @param nodePositions The reference positions of all nodes.
 * @param blockLocalDofNumber Block-local DOF numbering.
 * @param dofRankOffset The rank offset for DOFs.
 * @param localRhs The local RHS vector to assemble into.
 */
GEOS_HOST_DEVICE
inline
void assembleFaceTraction( real64 const ( &traction )[3],
                           localIndex const numNodes,
                           real64 const faceArea,
                           ArrayOfArraysView< localIndex const > const & faceToNodeMap,
                           localIndex const kf,
                           arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePositions,
                           arrayView1d< globalIndex const > const & blockLocalDofNumber,
                           globalIndex const dofRankOffset,
                           arrayView1d< real64 > const & localRhs )
{
  using TriFace = finiteElement::H1_TriangleFace_Lagrange1_Gauss1_impl;
  using QuadFace = finiteElement::H1_QuadrilateralFace_Lagrange1_GaussLegendre2_impl;

  if( numNodes == 3 )
  {
    // Get the permutation mapping from FE node ordering to mesh face-local node ordering.
    int permutation[TriFace::numNodes];
    TriFace::getPermutation( permutation );

    // Gather coordinates in FE ordering: xFace[feNode] = coords of mesh node permutation[feNode]
    real64 xFace[TriFace::numNodes][3];
    for( localIndex a = 0; a < TriFace::numNodes; ++a )
    {
      localIndex const nodeIdx = faceToNodeMap( kf, permutation[a] );
      xFace[a][0] = nodePositions( nodeIdx, 0 );
      xFace[a][1] = nodePositions( nodeIdx, 1 );
      xFace[a][2] = nodePositions( nodeIdx, 2 );
    }
    integrateFaceTraction< TriFace >( traction, xFace, permutation, faceToNodeMap, kf,
                                      blockLocalDofNumber, dofRankOffset, localRhs );
  }
  else if( numNodes == 4 )
  {
    // Get the permutation mapping from FE node ordering to mesh face-local node ordering.
    // For quads: FE uses Z-ordering (tensor product), mesh uses CCW ordering.
    // getPermutation returns {0, 1, 3, 2} meaning FE node 2 → mesh node 3 and vice versa.
    int permutation[QuadFace::numNodes];
    QuadFace::getPermutation( permutation );

    // Gather coordinates in FE ordering: xFace[feNode] = coords of mesh node permutation[feNode]
    real64 xFace[QuadFace::numNodes][3];
    for( localIndex a = 0; a < QuadFace::numNodes; ++a )
    {
      localIndex const nodeIdx = faceToNodeMap( kf, permutation[a] );
      xFace[a][0] = nodePositions( nodeIdx, 0 );
      xFace[a][1] = nodePositions( nodeIdx, 1 );
      xFace[a][2] = nodePositions( nodeIdx, 2 );
    }
    integrateFaceTraction< QuadFace >( traction, xFace, permutation, faceToNodeMap, kf,
                                       blockLocalDofNumber, dofRankOffset, localRhs );
  }
  else
  {
    integrateFaceTractionPolygon( traction, faceArea, numNodes, faceToNodeMap, kf,
                                  blockLocalDofNumber, dofRankOffset, localRhs );
  }
}


void TractionBoundaryCondition::launch( real64 const time,
                                        arrayView1d< globalIndex const > const blockLocalDofNumber,
                                        globalIndex const dofRankOffset,
                                        FaceManager const & faceManager,
                                        arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const nodePositions,
                                        SortedArrayView< localIndex const > const & targetSet,
                                        arrayView1d< real64 > const & localRhs ) const
{
  arrayView1d< real64 const > const faceArea  = faceManager.faceArea();
  arrayView2d< real64 const > const faceNormal  = faceManager.faceNormal();
  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  FunctionManager const & functionManager = FunctionManager::getInstance();

  string const & functionName = this->getFunctionName();

  globalIndex_array nodeDOF;
  real64_array nodeRHS;

  R1Tensor const direction = this->getDirection();
  TractionType const tractionType = m_tractionType;
  R2SymTensor inputStress = m_inputStress;

  real64 tractionMagnitude0;
  array1d< real64 > tractionMagnitudeArray( targetSet.size() );
  bool spatialFunction = false;

  if( functionName.empty() )
  {
    tractionMagnitude0 = this->getScale();
  }
  else
  {
    FunctionBase const & function = functionManager.getGroup< FunctionBase >( functionName );
    if( function.isFunctionOfTime() == 2 )
    {
      tractionMagnitude0 = this->getScale() * function.evaluate( &time );
    }
    else
    {
      tractionMagnitudeArray.setName( getName() + " function results" );
      function.evaluate( faceManager, time, targetSet, tractionMagnitudeArray );
      spatialFunction = true;
      tractionMagnitude0 = 1e99;
    }
  }

  arrayView1d< real64 const > const tractionMagnitudeArrayView = tractionMagnitudeArray;

  {
    integer const nodalScaleFlag = m_nodalScaleFlag;
    auto const scaleSet = m_scaleSet.toViewConst();

    forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const i )
    {
      localIndex const kf = targetSet[ i ];
      localIndex const numNodes = faceToNodeMap.sizeOfArray( kf );

      real64 const nodalScale = ( nodalScaleFlag == 1 ) ? scaleSet[i] : 1.0;

      real64 const tractionMagnitude = spatialFunction ? tractionMagnitudeArrayView[i] : (tractionMagnitude0 * nodalScale);

      real64 traction[3] = { 0 };
      if( tractionType == TractionType::vector )
      {
        traction[0] = tractionMagnitude * direction[0];
        traction[1] = tractionMagnitude * direction[1];
        traction[2] = tractionMagnitude * direction[2];
      }
      else
      {
        real64 const temp[3] = { tractionMagnitude * faceNormal( kf, 0 ),
                                 tractionMagnitude * faceNormal( kf, 1 ),
                                 tractionMagnitude * faceNormal( kf, 2 ) };

        if( tractionType == TractionType::normal )
        {
          traction[0] = temp[0];
          traction[1] = temp[1];
          traction[2] = temp[2];
        }
        else if( tractionType == TractionType::stress )
        {
          traction[0] = inputStress[0] * temp[0] + inputStress[5] * temp[1] + inputStress[4] * temp[2];
          traction[1] = inputStress[5] * temp[0] + inputStress[1] * temp[1] + inputStress[3] * temp[2];
          traction[2] = inputStress[4] * temp[0] + inputStress[3] * temp[1] + inputStress[2] * temp[2];
        }
      }

      // Dispatch face traction integration to the appropriate FE formulation
      assembleFaceTraction( traction, numNodes, faceArea[kf], faceToNodeMap, kf,
                            nodePositions, blockLocalDofNumber, dofRankOffset, localRhs );
    } );
  }
}

void TractionBoundaryCondition::reinitScaleSet( FaceManager const & faceManager,
                                                SortedArrayView< localIndex const > const & targetSet,
                                                arrayView1d< real64 const > const nodalScaleSet )
{
  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  localIndex const faceSize = targetSet.size();

  m_scaleSet.resize( faceSize );

  auto const scaleSet = m_scaleSet.toView();

  // Loop over targetSet to assign damage values to m_scaleSet
  forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const i )
  {
    localIndex const kf = targetSet[ i ];
    localIndex const numNodes = faceToNodeMap.sizeOfArray( kf );

    real64 faceScale = 0.0;

    for( localIndex a=0; a<numNodes; ++a )
    {
      faceScale  += nodalScaleSet[ faceToNodeMap( kf, a ) ];
    }

    scaleSet[i] = LvArray::math::min( 1.0, (1.0 - faceScale/numNodes)*(1.0 - faceScale/numNodes) );
  } );
}

REGISTER_CATALOG_ENTRY( FieldSpecificationABC, TractionBoundaryCondition, string const &, Group * const )


} /* namespace geos */
