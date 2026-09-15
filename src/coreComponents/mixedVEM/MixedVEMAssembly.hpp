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
 * @file MixedVEMAssembly.hpp
 *
 * Mesh adaptors and scatter helpers for the mixed VEM element operators.
 *
 * The face unknowns are the moments (5)-(6) in the frame of f, shared by both neighbours; the
 * outward sign s_{E,f} is folded into B_E and P_E.
 */

#ifndef GEOS_MIXEDVEM_MIXEDVEMASSEMBLY_HPP_
#define GEOS_MIXEDVEM_MIXEDVEMASSEMBLY_HPP_

#include "common/DataLayouts.hpp"
#include "mixedVEM/HybridMixedVEM.hpp"

namespace geos
{

namespace mixedVEM
{

/**
 * @struct IndexedNodeCoordinates
 * @brief Adaptor presenting an indirect node list as a (node, dim) -> real64 callable.
 */
struct IndexedNodeCoordinates
{
  /// The reference positions of all mesh nodes.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > nodePositions;

  /// The nodes of the entity being read, in order.
  arraySlice1d< localIndex const > nodeList;

  /**
   * @brief Coordinate accessor.
   * @param[in] k the local node index
   * @param[in] i the dimension
   * @return the coordinate
   */
  GEOS_HOST_DEVICE
  real64 operator()( integer const k, integer const i ) const
  { return nodePositions( nodeList[k], i ); }
};

/**
 * @brief Build the oriented face geometry and the moments of one element.
 * @tparam FACE_NODES the type of the face to node map
 * @param[in] nodePositions the reference positions of the mesh nodes
 * @param[in] faceToNodes the map from a face to its nodes, in loop order
 * @param[in] faceNormals the face-intrinsic unit normals
 * @param[in] elemToFaces the faces of the element
 * @param[in] numFaces the number of faces of the element
 * @param[in] elemCenter the point x_E
 * @param[out] faceGeom the geometry of each face, oriented for this element
 * @param[out] moments the element moments
 */
template< typename FACE_NODES >
GEOS_HOST_DEVICE
void buildElementGeometry( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePositions,
                           FACE_NODES const & faceToNodes,
                           arrayView2d< real64 const > const & faceNormals,
                           arraySlice1d< localIndex const > const & elemToFaces,
                           integer const numFaces,
                           real64 const (&elemCenter)[3],
                           FaceGeometry * const faceGeom,
                           ElementMoments & moments )
{
  initialize( moments );

  for( integer lf = 0; lf < numFaces; ++lf )
  {
    localIndex const faceIndex = elemToFaces[lf];

    IndexedNodeCoordinates const X { nodePositions, faceToNodes[faceIndex] };
    integer const numNodes = LvArray::integerConversion< integer >( faceToNodes.sizeOfArray( faceIndex ) );

    real64 const normal[3] = { faceNormals( faceIndex, 0 ),
                               faceNormals( faceIndex, 1 ),
                               faceNormals( faceIndex, 2 ) };

    computeFaceGeometry( X, numNodes, normal, faceGeom[lf] );
    orientFace( elemCenter, faceGeom[lf] );

    accumulateElementMoments( X, numNodes, faceGeom[lf], elemCenter, moments );
  }
}

/**
 * @brief Collect the global stress degrees of freedom of an element.
 * @param[in] elemToFaces the faces of the element
 * @param[in] numFaces the number of faces of the element
 * @param[in] faceDofNumber the first stress degree of freedom of each face
 * @param[out] dofIndices the 6 numFaces global indices, grouped by face
 *
 * The layout matches the column ordering of B_E and M_E, so the scatter below is a
 * plain index lookup with no permutation.
 */
GEOS_HOST_DEVICE
inline void gatherStressDofIndices( arraySlice1d< localIndex const > const & elemToFaces,
                                    integer const numFaces,
                                    arrayView1d< globalIndex const > const & faceDofNumber,
                                    globalIndex * const dofIndices )
{
  for( integer lf = 0; lf < numFaces; ++lf )
  {
    globalIndex const offset = faceDofNumber[elemToFaces[lf]];
    for( integer j = 0; j < NUM_FACE_DOF; ++j )
    {
      dofIndices[NUM_FACE_DOF * lf + j] = offset + j;
    }
  }
}

/**
 * @brief Scatter one element block into the global saddle-point matrix.
 * @tparam ATOMIC the atomic policy used by the scatter
 * @param[in] localMatrix the global matrix owned by this rank
 * @param[in] rankOffset the first global row owned by this rank
 * @param[in] stressDofIndices the global stress degrees of freedom, grouped by face
 * @param[in] numStressDof the number of stress degrees of freedom, 6 numFaces
 * @param[in] dispDofIndices the six global displacement degrees of freedom
 * @param[in] complianceMatrix the matrix M_E
 * @param[in] divergence the matrix B_E
 * @param[out] rowBuffer a scratch buffer of length @p numStressDof
 *
 * Rows of M_E and of B_E are contiguous, so each is scattered with one call per row;
 * the B_E^T block is the only strided read and is packed into a six entry buffer.
 *
 * The RM(E) rows carry -B_E, the second equation of (19), so the Schur complement B M^{-1} B^T is
 * positive definite.
 */
template< typename ATOMIC >
GEOS_HOST_DEVICE
void addElementToMatrix( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                         globalIndex const rankOffset,
                         globalIndex const * const stressDofIndices,
                         integer const numStressDof,
                         globalIndex const * const dispDofIndices,
                         MatrixSliceConst const & complianceMatrix,
                         MatrixSliceConst const & divergence,
                         real64 * const rowBuffer )
{
  for( integer i = 0; i < numStressDof; ++i )
  {
    localIndex const localRow = LvArray::integerConversion< localIndex >( stressDofIndices[i] - rankOffset );
    if( localRow < 0 || localRow >= localMatrix.numRows() )
    {
      continue;
    }

    localMatrix.template addToRowBinarySearchUnsorted< ATOMIC >( localRow,
                                                                 stressDofIndices,
                                                                 &complianceMatrix( i, 0 ),
                                                                 numStressDof );

    real64 divergenceColumn[NUM_RM_DOF];
    for( integer k = 0; k < NUM_RM_DOF; ++k )
    {
      divergenceColumn[k] = divergence( k, i );
    }

    localMatrix.template addToRowBinarySearchUnsorted< ATOMIC >( localRow,
                                                                 dispDofIndices,
                                                                 divergenceColumn,
                                                                 NUM_RM_DOF );
  }

  for( integer k = 0; k < NUM_RM_DOF; ++k )
  {
    localIndex const localRow = LvArray::integerConversion< localIndex >( dispDofIndices[k] - rankOffset );
    if( localRow < 0 || localRow >= localMatrix.numRows() )
    {
      continue;
    }

    for( integer j = 0; j < numStressDof; ++j )
    {
      rowBuffer[j] = -divergence( k, j );
    }

    localMatrix.template addToRowBinarySearchUnsorted< ATOMIC >( localRow,
                                                                 stressDofIndices,
                                                                 rowBuffer,
                                                                 numStressDof );
  }
}

/**
 * @struct FaceConstraint
 * @brief Which traction modes of a boundary face are prescribed, and with what data.
 */
struct FaceConstraint
{
  /// True when the mode is an essential unknown fixed by a prescribed traction.
  bool essential[NUM_FACE_DOF];

  /// Value the degree of freedom takes when @p essential.
  real64 value[NUM_FACE_DOF];

  /// Right hand side contribution when the mode is free, from the prescribed displacement.
  real64 load[NUM_FACE_DOF];
};

/**
 * @brief Decide the boundary role of each traction mode of a face, component by component.
 * @param[in] geom the face geometry, oriented for the adjacent element
 * @param[in] mask bit k set when displacement component k is prescribed
 * @param[in] displacement the prescribed displacement, components that are not masked are ignored
 * @param[in] traction the prescribed traction, components that are masked are ignored
 * @param[out] constraint the resulting per mode roles and data
 *
 * Constant modes psi_k follow component k: a prescribed displacement is natural and loads the right
 * hand side, a prescribed traction psi is essential and fixes the mode at |f| s_{E,f} psi_k.
 * Rotational mode 3+k is essential when e_k ^ span{t1, t2} carries no prescribed displacement.
 */
GEOS_HOST_DEVICE
inline void computeFaceConstraint( FaceGeometry const & geom,
                                   integer const mask,
                                   real64 const (&displacement)[3],
                                   real64 const (&traction)[3],
                                   FaceConstraint & constraint )
{
  // a direction is carried by a prescribed traction when it has no component along a
  // masked axis
  auto inTractionSpan = [mask]( real64 const (&v)[3] )
  {
    for( integer k = 0; k < 3; ++k )
    {
      if( ( ( mask >> k ) & 1 ) && LvArray::math::abs( v[k] ) > 1e-10 )
      {
        return false;
      }
    }
    return true;
  };

  for( integer k = 0; k < NUM_FACE_DOF; ++k )
  {
    constraint.essential[k] = false;
    constraint.value[k] = 0.0;
    constraint.load[k] = 0.0;
  }

  for( integer k = 0; k < 3; ++k )
  {
    if( ( mask >> k ) & 1 )
    {
      constraint.load[k] = geom.outwardSign * displacement[k];
    }
    else
    {
      constraint.essential[k] = true;
      constraint.value[k] = geom.outwardSign * geom.area * traction[k];
    }
  }

  // mode 3+k takes its values in e_k ^ span{t1, t2}, so it is a traction mode exactly when
  // that plane avoids every direction whose displacement is prescribed
  for( integer k = 0; k < 3; ++k )
  {
    real64 v1[3], v2[3];
    axisCross( k, geom.t1, v1 );
    axisCross( k, geom.t2, v2 );

    constraint.essential[3 + k] = inTractionSpan( v1 ) && inTractionSpan( v2 );
  }
}

/**
 * @brief Collect the interface multiplier degrees of freedom of an element.
 * @param[in] elemToFaces the faces of the element
 * @param[in] numFaces the number of faces of the element
 * @param[in] faceDofNumber the first multiplier degree of freedom of each face, negative
 *   on faces carrying no unknown, that is the Dirichlet part of the boundary
 * @param[out] dofIndices the 6 numFaces global indices, negative where eliminated
 */
GEOS_HOST_DEVICE
inline void gatherMultiplierDofIndices( arraySlice1d< localIndex const > const & elemToFaces,
                                        integer const numFaces,
                                        arrayView1d< globalIndex const > const & faceDofNumber,
                                        globalIndex * const dofIndices )
{
  for( integer lf = 0; lf < numFaces; ++lf )
  {
    globalIndex const offset = faceDofNumber[elemToFaces[lf]];

    for( integer j = 0; j < NUM_FACE_DOF; ++j )
    {
      dofIndices[NUM_FACE_DOF * lf + j] = ( offset < 0 ) ? -1 : offset + j;
    }
  }
}

/**
 * @brief Scatter one element contribution C_E S_E C_E^T into the interface operator.
 * @tparam ATOMIC the atomic policy used by the scatter
 * @param[in] localMatrix the global interface matrix owned by this rank
 * @param[in] rankOffset the first global row owned by this rank
 * @param[in] dofIndices the multiplier degrees of freedom, negative where eliminated
 * @param[in] numStressDof the number of stress degrees of freedom, 6 numFaces
 * @param[in] contribution the matrix C_E S_E C_E^T
 * @param[out] packedColumns a scratch buffer of length @p numStressDof
 * @param[out] packedOrder a scratch buffer of length @p numStressDof
 * @param[out] packedValues a scratch buffer of length @p numStressDof
 *
 * H = sum_E C_E S_E C_E^T is symmetric positive semidefinite, and definite once the
 * Dirichlet multipliers are dropped, so the assembled system is amenable to conjugate
 * gradients rather than to an indefinite solver.
 */
template< typename ATOMIC >
GEOS_HOST_DEVICE
void addElementToInterfaceMatrix( CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                  globalIndex const rankOffset,
                                  globalIndex const * const dofIndices,
                                  integer const numStressDof,
                                  MatrixSliceConst const & contribution,
                                  globalIndex * const packedColumns,
                                  real64 * const packedValues,
                                  integer * const packedOrder )
{
  // every row of the block has the same column set, so pack and sort it once
  integer numPacked = 0;
  for( integer j = 0; j < numStressDof; ++j )
  {
    if( dofIndices[j] < 0 )
    {
      continue;
    }

    globalIndex const column = dofIndices[j];
    integer p = numPacked++;
    while( p > 0 && packedColumns[p - 1] > column )
    {
      packedColumns[p] = packedColumns[p - 1];
      packedOrder[p] = packedOrder[p - 1];
      --p;
    }
    packedColumns[p] = column;
    packedOrder[p] = j;
  }

  for( integer i = 0; i < numStressDof; ++i )
  {
    if( dofIndices[i] < 0 )
    {
      continue;
    }

    localIndex const localRow = LvArray::integerConversion< localIndex >( dofIndices[i] - rankOffset );
    if( localRow < 0 || localRow >= localMatrix.numRows() )
    {
      continue;
    }

    for( integer p = 0; p < numPacked; ++p )
    {
      packedValues[p] = contribution( i, packedOrder[p] );
    }

    // sorted columns let the row be merged in one pass instead of a search per entry
    localMatrix.template addToRow< ATOMIC >( localRow, packedColumns, packedValues, numPacked );
  }
}

} // namespace mixedVEM

} // namespace geos

#endif // GEOS_MIXEDVEM_MIXEDVEMASSEMBLY_HPP_
