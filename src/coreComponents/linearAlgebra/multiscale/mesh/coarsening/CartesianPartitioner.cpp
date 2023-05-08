/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CartesianPartitioner.cpp
 */

#include "CartesianPartitioner.hpp"

#include "common/MpiWrapper.hpp"
#include "linearAlgebra/multiscale/mesh/MeshData.hpp"

namespace geos
{
namespace multiscale
{

localIndex CartesianPartitioner::generate( MeshLevel const & mesh,
                                           arrayView1d< localIndex > const & partition )
{
  GEOS_MARK_FUNCTION;

  MeshObjectManager const & cellManager = mesh.cellManager();

  GEOS_THROW_IF_NE_MSG( m_params.ratio.size(), 3,
                        "Cartesian partitioning requires 3 coarsening ratio values",
                        InputError );
  GEOS_THROW_IF( !cellManager.hasField< fields::StructuredIndex >(),
                 "Cartesian partitioning is only compatible with InternalMesh",
                 InputError );

  arrayView2d< integer const > const cartIndex =
    cellManager.getField< fields::StructuredIndex >().toViewConst();

  // Find the range of locally owned cell cartesian indices
  integer constexpr minIdx = std::numeric_limits< integer >::min();
  integer constexpr maxIdx = std::numeric_limits< integer >::max();

  RAJA::ReduceMin< parallelHostReduce, integer > lo0( maxIdx ), lo1( maxIdx ), lo2( maxIdx );
  RAJA::ReduceMax< parallelHostReduce, integer > hi0( minIdx ), hi1( minIdx ), hi2( minIdx );

  forAll< parallelHostPolicy >( cellManager.numOwnedObjects(), [=]( localIndex const i )
  {
    lo0.min( cartIndex[i][0] ); lo1.min( cartIndex[i][1] ); lo2.min( cartIndex[i][2] );
    hi0.max( cartIndex[i][0] ); hi1.max( cartIndex[i][1] ); hi2.max( cartIndex[i][2] );
  } );

  integer loCartIndex[3] = { lo0.get(), lo1.get(), lo2.get() };
  integer hiCartIndex[3] = { hi0.get(), hi1.get(), hi2.get() };

  // Special treatment for ranks that don't have a piece of the mesh
  if( cellManager.numOwnedObjects() == 0 )
  {
    for( int dim = 0; dim < 3; ++dim )
    {
      loCartIndex[dim] = 0;
      hiCartIndex[dim] = -1;
    }
  }

  // When semi-coarsening is enabled, only coarsen in z-direction,
  // until every rank has reduced its mesh to 1 cell in z-direction.
  if( m_params.structured.semicoarsening )
  {
    integer const numCellsZ = hiCartIndex[2] - loCartIndex[2] + 1;
    integer const maxNumCellsZ = MpiWrapper::max( numCellsZ );
    if( maxNumCellsZ > 1 )
    {
      m_params.ratio[0] = 1.0;
      m_params.ratio[1] = 1.0;
    }
  }

  // Compute cartesian sizes of coarse grid
  integer ratio[3]{ 1, 1, 1 };
  for( int dim = 0; dim < 3; ++dim )
  {
    integer const numCellsDim = hiCartIndex[dim] - loCartIndex[dim] + 1;
    m_numPart[dim] = static_cast< integer >( std::ceil( numCellsDim / m_params.ratio[dim] ) ); // rounded up
    if( m_numPart[dim] > 0 )
    {
      ratio[dim] = ( numCellsDim + m_numPart[dim] - 1 ) / m_numPart[dim]; // rounded up
    }
  }

  localIndex const partStride[3] = { 1, m_numPart[0], m_numPart[0] * m_numPart[1] };

  // Compute cartesian coarse cell indices
  forAll< parallelHostPolicy >( cellManager.numOwnedObjects(), [=, &loCartIndex]( localIndex const i )
  {
    localIndex part = 0;
    for( int dim = 0; dim < 3; ++dim )
    {
      integer const cartIndexCoarse = ( cartIndex[i][dim] - loCartIndex[dim] ) / ratio[dim];
      part += cartIndexCoarse * partStride[dim];
    }
    partition[i] = part;
  } );

  return m_numPart[0] * m_numPart[1] * m_numPart[2];
}

void CartesianPartitioner::setCoarseData( multiscale::MeshLevel & coarseMesh ) const
{
  array2d< integer > & cartIndex =
    coarseMesh.cellManager().registerField< fields::StructuredIndex >( {} ).reference();
  cartIndex.resizeDimension< 1 >( 3 );

  localIndex const cartPartStride[3] = { 1, m_numPart[0], m_numPart[0] * m_numPart[1] };

  forAll< parallelHostPolicy >( coarseMesh.cellManager().numOwnedObjects(),
                                [cartIndex = cartIndex.toView(), cartPartStride]( localIndex const i )
  {
    localIndex cellIndex = i;
    for( int dim = 2; dim >= 0; --dim )
    {
      cartIndex[i][dim] = LvArray::integerConversion< integer >( cellIndex / cartPartStride[dim] );
      cellIndex = cellIndex % cartPartStride[dim];
    }
  } );
}

} // namespace multiscale
} // namespace geos
