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
 * @file MultivariableTableFunction.cpp
 */

#include "MultivariableTableFunction.hpp"

#include "common/DataTypes.hpp"
#include <algorithm>

namespace geosx
{

using namespace dataRepository;

MultivariableTableFunction::MultivariableTableFunction( const string & name,
                                                        Group * const parent ):
  FunctionBase( name, parent )
{}

void MultivariableTableFunction::initializeFunctionFromFile( string const & filename )
{
  std::ifstream file( filename.c_str() );
  GEOSX_THROW_IF( !file, catalogName() << " " << getName() << ": could not read input file " << filename, InputError );

  integer numDims, numOps;
  globalIndex numPointsTotal = 1;
  real64_array axisMinimums, axisMaximums;
  integer_array axisPoints;

  // 1. Read numDims and numOps

  GEOSX_THROW_IF( file.eof(), "Table file is shorter than expected", std::runtime_error );
  file >> numDims >> numOps;

  // assume no more than 10 dimensions
  GEOSX_ASSERT_GE( numDims, 1 );
  GEOSX_ASSERT_GT( 10, numDims );

  // assume no more than 100 operators
  GEOSX_ASSERT_GE( numOps, 1 );
  GEOSX_ASSERT_GT( 100, numDims );

  axisMinimums.resize( numDims );
  axisMaximums.resize( numDims );
  axisPoints.resize( numDims );

  // 2. Read axis parameters

  for( integer i = 0; i < numDims; i++ )
  {
    GEOSX_THROW_IF( file.eof(), "Table file is shorter than expected", std::runtime_error );
    file >> axisPoints[i] >> axisMinimums[i] >> axisMaximums[i];
    GEOSX_ASSERT_GT( axisPoints[i], 1 );
    numPointsTotal *= axisPoints[i];
  }
  m_pointData.resize( numPointsTotal * numOps );

  // 3. Read table values
  for( auto i = 0; i < numPointsTotal; i++ )
  {
    for( auto j = 0; j < numOps; j++ )
    {
      GEOSX_THROW_IF( file.eof(), "Table file is shorter than expected", std::runtime_error );
      file >> m_pointData[i * numOps + j];
    }
  }

  file.close();

  setTableCoordinates( numDims, numOps, axisMinimums, axisMaximums, axisPoints );
  initializeFunction();
}


void MultivariableTableFunction::setTableCoordinates( integer const numDims,
                                                      integer const numOps,
                                                      real64_array const & axisMinimums,
                                                      real64_array const & axisMaximums,
                                                      integer_array const & axisPoints )
{
  m_numDims = numDims;
  m_numOps = numOps;
  m_numVerts =  1 << numDims;

  m_axisMinimums = axisMinimums;
  m_axisMaximums = axisMaximums;
  m_axisPoints = axisPoints;
}

void MultivariableTableFunction::setTableValues( real64_array values )
{
  m_pointData = std::move( values );
}
void MultivariableTableFunction::getHypercubePoints( globalIndex const hypercubeIndex, globalIndex_array & hypercubePoints ) const
{
  auto remainder = hypercubeIndex;
  auto pwr = m_numVerts;

  for( auto j = 0; j < m_numVerts; ++j )
  {
    hypercubePoints[j] = 0;
  }

  for( auto i = 0; i < m_numDims; ++i )
  {

    integer axis_idx = remainder / m_axisHypercubeMults[i];
    remainder = remainder % m_axisHypercubeMults[i];

    pwr /= 2;

    for( auto j = 0; j < m_numVerts; ++j )
    {
      auto zero_or_one = (j / pwr) % 2;
      hypercubePoints[j] += (axis_idx + zero_or_one) * m_axisPointMults[i];
    }
  }
}


void MultivariableTableFunction::initializeFunction()
{
  // check input


  GEOSX_ASSERT_EQ( m_numDims, m_axisMinimums.size());
  GEOSX_ASSERT_EQ( m_numDims, m_axisMaximums.size());
  GEOSX_ASSERT_EQ( m_numDims, m_axisPoints.size());

  m_axisSteps.resize( m_numDims );
  m_axisStepInvs.resize( m_numDims );
  m_axisPointMults.resize( m_numDims );
  m_axisHypercubeMults.resize( m_numDims );

  // compute service data

  for( integer dim = 0; dim < m_numDims; dim++ )
  {
    m_axisSteps[dim] = (m_axisMaximums[dim] - m_axisMinimums[dim]) / (m_axisPoints[dim] - 1);
    m_axisStepInvs[dim] = 1 / m_axisSteps[dim];
  }

  m_axisPointMults[m_numDims - 1] = 1;
  m_axisHypercubeMults[m_numDims - 1] = 1;
  for( integer dim = m_numDims - 2; dim >= 0; --dim )
  {
    m_axisPointMults[dim] = m_axisPointMults[dim + 1] * m_axisPoints[dim + 1];
    m_axisHypercubeMults[dim] = m_axisHypercubeMults[dim + 1] * (m_axisPoints[dim + 1] - 1);
  }

  // check for point index overflow
  // fp type is intentional - to prevent overflow during computation and detect it later
  real64 numTablePoints = 1.0;
  globalIndex numTableHypercubes = 1;
  for( int dim = 0; dim < m_numDims; dim++ )
  {
    numTablePoints *= m_axisPoints[dim];
    numTableHypercubes *= m_axisPoints[dim] - 1;
  }

  GEOSX_ASSERT_GT_MSG( std::numeric_limits< globalIndex >::max(), numTablePoints, "MultivariableTableFunction: too many points to index, reduce axisPPoints values" );

  // check is point data size is correct
  GEOSX_ASSERT_EQ( globalIndex( numTablePoints ) * m_numOps, m_pointData.size());


  // initialize hypercube data storage


  m_hypercubeData.resize( numTableHypercubes * m_numVerts * m_numOps );
  globalIndex_array points( m_numVerts );

  // fill each hypercube directly on the device with corresponding data from m_pointData
  for( auto i = 0; i < numTableHypercubes; i++ )
  {
    getHypercubePoints( i, points );

    for( auto j = 0; j < m_numVerts; ++j )
    {
      std::copy( m_pointData.begin() + points[j] * m_numOps,
                 m_pointData.begin() + (points[j] + 1) * m_numOps,
                 m_hypercubeData.begin() + m_numOps * (i * m_numVerts + j));
    }
  }

}

REGISTER_CATALOG_ENTRY( FunctionBase, MultivariableTableFunction, string const &, Group * const )

} /* namespace ANST */
