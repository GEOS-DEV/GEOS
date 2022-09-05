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
 * @file SolidMechanicsStatistics.cpp
 */

#include "SolidMechanicsStatistics.hpp"

#include "common/MpiWrapper.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"

namespace geosx
{

using namespace constitutive;
using namespace dataRepository;

SolidMechanicsStatistics::SolidMechanicsStatistics( const string & name,
                                                    Group * const parent ):
  Base( name, parent )
{}

void SolidMechanicsStatistics::registerDataOnMesh( Group & meshBodies )
{
  // the fields have to be registered in "registerDataOnMesh" (and not later)
  // otherwise they cannot be targeted by TimeHistory

  // for now, this guard is needed to avoid breaking the xml schema generation
  if( m_solver == nullptr )
  {
    return;
  }

  m_solver->forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                              MeshLevel & mesh,
                                                              arrayView1d< string const > const & )
  {
    NodeManager & nodeManager = mesh.getNodeManager();
    nodeManager.registerWrapper< NodeStatistics >( viewKeyStruct::nodeStatisticsString() ).
      setRestartFlags( RestartFlags::NO_WRITE );
    nodeManager.excludeWrappersFromPacking( { viewKeyStruct::nodeStatisticsString() } );
    NodeStatistics & nodeStatistics = nodeManager.getReference< NodeStatistics >( viewKeyStruct::nodeStatisticsString() );

    nodeStatistics.minDisplacement.resizeDimension< 0 >( 3 );
    nodeStatistics.maxDisplacement.resizeDimension< 0 >( 3 );
  } );
}

bool SolidMechanicsStatistics::execute( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                        real64 const GEOSX_UNUSED_PARAM( dt ),
                                        integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                                        integer const GEOSX_UNUSED_PARAM( eventCounter ),
                                        real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                                        DomainPartition & domain )
{
  m_solver->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                          MeshLevel & mesh,
                                                                          arrayView1d< string const > const & )
  {
    computeNodeStatistics( mesh );
  } );
  return false;
}

void SolidMechanicsStatistics::computeNodeStatistics( MeshLevel & mesh ) const
{
  GEOSX_MARK_FUNCTION;

  // Step 1: increment the min/max quantities

  NodeManager & nodeManager = mesh.getNodeManager();
  arrayView1d< integer const > const ghostRank = nodeManager.ghostRank();
  arrayView2d< real64, nodes::TOTAL_DISPLACEMENT_USD > const & u = nodeManager.totalDisplacement();

  RAJA::ReduceMax< parallelDeviceReduce, real64 > maxDispX( -LvArray::NumericLimits< real64 >::max );
  RAJA::ReduceMax< parallelDeviceReduce, real64 > maxDispY( -LvArray::NumericLimits< real64 >::max );
  RAJA::ReduceMax< parallelDeviceReduce, real64 > maxDispZ( -LvArray::NumericLimits< real64 >::max );
  RAJA::ReduceMin< parallelDeviceReduce, real64 > minDispX( LvArray::NumericLimits< real64 >::max );
  RAJA::ReduceMin< parallelDeviceReduce, real64 > minDispY( LvArray::NumericLimits< real64 >::max );
  RAJA::ReduceMin< parallelDeviceReduce, real64 > minDispZ( LvArray::NumericLimits< real64 >::max );

  forAll< parallelDevicePolicy<> >( nodeManager.size(), [u,
                                                         ghostRank,
                                                         maxDispX,
                                                         maxDispY,
                                                         maxDispZ,
                                                         minDispX,
                                                         minDispY,
                                                         minDispZ] GEOSX_HOST_DEVICE ( localIndex const a )
  {
    if( ghostRank[a] < 0 )
    {
      maxDispX.max( u[a][0] );
      maxDispY.max( u[a][1] );
      maxDispZ.max( u[a][2] );
      minDispX.min( u[a][0] );
      minDispY.min( u[a][1] );
      minDispZ.min( u[a][2] );
    }
  } );

  // Step 2: synchronize the results over the MPI ranks

  NodeStatistics & nodeStatistics = nodeManager.getReference< NodeStatistics >( viewKeyStruct::nodeStatisticsString() );

  nodeStatistics.maxDisplacement[0] = maxDispX.get();
  nodeStatistics.maxDisplacement[1] = maxDispY.get();
  nodeStatistics.maxDisplacement[2] = maxDispZ.get();
  nodeStatistics.minDisplacement[0] = minDispX.get();
  nodeStatistics.minDisplacement[1] = minDispY.get();
  nodeStatistics.minDisplacement[2] = minDispZ.get();

  MpiWrapper::allReduce( nodeStatistics.maxDisplacement.data(),
                         nodeStatistics.maxDisplacement.data(),
                         3,
                         MpiWrapper::getMpiOp( MpiWrapper::Reduction::Max ),
                         MPI_COMM_GEOSX );

  MpiWrapper::allReduce( nodeStatistics.minDisplacement.data(),
                         nodeStatistics.minDisplacement.data(),
                         3,
                         MpiWrapper::getMpiOp( MpiWrapper::Reduction::Min ),
                         MPI_COMM_GEOSX );

  GEOSX_LOG_LEVEL_RANK_0( 1, getName() << ": Min displacement (X, Y, Z): "
                                       << nodeStatistics.minDisplacement[0] << ", "
                                       << nodeStatistics.minDisplacement[1] << ", "
                                       << nodeStatistics.minDisplacement[2] << " m" );
  GEOSX_LOG_LEVEL_RANK_0( 1, getName() << ": Max displacement (X, Y, Z): "
                                       << nodeStatistics.maxDisplacement[0] << ", "
                                       << nodeStatistics.maxDisplacement[1] << ", "
                                       << nodeStatistics.maxDisplacement[2] << " m" );
}

REGISTER_CATALOG_ENTRY( TaskBase,
                        SolidMechanicsStatistics,
                        string const &, dataRepository::Group * const )

} /* namespace geosx */
