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
 * @file GravityFE.cpp
 */

#include "GravityFE.hpp"
#include "GravityFEKernel.hpp"

#include "finiteElement/FiniteElementDiscretization.hpp"
#include "mesh/DomainPartition.hpp"
#include "common/MpiWrapper.hpp"


namespace geos
{

using namespace constitutive;
using namespace dataRepository;

constexpr localIndex MAX_SUPPORT_POINTS = 8;



GravityFE::GravityFE( const std::string & name,
                      Group * const parent ):
  GravitySolverBase( name, parent )
{}

GravityFE::~GravityFE()
{}


void GravityFE::initializePreSubGroups()
{
  GravitySolverBase::initializePreSubGroups();

  DomainPartition & domain = getGroupByPath< DomainPartition >( "/Problem/domain" );
  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteElementDiscretizationManager const &
  feDiscretizationManager = numericalMethodManager.getFiniteElementDiscretizationManager();

  FiniteElementDiscretization const * const
  feDiscretization = feDiscretizationManager.getGroupPointer< FiniteElementDiscretization >( m_discretizationName );
  GEOS_THROW_IF( feDiscretization == nullptr,
                 getName() << ": FE discretization not found: " << m_discretizationName,
                 InputError );
  GEOS_THROW_IF( feDiscretization->getOrder() >1,
                 getName() << ": FE discretization order should be 1, but we get: " << feDiscretization->getOrder(),
                 InputError );
}


void GravityFE::registerDataOnMesh( Group & meshBodies )
{
  forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                    MeshLevel & mesh,
                                                    string_array const & )
  {
    NodeManager & nodeManager = mesh.getNodeManager();
    nodeManager.registerField< fields::VolumeIntegral >( this->getName() );

    ElementRegionManager & elemManager = mesh.getElemManager();
    elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
    {
      subRegion.registerField< fields::MediumDensity >( this->getName() );
      subRegion.registerField< fields::Adjoint >( this->getName() );
      subRegion.registerField< fields::VolumeIntegral2d >( this->getName() );

      // Assume the maximum number of points per cell is MAX_SUPPORT_POINTS...
      subRegion.getField< fields::VolumeIntegral2d >().resizeDimension< 1 >( MAX_SUPPORT_POINTS );
    } );
  } );
}


void GravityFE::postInputInitialization()
{
  GravitySolverBase::postInputInitialization();
}


void GravityFE::initializePostInitialConditionsPreSubGroups()
{
  GravitySolverBase::initializePostInitialConditionsPreSubGroups();
}


real64 GravityFE::explicitStepModeling( real64 const & time_n,
                                        real64 const & dt,
                                        integer const cycleNumber,
                                        DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( time_n, cycleNumber );

  array1d< real64 > localGzAtStations( m_stationCoordinates.size( 0 ) );
  localGzAtStations.setValues< parallelHostPolicy >( 0. );
  auto localGzAtStationsView = localGzAtStations.toView();

  forDiscretizationOnMeshTargets( domain.getMeshBodies(),
                                  [&] ( string const &,
                                        MeshLevel & mesh,
                                        string_array const & regionNames )
  {
    // Step 1: Precompute volumeIntegral once and for all.
    NodeManager & nodeManager = mesh.getNodeManager();
    auto const nodeGhostRank = nodeManager.ghostRank().toViewConst();

    // VolumeIntegral matrix to be computed in this function.
    auto volumeIntegral = nodeManager.getField< fields::VolumeIntegral >().toView();


    volumeIntegral.setValues< parallelHostPolicy >( 0. );

    // Loop over all sub-regions in regions of type "CellElements".
    mesh.getElemManager().forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                          CellElementSubRegion & elementSubRegion )
    {
      arrayView1d< real64 > const density = elementSubRegion.getReference< array1d< real64 > >( fields::MediumDensity::key());
      if( this->getLogLevel()>3 )
      {
        for( localIndex i=0; i<elementSubRegion.size(); ++i )
        {
          GEOS_LOG( "GravityFE: Cell["<<i<<" / "<< elementSubRegion.size()<<", "<< elementSubRegion.getElementType()<<"]: density= " << density[i] );
        }
      }

      // Compute volume integral at each mesh node.
      arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes = elementSubRegion.nodeList();
      arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition().toViewConst();

      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

      finiteElement::dispatchlowOrder3D( fe, [&] ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        gravityFEKernel::
          ForwardVolumeIntegralKernel< FE_TYPE > kernel( finiteElement );
        kernel.template launch< EXEC_POLICY, ATOMIC_POLICY >
          ( elementSubRegion.size(),
          X,
          elemsToNodes,
          density,
          volumeIntegral );
      } );

      // Zero out ghost node.
      forAll< parallelDevicePolicy<> >( nodeManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
      {
        if( nodeGhostRank[a]>=0 )
        {
          volumeIntegral[a]=0.;
        }
      } );
    } );  // loop cellElements


    // Step #2: Compute contribution to all stations.
    // Make volumeIntegral const
    arrayView1d< real64 const > const volumeIntegralConst = volumeIntegral.toViewConst();
  
    // Hardcode lowest order... i.e. take advantage that mesh vertices and quadrature nodes are collocated.
    // Warning: this code will fail for higher order.
    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const nodePosition =
      nodeManager.referencePosition().toViewConst();
    auto const station = m_stationCoordinates.toViewConst();

    // Constants
    constexpr real64 GRAVITATIONAL_CONSTANT_LOCAL = GRAVITATIONAL_CONSTANT;
    constexpr real64 EPSILON = 1e-15;


    for( localIndex iStation=0; iStation<m_stationCoordinates.size( 0 ); ++iStation )
    {
      // Deal with one station.
      const real64 cx = station( iStation, 0 );
      const real64 cy = station( iStation, 1 );
      const real64 cz = station( iStation, 2 );

      RAJA::ReduceSum< parallelDeviceReduce, real64 > gz( 0.0 );

      forAll< parallelDevicePolicy<> >( nodeManager.size(), [=] GEOS_HOST_DEVICE ( localIndex const a )
      {
        const real64 dx = nodePosition( a, 0 ) - cx;
        const real64 dy = nodePosition( a, 1 ) - cy;
        const real64 dz = nodePosition( a, 2 ) - cz;
        const real64 r2 = dx*dx + dy*dy + dz*dz;
        real64 inv_r3 = 0.0;
        if( r2 > EPSILON )
        {
          inv_r3 = 1.0 / ( r2 * std::sqrt( r2 ) );
        }
        gz += GRAVITATIONAL_CONSTANT_LOCAL * volumeIntegralConst[a] * dz * inv_r3;
      } );

      localGzAtStationsView[iStation]=gz.get();
    } // Loop station
  } );   // Loop mesh


  // Brute force reduce...
  MpiWrapper::allReduce( localGzAtStations, localGzAtStations, MpiWrapper::Reduction::Sum, MPI_COMM_GEOS );

  arrayView1d< real64 > gzAtStations = m_gzAtStations.toView();
  for( localIndex i = 0; i < gzAtStations.size(); ++i )
  {
    gzAtStations[i] = localGzAtStations[i];
  }

  if( this->getLogLevel()>1 )
  {
    for( localIndex iStation = 0; iStation < m_stationCoordinates.size( 0 ); ++iStation )
    {
      auto const & coords = m_stationCoordinates[iStation];
      std::ostringstream logStream;
      logStream << std::fixed << std::setprecision( 2 );
      logStream << "GravityFE: station[" << std::setw( 5 ) << iStation << "] "
                << std::setw( 15 ) << coords[0] << " "
                << std::setw( 15 ) << coords[1] << " "
                << std::setw( 10 )  << coords[2] << " "
                << std::scientific << std::setprecision( 6 )
                << std::setw( 14 ) << gzAtStations[iStation];
      GEOS_LOG_RANK_0( logStream.str() );
    }
  }

  return dt;
}


real64 GravityFE::explicitStepAdjoint( real64 const & time_n,
                                       real64 const & dt,
                                       integer const cycleNumber,
                                       DomainPartition & domain )
{
  GEOS_MARK_FUNCTION;
  GEOS_UNUSED_VAR( time_n, cycleNumber );

  array1d< real64 > localGzAtStations( m_stationCoordinates.size( 0 ) );
  localGzAtStations.setValues< parallelHostPolicy >( 0. );

  forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&] ( string const &,
                                                                MeshLevel & mesh,
                                                                string_array const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();
    NodeManager const & nodeManager = mesh.getNodeManager();

    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const nodePosition =
      nodeManager.referencePosition().toViewConst();

    elemManager.forElementSubRegions< CellElementSubRegion >( regionNames, [&]( localIndex const,
                                                                                CellElementSubRegion & elementSubRegion )
    {
      // Const views for input data
      arrayView1d< real64 const > const residue = m_residue.toViewConst();
      auto const station = m_stationCoordinates.toViewConst();
      arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes = elementSubRegion.nodeList();
      arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition().toViewConst();

      // Output arrays
      arrayView2d< real64 > volumeIntegral2d = elementSubRegion.getReference< array2d< real64 > >( fields::VolumeIntegral2d::key()).toView();
      arrayView1d< real64 > adjoint = elementSubRegion.getReference< array1d< real64 > >( fields::Adjoint::key()).toView();

      volumeIntegral2d.setValues< parallelHostPolicy >( 0.0 );
      adjoint.zero();

      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );

      localIndex const numSupportPoints = fe.getNumSupportPoints();
      localIndex const numStations = m_stationCoordinates.size( 0 );
      GEOS_ERROR_IF_GT_MSG( numSupportPoints, MAX_SUPPORT_POINTS, GEOS_FMT( "Maximum number of SupportPoints is {}", MAX_SUPPORT_POINTS ) );

      finiteElement::dispatchlowOrder3D( fe, [&] ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        gravityFEKernel::
          AdjointVolumeIntegralKernel< FE_TYPE > kernel( finiteElement );
        kernel.template launch< EXEC_POLICY, ATOMIC_POLICY >
          ( elementSubRegion.size(),
          X,
          elemsToNodes,
          volumeIntegral2d );
      } );

      // Make volumeIntegral2d const
      arrayView2d< real64 const > const volumeIntegral2dConst = volumeIntegral2d.toViewConst();

      // Constants
      constexpr real64 GRAVITATIONAL_CONSTANT_LOCAL = GRAVITATIONAL_CONSTANT;
      constexpr real64 EPSILON = 1e-15;

      // Optimized memory access: compute all stations for each element in one kernel
      forAll< EXEC_POLICY >( elementSubRegion.size(), [=] GEOS_HOST_DEVICE ( localIndex const k )
      {
        adjoint[k] = 0.0; // Initialize once
        
        for( localIndex iStation = 0; iStation < numStations; ++iStation )
        {
          const real64 cx = station( iStation, 0 );
          const real64 cy = station( iStation, 1 );
          const real64 cz = station( iStation, 2 );
          const real64 res = residue[iStation];
          
          for( localIndex iLoc = 0; iLoc < numSupportPoints; ++iLoc )
          {
            const localIndex a = elemsToNodes( k, iLoc );
            const real64 dx = nodePosition( a, 0 ) - cx;
            const real64 dy = nodePosition( a, 1 ) - cy;
            const real64 dz = nodePosition( a, 2 ) - cz;
            const real64 r2 = dx*dx + dy*dy + dz*dz;
            real64 inv_r3 = 0.0;
            if( r2 > EPSILON )
            {
              inv_r3 = 1.0 / ( r2 * std::sqrt( r2 ) );
            }
            
            adjoint[k] += GRAVITATIONAL_CONSTANT_LOCAL * volumeIntegral2dConst( k, iLoc ) * res * dz * inv_r3;
          }
        }
      } );  // Loop elem

      adjoint.move( LvArray::MemorySpace::host, true );

      if( this->getLogLevel() > 2 )
      {
        for( localIndex i=0; i<elementSubRegion.size(); ++i )
        {
          GEOS_LOG( "GravityFE: adjoint [" << i << "]= "<< adjoint[i] );
        }
      }

    } ); // Loop subregion
  } );   // Loop mesh

  return dt;
}

REGISTER_CATALOG_ENTRY( PhysicsSolverBase, GravityFE, string const &, dataRepository::Group * const )

} // namespace geos
