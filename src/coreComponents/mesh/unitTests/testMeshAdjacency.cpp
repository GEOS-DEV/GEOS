/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/** @file testMeshAdjacency.cpp */
#include "common/MpiWrapper.hpp"
#include "mesh/MeshLevel.hpp"
#include "mesh/FaceElementSubRegion.hpp"
#include "mesh/SurfaceElementRegion.hpp"

#include <gtest/gtest.h>
#include <conduit.hpp>
#include <algorithm>
#include <cstdint>
#include <initializer_list>
#include <utility>

using namespace geos;

namespace
{
class CollocatedAdjacency : public ::testing::Test
{
protected:
  conduit::Node repository;
  dataRepository::Group root{ "Problem", repository };
  MeshLevel mesh{ "Level0", &root };
  globalIndex const base = sizeof( globalIndex ) == 8 ? INT64_C( 9007199254741001 ) : 10001;

  CollocatedAdjacency()
  {
    auto & nodes = mesh.getNodeManager();
    nodes.resize( 4 );
    for( localIndex n = 0; n < nodes.size(); ++n ) nodes.localToGlobalMap()[n] = base + n;
    static_cast< ObjectManagerBase & >( nodes ).constructGlobalToLocalMap();
  }

  void addFracture( string const & name, std::initializer_list< stdVector< globalIndex > > values )
  {
    auto & elements = mesh.getElemManager();
    elements.createChild( "SurfaceElementRegion", name );
    auto & region = elements.getRegion< SurfaceElementRegion >( name );
    auto & faces = region.createElementSubRegion< FaceElementSubRegion >( "faces" );
    faces.resize( 1 );
    faces.registerWrapper< array1d< localIndex > >( ObjectManagerBase::viewKeyStruct::adjacencyListString() );
    auto & buckets = faces.getReference< ArrayOfArrays< array1d< globalIndex > > >(
      FaceElementSubRegion::viewKeyStruct::elem2dToCollocatedNodesBucketsString() );
    buckets.resize( 1 );
    for( auto const & valuesInBucket : values )
    {
      array1d< globalIndex > bucket;
      bucket.resize( valuesInBucket.size() );
      std::copy( valuesInBucket.begin(), valuesInBucket.end(), bucket.begin() );
      buckets.emplaceBack( 0, std::move( bucket ) );
    }
  }

  void check( stdVector< localIndex > const & seeds, integer depth, stdVector< localIndex > const & expected )
  {
    array1d< localIndex > nodes, edges, faces;
    auto elements = mesh.getElemManager().constructReferenceAccessor< array1d< localIndex > >(
      ObjectManagerBase::viewKeyStruct::adjacencyListString() );
    array1d< localIndex > seedNodes;
    seedNodes.resize( seeds.size() );
    std::copy( seeds.begin(), seeds.end(), seedNodes.begin() );
    mesh.generateAdjacencyLists( seedNodes.toViewConst(), nodes, edges, faces, elements, depth );
    EXPECT_EQ( ( stdVector< localIndex >( nodes.begin(), nodes.end() ) ), expected );
    EXPECT_EQ( edges.size(), 0 );
    EXPECT_EQ( faces.size(), 0 );
  }
};
} // namespace

TEST_F( CollocatedAdjacency, OneHopPerSubregionAndDepth )
{
  addFracture( "fault", { { base, base + 1 }, { base + 1, base + 2 } } );
  check( { 0 }, 0, { 0 } );
  check( { 0 }, 1, { 0, 1 } );
  check( { 0 }, 2, { 0, 1, 2 } );
  check( { 2 }, 1, { 1, 2 } );
  check( { 3 }, 2, { 3 } );
  check( {}, 2, {} );
}

TEST_F( CollocatedAdjacency, JunctionsDuplicatesAndRemoteNodes )
{
  addFracture( "fault", { { base, base + 1, base + 2, base + 9 },
                 { base, base + 1, base + 2, base + 9 }, { base + 8, base + 9 }, {} } );
  check( { 0 }, 1, { 0, 1, 2 } );
  check( { 1 }, 1, { 0, 1, 2 } );
  check( { 2, 3 }, 1, { 0, 1, 2, 3 } );
}

TEST_F( CollocatedAdjacency, NamedSubregionsPreserveSequentialExpansion )
{
  addFracture( "first", { { base, base + 1 } } );
  addFracture( "second", { { base + 1, base + 2 } } );
  check( { 0 }, 1, { 0, 1, 2 } );
}

int main( int argc, char * * argv )
{
  MpiWrapper::init( &argc, &argv );
  MPI_COMM_GEOS = MpiWrapper::commDup( MPI_COMM_WORLD );
  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();
  MpiWrapper::finalize();
  return result;
}
