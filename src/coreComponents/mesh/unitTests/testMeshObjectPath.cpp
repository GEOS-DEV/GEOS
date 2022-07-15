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
 * @file testMeshObjectPath.cpp
 */

#define MESH_OBJECT_PATH_PRIVATE_FUNCTION_UNIT_TESTING
#include "../MeshObjectPath.hpp"
#include "../MeshBody.hpp"
#include "dataRepository/Group.hpp"

#include <gtest/gtest.h>

namespace geosx
{
using namespace dataRepository;

class TestMesh
{
public:


  static TestMesh & getTestMesh()
  {
    static TestMesh testMesh;
    return testMesh;
  }


  MeshObjectPath::permutationMapType const & pathPermutations() const
  {
    return m_pathPermutations;
  }

  Group const & meshBodies() const
  {
    return m_meshBodies;
  }

  Group & meshBodies()
  {
    return m_meshBodies;
  }

private:
  TestMesh():
    m_pathPermutations(),
    m_rootNode(),
    m_meshBodies( "meshBodies", m_rootNode )
  {
    createTestMesh();
  }

  void createTestMesh();

  MeshObjectPath::permutationMapType m_pathPermutations;
  conduit::Node m_rootNode;
  Group m_meshBodies;
};

void TestMesh::createTestMesh()
{

  m_pathPermutations["body0"]["level0"]["region1"]={"subreg0", "subreg1"};
  m_pathPermutations["body0"]["level0"]["region0"]={"subreg0", "subreg1"};
  m_pathPermutations["body0"]["level1"]["region0"]={"subreg0", "subreg1"};
  m_pathPermutations["body0"]["level1"]["region1"]={"subreg0", "subreg1"};
  m_pathPermutations["body1"]["level0"]["region0"]={"subreg0", "subreg1"};
  m_pathPermutations["body1"]["level0"]["region1"]={"subreg0", "subreg1"};
  m_pathPermutations["body1"]["level1"]["region0"]={"subreg0", "subreg1"};
  m_pathPermutations["body1"]["level1"]["region1"]={"subreg0", "subreg1"};
  m_pathPermutations["body1"]["level1"]["region2"]={"subreg0", "subreg2"};
  m_pathPermutations["body3"]["level0"]["region0"]={"subreg0", "subreg1"};
  m_pathPermutations["body3"]["level0"]["region1"]={"subreg0", "subreg1"};
  m_pathPermutations["body3"]["level2"]["region0"]={"subreg0", "subreg1"};
  m_pathPermutations["body3"]["level2"]["region1"]={"subreg0", "subreg1"};


  for( auto const & meshBodyPair : m_pathPermutations )
  {
    MeshBody & meshBody = m_meshBodies.registerGroup< MeshBody >( meshBodyPair.first );
    for( auto const & meshLevelPair : meshBodyPair.second )
    {
      MeshLevel & level = meshBody.createMeshLevel( meshLevelPair.first );
      ElementRegionManager & elemRegMan = level.getElemManager();
      for( auto const & elemRegionPair : meshLevelPair.second )
      {
        elemRegMan.createChild( "CellElementRegion", elemRegionPair.first );
        CellElementRegion & elemRegion = elemRegMan.getRegion< CellElementRegion >( elemRegionPair.first );
        for( auto const & elemSubRegionName : elemRegionPair.second )
        {
          elemRegion.createElementSubRegion< CellElementSubRegion >( elemSubRegionName );
        }
      }
    }
  }
}

TEST( testMeshObjectPath, checkObjectTypeConsistency )
{
  TestMesh & testMesh = TestMesh::getTestMesh();
  Group const & meshBodies = testMesh.meshBodies();

  {
    string const path = "nodeManager";
    MeshObjectPath meshObjectPath( path, meshBodies );
    meshObjectPath.testCheckObjectTypeConsistency< NodeManager >();
    EXPECT_DEATH_IF_SUPPORTED( meshObjectPath.testCheckObjectTypeConsistency< EdgeManager >(); , ".*" );
    EXPECT_DEATH_IF_SUPPORTED( meshObjectPath.testCheckObjectTypeConsistency< FaceManager >(); , ".*" );
    EXPECT_DEATH_IF_SUPPORTED( meshObjectPath.testCheckObjectTypeConsistency< ElementRegionManager >(); , ".*" );
  }
}

TEST( testMeshObjectPath, fillPathTokens )
{
  TestMesh & testMesh = TestMesh::getTestMesh();
  Group const & meshBodies = testMesh.meshBodies();

  {
    string const path = "ElementRegions";
    MeshObjectPath meshObjectPath( path, meshBodies );
    auto pathTokens = meshObjectPath.testFillPathTokens( path, meshBodies );

    EXPECT_TRUE( pathTokens.size() == 5 );
    EXPECT_TRUE( pathTokens[0] == "*" );
    EXPECT_TRUE( pathTokens[1] == "*" );
    EXPECT_TRUE( pathTokens[2] == "ElementRegions" );
    EXPECT_TRUE( pathTokens[3] == "*" );
    EXPECT_TRUE( pathTokens[4] == "*" );
  }

  {
    string const path = "nodeManager";
    MeshObjectPath meshObjectPath( path, meshBodies );
    auto pathTokens = meshObjectPath.testFillPathTokens( path, meshBodies );

    EXPECT_TRUE( pathTokens.size() == 3 );
    EXPECT_TRUE( pathTokens[0] == "*" );
    EXPECT_TRUE( pathTokens[1] == "*" );
    EXPECT_TRUE( pathTokens[2] == "nodeManager" );
  }

  {
    string const path = "edgeManager";
    MeshObjectPath meshObjectPath( path, meshBodies );
    auto pathTokens = meshObjectPath.testFillPathTokens( path, meshBodies );

    EXPECT_TRUE( pathTokens.size() == 3 );
    EXPECT_TRUE( pathTokens[0] == "*" );
    EXPECT_TRUE( pathTokens[1] == "*" );
    EXPECT_TRUE( pathTokens[2] == "edgeManager" );
  }

  {
    string const path = "faceManager";
    MeshObjectPath meshObjectPath( path, meshBodies );
    auto pathTokens = meshObjectPath.testFillPathTokens( path, meshBodies );

    EXPECT_TRUE( pathTokens.size() == 3 );
    EXPECT_TRUE( pathTokens[0] == "*" );
    EXPECT_TRUE( pathTokens[1] == "*" );
    EXPECT_TRUE( pathTokens[2] == "faceManager" );
  }

}

TEST( testMeshObjectPath, ExpandPathTokens )
{}


TEST( testMeshObjectPath, fullPathExpansion )
{
  TestMesh & testMesh = TestMesh::getTestMesh();
  Group const & meshBodies = testMesh.meshBodies();
  MeshObjectPath::permutationMapType const & pathPermutations = testMesh.pathPermutations();

  {
    string const path = "*/*/ElementRegions";
    MeshObjectPath meshObjectPath( path, meshBodies );
    EXPECT_TRUE( meshObjectPath.pathPermutations() == pathPermutations );
  }

  MeshObjectPath::permutationMapType pathPermutationSub = pathPermutations;

  for( auto & meshBodyPair : pathPermutationSub )
  {
    for( auto & meshLevelPair : meshBodyPair.second )
    {
      meshLevelPair.second.clear();
    }
  }

  {
    string const path = "nodeManager";
    MeshObjectPath meshObjectPath( path, meshBodies );
    EXPECT_TRUE( meshObjectPath.pathPermutations() == pathPermutationSub );
  }
  {
    string const path = "edgeManager";
    MeshObjectPath meshObjectPath( path, meshBodies );
    EXPECT_TRUE( meshObjectPath.pathPermutations() == pathPermutationSub );
  }
  {
    string const path = "faceManager";
    MeshObjectPath meshObjectPath( path, meshBodies );
    EXPECT_TRUE( meshObjectPath.pathPermutations() == pathPermutationSub );
  }
}


TEST( testMeshObjectPath, invalidMeshBody )
{
  TestMesh & testMesh = TestMesh::getTestMesh();
  Group const & meshBodies = testMesh.meshBodies();
  {
    string const path = "*/level2/ElementRegions";
    EXPECT_DEATH_IF_SUPPORTED( MeshObjectPath meshObjectPath( path, meshBodies ), ".*" );
  }
}


TEST( testMeshObjectPath, invalidMeshLevel )
{
  TestMesh & testMesh = TestMesh::getTestMesh();
  Group const & meshBodies = testMesh.meshBodies();
  {
    string const path = "*/*/ElementRegions/{region2}";
    EXPECT_DEATH_IF_SUPPORTED( MeshObjectPath meshObjectPath( path, meshBodies ), ".*" );
  }
}

TEST( testMeshObjectPath, invalidMeshRegion )
{
  TestMesh & testMesh = TestMesh::getTestMesh();
  Group const & meshBodies = testMesh.meshBodies();
  {
    string const path = "*/*/ElementRegions/*/subreg2";
    EXPECT_DEATH_IF_SUPPORTED( MeshObjectPath meshObjectPath( path, meshBodies ), ".*" );
  }
}


void checkSubRegionNames( std::vector< string > const & names )
{
  EXPECT_TRUE( names[0] == "subreg0" );
  EXPECT_TRUE( names[1] == "subreg1" );
  EXPECT_TRUE( names[2] == "subreg0" );
  EXPECT_TRUE( names[3] == "subreg1" );
  EXPECT_TRUE( names[4] == "subreg0" );
  EXPECT_TRUE( names[5] == "subreg1" );
  EXPECT_TRUE( names[6] == "subreg0" );
  EXPECT_TRUE( names[7] == "subreg1" );
  EXPECT_TRUE( names[8] == "subreg0" );
  EXPECT_TRUE( names[9] == "subreg1" );
  EXPECT_TRUE( names[10] == "subreg0" );
  EXPECT_TRUE( names[11] == "subreg1" );
  EXPECT_TRUE( names[12] == "subreg0" );
  EXPECT_TRUE( names[13] == "subreg1" );
  EXPECT_TRUE( names[14] == "subreg0" );
  EXPECT_TRUE( names[15] == "subreg1" );
  EXPECT_TRUE( names[16] == "subreg0" );
  EXPECT_TRUE( names[17] == "subreg2" );
  EXPECT_TRUE( names[18] == "subreg0" );
  EXPECT_TRUE( names[19] == "subreg1" );
  EXPECT_TRUE( names[20] == "subreg0" );
  EXPECT_TRUE( names[21] == "subreg1" );
  EXPECT_TRUE( names[22] == "subreg0" );
  EXPECT_TRUE( names[23] == "subreg1" );
  EXPECT_TRUE( names[24] == "subreg0" );
  EXPECT_TRUE( names[25] == "subreg1" );
}

void checkRegionNames( std::vector< string > const & names )
{
  EXPECT_TRUE( names[0] == "region0" );
  EXPECT_TRUE( names[1] == "region1" );
  EXPECT_TRUE( names[2] == "region0" );
  EXPECT_TRUE( names[3] == "region1" );
  EXPECT_TRUE( names[4] == "region0" );
  EXPECT_TRUE( names[5] == "region1" );
  EXPECT_TRUE( names[6] == "region0" );
  EXPECT_TRUE( names[7] == "region1" );
  EXPECT_TRUE( names[8] == "region2" );
  EXPECT_TRUE( names[9] == "region0" );
  EXPECT_TRUE( names[10] == "region1" );
  EXPECT_TRUE( names[11] == "region0" );
  EXPECT_TRUE( names[12] == "region1" );
}

TEST( testMeshObjectPath, forObjectsInPath )
{
  TestMesh & testMesh = TestMesh::getTestMesh();
  Group & meshBodies = testMesh.meshBodies();
  Group const & meshBodiesConst = meshBodies;
  string const path = "*/*/ElementRegions";
  MeshObjectPath meshObjectPath( path, meshBodiesConst );

  {
    std::vector< string > names;
    meshObjectPath.forObjectsInPath< CellElementSubRegion >( meshBodiesConst,
                                                             [&]( ElementSubRegionBase const & elemSubRegionBase )
    {
      names.push_back( elemSubRegionBase.getName() );
    } );
    checkSubRegionNames( names );
  }

  {
    std::vector< string > names;
    meshObjectPath.forObjectsInPath< CellElementSubRegion >( meshBodies,
                                                             [&]( ElementSubRegionBase & elemSubRegionBase )
    {
      names.push_back( elemSubRegionBase.getName() );
    } );
    checkSubRegionNames( names );
  }

  {
    std::vector< string > names;
    meshObjectPath.forObjectsInPath< CellElementRegion >( meshBodiesConst,
                                                          [&]( CellElementRegion const & elemRegionBase )
    {
      names.push_back( elemRegionBase.getName() );
    } );
    checkRegionNames( names );
  }
}


} /* namespace geosx */
