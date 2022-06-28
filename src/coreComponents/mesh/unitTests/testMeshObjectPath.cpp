/*
 * testMeshPath.cpp
 *
 *  Created on: Jun 21, 2022
 *      Author: settgast
 */

#include "../MeshObjectPath.hpp"
#include "../MeshBody.hpp"
#include "dataRepository/Group.hpp"

#include <gtest/gtest.h>

namespace geosx
{
using namespace dataRepository;


//void shortPathCheck( string const & path )
//{
//  conduit::Node rootNode;
//  Group meshBodies( "meshBodies", rootNode );
//  MeshBody & meshBody1 = meshBodies.registerGroup<MeshBody>( "body1" );
//  MeshBody & meshBody2 = meshBodies.registerGroup<MeshBody>( "body2" );
//  meshBody1.createMeshLevel("level1");
//  meshBody1.createMeshLevel("level2");
//  meshBody2.createMeshLevel("level3");
//  meshBody2.createMeshLevel("level4");
//
//  {
//    MeshObjectPath meshObjectPath( path, meshBodies );
//    auto const & meshBodyLevelMap = meshObjectPath.meshBodyLevelMap();
////  auto const & regionSubRegionMap = meshObjectPath.regionSubRegionMap();
//
//    EXPECT_STREQ( meshObjectPath.getObjectType().c_str(), path.c_str() );
//
//    EXPECT_EQ( meshBodyLevelMap.size(), 2 );
//
//    auto iter = meshBodyLevelMap.begin();
//    EXPECT_STREQ( iter->first.c_str(), "body1" );
//    EXPECT_EQ( iter->second.size(), 2 );
//    EXPECT_STREQ( iter->second[0].c_str(), "level1" );
//    EXPECT_STREQ( iter->second[1].c_str(), "level2" );
//
//    ++iter;
//    EXPECT_STREQ( iter->first.c_str(), "body2" );
//    EXPECT_EQ( iter->second.size(), 2 );
//    EXPECT_STREQ( iter->second[0].c_str(), "level3" );
//    EXPECT_STREQ( iter->second[1].c_str(), "level4" );
//  }
//}
//
//TEST( testMeshObjectPath, shortPathInputNode )
//{
//  shortPathCheck( MeshLevel::groupStructKeys::nodeManagerString() );
//}
//TEST( testMeshObjectPath, shortPathInputEdge )
//{
//  shortPathCheck( MeshLevel::groupStructKeys::edgeManagerString() );
//}
//TEST( testMeshObjectPath, shortPathInputFace )
//{
//  shortPathCheck( MeshLevel::groupStructKeys::faceManagerString() );
//}



TEST( testMeshObjectPath, fullPathInputNode )
{
  conduit::Node rootNode;
  Group meshBodies( "meshBodies", rootNode );
  MeshBody & meshBody0 = meshBodies.registerGroup<MeshBody>( "body0" );
  MeshBody & meshBody1 = meshBodies.registerGroup<MeshBody>( "body1" );
  MeshLevel & level00 = meshBody0.createMeshLevel("level0");
  MeshLevel & level01 = meshBody0.createMeshLevel("level1");
  MeshLevel & level10 = meshBody1.createMeshLevel("level0");
  MeshLevel & level11 = meshBody1.createMeshLevel("level1");

  ElementRegionManager & elemRegMan00 = level00.getElemManager();

  elemRegMan00.createChild( "CellElementRegion", "Region0" );

  CellElementRegion & elemRegion0 = elemRegMan00.getRegion<CellElementRegion>("Region0");
  elemRegion0.createElementSubRegion<CellElementSubRegion>("subReg0");

  {
    string const path = "ElementRegions";
    MeshObjectPath meshObjectPath( path, meshBodies );
//    std::vector<string> splitPath = meshObjectPath.getPath();
//    EXPECT_STREQ( splitPath[0].c_str(), "body0" );
//    EXPECT_STREQ( splitPath[1].c_str(), "level0" );
//    EXPECT_STREQ( splitPath[2].c_str(), "nodeManager" );
  }



}

//TEST( testMeshObjectPath, fullPathInputElement )
//{
//  conduit::Node rootNode;
//  Group meshBodies( "meshBodies", rootNode );
//  MeshBody & meshBody = meshBodies.registerGroup<MeshBody>( "body" );
//  MeshLevel & meshLevel = meshBody.createMeshLevel("level");
//
//  string const path = "body/level/nodeManager";
//  MeshObjectPath meshObjectPath( path, meshBodies );
//  std::vector<string> splitPath = meshObjectPath.getPath();
//
//  EXPECT_STREQ( splitPath[0].c_str(), "body" );
//  EXPECT_STREQ( splitPath[1].c_str(), "level" );
//  EXPECT_STREQ( splitPath[2].c_str(), "nodeManager" );
//}






} /* namespace geosx */
