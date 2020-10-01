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


// Source includes
#include "managers/initialization.hpp"
#include "dataRepository/xmlWrapper.hpp"
#include "tests/graphFileNames.hpp"
#include "mesh/GraphManager.hpp"
#include "mesh/GraphBase.hpp"
#include "mesh/GraphVertex.hpp"
#include "mesh/GraphFromText.hpp"




// TPL includes
#include <gtest/gtest.h>


using namespace geosx;
using namespace geosx::dataRepository;

void TestGraphModif( string const & inputStringGraph)
{
  GraphManager graphManager( "graph", nullptr );

  // Load the graph
  xmlWrapper::xmlDocument xmlDocument;
  xmlDocument.load_buffer( inputStringGraph.c_str(), inputStringGraph.size() );

  
  xmlWrapper::xmlNode xmlGraphNode = xmlDocument.child( "Graph" );
  
  graphManager.ProcessInputFileRecursive( xmlGraphNode );
  graphManager.PostProcessInputRecursive();
    
  graphManager.GenerateGraphs( );
  for( auto & graph : graphManager.GetSubGroups() )
  {
    GraphFromText * gb = Group::group_cast< GraphFromText * >( graph.second );
    gb->AddVertex(0,0,0);
    gb->AddVertex(0,0,1);
    gb->AddVertex(0,0,2);
    gb->AddVertex(0,0,3);
    gb->AddVertex(0,0,4);
    gb->AddVertex(0,0,5);
    array1d<std::shared_ptr<GraphVertex>> vertices = gb->getVertices();
    gb->AddEdge(0,vertices[0],vertices[1],0);
    gb->AddEdge(1,vertices[1],vertices[2],0);
    gb->AddEdge(2,vertices[2],vertices[3],0);
    gb->AddEdge(3,vertices[3],vertices[4],0);
    gb->AddEdge(4,vertices[4],vertices[5],0);
    gb->AddEdge(5,vertices[5],vertices[0],0);
    gb->AddEdge(6,vertices[0],vertices[2],0);
    gb->AddEdge(7,vertices[2],vertices[5],0);
    gb->AddEdge(8,vertices[5],vertices[3],0);
    gb->AddEdge(9,vertices[3],vertices[1],0);
    gb->AddEdge(10,vertices[1],vertices[4],0);
    gb->AddEdge(11,vertices[4],vertices[0],0);

    std::map<std::shared_ptr<GraphVertex>, array1d<GraphEdge*>> map = gb->getVertexWithEdgesMap();
    for ( auto & vertex : map )
    {
      array1d<GraphEdge*> edges = vertex.second;
      GEOSX_ASSERT_EQ(edges.size(),4);
    }
    gb->RemoveVertex(vertices[1]);
    gb->RemoveVertex(vertices[3]);
    gb->RemoveVertex(vertices[5]);
    gb->AddEdge(12,vertices[2],vertices[4],0);
    map = gb->getVertexWithEdgesMap();

    for ( auto & vertex : map )
    {
      array1d<GraphEdge*> edges = vertex.second;
      GEOSX_ASSERT_EQ(edges.size(), 2);
    }

  }

  
}

void TestGraphGenerate( string const & inputStringGraph)
{
  GraphManager graphManager( "graph", nullptr );

  // Load the graph
  xmlWrapper::xmlDocument xmlDocument;
  xmlDocument.load_buffer( inputStringGraph.c_str(), inputStringGraph.size() );

  
  xmlWrapper::xmlNode xmlGraphNode = xmlDocument.child( "Graph" );
  
  graphManager.ProcessInputFileRecursive( xmlGraphNode );
  graphManager.PostProcessInputRecursive();
    
  graphManager.GenerateGraphs( );
  for( auto & graph : graphManager.GetSubGroups() )
  {
    GraphFromText * gb = Group::group_cast< GraphFromText * >( graph.second );
    std::map<std::shared_ptr<GraphVertex>, array1d<GraphEdge*>> map = gb->getVertexWithEdgesMap();
    std::cout<<gb->getFile();
    for ( auto & vertex : map )
    {
      array1d<GraphEdge*> edges = vertex.second;
      GEOSX_ASSERT_EQ(edges.size(),2);
    }
  }
}

TEST( Graph, GraphModification )
{
  
  GraphManager graphManager( "graph", nullptr );

  std::stringstream inputStreamGraph;
  inputStreamGraph <<
    "<?xml version=\"1.0\" ?>" <<
    "  <Graph xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:noNamespaceSchemaLocation=\"geos_v0.0.xsd\">" <<
    "  <GraphFromText name=\"BasicGraph\" " <<
    "  file=\"" <<"Foo.txt"<< "\"/>"<<
    "</Graph>";
  const string inputStringGraph = inputStreamGraph.str();

  TestGraphModif( inputStringGraph);
}

TEST( Graph, GraphGeneration )
{
  
  GraphManager graphManager( "graph", nullptr );

  std::stringstream inputStreamGraph;
  inputStreamGraph <<
    "<?xml version=\"1.0\" ?>" <<
    "  <Graph xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:noNamespaceSchemaLocation=\"geos_v0.0.xsd\">" <<
    "  <GraphFromText name=\"BasicGraph\" " <<
    "  file=\"" <<txtFilePath.c_str()<< "\"/>"<<
    "</Graph>";
  const string inputStringGraph = inputStreamGraph.str();

  TestGraphGenerate( inputStringGraph);
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geosx::basicSetup( argc, argv );

  int const result = RUN_ALL_TESTS();

  geosx::basicCleanup();

  return result;
}
