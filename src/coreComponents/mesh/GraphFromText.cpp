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
 * @file GraphFromText.cpp
 */
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "common/Path.hpp"
#include "GraphFromText.hpp"
#include "managers/Functions/FunctionManager.hpp"
#include "mesh/GraphEdge.hpp"
#include "mesh/GraphVertex.hpp"
#include "mesh/MeshLevel.hpp"
#include "common/Logger.hpp"

namespace geosx
{

using namespace dataRepository;


GraphFromText::GraphFromText( const std::string & name,
                              Group * const parent ):
  GraphBase( name, parent ),m_file("")
  {
     setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

    // This enables logLevel filtering
    enableLogLevelInput();

    registerWrapper( viewKeyStruct::fileString, &m_file )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of the graph file." );
    
  }


GraphFromText::~GraphFromText()
{
  for (long unsigned int i=0;i<m_edges.size();i++)
  {
    delete(m_edges[i]);
  }
  for(map<GraphVertex*,std::vector<GraphEdge*>>::iterator it = m_vertexWithEdgesMap.begin(); it != m_vertexWithEdgesMap.end(); ++it) 
  {
    delete(it->first);
  }

}


void GraphFromText::AddEdge(localIndex ind, GraphVertex* v1, GraphVertex* v2, real64 transm)
{
  GraphEdge* e;
  e=new GraphEdge(ind, v1, v2, transm);
  m_edges.push_back(e);
  m_vertexWithEdgesMap[v1].push_back(e);
  m_vertexWithEdgesMap[v2].push_back(e);
}

void GraphFromText::RemoveEdge(localIndex ind)
{
  forAll< serialPolicy >( m_edges.size(), [=]( localIndex const i )
  {
    if (m_edges[i]->getEdgeIndex()==ind)
    {
      GraphEdge* edge = m_edges[i];
      //Delete the edge from the edges linked to the Vertex1
      forAll< serialPolicy >( m_vertexWithEdgesMap[edge->getVertex1()].size(), [=]( localIndex const j )
      {
        if (m_vertexWithEdgesMap[edge->getVertex1()][j]->getEdgeIndex()==ind)
        {
          int convertedJ = LvArray::integerConversion< int >(j);
          m_vertexWithEdgesMap[edge->getVertex1()].erase(m_vertexWithEdgesMap[edge->getVertex1()].begin()+convertedJ);
        }
      });
      //Delete the edge from the edges linked to the Vertex2
      forAll< serialPolicy >( m_vertexWithEdgesMap[edge->getVertex2()].size(), [=]( localIndex const j )
      {
        if (m_vertexWithEdgesMap[edge->getVertex2()][j]->getEdgeIndex()==ind)
        {
          int convertedJ = LvArray::integerConversion< int >(j);
          m_vertexWithEdgesMap[edge->getVertex2()].erase(m_vertexWithEdgesMap[edge->getVertex2()].begin()+convertedJ);
        }
      });
      //Delete the edge from the edge list
      int convertedI = LvArray::integerConversion< int >(i);
      m_edges.erase(m_edges.begin()+convertedI);
      delete(edge);
    }
  });

}

void GraphFromText::AddVertex(localIndex er, localIndex esr, globalIndex ei)
{
  GraphVertex* v=new GraphVertex(er,esr,ei);
  m_vertices.push_back(v);
  std::vector<GraphEdge*> edges;
  m_vertexWithEdgesMap.insert({v,edges});
}

void GraphFromText::RemoveVertex(localIndex er, localIndex esr, globalIndex ei)
{
  GraphVertex* v = getVertexWithGlobalIndex(er,esr,ei);
  std::vector<localIndex> indexes;
  localIndex loc = -1;
  for(long unsigned int i = 0; i<m_vertices.size(); i++)
  {
    if (m_vertices[i]->getRegionIndex()==er && m_vertices[i]->getSubRegionIndex()==esr && m_vertices[i]->getGlobalVertexIndex()==ei)
    {
      loc = i;
    }
  }
  if (loc != -1)
  {
    m_vertices.erase(m_vertices.begin()+loc);    
  }
  for(long unsigned int i = 0; i<m_vertexWithEdgesMap[v].size(); i++)
  {
    indexes.push_back(m_vertexWithEdgesMap[v][i]->getEdgeIndex());
  }

  forAll< serialPolicy >( indexes.size(), [=]( localIndex const i )
  {
    //Deleting all edges linked to the vertex to delete (as they will have a void target on one extremity)
    //std::cout<<indexes[i]<<"\n";
    RemoveEdge(indexes[i]);
  });
  m_vertexWithEdgesMap.erase(v);
  delete(v);
  
}

void GraphFromText::RemoveVertex(GraphVertex* vertex)
{
  std::vector<localIndex> indexes;
  localIndex loc = -1;
  for(long unsigned int i = 0; i<m_vertices.size(); i++)
  {
    if (m_vertices[i]->getRegionIndex()==vertex->getRegionIndex() && m_vertices[i]->getSubRegionIndex()==vertex->getSubRegionIndex() && m_vertices[i]->getGlobalVertexIndex()==vertex->getGlobalVertexIndex())
    {
      loc = i;
    }
  }
  if (loc != -1)
  {
    m_vertices.erase(m_vertices.begin()+loc);    
  }

  for(long unsigned int i = 0; i<m_vertexWithEdgesMap[vertex].size(); i++)
  {
    indexes.push_back(m_vertexWithEdgesMap[vertex][i]->getEdgeIndex());
  }

  forAll< serialPolicy >( indexes.size(), [=]( localIndex const i )
  {
    //Deleting all edges linked to the vertex to delete (as they will have a void target on one extremity)
    //std::cout<<indexes[i]<<"\n";
    RemoveEdge(indexes[i]);
  });
  m_vertexWithEdgesMap.erase(vertex);
  delete(vertex);
  
}


GraphVertex* GraphFromText::getVertexWithGlobalIndex(localIndex er, localIndex esr, globalIndex ei)
{
  for(long unsigned int i=0; i < m_vertices.size(); i++)
  {
    if (m_vertices[i]->getGlobalVertexIndex()==ei && m_vertices[i]->getRegionIndex()==er && m_vertices[i]->getSubRegionIndex()==esr  )
    {
      return m_vertices[i];
    }
  }
  return nullptr;
}



void GraphFromText::GenerateGraph()
{
  std::ifstream infile(m_file);
  std::string line;
  //std::vector<GraphEdge*> edges;
  while (std::getline(infile, line))
  {
    std::istringstream iss(line);
    int er1, esr1, ei1, er2, esr2, ei2;
    real64 transm;
    //The input file is written with thr triple index (regionIndex, subRegionIndex,globalVertexIndex) of the two connected point and associated transmissibility
    if (!(iss >>er1 >> esr1 >> ei1 >> er2 >> esr2 >> ei2 >> transm)) { break; } // error
    int place_1=-1;
    int place_2=-1;
    long unsigned int i=0;
    //We search if any of the two point is already created
    //If it is the case we recover its position in the place variables
    while (i<m_vertices.size() && (place_1==-1 || place_2==-1))
    {
      if (m_vertices[i]->getGlobalVertexIndex()==ei1 && m_vertices[i]->getRegionIndex()==er1 && m_vertices[i]->getSubRegionIndex()==esr1 )
      {
        place_1=i;
      }
      if (m_vertices[i]->getGlobalVertexIndex()==ei2 && m_vertices[i]->getRegionIndex()==er2 && m_vertices[i]->getSubRegionIndex()==esr2 )
      {
        place_2=i;
      }
      ++i;
    }
    //Create or recover the two vertices
    //Add them to the vertex list and create an entry in the map for them
    GraphVertex* v1;
    GraphVertex* v2;
    if(place_1==-1)
    {
      v1=new GraphVertex(er1,esr1,ei1);
      m_vertices.push_back(v1);
      std::vector<GraphEdge*> edges;
      m_vertexWithEdgesMap.insert({v1,edges});
    }
    else
    {
      v1=m_vertices[place_1];
    }
    if(place_2==-1)
    {
      v2=new GraphVertex(er2,esr2,ei2);
      m_vertices.push_back(v2);
      std::vector<GraphEdge*> edges;
      m_vertexWithEdgesMap.insert({v2,edges});

    }
    else
    {
      v2=m_vertices[place_2];
    }
    //Create the edge, add it to the edge list
    GraphEdge* e= new GraphEdge(m_edges.size(), v1, v2,transm);
    m_edges.push_back(e);
    //Add the edge to the list of linked edge for his two vertices
    m_vertexWithEdgesMap[v1].push_back(e);
    m_vertexWithEdgesMap[v2].push_back(e);
    }
  /*
  this->RemoveEdge(LvArray::integerConversion< localIndex >(2690));
  this->RemoveEdge(LvArray::integerConversion< localIndex >(2688));
  this->RemoveEdge(LvArray::integerConversion< localIndex >(2671));
  this->RemoveEdge(LvArray::integerConversion< localIndex >(2557));
  this->RemoveEdge(LvArray::integerConversion< localIndex >(2555));
  this->RemoveEdge(LvArray::integerConversion< localIndex >(2538));
  this->RemoveEdge(LvArray::integerConversion< localIndex >(2500));
  this->RemoveEdge(LvArray::integerConversion< localIndex >(2297));
  GraphVertex* v819 = this->getVertexWithIndex(LvArray::integerConversion< localIndex >(819));
  GraphVertex* v919 = this->getVertexWithIndex(LvArray::integerConversion< localIndex >(919));
  GraphVertex* v909 = this->getVertexWithIndex(LvArray::integerConversion< localIndex >(909));
  GraphVertex* v918 = this->getVertexWithIndex(LvArray::integerConversion< localIndex >(918));
  GraphVertex* v929 = this->getVertexWithIndex(LvArray::integerConversion< localIndex >(929));
  GraphVertex* v889 = this->getVertexWithIndex(LvArray::integerConversion< localIndex >(889));
  GraphVertex* v989 = this->getVertexWithIndex(LvArray::integerConversion< localIndex >(989));
  GraphVertex* v979 = this->getVertexWithIndex(LvArray::integerConversion< localIndex >(979));
  GraphVertex* v988 = this->getVertexWithIndex(LvArray::integerConversion< localIndex >(988));
  GraphVertex* v999 = this->getVertexWithIndex(LvArray::integerConversion< localIndex >(999));
  this->AddEdge(LvArray::integerConversion< localIndex >(2690),v999,v919);
  this->AddEdge(LvArray::integerConversion< localIndex >(2688),v988,v919);
  this->AddEdge(LvArray::integerConversion< localIndex >(2671),v979,v919);
  this->AddEdge(LvArray::integerConversion< localIndex >(2557),v929,v989);
  this->AddEdge(LvArray::integerConversion< localIndex >(2555),v918,v989);
  this->AddEdge(LvArray::integerConversion< localIndex >(2538),v909,v989);
  this->AddEdge(LvArray::integerConversion< localIndex >(2500),v889,v919);
  this->AddEdge(LvArray::integerConversion< localIndex >(2297),v819,v989);
  this->RemoveVertex(919);
  this->RemoveVertex(989);
  */
  /*
  for (long unsigned int i=400; i<500; ++i)
  {
    std::cout<<i<<"\n";
    this->RemoveVertex(0,0,i);
  }
  */
}

void GraphFromText::PartitionGraph(const MeshLevel & mesh)
{
  ElementRegionManager const & elemManager = *mesh.getElemManager();
  ElementRegionManager::ElementViewAccessor< arrayView1d< integer const > > const elemGhostRank =
    elemManager.ConstructArrayViewAccessor< integer, 1 >( ObjectManagerBase::viewKeyStruct::ghostRankString );
  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > const localToGlobalMap =
    elemManager.ConstructArrayViewAccessor< globalIndex, 1 >( ObjectManagerBase::viewKeyStruct::localToGlobalMapString );
  //Create and populate a vector containing the vertices to delete as they are not present in the partition
  std::vector <GraphVertex*> vertexToDelete;
  for( long unsigned int i=0; i<m_vertices.size(); i++)
  {
    vertexToDelete.push_back(m_vertices[i]);
  }
  for (long int i = 0; i<elemGhostRank[0][0].size(); i++)
  {
    int vertexToKeep = -1;
    for (long unsigned int j = 0; j<vertexToDelete.size();j++)
    {
            //For each vertex appearing in the elemGhostRank matrix, hence owned or ghosted, we suppress it from the list of point to delete
            if (localToGlobalMap[0][0][i] ==vertexToDelete[j]->getGlobalVertexIndex())
      {
        vertexToKeep=j;
        //The values are copied in the vertex object to avoid requesting the local to global conversion again
        vertexToDelete[j]->setGhostIndex(elemGhostRank[0][0][i]);
        vertexToDelete[j]->setLocalIndex(i);
        //std::cout<<vertexToDelete[j]->getLocalVertexIndex()<<" ";
        //std::cout<<vertexToDelete[j]->getGhostIndex()<<"\n";
      }
    }
    //If we found a corresponding vertex we suppress it from the list to vertex to delete
    if (vertexToKeep != -1)
    {
      vertexToDelete.erase(vertexToDelete.begin()+vertexToKeep);
    }
  }
  //Every point that have not been removed from the list is not present in the partition and is deleted
  for ( long unsigned int i=0; i<vertexToDelete.size(); i++)
  {
    int to_delete = -1;
    for (long unsigned int j = 0; j<m_vertices.size();j++)
    {
      if (m_vertices[j]->getGlobalVertexIndex() == vertexToDelete[i]->getGlobalVertexIndex())
      {
        to_delete=j;          
      }
    }
    if (to_delete != -1)
    {
      RemoveVertex(m_vertices[to_delete]);
    }
  }
  
}


REGISTER_CATALOG_ENTRY( GraphBase, GraphFromText, std::string const &, Group * const )
} /* namespace geosx */
