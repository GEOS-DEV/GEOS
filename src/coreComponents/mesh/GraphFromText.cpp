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
    registerWrapper( viewKeyStruct::meshString, &m_meshString )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Name of the mesh to link to the graph." );

    
  }


GraphFromText::~GraphFromText()
{
  for (long unsigned int i=0;i<m_edges.size();i++)
  {
    delete(m_edges[i]);
  }
  for(map<GraphVertex*,std::vector<GraphEdge*>>::iterator it = m_vertex_edges.begin(); it != m_vertex_edges.end(); ++it) 
  {
    delete(it->first);
  }

}

void GraphFromText::GetTargetReferences()
{
  if( !m_meshString.empty())
  {
    Group * tmp = this->GetGroupByPath( m_meshString );
    m_mesh = tmp;
    GEOSX_ERROR_IF( m_mesh == nullptr, "The event " << m_meshString << " does not exist or it is not executable." );
  }
}

void GraphFromText::AddEdge(localIndex ind, GraphVertex* v1, GraphVertex* v2)
{
  GraphEdge* e;
  e=new GraphEdge(ind, v1, v2, 0);
  m_edges.push_back(e);
  m_vertex_edges[v1].push_back(e);
  m_vertex_edges[v2].push_back(e);
}

void GraphFromText::RemoveEdge(localIndex ind)
{
  forAll< serialPolicy >( m_edges.size(), [=]( localIndex const i )
  {
    if (m_edges[i]->getIndice()==ind)
    {
      GraphEdge* edge = m_edges[i];
      forAll< serialPolicy >( m_vertex_edges[edge->getN1()].size(), [=]( localIndex const j )
      {
        if (m_vertex_edges[edge->getN1()][j]->getIndice()==ind)
        {
          int test = LvArray::integerConversion< int >(j);
          m_vertex_edges[edge->getN1()].erase(m_vertex_edges[edge->getN1()].begin()+test);
        }
      });
      forAll< serialPolicy >( m_vertex_edges[edge->getN2()].size(), [=]( localIndex const j )
      {
        if (m_vertex_edges[edge->getN2()][j]->getIndice()==ind)
        {
          int test = LvArray::integerConversion< int >(j);
          m_vertex_edges[edge->getN2()].erase(m_vertex_edges[edge->getN2()].begin()+test);
        }
      });


      int test = LvArray::integerConversion< int >(i);
      m_edges.erase(m_edges.begin()+test);
      delete(edge);
    }
  });

}

void GraphFromText::AddVertex(localIndex ind)
{
  GraphVertex* v=new GraphVertex(0,0,ind);
  m_vertices.push_back(v);
  std::vector<GraphEdge*> edges;
  m_vertex_edges.insert({v,edges});
}

void GraphFromText::RemoveVertex(localIndex ind)
{
  GraphVertex* v = getVertexWithIndex(ind);
  std::vector<localIndex> indexes;
  for(long unsigned int i = 0; i<m_vertex_edges[v].size(); i++)
  {
    indexes.push_back(m_vertex_edges[v][i]->getIndice());
  }

  forAll< serialPolicy >( indexes.size(), [=]( localIndex const i )
  {
    std::cout<<indexes[i]<<"\n";
    RemoveEdge(indexes[i]);
  });
  m_vertex_edges.erase(v);
  delete(v);
  
}

GraphVertex* GraphFromText::getVertexWithIndex(localIndex ind)
{
  for(long unsigned int i=0; i < m_vertices.size(); i++)
  {
    if (m_vertices[i]->getIndice()==ind)
    {
      return m_vertices[i];
    }
  }
  return nullptr;
}



void GraphFromText::GenerateGraph()
{
  GetTargetReferences();
  std::ifstream infile(m_file);
  std::string line;
  //std::vector<GraphEdge*> edges;
  while (std::getline(infile, line))
  {
    std::istringstream iss(line);
    int er1, esr1, ei1, er2, esr2, ei2;
    real64 transm;
    if (!(iss >>er1 >> esr1 >> ei1 >> er2 >> esr2 >> ei2 >> transm)) { break; } // error
    //std::cout<<a<<" "<<b<<" "<<c<<"\n" ;
    int place_1=-1;
    int place_2=-1;
    long unsigned int i=0;
    while (i<m_vertices.size() && (place_1==-1 || place_2==-1))
    {
      if (m_vertices[i]->getIndice()==ei1)
      {
        place_1=i;
      }
      if (m_vertices[i]->getIndice()==ei2)
      {
        place_2=i;
      }
      ++i;
    }
    GraphVertex* v1;
    GraphVertex* v2;
    if(place_1==-1)
    {
      v1=new GraphVertex(er1,esr1,ei1);
      m_vertices.push_back(v1);
      std::vector<GraphEdge*> edges;
      m_vertex_edges.insert({v1,edges});
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
      m_vertex_edges.insert({v2,edges});

    }
    else
    {
      v2=m_vertices[place_2];
    }
    GraphEdge* e1= new GraphEdge(m_edges.size(), v1, v2,transm);
    m_vertex_edges[v1].push_back(e1);
    m_vertex_edges[v2].push_back(e1);
    m_edges.push_back(e1);
  }
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
}

REGISTER_CATALOG_ENTRY( GraphBase, GraphFromText, std::string const &, Group * const )
} /* namespace geosx */
