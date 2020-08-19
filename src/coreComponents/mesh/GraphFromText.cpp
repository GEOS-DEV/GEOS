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
  for (long unsigned int i=0;i<m_vertices.size();i++)
  {
    delete(m_vertices[i]);
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
  e=new GraphEdge(ind, v1, v2);
  m_edges.push_back(e);
}

void GraphFromText::RemoveEdge(localIndex ind)
{
  forAll< serialPolicy >( m_edges.size(), [=]( localIndex const i )
  {
    if (m_edges[i]->getIndice()==ind)
    {
      int test = LvArray::integerConversion< int >(i);
      m_edges.erase(m_edges.begin()+test);
    }
  });

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
    int a, b, c;
    if (!(iss >> a >> b >> c)) { break; } // error
    std::cout<<a<<" "<<b<<" "<<c<<"\n" ;
    int place_b=-1;
    int place_c=-1;
    long unsigned int i=0;
    while (i<m_vertices.size() && (place_b==-1 || place_c==-1))
    {
      if (m_vertices[i]->getIndice()==b)
      {
        place_b=i;
      }
      if (m_vertices[i]->getIndice()==c)
      {
        place_c=i;
      }
      ++i;
    }
    GraphVertex* v1;
    GraphVertex* v2;
    if(place_b==-1)
    {
      v1=new GraphVertex(b);
      m_vertices.push_back(v1);
    }
    else
    {
      v1=m_vertices[place_b];
    }
    if(place_c==-1)
    {
      v2=new GraphVertex(c);
      m_vertices.push_back(v2);
    }
    else
    {
      v2=m_vertices[place_c];
    }
    GraphEdge* e1= new GraphEdge(a, v1, v2);
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
}

REGISTER_CATALOG_ENTRY( GraphBase, GraphFromText, std::string const &, Group * const )
} /* namespace geosx */
