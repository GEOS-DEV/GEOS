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
#include "mesh/GraphVertexFace.hpp"
#include "mesh/GraphVertexPoint.hpp"
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
  for (localIndex i=0;i<m_edges.size();i++)
  {
    delete(m_edges[i]);
  }
  for (localIndex i=0;i<m_boundaryEdges.size();i++)
  {
    delete(m_boundaryEdges[i]);
  }

  for(map<std::shared_ptr<GraphVertex>,array1d<GraphEdge*>>::iterator it = m_vertexWithEdgesMap.begin(); it != m_vertexWithEdgesMap.end(); ++it) 
  {
    //delete(it->first);
  }

}


void GraphFromText::AddEdge(localIndex ind, std::shared_ptr<GraphVertex> v1, std::shared_ptr<GraphVertex> v2, real64 transm)
{
  GraphEdge* e;
  e=new GraphEdge(ind, v1, v2, transm);
  m_edges.emplace(m_edges.size(),e);
  m_vertexWithEdgesMap[v1].emplace(m_vertexWithEdgesMap[v1].size(),e);
  m_vertexWithEdgesMap[v2].emplace(m_vertexWithEdgesMap[v2].size(), e);
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
          m_vertexWithEdgesMap[edge->getVertex1()].erase(convertedJ);
        }
      });
      //Delete the edge from the edges linked to the Vertex2
      forAll< serialPolicy >( m_vertexWithEdgesMap[edge->getVertex2()].size(), [=]( localIndex const j )
      {
        if (m_vertexWithEdgesMap[edge->getVertex2()][j]->getEdgeIndex()==ind)
        {
          int convertedJ = LvArray::integerConversion< int >(j);
          m_vertexWithEdgesMap[edge->getVertex2()].erase(convertedJ);
        }
      });
      //Delete the edge from the edge list
      int convertedI = LvArray::integerConversion< int >(i);
      m_edges.erase(convertedI);
      delete(edge);
    }
  });


  forAll< serialPolicy >( m_boundaryEdges.size(), [=]( localIndex const i )
  {
    if (m_boundaryEdges[i]->getEdgeIndex()==ind)
    {
      GraphEdge* edge = m_boundaryEdges[i];
      //Delete the edge from the edges linked to the Vertex1
      forAll< serialPolicy >( m_vertexWithEdgesMap[edge->getVertex1()].size(), [=]( localIndex const j )
      {
        if (m_vertexWithEdgesMap[edge->getVertex1()][j]->getEdgeIndex()==ind)
        {
          int convertedJ = LvArray::integerConversion< int >(j);
          m_vertexWithEdgesMap[edge->getVertex1()].erase(convertedJ);
        }
      });
      //Delete the edge from the edges linked to the Vertex2
      forAll< serialPolicy >( m_vertexWithEdgesMap[edge->getVertex2()].size(), [=]( localIndex const j )
      {
        if (m_vertexWithEdgesMap[edge->getVertex2()][j]->getEdgeIndex()==ind)
        {
          int convertedJ = LvArray::integerConversion< int >(j);
          m_vertexWithEdgesMap[edge->getVertex2()].erase(convertedJ);
        }
      });
      //Delete the edge from the edge list
      int convertedI = LvArray::integerConversion< int >(i);
      m_boundaryEdges.erase(convertedI);
      delete(edge);
    }
  });


}

void GraphFromText::AddVertex(localIndex er, localIndex esr, globalIndex ei)
{
  std::shared_ptr<GraphVertex> v=std::make_shared<GraphVertex>(GraphVertex(er,esr,ei));
  m_vertices.emplace(m_vertices.size(), v);
  array1d<GraphEdge*> edges;
  m_vertexWithEdgesMap.insert({v,edges});
}

void GraphFromText::RemoveVertex(localIndex er, localIndex esr, globalIndex ei)
{
  std::shared_ptr<GraphVertex> v = getVertexWithGlobalIndex(er,esr,ei);
  array1d<localIndex> indexes;
  localIndex loc = -1;
  for(localIndex i = 0; i<m_vertices.size(); i++)
  {
    if (m_vertices[i]->getRegionIndex()==er && m_vertices[i]->getSubRegionIndex()==esr && m_vertices[i]->getGlobalVertexIndex()==ei)
    {
      loc = i;
    }
  }
  if (loc != -1)
  {
    m_vertices.erase(loc);    
  }
  for(localIndex i = 0; i<m_vertexWithEdgesMap[v].size(); i++)
  {
    indexes.emplace(indexes.size(), m_vertexWithEdgesMap[v][i]->getEdgeIndex());
  }

  forAll< serialPolicy >( indexes.size(), [=]( localIndex const i )
  {
    //Deleting all edges linked to the vertex to delete (as they will have a void target on one extremity)
    //std::cout<<indexes[i]<<"\n";
    RemoveEdge(indexes[i]);
  });
  m_vertexWithEdgesMap.erase(v);
  //delete(v);
  
}

void GraphFromText::RemoveVertex(std::shared_ptr<GraphVertex> vertex)
{
  //std::cout<< vertex->getGlobalVertexIndex() << "/n";
  array1d<localIndex> indexes;
  localIndex loc = -1;
  for(localIndex i = 0; i<m_vertices.size(); i++)
  {
    if (m_vertices[i]->getRegionIndex()==vertex->getRegionIndex() && m_vertices[i]->getSubRegionIndex()==vertex->getSubRegionIndex() && m_vertices[i]->getGlobalVertexIndex()==vertex->getGlobalVertexIndex())
    {
      loc = i;
    }
  }
  if (loc != -1)
  {
    m_vertices.erase(loc);    
  }

  for(localIndex i = 0; i<m_vertexWithEdgesMap[vertex].size(); i++)
  {
    indexes.emplace(indexes.size(), m_vertexWithEdgesMap[vertex][i]->getEdgeIndex());
  }

  forAll< serialPolicy >( indexes.size(), [=]( localIndex const i )
  {
    //Deleting all edges linked to the vertex to delete (as they will have a void target on one extremity)
    //std::cout<<indexes[i]<<"\n";
    RemoveEdge(indexes[i]);
  });
  m_vertexWithEdgesMap.erase(vertex);
  //delete(vertex);
  
}


std::shared_ptr<GraphVertex> GraphFromText::getVertexWithGlobalIndex(localIndex er, localIndex esr, globalIndex ei)
{
  for(localIndex i=0; i < m_vertices.size(); i++)
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
  std::string type;
  int size;
  int count = 0;
  infile >> type >> size;
  //std::cout<<type<<" "<<size<<std::endl;
  std::getline(infile, line);
  //array1d<GraphEdge*> edges;
  while (std::getline(infile, line) && count< size)
  {
    std::istringstream iss(line);
    int er1, esr1, ei1, er2, esr2, ei2;
    real64 transm;
    //The input file is written with thr triple index (regionIndex, subRegionIndex,globalVertexIndex) of the two connected point and associated transmissibility
    if (!(iss >>er1 >> esr1 >> ei1 >> er2 >> esr2 >> ei2 >> transm)) { break; } // error
    count++;
    int place_1=-1;
    int place_2=-1;
    localIndex i=0;
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
    std::shared_ptr<GraphVertex> v1;
    std::shared_ptr<GraphVertex> v2;
    if(place_1==-1)
    {
      v1=std::make_shared<GraphVertex>(GraphVertex(er1,esr1,ei1));
      m_vertices.emplace(m_vertices.size(), v1);
      array1d<GraphEdge*> edges;
      m_vertexWithEdgesMap.insert({v1,edges});
    }
    else
    {
      v1=m_vertices[place_1];
    }
    if(place_2==-1)
    {
      v2=std::make_shared<GraphVertex>(GraphVertex(er2,esr2,ei2));
      m_vertices.emplace(m_vertices.size(), v2);
      array1d<GraphEdge*> edges;
      m_vertexWithEdgesMap.insert({v2,edges});

    }
    else
    {
      v2=m_vertices[place_2];
    }
    //Create the edge, add it to the edge list
    GraphEdge* e= new GraphEdge(m_edges.size(), v1, v2,transm);
    m_edges.emplace(m_edges.size(), e);
    //Add the edge to the list of linked edge for his two vertices
    m_vertexWithEdgesMap[v1].emplace(m_vertexWithEdgesMap[v1].size(), e);
    m_vertexWithEdgesMap[v2].emplace(m_vertexWithEdgesMap[v2].size(), e);
    }
  infile >> type >> size;
  std::getline(infile, line);
  count = 0;
  //array1d<GraphEdge*> edges;
  std::cout<<"test"<<std::endl;
  while (std::getline(infile, line) && count< size)
  {
    std::istringstream iss(line);
    int er1, esr1, ei1, er2, esr2, ei2;
    real64 transm;
    //The input file is written with thr triple index (regionIndex, subRegionIndex,globalVertexIndex) of the two connected point and associated transmissibility
    if (!(iss >>er1 >> esr1 >> ei1 >> er2 >> esr2 >> ei2 >> transm)) { break; } // error
    count++;
    int place_1=-1;
    int place_2=-1;
    localIndex i=0;
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
    std::shared_ptr<GraphVertex> v1;
    std::shared_ptr<GraphVertexFace> v2;
    if(place_1==-1)
    {
      v1=std::make_shared<GraphVertex>(GraphVertex(er1,esr1,ei1));
      m_vertices.emplace(m_vertices.size(), v1);
      array1d<GraphEdge*> edges;
      m_vertexWithEdgesMap.insert({v1,edges});
    }
    else
    {
      v1=m_vertices[place_1];
    }
    if(place_2==-1)
    {
      //std::cout<<er2 << " " << esr2 << " " << ei2 << "\n";
      v2=std::make_shared<GraphVertexFace>(GraphVertexFace(er2,esr2,ei2));
      m_vertices.emplace(m_vertices.size(), v2);
      array1d<GraphEdge*> edges;
      m_vertexWithEdgesMap.insert({v2,edges});
      

    }
    else
    {
      v2= std::dynamic_pointer_cast<GraphVertexFace>(m_vertices[place_2]);
    }
    //Create the edge, add it to the edge list
    GraphEdge* e= new GraphEdge(m_edges.size()+m_boundaryEdges.size(), v1, v2,transm);
    m_boundaryEdges.emplace(m_boundaryEdges.size(), e);
    //Add the edge to the list of linked edge for his two vertices
    m_vertexWithEdgesMap[v1].emplace(m_vertexWithEdgesMap[v1].size(), e);
    m_vertexWithEdgesMap[v2].emplace(m_vertexWithEdgesMap[v2].size(), e);

  } 
}

void GraphFromText::PartitionGraph(const MeshLevel & mesh)
{
  /*
  for( int i=0; i<1000; i++)
  {
    m_vertices[i]->setLocalIndex(m_vertices[i]->getGlobalVertexIndex());
    m_vertices[i]->setGhostIndex(-1);
  }
  */
  ElementRegionManager const & elemManager = *mesh.getElemManager();
  ElementRegionManager::ElementViewAccessor< arrayView1d< integer const > > const elemGhostRank =
    elemManager.ConstructArrayViewAccessor< integer, 1 >( ObjectManagerBase::viewKeyStruct::ghostRankString );
  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > const localToGlobalMap =
    elemManager.ConstructArrayViewAccessor< globalIndex, 1 >( ObjectManagerBase::viewKeyStruct::localToGlobalMapString );
  //Create and populate a vector containing the vertices to delete as they are not present in the partition
  array1d <std::shared_ptr<GraphVertex>> vertexToDelete;
  array1d <std::shared_ptr<GraphVertex>> vertexToDeleteLater;
  for( localIndex i=0; i<m_vertices.size(); i++)
  { 
    if (m_vertices[i]->getRegionIndex()!=-1)
    {
      vertexToDelete.emplace(vertexToDelete.size(), m_vertices[i]);
    }
    else
    {
      vertexToDeleteLater.emplace(vertexToDeleteLater.size(), m_vertices[i]);
    }
  }

  for (long int g = 0; g<elemGhostRank.size(); g++)
  {
    for (long int h = 0; h<elemGhostRank[g].size(); h++)
    {
      for (long int i = 0; i<elemGhostRank[g][h].size(); i++)
      {
        int vertexToKeep = -1;
        for (localIndex j = 0; j<vertexToDelete.size();j++)
        {
          //For each vertex appearing in the elemGhostRank matrix, hence owned or ghosted, we suppress it from the list of point to delete
          //std::cout << localToGlobalMap[g][h][i]<< " " << vertexToDelete[j]->getGlobalVertexIndex() << "\n";
          if (localToGlobalMap[g][h][i] ==vertexToDelete[j]->getGlobalVertexIndex())
          {
            vertexToKeep=j;
            //The values are copied in the vertex object to avoid requesting the local to global conversion again
            vertexToDelete[j]->setGhostIndex(elemGhostRank[g][h][i]);
            vertexToDelete[j]->setLocalIndex(i);
            //std::cout<<vertexToDelete[j]->getLocalVertexIndex()<<" ";
            //std::cout<<vertexToDelete[j]->getGhostIndex()<<"\n";
          }
        }
        //If we found a corresponding vertex we suppress it from the list to vertex to delete
        if (vertexToKeep != -1)
        {
          vertexToDelete.erase(vertexToKeep);
        }
      }
    }
  }

  //Every point that have not been removed from the list is not present in the partition and is deleted
  std::cout<<vertexToDelete.size()<<std::endl;
  for ( localIndex i=0; i<vertexToDelete.size(); i++)
  {
    int to_delete = -1;
    for (localIndex j = 0; j<m_vertices.size();j++)
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
  
  RemapFace(mesh);
  FaceManager const & faceManager = *mesh.getFaceManager();
  arrayView1d< integer const > const & faceGhostRank = faceManager.ghostRank(); 
  for (localIndex i=0; i<vertexToDeleteLater.size();i++)
  {
    std::shared_ptr<GraphVertexFace> face = std::dynamic_pointer_cast<GraphVertexFace>(vertexToDeleteLater[i]);
    if (face->getCorrespondingId()==-1)
    {
      RemoveVertex(face);
    }
    else
    {
      face->setGhostIndex(faceGhostRank[face->getCorrespondingId()]);
    }
  } 
}

void GraphFromText::RemapFace(const MeshLevel & mesh)
{
  FaceManager const & faceManager = *mesh.getFaceManager();
  arrayView2d< localIndex const > const & elemList = faceManager.elementList();
  //arrayView1d< integer const > const & faceGhostRank = faceManager.ghostRank(); 

  for(localIndex h = 0; h < elemList.size(); h++)
  {

  
    //GEOSX_ERROR_IF_GT( elemList[ h ][ 0 ], 1e9 );
    if (h>=2000)
    {
      std::cout<<"Error "<<h<<" "<<elemList[h][0]<<" "<<"\n";
    }
    if (h==0)
    {
      std::cout<<h<<" "<<elemList[h][0]<<"\n";
    }

    // Filter in boundary faces
    if( elemList[h][1] < 0 && elemList[h][0]< 1e9)
    {
      bool found = false;
            for (localIndex i = 0; i < m_boundaryEdges.size(); i++)
      { 
        if (elemList[h][0] == m_boundaryEdges[i]->getVertex1()->getLocalVertexIndex() && !found)
        { 
          std::shared_ptr<GraphVertexFace> face;
          face = std::dynamic_pointer_cast<GraphVertexFace>(m_boundaryEdges[i]->getVertex2());
          if (face==0) {std::cout << "Null pointer on second type-cast.\n";}
          
          if (face->getCorrespondingId() == -1)
          {
            //std::cout<< h << " " << m_boundaryEdges[i]->getEdgeIndex() << "\n";
            //std::cout<< elemList[h][0] << " " << m_boundaryEdges[i]->getVertex1()->getLocalVertexIndex() << " " << m_boundaryEdges[i]->getVertex2()->getLocalVertexIndex() << " " << m_boundaryEdges[i]->getVertex2()->getSubRegionIndex() << " " << m_boundaryEdges[i]->getVertex2()->getRegionIndex() << "\n\n";
            face->setCorrespondingId(h);
            found = true;
            //std::cout<< face->getCorrespondingId()<< "\n";
          }
        }
      }
    }
  }
}


REGISTER_CATALOG_ENTRY( GraphBase, GraphFromText, std::string const &, Group * const )
} /* namespace geosx */
