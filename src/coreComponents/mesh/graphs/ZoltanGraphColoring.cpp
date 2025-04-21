/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License: LGPL-2.1-only
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
 * @file ZoltanGraphColoring.cpp
 */

#include "ZoltanGraphColoring.hpp"
#include "GraphToolsMPI.hpp"
#include <algorithm>


#if 0

#define GEOS_ZOLTAN_CHECK( call ) \
  do { \
    auto const ierr = call; \
    GEOS_ERROR_IF_NE_MSG( ierr, ZOLTAN_OK, "Error in call to:\n" << #call ); \
  } while( false )
#endif

namespace geos
{
namespace graph
{

ZoltanGraphColoring::ZoltanGraphColoring( MPI_Comm comm ): GraphColoringBase( comm )
{
  float version;
  Zoltan_Initialize( 0, nullptr, &version );

  m_zz = new Zoltan( comm );
  if( !m_zz )
  {
    exit( 0 );
  }

  m_zz->Set_Param( "DEBUG_LEVEL", "0" );
  m_zz->Set_Param( "NUM_GID_ENTRIES", "1" );
  m_zz->Set_Param( "NUM_LID_ENTRIES", "1" );
  m_zz->Set_Param( "OBJ_WEIGHT_DIM", "0" );
  m_zz->Set_Param( "COLORING_PROBLEM", "distance-1" );
}

//Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "1"); /* weights are 1 float */
/*		  { "SUPERSTEP_SIZE",     NULL, "INT", 0},
      { "COMM_PATTERN",       NULL, "CHAR", 0 },
      { "VERTEX_VISIT_ORDER", NULL, "CHAR", 0 },
      { "COLORING_METHOD",    NULL, "CHAR", 0},
 */


ZoltanGraphColoring::~ZoltanGraphColoring()
{
  delete m_zz;
}


int ZoltanGraphColoring::colorGraph( const std::vector< camp::idx_t > & localAdjncy )
{
  std::vector< camp::idx_t > localXadj = createXadjFromAdjncy( localAdjncy, m_comm );
  std::vector< int > colors = colorGraph( localXadj, localAdjncy );
  return colors[0];
}


std::vector< int > ZoltanGraphColoring::colorGraph( const std::vector< camp::idx_t > & xadj,
                                                    const std::vector< camp::idx_t > & adjncy )
{
  //int const rank = MpiWrapper::commRank( m_comm );
  int rank;
  MPI_Comm_rank( m_comm, &rank );


  ZoltanGraph graph;
  graph.m_xadj.assign( xadj.begin(), xadj.end());
  graph.m_adjncy.assign( adjncy.begin(), adjncy.end());
  graph.m_numVertices = xadj.size() - 1;
  graph.m_rank = rank;

  std::vector< int > vertexGID = createVertexGlobalID( xadj, m_comm );
  graph.m_vertexGID.resize( graph.m_numVertices );
  for( int i = 0; i < graph.m_numVertices; i++ )
  {
    graph.m_vertexGID[i] = vertexGID[i];
  }

  /* set call backs */
  m_zz->Set_Num_Obj_Fn( getNumberOfVertices, &graph );
  m_zz->Set_Obj_List_Fn( getVertexList, &graph );
  m_zz->Set_Num_Edges_Multi_Fn( getNumEdgesList, &graph );
  m_zz->Set_Edge_List_Multi_Fn( getEdgeList, &graph );


  int numGidEntries = 1;
  int numReqObjs = graph.m_numVertices;
  ZOLTAN_ID_PTR reqObjs = new ZOLTAN_ID_TYPE[graph.m_numVertices];
  for( int i = 0; i < graph.m_numVertices; i++ )
  {
    reqObjs[i] = graph.m_vertexGID[i];
  }

  int * color = new int[graph.m_numVertices];
  std::fill( color, color + graph.m_numVertices, -1 );

  //GEOS_ZOLTAN_CHECK( m_zz->Color( numGidEntries, numReqObjs, reqObjs, color ));
  m_zz->Color( numGidEntries, numReqObjs, reqObjs, color );
  std::vector< int > coloringVector;
  coloringVector.assign( color, color + graph.m_numVertices );

  delete[] reqObjs;
  delete[] color;

  return coloringVector;
}

int ZoltanGraphColoring::getNumberOfVertices( void * data, int * ierr )
{
  ZoltanGraph * graph = static_cast< ZoltanGraph * >(data);
  *ierr = ZOLTAN_OK;
  return graph->m_numVertices;
}

void ZoltanGraphColoring::getVertexList( void * data, int GEOS_UNUSED_PARAM( sizeGID ), int GEOS_UNUSED_PARAM( sizeLID ),
                                         ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                                         int GEOS_UNUSED_PARAM( wgtDim ), float * GEOS_UNUSED_PARAM( objWgts ), int * ierr )
{
  ZoltanGraph * graph = static_cast< ZoltanGraph * >(data);
  *ierr = ZOLTAN_OK;

  for( int i = 0; i < graph->m_numVertices; ++i )
  {
    globalID[i] = graph->m_vertexGID[i];
    localID[i] = i;
  }
}

void ZoltanGraphColoring::getNumEdgesList( void * data, int sizeGID, int sizeLID,
                                           int numObj, ZOLTAN_ID_PTR GEOS_UNUSED_PARAM( globalID ), ZOLTAN_ID_PTR localID,
                                           int * numEdges, int * ierr )
{
  ZoltanGraph * graph = static_cast< ZoltanGraph * >(data);
  *ierr = ZOLTAN_OK;

  if( sizeGID != 1 || sizeLID != 1 || numObj != graph->m_numVertices )
  {
    *ierr = ZOLTAN_FATAL;
    return;
  }

  for( int i = 0; i < numObj; ++i )
  {
    int idx = localID[i];
    numEdges[i] = graph->m_xadj[idx + 1] - graph->m_xadj[idx];
  }
}

void ZoltanGraphColoring::getEdgeList( void * data, int sizeGID, int sizeLID,
                                       int numObj, ZOLTAN_ID_PTR GEOS_UNUSED_PARAM( globalID ), ZOLTAN_ID_PTR localID,
                                       int * numEdges, ZOLTAN_ID_PTR adjncy, int * nborProc,
                                       int wgtDim, float * GEOS_UNUSED_PARAM( ewgts ), int * ierr )
{
  ZoltanGraph * graph = static_cast< ZoltanGraph * >(data);
  *ierr = ZOLTAN_OK;

  if( sizeGID != 1 || sizeLID != 1 || numObj != graph->m_numVertices || wgtDim != 0 )
  {
    *ierr = ZOLTAN_FATAL;
    return;
  }

  ZOLTAN_ID_TYPE * nextNbor = adjncy;
  int * nextProc = nborProc;

  for( int i = 0; i < numObj; ++i )
  {
    int from = graph->m_xadj[localID[i]] - graph->m_xadj[localID[0]];
    int to = graph->m_xadj[localID[i] + 1] - graph->m_xadj[localID[0]];

    if((to - from) != numEdges[i] )
    {
      *ierr = ZOLTAN_FATAL;
      return;
    }

    for( int j = from; j < to; ++j )
    {
      *nextNbor++ = graph->m_adjncy[j];
      *nextProc++ = graph->m_rank;
    }
  }
}


size_t ZoltanGraphColoring::getNumberOfColors( const int color ) const
{
  return GraphColoringBase::getNumberOfColors( color, m_comm );
}

size_t ZoltanGraphColoring::getNumberOfColors( const std::vector< int > & colors ) const
{
  return GraphColoringBase::getNumberOfColors( colors, m_comm );
}


bool ZoltanGraphColoring::isColoringValid( const std::vector< camp::idx_t > & adjncy, const int color ) const
{
  return GraphColoringBase::isColoringValid( adjncy, color, m_comm );
}

} // namespace graph
} // namespace geos
