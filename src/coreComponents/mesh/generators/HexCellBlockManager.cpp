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

#include "HexCellBlockManager.hpp"
#include "CellBlockUtilities.hpp" // TODO Remove
#include "CellBlock.hpp"

#include "dataRepository/Group.hpp"

#include <algorithm>
#include <numeric>

#pragma clang diagnostic ignored "-Wunused" // TODO Remove this - too much of a pain when developing new features
#pragma clang diagnostic ignored "-Wcomment"
#pragma clang diagnostic ignored "-Wsign-compare"


namespace geosx
{
using namespace dataRepository;

//Typedefs for the awful datastructures 
typedef array2d<localIndex, cells::NODE_MAP_PERMUTATION> CellVertexIndices;

typedef localIndex LocalVertexIndex;
typedef globalIndex GlobalVertexIndex;

typedef localIndex LocalCellIndex;
typedef globalIndex GlobalCellIndex;

typedef localIndex LocalFaceIndex;  // A face in a cell 6*cell + f
typedef localIndex LocalEdgeIndex;  // An edge in a cell 12*cell + e

typedef localIndex CellBlockIndex;

typedef localIndex SizeOfStuff;


// Structure for face computation
// Possible optimisation - using a similar strategies to edges
struct LocalIndexQuadruplet{
  localIndex v[4];
};
struct {
  bool operator()( LocalIndexQuadruplet const & a, LocalIndexQuadruplet const & b) const {
    if(a.v[0] != b.v[0]) return a.v[0] < b.v[0];
    if(a.v[1] != b.v[1]) return a.v[1] < b.v[1];
    if(a.v[2] != b.v[2]) return a.v[2] < b.v[2];
    return false;
  }
} sortV0V1V2;


// Storage of mesh edges
typedef std::pair< LocalVertexIndex, LocalEdgeIndex > EdgeInfo;

struct {
  bool operator() ( EdgeInfo const & a, EdgeInfo const & b ) const
  {
    return a.first < b.first;
  }
} compareEdges;

struct {
  bool operator() ( EdgeInfo const & a, EdgeInfo const & b ) const
  {
    return a.first == b.first;
  }
} equalEdges;

// Use explicit aliases for indices
// HEX
typedef unsigned int HexVertexIndex;      // a vertex in a hex 0 to 8
typedef unsigned int HexFacetVertexIndex; // a vertex in a hex facet 0 to 4
typedef unsigned int HexEdgeIndex;        // a edge in a hex 0 to 12
typedef unsigned int HexFacetIndex;       // a facet in a hex 0 to 6

/*  Hexahedron template
*
*   7----------6
*   |\         |\
*   | \        | \
*   |  \       |  \
*   |   4------+---5
*   |   |      |   |
*   3---+------2   |
*    \  |       \  |
*     \ |        \ |
*      \|         \|
*       0----------1
*/
struct Hex
{
static const unsigned int nbVertices = 8;
static const unsigned int nbEdges = 12;
static const unsigned int nbFacets = 6;
static const unsigned int nbTriFacets = 0;
static const unsigned int nbQuadFacets = 6;
static const unsigned int nbEdgesPerFacet = 4;


static constexpr HexVertexIndex facetVertex[6][4] = {
    {0, 1, 2, 3}, {4, 5, 6, 7}, {0, 1, 5, 4}, {1, 2, 6, 5}, {3, 2, 6, 7}, {0, 3, 7, 4}};
//      0            1               2             3              4             5
static constexpr HexVertexIndex edgeVertex[12][2]{
    {0, 1}, {0, 3}, {0, 4}, {1, 2}, {1, 5}, {2, 3}, {2, 6}, {3, 7}, {4, 5}, {4, 7}, {5, 6}, {6, 7}};

// TODO to check
static constexpr HexVertexIndex edgeAdjacentFacet[12][2]{
    {0, 2}, {0, 5}, {2, 5}, {0, 3}, {2, 3}, {0, 4}, {3, 4}, {5, 4}, {2, 1}, {1, 5}, {1, 3}, {1, 4}};

static constexpr HexVertexIndex facetEdgesVertex[6][4][2] = {
    {{0, 1}, {1, 2}, {2, 3}, {3, 0}}, {{4, 5}, {5, 6}, {6, 7}, {7, 4}},
    {{0, 1}, {1, 5}, {5, 4}, {4, 0}}, {{1, 2}, {2, 6}, {6, 5}, {5, 1}}, 
    {{3, 2}, {2, 6}, {6, 7}, {7, 3}}, {{0, 3}, {3, 7}, {7, 4}, {4, 0}}
};

static constexpr HexFacetIndex vertexAdjacentFacet[8][3] = {
    {0, 2, 5}, {0, 2, 3}, {0, 3, 4}, {0, 4, 5}, {1, 2, 5}, {1, 2, 3}, {1, 3, 4}, {1, 4, 5}};

static constexpr HexVertexIndex orientedFacetVertex[6][4] = {
    {0, 1, 2, 3}, {4, 7, 6, 5}, {0, 4, 5, 1}, {1, 5, 6, 2}, {3, 2, 6, 7}, {0, 3, 7, 4}};
};


class HexMeshConnectivityBuilder
{
public:

  HexMeshConnectivityBuilder( HexCellBlockManager & cellBlockManager );
  ~HexMeshConnectivityBuilder() = default;

  ArrayOfSets<localIndex>          getNodeToEdges   ();
  ArrayOfSets<localIndex>          getNodeToFaces   ();
  ArrayOfArrays<localIndex>        getNodeToElements();

  array2d<geosx::localIndex>       getEdgeToNodes   ();
  ArrayOfSets<geosx::localIndex>   getEdgeToFaces   ();
  
  ArrayOfArrays<localIndex>        getFaceToNodes   ();
  array2d<localIndex>              getFaceToElements();
  ArrayOfArrays<geosx::localIndex> getFaceToEdges   ();

  // Compute functions for all maps
  // Do we really need to  implement it
  // What is the best order to compute this
  void computeAllMaps() {
    computeVertexGlobalToLocal();
    computeCellNeighbors();
    computeEdges();
  };

  bool isRegularMesh() { return false; }

  // Why is this here? The information is on the CellBlocks
  bool computeElementsToFaces(); 
  bool computeElementsToEdges();

private:
  bool computeVertexGlobalToLocal(); // To remove ?
  bool computeCellNeighbors();
  bool computeEdges(); 
  
  bool computeNodesToEdges   ();
  bool computeNodesToFaces   ();
  bool computeNodesToElements();

  bool computeEdgesToFaces   ();
  bool computeEdgesToElements() { return false; } // Not implemented - Not used priori used by no one
  
  bool computeFacesToNodes   ();
  bool computeFacesToElements();
  bool computeFacesToEdges   ();

  std::pair< CellBlockIndex, LocalCellIndex> 
  getBlockCellFromManagerCell( LocalCellIndex cellId ) const;

  unsigned int nbCellBlocks() {
    return m_cellBlocks.size();
  } 
  CellBlock const & getCellBlock( unsigned int id ) const {
    return *(m_cellBlocks[id]);
  }
  CellBlock & getCellBlock( unsigned int id ){
    return *(m_cellBlocks[id]);
  }

private:
// HexCellBlockManager for which we build the mapping
HexCellBlockManager const &  m_cellBlockManager;

SizeOfStuff nbNodes = 0;
SizeOfStuff nbEdges = 0;
SizeOfStuff nbFaces = 0;
SizeOfStuff nbElements = 0;

std::vector< CellBlock * > m_cellBlocks;
std::vector< LocalCellIndex > m_blockCellIndexOffset;

// Do we need it ?
// Or areall indices local and we are good
std::unordered_map< GlobalVertexIndex, LocalVertexIndex > vertexGlobalToLocal; 


// Cell adjacencies - for all cell blocks
// I do not understand LvArray - and even less multidimensional LvArrays
// The advantages of multidimensional tables are a mystery to me
std::vector< LocalFaceIndex > m_cellNeighbors;
std::vector< localIndex > m_uniqueFaces;
std::vector< bool > m_isBoundaryFace;

std::vector< EdgeInfo > m_edges;
std::vector< localIndex > m_uniqueEdgeIds;

// Some of the mappings that this class is responsible of building
// TODO Why are storage strategies different? 
ArrayOfSets<localIndex> m_nodeToEdges;
ArrayOfSets<localIndex> m_nodeToFaces;
ArrayOfArrays<localIndex> m_nodeToElements;

ArrayOfSets<localIndex> m_edgeToFaces;
array2d<localIndex> m_edgeToNodes;

ArrayOfArrays<localIndex> m_faceToNodes;
ArrayOfArrays<localIndex> m_faceToEdges;
array2d<localIndex> m_faceToElements;

};

// Be able to get out and return an error message if 3 cells have the same facet
// Flagging boundary facet if the MPI rank should be possible at this point

HexMeshConnectivityBuilder::HexMeshConnectivityBuilder( HexCellBlockManager & cellBlockManager ):
  m_cellBlockManager(cellBlockManager)
{
  nbNodes = cellBlockManager.numNodes();

  Group & group =  cellBlockManager.getCellBlocks();
  localIndex nbBlocks = group.numSubGroups();

  m_cellBlocks.resize(nbBlocks);
  m_blockCellIndexOffset.resize(nbBlocks);

  localIndex prevSum = 0;
  for (int i = 0; i < nbBlocks; ++i)
  {
    // TODO - Check that we indeed have CellBlock instantiation 
    m_cellBlocks[i] = &(group.getGroup< CellBlock >( i ));
    m_blockCellIndexOffset[i] = prevSum + m_cellBlocks[i]->numElements();
  }

  nbElements =  m_blockCellIndexOffset[nbBlocks];
}

ArrayOfArrays<localIndex> HexMeshConnectivityBuilder::getFaceToNodes()
{
  return m_faceToNodes;
}
array2d<geosx::localIndex> HexMeshConnectivityBuilder::getEdgeToNodes()
{
  return m_edgeToNodes;
}
ArrayOfSets<localIndex> HexMeshConnectivityBuilder::getNodeToEdges()
{
  return m_nodeToEdges;
}
ArrayOfSets<localIndex> HexMeshConnectivityBuilder::getNodeToFaces()
{
  return m_nodeToFaces;
}
ArrayOfArrays<localIndex> HexMeshConnectivityBuilder::getNodeToElements()
{
  return m_nodeToElements;
}
array2d<localIndex> HexMeshConnectivityBuilder::getFaceToElements()
{
  return m_faceToElements;
}
ArrayOfArrays<geosx::localIndex> HexMeshConnectivityBuilder::getFaceToEdges()
{
  return m_faceToEdges;
}
ArrayOfSets<geosx::localIndex> HexMeshConnectivityBuilder::getEdgeToFaces()
{
  return m_edgeToFaces;
}

std::pair< CellBlockIndex, LocalCellIndex> 
HexMeshConnectivityBuilder::getBlockCellFromManagerCell( LocalCellIndex cellId ) const
{
  CellBlockIndex b = 0;
  while( cellId > m_blockCellIndexOffset[b] )
  {
    b++;
  }
  assert(b < nbCellBlocks()); // debug paranoia

  localIndex offset = b > 0 ? m_blockCellIndexOffset[b-1] : 0;
  LocalCellIndex c = cellId - offset;

  return std::pair< CellBlockIndex, LocalCellIndex> (b, c);
}

bool HexMeshConnectivityBuilder::computeVertexGlobalToLocal()
{
  vertexGlobalToLocal.reserve(nbNodes);
  array1d<GlobalVertexIndex> globalIds = m_cellBlockManager.getNodeLocalToGlobal();
  for (LocalVertexIndex i = 0; i < nbNodes; ++i)
  {
    GlobalVertexIndex id = globalIds[i];
    vertexGlobalToLocal[id] = i;
  }
  return true;
}

/* Compute and store the cell neighors for each face
 */
bool HexMeshConnectivityBuilder::computeCellNeighbors()
{
  // 1 - Allocate
  SizeOfStuff nbTotalFaces = nbElements * Hex::nbFacets;
  m_cellNeighbors.resize( nbTotalFaces);
  
  // To collect and sort the facets 
  std::vector< LocalIndexQuadruplet > allFaces ( nbTotalFaces );

  // 2 - Fill 
  localIndex curFace = 0;
  for (int i = 0; i < nbCellBlocks(); ++i)
  {
    CellBlock const & block = getCellBlock(i);
    CellVertexIndices const & cells = block.getElemToNodes();

    for (int j = 0; j < cells.size(); ++j)
    {
      for (int f = 0; f < Hex::nbFacets; ++f)
      {
        // TODO Let's bet that we get VertexIds locally to the partition 
        LocalVertexIndex v0 = cells( j, Hex::facetVertex[f][0] );
        LocalVertexIndex v1 = cells( j, Hex::facetVertex[f][1] );
        LocalVertexIndex v2 = cells( j, Hex::facetVertex[f][2] );
        LocalVertexIndex v3 = cells( j, Hex::facetVertex[f][3] );

        // Sort the vertices 
        if( v0 > v1 ) std::swap( v0, v1 );
        if( v2 > v3 ) std::swap( v2, v3 );
        if( v0 > v2 ) std::swap( v0, v2 );
        if( v1 > v3 ) std::swap( v1, v3 );
        if( v1 > v2 ) std::swap( v1, v2 );

        allFaces[curFace].v[0] = v0;
        allFaces[curFace].v[1] = v1;
        allFaces[curFace].v[2] = v2;
        // If the mesh is valid, then if 2 quad faces share 3 vertices they are the same
        // Last slot is used to identify the facet from its cell 
        allFaces[curFace].v[3]= Hex::nbFacets * j + f; 
        
        curFace++;
      }
    }
  }

  // 3 - Sort 
  // TODO Definitely not the fastest we can do - Use of HXTSort if possible? LvArray is fast? 
  std::sort( allFaces.begin(), allFaces.end(), sortV0V1V2);

  // That is an overallocation - about twice as big as needed
  m_uniqueFaces.resize( allFaces.size() ); 
  m_isBoundaryFace.resize( allFaces.size(), false ); 


  // 4 - Counting + Set Cell Adjacencies
  curFace = 0;
  int i = 0;
  for( ; i+1 < nbTotalFaces; ++i )
  {
    LocalIndexQuadruplet f = allFaces[i];
    LocalIndexQuadruplet f1 = allFaces[i+1];

    // If two successive faces are the same - this is an interior facet
    if( f.v[0]==f1.v[0] && f.v[1]==f1.v[1] && f.v[2]==f1.v[2] )
    {
      m_cellNeighbors[f.v[3]] = f1.v[3];
      m_cellNeighbors[f1.v[3]] = f.v[3];
      m_uniqueFaces[curFace] = f.v[3];
      curFace++;
      ++i;
    }
    // If not this is a boundary face
    else
    {
      m_cellNeighbors[f.v[3]] = -1; //Replace by NO_ID - maxvalue of  something
      m_uniqueFaces[curFace] = f.v[3];
      m_isBoundaryFace[curFace] = true;
      curFace++;
    }
  }
  // Last facet is a boundary facet 
  if (i < nbTotalFaces) 
  {
    m_uniqueFaces[curFace] = allFaces[i].v[3];
    m_isBoundaryFace[curFace] = true;
    curFace++;
  }
  // Set the number of faces 
  nbFaces = curFace; 
  // Cut off the non filled values 
  m_uniqueFaces.resize( nbFaces );
  m_isBoundaryFace.resize( nbFaces );

  return true;
}

// Compute Cell Neighbors MUST have been called beforehand
bool HexMeshConnectivityBuilder::computeFacesToNodes()
{
  // 1 - Allocate - No overallocation
  m_faceToNodes.resize(0);

  // LvArray allocation is an unclear business
  m_faceToNodes.reserve( nbFaces );
  m_faceToNodes.reserveValues( 4 * nbFaces );
  for (int i = 0; i < nbFaces; ++i)
  {
    m_faceToNodes.appendArray( 4 );
  }

  // 2 - Fill FaceToNode  -- Could be avoided and done when required from adjacencies
  for( int curFace = 0; curFace < m_uniqueFaces.size(); ++curFace)
  {
    localIndex f = m_uniqueFaces[curFace];
    LocalCellIndex c = f / Hex::nbFacets;

    // We want the nodes of face f % 6 in cell c 
    auto blockCell = getBlockCellFromManagerCell(c) ;
    localIndex blockIndex = blockCell.first;
    localIndex cellInBlock = blockCell.second;

    CellBlock const & cB = getCellBlock( blockIndex );
    HexFacetIndex faceToStore = f % Hex::nbFacets;

    assert(cellInBlock >= 0); // debug 

    // Maybe the oriented facets would be useful - if they are the ones recomputed 
    // later one by the FaceManager
    m_faceToNodes[curFace][0] = cB.getElementNode(cellInBlock, Hex::facetVertex[faceToStore][0]);
    m_faceToNodes[curFace][1] = cB.getElementNode(cellInBlock, Hex::facetVertex[faceToStore][1]);
    m_faceToNodes[curFace][2] = cB.getElementNode(cellInBlock, Hex::facetVertex[faceToStore][2]);
    m_faceToNodes[curFace][3] = cB.getElementNode(cellInBlock, Hex::facetVertex[faceToStore][3]);
  }

  return true;
}


bool HexMeshConnectivityBuilder::computeEdges()
{
  // Let's use a different strategy than from the faces 
  // Each edge is identified by a pair permitting to identify 
  // the 2 vertices: v0 * nbNodes + v1 - and the edge:  12 *idCell + e

  // 1 - Allocate 
  m_edges.resize( nbElements * Hex::nbEdges );
  
  // 2 - Get all edges
  localIndex cur = 0 ;
  for (int i = 0; i < nbCellBlocks(); ++i)
  {
    CellBlock const & block = getCellBlock(i);
    CellVertexIndices const & cells = block.getElemToNodes();

    for (int j = 0; j < cells.size(); ++j)
    {
      for( int e = 0; e < Hex::nbEdges; ++e)
      {
        LocalVertexIndex v0 = cells( j, Hex::edgeVertex[e][0]);
        LocalVertexIndex v1 = cells( j, Hex::edgeVertex[e][1]);
        if (v1 < v0 ) { std::swap(v0, v1); }
        LocalVertexIndex v = v0 * nbNodes + v1;
        LocalEdgeIndex id = j * Hex::nbEdges + e;

        m_edges[cur] = std::pair< LocalVertexIndex,LocalEdgeIndex > (v, id);
      }
    }
  }
  // 3 - Remove duplicates
  std::sort( m_edges.begin(), m_edges.end(), compareEdges );
  
  // This is no estimation if CellBlockManager manages a topological ball
  // Connected set of cells - boundary is a sphere (to check?) ie no holes
  // v - e + f - c = 1 
  SizeOfStuff guessNbEdges = nbNodes + nbFaces - nbElements;
  m_uniqueEdgeIds.reserve( guessNbEdges );

  // 4 - Count the edges and set the first one of duplicates
  int i = 0;
  while( i < m_edges.size() )
  {
    m_uniqueEdgeIds.push_back(i);
    int j = i+1;
    while( j < m_edges.size() && equalEdges(m_edges[i], m_edges[j]) )
    { 
      ++j; 
    }
    i = j;
  }
  // TODO Tests and asserts

  nbEdges = m_uniqueEdgeIds.size();

  return true;
}


bool HexMeshConnectivityBuilder::computeNodesToEdges()
{ 
  // 1 - Counting 
  std::vector<unsigned int> nbEdgesPerNode(nbNodes, 0);
  if (!isRegularMesh())
  {
    for( int i = 0; i < m_uniqueEdgeIds.size(); ++i)
    {
      EdgeInfo e = m_edges[ m_uniqueEdgeIds[i] ];
      LocalVertexIndex v0 = e.first / nbNodes;
      LocalVertexIndex v1 = e.first % nbNodes;
      nbEdgesPerNode[v0]++;
      nbEdgesPerNode[v1]++;
    }
  }
  else{
    nbEdgesPerNode.assign(nbNodes, 6);
  }
  localIndex valuesToReserve = std::accumulate(nbEdgesPerNode.begin(), nbEdgesPerNode.end(), 0);

  // 2 - Allocating 
  m_nodeToEdges.resize(0);
  m_nodeToEdges.reserve(nbNodes);
  m_nodeToEdges.reserveValues(valuesToReserve);
  for (localIndex n = 0; n < nbNodes; ++n)
  {
    m_nodeToEdges.appendSet( nbEdgesPerNode[n] );
  }
  
  // 3 - Filling
  for( int i = 0; i < m_uniqueEdgeIds.size(); ++i)
  {
      EdgeInfo e = m_edges[ m_uniqueEdgeIds[i] ];
      LocalVertexIndex v0 = e.first / nbNodes;
      LocalVertexIndex v1 = e.first % nbNodes;

    m_nodeToEdges.insertIntoSet(v0, i); 
    m_nodeToEdges.insertIntoSet(v1, i);
  }

  return true;
}

bool HexMeshConnectivityBuilder::computeNodesToFaces()
{
  if( m_faceToNodes.size()== 0) 
  {
    computeFacesToNodes();
  }
  // 1 - Counting
  std::vector<unsigned int> nbFacesPerNode(nbNodes, 0);
  if (!isRegularMesh())
  {
    for(int i = 0; i < m_faceToNodes.size(); ++i)
    {
      nbFacesPerNode[m_faceToNodes[i][0]]++;
      nbFacesPerNode[m_faceToNodes[i][1]]++;
      nbFacesPerNode[m_faceToNodes[i][2]]++;
      nbFacesPerNode[m_faceToNodes[i][3]]++;
    }
  }
  else
  {
    nbFacesPerNode.assign(nbNodes, 12);
  }
  localIndex valuesToReserve = std::accumulate(nbFacesPerNode.begin(), nbFacesPerNode.end(), 0);

  // 2 - Allocating 
  m_nodeToFaces.resize(0);
  m_nodeToFaces.reserve(nbNodes);
  m_nodeToFaces.reserveValues(valuesToReserve);
  for (localIndex n = 0; n < nbNodes; ++n)
  {
    m_nodeToFaces.appendSet( nbFacesPerNode[n] );
  }
  
  // 3 - Filling 
  for(int i = 0; i < m_faceToNodes.size(); ++i)
  {
    m_nodeToFaces.insertIntoSet( m_faceToNodes[i][0], i );
    m_nodeToFaces.insertIntoSet( m_faceToNodes[i][1], i );
    m_nodeToFaces.insertIntoSet( m_faceToNodes[i][2], i );
    m_nodeToFaces.insertIntoSet( m_faceToNodes[i][3], i );
  }

  return false;
}

// TODO We have a problem - where are the Element index valid ? 
bool HexMeshConnectivityBuilder::computeNodesToElements()
{
  // 1 -  Counting
  std::vector<unsigned int> nbElementsPerNode(nbNodes, 0);

  if (!isRegularMesh())
  {
    for (int i = 0; i < nbCellBlocks(); ++i)
    {
      CellBlock const & block = getCellBlock(i);
      CellVertexIndices const & cells = block.getElemToNodes();
    
      for (int j = 0; j < block.numElements(); ++j)
      {
        for (int v = 0; v < Hex::nbVertices; ++v)
        {
          nbElementsPerNode[ cells( j, v ) ]++;   
        }
      }
    }
  }
  nbElementsPerNode.assign(nbNodes, 8);
  localIndex nbValues = std::accumulate(nbElementsPerNode.begin(), nbElementsPerNode.end(), 0);

  //  2 - Allocating - No overallocation
  m_nodeToElements.reserve(nbNodes);
  m_nodeToElements.reserveValues(nbValues);
  for (localIndex v = 0; v < nbNodes; ++v)
  {
    m_nodeToElements.appendArray(0);
    m_nodeToElements.setCapacityOfArray(v, nbElementsPerNode[v]);
  }

  // 3 - Set the values
  for (int i = 0; i < nbCellBlocks(); ++i)
  {
    CellBlock const & block = getCellBlock(i);
    CellVertexIndices const & cells = block.getElemToNodes();

    for (int j = 0; j < block.numElements(); ++j)
    {
      for (int v = 0; v < Hex::nbVertices; ++v)
      {
        localIndex const nodeIndex = cells(j, v);
        m_nodeToElements.emplaceBack(nodeIndex, j);
      }
    }
  }
  return true;
}

// TODO We have a problem - where are the Element index valid ? 
bool HexMeshConnectivityBuilder::computeFacesToElements() 
{
  m_faceToElements.resize( nbFaces, 2);

  for( int curFace = 0; curFace < m_uniqueFaces.size(); ++curFace)
  {
    localIndex f = m_uniqueFaces[curFace];
    LocalFaceIndex neighbor = m_cellNeighbors[f];

    LocalCellIndex c = f/Hex::nbFacets;
    LocalCellIndex cNeighbor = -1;
    
    localIndex cellInBlock = getBlockCellFromManagerCell(c).second;
    localIndex neighborCellInMaybeNotSameBlock = -1;
    
    if (neighbor != -1) {
      cNeighbor = neighbor/Hex::nbFacets;
      neighborCellInMaybeNotSameBlock = getBlockCellFromManagerCell( cNeighbor ).second;
    }

    // TODO and here is ONE problem: mapping is toward Elements that may not be in the same CellBlock
    // Check what was done - Where is this used and what for?
    m_faceToElements( curFace, 0) = cellInBlock;
    m_faceToElements( curFace, 1) = neighborCellInMaybeNotSameBlock;
  }

  return false;
}

bool HexMeshConnectivityBuilder::computeFacesToEdges   () 
{
  // What is the point of this? What the hell is it useful for?
  // 1 - Allocate - No overallocation
  m_faceToEdges.resize(0);

  // LvArray allocation is an unclear business
  m_faceToEdges.reserve( nbFaces );
  m_faceToEdges.reserveValues( 4 * nbFaces );
  for (int i = 0; i < nbFaces; ++i)
  {
    m_faceToEdges.appendArray( Hex::nbEdgesPerFacet  ); // What is the default value ? 
  }

  // Pre-compute the mapping to have uniqueFaceIndex from allFaceIndex
  std::vector< localIndex > allFacesToUniqueFace( m_cellNeighbors.size() );
  for( int curFace = 0; curFace < m_uniqueFaces.size(); ++curFace)
  {
    localIndex f = m_uniqueFaces[curFace];
    allFacesToUniqueFace[f] = curFace;
    if( !m_isBoundaryFace[curFace] )
    {
      allFacesToUniqueFace[f+1] = curFace;
    }
  }  

  // Then what is the best way to get this useless map ?
  // No idea - let's try one - 
  // Iterate on the faces to find the edges - or iterate on the edges to find the faces ?

  for (int edgeIndex = 0; edgeIndex < m_uniqueEdgeIds.size(); ++edgeIndex)
  {
    int last = (edgeIndex+1 < m_uniqueEdgeIds.size()) ? m_uniqueEdgeIds[edgeIndex+1] : m_edges.size();

    //performance is bad but I'm lazy
    std::set<localIndex> faces;

    for( int i = m_uniqueEdgeIds[edgeIndex]; i < last; ++i)
    {
      LocalEdgeIndex id = m_edges[i].second;
      LocalCellIndex c = id / Hex::nbEdges;
      HexEdgeIndex e = id % Hex::nbEdges;

      HexFacetIndex f0 = Hex::edgeAdjacentFacet[e][0];
      HexFacetIndex f1 = Hex::edgeAdjacentFacet[e][1];

      // Now we need to find  the index of these faces 
      // The indices in the neighbor cell vector are known
      LocalFaceIndex face0 = Hex::nbFacets * c + f0;
      LocalFaceIndex face1 = Hex::nbFacets * c + f1;
      
      // Get the uniqueFaceId 
      localIndex uniqueFace0 = allFacesToUniqueFace[face0];
      localIndex uniqueFace1 = allFacesToUniqueFace[face1];

      faces.insert(uniqueFace0);
      faces.insert(uniqueFace1);
    }
    for( auto itr( faces.begin()); itr !=faces.end(); ++itr)
    {
      // Stupid but whatever -  this assumes arrays default values is ZERO  
      if     ( m_faceToEdges( *itr, 0) == 0 )  m_faceToEdges( *itr, 0 ) = edgeIndex ;
      else if( m_faceToEdges( *itr, 1) == 0 )  m_faceToEdges( *itr, 1 ) = edgeIndex ;
      else if( m_faceToEdges( *itr, 2) == 0 )  m_faceToEdges( *itr, 2 ) = edgeIndex ;
      else if( m_faceToEdges( *itr, 3) == 0 )  m_faceToEdges( *itr, 3 ) = edgeIndex ;
      else assert( false ); 
    }
  }

  // TODO Asserts to check that everything went well 
  // 4 edges exactly per facet  

  return false;
}

bool HexMeshConnectivityBuilder::computeEdgesToFaces   () 
{
  // What is the point of this? What the hell is it useful for?
  // We are going to do exactly the same loop than computeFacesToEdges

  // These two maps seem both as useless  - compute them at the same time ?
  // What is the point of this? What the hell is it useful for?

  // 1 - Allocate - Without counting  -- Should we take an overallocation shit
  // Probably for tet meshes one should count
  m_edgeToFaces.resize(0);

  // LvArray allocation is an unclear business
  m_edgeToFaces.reserve( nbEdges );
  m_edgeToFaces.reserveValues( 4 * nbEdges );
  for (int i = 0; i < nbEdges; ++i)
  {
    m_edgeToFaces.appendSet( 4 );
  }

  // Pre-compute the mapping to have uniqueFaceIndex from allFaceIndex
  std::vector< localIndex > allFacesToUniqueFace( m_cellNeighbors.size() );
  for( int curFace = 0; curFace < m_uniqueFaces.size(); ++curFace)
  {
    localIndex f = m_uniqueFaces[curFace];
    allFacesToUniqueFace[f] = curFace;
    if( !m_isBoundaryFace[curFace] )
    {
      allFacesToUniqueFace[f+1] = curFace;
    }
  }  

  // Then what is the best way to get this useless map ?
  // No idea - let's try one - 
  // Iterate on the faces to find the edges - or iterate on the edges to find the faces ?

  for (int edgeIndex = 0; edgeIndex < m_uniqueEdgeIds.size(); ++edgeIndex)
  {
    int last = (edgeIndex+1 < m_uniqueEdgeIds.size()) ? m_uniqueEdgeIds[edgeIndex+1] : m_edges.size();
    //performance is bad but I'm lazy
    std::set<localIndex> faces;

    for( int i = m_uniqueEdgeIds[edgeIndex]; i < last; ++i)
    {
      LocalEdgeIndex id = m_edges[i].second;
      LocalCellIndex c = id / Hex::nbEdges;
      HexEdgeIndex e = id % Hex::nbEdges;

      HexFacetIndex f0 = Hex::edgeAdjacentFacet[e][0];
      HexFacetIndex f1 = Hex::edgeAdjacentFacet[e][1];

      // Now we need to find  the index of these faces 
      // The indices in the neighbor cell vector are known
      LocalFaceIndex face0 = Hex::nbFacets * c + f0;
      LocalFaceIndex face1 = Hex::nbFacets * c + f1;
      
      // Get the uniqueFaceId 
      localIndex uniqueFace0 = allFacesToUniqueFace[face0];
      localIndex uniqueFace1 = allFacesToUniqueFace[face1];
    
      faces.insert(uniqueFace0);
      faces.insert(uniqueFace1);
    }
    for( auto itr( faces.begin()); itr !=faces.end(); ++itr)
    {
      m_edgeToFaces.insertIntoSet( edgeIndex, *itr );
    }
  }

  return false;
}

// We must have the CellNeighbors 
bool HexMeshConnectivityBuilder::computeElementsToFaces() 
{
  // 1 - Allocation is managed by the CellBlock 

  // 2 - Fill ElementToFace
  for( int curFace = 0; curFace < m_uniqueFaces.size(); ++curFace)
  {
    LocalFaceIndex id = m_uniqueFaces[curFace];

    LocalCellIndex c = id / Hex::nbFacets;
    HexFacetIndex face = id % Hex::nbFacets;
    auto blockCell = getBlockCellFromManagerCell(c) ;
    getCellBlock( blockCell.first ).setElementToFaces( blockCell.second, face, curFace );

    LocalFaceIndex idNeighbor = m_cellNeighbors[id];
    if( idNeighbor != -1)
    {
      LocalCellIndex cNeighbor = idNeighbor / Hex::nbFacets;
      HexFacetIndex faceNeighbor = idNeighbor % Hex::nbFacets;
      auto blockCellNeighbor = getBlockCellFromManagerCell(cNeighbor) ;
      getCellBlock( blockCellNeighbor.first ).setElementToFaces( blockCellNeighbor.second, faceNeighbor, curFace );
    }
  }
  return true;
}


bool HexMeshConnectivityBuilder::computeElementsToEdges() 
{
  for (int edgeIndex = 0; edgeIndex < m_uniqueEdgeIds.size(); ++edgeIndex)
  {
    int last = (edgeIndex+1 < m_uniqueEdgeIds.size()) ? m_uniqueEdgeIds[edgeIndex+1] : m_edges.size();
    for( int i = m_uniqueEdgeIds[edgeIndex]; i < last; ++i)
    {
      LocalEdgeIndex id = m_edges[i].second;
      LocalCellIndex c = id / Hex::nbEdges;
      HexEdgeIndex e = id % Hex::nbEdges;
      auto blockCell = getBlockCellFromManagerCell(c);
      getCellBlock(blockCell.first).setElementToEdges(blockCell.second, e, edgeIndex);
    }
  }
  return true;
}

/*********************************************************************************************************************/
/*********************************************************************************************************************/
/*********************************************************************************************************************/
/*********************************************************************************************************************/
/*********************************************************************************************************************/


HexCellBlockManager::HexCellBlockManager(string const &name, Group *const parent) : CellBlockManagerABC(name, parent),
                                                                                    m_nodesPositions(0, 3)
{
  this->registerGroup<Group>(viewKeyStruct::cellBlocks());
  m_theOneWhoDoesTheJob = new HexMeshConnectivityBuilder( *this );
}

HexCellBlockManager::~HexCellBlockManager()
{
  delete m_theOneWhoDoesTheJob;
  m_theOneWhoDoesTheJob = nullptr;
}

void HexCellBlockManager::resize(integer_array const &numElements,
                                 string_array const &regionNames)
{
  localIndex const numRegions = LvArray::integerConversion<localIndex>(regionNames.size());
  for (localIndex reg = 0; reg < numRegions; ++reg)
  {
    this->getCellBlock(regionNames[reg]).resize(numElements[reg]);
  }
}

Group *HexCellBlockManager::createChild(string const &GEOSX_UNUSED_PARAM(childKey), string const &GEOSX_UNUSED_PARAM(childName))
{
  return nullptr;
}

array2d<geosx::localIndex> HexCellBlockManager::getEdgeToNodes()
{
  return m_theOneWhoDoesTheJob->getEdgeToNodes();
}
ArrayOfSets<geosx::localIndex> HexCellBlockManager::getEdgeToFaces()
{
  return m_theOneWhoDoesTheJob->getEdgeToFaces();
}

ArrayOfArrays<localIndex> HexCellBlockManager::getFaceToNodes()
{
  return m_theOneWhoDoesTheJob->getFaceToNodes();
}
ArrayOfArrays<geosx::localIndex> HexCellBlockManager::getFaceToEdges()
{
  return m_theOneWhoDoesTheJob->getFaceToEdges();
}
array2d<localIndex> HexCellBlockManager::getFaceToElements()
{
  return m_theOneWhoDoesTheJob->getFaceToElements();
}
ArrayOfSets<localIndex> HexCellBlockManager::getNodeToEdges()
{
  return m_theOneWhoDoesTheJob->getNodeToEdges();
}
ArrayOfSets<localIndex> HexCellBlockManager::getNodeToFaces()
{
  return m_theOneWhoDoesTheJob->getNodeToFaces();
}
ArrayOfArrays<localIndex> HexCellBlockManager::getNodeToElements()
{
  return m_theOneWhoDoesTheJob->getNodeToElements();
}

void HexCellBlockManager::buildMaps()
{
  m_theOneWhoDoesTheJob->computeAllMaps();
}

void HexCellBlockManager::setNumNodes(localIndex numNodes)
{
  m_nodesPositions.resize(numNodes);
  m_nodeLocalToGlobal.resize(numNodes);
  m_nodeLocalToGlobal.setValues<serialPolicy>(-1);
}

} // namespace geosx
