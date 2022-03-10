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

#include <algorithm>
#include <numeric>

#pragma clang diagnostic ignored "-Wunused" // TODO Remove this - too much of a pain when developing new features
//#pragma clang diagnostic ignored "-Wcomment"
#pragma clang diagnostic ignored "-Wsign-compare"

namespace geosx
{

typedef localIndex  LocalVertexIndex;
typedef globalIndex GlobalVertexIndex;
typedef localIndex  LocalCellIndex;
typedef globalIndex GlobalCellIndex; 
typedef localIndex  LocalFaceIndex;  // A  face in a cell nbFaces * idCell + f
typedef localIndex  LocalEdgeIndex;  // An edge in a cell nbEdges * idCell + e
typedef localIndex  CellBlockIndex;
typedef localIndex  SizeOfStuff;


// Use explicit aliases for indices in a hexahedron
typedef unsigned int HexVertexIndex;      // a vertex in a hex 0 to 8
typedef unsigned int HexEdgeIndex;        // a edge in a hex 0 to 12
typedef unsigned int HexFacetIndex;       // a facet in a hex 0 to 6

static const int NO_ID = -1;

/*  Hexahedron template
*  WARNING - Hex vertex numbering in GEOSX differs from the one used 
*  by most mesh datastructurees
*
*   6----------7
*   |\         |\
*   | \        | \
*   |  \       |  \
*   |   4------+---5
*   |   |      |   |
*   2---+------3   |
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
static const unsigned int nbEdgesPerFacet = 4;

static constexpr HexVertexIndex facetVertex[6][4] = {
    {0, 1, 3, 2}, {4, 5, 7, 6}, {0, 1, 5, 4}, {1, 3, 7, 5}, {2, 3, 7, 6}, {0, 2, 6, 4}};
//         0             1           2            3             4               5          
static constexpr HexVertexIndex edgeVertex[12][2]{
    {0, 1}, {0, 2}, {0, 4}, {1, 3}, {1, 5}, {2, 3}, {2, 6}, {3, 7}, {4, 5}, {4, 6}, {5, 7}, {6, 7}};

static constexpr HexVertexIndex edgeAdjacentFacet[12][2]{
    {0, 2}, {0, 5}, {2, 5}, {0, 3}, {2, 3}, {0, 4}, {5, 4}, {3, 4}, {2, 1}, {1, 5}, {1, 3}, {1, 4}};

// This the ordering needed to compute consistent normals
// static constexpr HexVertexIndex orientedFacetVertex[6][4] = {
//    {0, 1, 3, 2}, {4, 6, 7, 5}, {0, 4, 5, 1}, {1, 5, 7, 3}, {2, 3, 7, 6}, {0, 2, 6, 4}};
};


// Structure for face computation
// TODO Optimize: replace (v0, v1, v2, f) by (v0*nbNodes + v1, v2, f)
struct FaceInfo{
  localIndex v[4];
};

struct {
  bool operator()( FaceInfo const & a, FaceInfo const & b) const {
    if(a.v[0] != b.v[0]) return a.v[0] < b.v[0];
    if(a.v[1] != b.v[1]) return a.v[1] < b.v[1];
    if(a.v[2] != b.v[2]) return a.v[2] < b.v[2];
    return false;
  }
} compareFaceInfo;

// Structure for edge computation and storage
//typedef std::pair< LocalVertexIndex, LocalEdgeIndex > EdgeInfo;
struct EdgeInfo
{
  EdgeInfo(LocalVertexIndex a, LocalEdgeIndex b){
   first = a;
   second = b;
  }
  EdgeInfo(){
   first = 0;
   second = 0;
  }

  LocalVertexIndex first;
  LocalEdgeIndex second;
};

struct {
  bool operator() ( EdgeInfo const & a, EdgeInfo const & b ) const
  {
    return a.first < b.first;
  }
} compareEdgeInfo;

struct {
  bool operator() ( EdgeInfo const & a, EdgeInfo const & b ) const
  {
    return a.first == b.first;
  }
} equalEdgeInfo;


void print(std::vector<EdgeInfo> const &in)
{
  for (int i = 0; i < in.size(); ++i)
  {

    std::cout << std::setw(5) << i
              << std::setw(5) << std::left << in[i].first 
              << std::setw(5) << std::left << in[i].second 
              << std::endl;
  }
  std::cout << std::endl
            << std::endl;
}

typedef array2d<localIndex, cells::NODE_MAP_PERMUTATION> CellVertexIndices;


// LvArrays are a mystery 
// TODO Switch to LvAraray once the code is tested  
typedef std::vector<std::vector<localIndex>> stdArrayOfArraysOfIndex;
typedef std::vector<localIndex> stdArrayOfIndex;


/***************************************************************************************/
/***************************************************************************************/

// DEBUGGING
void print(std::vector<localIndex> const &in)
{
  int count = 0 ;
  for (int i = 0; i < in.size(); ++i)
  {
    if(count == 10 )
    {
      std::cout << std::endl; 
      count = 0;
    }
    std::cout << std::setw(5) << std::left << in[i];
    count++;
  }
  std::cout << std::endl;
}
void print(std::vector<bool> const &in)
{
  int count = 0 ;
  for (int i = 0; i < in.size(); ++i)
  {
    if(count == 10 )
    {
      std::cout << std::endl; 
      count = 0;
    }
    std::cout << std::setw(5) << std::left << in[i];
    count++;
  }
  std::cout << std::endl    << std::endl;
}

void print( ArrayOfSets<localIndex> const & in )
{
  for( int i = 0; i < in.size(); ++i) {
    for(int j = 0; j < in.sizeOfSet(i);++j)
    {
      std::cout << std::setw(5) << std::left << in(i,j);
    }
    std::cout << std::endl;
  }
}

void print(stdArrayOfArraysOfIndex const &in)
{
  for (int i = 0; i < in.size(); ++i)
  {
    std::cout << std::setw(5) << std::left << i ;
    for (int j = 0; j < in[i].size(); ++j)
    {
      std::cout << std::setw(5) << std::right << in[i][j];
    }
    std::cout << std::endl;
  }
}

void print(array2d<localIndex> const &in)
{
  for (int i = 0; i < in.size(); ++i)
  {
      std::cout << std::setw(5) << std::left << in(i,0);
      std::cout << std::setw(5) << std::left << in(i,1);
  }
  std::cout << std::endl
            << std::endl;
}

/***************************************************************************************/
/***************************************************************************************/


/**
 * @class HexMeshConnectivityBuilder
 * @brief The HexMeshConnectivityBuilder to build the connectivity maps 
 * 
 * Initially designed for hexahedral meshes (unstructured)
 * TODO How do we reuse this for other type of cells 
 * Most of the code will be the same 
 * Is templating an option? Or derivation? 
 * 
 * TODO and here is ONE problem: all mapping is toward Elements are not safe 
 * since the elements may not be in the same CellBlock
 * Check what was done - Where is this used and what for?
 *
 * TODO Why are storage strategies different for the mappings ?
 * TODO Why multidimensional arrays? Isn't is more expensive? 
 * TODO Switch storage to LvArrays - What gains? 
 * TODO Implement specializationfor regular hex mesh
 *  
 */
class HexMeshConnectivityBuilder
{
public:
  HexMeshConnectivityBuilder() = default;
  ~HexMeshConnectivityBuilder() = default;

  bool initialize( HexCellBlockManager & cellBlockManager );

  stdArrayOfArraysOfIndex getFaceToNodes();
  stdArrayOfArraysOfIndex getEdgeToNodes();

  stdArrayOfArraysOfIndex getFaceToElements();

  stdArrayOfArraysOfIndex getNodeToEdges();
  stdArrayOfArraysOfIndex getNodeToFaces();

  stdArrayOfArraysOfIndex getNodeToElements();

  stdArrayOfArraysOfIndex getEdgeToFaces();

  stdArrayOfArraysOfIndex getFaceToEdges();

  // Why is this here? The information is on the CellBlocks
  bool computeElementsToFacesOfCellBlocks(); 
  bool computeElementsToEdgesOfCellBlocks();

  bool debuggingComputeAllMaps();

private:
  bool computeFaces();
  bool computeEdges(); 
  
  bool computeNodesToEdges( stdArrayOfArraysOfIndex & result ) const;
  bool computeNodesToFaces( stdArrayOfArraysOfIndex & result ) const ;
  bool computeNodesToElements( stdArrayOfArraysOfIndex & result ) const;

  bool computeEdgesToNodes( stdArrayOfArraysOfIndex& result ) const;
  bool computeEdgesToFaces( stdArrayOfArraysOfIndex & result ) const;

  bool computeFacesToNodes( stdArrayOfArraysOfIndex & result ) const;
  bool computeFacesToElements( stdArrayOfArraysOfIndex &  result ) const;
  bool computeFacesToEdges( stdArrayOfArraysOfIndex & result ) const;

  std::pair< CellBlockIndex, LocalCellIndex> 
  getBlockCellFromManagerCell( LocalCellIndex cellId ) const;

  unsigned int nbCellBlocks() const {
    return m_cellBlocks.size();
  } 
  CellBlock const & getCellBlock( unsigned int id ) const {
    return *(m_cellBlocks[id]);
  }
  CellBlock & getCellBlock( unsigned int id ){
    return *(m_cellBlocks[id]);
  }
  bool getEdges() {
    if( nbEdges == 0  ){
      return computeEdges();
    }
    else return true;
  }
  bool getFaces() {
    if( nbFaces == 0  ){
      return computeFaces();
    } 
    else return true;
  }
  // TODO Do we need it? It might prove useful to settle any issues 
  // arising with indices
  bool computeVertexGlobalToLocal( HexCellBlockManager const & cellBlockManager );

  void computeAllFacesToUniqueFace( std::vector< localIndex > & result ) const;

  void getFacesAroundEdge( int edgeIndex,
                           std::vector< localIndex > const & allFacesToUniqueFace, 
                           std::set <localIndex> & faces ) const;

private:
// Input 
SizeOfStuff nbNodes = 0;
SizeOfStuff nbElements = 0;   // in all CellBlocks

// Computed 
SizeOfStuff nbEdges = 0;
SizeOfStuff nbFaces = 0;

// Identification and number of the Cells hold by the cell blocks 
std::vector< CellBlock * > m_cellBlocks;              // nbBlocks
std::vector< LocalCellIndex > m_blockCellIndexOffset; // nbBlocks

// TODO Switch to LvArrays - For what gains? 
std::vector< LocalFaceIndex > m_allFacesToNeighbors;  // 6 * nbElements
std::vector< localIndex > m_uniqueFaces;              // nbFaces
std::vector< bool > m_isBoundaryFace;                 // nbFaces

// TODO The storage of Faces and Edges is dependant on the Cell Types  
// How do we manage the CellBlocks with different types of cells?
// Strategy 1: allocate the space for hexahedra and keep a lot of invalid stuff 
//             in these vectors - not the worst idea since these meshes should be hex-dominant
//             do we need to store the cell type?
// TODO For full tetrahedral meshes we need a dedicated implementation

std::vector< EdgeInfo > m_allEdges;                 // 12 * nbElements
std::vector< localIndex > m_uniqueEdges;            // nbEdges

std::unordered_map< GlobalVertexIndex, LocalVertexIndex > vertexGlobalToLocal; 
};


bool HexMeshConnectivityBuilder::initialize( HexCellBlockManager & cellBlockManager )
{
  nbNodes = cellBlockManager.numNodes();
  
  std::cout << "nbNodes " <<  nbNodes << std::endl;

  dataRepository::Group & group = cellBlockManager.getCellBlocks();
  localIndex nbBlocks = group.numSubGroups();
  
  std::cout << "nbBlocks " <<  nbBlocks << std::endl;

  m_cellBlocks.resize(nbBlocks);
  m_blockCellIndexOffset.resize(nbBlocks);

  localIndex prevSum = 0;
  for (int i = 0; i < nbBlocks; ++i)
  {
    // TODO - Check that we indeed have CellBlock instantiation 
    CellBlock & cur = group.getGroup< CellBlock >( i );
    std::cout << "nbBlockElements " <<  cur.numElements() << std::endl;

    m_cellBlocks[i] = &cur;
    m_blockCellIndexOffset[i] = prevSum + m_cellBlocks[i]->numElements();
  }

  nbElements =  m_blockCellIndexOffset[nbBlocks-1];
  std::cout << "nbElements " <<  nbElements << std::endl;

  return nbElements > 0;

  // Do we need it to solve potential inconsistencies / bugs?
  // computeVertexGlobalToLocal( cellBlockManager );
}

stdArrayOfArraysOfIndex HexMeshConnectivityBuilder::getFaceToNodes()
{
  getFaces();
  stdArrayOfArraysOfIndex result;
  computeFacesToNodes( result );
  return result;
}

stdArrayOfArraysOfIndex HexMeshConnectivityBuilder::getEdgeToNodes()
{
  getEdges();
  stdArrayOfArraysOfIndex result;
  computeEdgesToNodes( result );
  return result;
}

stdArrayOfArraysOfIndex HexMeshConnectivityBuilder::getNodeToEdges()
{
  getEdges();
  stdArrayOfArraysOfIndex result;
  computeNodesToEdges( result );
  return result;
}

stdArrayOfArraysOfIndex HexMeshConnectivityBuilder::getNodeToFaces()
{
  getFaces();
  stdArrayOfArraysOfIndex result;
  computeNodesToFaces( result );
  return result;
}

stdArrayOfArraysOfIndex HexMeshConnectivityBuilder::getNodeToElements()
{
  stdArrayOfArraysOfIndex result;
  computeNodesToElements( result );
  return result;
}

stdArrayOfArraysOfIndex HexMeshConnectivityBuilder::getFaceToElements()
{
  getFaces();
  stdArrayOfArraysOfIndex result;
  computeFacesToElements( result );
  return result;
}

stdArrayOfArraysOfIndex HexMeshConnectivityBuilder::getFaceToEdges()
{
  getFaces();
  getEdges();
  stdArrayOfArraysOfIndex result;
  computeFacesToEdges( result );
  return result;
}

stdArrayOfArraysOfIndex HexMeshConnectivityBuilder::getEdgeToFaces()
{
  getFaces();
  getEdges();
  stdArrayOfArraysOfIndex result;
  computeEdgesToFaces( result );
  return result;
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

bool HexMeshConnectivityBuilder::computeVertexGlobalToLocal(
   HexCellBlockManager const & cellBlockManager )
{
  vertexGlobalToLocal.reserve(nbNodes);
  array1d<GlobalVertexIndex> globalIds = cellBlockManager.getNodeLocalToGlobal();
  for (LocalVertexIndex i = 0; i < nbNodes; ++i)
  {
    GlobalVertexIndex id = globalIds[i];
    vertexGlobalToLocal[id] = i;
  }
  return true;
}


/**
 * @brief Compute all the faces of the cells of the cell blocks
 * Fills the face information 
 * TODO Be able to get out and return an error message if 3 cells have the same face
 * 
 * TODO CellType dependant 
 */
bool HexMeshConnectivityBuilder::computeFaces()
{
  // 1 - Allocate - 
  SizeOfStuff nbTotalFaces = nbElements * Hex::nbFacets;
  m_allFacesToNeighbors.resize( nbTotalFaces);
  
  // To collect and sort the facets 
  std::vector< FaceInfo > allFaces ( nbTotalFaces );

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
        // and that they are not local to the CellBlock 
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
  std::sort( allFaces.begin(), allFaces.end(), compareFaceInfo);
  

  // That is an overallocation, about twice as big as needed
  m_uniqueFaces.resize( allFaces.size() ); 
  m_isBoundaryFace.resize( allFaces.size(), false ); 

  // 4 - Counting + Set Cell Adjacencies
  curFace = 0;
  int i = 0;
  for( ; i+1 < nbTotalFaces; ++i )
  {
    FaceInfo f = allFaces[i];
    FaceInfo f1 = allFaces[i+1];

    // If two successive faces are the same - this is an interior facet
    if( f.v[0]==f1.v[0] && f.v[1]==f1.v[1] && f.v[2]==f1.v[2] )
    {
      m_allFacesToNeighbors[f.v[3]] = f1.v[3];
      m_allFacesToNeighbors[f1.v[3]] = f.v[3];
      m_uniqueFaces[curFace] = f.v[3];
      curFace++;
      ++i;
    }
    // If not this is a boundary face
    else
    {
      m_allFacesToNeighbors[f.v[3]] = NO_ID;
      m_uniqueFaces[curFace] = f.v[3];
      m_isBoundaryFace[curFace] = true;
      curFace++;
    }
  }
  // Last facet is a boundary facet 
  if (i < nbTotalFaces) 
  {
    FaceInfo f = allFaces[i];
    m_allFacesToNeighbors[f.v[3]] = NO_ID;
    m_uniqueFaces[curFace] = allFaces[i].v[3];
    m_isBoundaryFace[curFace] = true;
    curFace++;
  }
  // Set the number of faces 
  nbFaces = curFace; 
  // Remove the non-filled values 
  m_uniqueFaces.resize( nbFaces );
  m_isBoundaryFace.resize( nbFaces );

  return true;
}

/**
 * @brief Compute all the edges of the cells of the cell blocks
 * Fills the edge information 
 * 
 * TODO CellType dependant 
 */
bool HexMeshConnectivityBuilder::computeEdges()
{
  // 1 - Allocate 
  SizeOfStuff nbAllEdges = nbElements * Hex::nbEdges;
  m_allEdges.resize( nbAllEdges );
  
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
        assert (cur < nbAllEdges );

        LocalVertexIndex v0 = cells( j, Hex::edgeVertex[e][0]);
        LocalVertexIndex v1 = cells( j, Hex::edgeVertex[e][1]);

        if (v1 < v0 ) { std::swap(v0, v1); }
        LocalVertexIndex v = v0 * nbNodes + v1;
        LocalEdgeIndex id = j * Hex::nbEdges + e;

        m_allEdges[cur] = EdgeInfo(v, id);
        ++cur;
      }
    }
  }
  // 3 - Sort according to vertices 
  std::sort( m_allEdges.begin(), m_allEdges.end(), compareEdgeInfo );
  
  // 4 - Reserve space unique edges 
  // If a CellBlockManager manages a connected set of cells (any type of cell)
  // bounded by a sphere nbNodes - nbEdges + nbFaces - nbElements = 1 
  // For Hexes: 12 * nbCells = 4 * nbEdges + BoundaryEdges + Singularities - speculation 3.5 factor
  SizeOfStuff guessNbEdges = nbFaces != 0 ? nbNodes + nbFaces - nbElements : 3.5 * nbElements;
  m_uniqueEdges.reserve( guessNbEdges );

  // 4 - Get unique edges
  int i = 0;
  while( i < m_allEdges.size() )
  {
    m_uniqueEdges.push_back(i);
    int j = i+1;
    while( j < m_allEdges.size() && equalEdgeInfo(m_allEdges[i], m_allEdges[j]) )
    { 
      ++j; 
    }
    i = j;
  }
  // TODO Tests and asserts

  std::cout << "computed Edges " << m_uniqueEdges.size() << std::endl;

  nbEdges = m_uniqueEdges.size();
  return true;
}


bool HexMeshConnectivityBuilder::computeFacesToNodes( stdArrayOfArraysOfIndex & faceToNodes ) const
{
  // 1 - Allocate - No overallocation
  faceToNodes.resize(0);
  faceToNodes.resize( nbFaces, std::vector<localIndex>(4, 0) );

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
    faceToNodes[curFace][0] = cB.getElementNode(cellInBlock, Hex::facetVertex[faceToStore][0]);
    faceToNodes[curFace][1] = cB.getElementNode(cellInBlock, Hex::facetVertex[faceToStore][1]);
    faceToNodes[curFace][2] = cB.getElementNode(cellInBlock, Hex::facetVertex[faceToStore][2]);
    faceToNodes[curFace][3] = cB.getElementNode(cellInBlock, Hex::facetVertex[faceToStore][3]);
  }
  return true;
}


bool HexMeshConnectivityBuilder::computeEdgesToNodes ( stdArrayOfArraysOfIndex & edgeToNodes ) const
{
  edgeToNodes.resize(0);
  edgeToNodes.resize(nbEdges, std::vector<localIndex>(2, 0));

  for (int i = 0; i < m_uniqueEdges.size(); ++i)
  {
    EdgeInfo e = m_allEdges[m_uniqueEdges[i]];
    LocalVertexIndex v0 = e.first / nbNodes;
    LocalVertexIndex v1 = e.first % nbNodes;

    edgeToNodes[i][0] = v0;
    edgeToNodes[i][1] = v1;
  }
  return true;
}

bool HexMeshConnectivityBuilder::computeNodesToEdges( stdArrayOfArraysOfIndex & nodeToEdges ) const
{ 
  // 1 - Counting 
  std::vector<unsigned int> nbEdgesPerNode(nbNodes, 0);
  for( int i = 0; i < m_uniqueEdges.size(); ++i)
  {
    EdgeInfo e = m_allEdges[ m_uniqueEdges[i] ];
    LocalVertexIndex v0 = e.first / nbNodes;
    LocalVertexIndex v1 = e.first % nbNodes;
    nbEdgesPerNode[v0]++;
    nbEdgesPerNode[v1]++;
  }
  localIndex valuesToReserve = std::accumulate(nbEdgesPerNode.begin(), nbEdgesPerNode.end(), 0);

  // 2 - Allocating 
  nodeToEdges.resize(nbNodes);
  for (localIndex n = 0; n < nbNodes; ++n)
  {
    nodeToEdges[n].reserve( nbEdgesPerNode[n] );
  }
  
  // 3 - Filling
  for( int i = 0; i < m_uniqueEdges.size(); ++i)
  {
    EdgeInfo e = m_allEdges[ m_uniqueEdges[i] ];
    LocalVertexIndex v0 = e.first / nbNodes;
    LocalVertexIndex v1 = e.first % nbNodes;

    nodeToEdges[v0].push_back(i); 
    nodeToEdges[v1].push_back(i);
  }
  return true;
}


bool HexMeshConnectivityBuilder::computeNodesToFaces( stdArrayOfArraysOfIndex & nodeToFaces ) const
{
  stdArrayOfArraysOfIndex faceToNodes; 
  computeFacesToNodes( faceToNodes );

  // 1 - Counting
  std::vector<unsigned int> nbFacesPerNode(nbNodes, 0);
  for(int i = 0; i < faceToNodes.size(); ++i)
  {
    nbFacesPerNode[faceToNodes[i][0]]++;
    nbFacesPerNode[faceToNodes[i][1]]++;
    nbFacesPerNode[faceToNodes[i][2]]++;
    nbFacesPerNode[faceToNodes[i][3]]++;
  }

  localIndex valuesToReserve = std::accumulate(nbFacesPerNode.begin(), nbFacesPerNode.end(), 0);

  // 2 - Allocating 
  nodeToFaces.resize(0);
  nodeToFaces.resize(nbNodes);
  for (localIndex n = 0; n < nbNodes; ++n)
  {
    nodeToFaces[n].reserve( nbFacesPerNode[n] );
  }
  
  // 3 - Filling 
  for(int i = 0; i < faceToNodes.size(); ++i)
  {
    nodeToFaces[faceToNodes[i][0]].push_back(i);
    nodeToFaces[faceToNodes[i][1]].push_back(i);
    nodeToFaces[faceToNodes[i][2]].push_back(i);
    nodeToFaces[faceToNodes[i][3]].push_back(i);
  }

  return false;
}

bool HexMeshConnectivityBuilder::computeNodesToElements( stdArrayOfArraysOfIndex & nodeToElements ) const
{
  // 1 -  Counting
  std::vector<unsigned int> nbElementsPerNode(nbNodes, 0);

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
  localIndex nbValues = std::accumulate(nbElementsPerNode.begin(), nbElementsPerNode.end(), 0);

  //  2 - Allocating - No overallocation
  nodeToElements.resize(nbNodes);
  for (localIndex v = 0; v < nbNodes; ++v)
  {
    nodeToElements[v].reserve(nbElementsPerNode[v]);
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
        nodeToElements[nodeIndex].push_back(j);
      }
    }
  }
  return true;
}

// TODO We have a problem - where are the Element index valid ?
bool HexMeshConnectivityBuilder::computeFacesToElements( stdArrayOfArraysOfIndex & faceToElements ) const
{
  faceToElements.resize(0);
  faceToElements.resize( nbFaces, std::vector<localIndex>(2,0) );

  for( int curFace = 0; curFace < m_uniqueFaces.size(); ++curFace)
  {
    localIndex f = m_uniqueFaces[curFace];
    LocalFaceIndex neighbor = m_allFacesToNeighbors[f];

    LocalCellIndex c = f/Hex::nbFacets;
    LocalCellIndex cNeighbor = -1;
    
    localIndex cellInBlock = getBlockCellFromManagerCell(c).second;
    localIndex neighborCellInMaybeNotSameBlock = -1;
    
    if (neighbor != -1) {
      cNeighbor = neighbor/Hex::nbFacets;
      neighborCellInMaybeNotSameBlock = getBlockCellFromManagerCell( cNeighbor ).second;
    }

    faceToElements[curFace][0] = cellInBlock;
    faceToElements[curFace][1] = neighborCellInMaybeNotSameBlock;
  }
  return true;
}

void HexMeshConnectivityBuilder::computeAllFacesToUniqueFace(
   std::vector< localIndex > & allFacesToUniqueFace ) const
{
  allFacesToUniqueFace.resize( m_allFacesToNeighbors.size(), -1 );

  for( int curFace = 0; curFace < m_uniqueFaces.size(); ++curFace)
  {
    localIndex f = m_uniqueFaces[curFace];
    allFacesToUniqueFace[f] = curFace;
    if( !m_isBoundaryFace[curFace] )
    {
      localIndex f1 = m_allFacesToNeighbors[f];
      allFacesToUniqueFace[f1] = curFace;
    }
  }  
}

// Using a set means bad performances and annoying code 
// TODO Change set to vector
void HexMeshConnectivityBuilder::getFacesAroundEdge( 
  int edgeIndex, 
  std::vector< localIndex > const & allFacesToUniqueFace,
  std::set <localIndex> & faces ) const
{
  faces.clear();

  int first = m_uniqueEdges[edgeIndex];
  int last = m_allEdges.size();

  if( edgeIndex+1 < nbEdges )
  {
    last = m_uniqueEdges[edgeIndex+1];
  }

  // For each duplicata of this edge
  for( int i = first; i < last; ++i)
  {
    LocalEdgeIndex id = m_allEdges[i].second;
    // Get the cell
    LocalCellIndex c = id / Hex::nbEdges;
    // Get the hex facet incident to the edge
    HexEdgeIndex e = id % Hex::nbEdges;
    HexFacetIndex f0 = Hex::edgeAdjacentFacet[e][0];
    HexFacetIndex f1 = Hex::edgeAdjacentFacet[e][1];
    // Get the face index in allFaces 
    LocalFaceIndex face0 = Hex::nbFacets * c + f0;
    LocalFaceIndex face1 = Hex::nbFacets * c + f1;
    // Get the uniqueFaceId 
    localIndex uniqueFace0 = allFacesToUniqueFace[face0];
    localIndex uniqueFace1 = allFacesToUniqueFace[face1];

    faces.insert(uniqueFace0);
    faces.insert(uniqueFace1);
  }
}

// What is the point of this? What the hell is it useful for?
bool HexMeshConnectivityBuilder::computeFacesToEdges( stdArrayOfArraysOfIndex & faceToEdges) const 
{
  // 1 - Allocate
  faceToEdges.resize(nbFaces, std::vector<localIndex>( Hex::nbEdgesPerFacet, -1) );

  std::vector< localIndex > allFacesToUniqueFace;
  computeAllFacesToUniqueFace( allFacesToUniqueFace );

  for (int edgeIndex = 0; edgeIndex < m_uniqueEdges.size(); ++edgeIndex)
  {
    // Performance is not the best
    std::set<localIndex> faces;
    getFacesAroundEdge( edgeIndex, allFacesToUniqueFace, faces );

    for( auto itr( faces.begin()); itr !=faces.end(); ++itr)
    { 
      if     ( faceToEdges[*itr][0] == -1 )  faceToEdges[*itr][0] = edgeIndex ;
      else if( faceToEdges[*itr][1] == -1 )  faceToEdges[*itr][1] = edgeIndex ;
      else if( faceToEdges[*itr][2] == -1 )  faceToEdges[*itr][2] = edgeIndex ;
      else if( faceToEdges[*itr][3] == -1 )  faceToEdges[*itr][3] = edgeIndex ;
      else assert( false ); 
    }
  }
  // TODO Asserts to check that everything went well 
  // 4 edges exactly per facet  

  return true;
}

bool HexMeshConnectivityBuilder::computeEdgesToFaces(stdArrayOfArraysOfIndex & edgeToFaces ) const
{
  // 1 - Allocate - Without counting  -- Should we overallocate a bit ?
  // Probably for tet meshes one should count
  edgeToFaces.resize(0);
  edgeToFaces.resize( nbEdges );
  for (int i = 0; i < nbEdges; ++i)
  {
    edgeToFaces[i].reserve(4); 
  }

  std::vector< localIndex > allFacesToUniqueFace;
  computeAllFacesToUniqueFace( allFacesToUniqueFace );

  for (int edgeIndex = 0; edgeIndex < m_uniqueEdges.size(); ++edgeIndex)
  {
    std::set<localIndex> faces;
    getFacesAroundEdge( edgeIndex,  allFacesToUniqueFace, faces );
    
    for( auto itr( faces.begin()); itr !=faces.end(); ++itr)
    {
      edgeToFaces[edgeIndex].push_back(*itr);
    }
  }
  return true;
}

// Allocation is managed by the CellBlock 
bool HexMeshConnectivityBuilder::computeElementsToFacesOfCellBlocks() 
{
  getFaces();

  // Debugging - I cannot check LvArrays contents - I do not understand it
  std::vector < std::vector < std::vector <localIndex> > > store;
  store.resize( nbCellBlocks() );
  for (int i = 0; i < nbCellBlocks(); ++i)
  {
    store[i].resize(
      getCellBlock(i).getElemToNodes().size(0), std::vector<localIndex>( Hex::nbFacets, -1 )
    );
  }

  for( int curFace = 0; curFace < m_uniqueFaces.size(); ++curFace)
  {
    LocalFaceIndex id = m_uniqueFaces[curFace];

    LocalCellIndex c = id / Hex::nbFacets;
    HexFacetIndex face = id % Hex::nbFacets;
    auto blockCell = getBlockCellFromManagerCell(c);
    assert( blockCell.second == c); // DEBUG TO REMOVE

    store[blockCell.first][blockCell.second][face] = curFace;
    // getCellBlock( blockCell.first ).setElementToFaces( blockCell.second, face, curFace );

    LocalFaceIndex idNeighbor = m_allFacesToNeighbors[id];
    if( idNeighbor != -1)
    {
      LocalCellIndex cNeighbor = idNeighbor / Hex::nbFacets;
      HexFacetIndex faceNeighbor = idNeighbor % Hex::nbFacets;
      auto blockCellNeighbor = getBlockCellFromManagerCell(cNeighbor) ;
     // getCellBlock( blockCellNeighbor.first ).setElementToFaces( blockCellNeighbor.second, faceNeighbor, curFace );
     store[blockCellNeighbor.first][blockCellNeighbor.second][faceNeighbor] = curFace;
    }
  }

  // DEBUG
  print( store[0]);

  return true;
}


bool HexMeshConnectivityBuilder::computeElementsToEdgesOfCellBlocks() 
{
  getEdges();

  // Debugging - I cannot check LvArrays contents - I do not understand it
  std::vector < std::vector < std::vector <localIndex> > > store;
  store.resize( nbCellBlocks() );
  for (int i = 0; i < nbCellBlocks(); ++i)
  {
    store[i].resize(
      getCellBlock(i).getElemToNodes().size(0), std::vector<localIndex>( Hex::nbEdges, -1 )
    );
  }

  for (int edgeIndex = 0; edgeIndex < m_uniqueEdges.size(); ++edgeIndex)
  {
    int last = (edgeIndex+1 < m_uniqueEdges.size()) ? m_uniqueEdges[edgeIndex+1] : m_allEdges.size();
    for( int i = m_uniqueEdges[edgeIndex]; i < last; ++i)
    {
      LocalEdgeIndex id = m_allEdges[i].second;
      LocalCellIndex c = id / Hex::nbEdges;
      HexEdgeIndex e = id % Hex::nbEdges;
      auto blockCell = getBlockCellFromManagerCell(c);
      //getCellBlock(blockCell.first).setElementToEdges(blockCell.second, e, edgeIndex);
      store[blockCell.first][blockCell.second][e] = edgeIndex;
    }
  }

  // DEBUG
  print( store[0]);
  return true;
}

bool HexMeshConnectivityBuilder::debuggingComputeAllMaps()
{
  std::cout << std::endl
            << std::endl;

  std::cout << " Faces " << std::endl
            << std::endl;

  computeFaces();
  std::cout << " NB Faces : " << nbFaces << std::endl
            << std::endl;

  std::cout << " Face Info " << std::endl
            << std::endl;

  print(m_allFacesToNeighbors);
  std::cout << std::endl
            << std::endl;

  print(m_uniqueFaces);
  std::cout << std::endl
            << std::endl;

  print(m_isBoundaryFace);
  std::cout << std::endl
            << std::endl;

  std::cout << " Faces To Nodes " << std::endl
            << std::endl;
  stdArrayOfArraysOfIndex toto5 = getFaceToNodes();

  print(toto5);

  std::cout << " Edges " << std::endl
            << std::endl;

  computeEdges();

  std::cout << " NB Edges : " << nbEdges << std::endl
            << std::endl;

  std::cout << " Edge Info " << std::endl
            << std::endl;

  print(m_allEdges);
  print(m_uniqueEdges);

  std::cout << " Edges To Nodes " << std::endl
            << std::endl;
  stdArrayOfArraysOfIndex toto3 = getEdgeToNodes();
  print(toto3);

  std::cout << "Face To Element " << std::endl
            << std::endl;
  stdArrayOfArraysOfIndex toto6 = getFaceToElements();
  print(toto6);

  std::cout << "Node To Element " << std::endl
            << std::endl;

  stdArrayOfArraysOfIndex toto2 = getNodeToElements();
  print(toto2);

  std::cout << "Face To Edges " << std::endl
            << std::endl;

  stdArrayOfArraysOfIndex toto7 = getFaceToEdges();
  print(toto7);

  std::cout << "Node To Edges " << std::endl
            << std::endl;

  stdArrayOfArraysOfIndex toto = getNodeToEdges();
  print(toto);

  std::cout << "Node To Faces " << std::endl
            << std::endl;
  stdArrayOfArraysOfIndex toto1 = getNodeToFaces();
  print(toto1);

  std::cout << "Edges To Faces " << std::endl
            << std::endl;

  stdArrayOfArraysOfIndex toto4 = getEdgeToFaces();
  print(toto4);

  std::cout << "Cell to Edges " << std::endl
            << std::endl;


  computeElementsToEdgesOfCellBlocks();


  computeElementsToFacesOfCellBlocks();

  return true;
}


/*********************************************************************************************************************/
/*********************************************************************************************************************/

HexCellBlockManager::HexCellBlockManager(string const &name, 
                                         dataRepository::Group *const parent) 
  :CellBlockManagerABC(name, parent),
   m_nodesPositions(0, 3)
{
  this->registerGroup<dataRepository::Group>(viewKeyStruct::cellBlocks());
  m_theOneWhoDoesTheJob = new HexMeshConnectivityBuilder();
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

dataRepository::Group *HexCellBlockManager::createChild(
  string const &GEOSX_UNUSED_PARAM(childKey), 
  string const &GEOSX_UNUSED_PARAM(childName))
{
  return nullptr;
}

array2d<localIndex> HexCellBlockManager::getEdgeToNodes()
{
   stdArrayOfArraysOfIndex result = m_theOneWhoDoesTheJob->getEdgeToNodes();
   return array2d<localIndex> ();
}
ArrayOfSets<localIndex> HexCellBlockManager::getEdgeToFaces()
{
  stdArrayOfArraysOfIndex result = m_theOneWhoDoesTheJob->getEdgeToFaces();
  return ArrayOfSets<localIndex>();
}
ArrayOfArrays<localIndex> HexCellBlockManager::getFaceToNodes()
{
  stdArrayOfArraysOfIndex result = m_theOneWhoDoesTheJob->getFaceToNodes();
  return ArrayOfArrays<localIndex>(); 
}
ArrayOfArrays<localIndex> HexCellBlockManager::getFaceToEdges()
{
  stdArrayOfArraysOfIndex result = m_theOneWhoDoesTheJob->getFaceToEdges();
  return ArrayOfArrays<localIndex>(); 
}
array2d<localIndex> HexCellBlockManager::getFaceToElements()
{
  stdArrayOfArraysOfIndex result =  m_theOneWhoDoesTheJob->getFaceToElements();
  return array2d<localIndex> ();
}
ArrayOfSets<localIndex> HexCellBlockManager::getNodeToEdges()
{
   stdArrayOfArraysOfIndex result = m_theOneWhoDoesTheJob->getNodeToEdges();
  return ArrayOfSets<localIndex>();
}
ArrayOfSets<localIndex> HexCellBlockManager::getNodeToFaces()
{
   stdArrayOfArraysOfIndex result = m_theOneWhoDoesTheJob->getNodeToFaces();
  return ArrayOfSets<localIndex>();
}
ArrayOfArrays<localIndex> HexCellBlockManager::getNodeToElements()
{
  stdArrayOfArraysOfIndex result = m_theOneWhoDoesTheJob->getNodeToElements();
  return ArrayOfArrays<localIndex>();
}

// Lazy computation - This should not do anything 
// Who should be calling this mapping computations ? 
void HexCellBlockManager::buildMaps()
{
  m_theOneWhoDoesTheJob->initialize( *this ); 
  m_theOneWhoDoesTheJob->debuggingComputeAllMaps(); 

  std::cout << " All maps are computed " << std::endl;
  return ;
}

void HexCellBlockManager::setNumNodes(localIndex numNodes)
{
  std::cout << "coucou c'est moi " << std::endl;
  m_numNodes = numNodes;
  m_nodesPositions.resize(numNodes);
  // TODO why is this allocated here and filled somewhere else?
  m_nodeLocalToGlobal.resize(numNodes);
  m_nodeLocalToGlobal.setValues<serialPolicy>(-1);
}

} // namespace geosx
