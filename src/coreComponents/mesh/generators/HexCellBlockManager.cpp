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

//#pragma clang diagnostic ignored "-Wunused" // TODO Remove this - too much of a pain when developing new features
//#pragma clang diagnostic ignored "-Wcomment"
#pragma clang diagnostic ignored "-Wsign-compare"


// Hide the code from the outside of this file
// To hide the full implementation there a nasty reinterpret_cast from void* 
// could be used
namespace 
{
  using namespace geosx;

typedef localIndex  LocalVertexIndex;
typedef globalIndex GlobalVertexIndex;
typedef localIndex  LocalCellIndex;
typedef globalIndex GlobalCellIndex; 
typedef localIndex  LocalFaceIndex;  // A  face in a cell nbFaces * idCell + f
typedef localIndex  LocalEdgeIndex;  // An edge in a cell nbEdges * idCell + e
typedef localIndex  CellBlockIndex;
typedef localIndex  SizeOfStuff;


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
static const unsigned int nbEdgesPerFacet = 4;

static constexpr HexVertexIndex facetVertex[6][4] = {
    {0, 1, 2, 3}, {4, 5, 6, 7}, {0, 1, 5, 4}, {1, 2, 6, 5}, {3, 2, 6, 7}, {0, 3, 7, 4}};

static constexpr HexVertexIndex edgeVertex[12][2]{
    {0, 1}, {0, 3}, {0, 4}, {1, 2}, {1, 5}, {2, 3}, {2, 6}, {3, 7}, {4, 5}, {4, 7}, {5, 6}, {6, 7}};

static constexpr HexVertexIndex edgeAdjacentFacet[12][2]{
    {0, 2}, {0, 5}, {2, 5}, {0, 3}, {2, 3}, {0, 4}, {3, 4}, {5, 4}, {2, 1}, {1, 5}, {1, 3}, {1, 4}};

// This the ordering needed to compute consistent normals
// static constexpr HexVertexIndex orientedFacetVertex[6][4] = {
//    {0, 1, 2, 3}, {4, 7, 6, 5}, {0, 4, 5, 1}, {1, 5, 6, 2}, {3, 2, 6, 7}, {0, 3, 7, 4}};
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
typedef std::pair< LocalVertexIndex, LocalEdgeIndex > EdgeInfo;

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

} // anonymous namespace

namespace geosx
{

typedef array2d<localIndex, cells::NODE_MAP_PERMUTATION> CellVertexIndices;

/* The class to build the connectivity maps 
*
* TODO and here is ONE problem: all mapping is toward Elements are not safe 
* since the elements may not be in the same CellBlock
* Check what was done - Where is this used and what for?
*
* TODO Why are storage strategies different for the mappings ?
* TODO Why multidimensional arrays? Isn't is more expensive? 
*/
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

  // TODO Implement specializationfor regular hex mesh
  bool isRegularMesh() const { return false; }

  // Why is this here? The information is on the CellBlocks
  bool computeElementsToFacesOfCellBlocks(); 
  bool computeElementsToEdgesOfCellBlocks();

private:
  bool computeFaces();
  bool computeEdges(); 
  
  bool computeNodesToEdges   ( ArrayOfSets<localIndex> & result ) const;
  bool computeNodesToFaces   ( ArrayOfSets<localIndex> & result ) const ;
  bool computeNodesToElements( ArrayOfArrays<localIndex> & result ) const;

  bool computeEdgesToNodes ( array2d<geosx::localIndex> & result ) const;
  bool computeEdgesToFaces ( ArrayOfSets<localIndex> & result ) const;
  bool computeEdgesToElements() { return false; } // Not implemented

  bool computeFacesToNodes   ( ArrayOfArrays<localIndex> & result ) const;
  bool computeFacesToElements( array2d<localIndex> &  result ) const;
  bool computeFacesToEdges   ( ArrayOfArrays<localIndex> & result ) const;

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
    if( nbEdges != 0  ){
      return computeEdges();
    }
    else return true;
  }
  bool getFaces() {
    if( nbFaces != 0  ){
      return computeFaces();
    } 
    else return true;
  }
  // To remove? 
  bool computeVertexGlobalToLocal( HexCellBlockManager const & cellBlockManager );

  void computeAllFacesToUniqueFace( std::vector< localIndex > & result ) const;

  void getFacesAroundEdge( int edgeIndex,
                           std::vector< localIndex > const & allFacesToUniqueFace, 
                           std::set <localIndex> & faces ) const;

private:
// Input 
SizeOfStuff nbNodes = 0;
SizeOfStuff nbElements = 0;

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

std::vector< EdgeInfo > m_allEdges;                   // 12 * nbElements
std::vector< localIndex > m_uniqueEdgeIds;            // nbEdges

// Do we need it ? Or are all indices local and we are good
std::unordered_map< GlobalVertexIndex, LocalVertexIndex > vertexGlobalToLocal; 

};


HexMeshConnectivityBuilder::HexMeshConnectivityBuilder( HexCellBlockManager & cellBlockManager )
{
  nbNodes = cellBlockManager.numNodes();

  dataRepository::Group & group =  cellBlockManager.getCellBlocks();
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

  // Do we need it ? Or are all indices local and we are good
  // computeVertexGlobalToLocal( cellBlockManager );
}

ArrayOfArrays<localIndex> HexMeshConnectivityBuilder::getFaceToNodes()
{
  getFaces();
  ArrayOfArrays<localIndex> result; 
  computeFacesToNodes( result );
  return result;
}

array2d<geosx::localIndex> HexMeshConnectivityBuilder::getEdgeToNodes()
{
  getEdges();
  array2d<localIndex> result;
  computeEdgesToNodes( result );
  return result;
}

ArrayOfSets<localIndex> HexMeshConnectivityBuilder::getNodeToEdges()
{
  getEdges();
  ArrayOfSets<localIndex> result;
  computeNodesToEdges( result );
  return result;
}

ArrayOfSets<localIndex> HexMeshConnectivityBuilder::getNodeToFaces()
{
  getFaces();
  ArrayOfSets<localIndex> result;
  computeNodesToFaces( result );
  return result;
}

ArrayOfArrays<localIndex> HexMeshConnectivityBuilder::getNodeToElements()
{
  ArrayOfArrays<localIndex> result;
  computeNodesToElements( result );
  return result;
}

array2d<localIndex> HexMeshConnectivityBuilder::getFaceToElements()
{
  getFaces();
  array2d<localIndex> result;
  computeFacesToElements( result );
  return result;
}

ArrayOfArrays<geosx::localIndex> HexMeshConnectivityBuilder::getFaceToEdges()
{
  getFaces();
  getEdges();
  ArrayOfArrays<localIndex> result;
  computeFacesToEdges( result );
  return result;
}
ArrayOfSets<geosx::localIndex> HexMeshConnectivityBuilder::getEdgeToFaces()
{
  getFaces();
  getEdges();
  ArrayOfSets<geosx::localIndex> result;
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

/* Compute and store the cell neighors for each face
 */
// TODO Be able to get out and return an error message if 3 cells have the same facet
// Flagging boundary facet if the MPI rank should be possible at this point
bool HexMeshConnectivityBuilder::computeFaces()
{
  // 1 - Allocate
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

  // That is an overallocation - about twice as big as needed
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
      m_allFacesToNeighbors[f.v[3]] = -1; //Replace by NO_ID - maxvalue of  something
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


bool HexMeshConnectivityBuilder::computeEdges()
{
  // 1 - Allocate 
  m_allEdges.resize( nbElements * Hex::nbEdges );
  
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

        m_allEdges[cur] = std::pair< LocalVertexIndex,LocalEdgeIndex > (v, id);
        ++cur;
      }
    }
  }
  // 3 - Sort according to vertices 
  std::sort( m_allEdges.begin(), m_allEdges.end(), compareEdgeInfo );
  
  // 4 - Reserve space uniqueedges 
  // If a CellBlockManager manages a connected set of cells (any type of cell)
  // bounded by a sphere nbNodes - nbEdges + nbFaces - nbElements = 1 
  // For Hexes: 12 * nbCells = 4 * nbEdges + BoundaryEdges + Singularities - speculation 3.5 factor
  SizeOfStuff guessNbEdges = nbFaces != 0 ? nbNodes + nbFaces - nbElements : 3.5 * nbElements;
  m_uniqueEdgeIds.reserve( guessNbEdges );

  // 4 - Get unique edges
  int i = 0;
  while( i < m_allEdges.size() )
  {
    m_uniqueEdgeIds.push_back(i);
    int j = i+1;
    while( j < m_allEdges.size() && equalEdgeInfo(m_allEdges[i], m_allEdges[j]) )
    { 
      ++j; 
    }
    i = j;
  }
  // TODO Tests and asserts

  nbEdges = m_uniqueEdgeIds.size();

  return true;
}

bool HexMeshConnectivityBuilder::computeFacesToNodes( ArrayOfArrays<localIndex> & faceToNodes ) const
{
  assert (nbFaces > 0 );
  // 1 - Allocate - No overallocation
  faceToNodes.resize(0);

  
  faceToNodes.reserve( nbFaces );
  faceToNodes.reserveValues( 4 * nbFaces );
  for (int i = 0; i < nbFaces; ++i)
  {
    faceToNodes.appendArray( 4 );
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
    faceToNodes[curFace][0] = cB.getElementNode(cellInBlock, Hex::facetVertex[faceToStore][0]);
    faceToNodes[curFace][1] = cB.getElementNode(cellInBlock, Hex::facetVertex[faceToStore][1]);
    faceToNodes[curFace][2] = cB.getElementNode(cellInBlock, Hex::facetVertex[faceToStore][2]);
    faceToNodes[curFace][3] = cB.getElementNode(cellInBlock, Hex::facetVertex[faceToStore][3]);
  }

  return true;
}


bool HexMeshConnectivityBuilder::computeEdgesToNodes ( array2d<geosx::localIndex> & edgeToNodes ) const
{
  assert(nbEdges > 0);

  edgeToNodes.resize(nbEdges, 2);

  for (int i = 0; i < m_uniqueEdgeIds.size(); ++i)
  {
    EdgeInfo e = m_allEdges[m_uniqueEdgeIds[i]];
    LocalVertexIndex v0 = e.first / nbNodes;
    LocalVertexIndex v1 = e.first % nbNodes;

    edgeToNodes(i, 0) = v0;
    edgeToNodes(i, 1) = v1;
  }
  return true;
}

bool HexMeshConnectivityBuilder::computeNodesToEdges( ArrayOfSets<localIndex> & nodeToEdges ) const
{ 
  // 1 - Counting 
  std::vector<unsigned int> nbEdgesPerNode(nbNodes, 0);
  if (!isRegularMesh())
  {
    for( int i = 0; i < m_uniqueEdgeIds.size(); ++i)
    {
      EdgeInfo e = m_allEdges[ m_uniqueEdgeIds[i] ];
      LocalVertexIndex v0 = e.first / nbNodes;
      LocalVertexIndex v1 = e.first % nbNodes;
      nbEdgesPerNode[v0]++;
      nbEdgesPerNode[v1]++;
    }
  }
  else 
  {
    nbEdgesPerNode.assign(nbNodes, 6); // because this is a Hex mesh
  }
  localIndex valuesToReserve = std::accumulate(nbEdgesPerNode.begin(), nbEdgesPerNode.end(), 0);

  // 2 - Allocating 
  nodeToEdges.resize(0);
  nodeToEdges.reserve(nbNodes);
  nodeToEdges.reserveValues(valuesToReserve);
  for (localIndex n = 0; n < nbNodes; ++n)
  {
    nodeToEdges.appendSet( nbEdgesPerNode[n] );
  }
  
  // 3 - Filling
  for( int i = 0; i < m_uniqueEdgeIds.size(); ++i)
  {
    EdgeInfo e = m_allEdges[ m_uniqueEdgeIds[i] ];
    LocalVertexIndex v0 = e.first / nbNodes;
    LocalVertexIndex v1 = e.first % nbNodes;

    nodeToEdges.insertIntoSet(v0, i); 
    nodeToEdges.insertIntoSet(v1, i);
  }

  return true;
}

bool HexMeshConnectivityBuilder::computeNodesToFaces( ArrayOfSets<localIndex> & nodeToFaces ) const
{
  assert( nbFaces > 0);

  ArrayOfArrays<localIndex> faceToNodes; 
  computeFacesToNodes( faceToNodes );

  // 1 - Counting
  std::vector<unsigned int> nbFacesPerNode(nbNodes, 0);
  if (!isRegularMesh())
  {
    for(int i = 0; i < faceToNodes.size(); ++i)
    {
      nbFacesPerNode[faceToNodes[i][0]]++;
      nbFacesPerNode[faceToNodes[i][1]]++;
      nbFacesPerNode[faceToNodes[i][2]]++;
      nbFacesPerNode[faceToNodes[i][3]]++;
    }
  }
  else
  {
    nbFacesPerNode.assign(nbNodes, 12);
  }
  localIndex valuesToReserve = std::accumulate(nbFacesPerNode.begin(), nbFacesPerNode.end(), 0);

  // 2 - Allocating 
  nodeToFaces.resize(0);
  nodeToFaces.reserve(nbNodes);
  nodeToFaces.reserveValues(valuesToReserve);
  for (localIndex n = 0; n < nbNodes; ++n)
  {
    nodeToFaces.appendSet( nbFacesPerNode[n] );
  }
  
  // 3 - Filling 
  for(int i = 0; i < faceToNodes.size(); ++i)
  {
    nodeToFaces.insertIntoSet( faceToNodes[i][0], i );
    nodeToFaces.insertIntoSet( faceToNodes[i][1], i );
    nodeToFaces.insertIntoSet( faceToNodes[i][2], i );
    nodeToFaces.insertIntoSet( faceToNodes[i][3], i );
  }

  return false;
}

bool HexMeshConnectivityBuilder::computeNodesToElements( ArrayOfArrays<localIndex> & nodeToElements ) const
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
  else 
  {
    nbElementsPerNode.assign(nbNodes, 8);
  }
  localIndex nbValues = std::accumulate(nbElementsPerNode.begin(), nbElementsPerNode.end(), 0);

  //  2 - Allocating - No overallocation
  nodeToElements.reserve(nbNodes);
  nodeToElements.reserveValues(nbValues);
  for (localIndex v = 0; v < nbNodes; ++v)
  {
    nodeToElements.appendArray(0);
    nodeToElements.setCapacityOfArray(v, nbElementsPerNode[v]);
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
        nodeToElements.emplaceBack(nodeIndex, j);
      }
    }
  }
  return true;
}

// TODO We have a problem - where are the Element index valid ? 
bool HexMeshConnectivityBuilder::computeFacesToElements( array2d<localIndex> & faceToElements ) const
{
  faceToElements.resize( nbFaces, 2);

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

    faceToElements( curFace, 0) = cellInBlock;
    faceToElements( curFace, 1) = neighborCellInMaybeNotSameBlock;
  }
  return true;
}

void HexMeshConnectivityBuilder::computeAllFacesToUniqueFace(
   std::vector< localIndex > & allFacesToUniqueFace ) const
{
  allFacesToUniqueFace.resize( m_allFacesToNeighbors.size() );
  for( int curFace = 0; curFace < m_uniqueFaces.size(); ++curFace)
  {
    localIndex f = m_uniqueFaces[curFace];
    allFacesToUniqueFace[f] = curFace;
    if( !m_isBoundaryFace[curFace] )
    {
      allFacesToUniqueFace[f+1] = curFace;
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

  int first = m_uniqueEdgeIds[edgeIndex];
  int last = m_allEdges.size();
  if( edgeIndex+1 < nbEdges )
  {
    last = m_uniqueEdgeIds[edgeIndex+1];
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
bool HexMeshConnectivityBuilder::computeFacesToEdges( ArrayOfArrays<localIndex> & faceToEdges) const 
{
  // 1 - Allocate
  faceToEdges.resize(0);
  faceToEdges.reserve( nbFaces );
  faceToEdges.reserveValues( Hex::nbEdgesPerFacet * nbFaces );
  for (int i = 0; i < nbFaces; ++i)
  {
    faceToEdges.appendArray( Hex::nbEdgesPerFacet ); // What is the default value? Can we set it?
  }

  std::vector< localIndex > allFacesToUniqueFace;
  computeAllFacesToUniqueFace( allFacesToUniqueFace );
  
  for (int edgeIndex = 0; edgeIndex < m_uniqueEdgeIds.size(); ++edgeIndex)
  {
    //performance is bad but I'm lazy
    std::set<localIndex> faces;
    getFacesAroundEdge( edgeIndex, allFacesToUniqueFace, faces );

    for( auto itr( faces.begin()); itr !=faces.end(); ++itr)
    {
      // Not smart - WARNING this assumes arrays default values is ZERO  
      if     ( faceToEdges( *itr, 0) == 0 )  faceToEdges( *itr, 0 ) = edgeIndex ;
      else if( faceToEdges( *itr, 1) == 0 )  faceToEdges( *itr, 1 ) = edgeIndex ;
      else if( faceToEdges( *itr, 2) == 0 )  faceToEdges( *itr, 2 ) = edgeIndex ;
      else if( faceToEdges( *itr, 3) == 0 )  faceToEdges( *itr, 3 ) = edgeIndex ;
      else assert( false ); 
    }
  }
  // TODO Asserts to check that everything went well 
  // 4 edges exactly per facet  

  return true;
}

bool HexMeshConnectivityBuilder::computeEdgesToFaces(ArrayOfSets<localIndex> & edgeToFaces ) const
{
  // 1 - Allocate - Without counting  -- Should we overallocate a bit ?
  // Probably for tet meshes one should count
  edgeToFaces.resize(0);
  edgeToFaces.reserve( nbEdges );
  edgeToFaces.reserveValues( 4 * nbEdges );
  for (int i = 0; i < nbEdges; ++i)
  {
    edgeToFaces.appendSet( 4 );
  }

  std::vector< localIndex > allFacesToUniqueFace;
  computeAllFacesToUniqueFace( allFacesToUniqueFace );

  for (int edgeIndex = 0; edgeIndex < m_uniqueEdgeIds.size(); ++edgeIndex)
  {
    std::set<localIndex> faces;
    getFacesAroundEdge( edgeIndex,  allFacesToUniqueFace, faces );
    
    for( auto itr( faces.begin()); itr !=faces.end(); ++itr)
    {
      edgeToFaces.insertIntoSet( edgeIndex, *itr );
    }
  }
  return true;
}

// We must have the CellNeighbors 
// Allocation is managed by the CellBlock 
bool HexMeshConnectivityBuilder::computeElementsToFacesOfCellBlocks() 
{
  getFaces();

  for( int curFace = 0; curFace < m_uniqueFaces.size(); ++curFace)
  {
    LocalFaceIndex id = m_uniqueFaces[curFace];

    LocalCellIndex c = id / Hex::nbFacets;
    HexFacetIndex face = id % Hex::nbFacets;
    auto blockCell = getBlockCellFromManagerCell(c) ;
    getCellBlock( blockCell.first ).setElementToFaces( blockCell.second, face, curFace );

    LocalFaceIndex idNeighbor = m_allFacesToNeighbors[id];
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


bool HexMeshConnectivityBuilder::computeElementsToEdgesOfCellBlocks() 
{
  getEdges();

  for (int edgeIndex = 0; edgeIndex < m_uniqueEdgeIds.size(); ++edgeIndex)
  {
    int last = (edgeIndex+1 < m_uniqueEdgeIds.size()) ? m_uniqueEdgeIds[edgeIndex+1] : m_allEdges.size();
    for( int i = m_uniqueEdgeIds[edgeIndex]; i < last; ++i)
    {
      LocalEdgeIndex id = m_allEdges[i].second;
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


HexCellBlockManager::HexCellBlockManager(string const &name, 
                                         dataRepository::Group *const parent) 
  :CellBlockManagerABC(name, parent),
   m_nodesPositions(0, 3)
{
  this->registerGroup<dataRepository::Group>(viewKeyStruct::cellBlocks());
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

dataRepository::Group *HexCellBlockManager::createChild(
  string const &GEOSX_UNUSED_PARAM(childKey), 
  string const &GEOSX_UNUSED_PARAM(childName))
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
  // Lazy computation - This should not do anything 
  // Who should be calling this mapping computations ? 
  m_theOneWhoDoesTheJob->computeElementsToFacesOfCellBlocks(); 
  m_theOneWhoDoesTheJob->computeElementsToEdgesOfCellBlocks();
  return ;
}

void HexCellBlockManager::setNumNodes(localIndex numNodes)
{
  m_nodesPositions.resize(numNodes);
  // TODO why is this allocated here and filled somewhere else?
  m_nodeLocalToGlobal.resize(numNodes);
  m_nodeLocalToGlobal.setValues<serialPolicy>(-1);
}

} // namespace geosx
