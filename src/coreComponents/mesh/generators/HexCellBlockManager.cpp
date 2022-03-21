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
#include <cassert>


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

typedef array2d<localIndex, cells::NODE_MAP_PERMUTATION> CellVertexIndices;

// Use explicit aliases for indices in a hexahedron
typedef unsigned int HexVertexIndex;      // a vertex in a hex 0 to 8
typedef unsigned int HexEdgeIndex;        // a edge in a hex 0 to 12
typedef unsigned int HexFacetIndex;       // a facet in a hex 0 to 6

static const int NO_ID = -1;
static const int NO_FACE = -10;

/*  Hexahedron template
*  WARNING - Hex vertex numbering in GEOSX differs from the one used 
*  by most mesh datastructures - There are further variations in GEOSX itself 
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
static const unsigned int nbNodesPerFacet = 4;

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



/***************************************************************************************/
/***************************************************************************************/

/* Debugging functionalites */ 

void print(std::vector<localIndex> const &in)
{
  int count = 0 ;
  for (unsigned int i = 0; i < in.size(); ++i)
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
  for (unsigned int i = 0; i < in.size(); ++i)
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
  for( unsigned int i = 0; i < in.size(); ++i) {
    for(unsigned int j = 0; j < in.sizeOfSet(i);++j)
    {
      std::cout << std::setw(5) << std::left << in(i,j);
    }
    std::cout << std::endl;
  }
}

void print(array2d<localIndex> const &in)
{
  for (unsigned int i = 0; i < in.size(0); ++i)
  {
      std::cout << std::setw(5) << std::left << in(i,0);
      std::cout << std::setw(5) << std::left << in(i,1);
  }
  std::cout << std::endl
            << std::endl;
}
void print(std::vector<EdgeInfo> const &in)
{
  for (unsigned int  i = 0; i < in.size(); ++i)
  {
    std::cout << std::setw(5) << i
              << std::setw(5) << std::left << in[i].first 
              << std::setw(5) << std::left << in[i].second 
              << std::endl;
  }
  std::cout << std::endl
            << std::endl;
}

/***************************************************************************************/
/***************************************************************************************/


/**
 * @class MeshConnectivityBuilder
 * @brief The MeshConnectivityBuilder to build the connectivity maps 
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
 * TODO Implement specialization for regular hex mesh
 * 
 * TODO The storage of Faces and Edges is dependant on the Cell Types  
 * How do we manage the CellBlocks with different types of cells?
 * Strategy 1: allocate the space for hexahedra and keep a lot of invalid stuff 
 *             in these vectors - not the worst idea since these meshes should be hex-dominant
 *             do we need to store the cell type?
 * TODO For full tetrahedral meshes we need a dedicated implementation
 * 
 * Options : Store if the face is a triangle or quad ? 
 *  Store a  NO_FACE  value if  same storage space for hybrid meshes
 *  
 */
class MeshConnectivityBuilder
{
public:
  MeshConnectivityBuilder(CellBlockManagerBase & cellBlockManager);
  MeshConnectivityBuilder( const MeshConnectivityBuilder & ) = delete;
  MeshConnectivityBuilder & operator=( const MeshConnectivityBuilder & ) = delete;
  virtual ~MeshConnectivityBuilder() = default;

  localIndex numEdges() const {
    return nbEdges;
  }
  localIndex numFaces() const {
    return nbFaces;
  }

  // Dependent on Cell Type
  virtual bool computeFaces() = 0;
  virtual bool computeEdges() = 0;

  // Independent on Cell Type
  bool computeNodesToEdges( ArrayOfSets<localIndex> & result ) const; 
  bool computeEdgesToNodes( array2d<localIndex>& result ) const;  
  bool computeNodesToElements( ArrayOfArrays<localIndex>  & result ) const;

  // Dependent on Cell Type
  virtual bool computeNodesToFaces( ArrayOfSets<localIndex> & result ) const = 0 ;
  virtual bool computeEdgesToFaces( ArrayOfSets<localIndex> & result ) const  = 0;
  virtual bool computeFacesToNodes( ArrayOfArrays<localIndex> & result ) const = 0;
  virtual bool computeFacesToElements( array2d<localIndex> &  result ) const  = 0;
  virtual bool computeFacesToEdges( ArrayOfArrays<localIndex> & result ) const =0 ;

  // Why is this here? The information is on the CellBlocks
  virtual bool computeElementsToFacesOfCellBlocks() = 0; 
  virtual bool computeElementsToEdgesOfCellBlocks() = 0;

  bool debuggingComputeAllMaps() const;

  void printDebugInformation() const;

protected:
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
 
  void computeAllFacesToUniqueFace( std::vector< localIndex > & result ) const;

  virtual void getFacesAroundEdge( int edgeIndex,
                           std::vector< localIndex > const & allFacesToUniqueFace, 
                           std::set <localIndex> & faces ) const = 0;

protected:
/// Number of vertices 
SizeOfStuff nbNodes = 0;
/// All Elements in all CellBlocks
SizeOfStuff nbElements = 0;  

/// Number of unique edges - Size of m_uniqueEdges 
SizeOfStuff nbEdges = 0;
/// Number of unique faces - Size of m_uniqueFaces and m_isBoundaryFace
SizeOfStuff nbFaces = 0;

/// Cell Blocks on which the class operates - Size of nbBlocks
std::vector< CellBlock * > m_cellBlocks;

/// Offset for the numbering of all the cells - First value is the number of cells of block 0
/// Size of nbBlocks
std::vector< LocalCellIndex > m_blockCellIndexOffset;

/* Storage of a minimal set of information to iterate through
 * the faces while storing to which face of which cell they belong and which
 * is the neighbor face is the neighbor cell.
 * Use the numbering of cells managed by this class max is nbElements
 * Each face of each cell is encoded by 6 * cellIndex + faceIndexInCell
 * 
 * TODO Implement for tetrahedra (6 becomes 4)
 * TODO Define the strategy for hybrid FE meshes (hex, prism, pyramids, tets)
 */
std::vector< LocalFaceIndex > m_allFacesToNeighbors;  // 6 * nbElements
std::vector< localIndex > m_uniqueFaces;              // nbFaces
std::vector< bool > m_isBoundaryFace;                 // nbFaces

/* Storage of a minimal set of information to iterate through
 * the edges faces while storing to which face of which cell they belong to
 */
std::vector< EdgeInfo > m_allEdges;                 // 12 * nbElements
std::vector< localIndex > m_uniqueEdges;            // nbEdges

};

class HexMeshConnectivityBuilder : public MeshConnectivityBuilder
{
public:
  HexMeshConnectivityBuilder(CellBlockManagerBase & cellBlockManager):
    MeshConnectivityBuilder(cellBlockManager){};
  HexMeshConnectivityBuilder( const HexMeshConnectivityBuilder & ) = delete;
  HexMeshConnectivityBuilder & operator=( const HexMeshConnectivityBuilder & ) = delete;
  ~HexMeshConnectivityBuilder() = default;

  bool computeFaces() override;
  bool computeEdges() override;
  
  virtual bool computeElementsToFacesOfCellBlocks() override;
  virtual bool computeElementsToEdgesOfCellBlocks() override;

  bool computeNodesToFaces(ArrayOfSets<localIndex> &result) const override;
  bool computeEdgesToFaces(ArrayOfSets<localIndex> &result) const override;
  bool computeFacesToNodes(ArrayOfArrays<localIndex> &result) const override;
  bool computeFacesToElements(array2d<localIndex> &result) const override;
  bool computeFacesToEdges(ArrayOfArrays<localIndex> &result) const override;

protected:
  void getFacesAroundEdge(int edgeIndex,
                          std::vector<localIndex> const &allFacesToUniqueFace,
                          std::set<localIndex> &faces) const override;
};



MeshConnectivityBuilder::MeshConnectivityBuilder( CellBlockManagerBase & cellBlockManager )
{
  nbNodes = cellBlockManager.numNodes();
  
  dataRepository::Group & group = cellBlockManager.getCellBlocks();
  localIndex nbBlocks = group.numSubGroups();
  
  m_cellBlocks.resize(nbBlocks);
  m_blockCellIndexOffset.resize(nbBlocks);

  localIndex prevSum = 0;
  for (int i = 0; i < nbBlocks; ++i)
  {
    // TODO - Check that we indeed have CellBlock instantiation 
    CellBlock & cur = group.getGroup< CellBlock >( i );
    m_cellBlocks[i] = &cur;
    m_blockCellIndexOffset[i] = prevSum + m_cellBlocks[i]->numElements();
  }
  assert( nbBlocks > 0);

  nbElements =  m_blockCellIndexOffset[nbBlocks-1];
}

// TODO convert to log output
void MeshConnectivityBuilder::printDebugInformation() const
{ 
   std::cout << std::endl
            << std::endl;
  
  std::cout << " Number of blocks : " << m_cellBlocks.size() << std::endl
            << std::endl;
  print(m_blockCellIndexOffset);  

  std::cout << std::endl
            << std::endl;

  std::cout << " Total number of elements : " << nbElements << std::endl
            << std::endl;


  std::cout << " Number of unique faces : " << nbFaces << std::endl
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

  std::cout << " NB Edges : " << nbEdges << std::endl
            << std::endl;

  std::cout << " Edge Info " << std::endl
            << std::endl;

  print(m_allEdges);
  std::cout << std::endl
          << std::endl;
  
  print(m_uniqueEdges);
  std::cout << std::endl
            << std::endl;

  return;
}

std::pair< CellBlockIndex, LocalCellIndex> 
MeshConnectivityBuilder::getBlockCellFromManagerCell( LocalCellIndex cellId ) const
{
  CellBlockIndex b = 0;
  while( cellId > m_blockCellIndexOffset[b] )
  {
    b++;
  }
  localIndex offset = b > 0 ? m_blockCellIndexOffset[b-1] : 0;
  LocalCellIndex c = cellId - offset;

  return std::pair< CellBlockIndex, LocalCellIndex> (b, c);
}


void MeshConnectivityBuilder::computeAllFacesToUniqueFace(
   std::vector< localIndex > & allFacesToUniqueFace ) const
{
  allFacesToUniqueFace.resize( m_allFacesToNeighbors.size(), NO_FACE );

  for( unsigned int curFace = 0; curFace < m_uniqueFaces.size(); ++curFace)
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

// Cell type independent
bool MeshConnectivityBuilder::computeEdgesToNodes ( array2d<localIndex> & edgeToNodes ) const
{
  edgeToNodes.resize(nbEdges, 2);
  // Initialize contents
  edgeToNodes.setValues< serialPolicy >( NO_ID );
 
  for (unsigned int i = 0; i < m_uniqueEdges.size(); ++i)
  {
    EdgeInfo e = m_allEdges[m_uniqueEdges[i]];
    LocalVertexIndex v0 = e.first / nbNodes;
    LocalVertexIndex v1 = e.first % nbNodes;

    edgeToNodes[i][0] = v0;
    edgeToNodes[i][1] = v1;
  }
  return true;
}

// Cell type independent - There is no need for this to be an ArrayOfSets
// We fill it with unique (maybe sorted) values
bool MeshConnectivityBuilder::computeNodesToEdges( ArrayOfSets<localIndex> & nodeToEdges ) const
{ 
  // 1 - Counting 
  // TODO Can be skipped for hexahedral meshes - 6 for regular nodes - 10 tops for singular nodes
  std::vector<unsigned int> nbEdgesPerNode(nbNodes, 0);
  for( unsigned int i = 0; i < m_uniqueEdges.size(); ++i)
  {
    EdgeInfo e = m_allEdges[ m_uniqueEdges[i] ];
    LocalVertexIndex v0 = e.first / nbNodes;
    LocalVertexIndex v1 = e.first % nbNodes;
    nbEdgesPerNode[v0]++;
    nbEdgesPerNode[v1]++;
  }
  localIndex valuesToReserve = std::accumulate(nbEdgesPerNode.begin(), nbEdgesPerNode.end(), 0);

  // 2 - Allocating 
  nodeToEdges.resize( 0 );
  nodeToEdges.reserve( nbNodes );
  nodeToEdges.reserveValues( valuesToReserve );

  // Append and set the capacity of the individual arrays.
  for (localIndex n = 0; n < nbNodes; ++n)
  {
    nodeToEdges.appendSet( nbEdgesPerNode[n] );
  }
  
  // 3 - Filling
  for(unsigned  int i = 0; i < m_uniqueEdges.size(); ++i)
  {
    EdgeInfo e = m_allEdges[ m_uniqueEdges[i] ];
    LocalVertexIndex v0 = e.first / nbNodes;
    LocalVertexIndex v1 = e.first % nbNodes;

    nodeToEdges.insertIntoSet( v0, i); 
    nodeToEdges.insertIntoSet( v1, i);
  }
  return true;
}

// Cell type independent
bool MeshConnectivityBuilder::computeNodesToElements( ArrayOfArrays<localIndex> & nodeToElements ) const
{
  // 1 -  Counting
  // TODO Can be skipped for hexahedral meshes - 8 for regular nodes - 12 tops for singular nodes
  std::vector<unsigned int> nbElementsPerNode(nbNodes, 0);
  for (unsigned int i = 0; i < nbCellBlocks(); ++i)
  {
    CellBlock const & block = getCellBlock(i);
    CellVertexIndices const & cells = block.getElemToNodes();
    localIndex nbVertices = block.numNodesPerElement();

    for (int j = 0; j < block.numElements(); ++j)
    {
      for (unsigned int v = 0; v < nbVertices; ++v)
      {
        nbElementsPerNode[ cells( j, v ) ]++;   
      }
    }
  }
  localIndex nbValues = std::accumulate(nbElementsPerNode.begin(), nbElementsPerNode.end(), 0);

  //  2 - Allocating - No overallocation
  nodeToElements.resize(0);
  nodeToElements.resize(nbNodes);
  nodeToElements.reserveValues( nbValues ); // Does this accelerate allocation?

  for( int i = 0; i < nbNodes; ++i )
  {
    nodeToElements.setCapacityOfArray(i, nbElementsPerNode[i]);
  }

  // 3 - Set the values
  for (unsigned int i = 0; i < nbCellBlocks(); ++i)
  {
    CellBlock const & block = getCellBlock(i);
    CellVertexIndices const & cells = block.getElemToNodes();
    localIndex nbVertices = block.numNodesPerElement();

    for (int j = 0; j < block.numElements(); ++j)
    {
      for (unsigned int v = 0; v < nbVertices; ++v)
      {
        localIndex const nodeIndex = cells(j, v);
        nodeToElements.emplaceBack(nodeIndex, j);
      }
    }
  }
  return true;
}


/**
 * @brief Compute all the faces of the cells of the cell blocks
 * Fills the face information 
 * TODO Be able to get out and return an error message if 3 cells have the same face
 * 
 * TODO CellType dependant 
 * TODO What management for triangles? Fill with NO_VERTEX the 4th vertex and 
 * have consistent cell descriptions?
 * The block knows the cell type, we could use this? 
 * Maybe sort CellBlocks by cell type, and template things
 */
bool HexMeshConnectivityBuilder::computeFaces()
{
  m_allFacesToNeighbors.resize( 0);
  m_uniqueFaces.resize( 0 ); 
  m_isBoundaryFace.resize( 0 ); 
  
  // 1 - Allocate - 
  SizeOfStuff nbTotalFaces = nbElements * Hex::nbFacets;
  m_allFacesToNeighbors.resize( nbTotalFaces);
  
  // To collect and sort the facets 
  std::vector< FaceInfo > allFaces ( nbTotalFaces );

  // 2 - Fill 
  localIndex curFace = 0;
  for (unsigned int i = 0; i < nbCellBlocks(); ++i)
  {
    CellBlock const & block = getCellBlock(i);
    CellVertexIndices const & cells = block.getElemToNodes();

    for (int j = 0; j < cells.size(0); ++j)
    {
      for (unsigned int f = 0; f < Hex::nbFacets; ++f)
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
        // Last slot is used to identify the facet in its cell 
        allFaces[curFace].v[3]= Hex::nbFacets * j + f; 
        
        curFace++;
      }
    }
  }

  // 3 - Sort 
  // TODO Definitely not the fastest we can do - Use of HXTSort if possible? Has LvArray faster sort? 
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
  m_allEdges.resize( 0 );
  m_uniqueEdges.resize( 0 );

  // 1 - Allocate 
  SizeOfStuff nbAllEdges = nbElements * Hex::nbEdges;
  m_allEdges.resize( nbAllEdges );
  
  // 2 - Get all edges
  localIndex cur = 0 ;
  for (unsigned int i = 0; i < nbCellBlocks(); ++i)
  {
    CellBlock const & block = getCellBlock(i);
    CellVertexIndices const & cells = block.getElemToNodes();

    for (int j = 0; j < cells.size(0); ++j)
    {
      for(unsigned int e = 0; e < Hex::nbEdges; ++e)
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
  unsigned int i = 0;
  while( i < m_allEdges.size() )
  {
    m_uniqueEdges.push_back(i);
    unsigned int j = i+1;
    while( j < m_allEdges.size() && equalEdgeInfo(m_allEdges[i], m_allEdges[j]) )
    { 
      ++j; 
    }
    i = j;
  }
  nbEdges = m_uniqueEdges.size();
  return true;
}


bool HexMeshConnectivityBuilder::computeFacesToNodes( ArrayOfArrays<localIndex> & faceToNodes ) const
{
  // 1 - Allocate - No overallocation
  faceToNodes.resize(0);
  faceToNodes.resize( nbFaces, Hex::nbNodesPerFacet );

  // TODO Check if this resize or reserve the space for the inner arrays
  // In doubt resize inner arrays
  for (int i = 0; i < nbFaces; ++i)
  {
    faceToNodes.resizeArray(i, Hex::nbNodesPerFacet);
  }

  // 2 - Fill FaceToNode  -- Could be avoided and done when required from adjacencies
  for( unsigned int curFace = 0; curFace < m_uniqueFaces.size(); ++curFace)
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


bool HexMeshConnectivityBuilder::computeNodesToFaces( ArrayOfSets<localIndex> & nodeToFaces ) const
{
  ArrayOfArrays<localIndex> faceToNodes; 
  computeFacesToNodes( faceToNodes );

  // 1 - Counting  - Quite unecessary for hexahedral meshes 
  //( 12 for regular nodes  - TODO check a good max for singular node )
  std::vector<unsigned int> nbFacesPerNode(nbNodes, 0);
  for(unsigned int i = 0; i < faceToNodes.size(); ++i)
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
  nodeToFaces.reserveValues( valuesToReserve );
  for (localIndex n = 0; n < nbNodes; ++n)
  {
    nodeToFaces.appendSet( nbFacesPerNode[n] );
  }
  
  // 3 - Filling 
  for(unsigned int i = 0; i < faceToNodes.size(); ++i)
  {
    nodeToFaces.insertIntoSet( faceToNodes[i][0],i );
    nodeToFaces.insertIntoSet( faceToNodes[i][1],i );
    nodeToFaces.insertIntoSet( faceToNodes[i][2],i );
    nodeToFaces.insertIntoSet( faceToNodes[i][3],i );
  }

  return false;
}


// TODO We have a problem - where are the Element index valid ?
bool HexMeshConnectivityBuilder::computeFacesToElements( array2d<localIndex> & faceToElements ) const
{
  faceToElements.resize(0);
  faceToElements.resize(nbFaces, 2);
  // Initialize contents to NO_ID 
  faceToElements.setValues< serialPolicy >( NO_ID );

  for(unsigned int curFace = 0; curFace < m_uniqueFaces.size(); ++curFace)
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

  // For each duplicate of this edge
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


bool HexMeshConnectivityBuilder::computeFacesToEdges( ArrayOfArrays<localIndex> & faceToEdges) const 
{
  // 1 - Allocate
  faceToEdges.resize(0);
  // TODO Check if this resize or reserve the space for the inner arrays
  faceToEdges.resize( nbFaces, Hex::nbEdgesPerFacet); 
  // In doubt resize inner arrays
  for (int i = 0; i < nbFaces; ++i)
  {
    faceToEdges.resizeArray(i, Hex::nbEdgesPerFacet );
    // Initialize the values - Is there a simpler way to do this?
    faceToEdges[i][0] = NO_ID;
    faceToEdges[i][1] = NO_ID;
    faceToEdges[i][2] = NO_ID;
    faceToEdges[i][3] = NO_ID;
  }

  std::vector< localIndex > allFacesToUniqueFace;
  computeAllFacesToUniqueFace( allFacesToUniqueFace );

  for (unsigned int edgeIndex = 0; edgeIndex < m_uniqueEdges.size(); ++edgeIndex)
  {
    std::set<localIndex> faces;
    getFacesAroundEdge( edgeIndex, allFacesToUniqueFace, faces );

    for( auto itr( faces.begin()); itr !=faces.end(); ++itr)
    { 
      if     ( faceToEdges[*itr][0] == NO_ID )  faceToEdges[*itr][0] = edgeIndex ;
      else if( faceToEdges[*itr][1] == NO_ID )  faceToEdges[*itr][1] = edgeIndex ;
      else if( faceToEdges[*itr][2] == NO_ID )  faceToEdges[*itr][2] = edgeIndex ;
      else if( faceToEdges[*itr][3] == NO_ID )  faceToEdges[*itr][3] = edgeIndex ;
      else assert( false ); 
    }
  }
  return true;
}

bool HexMeshConnectivityBuilder::computeEdgesToFaces(ArrayOfSets<localIndex> & edgeToFaces ) const
{
  // 1 - Allocate - Without counting  -- Should we overallocate a bit ?
  // TODO Count for tet meshes one should count
  int numberOfFacesAroundEdge = 5; // In a hexahedral mesh this is enough
  edgeToFaces.resize( 0 );
  edgeToFaces.reserve( nbEdges );
  edgeToFaces.reserveValues( nbEdges * numberOfFacesAroundEdge );

  // Append and set the capacity of the individual arrays.
  for (int i = 0; i < nbEdges; ++i)
  {
    edgeToFaces.appendSet(numberOfFacesAroundEdge);
  }

  // 2 - Get the mapping
  std::vector< localIndex > allFacesToUniqueFace;
  computeAllFacesToUniqueFace( allFacesToUniqueFace );

  for (unsigned int edgeIndex = 0; edgeIndex < m_uniqueEdges.size(); ++edgeIndex)
  {
    std::set<localIndex> faces;
    getFacesAroundEdge( edgeIndex, allFacesToUniqueFace, faces );
    
    for( auto itr( faces.begin()); itr !=faces.end(); ++itr)
    {
      edgeToFaces.insertIntoSet(edgeIndex, *itr);
    }
  }
  return true;
}

// Allocation is managed by the CellBlock 
bool HexMeshConnectivityBuilder::computeElementsToFacesOfCellBlocks() 
{
  for(unsigned int curFace = 0; curFace < m_uniqueFaces.size(); ++curFace)
  {
    LocalFaceIndex id = m_uniqueFaces[curFace];

    LocalCellIndex c = id / Hex::nbFacets;
    HexFacetIndex face = id % Hex::nbFacets;
    auto blockCell = getBlockCellFromManagerCell(c);

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


// Allocation is managed by the CellBlock 
bool HexMeshConnectivityBuilder::computeElementsToEdgesOfCellBlocks() 
{
  for (unsigned int edgeIndex = 0; edgeIndex < m_uniqueEdges.size(); ++edgeIndex)
  {
    unsigned int first = m_uniqueEdges[edgeIndex];
    unsigned int last = m_allEdges.size();

    if( edgeIndex+1 < m_uniqueEdges.size() ) {
      last = m_uniqueEdges[edgeIndex+1];
    }

    for(unsigned int i = first ; i < last; ++i)
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

HexCellBlockManager::HexCellBlockManager(string const &name, 
                                         dataRepository::Group *const parent) 
  :CellBlockManagerBase(name, parent)
{
  m_theOneWhoDoesTheJob = nullptr; 
}

HexCellBlockManager::~HexCellBlockManager()
{
  delete m_theOneWhoDoesTheJob;
  m_theOneWhoDoesTheJob = nullptr;
}

void HexCellBlockManager::buildMaps()
{
  // Create the MeshConnectivityBuilder after the CellBlocks are 
  // filled and cell types are known
  // TODO Change the instanciated class depending on the incoming types of cells
  m_theOneWhoDoesTheJob = new HexMeshConnectivityBuilder( *this );

  // If not called here - the number of faces and edges are not 
  // available resize EdgeManager and FaceManager.
  m_theOneWhoDoesTheJob->computeFaces();
  m_theOneWhoDoesTheJob->computeEdges();
  
  m_theOneWhoDoesTheJob->printDebugInformation(); 

  // Otherwise this are not called. Who should?
  m_theOneWhoDoesTheJob->computeElementsToFacesOfCellBlocks(); 
  m_theOneWhoDoesTheJob->computeElementsToEdgesOfCellBlocks();

  return ;
}

localIndex HexCellBlockManager::numEdges() const
{
  return m_theOneWhoDoesTheJob->numEdges();
}
localIndex HexCellBlockManager::numFaces() const
{
  return m_theOneWhoDoesTheJob->numFaces();
}

array2d<localIndex> HexCellBlockManager::getEdgeToNodes() const
{
  array2d<localIndex> result;
  m_theOneWhoDoesTheJob->computeEdgesToNodes(result);
  return result;
}

ArrayOfSets<localIndex> HexCellBlockManager::getEdgeToFaces() const
{
  ArrayOfSets<localIndex> result;
  m_theOneWhoDoesTheJob->computeEdgesToFaces(result);
  return result;
}
ArrayOfArrays<localIndex> HexCellBlockManager::getFaceToNodes() const
{
  ArrayOfArrays<localIndex> result;
  m_theOneWhoDoesTheJob->computeFacesToNodes(result);
  return result;
}
ArrayOfArrays<localIndex> HexCellBlockManager::getFaceToEdges() const
{
  ArrayOfArrays<localIndex> result;
  m_theOneWhoDoesTheJob->computeFacesToEdges(result);
  return result;
}
array2d<localIndex> HexCellBlockManager::getFaceToElements() const
{
  array2d<localIndex> result;
  m_theOneWhoDoesTheJob->computeFacesToElements(result);
  return result;
}
ArrayOfSets<localIndex> HexCellBlockManager::getNodeToEdges() const
{
  ArrayOfSets<localIndex> result;
  m_theOneWhoDoesTheJob->computeNodesToEdges(result);
  return result;
}
ArrayOfSets<localIndex> HexCellBlockManager::getNodeToFaces() const
{
  ArrayOfSets<localIndex> result;
  m_theOneWhoDoesTheJob->computeNodesToFaces(result);
  return result;
}
ArrayOfArrays<localIndex> HexCellBlockManager::getNodeToElements() const
{
  ArrayOfArrays<localIndex> result;
  m_theOneWhoDoesTheJob->computeNodesToElements(result);
  return result;
}

} // namespace geosx
