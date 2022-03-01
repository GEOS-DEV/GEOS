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

static constexpr HexVertexIndex facetVertex[6][4] = {
    {0, 1, 2, 3}, {4, 5, 6, 7}, {0, 1, 5, 4}, {1, 2, 6, 5}, {3, 2, 6, 7}, {0, 3, 7, 4}};

static constexpr HexVertexIndex edgeVertex[12][2]{
    {0, 1}, {0, 3}, {0, 4}, {1, 2}, {1, 5}, {2, 3}, {2, 6}, {3, 7}, {4, 5}, {4, 7}, {5, 6}, {6, 7}};

static constexpr HexVertexIndex quadFacetTriangleVertex[24][3] = {
    {0, 1, 2}, {0, 2, 3}, {0, 1, 3}, {1, 2, 3}, {4, 5, 6}, {4, 6, 7}, {4, 5, 7}, {5, 6, 7}, {0, 1, 5}, {0, 5, 4}, {0, 1, 4}, {1, 5, 4}, {1, 2, 6}, {1, 6, 5}, {1, 2, 5}, {2, 6, 5}, {3, 2, 6}, {3, 6, 7}, {3, 2, 7}, {2, 6, 7}, {0, 3, 7}, {0, 7, 4}, {3, 7, 4}, {3, 4, 0}};

static constexpr HexFacetIndex vertexAdjacentFacet[8][3] = {
    {0, 2, 5}, {0, 2, 3}, {0, 3, 4}, {0, 4, 5}, {1, 2, 5}, {1, 2, 3}, {1, 3, 4}, {1, 4, 5}};

static constexpr HexVertexIndex vertexAdjacentVertex[8][3] = {
    {1, 3, 4}, {2, 0, 5}, {3, 1, 6}, {0, 2, 7}, {7, 5, 0}, {4, 6, 1}, {5, 7, 2}, {6, 4, 3}};

static constexpr HexVertexIndex orientedFacetVertex[6][4] = {
    {0, 1, 2, 3}, {4, 7, 6, 5}, {0, 4, 5, 1}, {1, 5, 6, 2}, {3, 2, 6, 7}, {0, 3, 7, 4}};
};


class HexMeshConnectivityBuilder
{
public:

  HexMeshConnectivityBuilder( HexCellBlockManager & cellBlockManager ):
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
      m_blockCellIndexOffset[i] = prevSum + m_cellBlocks[i]->getElemToNodes().size();
    }
  };

  ~HexMeshConnectivityBuilder() = default;

  ArrayOfArrays<localIndex>        getFaceToNodes   () { return m_faceToNodes;    }
  array2d<geosx::localIndex>       getEdgeToNodes   () { return m_edgeToNodes;    }
  ArrayOfSets<localIndex>          getNodeToEdges   () { return m_nodeToEdges;    }
  ArrayOfSets<localIndex>          getNodeToFaces   () { return m_nodeToFaces;    }
  ArrayOfArrays<localIndex>        getNodeToElements() { return m_nodeToElements; }
  array2d<localIndex>              getFaceToElements() { return m_faceToElements; }
  ArrayOfArrays<geosx::localIndex> getFaceToEdges   () { return m_faceToEdges;    }
  ArrayOfSets<geosx::localIndex>   getEdgeToFaces   () { return m_edgeToFaces;    }

  // Compute functions for all maps
  // Do we really need to  implement it
  // What is the best order to compute this
  void computeAllMaps() {
    computeVertexGlobalToLocal();
    computeNbElements();
    computeCellNeighbors();
    computeEdges();
  };

  bool isRegularMesh() { return false; }


private:
  bool doneFaces          () { return false; }
  bool doneEdges          () { return false; }
  bool doneNodesToEdges   () { return false; }
  bool doneNodesToFaces   () { return false; }
  bool doneNodesToElements() { return false; }
  bool doneFacesToElements() { return false; }
  bool doneFacesToEdges   () { return false; }
  bool doneEdgesToFaces   () { return false; }
// bool doneEdgesToElements() { return false; } // A priori used by no one
  bool doneElementsToFaces() { return false; } // Here ? yes probably it depends on computation done here
  bool doneElementsToEdges() { return false; }

  // 
  bool computeVertexGlobalToLocal();
  bool computeNbElements();
  bool computeCellNeighbors();
  bool computeEdges          (); 
  
  
  bool computeNodesToEdges   ();
  bool computeNodesToFaces   ();
  bool computeNodesToElements();
  bool computeFacesToNodes   ();
  bool computeFacesToElements();
  bool computeFacesToEdges   ();
  bool computeEdgesToFaces   ();
//  bool computeEdgesToElements(); // A priori used by no one
  bool computeElementsToFaces(); // Done 
  bool computeElementsToEdges();

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
// Bool functions to know which map is already computed

CellBlockIndex getBlockFromCell( LocalCellIndex cellId ) const
{
  CellBlockIndex result = 0;
  while( cellId > m_blockCellIndexOffset[result] )
  {
    result++;
  }
  assert(result < nbCellBlocks()); // debug paranoia
  return result;
}



private:
// HexCellBlockManager for which we build the mapping
HexCellBlockManager const &  m_cellBlockManager;

LocalVertexIndex nbNodes = 0;
localIndex nbEdges = 0;
localIndex nbFaces = 0;
localIndex nbElements = 0;

std::vector< CellBlock * > m_cellBlocks;
std::vector< LocalCellIndex > m_blockCellIndexOffset;

// Do we need it ?
// Or areall indices local and we are good
std::unordered_map< GlobalVertexIndex, LocalVertexIndex > vertexGlobalToLocal; 


// Cell adjacencies - for all cell blocks
// I do not understand LvArray - and even less multidimensional LvArrays
// The advantages of multidimensional tables area mystery to me
array1d< LocalFaceIndex > m_cellNeighbors;

std::vector< EdgeInfo > m_edges;

// The mappings

// The cells/Elements are stored by the CellBlocks

// The mappings that this class is responsible of building
// TODO Why are storage strategies different? No real reason
// Do we really need to store them ? Probably not necessary for all of them
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

// The cells are stored by the CellBlocks

// TODO Check that there is only 1 CellBlock if there is only 1 type of cells -  Hex -
//      Looks likes there is 1 CellBlock per Region ?
//      What is the most efficient way to get a numbering of cells at the CellBlockManager level ?
// TODO Fast accessors to cells - in the CellBlocks - using the numbering

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

bool HexMeshConnectivityBuilder::computeNbElements()
{
  nbElements = 0;
  for (int i = 0; i < nbCellBlocks(); ++i)
  {
    CellBlock const & block = getCellBlock(i);
    nbElements += block.numElements();
  }
  return true;
}


bool HexMeshConnectivityBuilder::computeCellNeighbors()
{
  // 1 - Allocate
  SizeOfStuff nbTotalFaces = nbElements * Hex::nbFacets;

  // One goal is to compute the neighbors 
  m_cellNeighbors.resize( nbTotalFaces);
  
  // A big vector in which we put all facets
  // We use the LocalVertexIndex because the range is known go from 0 to nbNodes 
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
  // TODO Definitely not the fastest we can do - Use of HXTSort doable? 
  std::sort( allFaces.begin(), allFaces.end(), sortV0V1V2);

  // 4 - Counting + Set Cell Adjacencies
  SizeOfStuff nbBoundaryFaces = 0;
  SizeOfStuff nbInteriorFaces = 0;
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
      nbInteriorFaces++;
      ++i;
    }
    // If not this is a boundary face
    else
    {
      m_cellNeighbors[f.v[3]] = -1; //Replace by NO_ID - maxvalue of  something
      nbBoundaryFaces++;
    }
  }
  // Last facet is a boundary facet 
  if (i < nbTotalFaces) 
  {
    nbBoundaryFaces++;
  }
  // Set the number of faces 
  nbFaces = nbBoundaryFaces + nbInteriorFaces; 
  
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
  localIndex curFace = 0;
  for( int f = 0; f < m_cellNeighbors.size(); ++f)
  {
    LocalCellIndex c = f / Hex::nbFacets;
    LocalFaceIndex f1 = m_cellNeighbors[f];

    LocalCellIndex c1 = f1 != -1 ? f1 / Hex::nbFacets : c;
    if( (f1 != -1 && c < c1 ) || (f1 == -1) )
    {
      // We want the nodes of face f % 6 in cell c 
      localIndex blockIndex = getBlockFromCell( c );
      localIndex offset = blockIndex > 0 ? m_blockCellIndexOffset[blockIndex-1] : 0 ;
      
      CellBlock const & cB = getCellBlock( blockIndex );
      localIndex cellInBlock = c - offset;
      HexFacetIndex faceToStore = f % Hex::nbFacets;

      assert(cellInBlock >= 0); // debug 

      // Maybe the oriented facets would be useful - if they are the ones recomputed 
      // later one by the FaceManager
      m_faceToNodes[curFace][0] = cB.getElementNode(cellInBlock, Hex::facetVertex[faceToStore][0]);
      m_faceToNodes[curFace][1] = cB.getElementNode(cellInBlock, Hex::facetVertex[faceToStore][1]);
      m_faceToNodes[curFace][2] = cB.getElementNode(cellInBlock, Hex::facetVertex[faceToStore][2]);
      m_faceToNodes[curFace][3] = cB.getElementNode(cellInBlock, Hex::facetVertex[faceToStore][3]);

      ++curFace; 
    }
  }

  return true;
}



bool HexMeshConnectivityBuilder::computeEdges()
{
  // Let's use a different strategy than from the faces 
  // We are going for each edge to store a pair permitting to identify 
  // the 2 vertices: v0 * nbNodes + v1 - and the identification of the edge:  12 *idCell + e

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
  auto newEnd = std::unique( m_edges.begin(), m_edges.end(), equalEdges );  
  m_edges.erase( newEnd, m_edges.end());

  nbEdges = m_edges.size();
  return false;
}


bool HexMeshConnectivityBuilder::computeNodesToEdges()
{
  return false;
}

bool HexMeshConnectivityBuilder::computeNodesToFaces()
{
  return false;
}

bool HexMeshConnectivityBuilder::computeNodesToElements()
{
  // 1 -  Counting
  std::vector<unsigned int> nbElementsPerNode(nbNodes, 8);

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
  unsigned int nbValues = std::accumulate(nbElementsPerNode.begin(), nbElementsPerNode.end(), 0);

  //  2 - Allocating - No overallocation
  m_nodeToElements.reserve(nbNodes);
  m_nodeToElements.reserveValues(nbValues);

  // Set the adequate array for each node
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

bool computeFacesToElements() 
{
  return false;
}
bool computeFacesToEdges   () 
{
  return false;
}
bool computeEdgesToFaces   () 
{
  return false;
}
bool computeEdgesToElements() 
{
  return false;
}

// We must have the CellNeighbors 
// Same loop than FaceToNode - more or less - 
// TODO ?  Assign a FacetIndex to each of the faces and flag those to ignore in CellNeighbors
bool HexMeshConnectivityBuilder::computeElementsToFaces() 
{
  // 1 - Allocation is managed by the CellBlock 

  // 2 - Fill ElementToFace
  localIndex curFace = 0;
  for( int id = 0; id < m_cellNeighbors.size(); ++id)
  {
    LocalCellIndex c = id / Hex::nbFacets;
    LocalFaceIndex idNeighbor = m_cellNeighbors[id];
    
    if( id < idNeighbor || idNeighbor == -1 )
    {
      localIndex blockIndex = getBlockFromCell( c );
      localIndex offset = blockIndex > 0 ? m_blockCellIndexOffset[blockIndex-1] : 0 ;

      localIndex cellInBlock = c - offset;
      HexFacetIndex face = id % Hex::nbFacets;
      getCellBlock( blockIndex ).setElementToFaces( cellInBlock, face, curFace );
      ++curFace; 
    }
  }
  return true;
}


bool computeElementsToEdges() 
{
  return false;
}


# ifdef CODE_TOREMOVE

bool computeNodesToEdges()
{

  ArrayOfArrays<localIndex> toEdgesTemp(m_numNodes, maxEdgesPerNode());
  RAJA::ReduceSum<parallelHostReduce, localIndex> totalNodeEdges = 0;

  forAll<parallelHostPolicy>(m_numEdges, [&](localIndex const edgeID)
                             {
      toEdgesTemp.emplaceBackAtomic< parallelHostAtomic >( m_edgeToNodes( edgeID, 0 ), edgeID );
      toEdgesTemp.emplaceBackAtomic< parallelHostAtomic >( m_edgeToNodes( edgeID, 1 ), edgeID );
      totalNodeEdges += 2; });

  // Resize the node to edge map.
  m_nodeToEdges.resize(0);

  // Reserve space for the number of current nodes plus some extra.
  double const overAllocationFactor = 0.3;
  localIndex const entriesToReserve = (1 + overAllocationFactor) * m_numNodes;
  m_nodeToEdges.reserve(entriesToReserve);

  // Reserve space for the total number of face nodes + extra space for existing faces + even more space for new faces.
  localIndex const valuesToReserve = totalNodeEdges.get() + m_numNodes * getEdgeMapOverallocation() * (1 + 2 * overAllocationFactor);
  m_nodeToEdges.reserveValues(valuesToReserve);

  // Append the individual sets.
  for (localIndex nodeID = 0; nodeID < m_numNodes; ++nodeID)
  {
    m_nodeToEdges.appendSet(toEdgesTemp.sizeOfArray(nodeID) + getEdgeMapOverallocation());
  }

  ArrayOfSetsView<localIndex> const &toEdgesView = m_nodeToEdges.toView(); // FIXME why a view here?
  forAll<parallelHostPolicy>(m_numNodes, [&](localIndex const nodeID)
                             {
localIndex * const edges = toEdgesTemp[ nodeID ];
localIndex const numNodeEdges = toEdgesTemp.sizeOfArray( nodeID );
localIndex const numUniqueEdges = LvArray::sortedArrayManipulation::makeSortedUnique( edges, edges + numNodeEdges );
toEdgesView.insertIntoSet( nodeID, edges, edges + numUniqueEdges ); });

  return true;
}

bool computeNodesToFaces()
{
  localIndex const numFaces = m_faceToNodes.size();
  ArrayOfArrays<localIndex> nodeToFacesTemp(m_numNodes, maxFacesPerNode());
  RAJA::ReduceSum<parallelHostReduce, localIndex> totalNodeFaces = 0;

  forAll<parallelHostPolicy>(numFaces, [&](localIndex const faceID)
                             {
localIndex const numFaceNodes = m_faceToNodes.sizeOfArray( faceID );
totalNodeFaces += numFaceNodes;
for( localIndex a = 0; a < numFaceNodes; ++a )
{
  nodeToFacesTemp.emplaceBackAtomic< parallelHostAtomic >( m_faceToNodes( faceID, a ), faceID );
} });

  // Resize the node to face map.
  m_nodeToFaces.resize(0);

  // Reserve space for the number of nodes faces plus some extra.
  double const overAllocationFactor = 0.3;
  localIndex const entriesToReserve = (1 + overAllocationFactor) * m_numNodes;
  m_nodeToFaces.reserve(entriesToReserve);

  // Reserve space for the total number of node faces + extra space for existing nodes + even more space for new nodes.
  localIndex const valuesToReserve = totalNodeFaces.get() + m_numNodes * nodeMapExtraSpacePerFace() * (1 + 2 * overAllocationFactor);
  m_nodeToFaces.reserveValues(valuesToReserve);

  // Append the individual arrays.
  for (localIndex nodeID = 0; nodeID < m_numNodes; ++nodeID)
  {
    m_nodeToFaces.appendSet(nodeToFacesTemp.sizeOfArray(nodeID) + getFaceMapOverallocation());
  }

  forAll<parallelHostPolicy>(m_numNodes, [&](localIndex const nodeID)
                             {
localIndex * const faces = nodeToFacesTemp[ nodeID ];
localIndex const numNodeFaces = nodeToFacesTemp.sizeOfArray( nodeID );
localIndex const numUniqueFaces = LvArray::sortedArrayManipulation::makeSortedUnique( faces, faces + numNodeFaces );
m_nodeToFaces.insertIntoSet( nodeID, faces, faces + numUniqueFaces ); });

  return true;
}



/**
* @brief Fills the element to edges mappings of all the cells provided through @p cellBlocks.
* @param faceToEdges We need the face to edges mapping to get some edge index.
* @param cellBlocks The cell blocks for which the mappings will be constructed.
*/
void fillElementToEdgesOfCellBlocks(ArrayOfArrays<localIndex> const &faceToEdges,
                                  Group &cellBlocks)
{
for (localIndex iCellBlock = 0; iCellBlock < cellBlocks.numSubGroups(); ++iCellBlock)
{
  CellBlock &cellBlock = cellBlocks.getGroup<CellBlock>(iCellBlock);
  arrayView2d<localIndex const> const cellToFaces = cellBlock.getElemToFacesConstView();

  // We build the edges of each face of each cell,
  // so we can construct the cell to edges mapping.
  // Another implementation (used in other contexts) would use some edge signature
  // Some specific care is required not to insert edges twice (faces share edges).
  // to remove the duplicates.

  // Loop over the cells
  for (localIndex kc = 0; kc < cellBlock.numElements(); kc++)
  {
    int count = 0;
    for (localIndex kf = 0; kf < cellBlock.numFacesPerElement(); kf++)
    {
      // Loop over edges of each face
      localIndex const faceIndex = cellToFaces[kc][kf];
      for (localIndex ke = 0; ke < faceToEdges.sizeOfArray(faceIndex); ke++)
      {
        bool isUnique = true;
        localIndex edgeIndex = faceToEdges[faceIndex][ke];

        // Loop over edges that have already been added to the element.
        for (localIndex kec = 0; kec < count + 1; kec++)
        {
          // make sure that the edge has not been counted yet
          if (cellBlock.hasElementToEdges(kc, kec, edgeIndex))
          {
            isUnique = false;
            break;
          }
        }
        if (isUnique)
        {
          cellBlock.setElementToEdges(kc, count, edgeIndex);
          count++;
        }
      } // end edge loop
    }   // end face loop
  }     // end cell loop
}
}

} // ANONYMOUS

#endif



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
