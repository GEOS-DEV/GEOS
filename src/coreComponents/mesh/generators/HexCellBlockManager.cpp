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

typedef localIndex CellBlockIndex;

typedef localIndex SizeOfStuff;


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

  HexMeshConnectivityBuilder( HexCellBlockManager const & cellBlockManager ):
    m_cellBlockManager(cellBlockManager),
    m_cellBlocks( m_cellBlockManager.getCellBlocks() )
  {

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
  bool doneEdgesToElements() { return false; } // A priori used by no one
  bool doneElementsToFaces() { return false; } // Here ? yes probably it depends on computation done here
  bool doneElementsToEdges() { return false; }

  bool computeVertexGlobalToLocal(); // Done
  bool computeNbElements();          // Done
  bool computeCellNeighbors();

  bool computeEdges          ();
  bool computeNodesToEdges   ();
  bool computeNodesToFaces   ();
  bool computeNodesToElements();  //  Done 
  bool computeFacesToNodes   ();
  bool computeFacesToElements();
  bool computeFacesToEdges   ();
  bool computeEdgesToFaces   ();
  bool computeEdgesToElements(); // A priori used by no one
  bool computeElementsToFaces(); // Here ? yes probably it depends on computation done here
  bool computeElementsToEdges();

  unsigned int nbCellBlocks() {
    return m_cellBlocks.numSubGroups();
  } 
  CellBlock const & getCellBlock( unsigned int id ){
    // Maybe need for dynamic cast
    return m_cellBlocks.getGroup<CellBlock>( id );
  }

private:
// Bool functions to know which map is already computed



private:
// HexCellBlockManager for which we build the mapping
HexCellBlockManager const &  m_cellBlockManager;
Group const & m_cellBlocks;

LocalVertexIndex nbNodes = 0;
localIndex nbEdges = 0;
localIndex nbFaces = 0;
localIndex nbElements = 0;

// Do we need it ?
// Or areall indices local and we are good
std::unordered_map< GlobalVertexIndex, LocalVertexIndex > vertexGlobalToLocal; 


// Cell adjacencies - for all cell blocks
// I do not understand LvArray - and even less multidimensional LvArrays
// The advantages of multidimensional tables area mystery to me
array1d< LocalFaceIndex > m_cellNeighbors;

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


// We need a fast a smart way to sort the facets
// while keeping information on where this facets is coming from
// A similar strategy is possible for tet - hex - mixed element meshes

// A facet is a set of sorted nodes - maybe *n or *n*n - to reduce number of sorts
// Use sorting with adequate functors

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

  // OOOOOKKK we are going to need this
  std::vector< CellBlock const * > blocks ( nbCellBlocks(), nullptr );
  std::vector< LocalCellIndex > offsets ( nbCellBlocks() );
  LocalCellIndex prevSum = 0;
  for (int i = 0; i < nbCellBlocks(); ++i)
  {
    blocks[i] = &getCellBlock(i);
    // TODO Accessors at CellBlock Level - I doubt what happens here 
    offsets[i] = prevSum + blocks[i]->getElemToNodes().size();
  }
  
  // 2 - Fill FaceToNode  -- Could be avoided and done when required from adjacencies
  localIndex curFace = 0;
  for( int f= 0; f < m_cellNeighbors.size(); ++f)
  {
    LocalCellIndex c = f / Hex::nbFacets;
    LocalFaceIndex f1 = m_cellNeighbors[f];

    LocalCellIndex c1 = f1 != -1 ? f1 / Hex::nbFacets : c;
    if( (f1 != -1 && c < c1 ) || (f1 == -1) )
    {
      // We want the nodes of face f % 6 in cell c 
      // TODO Add accessor and functions to ease and clarify this 
      localIndex blockIndex = 0;
      localIndex offset = 0;
      while( c > offsets[blockIndex] )
      {
        offset = offsets[blockIndex];
        blockIndex++;
        assert(blockIndex <nbCellBlocks()); // debug paranoia
      }
      
      CellBlock const * cB = blocks[blockIndex];
      localIndex cellInBlock = c - offset;
      HexFacetIndex faceToStore = f % Hex::nbFacets;

      assert(cellInBlock >= 0); // debug 

      // Maybe the oriented facets would be useful - if they are the ones recomputed 
      // later one by the FaceManager
      m_faceToNodes[curFace][0] = cB->getElementNode(cellInBlock, Hex::facetVertex[faceToStore][0]);
      m_faceToNodes[curFace][1] = cB->getElementNode(cellInBlock, Hex::facetVertex[faceToStore][1]);
      m_faceToNodes[curFace][2] = cB->getElementNode(cellInBlock, Hex::facetVertex[faceToStore][2]);
      m_faceToNodes[curFace][3] = cB->getElementNode(cellInBlock, Hex::facetVertex[faceToStore][3]);

      ++curFace; 
    }
  }

  return true;
}


bool HexMeshConnectivityBuilder::computeEdges()
{
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
bool computeElementsToFaces() 
{
  return false;
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





namespace
{

using namespace geosx;
using geosx::dataRepository::Group;


/**
* @brief Convenience programming structure that holds the nodes for a given Face. No information about which cell though.
*
* This structure holds `<` and `==` operators such that equal instances in a sorted array
* are meant to be consecutive. (see `std::unique` for example).
* Two instances with same nodes in different order are considered equal.
* Such a case often happen since the same surface is shared by two adjacent elements
* in conformal meshes.
*/
struct NodesAndElementOfFace
{
NodesAndElementOfFace(array1d<localIndex> nodes_, localIndex element_, localIndex iCellBlock_, localIndex iFace_) : nodes(nodes_),
                                                                                                                    element(element_),
                                                                                                                    iCellBlock(iCellBlock_),
                                                                                                                    iFace(iFace_),
                                                                                                                    sortedNodes(nodes_)
{
  std::sort(sortedNodes.begin(), sortedNodes.end());
}

/**
* @brief Imposes an ordering on NodesAndElementOfFace.
* @param [in] rhs the NodesAndElementOfFace to compare against.
* @return a boolean.
*/
bool operator<(NodesAndElementOfFace const &rhs) const
{
  // Using some standard comparison like vector::operator<.
  // Two subsequent NodesAndElementOfFace may still be equal
  // We are consistent with operator==, which is what we require.
  return std::lexicographical_compare(sortedNodes.begin(), sortedNodes.end(),
                                      rhs.sortedNodes.begin(), rhs.sortedNodes.end());
}

/**
* @brief Two NodesAndElementOfFace instances are considered equal if they share the same node set, whatever the order.
* @param [in] rhs the NodesAndElementOfFace to compare against.
* @return a boolean.
*/
bool operator==(NodesAndElementOfFace const &rhs) const
{
  // Comparing term by term like STL does.
  return (sortedNodes.size() == rhs.sortedNodes.size() && std::equal(sortedNodes.begin(), sortedNodes.end(), rhs.sortedNodes.begin()));
}

/// The list of nodes describing the face.
array1d<localIndex> nodes;

/**
* @brief The element to which this face belongs.
*
* Each face may belong to multiple elements.
* But during the identification process (we loop on each face of each element),
* we store the cell we are iterating on.
* The we'll be able to identify the duplicated faces because we also have the nodes.
*/
localIndex element;

/**
* @brief Cell block index
*
* During the process, we need to know form which cell block the instance was created.
*/
localIndex iCellBlock;

/**
* @brief Face index
*
* During the process, we need to know what was the face index when this instance was created.
*/
localIndex iFace;

private:
/// Sorted nodes describing the face; only for comparison/sorting reasons.
array1d<localIndex> sortedNodes;
};

/**
* @brief Return the total number of unique faces and fill in the uniqueFaceOffsets array.
* @param [in] lowestNodeToFaces An array of size numNodes of arrays of NodesAndElementOfFace associated with each node.
* @param [out] uniqueFaceOffsets An array of size numNodes + 1. After this function returns node i contains
*              faces with IDs ranging from uniqueFaceOffsets[ i ] to uniqueFaceOffsets[ i + 1 ] - 1.
* @return return total number of faces
*/
localIndex calculateTotalNumberOfFaces(ArrayOfArraysView<NodesAndElementOfFace const> const &lowestNodeToFaces,
                                      arrayView1d<localIndex> const &uniqueFaceOffsets)
{
localIndex const numNodes = lowestNodeToFaces.size();
GEOSX_ERROR_IF_NE(numNodes, uniqueFaceOffsets.size() - 1);
uniqueFaceOffsets.setValues<serialPolicy>(0.);

// Loop over all the nodes.
forAll<parallelHostPolicy>(numNodes, [&](localIndex const nodeID)
                            {
localIndex const numFaces = lowestNodeToFaces.sizeOfArray( nodeID );
// Since lowestNodeToFaces[ nodeID ] is sorted we can compare subsequent entries
// count up the unique entries. Since each face can appear at most twice if we find a match we
// can skip the next entry as well.
localIndex & numUniqueFaces = uniqueFaceOffsets[ nodeID + 1 ];
for( localIndex j = 0; j < numFaces; ++j )
{
  ++numUniqueFaces;
  if( j < numFaces - 1 )
  {
    if( lowestNodeToFaces( nodeID, j ) == lowestNodeToFaces( nodeID, j + 1 ) )
    {
      ++j;
    }
  }
} });
// At this point uniqueFaceOffsets[ i ] holds the number of unique face associated with node i - 1.
// Perform an inplace prefix-sum to get the unique face offset.
RAJA::inclusive_scan_inplace<parallelHostPolicy>(uniqueFaceOffsets.begin(), uniqueFaceOffsets.end());
return uniqueFaceOffsets.back();
}

/**
* @brief Copies the nodes from @p nodesAndElementOfFace into @p faceToNodes[@p faceID ].
* @param [in] faceID The face index.
* @param [in] nodesAndElementOfFace The nodes and element for @p faceID.
* @param [out] faceToNodes No input data is used.
*/
void insertFaceToNodesEntry(localIndex const faceID,
                          NodesAndElementOfFace const &nodesAndElementOfFace,
                          ArrayOfArrays<localIndex> &faceToNodes)
{
localIndex const numFaceNodes = nodesAndElementOfFace.nodes.size();
// FIXME The size should be OK because it's been allocated previously.
for (localIndex i = 0; i < numFaceNodes; ++i)
{
  faceToNodes[faceID][i] = nodesAndElementOfFace.nodes[i];
}
GEOSX_ASSERT_EQ(numFaceNodes, faceToNodes.sizeOfArray(faceID));
GEOSX_DEBUG_VAR(numFaceNodes);
}

/**
* @brief Fills the face to nodes map and face to element maps
* @param [in] lowestNodeToFaces and array of size numNodes of arrays of NodesAndElementOfFace associated with each node.
* @param [in] uniqueFaceOffsets an array containing the unique ID of the first face associated with each node.
* @param [inout] faceToElements the face to element map.
* @param [inout] faceToNodes the face to node map.
*/
void populateFaceMaps(ArrayOfArraysView<NodesAndElementOfFace const> const &lowestNodeToFaces,
                    arrayView1d<localIndex const> const &uniqueFaceOffsets,
                    ArrayOfArrays<localIndex> &faceToNodes,
                    arrayView2d<localIndex> &faceToElements)
{
GEOSX_MARK_FUNCTION;

localIndex const numNodes = lowestNodeToFaces.size();
localIndex const numUniqueFaces = uniqueFaceOffsets.back();
GEOSX_ERROR_IF_NE(numNodes, uniqueFaceOffsets.size() - 1);
GEOSX_ERROR_IF_NE(numUniqueFaces, faceToNodes.size());
GEOSX_ERROR_IF_NE(numUniqueFaces, faceToElements.size(0));
GEOSX_ERROR_IF_NE(2, faceToElements.size(1));

// loop over all the nodes.
forAll<parallelHostPolicy>(numNodes, [&](localIndex const nodeID)
                            {
localIndex const numFaces = lowestNodeToFaces.sizeOfArray( nodeID );
// loop over all the `NodesAndElementOfFace` associated with the node.
for( localIndex j = 0, curFaceID = uniqueFaceOffsets[nodeID]; j < numFaces; ++j, ++curFaceID )
{
  // If two subsequent `NodesAndElementOfFace` compare equal then they describe an interior face.
  // It must therefore be considered "twice".
  // Otherwise it's a boundary face, and "once" is enough.
  NodesAndElementOfFace const & f0 = lowestNodeToFaces( nodeID, j );
  insertFaceToNodesEntry( curFaceID, f0, faceToNodes );
  faceToElements( curFaceID, 0 ) = f0.element;
  faceToElements( curFaceID, 1 ) = -1; // TODO Make a constant

  // This is where we check for the two subsequent faces (when they exist).
  if( j < numFaces - 1 )
  {
    NodesAndElementOfFace const & f1 = lowestNodeToFaces( nodeID, j + 1 );
    if( f0 == f1 )
    {
      faceToElements( curFaceID, 1 ) = f1.element;
      ++j;
    }
  }
} });
}

/**
* @brief Resize the face maps
* @param [in] lowestNodeToFaces and array of size numNodes of arrays of NodesAndElementOfFace associated with each node.
* @param [in] uniqueFaceOffsets an containing the unique face IDs for each node in lowestNodeToFaces.
* @param [out] faceToNodeMap the map from faces to nodes. This function resizes the array appropriately.
* @param [out] faceToElemMap the map from faces to elements. This function resizes the array appropriately.
*/
void resizeFaceMaps(ArrayOfArraysView<NodesAndElementOfFace const> const &lowestNodeToFaces,
                  arrayView1d<localIndex const> const &uniqueFaceOffsets,
                  ArrayOfArrays<localIndex> &faceToNodeMap,
                  ArrayOfArrays<localIndex> &faceToEdgesMap,
                  array2d<localIndex> &faceToElemMap)
{
GEOSX_MARK_FUNCTION;

localIndex const numNodes = lowestNodeToFaces.size();
localIndex const numUniqueFaces = uniqueFaceOffsets.back();
array1d<localIndex> numNodesPerFace(numUniqueFaces);
RAJA::ReduceSum<parallelHostReduce, localIndex> totalFaceNodes(0.0);

// loop over all the nodes.
forAll<parallelHostPolicy>(numNodes, [&](localIndex const nodeID)
                            {
localIndex const numFaces = lowestNodeToFaces.sizeOfArray( nodeID );
// loop over all the NodesAndElementOfFace associated with the node
for( localIndex j = 0, curFaceID = uniqueFaceOffsets[ nodeID ]; j < numFaces; ++j, ++curFaceID )
{
  const NodesAndElementOfFace & f0 = lowestNodeToFaces( nodeID, j );
  numNodesPerFace[curFaceID] = f0.nodes.size();
  totalFaceNodes += numNodesPerFace[curFaceID];

  if( ( j < numFaces - 1 ) and ( f0 == lowestNodeToFaces( nodeID, j + 1 ) ) )
  {
    ++j;
  }
} });

// Resize the face to node map.
faceToNodeMap.resize(0);

// Reserve space for the number of current faces plus some extra.
double const overAllocationFactor = 0.3;
localIndex const entriesToReserve = (1 + overAllocationFactor) * numUniqueFaces; // TODO why this allocation factor
faceToNodeMap.reserve(entriesToReserve);

// Reserve space for the total number of face nodes + extra space for existing faces + even more space for new faces.
localIndex const valuesToReserve = totalFaceNodes.get() + numUniqueFaces * HexCellBlockManager::nodeMapExtraSpacePerFace() * (1 + 2 * overAllocationFactor);
faceToNodeMap.reserveValues(valuesToReserve);
faceToNodeMap.reserveValues(2 * entriesToReserve); // TODO I don"t undertand anythin about LVARRAY :@

// Append the individual arrays.
for (localIndex faceID = 0; faceID < numUniqueFaces; ++faceID)
{
  faceToNodeMap.appendArray(numNodesPerFace[faceID]);
  faceToNodeMap.setCapacityOfArray(faceToNodeMap.size() - 1,
                                    numNodesPerFace[faceID] + HexCellBlockManager::nodeMapExtraSpacePerFace());
}

// Each face may belong to _maximum_ 2 elements.
// If it belongs to only one, we put `-1` for the undefined value.
faceToElemMap.resize(numUniqueFaces, 2);

// TODO I'm really not sure.
faceToEdgesMap.resize(numUniqueFaces, 2 * HexCellBlockManager::edgeMapExtraSpacePerFace());
}

/**
* @brief Populate the lowestNodeToFaces map.
* @param [in] numNodes Number of nodes
* @param [in] cellBlocks The cell blocks on which we need to operate.
*
* For each face of each element of each cell blocks,
* this function stores some information of the faces (@see NodesAndElementOfFace)
* in a node to faces (information) mapping.
* The key of this mapping is the lowest node index of the face.
* E.g. faces {3, 5, 6, 2} and {4, 2, 9, 7} will both be stored in "bucket" of node 2.
* Also, bucket of faces information are sorted (@see NodesAndElementOfFace) to make specific computations possible.
*/
ArrayOfArrays<NodesAndElementOfFace> createLowestNodeToFaces(localIndex numNodes, const Group &cellBlocks)
{
ArrayOfArrays<NodesAndElementOfFace> lowestNodeToFaces(numNodes, 2 * HexCellBlockManager::maxFacesPerNode());

// The function is a bit simplified and is not run in parallel anymore.
// Can be improved.
for (localIndex iCellBlock = 0; iCellBlock < cellBlocks.numSubGroups(); ++iCellBlock)
{
  const CellBlock &cb = cellBlocks.getGroup<CellBlock>(iCellBlock);
  localIndex const numFacesPerElement = cb.numFacesPerElement();
  localIndex const numElements = cb.numElements();

  for (localIndex iElement = 0; iElement < numElements; ++iElement)
  {
    // Looping on the faces of the cell
    for (localIndex iFace = 0; iFace < numFacesPerElement; ++iFace)
    {
      // Get all the nodes of the face
      array1d<localIndex> nodesInFace;
      cb.getFaceNodes(iElement, iFace, nodesInFace);
      // Fill the result with the collected data.
      localIndex const &lowestNode = *std::min_element(nodesInFace.begin(), nodesInFace.end());
      lowestNodeToFaces.emplaceBack(lowestNode, nodesInFace, iElement, iCellBlock, iFace);
    }
  }
}

// Loop over all the nodes and sort the associated faces.
forAll<parallelHostPolicy>(numNodes, [&](localIndex const nodeID)
                            {
NodesAndElementOfFace * const faces = lowestNodeToFaces[ nodeID ];
std::sort( faces, faces + lowestNodeToFaces.sizeOfArray( nodeID ) ); });

return lowestNodeToFaces;
}

/**
* @brief Filling the elements to faces maps in the cell blocks.
* @param lowestNodeToFaces The lowest node to faces information array.
* @param uniqueFaceOffsets The unique face offsets.
* @param cellBlocks The cell blocks for which we need to compute the element to faces mappings.
*
* @note @p lowestNodeToFaces and @p uniqueFaceOffsets are better described in the documentations of the functions that build them.
*/
void fillElementToFacesOfCellBlocks(ArrayOfArrays<NodesAndElementOfFace> const &lowestNodeToFaces,
                                  array1d<localIndex> const &uniqueFaceOffsets,
                                  Group &cellBlocks)
{
localIndex const numNodes = lowestNodeToFaces.size();
for (localIndex nodeID = 0; nodeID < numNodes; ++nodeID)
{
  localIndex const numFaces = lowestNodeToFaces.sizeOfArray(nodeID);
  for (localIndex j = 0, curFaceID = uniqueFaceOffsets[nodeID]; j < numFaces; ++j, ++curFaceID)
  {
    const NodesAndElementOfFace &f0 = lowestNodeToFaces(nodeID, j);
    CellBlock &cb0 = cellBlocks.getGroup<CellBlock>(f0.iCellBlock);
    cb0.setElementToFaces(f0.element, f0.iFace, curFaceID);

    // If the following face exists and is identical to the current one,
    // Then we insert the following face with the same face id
    // and remove it from the next iteration (it's already inserted).
    if (j < numFaces - 1)
    {
      const NodesAndElementOfFace &f1 = lowestNodeToFaces(nodeID, j + 1);
      if (f0 == f1)
      {
        CellBlock &cb1 = cellBlocks.getGroup<CellBlock>(f1.iCellBlock);
        cb1.setElementToFaces(f1.element, f1.iFace, curFaceID);
        ++j;
      }
    }
  }
}
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
m_numNodes = numNodes;
m_nodesPositions.resize(m_numNodes);
m_nodeLocalToGlobal.resize(m_numNodes);
m_nodeLocalToGlobal.setValues<serialPolicy>(-1);
}

} // namespace geosx
