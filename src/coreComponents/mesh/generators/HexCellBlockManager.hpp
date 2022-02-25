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

/**
 * @file HexCellBlockManager.hpp
 */

#ifndef GEOSX_MESH_HEXCELLBLOCKMANAGER_H_
#define GEOSX_MESH_HEXCELLBLOCKMANAGER_H_

#include "mesh/generators/CellBlockManagerABC.hpp"
#include "mesh/generators/CellBlock.hpp"

#include "common/DataTypes.hpp"


#include <cassert>

namespace geosx
{


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
struct Hex {
  static const unsigned int nbVertices = 8;
  static const unsigned int nbEdges = 12;
  static const unsigned int nbFacets = 6;
  static const unsigned int nbTriFacets = 0;
  static const unsigned int nbQuadFacets = 6;

  static constexpr HexVertexIndex facetVertex[6][4] = {
    { 0,1,2,3 },{ 4,5,6,7 },{ 0,1,5,4 },
    { 1,2,6,5 },{ 3,2,6,7 },{ 0,3,7,4 } };

  static constexpr HexVertexIndex edgeVertex[12][2]{
    { 0,1 },{ 0,3 },{ 0,4 },{ 1,2 },{ 1,5 },{ 2,3 },
    { 2,6 },{ 3,7 },{ 4,5 },{ 4,7 },{ 5,6 },{ 6,7 } };

  static constexpr HexVertexIndex quadFacetTriangleVertex[24][3] = {
    { 0,1,2 },{ 0,2,3 },{ 0,1,3 },{ 1,2,3 },
    { 4,5,6 },{ 4,6,7 },{ 4,5,7 },{ 5,6,7 },
    { 0,1,5 },{ 0,5,4 },{ 0,1,4 },{ 1,5,4 },
    { 1,2,6 },{ 1,6,5 },{ 1,2,5 },{ 2,6,5 },
    { 3,2,6 },{ 3,6,7 },{ 3,2,7 },{ 2,6,7 },
    { 0,3,7 },{ 0,7,4 },{ 3,7,4 },{ 3,4,0 } };

  static constexpr HexFacetIndex vertexAdjacentFacet[8][3] = {
    { 0,2,5 },{ 0,2,3 },{ 0,3,4 },{ 0,4,5 },
    { 1,2,5 },{ 1,2,3 },{ 1,3,4 },{ 1,4,5 } };

  static constexpr HexVertexIndex vertexAdjacentVertex[8][3] = {
    { 1,3,4 },{ 2,0,5 },{ 3,1,6 },{ 0,2,7 },
    { 7,5,0 },{ 4,6,1 },{ 5,7,2 },{ 6,4,3 } };

  static constexpr HexVertexIndex orientedFacetVertex[6][4] = {
    { 0,1,2,3 },{ 4,7,6,5 },{ 0,4,5,1 },
    { 1,5,6,2 },{ 3,2,6,7 },{ 0,3,7,4 } };
};




/**
 * @class HexCellBlockManager
 * @brief The HexCellBlockManager specializes CellBlockManagerABC for full hexahedral meshes  
 * The hexahedral mesh may be structured or unstructured.
 */
class HexCellBlockManager : public CellBlockManagerABC
{
public:

  /**
   * @brief Constructor for HexCellBlockManager object.
   * @param name name of this instantiation of HexCellBlockManager
   * @param parent pointer to the parent Group of this instantiation of HexCellBlockManager
   */
  HexCellBlockManager( string const & name, Group * const parent );
  HexCellBlockManager( const HexCellBlockManager & ) = delete;
  HexCellBlockManager & operator=( const HexCellBlockManager & ) = delete;
  ~HexCellBlockManager() override = default;

  // ------  Group functions  ---- 
  virtual Group * createChild( string const & childKey, string const & childName ) override;

  using Group::resize;

  /**
   * @brief Set the number of elements for a set of element regions.
   * @param[in] numElements list of the new element numbers
   * @param[in] regionNames list of the element region names
   */
  void resize( integer_array const & numElements,
               string_array const & regionNames );


  // ------ CellBlock management functions

  // TODO Think about putting these at the ABC level ? 
  /**
   * @brief Get cell block by name.
   * @param[in] name Name of the cell block.
   * @return Reference to the cell block instance.
   */
  CellBlock & getCellBlock( string const & name )
  { return this->getGroup( viewKeyStruct::cellBlocks() ).getGroup< CellBlock >( name ); }

  const Group & getCellBlocks() const override
  { return this->getGroup( viewKeyStruct::cellBlocks() ); }

  Group & getCellBlocks() override
  { return this->getGroup( viewKeyStruct::cellBlocks() ); }
  
  /**
   * @brief Registers and returns a cell block of name @p name.
   * @param name The name of the created cell block.
   * @return A reference to the new cell block. The HexCellBlockManager owns this new instance.
   */
  CellBlock & registerCellBlock( string name )
  { return this->getCellBlocks().registerGroup< CellBlock >( name ); }

  /**
   * @brief Launch kernel function over all the sub-regions
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   * 
   * //TODO Renaming: these are CellBlocks and not SubRegions
   */
  template< typename LAMBDA >
  void forElementSubRegions( LAMBDA lambda )
  {
    this->getGroup( viewKeyStruct::cellBlocks() ).forSubGroups< CellBlock >( lambda );
  }

  // ----------- 

  // TODO Get rid of this
  static constexpr int maxFacesPerNode()
  { 
    assert(false);
    return -1; 
  }
  // TODO Get rid of this 
  static constexpr localIndex getFaceMapOverallocation()
  {
    assert(false);
    return -1; 
  }

  localIndex numNodes() const override
  { return m_numNodes; }

  localIndex numEdges() const override
  { return m_numEdges; }

  localIndex numFaces() const override
  { return m_numFaces; }


  /**
   * @brief Compute all the mappings
   */
  void buildMaps();

  bool computeFacesToNodes()    { return false; }
  bool computeEdgesToNodes()    { return false; }

  bool computeNodesToEdges();
  bool computeNodesToFaces(); 
  bool computeNodesToElements() { return false; }

  bool computeFacesToElements() { return false; }
  bool computeFacesToEdges()    { return false; }
  
  bool computeEdgesToFaces()    { return false; }
  bool computeEdgesToElements() { return false; }  // A priori used by no one

  
  array2d< geosx::localIndex >     getEdgeToNodes() const override { return m_edgeToNodes; }
  ArrayOfSets< geosx::localIndex > getEdgeToFaces() const override { return m_edgeToFaces; }


  ArrayOfArrays< localIndex >        getFaceToNodes() const override    { return m_faceToNodes; }
  ArrayOfArrays< geosx::localIndex > getFaceToEdges() const override    { return m_faceToEdges; }
  array2d< localIndex >              getFaceToElements() const override { return m_faceToElements; }

  
  ArrayOfSets< localIndex >   getNodeToEdges() const override { return m_nodeToEdges; }
  ArrayOfSets< localIndex >   getNodeToFaces() const override { return m_nodeToFaces; }
  ArrayOfArrays< localIndex > getNodeToElements() const override;

  /**
   * @brief Defines the number of nodes and resizes some underlying arrays appropriately.
   * @param[in] numNodes The number of nodes.
   *
   * The nodes coordinates and nodes local to global mappings get resized to @p numNodes.
   */
  // TODO Improve doc. Is it per domain, are there duplicated nodes because of subregions?
  void setNumNodes( localIndex numNodes ); 


  // TODO What is the point of implementing the abstract accessor
  // and having the same name for functions giving open-access 

  /**
   * @brief Returns a view to the vector holding the nodes coordinates
   * @return The reference
   */
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > getNodesPositions()
  { return m_nodesPositions.toView(); }

  array2d< real64, nodes::REFERENCE_POSITION_PERM > getNodesPositions() const override
  { return m_nodesPositions; }

  array1d< globalIndex > getNodeLocalToGlobal() const override
  { return m_nodeLocalToGlobal; }
  
  arrayView1d< globalIndex > getNodeLocalToGlobal()
  { return m_nodeLocalToGlobal.toView(); }
  
  /**
   * @brief Returns a mutable reference to the node sets.
   * @return A reference to the mapping.
   */
  std::map< string, SortedArray< localIndex > > & getNodeSets()
  { return m_nodeSets; }
  
  std::map< string, SortedArray< localIndex > > const & getNodeSets() const override
  { return m_nodeSets; }


private:

  struct viewKeyStruct
  {
    /// Cell blocks key
    static constexpr char const * cellBlocks() { return "cellBlocks"; }
  };

  /**
   * @brief Get cell block at index @p iCellBlock.
   * @param[in] iCellBlock The cell block index.
   */
  const CellBlockABC & getCellBlock( localIndex iCellBlock ) const
  { return this->getCellBlocks().getGroup< const CellBlockABC >( iCellBlock ); }


  /**
   * @brief Returns the number of cells blocks
   */
  localIndex numCellBlocks() const
  { return this->getCellBlocks().numSubGroups(); }


  /**
   * @brief Trigger the face to nodes, edges and elements mappings.
   */
  void buildFaceMaps();

private:
  // The numbers of things we are dealing with 
  localIndex m_numNodes = 0;
  localIndex m_numEdges = 0;
  localIndex m_numFaces = 0;
  localIndex m_numElements = 0;

  // The nodes that this HexCellBlockManager is responsible of
  // with their global index in the full mesh
  // Typically one CellBlockManager per MPI rank
  array2d< real64, nodes::REFERENCE_POSITION_PERM > m_nodesPositions;
  array1d< globalIndex >  m_nodeLocalToGlobal;

  // The cells/Elements are stored by the CellBlocks

  // The mappings that this class is responsible of building
  // TODO Why are storage strategies different is beyond my understanding

  ArrayOfSets< localIndex >   m_nodeToEdges;
  ArrayOfSets< localIndex >   m_nodeToFaces;
  ArrayOfArrays< localIndex > m_nodeToElements;

  ArrayOfSets< localIndex > m_edgeToFaces;
  array2d< localIndex >     m_edgeToNodes;

  ArrayOfArrays< localIndex > m_faceToNodes;
  ArrayOfArrays< localIndex > m_faceToEdges;
  array2d< localIndex >       m_faceToElements;

  // TODO  What are these ? 
  std::map< string, SortedArray< localIndex > > m_nodeSets;

};

}
#endif /* GEOSX_MESH_CELLBLOCKMANAGER_H_ */
