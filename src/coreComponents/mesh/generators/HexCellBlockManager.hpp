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
  
class HexMeshConnectivityBuilder;

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
  ~HexCellBlockManager() override;

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

  // TODO Think about putting these at the ABC level? Or we cannot inherit the Group related stuff properly?  
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
   * // TODO Renaming: these are CellBlocks and not SubRegions
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

  localIndex numElements() const
  { return m_numElements; }
  
  array2d< geosx::localIndex >     getEdgeToNodes() override; 
  ArrayOfSets< geosx::localIndex > getEdgeToFaces() override; 

  ArrayOfArrays< localIndex >        getFaceToNodes() override;
  ArrayOfArrays< geosx::localIndex > getFaceToEdges() override;

  // TODO We have a problem - where are the Element index valid ? 
  array2d< localIndex >              getFaceToElements() override;

  ArrayOfSets< localIndex >   getNodeToEdges() override;
  ArrayOfSets< localIndex >   getNodeToFaces() override;
  
  // TODO We have a problem - where are the Element index valid ? 
  ArrayOfArrays< localIndex > getNodeToElements() override;

    /**
   * @brief Compute all possible maps and more
   */
  void buildMaps();

  /**
   * @brief Defines the number of nodes and resizes some underlying arrays appropriately.
   * @param[in] numNodes The number of nodes.
   *
   * The nodes coordinates and nodes local to global mappings get resized to @p numNodes.
   */
  // TODO Improve doc. Is it per domain, are there duplicated nodes because of subregions?
  void setNumNodes( localIndex numNodes ); 


  // TODO What is the point of implementing the abstract accessor
  // and having the same name for functions giving open-bar access 

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

  // TODO Can't these functions be only at the ABC level ? 
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

private:
  HexMeshConnectivityBuilder * m_theOneWhoDoesTheJob;

  // The numbers of things we are dealing with 
  localIndex m_numNodes = 0;
  localIndex m_numEdges = 0;
  localIndex m_numFaces = 0;
  localIndex m_numElements = 0;

  // The nodes that this HexCellBlockManager is responsible of
  // with their global index in the full mesh
  // Typically one CellBlockManager per MPI rank  (if unique MeshLevel and unique MeshBody?)
  array2d< real64, nodes::REFERENCE_POSITION_PERM > m_nodesPositions;

  // TODO Why is this allocated here - initialized to -1
  // but not computed by this class 
  // before the transfer to NodeManager - SubElementRegionManager - 

  // This is filled by the MeshGenerators 
  array1d< globalIndex >  m_nodeLocalToGlobal;

  // TODO  What are these ?
  // Who computes them ?
  std::map<string, SortedArray<localIndex>> m_nodeSets;

};

}
#endif
