/*
* ------------------------------------------------------------------------------------------------------------
* SPDX-License-Identifier: LGPL-2.1-only
*
* Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
* Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
* Copyright (c) 2018-2020 TotalEnergies
* Copyright (c) 2020-     GEOSX Contributors
* All rights reserved
*
* See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
* ------------------------------------------------------------------------------------------------------------
*/

#ifndef GEOSX_CELLBLOCKMANAGERBASE_HPP
#define GEOSX_CELLBLOCKMANAGERBASE_HPP

#include "mesh/generators/CellBlockManagerABC.hpp"
#include "mesh/generators/CellBlock.hpp"

#include <map>

namespace geosx
{

/**
 * @brief Base class for a CellBlockManager for unstructured meshes
 * that manages a set of CellBlocks where the Nodes position and Elements are filled.
 * This class is derived to specialize the CellBlockManager per type of mesh
 * Hex - Tet (not implemented) - Hybrid (not implemented)
 * No overallocation
 */
class CellBlockManagerBase : public CellBlockManagerABC
{
public:
  /**
   * @brief Constructor
   * @param name The name of this Group.
   * @param parent The parent Group.
   */
  CellBlockManagerBase(string const &name, Group *const parent)
      : CellBlockManagerABC(name, parent),
        m_nodesPositions(0, 3)
  {
    this->registerGroup<dataRepository::Group>(viewKeyStruct::cellBlocks());
  }

  virtual ~CellBlockManagerBase() = default;

  // ------  Group functions  ----
  virtual Group *createChild(string const & GEOSX_UNUSED_PARAM( childKey ),
                             string const & GEOSX_UNUSED_PARAM( childName )) override
  {
    return nullptr;
  }

  using Group::resize;

  /**
   * @brief Set the number of elements for a set of element regions.
   * @param[in] numElements list of the new element numbers
   * @param[in] regionNames list of the element region names
   */
  void resize(integer_array const &numElements,
              string_array const &regionNames)
  {
    localIndex const numRegions = LvArray::integerConversion<localIndex>(regionNames.size());
    for (localIndex reg = 0; reg < numRegions; ++reg)
    {
      this->getCellBlock(regionNames[reg]).resize(numElements[reg]);
    }
  }

  // ------ CellBlock management functions

  // TODO Think about putting these at the ABC level?
  // Or can we not inherit the Group related stuff properly?
  /**
   * @brief Get cell block by name.
   * @param[in] name Name of the cell block.
   * @return Reference to the cell block instance.
   */
  CellBlock &getCellBlock(string const &name)
  {
    return this->getGroup(viewKeyStruct::cellBlocks()).getGroup<CellBlock>(name);
  }

  const Group &getCellBlocks() const override
  {
    return this->getGroup(viewKeyStruct::cellBlocks());
  }

  Group &getCellBlocks() override
  {
    return this->getGroup(viewKeyStruct::cellBlocks());
  }

  /**
   * @brief Registers and returns a cell block of name @p name.
   * @param name The name of the created cell block.
   * @return A reference to the new cell block. The HexCellBlockManager owns this new instance.
   */
  CellBlock &registerCellBlock(string name)
  {
    return this->getCellBlocks().registerGroup<CellBlock>(name);
  }

  /**
   * @brief Launch kernel function over all the sub-regions
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   *
   * // TODO Renaming: these are CellBlocks and not SubRegions
   */
  template <typename LAMBDA>
  void forElementSubRegions(LAMBDA lambda)
  {
    this->getGroup(viewKeyStruct::cellBlocks()).forSubGroups<CellBlock>(lambda);
  }

  /**
   * @brief Total number of nodes across all the cell blocks.
   * @return The total number of nodes.
   *
   * Nodes shared by multiple cell blocks are counted only once.
   */
  localIndex numNodes() const override
  {
    return m_numNodes;
  }

  virtual localIndex numEdges() const override = 0;
  virtual localIndex numFaces() const override = 0;

  virtual ArrayOfSets<localIndex> getNodeToEdges() override = 0;
  virtual ArrayOfSets<localIndex> getNodeToFaces() override = 0;
  virtual ArrayOfArrays<localIndex> getNodeToElements() override  = 0;
  virtual array2d<geosx::localIndex> getEdgeToNodes() override  = 0;
  virtual ArrayOfSets<geosx::localIndex> getEdgeToFaces() override = 0;
  virtual ArrayOfArrays<localIndex> getFaceToNodes() override = 0;
  virtual ArrayOfArrays<geosx::localIndex> getFaceToEdges() override = 0;
  virtual array2d<localIndex> getFaceToElements() override = 0;

  /**
   * @brief Defines the number of nodes and resizes some underlying arrays appropriately.
   * @param[in] numNodes The number of nodes.
   *
   * The nodes coordinates and nodes local to global mappings get resized to @p numNodes.
   */
  // TODO Improve doc. Is it per domain, are there duplicated nodes because of subregions?
  void setNumNodes(localIndex numNodes) 
  {
    std::cout << "coucou c'est moi je suis la base  " << std::endl;
    m_numNodes = numNodes;
    m_nodesPositions.resize(numNodes);
    // TODO why is this allocated here and filled somewhere else?
    m_nodeLocalToGlobal.resize(numNodes);
    m_nodeLocalToGlobal.setValues<serialPolicy>(-1);
  }

  // TODO What is the point of implementing the abstract accessor
  // and having the same name for functions giving open-bar access

  /**
   * @brief Returns a view to the vector holding the nodes coordinates
   * @return The reference
   */
  arrayView2d<real64, nodes::REFERENCE_POSITION_USD> getNodesPositions()
  {
    return m_nodesPositions.toView();
  }

  array2d<real64, nodes::REFERENCE_POSITION_PERM> getNodesPositions() const override
  {
    return m_nodesPositions;
  }

  array1d<globalIndex> getNodeLocalToGlobal() const override
  {
    return m_nodeLocalToGlobal;
  }

  arrayView1d<globalIndex> getNodeLocalToGlobal()
  {
    return m_nodeLocalToGlobal.toView();
  }

  /**
   * @brief Returns a mutable reference to the node sets.
   * @return A reference to the mapping.
   */
  std::map<string, SortedArray<localIndex>> &getNodeSets()
  {
    return m_nodeSets;
  }

  std::map<string, SortedArray<localIndex>> const &getNodeSets() const override
  {
    return m_nodeSets;
  }

protected:
  // TODO Can't these functions be only at the ABC level ?
  struct viewKeyStruct
  {
    /// Cell blocks key
    static constexpr char const *cellBlocks() { return "cellBlocks"; }
  };

  /**
   * @brief Get cell block at index @p iCellBlock.
   * @param[in] iCellBlock The cell block index.
   */
  const CellBlockABC &getCellBlock(localIndex iCellBlock) const
  {
    return this->getCellBlocks().getGroup<const CellBlockABC>(iCellBlock);
  }

  /**
   * @brief Returns the number of cells blocks
   */
  localIndex numCellBlocks() const
  {
    return this->getCellBlocks().numSubGroups();
  }

protected:
  localIndex m_numNodes = 0;

  // The nodes that this HexCellBlockManager is responsible of
  // with their global index in the full mesh
  // Typically one CellBlockManager per MPI rank  (if unique MeshLevel and unique MeshBody?)
  array2d<real64, nodes::REFERENCE_POSITION_PERM> m_nodesPositions;

  // TODO Why is this allocated here - initialized to -1
  // but not computed by this class
  // before the transfer to NodeManager - SubElementRegionManager -

  // This is filled by the MeshGenerators
  array1d<globalIndex> m_nodeLocalToGlobal;

  // TODO  What are these ?
  // Who computes them ?
  std::map<string, SortedArray<localIndex>> m_nodeSets;
};

}
#endif // include guard
