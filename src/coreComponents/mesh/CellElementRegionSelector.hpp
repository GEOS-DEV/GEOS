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
 * @file CellElementRegionSelector.hpp
 */

#ifndef GEOS_MESH_CELLELEMENTREGIONSELECTOR_HPP_
#define GEOS_MESH_CELLELEMENTREGIONSELECTOR_HPP_

#include "mesh/CellElementRegion.hpp"

namespace geos
{

/**
 * @class CellElementRegionSelector
 *
 * CellElementRegionSelector allows a CellElementRegion to safely select cell-blocks according to the user input.
 * It handles all user input checks and throws an exception in the event of inconsistency.
 */
class CellElementRegionSelector
{
public:

  /**
   * @brief Construct a new CellElementRegionSelector.
   * @param cellBlocks a Group containing all the available cell-blocks.
   */
  CellElementRegionSelector( dataRepository::Group const & cellBlocks );

  /**
   * @brief Select the mesh cell-blocks for the specified region following the user inputs.
   * @throw an InputError if the user setting are inconsistant.
   * @param region the region for which we want to select the cell-blocks.
   * @return the selected cell-blocks names.
   */
  std::set< string > selectRegionCellBlocks( CellElementRegion const & region );

  /**
   * @throw An InputError if region cell-blocks selections is inconsistent:
   *        - cellBlock is in more than one region,
   *        - orphan cellBlock
   * @todo For now, multiple regions per cell is not supported (ElementRegionManager::getCellBlockToSubRegionMap()).
   *       We could refactor the CellElementRegion & Mesh classes so regions are mapped to cell-blocks IN the mesh (and potentially
   *       to multiple regions per cell). So, for external meshes, the cell-blocks would no longer be exposed to the final user.
   */
  void checkSelectionConsistency() const;

private:

  /// @brief The Group containing all the available cell-blocks.
  dataRepository::Group const & m_cellBlocks;

  /// @brief A map that link every cell-block name to the CellElementRegion(s) that references it (0 -> n).
  std::map< string, std::vector< CellElementRegion const * > > m_cellBlocksOwners;

  /// @brief A map that link every region attribute values to the CellElementRegion(s) that references it (0 -> n).
  std::map< string, std::vector< CellElementRegion const * > > m_regionAttributeOwners;


  std::set< string > getFNMatchPatterns( CellElementRegion const & region,
                                         std::set< integer > const & requestedAttributeValues,
                                         std::set< string > const & requestedMatchPatterns ) const;

  std::set< string > getFNMatchSelection( CellElementRegion const & region,
                                          std::set< string > const & requestedMatchPatterns ) const;

  std::set< string > getOneByOneSelection( CellElementRegion const & region,
                                           std::set< string > const & requestedCellBlockNames ) const;

  /**
   * @brief Select the specified cell-blocks & region for the specified region.
   * @param region The region for which when want to select cell-blocks
   * @param attributeValues The attribute values we want to select (can be empty).
   * @param cellBlockNames The cell-block names we want to select (can be empty).
   */
  void selectRegionCellBlocks( CellElementRegion const & region,
                               std::set< integer > const & attributeValues,
                               std::set< string > const & cellBlockNames );

};

} /* namespace geos */

#endif /* GEOS_MESH_CELLELEMENTREGIONSELECTOR_HPP_ */
