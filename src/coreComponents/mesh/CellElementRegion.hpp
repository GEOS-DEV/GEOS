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
 * @file CellElementRegion.hpp
 */

#ifndef GEOS_MESH_CELLELEMENTREGION_HPP_
#define GEOS_MESH_CELLELEMENTREGION_HPP_

#include "ElementRegionBase.hpp"

namespace geos
{

/**
 * @class CellElementRegion
 *
 * The CellElementRegion class contains the functionality to support the concept of a CellElementRegion in the element
 * hierarchy. CellElementRegion derives from ElementRegionBase and has an entry in the ObjectManagerBase catalog.
 */
class CellElementRegion : public ElementRegionBase
{
public:

  /**
   * @name Constructor / Destructor
   */
  ///@{


  /**
   * @brief Constructor.
   * @param name the name of the object in the data hierarchy.
   * @param parent a pointer to the parent group in the data hierarchy.
   */
  CellElementRegion( string const & name, Group * const parent );

  /**
   * @brief Deleted default constructor.
   */
  CellElementRegion() = delete;

  /**
   * @brief Destructor.
   */
  virtual ~CellElementRegion() override;

  /**
   * @name Static factory catalog functions
   */
  ///@{

  /**
   * @brief The key name for the FaceElementRegion in the object catalog.
   * @return A string containing the key name.
   */
  static string catalogName()
  { return "CellElementRegion"; }

  /**
   * @copydoc catalogName()
   */
  virtual string getCatalogName() const override final
  { return catalogName(); }

  ///@}

  /**
   * @name Generation of the cell element mesh region
   */
  ///@{


  /**
   * @brief Getter for m_cellBlockNames
   * @return The array of cell block names.
   */
  arrayView1d< string const > getCellBlockNames() const
  {
    return m_cellBlockNames.toViewConst();
  }

  /**
   * @brief Add a cellBlockRegion name to the list.
   * @param cellBlockName string containing the cell block region name.
   */
  void addCellBlockName( string const & cellBlockName )
  {
    m_cellBlockNames.emplace_back( cellBlockName );
  }

  /**
   * @brief Add an array cellBlockRegion name to the list.
   * @param cellBlockNames array of string containing the cell block region names.
   */
  void addCellBlockNames( arrayView1d< string const > const & cellBlockNames )
  {
    for( auto const & name: cellBlockNames )
    {
      m_cellBlockNames.emplace_back( name );
    }
  }

  /**
   * @brief register every entry of m_cellBlockNames that is in the provided list.
   * After the call, the content of m_cellBlockNames no longer contain duplicates, and takes into account
   * m_cellBlockAttributeValues and m_cellBlockMatchPatterns.
   * @param cellBlocks the list of available cellBlocks.
   */
  virtual void generateMesh( Group const & cellBlocks ) override;

  ///@}

  /**
   * @brief A struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : public ElementRegionBase::viewKeyStruct
  {
    /// @return String key for the coarsening ratio
    static constexpr char const * coarseningRatioString() {return "coarseningRatio"; }

    /// @return String key for the cell block names
    static constexpr char const * sourceCellBlockNamesString() {return "cellBlocks"; }

    /// @return String key for the cell block names
    static constexpr char const * cellBlockAttributeValuesString() {return "cellBlockAttributeValues"; }

    /// @return String key for the cell block matches
    static constexpr char const * cellBlockMatchPatternsString() {return "cellBlocksMatch"; }
  };

  /**
   * @param cellBlockName the cellBlock name
   * @return if we are in the form "1_tetrahedra", we return the region attribute value.
   * Elsewise, we return an empty string.
   */
  static string getCellBlockAttributeValue( string_view cellBlockName );

private:

  integer_array m_cellBlockAttributeValues;

  string_array m_cellBlockMatchPatterns;

  // Cell block names.
  string_array m_cellBlockNames;

  // Coarsening ratio
  real64 m_coarseningRatio;


  /**
   * @brief Register a CellElementSubRegion in the element sub region group (named after viewKeyStruct::elementSubRegions()).
   * The sub-regions are a set of the region that groups the same primitives.
   * @param cellBlock The cellBlock that the sub region will be made of.
   */
  void registerSubRegion( CellBlockABC const & cellBlock );

  /**
   * @return all cellBlock names entries from m_cellBlockAttributeValues,
   * m_cellBlockMatchPatterns and m_cellBlockNames.
   * @param cellBlocks the input mesh cellBlock list
   */
  std::set< string > computeSelectedCellBlocks( Group const & cellBlocks ) const;

};

} /* namespace geos */

#endif /* GEOS_MESH_CELLELEMENTREGION_HPP_ */
