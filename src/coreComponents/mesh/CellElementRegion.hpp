/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CellElementRegion.hpp
 *
 */

#ifndef GEOSX_MESH_CELLELEMENTREGION_HPP_
#define GEOSX_MESH_CELLELEMENTREGION_HPP_

#include "ElementRegionBase.hpp"

namespace geosx
{
class EdgeManager;
class EmbeddedSurfaceGenerator;

/**
 * @class CellElementRegion
 *
 * The CellElementRegion class contains the functionality to support the concept of a CellElementRegion in the element
 * hierarchy. CellElementRegion derives from ElementRegionBase and has an entry in the ObjectManagerBase catalog.
 *
 *
 */
class CellElementRegion : public ElementRegionBase
{
public:
  /**
   * @brief constructor
   * @param name The name of the object in the data hierarchy.
   * @param parent Pointer to the parent group in the data hierarchy.
   */
  CellElementRegion( string const & name, Group * const parent );

  CellElementRegion() = delete;
  virtual ~CellElementRegion() override;

  /**
   * @brief The key name for the FaceElementRegion in the object catalog.
   * @return A string containing the key name.
   */
  static const string CatalogName()
  { return "CellElementRegion"; }

  virtual const string getCatalogName() const override final
  { return CellElementRegion::CatalogName(); }


  void AddCellBlockName( string const & cellBlockName )
  {
    m_cellBlockNames.push_back( cellBlockName );
  }

  virtual void GenerateMesh( Group * const cellBlocks ) override;

  void GenerateAggregates( FaceManager const * const faceManager, NodeManager const * const NodeManager );

  struct viewKeyStruct : public ElementRegionBase::viewKeyStruct
  {
    static constexpr auto coarseningRatioString = "coarseningRatio";
    static constexpr auto sourceCellBlockNames = "cellBlocks";
  };


private:
  string_array m_cellBlockNames;
  real64 m_coarseningRatio;

};

} /* namespace geosx */

#endif /* GEOSX_MESH_CELLELEMENTREGION_HPP_ */
