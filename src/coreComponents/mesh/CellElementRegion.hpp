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
 * @file ElementRegionBase.hpp
 *
 */

#ifndef CORECOMPONENTS_MESH_CELLELEMENTREGION_HPP_
#define CORECOMPONENTS_MESH_CELLELEMENTREGION_HPP_

#include "ElementRegionBase.hpp"

namespace geosx
{
class EdgeManager;
class EmbeddedSurfaceGenerator;

/**
 * @class ElementRegionBase
 *
 * The FaceElementRegion class contains the functionality to support the concept of a FaceElementRegion in the element
 * hierarchy. FaceElementRegion derives from ElementRegion and has an entry in the ObjectManagerBase catalog.
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

  virtual void GenerateMesh( Group const * const cellBlocks ) override;

  void GenerateAggregates( FaceManager const * const faceManager, NodeManager const * const NodeManager );

  void GenerateEmbeddedSurfaces( FaceManager              const * const faceManager,
                                 EdgeManager              const * const edgeManager,
                                 NodeManager              const * const nodeManager,
                                 EmbeddedSurfaceGenerator const * const embeddedSurface);

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

#endif /* CORECOMPONENTS_MESH_CELLELEMENTREGION_HPP_ */
