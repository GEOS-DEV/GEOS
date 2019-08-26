/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file FaceElementRegion.hpp
 *
 */

#ifndef CORECOMPONENTS_MESH_CELLELEMENTREGION_HPP_
#define CORECOMPONENTS_MESH_CELLELEMENTREGION_HPP_

#include "ElementRegionBase.hpp"

namespace geosx
{


/**
 * @class FaceElementRegion
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
  CellElementRegion( string const & name, ManagedGroup * const parent );

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

  virtual void GenerateMesh( ManagedGroup const * const cellBlocks ) override;

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

#endif /* CORECOMPONENTS_MESH_CELLELEMENTREGION_HPP_ */
