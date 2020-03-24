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
 * @file EmbeddedSurfaceRegion.hpp
 *
 */

#ifndef GEOSX_MESH_EMBEDDEDSURFACEREGION_HPP_
#define GEOSX_MESH_EMBEDDEDSURFACEREGION_HPP_

#include "ElementRegionBase.hpp"

namespace geosx
{

class EdgeManager;

/**
 * @class EmbeddedSurfaceRegion
 *
 * The EmbeddedSurfaceRegion class contains the functionality to support the concept of a EmbeddedSurfaceRegion in the
 * element
 * hierarchy. EmbeddedSurfaceRegion derives from ElementRegion and has an entry in the ObjectManagerBase catalog.
 *
 *
 */
class EmbeddedSurfaceRegion : public ElementRegionBase
{
public:
  /**
   * @brief constructor
   * @param name The name of the object in the data hierarchy.
   * @param parent Pointer to the parent group in the data hierarchy.
   */

  EmbeddedSurfaceRegion( string const & name, Group * const parent );

  EmbeddedSurfaceRegion() = delete;
  virtual ~EmbeddedSurfaceRegion() override;

  /**
   * @brief The key name for the EmbeddedSurfaceRegion in the object catalog.
   * @return A string containing the key name.
   */
  static const string CatalogName()
  { return "EmbeddedSurfaceElementRegion"; }

  virtual const string getCatalogName() const override final
  { return EmbeddedSurfaceRegion::CatalogName(); }

  struct viewKeyStruct : public ElementRegionBase::viewKeyStruct
  {
    static constexpr auto fractureSetString = "fractureSet";
  };


private:

};

} /* namespace geosx */

#endif /* CORECOMPONENTS_MESH_EMBEDDEDSURFACEREGION_HPP_ */
