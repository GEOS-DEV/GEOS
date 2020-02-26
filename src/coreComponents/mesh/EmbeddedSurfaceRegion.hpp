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
 * The EmbeddedSurfaceRegion class contains the functionality to support the concept of a EmbeddedSurfaceRegion in the element
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
