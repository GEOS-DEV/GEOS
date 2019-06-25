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
 * @file AggregateElementRegion.hpp
 *
 */

#ifndef CORECOMPONENTS_MESH_AGGREGATEELEMENTREGION_HPP_
#define CORECOMPONENTS_MESH_AGGREGATEELEMENTREGION_HPP_

#include "ElementRegion.hpp"

namespace geosx
{
class AggregateElementRegion : public ElementRegion
{
public:
  /**
   * @brief constructor
   * @param name The name of the object in the data hierarchy.
   * @param parent Pointer to the parent group in the data hierarchy.
   */
  AggregateElementRegion( string const & name, ManagedGroup * const parent );

  AggregateElementRegion() = delete;
  virtual ~AggregateElementRegion() override;

  /**
   * @brief The key name for the AggregateElementRegion in the object catalog.
   * @return A string containing the key name.
   */
  static const string CatalogName()
  {
    return "AggregateElementRegion";
  }

  virtual const string getCatalogName() const override final
  {
    return AggregateElementRegion::CatalogName();
  }

  virtual void GenerateMesh( ManagedGroup const * ) override {}

  void Generate( FaceManager const * const faceManager,
                 ElementRegion const * const elementRegion,
                 real64 coarseningRatio );

  struct viewKeyStruct : public ElementRegion::viewKeyStruct
  {
  };
};
}

#endif
