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
 * @file Basis.h
 * @author white230
 */

#ifndef BASIS_H
#define BASIS_H

#include "dataRepository/ManagedGroup.hpp"
#include "common/DataTypes.hpp"
#include "ObjectCatalog.hpp"

/**
 * Pure virtual base class representing a space
 * of finite element basis functions (i.e the set
 * of all basis functions defined on the parent
 * cell).
 */
namespace geosx
{


class BasisBase : public dataRepository::ManagedGroup
{
public:
  /// Main constructor
  explicit BasisBase( std::string const & name,
                      ManagedGroup * const parent );

  /// Destructor
  virtual ~BasisBase() override;

  // Catalog name interface
  static string CatalogName() { return "BasisBase"; }

  using CatalogInterface = cxx_utilities::CatalogInterface< BasisBase, std::string const &, ManagedGroup * const >;
  static CatalogInterface::CatalogType& GetCatalog();

  virtual int size() const = 0;

  virtual double value( const int index,
                        const R1Tensor &point ) const = 0;

  virtual R1Tensor gradient( const int index,
                             const R1Tensor &point ) const = 0;

  virtual R1Tensor support_point( const int index ) = 0;

};
}
#endif
