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

#ifndef QUADRATURE_H
#define QUADRATURE_H

/**
 * @file Quadrature.h
 * @author white230
 */

#include "dataRepository/ManagedGroup.hpp"
#include "common/DataTypes.hpp"
#include "ObjectCatalog.hpp"

#include <cassert>
#include "fileIO/xmlWrapper.hpp"

/*
 * Pure virtual base class representing a generic quadrature object.
 */
namespace geosx
{

class QuadratureBase : public dataRepository::ManagedGroup
{
public:
  /// Main constructor
  explicit QuadratureBase( std::string const & name,
                           ManagedGroup * const parent );

  /// Destructor
  virtual ~QuadratureBase() override;

  // Catalog name interface
  static string CatalogName() { return "QuadratureBase"; }

  QuadratureBase() = default;
  QuadratureBase( QuadratureBase const & ) = default;
  QuadratureBase( QuadratureBase &&) = default;
  QuadratureBase& operator=( QuadratureBase const & ) = default;
  QuadratureBase& operator=( QuadratureBase&& ) = default;  

  using CatalogInterface = cxx_utilities::CatalogInterface< QuadratureBase, std::string const &, ManagedGroup * const >;
  static CatalogInterface::CatalogType& GetCatalog();

  virtual int size() const = 0;
  virtual R1Tensor integration_point( const int index ) const = 0;
  virtual double integration_weight( const int index ) const = 0;

private:

};
}

#endif
