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

#ifndef QUADRATURE_H
#define QUADRATURE_H

/**
 * @file Quadrature.h
 */

#include "common/DataTypes.hpp"
#include "dataRepository/Group.hpp"
#include "dataRepository/ObjectCatalog.hpp"
#include "dataRepository/xmlWrapper.hpp"

#include <cassert>

/*
 * Pure virtual base class representing a generic quadrature object.
 */
namespace geosx
{

class QuadratureBase : public dataRepository::Group
{
public:
  /// Main constructor
  explicit QuadratureBase( std::string const & name,
                           Group * const parent );

  /// Destructor
  virtual ~QuadratureBase() override;

  // Catalog name interface
  static string CatalogName() { return "QuadratureBase"; }

  using CatalogInterface = dataRepository::CatalogInterface< QuadratureBase, std::string const &, Group * const >;
  static CatalogInterface::CatalogType & GetCatalog();

  virtual int size() const = 0;
  virtual R1Tensor integration_point( const int index ) const = 0;
  virtual double integration_weight( const int index ) const = 0;

private:

};
}

#endif
