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
 * @file BasisBase.hpp
 */

#ifndef BASIS_H
#define BASIS_H

#include "dataRepository/Group.hpp"
#include "common/DataTypes.hpp"
#include "dataRepository/ObjectCatalog.hpp"

/**
 * Pure virtual base class representing a space
 * of finite element basis functions (i.e the set
 * of all basis functions defined on the parent
 * cell).
 */
namespace geosx
{


class BasisBase : public dataRepository::Group
{
public:
  /// Main constructor
  explicit BasisBase( std::string const & name,
                      Group * const parent );

  /// Destructor
  virtual ~BasisBase() override;

  // Catalog name interface
  static string CatalogName() { return "BasisBase"; }

  using CatalogInterface = dataRepository::CatalogInterface< BasisBase, std::string const &, Group * const >;
  static CatalogInterface::CatalogType & GetCatalog();

  virtual int size() const = 0;

  virtual double value( const int index,
                        const R1Tensor & point ) const = 0;

  virtual R1Tensor gradient( const int index,
                             const R1Tensor & point ) const = 0;

  virtual R1Tensor support_point( const int index ) = 0;

};
}
#endif
