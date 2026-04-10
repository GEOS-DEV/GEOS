/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "FieldSpecificationABC.hpp"

namespace geos
{
using namespace dataRepository;

FieldSpecificationABC::FieldSpecificationABC( string const & name, Group * parent ):
  Group( name, parent )
{}

FieldSpecificationABC::~FieldSpecificationABC()
{}

FieldSpecificationABC::CatalogInterface::CatalogType &
FieldSpecificationABC::getCatalog()
{
  static FieldSpecificationABC::CatalogInterface::CatalogType catalog;
  return catalog;
}

}
