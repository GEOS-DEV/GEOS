/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SimpleGeometricObjectBase.cpp
 */

#include "SimpleGeometricObjectBase.hpp"

namespace geosx
{

SimpleGeometricObjectBase::SimpleGeometricObjectBase( string const & name,
                                                      Group * const parent ):
  Group( name, parent )
{
  setInputFlags( dataRepository::InputFlags::OPTIONAL_NONUNIQUE );
}


SimpleGeometricObjectBase::~SimpleGeometricObjectBase()
{}


SimpleGeometricObjectBase::CatalogInterface::CatalogType & SimpleGeometricObjectBase::getCatalog()
{
  static SimpleGeometricObjectBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

} /// namespace geosx
