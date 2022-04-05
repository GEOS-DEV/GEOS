/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "NullModel.hpp"

namespace geosx
{
namespace constitutive
{

NullModel::NullModel( string const & name,
                      Group * const parent ):
  ConstitutiveBase( name, parent )
{}

NullModel::~NullModel()
{}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, NullModel, string const &, dataRepository::Group * const )

} // constitutive
} /* namespace geosx */
