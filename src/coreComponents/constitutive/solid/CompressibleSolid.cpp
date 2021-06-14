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
 * @file CompressibleSolid.cpp
 */

#include "CompressibleSolid.hpp"
#include "porosity/PressurePorosity.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

template< typename PORO_TYPE >
CompressibleSolid< PORO_TYPE >::CompressibleSolid( string const & name, Group * const parent ):
  CoupledSolid< NullModel, PORO_TYPE >( name, parent )
{}

template< typename PORO_TYPE >
CompressibleSolid< PORO_TYPE >::~CompressibleSolid()
{}

// Register all CoupleSolid model types.
typedef CompressibleSolid< PressurePorosity > CompressibleRock;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, CompressibleRock, string const &, Group * const )

}
} /* namespace geosx */
