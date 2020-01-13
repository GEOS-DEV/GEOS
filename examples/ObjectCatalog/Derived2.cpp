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
#include "Derived2.hpp"


Derived2::Derived2( int junk, double const & junk2, Parameter& param ):
  Base( junk, junk2, param )
{
  GEOSX_LOG( "calling Derived2 constructor with arguments ("<<junk<<" "<<junk2<<")" );
}

Derived2::~Derived2()
{
  GEOSX_LOG( "calling Derived2 destructor" );
}

REGISTER_CATALOG_ENTRY( Base, Derived2, int, double const &, Parameter& )
