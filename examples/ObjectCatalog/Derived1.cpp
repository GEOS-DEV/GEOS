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

#include "Derived1.hpp"


Derived1::Derived1( int junk, double const & junk2, Parameter& param ):
  Base( junk, junk2, param )
{
  GEOSX_LOG( "calling Derived1 constructor with arguments ("<<junk<<" "<<junk2<<")" );
}

Derived1::~Derived1()
{
  GEOSX_LOG( "calling Derived1 destructor" );
}

REGISTER_CATALOG_ENTRY( Base, Derived1, int, double const &, Parameter& )
