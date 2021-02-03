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



#include <iostream>

#include "common/Logger.hpp"
#include "Base.hpp"
#include "Derived1.hpp"

int main( int argc, char *argv[] )
{
  GEOSX_LOG( "EXECUTING MAIN" );
  int junk = 1;
  double junk2 = 3.14;
  double junk3 = 2*3.14;
  Parameter param;


  GEOSX_LOG( "Attempting to create a Derived1 object" );
  std::unique_ptr<Base> derived1 = Base::CatalogInterface::Factory( "derived1", junk, junk2, param );
  GEOSX_LOG( "Attempting to create a Derived2 object" );
  std::unique_ptr<Base> derived2 = Base::CatalogInterface::Factory( "derived2", junk, junk3, param );

  Base::CatalogInterface::catalog_cast<Derived1>( *(derived2.get()));
  GEOSX_LOG( "EXITING MAIN" );
}
