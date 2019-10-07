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

#include <gtest/gtest.h>

#define OBJECTCATALOGVERBOSE 2
#include "dataRepository/ObjectCatalog.hpp"
#include "Logger.hpp"
#include <sys/time.h>
#include <stdint.h>
#include <string>

using namespace cxx_utilities;
using namespace geosx;
class Base
{
public:
  Base( int & junk, double const & junk2 )
  {
    GEOS_LOG( "calling Base constructor with arguments ("<<junk<<" "<<junk2<<")" );
  }

  virtual ~Base()
  {
    GEOS_LOG( "calling Base destructor" );
  }

  using CatalogInterface = dataRepository::CatalogInterface< Base, int &, double const & >;
  static CatalogInterface::CatalogType & GetCatalog()
  {
    static CatalogInterface::CatalogType catalog;
    return catalog;
  }

  virtual std::string getName() = 0;
};

class Derived1 : public Base
{
public:
  Derived1( int & junk, double const & junk2 ):
    Base( junk, junk2 )
  {
    GEOS_LOG( "calling Derived1 constructor with arguments ("<<junk<<" "<<junk2<<")" );
  }

  ~Derived1()
  {
    GEOS_LOG( "calling Derived1 destructor" );
  }
  static std::string CatalogName() { return "derived1"; }
  std::string getName() { return CatalogName(); }

};
REGISTER_CATALOG_ENTRY( Base, Derived1, int &, double const & )

class Derived2 : public Base
{
public:
  Derived2( int & junk, double const & junk2 ):
    Base( junk, junk2 )
  {
    GEOS_LOG( "calling Derived2 constructor with arguments ("<<junk<<" "<<junk2<<")" );
  }

  ~Derived2()
  {
    GEOS_LOG( "calling Derived2 destructor" );
  }
  static std::string CatalogName() { return "derived2"; }
  std::string getName() { return CatalogName(); }

};
REGISTER_CATALOG_ENTRY( Base, Derived2, int &, double const & )


TEST( testObjectCatalog, testRegistration )
{
  GEOS_LOG( "EXECUTING MAIN" );
  int junk = 1;
  double junk2 = 3.14;
  std::unique_ptr< Base > derived1 = Base::CatalogInterface::Factory( "derived1", junk, junk2 );
  std::unique_ptr< Base > derived2 = Base::CatalogInterface::Factory( "derived2", junk, junk2 );

  assert( derived1->getName() == Derived1::CatalogName() );
  assert( derived2->getName() == Derived2::CatalogName() );
  GEOS_LOG( "EXITING MAIN" );
}
