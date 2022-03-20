
#include "dataRepository/Group.hpp"

#include <conduit.hpp>
#include <conduit_relay_io.hpp>

#include <gtest/gtest.h>

using namespace geosx::dataRepository;

void testConduit( bool const removeOne )
{
  std::string const filename = GEOSX_FMT( "conduit-{}.hdf5", removeOne );

  // write data
  {
    conduit::Node root;
    conduit::Node & group = root["group"];
    for( int i = 0; i < 3; ++i )
    {
      conduit::Node & subgroup = group[std::to_string( i )];
      subgroup["data"] = i;
    }
    if( removeOne ) group.remove( "0" );
    conduit::relay::io::save( root, filename, "hdf5" );
  }

  // read data
  {
    conduit::Node root;
    conduit::relay::io::load( filename, "hdf5", root );
    conduit::Node const & group = root["group"];
    for( int i = 1; i < 3; ++i )
    {
      conduit::Node const & subgroup = group[ std::to_string( i ) ];
      int const value = subgroup["data"].value();
      EXPECT_EQ( value, i );
    }
  }
}

TEST( conduit, noRemoval )
{
  testConduit( false );
}

TEST( conduit, withRemoval )
{
  testConduit( true );
}

void testDataRepository( bool const removeOne )
{
  std::string const filename = GEOSX_FMT( "dataRepository-{}.hdf5", removeOne );

  // write data
  {
    conduit::Node root;
    Group group( "group", root );
    for( int i = 0; i < 3; ++i )
    {
      Group & subgroup = group.registerGroup( std::to_string( i ) );
      subgroup.registerWrapper< int >( "data" )
        .setSizedFromParent( false )
        .setRestartFlags( RestartFlags::WRITE_AND_READ )
        .reference() = i;
    }
    if( removeOne ) group.deregisterGroup( "0" );
    group.prepareToWrite();
    conduit::relay::io::save( root, filename, "hdf5" );
    group.finishWriting();
  }

  // read data
  {
    conduit::Node root;
    conduit::relay::io::load( filename, "hdf5", root );
    Group group( "group", root );
    for( int i = 1; i < 3; ++i )
    {
      Group & subgroup = group.registerGroup( std::to_string( i ) );
      subgroup.registerWrapper< int >( "data" )
        .setSizedFromParent( false )
        .setRestartFlags( RestartFlags::WRITE_AND_READ );
    }
    group.loadFromConduit();
    for( int i = 1; i < 3; ++i )
    {
      Group const & subgroup = group.getGroup( std::to_string( i ) );
      int const value = subgroup.getReference< int >( "data" );
      EXPECT_EQ( value, i );
    }
  }
}

TEST( dataRepository, noRemoval )
{
  testDataRepository( false );
}

TEST( dataRepository, withRemoval )
{
  testDataRepository( true );
}
