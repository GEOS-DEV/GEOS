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

// Source includes
#include "dataRepository/Group.hpp"
#include "dataRepository/Wrapper.hpp"

// TPL includes
#include <gtest/gtest.h>
#include <conduit.hpp>

using namespace geos;
using namespace dataRepository;

template< typename T >
class WrapperSetGet : public ::testing::Test
{
public:
  WrapperSetGet():
    m_node(),
    m_group( "root", m_node ),
    m_wrapper( "wrapper", m_group ),
    m_wrapperBase( m_wrapper )
  {}

  void testSizedFromParent( int const value )
  {
    {
      Wrapper< T > & rval = m_wrapper.setSizedFromParent( value );
      EXPECT_EQ( value, m_wrapper.sizedFromParent() );
      EXPECT_EQ( &rval, &m_wrapper );
    }

    {
      WrapperBase & rval = m_wrapperBase.setSizedFromParent( value );
      EXPECT_EQ( value, m_wrapperBase.sizedFromParent() );
      EXPECT_EQ( &rval, &m_wrapperBase );
    }
  }

  void testRestartFlags( RestartFlags const value )
  {
    {
      Wrapper< T > & rval = m_wrapper.setRestartFlags( value );
      EXPECT_EQ( value, m_wrapper.getRestartFlags() );
      EXPECT_EQ( &rval, &m_wrapper );
    }

    {
      WrapperBase & rval = m_wrapperBase.setRestartFlags( value );
      EXPECT_EQ( value, m_wrapperBase.getRestartFlags() );
      EXPECT_EQ( &rval, &m_wrapperBase );
    }
  }

  void testPlotLevel( PlotLevel const value )
  {
    {
      Wrapper< T > & rval = m_wrapper.setPlotLevel( value );
      EXPECT_EQ( value, m_wrapper.getPlotLevel() );
      EXPECT_EQ( &rval, &m_wrapper );
    }

    {
      WrapperBase & rval = m_wrapperBase.setPlotLevel( value );
      EXPECT_EQ( value, m_wrapperBase.getPlotLevel() );
      EXPECT_EQ( &rval, &m_wrapperBase );
    }
  }

  void testInputFlags( InputFlags const value )
  {
    {
      Wrapper< T > & rval = m_wrapper.setInputFlag( value );
      EXPECT_EQ( value, m_wrapper.getInputFlag() );
      EXPECT_EQ( &rval, &m_wrapper );
    }

    {
      WrapperBase & rval = m_wrapperBase.setInputFlag( value );
      EXPECT_EQ( value, m_wrapperBase.getInputFlag() );
      EXPECT_EQ( &rval, &m_wrapperBase );
    }
  }

  void testDescription( string const & value )
  {
    {
      Wrapper< T > & rval = m_wrapper.setDescription( value );
      EXPECT_EQ( value, m_wrapper.getDescription() );
      EXPECT_EQ( &rval, &m_wrapper );
    }

    {
      WrapperBase & rval = m_wrapperBase.setDescription( value );
      EXPECT_EQ( value, m_wrapperBase.getDescription() );
      EXPECT_EQ( &rval, &m_wrapperBase );
    }
  }

private:
  conduit::Node m_node;
  Group m_group;
  Wrapper< T > m_wrapper;
  WrapperBase & m_wrapperBase;
};

using WrapperSetGetTypes = ::testing::Types< int, array1d< real64 >, void *, std::function< void (void) > >;

TYPED_TEST_SUITE( WrapperSetGet, WrapperSetGetTypes, );

TYPED_TEST( WrapperSetGet, SizedFromParent )
{
  this->testSizedFromParent( true );
  this->testSizedFromParent( false );
}

TYPED_TEST( WrapperSetGet, RestartFlags )
{
  this->testRestartFlags( RestartFlags::NO_WRITE );
  this->testRestartFlags( RestartFlags::WRITE_AND_READ );
}

TYPED_TEST( WrapperSetGet, PlotLevel )
{
  this->testPlotLevel( PlotLevel::LEVEL_0 );
  this->testPlotLevel( PlotLevel::LEVEL_1 );
}

TYPED_TEST( WrapperSetGet, InputFlag )
{
  this->testInputFlags( InputFlags::OPTIONAL );
  this->testInputFlags( InputFlags::REQUIRED );
}

TYPED_TEST( WrapperSetGet, Description )
{
  this->testDescription( "First description." );
  this->testDescription( "Second description." );
}

class WrapperLimits : public ::testing::Test
{
protected:
  WrapperLimits():
    m_node(),
    m_group( "root", m_node )
  {}

  template< typename T >
  Wrapper< T > & makeWrapper( string const & name )
  {
    return m_group.template registerWrapper< T >( name );
  }

  conduit::Node m_node;
  Group m_group;
};

TEST_F( WrapperLimits, IsLimitableTrait )
{
  static_assert( is_limitable_v< integer >, "integer must be limitable" );
  static_assert( is_limitable_v< real64 >, "real64 must be limitable" );
  static_assert( is_limitable_v< array1d< integer > >, "array1d< integer > must be limitable" );
  static_assert( is_limitable_v< array2d< real64 > >, "array2d< real64 > must be limitable" );
  static_assert( is_limitable_v< array3d< integer > >, "array3d< integer > must be limitable" );

  static_assert( std::is_same< limit_value_type_t< real64 >, real64 >::value, "" );
  static_assert( std::is_same< limit_value_type_t< array1d< real64 > >, real64 >::value, "" );
  static_assert( std::is_same< limit_value_type_t< array2d< integer > >, integer >::value, "" );
}

TEST_F( WrapperLimits, ScalarSetGet )
{
  auto & w = makeWrapper< real64 >( "scalar" );
  w.setLimits( inclusive( 0.0 ), exclusive( 1.0 ), LimitsMode::Error );

  ASSERT_TRUE( w.getMinValue().has_value() );
  ASSERT_TRUE( w.getMaxValue().has_value() );
  EXPECT_DOUBLE_EQ( w.getMinValue()->value, 0.0 );
  EXPECT_TRUE( w.getMinValue()->isInclusive );
  EXPECT_DOUBLE_EQ( w.getMaxValue()->value, 1.0 );
  EXPECT_FALSE( w.getMaxValue()->isInclusive );
  EXPECT_EQ( w.getLimitsMode(), LimitsMode::Error );
}

TEST_F( WrapperLimits, Array1dSetGet )
{
  auto & w = makeWrapper< array1d< real64 > >( "array1d" );
  w.setLimits( 0.0, 1.0, LimitsMode::Error );

  ASSERT_TRUE( w.getMinValue().has_value() );
  ASSERT_TRUE( w.getMaxValue().has_value() );
  EXPECT_DOUBLE_EQ( w.getMinValue()->value, 0.0 );
  EXPECT_DOUBLE_EQ( w.getMaxValue()->value, 1.0 );
  EXPECT_EQ( w.getLimitsMode(), LimitsMode::Error );
}

TEST_F( WrapperLimits, ScalarValidateInRange )
{
  auto & w = makeWrapper< real64 >( "scalar" );
  w.setLimits( 0.0, 1.0, LimitsMode::Error );
  w.reference() = 0.5;
  EXPECT_NO_THROW( w.validateLimits() );
}

TEST_F( WrapperLimits, ScalarValidateBelowMinThrows )
{
  auto & w = makeWrapper< real64 >( "scalar" );
  w.setLimits( 0.0, 1.0, LimitsMode::Error );
  w.reference() = -0.1;
  EXPECT_THROW( w.validateLimits(), InputError );
}

TEST_F( WrapperLimits, ScalarValidateAboveMaxThrows )
{
  auto & w = makeWrapper< real64 >( "scalar" );
  w.setLimits( 0.0, 1.0, LimitsMode::Error );
  w.reference() = 1.1;
  EXPECT_THROW( w.validateLimits(), InputError );
}

TEST_F( WrapperLimits, ScalarValidateInclusiveBoundary )
{
  auto & w = makeWrapper< real64 >( "scalar" );
  w.setLimits( inclusive( 0.0 ), inclusive( 1.0 ), LimitsMode::Error );

  w.reference() = 0.0;
  EXPECT_NO_THROW( w.validateLimits() );

  w.reference() = 1.0;
  EXPECT_NO_THROW( w.validateLimits() );
}

TEST_F( WrapperLimits, ScalarValidateExclusiveBoundary )
{
  auto & w = makeWrapper< real64 >( "scalar" );
  w.setLimits( exclusive( 0.0 ), exclusive( 1.0 ), LimitsMode::Error );

  w.reference() = 0.0;
  EXPECT_THROW( w.validateLimits(), InputError );

  w.reference() = 0.5;
  EXPECT_NO_THROW( w.validateLimits() );

  w.reference() = 1.0;
  EXPECT_THROW( w.validateLimits(), InputError );
}

TEST_F( WrapperLimits, Array1dValidateAllInRange )
{
  auto & w = makeWrapper< array1d< real64 > >( "array1d" );
  w.setLimits( 0.0, 1.0, LimitsMode::Error );

  array1d< real64 > & data = w.reference();
  data.resize( 4 );
  data[ 0 ] = 0.0;
  data[ 1 ] = 0.25;
  data[ 2 ] = 0.75;
  data[ 3 ] = 1.0;

  EXPECT_NO_THROW( w.validateLimits() );
}

TEST_F( WrapperLimits, Array1dValidateOutOfRange )
{
  auto & w = makeWrapper< array1d< real64 > >( "array1d" );
  w.setLimits( 0.0, 1.0, LimitsMode::Error );

  array1d< real64 > & data = w.reference();
  data.resize( 4 );
  data[ 0 ] = 0.5;
  data[ 1 ] = 0.5;
  data[ 2 ] = 42.0;
  data[ 3 ] = 0.5;

  EXPECT_THROW( w.validateLimits(), InputError );
}

TEST_F( WrapperLimits, Array1dValidateEmpty )
{
  auto & w = makeWrapper< array1d< real64 > >( "array1d" );
  w.setLimits( 0.0, 1.0, LimitsMode::Error );

  EXPECT_NO_THROW( w.validateLimits() );
}

TEST_F( WrapperLimits, Array2dValidateAllInRange )
{
  auto & w = makeWrapper< array2d< real64 > >( "array2d" );
  w.setLimits( 0.0, 1.0, LimitsMode::Error );

  array2d< real64 > & data = w.reference();
  data.resize( 2, 3 );
  data( 0, 0 ) = 0.1;
  data( 0, 1 ) = 0.2;
  data( 0, 2 ) = 0.4;
  data( 1, 0 ) = 0.6;
  data( 1, 1 ) = 0.8;
  data( 1, 2 ) = 0.9;

  EXPECT_NO_THROW( w.validateLimits() );
}

TEST_F( WrapperLimits, Array2dValidateOutOfRange )
{
  auto & w = makeWrapper< array2d< real64 > >( "array2d" );
  w.setLimits( 0.0, 1.0, LimitsMode::Error );

  array2d< real64 > & data = w.reference();
  data.resize( 2, 3 );
  data( 0, 0 ) = 0.5;
  data( 0, 1 ) = 0.5;
  data( 0, 2 ) = 0.5;
  data( 1, 0 ) = 4000.0;
  data( 1, 1 ) = 0.5;
  data( 1, 2 ) = 0.5;

  EXPECT_THROW( w.validateLimits(), InputError );
}

TEST_F( WrapperLimits, Array2dValidateEmpty )
{
  auto & w = makeWrapper< array2d< real64 > >( "array2d" );
  w.setLimits( 0.0, 1.0, LimitsMode::Error );

  EXPECT_NO_THROW( w.validateLimits() );
}
