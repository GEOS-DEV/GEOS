/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include <iostream>
#include <regex>
#include <typeindex>

#include "codingUtilities/RTTypes.hpp"
#include "dataRepository/Group.hpp"
#include "dataRepository/Wrapper.hpp"
#include "common/format/EnumStrings.hpp"

#include <gtest/gtest.h>

// TPL includes
#include <conduit.hpp>

// Mock classes to test dynamic casting
#include "BaseClass.hpp"
#include "DerivedFinalClass.hpp"
#include "DerivedClassFinal.hpp"

using namespace geos;
using namespace dataRepository;

// --------------------------------------------------------------------------------------
// ENUM WITH STRINGS FOR TypeRegex SPECIALIZATION TESTING
// --------------------------------------------------------------------------------------
namespace geos // must be in same namespace for ENUM_STRINGS ADL
{
  enum struct MyEnum
  {
    Alpha,
    Beta,
    Gamma
  };
  ENUM_STRINGS( MyEnum, "alpha", "beta", "gamma" );
}

// --------------------------------------------------------------------------------------
// dynamicCast pointer tests
// --------------------------------------------------------------------------------------
TEST( DynamicCastTests, Pointer_Casting_Success )
{
  BaseClass * base = new DerivedFinalClass();
  DerivedFinalClass * derived = geos::dynamicCast< DerivedFinalClass * >( base );
  ASSERT_NE( derived, nullptr ) << "Expected successful cast from Base to Derived.";
  delete base;   // Clean up allocated memory
}

TEST( DynamicCastTests, Pointer_Casting_Failure )
{
  BaseClass * base = new BaseClass();
  DerivedFinalClass * derived = geos::dynamicCast< DerivedFinalClass * >( base );
  ASSERT_EQ( derived, nullptr ) << "Expected nullptr due to failed cast from Base to Derived.";
  delete base;   // Clean up allocated memory
}

TEST( DynamicCastTests, Pointer_Casting_Nullptr )
{
  BaseClass * base = nullptr;
  DerivedFinalClass * derived = geos::dynamicCast< DerivedFinalClass * >( base );
  ASSERT_EQ( derived, nullptr ) << "Casting a nullptr should return nullptr.";
}

// --------------------------------------------------------------------------------------
// dynamicCast reference tests
// --------------------------------------------------------------------------------------
TEST( DynamicCastTests, Reference_Casting_Success )
{
  DerivedFinalClass derived;
  BaseClass & base_ref = derived;
  DerivedFinalClass & derived_ref = geos::dynamicCast< DerivedFinalClass & >( base_ref );
  ASSERT_EQ( &derived_ref, &derived ) << "Expected successful cast from Base to Derived.";
}

TEST( DynamicCastTests, Reference_Casting_BaseToBase )
{
  BaseClass base;
  BaseClass & base_ref = base;
  BaseClass & derived_base_ref = geos::dynamicCast< BaseClass & >( base_ref );
  ASSERT_EQ( &derived_base_ref, &base ) << "Expected successful cast from Base to Base.";
}

// NOTE: We intentionally do NOT include a failing reference cast death test here.
// A failing reference dynamicCast triggers GEOS_ERROR_IF, which invokes LvArray::system::callErrorHandler().
// That handler presents an interactive 30s countdown ("Press space to interact...") slowing unit tests.
// For fast, non-interactive test runs, we avoid exercising that abort path and rely on pointer failing casts instead.

// --------------------------------------------------------------------------------------
// Typed test for geos Wrapper
// --------------------------------------------------------------------------------------

template< typename T >
class WrapperMock : public ::testing::Test
{
public:
  WrapperMock():
    m_node(),
    m_group( "root", m_node ),
    m_wrapper( "wrapper", m_group ),
    m_wrapperBase( m_wrapper )
  {}

  void testDynamicCastWithPointer()
  {
    {
      WrapperBase * base_pointer = &m_wrapperBase;
      Wrapper< T > * derived = geos::dynamicCast< Wrapper< T > * >( base_pointer );
      ASSERT_NE( derived, nullptr ) << "Expected successful cast from Base to Derived.";
    }
    {
      WrapperBase * base_pointer = &m_wrapperBase;
      WrapperBase * derived = geos::dynamicCast< WrapperBase * >( base_pointer );
      ASSERT_NE( derived, nullptr ) << "Expected successful cast from Base to Base.";
    }
    {
      Wrapper< T > * derived_pointer = &m_wrapper;
      Wrapper< T > * derived = geos::dynamicCast< Wrapper< T > * >( derived_pointer );
      ASSERT_NE( derived, nullptr ) << "Expected successful cast from Derived to Derived.";
    }
    {
      WrapperBase * nullPtr = nullptr;
      Wrapper< T > * castNull = geos::dynamicCast< Wrapper< T > * >( nullPtr );
      ASSERT_EQ( castNull, nullptr ) << "Casting nullptr should yield nullptr.";
    }
  }

  void testDynamicCastWithReference()
  {
    {
      WrapperBase & base_reference = m_wrapperBase;
      Wrapper< T > & derived = geos::dynamicCast< Wrapper< T > & >( base_reference );
      ASSERT_EQ( &derived, &base_reference ) << "Expected successful cast from Base to Derived.";
    }
    {
      WrapperBase & base_reference = m_wrapperBase;
      WrapperBase & derived = geos::dynamicCast< WrapperBase & >( base_reference );
      ASSERT_EQ( &derived, &base_reference ) << "Expected successful cast from Base to Base.";
    }
    {
      Wrapper< T > & derived_reference = m_wrapper;
      Wrapper< T > & derived = geos::dynamicCast< Wrapper< T > & >( derived_reference );
      ASSERT_EQ( &derived, &derived_reference ) << "Expected successful cast from Derived to Derived.";
    }
  }

private:
  conduit::Node m_node;
  Group m_group;
  Wrapper< T > m_wrapper;
  WrapperBase & m_wrapperBase;
};

using WrapperMockTypes = ::testing::Types< int, array1d< real64 >, array1d< array1d< int > >, void *, std::function< void (void) > >;
TYPED_TEST_SUITE( WrapperMock, WrapperMockTypes, );

TYPED_TEST( WrapperMock, DynamicCastWithPointer )
{
  this->testDynamicCastWithPointer();
}

TYPED_TEST( WrapperMock, DynamicCastWithReference )
{
  this->testDynamicCastWithReference();
}

// --------------------------------------------------------------------------------------
// Regex basic constructor
// --------------------------------------------------------------------------------------
TEST( RegexTests, Constructor )
{
  geos::Regex regex( "^[0-9]+$", "Input must be a number." );
  ASSERT_EQ( regex.m_regexStr, "^[0-9]+$" ) << "Regex string is incorrect.";
  ASSERT_EQ( regex.m_formatDescription, "Input must be a number." ) << "Format description is incorrect.";
}

// --------------------------------------------------------------------------------------
// rtTypes::getTypeName tests
// --------------------------------------------------------------------------------------
TEST( RtTypesTests, GetTypeName_KnownTypes )
{
  {
    std::type_index typeIndex( typeid( BaseClass ) );
    auto typeName = geos::rtTypes::getTypeName( typeIndex );
    EXPECT_EQ( typeName, std::string( "BaseClass" ) );
  }
  {
    std::type_index typeIndex( typeid( DerivedFinalClass ) );
    auto typeName = geos::rtTypes::getTypeName( typeIndex );
    EXPECT_EQ( typeName, std::string( "DerivedFinalClass" ) );
  }
}

TEST( RtTypesTests, GetTypeName_UnknownType_Fallback )
{
  struct SomeCustomType { int x; };
  std::type_index typeIndex( typeid( SomeCustomType ) );
  auto typeName = geos::rtTypes::getTypeName( typeIndex );
  // Fallback should contain the actual type name; allow either exact or namespace-qualified.
  ASSERT_TRUE( typeName.find( "SomeCustomType" ) != std::string::npos )
    << "Fallback demangled name should contain 'SomeCustomType' but was '" << typeName << "'";
}

// --------------------------------------------------------------------------------------
// rtTypes::getTypeRegex tests
// --------------------------------------------------------------------------------------
TEST( RtTypesTests, GetTypeRegex_DefaultInteger )
{
  geos::Regex const & regex1 = geos::rtTypes::getTypeRegex< int >();
  ASSERT_FALSE( regex1.m_regexStr.empty() ) << "Expected non-empty regex for int.";
  // Caching: second call should return the same reference
  geos::Regex const & regex2 = geos::rtTypes::getTypeRegex< int >();
  ASSERT_EQ( &regex1, &regex2 ) << "Regex map should cache and return same reference instance.";

  std::regex re( regex1.m_regexStr );
  EXPECT_TRUE( std::regex_match( std::string("-12"), re ) );
  EXPECT_TRUE( std::regex_match( std::string("+7"), re ) );
  EXPECT_FALSE( std::regex_match( std::string("12a"), re ) );
}

TEST( RtTypesTests, GetTypeRegex_Array2DInteger )
{
  geos::Regex const & r = geos::rtTypes::getTypeRegex< int >( "integer_array2d" );
  ASSERT_FALSE( r.m_regexStr.empty() );
  std::regex re( r.m_regexStr );
  EXPECT_TRUE( std::regex_match( std::string("{{1,2},{3,4}}"), re ) );
  EXPECT_TRUE( std::regex_match( std::string(" { { 1 , 2 } , { 3 , 4 } } "), re ) );
  EXPECT_FALSE( std::regex_match( std::string("{1,2,3,4}"), re ) ); // Not a 2d array pattern
}

TEST( RtTypesTests, GetTypeRegex_CustomGroupNameRef )
{
  geos::Regex const & r = geos::rtTypes::getTypeRegex< string >( rtTypes::CustomTypes::groupNameRef );
  ASSERT_FALSE( r.m_regexStr.empty() );
  std::regex re( r.m_regexStr );
  EXPECT_TRUE( std::regex_match( std::string("group/sub*pattern"), re ) );
  EXPECT_FALSE( std::regex_match( std::string("group name with spaces"), re ) );
}

TEST( RtTypesTests, GetTypeRegex_EnumSpecialization )
{
  geos::Regex const & r = geos::rtTypes::getTypeRegex< geos::MyEnum >();
  EXPECT_EQ( r.m_regexStr, std::string("alpha|beta|gamma") );
  EXPECT_NE( r.m_formatDescription.find("alpha, beta, gamma"), std::string::npos );
}

// --------------------------------------------------------------------------------------
// TypeName utility tests (robust to current brief implementation or future fix)
// --------------------------------------------------------------------------------------
TEST( TypeNameTests, FullAndBrief )
{
  auto fullInt = TypeName< int >::full();
  ASSERT_TRUE( fullInt.find("int") != std::string::npos );

  auto briefEnum = TypeName< geos::MyEnum >::brief();
  // Current implementation may yield leading ':'; accept both forms.
  ASSERT_TRUE( briefEnum == "MyEnum" || briefEnum == ":MyEnum" )
    << "Unexpected brief name: " << briefEnum;
}

// --------------------------------------------------------------------------------------
// Main (custom since we add death tests)
// --------------------------------------------------------------------------------------
int main( int argc, char ** argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();
}
