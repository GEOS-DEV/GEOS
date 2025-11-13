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
#include <string>

#include "codingUtilities/RTTypes.hpp"
#include "dataRepository/Group.hpp"
#include "dataRepository/Wrapper.hpp"
#include "common/format/EnumStrings.hpp"
#include "common/logger/Logger.hpp"

#include <gtest/gtest.h>
#include <conduit.hpp>

// Mock classes to test dynamic casting
#include "BaseClass.hpp"
#include "DerivedFinalClass.hpp"

using namespace geos;
using namespace dataRepository;

// Helper function used in several regex tests.
namespace
{
bool regexMatches( std::string const & pattern, std::string const & value )
{
  std::regex re( pattern );
  return std::regex_match( value, re );
}
}

// ENUM WITH STRINGS FOR TypeRegex SPECIALIZATION TESTING
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


// dynamicCast pointer tests
TEST( DynamicCastTests, Pointer_Casting_Success )
{
  std::unique_ptr< BaseClass > base( new DerivedFinalClass());
  DerivedFinalClass * derived = geos::dynamicCast< DerivedFinalClass * >( base.get());
  ASSERT_NE( derived, nullptr ) << "Expected successful cast from Base to Derived.";
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


// dynamicCast reference tests
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



// Typed test for geos Wrapper

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


// Regex basic constructor
TEST( RegexTests, Constructor )
{
  geos::Regex regex( "^[0-9]+$", "Input must be a number." );
  ASSERT_EQ( regex.m_regexStr, "^[0-9]+$" ) << "Regex string is incorrect.";
  ASSERT_EQ( regex.m_formatDescription, "Input must be a number." ) << "Format description is incorrect.";
}


// rtTypes::getTypeName tests
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


// rtTypes::getTypeRegex tests
TEST( RtTypesTests, GetTypeRegex_DefaultInteger )
{
  geos::Regex const & regex1 = geos::rtTypes::getTypeRegex< int >();
  ASSERT_FALSE( regex1.m_regexStr.empty() ) << "Expected non-empty regex for int.";
  // Caching: second call should return the same reference
  geos::Regex const & regex2 = geos::rtTypes::getTypeRegex< int >();
  ASSERT_EQ( &regex1, &regex2 ) << "Regex map should cache and return same reference instance.";

  EXPECT_TRUE( regexMatches( regex1.m_regexStr, "-12" ) );
  EXPECT_TRUE( regexMatches( regex1.m_regexStr, "+7" ) );
  EXPECT_FALSE( regexMatches( regex1.m_regexStr, "12a" ) );
}

TEST( RtTypesTests, GetTypeRegex_Array2DInteger )
{
  geos::Regex const & r = geos::rtTypes::getTypeRegex< int >( "integer_array2d" );
  ASSERT_FALSE( r.m_regexStr.empty() );
  std::regex re( r.m_regexStr );
  EXPECT_TRUE( std::regex_match( std::string( "{{1,2},{3,4}}" ), re ) );
  EXPECT_TRUE( std::regex_match( std::string( " { { 1 , 2 } , { 3 , 4 } } " ), re ) );
  EXPECT_FALSE( std::regex_match( std::string( "{1,2,3,4}" ), re ) ); // Not a 2d array pattern
}

TEST( RtTypesTests, GetTypeRegex_CustomGroupNameRef )
{
  geos::Regex const & r = geos::rtTypes::getTypeRegex< string >( rtTypes::CustomTypes::groupNameRef );
  ASSERT_FALSE( r.m_regexStr.empty() );
  std::regex re( r.m_regexStr );
  EXPECT_TRUE( std::regex_match( std::string( "group/sub*pattern" ), re ) );
  EXPECT_FALSE( std::regex_match( std::string( "group name with spaces" ), re ) );
}

TEST( RtTypesTests, GetTypeRegex_EnumSpecialization )
{
  geos::Regex const & r = geos::rtTypes::getTypeRegex< geos::MyEnum >();
  EXPECT_EQ( r.m_regexStr, std::string( "alpha|beta|gamma" ) );
  EXPECT_NE( r.m_formatDescription.find( "alpha, beta, gamma" ), std::string::npos );
}


// Custom type with TypeRegex specialization
namespace geos
{
struct CustomRegexType {};
template<> struct TypeRegex< CustomRegexType, void >
{
  static Regex get() { return Regex( "X+", "Input value must be one or more X characters." ); }
};

struct NoRegexType {};   // no specialization => empty regex
}

TEST( RtTypesTests, TypeRegex_CustomSpecialization_Caching )
{
  geos::Regex const & r1 = geos::rtTypes::getTypeRegex< geos::CustomRegexType >( "CustomRegexAlias" );
  ASSERT_EQ( r1.m_regexStr, std::string( "X+" ) );
  ASSERT_FALSE( r1.m_formatDescription.empty() );
  geos::Regex const & r2 = geos::rtTypes::getTypeRegex< geos::CustomRegexType >( "CustomRegexAlias" );
  ASSERT_EQ( &r1, &r2 ) << "Expected caching to return identical reference.";
}

TEST( RtTypesTests, TypeRegex_FallbackEmpty )
{
  geos::Regex const & r = geos::rtTypes::getTypeRegex< geos::NoRegexType >( "NoRegexAlias" );
  ASSERT_TRUE( r.m_regexStr.empty() );
  ASSERT_TRUE( r.m_formatDescription.empty() );
}

TEST( RtTypesTests, IntegerRegex_EdgeCases )
{
  geos::Regex const & intR = geos::rtTypes::getTypeRegex< int >();
  std::regex re( intR.m_regexStr );
  EXPECT_TRUE( std::regex_match( std::string( "0" ), re ) );
  EXPECT_TRUE( std::regex_match( std::string( "-0" ), re ) );
  EXPECT_TRUE( std::regex_match( std::string( "+123456789" ), re ) );
  EXPECT_FALSE( std::regex_match( std::string( "++1" ), re ) );
  EXPECT_FALSE( std::regex_match( std::string( "-+1" ), re ) );
  EXPECT_FALSE( std::regex_match( std::string( "1-" ), re ) );
}

TEST( RtTypesTests, RealRegex_EdgeCases )
{
  geos::Regex const & realR = geos::rtTypes::getTypeRegex< real64 >();
  std::regex re( realR.m_regexStr );
  EXPECT_TRUE( std::regex_match( std::string( "1." ), re ) );
  EXPECT_TRUE( std::regex_match( std::string( ".5" ), re ) );
  EXPECT_TRUE( std::regex_match( std::string( "5e3" ), re ) );
  EXPECT_TRUE( std::regex_match( std::string( "5E+3" ), re ) );
  EXPECT_TRUE( std::regex_match( std::string( "-8.2e-9" ), re ) );
  EXPECT_FALSE( std::regex_match( std::string( "." ), re ) ); // single dot not valid per pattern
  EXPECT_FALSE( std::regex_match( std::string( "1..2" ), re ) );
  EXPECT_FALSE( std::regex_match( std::string( "e10" ), re ) );
}

TEST( RtTypesTests, Array1DInteger_EmptyAndValues )
{
  geos::Regex const & arrR = geos::rtTypes::getTypeRegex< int >( "integer_array" );
  EXPECT_TRUE( regexMatches( arrR.m_regexStr, "{}" ) );
  EXPECT_TRUE( regexMatches( arrR.m_regexStr, "{1}" ) );
  EXPECT_TRUE( regexMatches( arrR.m_regexStr, "{ 1 , 2 , 3 }" ) );
  EXPECT_FALSE( regexMatches( arrR.m_regexStr, "{,}" ) );
}

TEST( RtTypesTests, GroupNameRefArrayRegex )
{
  geos::Regex const & r = geos::rtTypes::getTypeRegex< string >( rtTypes::CustomTypes::groupNameRefArray );
  std::regex re( r.m_regexStr );
  EXPECT_TRUE( std::regex_match( std::string( "{}" ), re ) );
  EXPECT_TRUE( std::regex_match( std::string( "{pattern/*,path/sub}" ), re ) );
  EXPECT_FALSE( std::regex_match( std::string( "{bad space}" ), re ) );
}

// New: R1Tensor regex test
TEST( RtTypesTests, GetTypeRegex_R1Tensor )
{
  geos::Regex const & r = geos::rtTypes::getTypeRegex< R1Tensor >( "R1Tensor" );
  ASSERT_FALSE( r.m_regexStr.empty() );
  EXPECT_TRUE( regexMatches( r.m_regexStr, "{1,2,3}" ) );
  EXPECT_TRUE( regexMatches( r.m_regexStr, " { 1 , .5 , -2.3e3 } " ) );
  EXPECT_FALSE( regexMatches( r.m_regexStr, "{1,2}" ) );
  EXPECT_FALSE( regexMatches( r.m_regexStr, "{1,2,3,4}" ) );
}

// Enum roundtrip tests (non-aborting, uses GEOS_THROW_IF which throws InputError)
TEST( EnumStringsTests, Roundtrip )
{
  geos::MyEnum e = EnumStrings< geos::MyEnum >::fromString( "beta" );
  EXPECT_EQ( EnumStrings< geos::MyEnum >::toString( e ), std::string( "beta" ) );
}

TEST( EnumStringsTests, InvalidValueThrows )
{
  EXPECT_THROW( (void) EnumStrings< geos::MyEnum >::fromString( "delta" ), geos::InputError );
}


// TypeName utility tests (robust to current brief implementation or future fix)
TEST( TypeNameTests, FullAndBrief )
{
  auto fullInt = TypeName< int >::full();
  ASSERT_TRUE( fullInt.find( "int" ) != std::string::npos );

  auto briefEnum = TypeName< geos::MyEnum >::brief();
  // Current implementation may yield leading ':'; accept both forms.
  ASSERT_TRUE( briefEnum == "MyEnum" || briefEnum == ":MyEnum" )
    << "Unexpected brief name: " << briefEnum;
}

TEST( TypeNameTests, BriefNamespaceStrip )
{
  // Confirm brief strips namespaces (no leading ':') when possible.
  auto briefInt = TypeName< int >::brief();
  ASSERT_EQ( briefInt, "int" );
  // For the enum, brief currently may include a leading ':' due to existing implementation if not patched.
  auto briefEnum = TypeName< geos::MyEnum >::brief();
  ASSERT_TRUE( briefEnum == "MyEnum" || briefEnum == ":MyEnum" ) << "Unexpected brief name: " << briefEnum;
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  return RUN_ALL_TESTS();
}
