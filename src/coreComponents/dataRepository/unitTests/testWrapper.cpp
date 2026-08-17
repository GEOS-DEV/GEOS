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
#include "common/TypeDispatch.hpp"
#include "dataRepository/Group.hpp"
#include "dataRepository/Wrapper.hpp"
#include "LvArray/src/system.hpp"

// TPL includes
#include <gtest/gtest.h>
#include <conduit.hpp>

// System includes
#include <cstdlib>
#include <stdexcept>
#include <tuple>

using namespace geos;
using namespace dataRepository;

namespace
{

class ScopedThrowingLvArrayErrorHandler
{
public:
  ScopedThrowingLvArrayErrorHandler()
  {
    LvArray::system::setErrorHandler( []()
    {
      throw std::runtime_error( "Expected WrapperStandardArrays error" );
    } );
  }

  ~ScopedThrowingLvArrayErrorHandler() noexcept
  {
    LvArray::system::setErrorHandler( []()
    {
      std::abort();
    } );
  }

  ScopedThrowingLvArrayErrorHandler( ScopedThrowingLvArrayErrorHandler const & ) = delete;
  ScopedThrowingLvArrayErrorHandler & operator=( ScopedThrowingLvArrayErrorHandler const & ) = delete;
};

template< typename ... ARRAY_TYPES >
void testStandardArrayVirtualInterface( camp::list< ARRAY_TYPES ... > )
{
  ScopedThrowingLvArrayErrorHandler const scopedErrorHandler;

  // Keep the arrays outside the Group so both the registered wrappers and their clones are
  // non-owning. This makes it safe to exercise clone() without sharing ownership of the data.
  std::tuple< ARRAY_TYPES ... > storage;
  conduit::Node rootNode;
  Group group( "Problem", rootNode );

  localIndex wrapperIndex = 0;
  std::apply( [&]( auto & ... value )
  {
    ( group.registerWrapper( "standardArray" + std::to_string( wrapperIndex++ ), &value ), ... );
  }, storage );

  xmlWrapper::xmlDocument emptyDocument;
  xmlWrapper::xmlNodePos const noPosition( emptyDocument,
                                           "",
                                           xmlWrapper::xmlDocument::npos,
                                           xmlWrapper::xmlDocument::npos,
                                           xmlWrapper::xmlDocument::npos );
  xmlWrapper::xmlNode const emptyNode;

  string_array const dimLabels = { "first", "second" };
  array1d< localIndex > packIndices( 1 );
  packIndices[0] = 0;
  parallelDeviceEvents events;
  localIndex numWrappersTested = 0;

  group.forWrappers( [&]( WrapperBase & wrapper )
  {
    ++numWrappersTested;

    EXPECT_FALSE( wrapper.getName().empty() );
    EXPECT_FALSE( wrapper.getPath().empty() );
    EXPECT_EQ( &wrapper.getParent(), &group );
    EXPECT_FALSE( wrapper.getSuccessfulReadFromInput() );
    EXPECT_FALSE( wrapper.processInputFile( emptyNode, noPosition ) );

    int const numDims = wrapper.numArrayDims();
    ASSERT_GT( numDims, 0 );
    stdVector< localIndex > const dims( numDims, 2 );
    wrapper.resize( numDims, dims.data() );

    EXPECT_GT( wrapper.size(), 0 );
    EXPECT_NE( wrapper.voidPointer(), nullptr );
    EXPECT_GT( wrapper.elementByteSize(), 0 );
    EXPECT_GT( wrapper.bytesAllocated(), 0 );

    wrapper.resize( 3 );
    wrapper.resize( numDims, dims.data() );
    wrapper.reserve( 3 );
    EXPECT_GT( wrapper.capacity(), 0 );

    // These validation failures occur before Array dimensions or storage are modified.
    EXPECT_THROW( wrapper.resize( numDims + 1, dims.data() ), std::runtime_error );
    EXPECT_THROW( wrapper.resize( localIndex( -1 ) ), std::runtime_error );
    EXPECT_THROW( wrapper.move( static_cast< LvArray::MemorySpace >( -1 ), false ), std::runtime_error );

    for( integer dim = 0; dim < numDims; ++dim )
    {
      EXPECT_EQ( &wrapper.setDimLabels( dim, dimLabels ), &wrapper );
      Span< string const > const labels = wrapper.getDimLabels( dim );
      ASSERT_EQ( labels.size(), dimLabels.size() );
      EXPECT_EQ( labels[0], dimLabels[0] );
      EXPECT_EQ( labels[1], dimLabels[1] );
    }

    EXPECT_THROW( wrapper.setDimLabels( -1, dimLabels ), std::runtime_error );
    EXPECT_THROW( wrapper.setDimLabels( numDims, dimLabels ), std::runtime_error );
    EXPECT_THROW( wrapper.getDimLabels( -1 ), std::runtime_error );
    EXPECT_THROW( wrapper.getDimLabels( numDims ), std::runtime_error );

    stdVector< string > const invalidComponentNames( wrapper.numArrayComp() + 1, "invalid" );
    conduit::Node invalidFields;
    conduit::Node invalidMCArray;
    EXPECT_THROW( wrapper.addBlueprintField( invalidFields,
                                             wrapper.getName(),
                                             "topology",
                                             invalidComponentNames ),
                  std::runtime_error );
    EXPECT_THROW( wrapper.populateMCArray( invalidMCArray, invalidComponentNames ), std::runtime_error );

    wrapper.resize( 0 );
    conduit::Node emptyFields;
    conduit::Node emptyMCArray;
    EXPECT_THROW( wrapper.addBlueprintField( emptyFields, wrapper.getName(), "topology" ), std::runtime_error );
    EXPECT_THROW( wrapper.populateMCArray( emptyMCArray ), std::runtime_error );
    wrapper.resize( numDims, dims.data() );

    // Restrict dispatch to the production type catalog so the test does not introduce unrelated
    // Wrapper<T> types solely for coverage.
    EXPECT_TRUE( types::dispatch( types::ListofTypeList< types::StandardArrays >{}, [&]( auto typeTuple )
    {
      using ArrayType = camp::first< decltype( typeTuple ) >;
      using ValueType = typename ArrayType::ValueType;

      Wrapper< ArrayType > & typedWrapper = Wrapper< ArrayType >::cast( wrapper );
      WrapperBase const & constBase = wrapper;
      EXPECT_EQ( &Wrapper< ArrayType >::cast( constBase ), &typedWrapper );

      ValueType const defaultValue = ValueType( 1 );
      typedWrapper.setDefaultValue( defaultValue );
      typedWrapper.setApplyDefaultValue( defaultValue );
      EXPECT_EQ( typedWrapper.getDefaultValue(), defaultValue );
      EXPECT_TRUE( typedWrapper.hasDefaultValue() );

      typedWrapper.setSizedFromParent( 1 );
      typedWrapper.setRestartFlags( RestartFlags::WRITE_AND_READ );
      typedWrapper.setPlotLevel( PlotLevel::LEVEL_1 );
      typedWrapper.setInputFlag( InputFlags::OPTIONAL );
      typedWrapper.setDescription( "standard array" );
      typedWrapper.appendDescription( " wrapper" );
      typedWrapper.setRegisteringObjects( "testWrapper" );
      typedWrapper.setRTTypeName( typedWrapper.getRTTypeName() );
      typedWrapper.setName();

      EXPECT_EQ( typedWrapper.getDescription(), "standard array wrapper" );
      EXPECT_EQ( typedWrapper.getRegisteringObjects().count( "testWrapper" ), 1 );
      EXPECT_EQ( typedWrapper.referenceAsView().size(), typedWrapper.reference().size() );

      ArrayType & array = typedWrapper.reference();
      for( localIndex valueIndex = 0; valueIndex < array.size(); ++valueIndex )
      {
        array.data()[valueIndex] = static_cast< ValueType >( valueIndex + 2 );
      }
    }, wrapper ) );

    EXPECT_FALSE( wrapper.getDefaultValueString().empty() );
    static_cast< void >( wrapper.getTypeRegex() );
    EXPECT_TRUE( wrapper.isPackable( false ) );
    static_cast< void >( wrapper.isPackable( true ) );

    HistoryMetadata const history = wrapper.getHistoryMetadata( 1 );
    EXPECT_EQ( history.getName(), wrapper.getName() );
    EXPECT_GT( history.getRank(), 0 );
    EXPECT_GT( history.size(), 0 );

    conduit::Node fields;
    wrapper.addBlueprintField( fields, wrapper.getName(), "topology" );
    EXPECT_GT( fields.number_of_children(), 0 );

    conduit::Node mcArray;
    wrapper.populateMCArray( mcArray );
    EXPECT_GT( mcArray.number_of_children(), 0 );

    std::unique_ptr< WrapperBase > average =
      wrapper.averageOverSecondDim( wrapper.getName() + "Average", group );
    ASSERT_NE( average, nullptr );
    EXPECT_GT( average->size(), 0 );

    std::unique_ptr< WrapperBase > clone = wrapper.clone( wrapper.getName() + "Clone", group );
    ASSERT_NE( clone, nullptr );
    EXPECT_TRUE( clone->getTypeId() == wrapper.getTypeId() );
    EXPECT_EQ( clone->voidPointer(), wrapper.voidPointer() );
    EXPECT_EQ( clone->bytesAllocated(), 0 );
    EXPECT_THROW( clone->copyWrapper( wrapper ), std::runtime_error );

    for( bool const withMetadata : { false, true } )
    {
      buffer_unit_type * unusedBuffer = nullptr;
      localIndex const packedSize = wrapper.pack< false >( unusedBuffer, withMetadata, false, events );
      ASSERT_GT( packedSize, 0 );

      buffer_type packed( packedSize );
      buffer_unit_type * packBuffer = packed.data();
      EXPECT_EQ( wrapper.pack< true >( packBuffer, withMetadata, false, events ), packedSize );
      EXPECT_EQ( packBuffer, packed.data() + packedSize );

      if( withMetadata )
      {
        buffer_unit_type const * mismatchedBuffer = packed.data();
        EXPECT_THROW( clone->unpack( mismatchedBuffer, true, false, events ), std::runtime_error );
      }

      EXPECT_TRUE( types::dispatch( types::ListofTypeList< types::StandardArrays >{}, [&]( auto typeTuple )
      {
        using ArrayType = camp::first< decltype( typeTuple ) >;
        using ValueType = typename ArrayType::ValueType;
        ArrayType & array = Wrapper< ArrayType >::cast( wrapper ).reference();
        for( localIndex valueIndex = 0; valueIndex < array.size(); ++valueIndex )
        {
          array.data()[valueIndex] = ValueType( 0 );
        }
      }, wrapper ) );

      buffer_unit_type const * unpackBuffer = packed.data();
      EXPECT_EQ( wrapper.unpack( unpackBuffer, withMetadata, false, events ), packedSize );
      EXPECT_EQ( unpackBuffer, packed.data() + packedSize );

      EXPECT_TRUE( types::dispatch( types::ListofTypeList< types::StandardArrays >{}, [&]( auto typeTuple )
      {
        using ArrayType = camp::first< decltype( typeTuple ) >;
        using ValueType = typename ArrayType::ValueType;
        ArrayType const & array = Wrapper< ArrayType >::cast( wrapper ).reference();
        for( localIndex valueIndex = 0; valueIndex < array.size(); ++valueIndex )
        {
          EXPECT_EQ( array.data()[valueIndex], static_cast< ValueType >( valueIndex + 2 ) );
        }
      }, wrapper ) );

      unusedBuffer = nullptr;
      localIndex const indexedPackedSize =
        wrapper.packByIndex< false >( unusedBuffer, packIndices.toViewConst(), withMetadata, false, events );
      ASSERT_GT( indexedPackedSize, 0 );

      buffer_type indexedPacked( indexedPackedSize );
      packBuffer = indexedPacked.data();
      EXPECT_EQ( wrapper.packByIndex< true >( packBuffer,
                                              packIndices.toViewConst(),
                                              withMetadata,
                                              false,
                                              events ),
                 indexedPackedSize );
      EXPECT_EQ( packBuffer, indexedPacked.data() + indexedPackedSize );

      if( withMetadata )
      {
        buffer_unit_type const * mismatchedBuffer = indexedPacked.data();
        EXPECT_THROW( clone->unpackByIndex( mismatchedBuffer,
                                            packIndices.toViewConst(),
                                            true,
                                            false,
                                            events,
                                            MPI_REPLACE ),
                      std::runtime_error );
      }

      EXPECT_TRUE( types::dispatch( types::ListofTypeList< types::StandardArrays >{}, [&]( auto typeTuple )
      {
        using ArrayType = camp::first< decltype( typeTuple ) >;
        using ValueType = typename ArrayType::ValueType;
        ArrayType & array = Wrapper< ArrayType >::cast( wrapper ).reference();
        LvArray::forValuesInSliceWithIndices( array[0], [&]( auto & value, auto const ... )
        {
          value = ValueType( 0 );
        } );
      }, wrapper ) );

      unpackBuffer = indexedPacked.data();
      EXPECT_EQ( wrapper.unpackByIndex( unpackBuffer,
                                        packIndices.toViewConst(),
                                        withMetadata,
                                        false,
                                        events,
                                        MPI_REPLACE ),
                 indexedPackedSize );
      EXPECT_EQ( unpackBuffer, indexedPacked.data() + indexedPackedSize );

      EXPECT_TRUE( types::dispatch( types::ListofTypeList< types::StandardArrays >{}, [&]( auto typeTuple )
      {
        using ArrayType = camp::first< decltype( typeTuple ) >;
        using ValueType = typename ArrayType::ValueType;
        ArrayType const & array = Wrapper< ArrayType >::cast( wrapper ).reference();
        for( localIndex valueIndex = 0; valueIndex < array.size(); ++valueIndex )
        {
          EXPECT_EQ( array.data()[valueIndex], static_cast< ValueType >( valueIndex + 2 ) );
        }
      }, wrapper ) );
    }

    wrapper.copy( 0, 1 );
    EXPECT_TRUE( types::dispatch( types::ListofTypeList< types::StandardArrays >{}, [&]( auto typeTuple )
    {
      using ArrayType = camp::first< decltype( typeTuple ) >;
      ArrayType const & array = Wrapper< ArrayType >::cast( wrapper ).reference();
      LvArray::forValuesInSliceWithIndices( array[0], [&]( auto const & expectedValue, auto const ... indices )
      {
        EXPECT_EQ( array( 1, indices ... ), expectedValue );
      } );
    }, wrapper ) );
    wrapper.move( hostMemorySpace, false );
    wrapper.move( hostMemorySpace, true );

    wrapper.setRestartFlags( RestartFlags::NO_WRITE );
    wrapper.registerToWrite();
    wrapper.finishWriting();
    wrapper.setRestartFlags( RestartFlags::WRITE );
    wrapper.registerToWrite();
    wrapper.finishWriting();
    EXPECT_FALSE( wrapper.loadFromConduit() );

    stdVector< camp::idx_t > invalidMetadata( numDims + 1, 0 );
    conduit::DataType const invalidMetadataType( conduitTypeInfo< camp::idx_t >::id,
                                                 invalidMetadata.size() );
    wrapper.setRestartFlags( RestartFlags::WRITE_AND_READ );

    wrapper.registerToWrite();
    group.getConduitNode()[wrapper.getName()]["__permutation__"].set( invalidMetadataType,
                                                                      invalidMetadata.data() );
    EXPECT_THROW( wrapper.loadFromConduit(), std::runtime_error );

    wrapper.registerToWrite();
    camp::idx_t * const permutation =
      group.getConduitNode()[wrapper.getName()]["__permutation__"].value();
    permutation[0] = numDims;
    EXPECT_THROW( wrapper.loadFromConduit(), std::runtime_error );

    wrapper.registerToWrite();
    group.getConduitNode()[wrapper.getName()]["__dimensions__"].set( invalidMetadataType,
                                                                     invalidMetadata.data() );
    EXPECT_THROW( wrapper.loadFromConduit(), std::runtime_error );

    wrapper.registerToWrite();
    buffer_unit_type invalidValue = 0;
    conduit::DataType const invalidValueType( conduitTypeInfo< buffer_unit_type >::id, 1 );
    group.getConduitNode()[wrapper.getName()]["__values__"].set( invalidValueType, &invalidValue );
    EXPECT_THROW( wrapper.loadFromConduit(), std::runtime_error );

    wrapper.registerToWrite();
    wrapper.finishWriting();

    EXPECT_THROW( wrapper.erase( std::set< localIndex >{} ), std::runtime_error );
    wrapper.erase( { 0 } );
  } );

  EXPECT_EQ( numWrappersTested, sizeof...( ARRAY_TYPES ) );
}

} // namespace

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

class WrapperLimitsTest : public ::testing::Test
{
protected:
  WrapperLimitsTest():
    m_node(),
    m_group( "Problem", m_node )
  {}

  template< typename T >
  Wrapper< T > & makeWrapper( string const & name )
  {
    return m_group.template registerWrapper< T >( name );
  }

  conduit::Node m_node;
  Group m_group;
};

TEST_F( WrapperLimitsTest, IsLimitableTrait )
{
  static_assert( wrapperLimits::is_limitable_v< integer >, "integer must be limitable" );
  static_assert( wrapperLimits::is_limitable_v< real64 >, "real64 must be limitable" );
  static_assert( wrapperLimits::is_limitable_v< array1d< integer > >, "array1d< integer > must be limitable" );
  static_assert( wrapperLimits::is_limitable_v< array2d< real64 > >, "array2d< real64 > must be limitable" );
  static_assert( wrapperLimits::is_limitable_v< array3d< integer > >, "array3d< integer > must be limitable" );

  static_assert( std::is_same< wrapperLimits::limit_value_type_t< real64 >, real64 >::value, "" );
  static_assert( std::is_same< wrapperLimits::limit_value_type_t< array1d< real64 > >, real64 >::value, "" );
  static_assert( std::is_same< wrapperLimits::limit_value_type_t< array2d< integer > >, integer >::value, "" );
}

TEST_F( WrapperLimitsTest, ScalarSetGet )
{
  auto & w = makeWrapper< real64 >( "scalar" );
  w.setLimits( wrapperLimits::inclusive( 0.0 ), wrapperLimits::exclusive( 1.0 ), wrapperLimits::LimitsMode::Error );

  ASSERT_TRUE( w.getMinBound().has_value() );
  ASSERT_TRUE( w.getMaxBound().has_value() );
  EXPECT_DOUBLE_EQ( w.getMinBound()->value, 0.0 );
  EXPECT_TRUE( w.getMinBound()->isInclusive );
  EXPECT_DOUBLE_EQ( w.getMaxBound()->value, 1.0 );
  EXPECT_FALSE( w.getMaxBound()->isInclusive );
  EXPECT_EQ( w.getLimitsMode(), wrapperLimits::LimitsMode::Error );
}

TEST_F( WrapperLimitsTest, Array1dSetGet )
{
  auto & w = makeWrapper< array1d< real64 > >( "array1d" );
  w.setLimits( 0.0, 1.0, wrapperLimits::LimitsMode::Error );

  ASSERT_TRUE( w.getMinBound().has_value() );
  ASSERT_TRUE( w.getMaxBound().has_value() );
  EXPECT_DOUBLE_EQ( w.getMinBound()->value, 0.0 );
  EXPECT_DOUBLE_EQ( w.getMaxBound()->value, 1.0 );
  EXPECT_EQ( w.getLimitsMode(), wrapperLimits::LimitsMode::Error );
}

TEST_F( WrapperLimitsTest, ScalarValidateInRange )
{
  auto & w = makeWrapper< real64 >( "scalar" );
  w.setLimits( 0.0, 1.0, wrapperLimits::LimitsMode::Error );
  w.reference() = 0.5;
  EXPECT_NO_THROW( w.validateLimits() );
}

TEST_F( WrapperLimitsTest, ScalarValidateBelowMinThrows )
{
  auto & w = makeWrapper< real64 >( "scalar" );
  w.setLimits( 0.0, 1.0, wrapperLimits::LimitsMode::Error );
  w.reference() = -0.1;
  EXPECT_THROW( w.validateLimits(), InputError );
}

TEST_F( WrapperLimitsTest, ScalarValidateAboveMaxThrows )
{
  auto & w = makeWrapper< real64 >( "scalar" );
  w.setLimits( 0.0, 1.0, wrapperLimits::LimitsMode::Error );
  w.reference() = 1.1;
  EXPECT_THROW( w.validateLimits(), InputError );
}

TEST_F( WrapperLimitsTest, ScalarValidateInclusiveBoundary )
{
  auto & w = makeWrapper< real64 >( "scalar" );
  w.setLimits( wrapperLimits::inclusive( 0.0 ), wrapperLimits::inclusive( 1.0 ), wrapperLimits::LimitsMode::Error );

  w.reference() = 0.0;
  EXPECT_NO_THROW( w.validateLimits() );

  w.reference() = 1.0;
  EXPECT_NO_THROW( w.validateLimits() );
}

TEST_F( WrapperLimitsTest, ScalarValidateExclusiveBoundary )
{
  auto & w = makeWrapper< real64 >( "scalar" );
  w.setLimits( wrapperLimits::exclusive( 0.0 ), wrapperLimits::exclusive( 1.0 ), wrapperLimits::LimitsMode::Error );

  w.reference() = 0.0;
  EXPECT_THROW( w.validateLimits(), InputError );

  w.reference() = 0.5;
  EXPECT_NO_THROW( w.validateLimits() );

  w.reference() = 1.0;
  EXPECT_THROW( w.validateLimits(), InputError );
}

TEST_F( WrapperLimitsTest, Array1dValidateAllInRange )
{
  auto & w = makeWrapper< array1d< real64 > >( "array1d" );
  w.setLimits( 0.0, 1.0, wrapperLimits::LimitsMode::Error );

  array1d< real64 > & data = w.reference();
  data.resize( 4 );
  data[ 0 ] = 0.0;
  data[ 1 ] = 0.25;
  data[ 2 ] = 0.75;
  data[ 3 ] = 1.0;

  EXPECT_NO_THROW( w.validateLimits() );
}

TEST_F( WrapperLimitsTest, Array1dValidateOutOfRange )
{
  auto & w = makeWrapper< array1d< real64 > >( "array1d" );
  w.setLimits( 0.0, 1.0, wrapperLimits::LimitsMode::Error );

  array1d< real64 > & data = w.reference();
  data.resize( 4 );
  data[ 0 ] = 0.5;
  data[ 1 ] = 0.5;
  data[ 2 ] = 42.0;
  data[ 3 ] = 0.5;

  EXPECT_THROW( w.validateLimits(), InputError );
}

TEST_F( WrapperLimitsTest, Array1dValidateEmpty )
{
  auto & w = makeWrapper< array1d< real64 > >( "array1d" );
  w.setLimits( 0.0, 1.0, wrapperLimits::LimitsMode::Error );

  EXPECT_NO_THROW( w.validateLimits() );
}

TEST_F( WrapperLimitsTest, Array2dValidateAllInRange )
{
  auto & w = makeWrapper< array2d< real64 > >( "array2d" );
  w.setLimits( 0.0, 1.0, wrapperLimits::LimitsMode::Error );

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

TEST_F( WrapperLimitsTest, Array2dValidateOutOfRange )
{
  auto & w = makeWrapper< array2d< real64 > >( "array2d" );
  w.setLimits( 0.0, 1.0, wrapperLimits::LimitsMode::Error );

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

TEST_F( WrapperLimitsTest, Array2dValidateEmpty )
{
  auto & w = makeWrapper< array2d< real64 > >( "array2d" );
  w.setLimits( 0.0, 1.0, wrapperLimits::LimitsMode::Error );

  EXPECT_NO_THROW( w.validateLimits() );
}

TEST( WrapperStandardArrays, VirtualInterface )
{
  testStandardArrayVirtualInterface( types::StandardArrays{} );
}
