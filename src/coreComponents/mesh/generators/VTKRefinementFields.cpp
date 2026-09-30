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

/** @file VTKRefinementFields.cpp */
#include "VTKRefinementFields.hpp"

#include <vtkArrayDispatch.h>
#include <vtkCellData.h>
#include <vtkDataArrayAccessor.h>
#include <vtkPointData.h>
#include <vtkBitArray.h>
#include <vtkStringArray.h>
#include <vtkStdString.h>
#include <vtkVariant.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <functional>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace geos::vtk::refinement
{
namespace
{
using Bytes = std::vector< unsigned char >;

void appendInteger( Bytes & bytes, std::uint64_t value, unsigned int width = 8 )
{
  for( unsigned int i = 0; i < width; ++i )
    bytes.push_back( static_cast< unsigned char >( value >> ( 8 * i ) ) );
}

void appendString( Bytes & bytes, char const * text )
{
  std::size_t const length = text == nullptr ? 0 : std::strlen( text );
  appendInteger( bytes, length );
  if( length )
    bytes.insert( bytes.end(), text, text + length );
}

template < typename Value >
using Bits = std::conditional_t<
    sizeof( Value ) == 1, std::uint8_t,
    std::conditional_t< sizeof( Value ) == 2, std::uint16_t, std::conditional_t< sizeof( Value ) == 4, std::uint32_t, std::uint64_t > > >;

struct PackTuple
{
  vtkIdType point;
  Bytes & bytes;
  template < typename Array > void operator()( Array * array ) const
  {
    vtkDataArrayAccessor< Array > access( array );
    using Value = typename decltype( access )::APIType;
    static_assert( sizeof( Value ) == sizeof( Bits< Value > ) );
    for( int c = 0; c < array->GetNumberOfComponents(); ++c )
    {
      Value const value = access.Get( point, c );
      Bits< Value > bits;
      std::memcpy( &bits, &value, sizeof( value ) );
      appendInteger( bytes, bits, sizeof( bits ) );
    }
  }
};

struct InstallTuple
{
  vtkIdType point;
  Bytes const & bytes;
  std::size_t & cursor;
  template < typename Array > void operator()( Array * array ) const
  {
    vtkDataArrayAccessor< Array > access( array );
    using Value = typename decltype( access )::APIType;
    for( int c = 0; c < array->GetNumberOfComponents(); ++c )
    {
      if( bytes.size() - cursor < sizeof( Value ) )
        throw std::invalid_argument( "Truncated refinement point-field tuple" );
      Bits< Value > bits{};
      for( unsigned int i = 0; i < sizeof( bits ); ++i )
        bits |= static_cast< Bits< Value > >( static_cast< std::uint64_t >( bytes[cursor++] ) << ( 8 * i ) );
      Value value;
      std::memcpy( &value, &bits, sizeof( value ) );
      access.Set( point, c, value );
    }
  }
};

void validateArrays( vtkDataSetAttributes & input, vtkIdType tuples )
{
  std::set< std::string > names;
  for( int i = 0; i < input.GetNumberOfArrays(); ++i )
  {
    vtkAbstractArray * array = input.GetAbstractArray( i );
    if( array->GetName() == nullptr || std::strlen( array->GetName() ) == 0 || !names.insert( array->GetName() ).second )
      throw std::invalid_argument( "Refinement requires unique, nonempty field-array names" );
    if( array->GetNumberOfTuples() != tuples || array->GetNumberOfComponents() < 1 )
      throw std::invalid_argument( "Refinement field tuple/component count mismatch: " + std::string( array->GetName() ) );
  }
}

template < typename Data > vtkSmartPointer< Data > allocate( Data & input, vtkIdType tuples, std::set< std::string > const & excluded )
{
  // VTK allocations index scalar values, not tuples. Prove component products
  // and addressable fixed storage before entering CopyAllocate/SetNumberOfTuples.
  if( tuples < 0 )
    throw std::invalid_argument( "Negative refinement field allocation count" );
  std::uint64_t totalBytes = 0;
  for( int i = 0; i < input.GetNumberOfArrays(); ++i )
  {
    auto * array = input.GetAbstractArray( i );
    if( excluded.count( array->GetName() ) || array == input.GetGlobalIds() || array == input.GetPedigreeIds() ||
        std::strcmp( array->GetName(), vtkDataSetAttributes::GhostArrayName() ) == 0 )
      continue;
    std::uint64_t const components = array->GetNumberOfComponents();
    std::uint64_t const count = static_cast< std::uint64_t >( tuples );
    if( count > static_cast< std::uint64_t >( std::numeric_limits< vtkIdType >::max() ) / components )
      throw std::overflow_error( "Refinement field tuple/component product exceeds vtkIdType" );
    std::uint64_t const values = count * components;
    std::uint64_t width = array->GetDataTypeSize();
    if( array->GetDataType() == VTK_STRING )
      width = sizeof( vtkStdString );
    if( array->GetDataType() == VTK_VARIANT )
      width = sizeof( vtkVariant );
    std::uint64_t bytes = values / 8 + ( values % 8 != 0 ); // Packed bits.
    if( width )
    {
      if( values > std::numeric_limits< std::size_t >::max() / width )
        throw std::overflow_error( "Refinement field storage exceeds addressable memory" );
      bytes = values * width;
    }
    else if( array->GetDataType() != VTK_BIT )
      throw std::invalid_argument( "Unsupported variable-sized refinement field storage" );
    if( bytes > std::numeric_limits< std::size_t >::max() - totalBytes )
      throw std::overflow_error( "Combined refinement field storage exceeds addressable memory" );
    totalBytes += bytes;
  }
  auto output = vtkSmartPointer< Data >::New();
  output->CopyAllOn();
  output->CopyGlobalIdsOff();
  output->CopyPedigreeIdsOff();
  output->CopyFieldOff( vtkDataSetAttributes::GhostArrayName() );
  for( int role : { vtkDataSetAttributes::GLOBALIDS, vtkDataSetAttributes::PEDIGREEIDS } )
    if( auto * array = input.GetAbstractAttribute( role ) )
      output->CopyFieldOff( array->GetName() );
  for( auto const & name : excluded )
    output->CopyFieldOff( name.c_str() );
  output->CopyAllocate( &input, tuples );
  for( int i = 0; i < output->GetNumberOfArrays(); ++i )
    output->GetAbstractArray( i )->SetNumberOfTuples( tuples );
  return output;
}

struct InterpolatePoints
{
  PointRegistry const & points;
  PointTransferPolicy policy;
  template < typename Source, typename Target > void operator()( Source * source, Target * target ) const
  {
    vtkDataArrayAccessor< Source > input( source );
    vtkDataArrayAccessor< Target > output( target );
    using Value = typename decltype( input )::APIType;
    if( policy != PointTransferPolicy::continuous && policy != PointTransferPolicy::equal && policy != PointTransferPolicy::nodeSet )
      throw std::invalid_argument( "Unknown refinement point transfer policy" );
    if( policy == PointTransferPolicy::continuous && !std::is_floating_point< Value >::value )
      throw std::invalid_argument( "Continuous refinement point fields must use floating point storage" );
    if( policy == PointTransferPolicy::nodeSet && ( !std::is_unsigned< Value >::value || source->GetNumberOfComponents() != 1 ) )
      throw std::invalid_argument( "Refinement node-set masks require unsigned scalar storage" );
    for( std::size_t i = 0; i < points.points().size(); ++i )
      for( int c = 0; c < source->GetNumberOfComponents(); ++c )
      {
        Connectivity const & support = points.points()[i].support;
        Value const first = input.Get( support.front(), c );
        Value value = first;
        if( i < static_cast< std::size_t >( points.originalSize() ) )
        {
          if constexpr( std::is_floating_point< Value >::value )
            if( policy == PointTransferPolicy::continuous && !std::isfinite( first ) )
              throw std::invalid_argument( "Nonfinite continuous refinement point field" );
          output.Set( i, c, first );
          continue;
        }
        if( policy == PointTransferPolicy::nodeSet )
        {
          // GEOS' nodeset reader treats exactly one as membership.
          if constexpr( std::is_unsigned< Value >::value )
          {
            bool member = true;
            for( vtkIdType corner : support )
              member = member && input.Get( corner, c ) == static_cast< Value >( 1 );
            value = static_cast< Value >( member );
          }
        }
        else if( policy == PointTransferPolicy::equal )
        {
          for( vtkIdType corner : support )
            if( !std::equal_to< Value >{}( input.Get( corner, c ), first ) )
              throw std::invalid_argument( "Differing categorical point values on refinement support" );
        }
        else
        {
          if constexpr( std::is_floating_point< Value >::value )
          {
            long double average = 0;
            for( vtkIdType corner : support )
            {
              Value const sample = input.Get( corner, c );
              if( !std::isfinite( sample ) )
                throw std::invalid_argument( "Nonfinite continuous refinement point field" );
              average += ( static_cast< long double >( sample ) - first ) / support.size();
            }
            value = static_cast< Value >( first + average );
            if( !std::isfinite( value ) )
              throw std::overflow_error( "Unrepresentable refinement point field" );
          }
        }
        output.Set( i, c, value );
      }
  }
};

struct ScaleExtensive
{
  std::vector< double > const & fractions;
  template < typename Array > void operator()( Array * array ) const
  {
    vtkDataArrayAccessor< Array > access( array );
    using Value = typename decltype( access )::APIType;
    if constexpr( !std::is_floating_point< Value >::value )
      throw std::invalid_argument( "Extensive refinement fields must use floating point storage" );
    else
      for( std::size_t i = 0; i < fractions.size(); ++i )
        for( int c = 0; c < array->GetNumberOfComponents(); ++c )
        {
          Value const value = static_cast< Value >( access.Get( i, c ) * fractions[i] );
          if( !std::isfinite( value ) )
            throw std::overflow_error( "Nonfinite or unrepresentable extensive refinement field" );
          access.Set( i, c, value );
        }
  }
};

} // namespace

vtkSmartPointer< vtkPointData > transferPointData( vtkPointData & input, PointRegistry const & points, TransferPolicies const & policies )
{
  validateArrays( input, points.originalSize() );
  for( auto const & policy : policies.pointArrays )
  {
    auto * array = input.GetAbstractArray( policy.first.c_str() );
    if( array == nullptr || array == input.GetGlobalIds() || array == input.GetPedigreeIds() ||
        policy.first == vtkDataSetAttributes::GhostArrayName() || policies.excludedPointArrays.count( policy.first ) )
      throw std::invalid_argument( "Unknown or relational refinement point-policy array: " + policy.first );
  }
  auto output = allocate( input, points.points().size(), policies.excludedPointArrays );
  for( int i = 0; i < output->GetNumberOfArrays(); ++i )
  {
    auto * target = vtkDataArray::SafeDownCast( output->GetAbstractArray( i ) );
    auto * source = input.GetArray( output->GetAbstractArray( i )->GetName() );
    if( target == nullptr || source == nullptr )
      throw std::invalid_argument( "Refinement point fields require supported numeric arrays" );
    auto const requested = policies.pointArrays.find( source->GetName() );
    PointTransferPolicy const policy =
        requested == policies.pointArrays.end()
            ? ( source->GetDataType() == VTK_FLOAT || source->GetDataType() == VTK_DOUBLE ? PointTransferPolicy::continuous
                                                                                          : PointTransferPolicy::equal )
            : requested->second;
    InterpolatePoints worker{ points, policy };
    if( !vtkArrayDispatch::Dispatch2BySameValueType< vtkArrayDispatch::AllTypes >::Execute( source, target, worker ) )
      throw std::invalid_argument( "Unsupported numeric refinement point-array implementation: " + std::string( source->GetName() ) );
  }
  return output;
}

vtkSmartPointer< vtkCellData > transferCellData( vtkCellData & input, vtkIdType parentCount, Connectivity const & parents,
                                                 std::vector< double > const & fractions, TransferPolicies const & policies )
{
  validateArrays( input, parentCount );
  if( parentCount < 0 || parents.size() != fractions.size() ||
      parents.size() > static_cast< std::size_t >( std::numeric_limits< vtkIdType >::max() ) )
    throw std::invalid_argument( "Refinement cell transfer count mismatch or overflow" );
  for( std::size_t i = 0; i < parents.size(); ++i )
    if( parents[i] < 0 || parents[i] >= parentCount || !std::isfinite( fractions[i] ) || fractions[i] <= 0 || fractions[i] > 1 )
      throw std::invalid_argument( "Invalid refinement cell parent or measured fraction" );
  std::map< vtkIdType, long double > sums;
  for( std::size_t i = 0; i < parents.size(); ++i )
    sums[parents[i]] += fractions[i];
  if( sums.size() != static_cast< std::size_t >( parentCount ) )
    throw std::invalid_argument( "Refinement cell transfer omits an active parent" );
  for( auto const & sum : sums )
    if( std::abs( sum.second - 1 ) > 1e-10L )
      throw std::invalid_argument( "Refinement measured child fractions do not sum to one" );
  auto output = allocate( input, parents.size(), policies.excludedCellArrays );
  for( std::size_t i = 0; i < parents.size(); ++i )
    output->CopyData( &input, parents[i], i );
  for( auto const & name : policies.extensiveCellArrays )
  {
    auto * array = output->GetArray( name.c_str() );
    if( array == nullptr )
      throw std::invalid_argument( "Unknown or nonnumeric extensive refinement cell array: " + name );
    ScaleExtensive worker{ fractions };
    if( !vtkArrayDispatch::Dispatch::Execute( array, worker ) )
      throw std::invalid_argument( "Unsupported extensive refinement array implementation: " + name );
  }
  return output;
}

PointFieldLayout::PointFieldLayout( vtkPointData & data ) : m_data( &data )
{
  validateArrays( data, data.GetNumberOfTuples() );
  for( int i = 0; i < data.GetNumberOfArrays(); ++i )
  {
    auto * array = data.GetArray( i );
    if( array == nullptr || array == data.GetGlobalIds() || array == data.GetPedigreeIds() )
      throw std::invalid_argument( "Point tuple layout must contain transferable numeric fields only" );
    m_arrays.push_back( array );
  }
  std::sort( m_arrays.begin(), m_arrays.end(),
             []( vtkDataArray * a, vtkDataArray * b ) { return std::strcmp( a->GetName(), b->GetName() ) < 0; } );
  appendInteger( m_schema, m_arrays.size() );
  for( auto * array : m_arrays )
  {
    appendString( m_schema, array->GetName() );
    appendInteger( m_schema, array->GetDataType() );
    appendInteger( m_schema, array->GetNumberOfComponents() );
    std::uint64_t roles = 0;
    for( int role = 0; role < vtkDataSetAttributes::NUM_ATTRIBUTES; ++role )
      if( data.GetAbstractAttribute( role ) == array )
        roles |= UINT64_C( 1 ) << role;
    appendInteger( m_schema, roles );
    for( int c = 0; c < array->GetNumberOfComponents(); ++c )
      appendString( m_schema, array->GetComponentName( c ) );
  }
}

std::vector< unsigned char > PointFieldLayout::pack( vtkIdType point ) const
{
  if( point < 0 || ( !m_arrays.empty() && point >= m_data->GetNumberOfTuples() ) )
    throw std::invalid_argument( "Refinement point-field tuple outside array" );
  Bytes bytes = m_schema;
  PackTuple worker{ point, bytes };
  for( auto * array : m_arrays )
    if( !vtkArrayDispatch::Dispatch::Execute( array, worker ) )
      throw std::invalid_argument( "Unsupported refinement point tuple array" );
  return bytes;
}

void PointFieldLayout::install( vtkIdType point, std::vector< unsigned char > const & tuple ) const
{
  if( point < 0 || ( !m_arrays.empty() && point >= m_data->GetNumberOfTuples() ) )
    throw std::invalid_argument( "Refinement point-field tuple outside array" );
  std::size_t expected = m_schema.size();
  for( auto * array : m_arrays )
  {
    std::size_t const width = static_cast< std::size_t >( array->GetDataTypeSize() ) * array->GetNumberOfComponents();
    if( width > std::numeric_limits< std::size_t >::max() - expected )
      throw std::overflow_error( "Refinement point tuple length overflow" );
    expected += width;
  }
  if( tuple.size() != expected || !std::equal( m_schema.begin(), m_schema.end(), tuple.begin() ) )
    throw std::invalid_argument( "Refinement point-field schema or tuple length mismatch" );
  std::size_t cursor = m_schema.size();
  InstallTuple worker{ point, tuple, cursor };
  for( auto * array : m_arrays )
    if( !vtkArrayDispatch::Dispatch::Execute( array, worker ) )
      throw std::invalid_argument( "Unsupported refinement point tuple array" );
}

CellFieldLayout::CellFieldLayout( vtkCellData & data ) : m_data( &data )
{
  validateArrays( data, data.GetNumberOfTuples() );
  for( int i = 0; i < data.GetNumberOfArrays(); ++i )
  {
    auto * array = data.GetAbstractArray( i );
    if( array == data.GetGlobalIds() || array == data.GetPedigreeIds() ||
        ( !vtkDataArray::SafeDownCast( array ) && !vtkStringArray::SafeDownCast( array ) ) )
      throw std::invalid_argument( "Cell tuple layout requires transferable numeric/string fields" );
    if( !vtkStringArray::SafeDownCast( array ) && !vtkBitArray::SafeDownCast( array ) &&
        ( array->GetDataTypeSize() <= 0 ||
          !vtkArrayDispatch::Dispatch::Execute( vtkDataArray::SafeDownCast( array ), []( auto * ) {} ) ) )
      throw std::invalid_argument( "Unsupported numeric refinement cell tuple storage" );
    m_arrays.push_back( array );
  }
  std::sort( m_arrays.begin(), m_arrays.end(), []( auto * a, auto * b ) { return std::strcmp( a->GetName(), b->GetName() ) < 0; } );
  appendInteger( m_schema, m_arrays.size() );
  for( auto * array : m_arrays )
  {
    appendString( m_schema, array->GetName() );
    appendInteger( m_schema, array->GetDataType() );
    appendInteger( m_schema, array->GetNumberOfComponents() );
    std::uint64_t roles = 0;
    for( int role = 0; role < vtkDataSetAttributes::NUM_ATTRIBUTES; ++role )
      if( data.GetAbstractAttribute( role ) == array )
        roles |= UINT64_C( 1 ) << role;
    appendInteger( m_schema, roles );
    for( int c = 0; c < array->GetNumberOfComponents(); ++c )
      appendString( m_schema, array->GetComponentName( c ) );
  }
}

std::vector< unsigned char > CellFieldLayout::pack( vtkIdType cell ) const
{
  if( cell < 0 || ( !m_arrays.empty() && cell >= m_data->GetNumberOfTuples() ) )
    throw std::invalid_argument( "Refinement cell tuple outside array" );
  Bytes bytes = m_schema;
  for( auto * array : m_arrays )
  {
    if( auto * strings = vtkStringArray::SafeDownCast( array ) )
      for( int c = 0; c < array->GetNumberOfComponents(); ++c )
      {
        auto const & value = strings->GetValue( cell * array->GetNumberOfComponents() + c );
        appendInteger( bytes, value.size() );
        bytes.insert( bytes.end(), value.begin(), value.end() );
      }
    else if( auto * bits = vtkBitArray::SafeDownCast( array ) )
      for( int c = 0; c < array->GetNumberOfComponents(); ++c )
        bytes.push_back( static_cast< unsigned char >( bits->GetValue( cell * array->GetNumberOfComponents() + c ) ) );
    else
    {
      PackTuple worker{ cell, bytes };
      if( !vtkArrayDispatch::Dispatch::Execute( vtkDataArray::SafeDownCast( array ), worker ) )
        throw std::invalid_argument( "Unsupported numeric refinement cell tuple array" );
    }
  }
  return bytes;
}

void CellFieldLayout::install( vtkIdType cell, std::vector< unsigned char > const & tuple ) const
{
  if( cell < 0 || ( !m_arrays.empty() && cell >= m_data->GetNumberOfTuples() ) )
    throw std::invalid_argument( "Refinement cell tuple outside array" );
  if( tuple.size() < m_schema.size() || !std::equal( m_schema.begin(), m_schema.end(), tuple.begin() ) )
    throw std::invalid_argument( "Refinement cell-field schema mismatch" );
  auto length = [&]( std::size_t & cursor )
  {
    if( tuple.size() - cursor < 8 )
      throw std::invalid_argument( "Truncated refinement string length" );
    std::uint64_t value = 0;
    for( int i = 0; i < 8; ++i )
      value |= static_cast< std::uint64_t >( tuple[cursor++] ) << ( 8 * i );
    if( value > tuple.size() - cursor )
      throw std::invalid_argument( "Truncated refinement string value" );
    return static_cast< std::size_t >( value );
  };
  // Validate every length before changing any tuple, including variable strings.
  std::size_t cursor = m_schema.size();
  for( auto * array : m_arrays )
  {
    if( vtkStringArray::SafeDownCast( array ) )
      for( int c = 0; c < array->GetNumberOfComponents(); ++c )
      {
        auto const n = length( cursor );
        cursor += n;
      }
    else
    {
      std::size_t const width = vtkBitArray::SafeDownCast( array ) ? 1 : array->GetDataTypeSize();
      if( static_cast< std::size_t >( array->GetNumberOfComponents() ) > ( tuple.size() - cursor ) / width )
        throw std::invalid_argument( "Truncated refinement numeric cell tuple" );
      if( vtkBitArray::SafeDownCast( array ) )
        for( int c = 0; c < array->GetNumberOfComponents(); ++c )
          if( tuple[cursor + c] > 1 )
            throw std::invalid_argument( "Invalid refinement bit value" );
      cursor += width * array->GetNumberOfComponents();
    }
  }
  if( cursor != tuple.size() )
    throw std::invalid_argument( "Trailing refinement cell tuple bytes" );
  cursor = m_schema.size();
  for( auto * array : m_arrays )
    if( auto * strings = vtkStringArray::SafeDownCast( array ) )
      for( int c = 0; c < array->GetNumberOfComponents(); ++c )
      {
        auto const n = length( cursor );
        strings->SetValue( cell * array->GetNumberOfComponents() + c,
                           std::string( reinterpret_cast< char const * >( tuple.data() + cursor ), n ) );
        cursor += n;
      }
    else if( auto * bits = vtkBitArray::SafeDownCast( array ) )
      for( int c = 0; c < array->GetNumberOfComponents(); ++c )
        bits->SetValue( cell * array->GetNumberOfComponents() + c, tuple[cursor++] );
    else
    {
      InstallTuple worker{ cell, tuple, cursor };
      if( !vtkArrayDispatch::Dispatch::Execute( vtkDataArray::SafeDownCast( array ), worker ) )
        throw std::invalid_argument( "Unsupported numeric refinement cell tuple array" );
    }
}

} // namespace geos::vtk::refinement
