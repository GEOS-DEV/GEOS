#include "BufferOpsDevice.hpp"

namespace geosx
{

namespace bufferOps
{

template< bool DO_PACKING, typename T >
GEOSX_HOST_DEVICE
typename std::enable_if< can_memcpy< T >, localIndex >::type
PackDevice( buffer_unit_type * & buffer,
            T const & var )
{
  localIndex const sizeOfPackedChars = sizeof(T);
  if( DO_PACKING )
  {
    memcpy( buffer, &var, sizeOfPackedChars );
    buffer += sizeOfPackedChars;
  }
  return sizeOfPackedChars;
}


template< bool DO_PACKING, typename T, typename INDEX_TYPE >
GEOSX_HOST_DEVICE
typename std::enable_if< can_memcpy< T >, localIndex >::type
PackRawPointerDevice( buffer_unit_type * & buffer,
                      T const * const restrict var,
                      INDEX_TYPE const length )
{
  localIndex const sizeOfPackedChars = length * sizeof(T);
  if( DO_PACKING )
  {
    memcpy( buffer, var, sizeOfPackedChars );
    buffer += sizeOfPackedChars;
  }
  return sizeOfPackedChars;
}

template< bool DO_PACKING, typename T >
GEOSX_HOST_DEVICE
typename std::enable_if< can_memcpy< T >, localIndex >::type
PackDevice( buffer_unit_type * & buffer,
            LvArray::ArraySlice< T, 1 > const & var )
{
  localIndex length = var.size();
  localIndex sizeOfPackedChars = PackDevice< DO_PACKING >( buffer, length );
  sizeOfPackedChars += length * sizeof(T);
  if( DO_PACKING )
  {
    for( localIndex ii = 0 ; ii < length ; ++ii )
    {
      memcpy( buffer, &var[ii], sizeof(T));
      buffer += sizeof(T);
    }
    buffer += length * sizeof(T);
  }
  return sizeOfPackedChars;
}

template< bool DO_PACKING, typename T, typename INDEX_TYPE >
GEOSX_HOST_DEVICE
typename std::enable_if< can_memcpy< T >, localIndex >::type
PackPointerDevice( buffer_unit_type * & buffer,
                   T const * const restrict var,
                   INDEX_TYPE const length )
{
  localIndex sizeOfPackedChars = PackDevice< DO_PACKING >( buffer, length );
  sizeOfPackedChars += length * sizeof(T);
  if( DO_PACKING )
  {
    memcpy( buffer, var, length*sizeof(T));
    buffer += length * sizeof(T);
  }
  return sizeOfPackedChars;
}

template< bool DO_PACKING, typename T, typename INDEX_TYPE >
GEOSX_HOST_DEVICE
typename std::enable_if< !can_memcpy< T >, localIndex >::type
PackPointerDevice( buffer_unit_type * & buffer,
                   T const * const restrict var,
                   INDEX_TYPE const length )
{
  localIndex sizeOfPackedChars = PackDevice< DO_PACKING >( buffer, length );
  for( INDEX_TYPE ii = 0 ; ii < length ; ++ii )
    sizeOfPackedChars += PackDevice< DO_PACKING >( buffer, var[ii] );
  return sizeOfPackedChars;
}

template< typename T >
GEOSX_HOST_DEVICE
typename std::enable_if< can_memcpy< T >, localIndex >::type
UnpackDevice( buffer_unit_type const * & buffer,
              T & var )
{
  localIndex const sizeOfUnpackedChars = sizeof(T);
  memcpy( &var, buffer, sizeOfUnpackedChars );
  buffer += sizeOfUnpackedChars;
  return sizeOfUnpackedChars;
}

template< typename T, int NDIM >
GEOSX_HOST_DEVICE
typename std::enable_if< can_memcpy< T >, localIndex >::type
UnpackDevice( buffer_unit_type const * & buffer,
              LvArray::ArraySlice< T, NDIM > const & var )
{
  localIndex length;
  localIndex sizeOfPackedChars = UnpackDevice( buffer, length );
  T const * const buffer_T = reinterpret_cast< T * >( buffer );
  for( localIndex ii = 0 ; ii < length ; ++ii )
    var[ii] = buffer_T[ ii ];
  buffer += length * sizeof(T);
  sizeOfPackedChars += length * sizeof(T);
  return sizeOfPackedChars;
}

template< typename T, typename INDEX_TYPE >
GEOSX_HOST_DEVICE
typename std::enable_if< can_memcpy< T >, localIndex >::type
UnpackPointerDevice( buffer_unit_type const * & buffer,
                     T * const restrict var,
                     INDEX_TYPE const expectedLength )
{
  INDEX_TYPE length = 0;
  localIndex sizeOfUnpackedChars = UnpackDevice( buffer, length );
  GEOSX_ASSERT_MSG( length == expectedLength, "expectedLength != length: " <<
                    expectedLength << " != " << length );
  GEOSX_DEBUG_VAR( expectedLength );
  memcpy( var, buffer, length * sizeof(T) );
  sizeOfUnpackedChars += length * sizeof(T);
  buffer += length * sizeof(T);
  return sizeOfUnpackedChars;
}

template< typename T, typename INDEX_TYPE >
GEOSX_HOST_DEVICE
typename std::enable_if< !can_memcpy< T >, localIndex >::type
UnpackPointerDevice( buffer_unit_type const * & buffer,
                     T * const restrict var,
                     INDEX_TYPE const expectedLength )
{
  localIndex sizeOfUnpackedChars = 0;
  INDEX_TYPE length = 0;
  UnpackDevice( buffer, length );
  GEOSX_ASSERT_MSG( length == expectedLength, "expectedLength != length: " <<
                    expectedLength << " != " << length );
  GEOSX_DEBUG_VAR( expectedLength );
  for( INDEX_TYPE ii = 0 ; ii < length ; ++ii )
    sizeOfUnpackedChars += UnpackDevice( buffer, var[ii] );
  return sizeOfUnpackedChars;
}

template< typename T, typename INDEX_TYPE >
GEOSX_HOST_DEVICE
typename std::enable_if< can_memcpy< T >, localIndex >::type
UnpackRawPointerDevice( buffer_unit_type const * & buffer,
                        T * const restrict var,
                        INDEX_TYPE const length )
{
  localIndex sizeOfUnpackedChars = 0;
  memcpy( var, buffer, length * sizeof(T) );
  sizeOfUnpackedChars += length * sizeof(T);
  buffer += length * sizeof(T);
  return sizeOfUnpackedChars;
}

template< typename POLICY, bool DO_PACKING, typename T, int NDIM, int UNIT_STRIDE_DIM, typename INDEX_TYPE >
typename std::enable_if< can_memcpy< T >, localIndex >::type
PackDevice( buffer_unit_type * & buffer,
            LvArray::ArrayView< T, NDIM, UNIT_STRIDE_DIM, INDEX_TYPE > const & var )
{
  localIndex packedSize = PackPointerDevice< DO_PACKING >( buffer, var.dims(), NDIM );
  packedSize += PackPointerDevice< DO_PACKING >( buffer, var.strides(), NDIM );
  if( DO_PACKING )
  {
    forall_in_range< POLICY >( 0, var.size(), GEOSX_HOST_DEVICE_LAMBDA( localIndex ii )
        {
          reinterpret_cast< std::remove_const_t< T > * >( buffer )[ ii ] = var.data()[ ii ];
        } );
  }
  packedSize += var.size() * sizeof(T);
  if( DO_PACKING )
  {
    buffer += var.size() * sizeof(T);
  }
  return packedSize;
}

template< typename POLICY, typename T, int NDIM, int UNIT_STRIDE_DIM, typename INDEX_TYPE >
typename std::enable_if< can_memcpy< T >, localIndex >::type
UnpackDevice( buffer_unit_type const * & buffer,
              LvArray::ArrayView< T, NDIM, UNIT_STRIDE_DIM, INDEX_TYPE > const & var )
{
  INDEX_TYPE dims[NDIM];
  localIndex packedSize = UnpackPointerDevice( buffer, dims, NDIM );
  INDEX_TYPE strides[NDIM];
  packedSize += UnpackPointerDevice( buffer, strides, NDIM );
  for( int dd = 0 ; dd < NDIM ; ++dd )
  {
    GEOSX_ASSERT_MSG( strides[dd] == var.strides()[dd], "Strides are inconsistent: " <<
                      strides[dd] << " != " << var.strides()[dd] );
  }
  forall_in_range< POLICY >( 0, var.size(), GEOSX_HOST_DEVICE_LAMBDA( localIndex ii )
      {
        var.data()[ ii ] = reinterpret_cast< const T * >( buffer )[ ii ];
      } );
  packedSize += var.size() * sizeof(T);
  buffer += var.size() * sizeof(T);
  return packedSize;
}

template< typename POLICY, bool DO_PACKING, typename T, int NDIM, typename T_INDICES >
typename std::enable_if< can_memcpy< T >, localIndex >::type
PackByIndexDevice ( buffer_unit_type * & buffer,
                    LvArray::ArrayView< T, NDIM > const & var,
                    const T_INDICES & indices )
{
  localIndex numIndices = indices.size();
  localIndex packedSize = PackPointerDevice< DO_PACKING >( buffer, var.strides(), NDIM );
  buffer_unit_type * const devBuffer = buffer;
  localIndex unitSize = var.strides()[0] * sizeof(T);
  if( DO_PACKING )
  {
    forall_in_range< POLICY >( 0, numIndices, GEOSX_HOST_DEVICE_LAMBDA( localIndex ii )
        {
          localIndex threadOffset = ii * unitSize;
          buffer_unit_type * threadBuffer = devBuffer + threadOffset;
          T const * const data = var.data( indices[ii] );
          PackRawPointerDevice< DO_PACKING >( threadBuffer, data, var.strides()[0] );
        } );
  }
  localIndex avPackedSize = numIndices * unitSize;
  packedSize += avPackedSize;
  if( DO_PACKING )
    buffer += avPackedSize;
  return packedSize;
}

template< typename POLICY, typename T, int NDIM, typename T_INDICES >
typename std::enable_if< can_memcpy< T >, localIndex >::type
UnpackByIndexDevice ( buffer_unit_type const * & buffer,
                      LvArray::ArrayView< T, NDIM > const & var,
                      T_INDICES const & indices )
{
  localIndex strides[NDIM];
  localIndex sizeOfPackedChars = UnpackPointerDevice( buffer, strides, NDIM );
  for( int dd = 0 ; dd < NDIM ; ++dd )
  {
    GEOSX_ASSERT_MSG( strides[dd] == var.strides()[dd], "Strides are inconsistent: " <<
                      strides[dd] << " != " << var.strides()[dd] );
  }
  localIndex numIndices = indices.size();
  buffer_unit_type const * devBuffer = buffer;
  localIndex unitSize = var.strides()[0] * sizeof(T);
  forall_in_range< POLICY >( 0, numIndices, GEOSX_HOST_DEVICE_LAMBDA( localIndex ii )
      {
        localIndex threadOffset = ii * unitSize;
        buffer_unit_type const * threadBuffer = devBuffer + threadOffset;
        T * const data = (var.data( indices[ii] ));
        UnpackRawPointerDevice( threadBuffer, data, var.strides()[0] );
      } );
  localIndex avSize = numIndices * unitSize;
  sizeOfPackedChars += avSize;
  buffer += avSize;
  return sizeOfPackedChars;
}

#define pd( DO_PACKING, TYPE, NDIM, UNIT_STRIDE_DIM ) \
  template localIndex PackDevice< parallelDevicePolicy< >, DO_PACKING, TYPE, NDIM, UNIT_STRIDE_DIM, localIndex > \
    ( buffer_unit_type * &buffer, \
    LvArray::ArrayView< TYPE, NDIM, UNIT_STRIDE_DIM, localIndex > const & var );

#define upd( TYPE, NDIM, UNIT_STRIDE_DIM ) \
  template localIndex UnpackDevice< parallelDevicePolicy< >, TYPE, NDIM, UNIT_STRIDE_DIM, localIndex > \
    ( buffer_unit_type const * & buffer, \
    LvArray::ArrayView< TYPE, NDIM, UNIT_STRIDE_DIM, localIndex > const & var );

#define pbid( DO_PACKING, TYPE, NDIM ) \
  template localIndex PackByIndexDevice< parallelDevicePolicy< >, DO_PACKING, TYPE, NDIM > \
    ( buffer_unit_type * &buffer, \
    LvArray::ArrayView< TYPE, NDIM > const & var, \
    LvArray::ArrayView< const localIndex, 1, 0, localIndex > const & indices );

#define upbid( TYPE, NDIM ) \
  template localIndex UnpackByIndexDevice< parallelDevicePolicy< >, TYPE, NDIM > \
    ( buffer_unit_type const * & buffer, \
    LvArray::ArrayView< TYPE, NDIM > const & var, \
    LvArray::ArrayView< const localIndex, 1, 0, localIndex > const & indices );

// 1d native types
pd( false, const int, 1, 0 )
pd( false, const localIndex, 1, 0 )
pd( false, const globalIndex, 1, 0 )
pd( false, const real32, 1, 0 )
pd( false, const real64, 1, 0 )
// 2d native types
pd( false, const int, 2, 0 )
pd( false, const localIndex, 2, 0 )
pd( false, const globalIndex, 2, 0 )
pd( false, const real32, 2, 0 )
pd( false, const real64, 2, 0 )

pd( false, const int, 2, 1 )
pd( false, const localIndex, 2, 1 )
pd( false, const globalIndex, 2, 1 )
pd( false, const real32, 2, 1 )
pd( false, const real64, 2, 1 )
// 3d native types
pd( false, const int, 3, 0 )
pd( false, const localIndex, 3, 0 )
pd( false, const globalIndex, 3, 0 )
pd( false, const real32, 3, 0 )
pd( false, const real64, 3, 0 )

pd( false, const int, 3, 1 )
pd( false, const localIndex, 3, 1 )
pd( false, const globalIndex, 3, 1 )
pd( false, const real32, 3, 1 )
pd( false, const real64, 3, 1 )

pd( false, const int, 3, 2 )
pd( false, const localIndex, 3, 2 )
pd( false, const globalIndex, 3, 2 )
pd( false, const real32, 3, 2 )
pd( false, const real64, 3, 2 )

// 4d native types
pd( false, const int, 4, 0 )
pd( false, const localIndex, 4, 0 )
pd( false, const globalIndex, 4, 0 )
pd( false, const real32, 4, 0 )
pd( false, const real64, 4, 0 )

pd( false, const int, 4, 1 )
pd( false, const localIndex, 4, 1 )
pd( false, const globalIndex, 4, 1 )
pd( false, const real32, 4, 1 )
pd( false, const real64, 4, 1 )

pd( false, const int, 4, 2 )
pd( false, const localIndex, 4, 2 )
pd( false, const globalIndex, 4, 2 )
pd( false, const real32, 4, 2 )
pd( false, const real64, 4, 2 )

pd( false, const int, 4, 3 )
pd( false, const localIndex, 4, 3 )
pd( false, const globalIndex, 4, 3 )
pd( false, const real32, 4, 3 )
pd( false, const real64, 4, 3 )
// 5d native types
pd( false, const int, 5, 0 )
pd( false, const localIndex, 5, 0 )
pd( false, const globalIndex, 5, 0 )
pd( false, const real32, 5, 0 )
pd( false, const real64, 5, 0 )

pd( false, const int, 5, 1 )
pd( false, const localIndex, 5, 1 )
pd( false, const globalIndex, 5, 1 )
pd( false, const real32, 5, 1 )
pd( false, const real64, 5, 1 )

pd( false, const int, 5, 2 )
pd( false, const localIndex, 5, 2 )
pd( false, const globalIndex, 5, 2 )
pd( false, const real32, 5, 2 )
pd( false, const real64, 5, 2 )

pd( false, const int, 5, 3 )
pd( false, const localIndex, 5, 3 )
pd( false, const globalIndex, 5, 3 )
pd( false, const real32, 5, 3 )
pd( false, const real64, 5, 3 )

pd( false, const int, 5, 4 )
pd( false, const localIndex, 5, 4 )
pd( false, const globalIndex, 5, 4 )
pd( false, const real32, 5, 4 )
pd( false, const real64, 5, 4 )
// 1d tensor-types
pd( false, const R1Tensor, 1, 0 )
pd( false, const R2Tensor, 1, 0 )
pd( false, const R2SymTensor, 1, 0 )
// 2d tensor-types
pd( false, const R1Tensor, 2, 0 )
pd( false, const R2Tensor, 2, 0 )
pd( false, const R2SymTensor, 2, 0 )

pd( false, const R1Tensor, 2, 1 )
pd( false, const R2Tensor, 2, 1 )
pd( false, const R2SymTensor, 2, 1 )

// 3d tensor-types
pd( false, const R1Tensor, 3, 0 )
pd( false, const R2Tensor, 3, 0 )
pd( false, const R2SymTensor, 3, 0 )

pd( false, const R1Tensor, 3, 1 )
pd( false, const R2Tensor, 3, 1 )
pd( false, const R2SymTensor, 3, 1 )

pd( false, const R1Tensor, 3, 2 )
pd( false, const R2Tensor, 3, 2 )
pd( false, const R2SymTensor, 3, 2 )
// 1d native types
pd( true, const int, 1, 0 )
pd( true, const localIndex, 1, 0 )
pd( true, const globalIndex, 1, 0 )
pd( true, const real32, 1, 0 )
pd( true, const real64, 1, 0 )
// 2d native types
pd( true, const int, 2, 0 )
pd( true, const localIndex, 2, 0 )
pd( true, const globalIndex, 2, 0 )
pd( true, const real32, 2, 0 )
pd( true, const real64, 2, 0 )

pd( true, const int, 2, 1 )
pd( true, const localIndex, 2, 1 )
pd( true, const globalIndex, 2, 1 )
pd( true, const real32, 2, 1 )
pd( true, const real64, 2, 1 )
// 3d native types
pd( true, const int, 3, 0 )
pd( true, const localIndex, 3, 0 )
pd( true, const globalIndex, 3, 0 )
pd( true, const real32, 3, 0 )
pd( true, const real64, 3, 0 )

pd( true, const int, 3, 1 )
pd( true, const localIndex, 3, 1 )
pd( true, const globalIndex, 3, 1 )
pd( true, const real32, 3, 1 )
pd( true, const real64, 3, 1 )

pd( true, const int, 3, 2 )
pd( true, const localIndex, 3, 2 )
pd( true, const globalIndex, 3, 2 )
pd( true, const real32, 3, 2 )
pd( true, const real64, 3, 2 )
// 4d native types
pd( true, const int, 4, 0 )
pd( true, const localIndex, 4, 0 )
pd( true, const globalIndex, 4, 0 )
pd( true, const real32, 4, 0 )
pd( true, const real64, 4, 0 )

pd( true, const int, 4, 1 )
pd( true, const localIndex, 4, 1 )
pd( true, const globalIndex, 4, 1 )
pd( true, const real32, 4, 1 )
pd( true, const real64, 4, 1 )

pd( true, const int, 4, 2 )
pd( true, const localIndex, 4, 2 )
pd( true, const globalIndex, 4, 2 )
pd( true, const real32, 4, 2 )
pd( true, const real64, 4, 2 )

pd( true, const int, 4, 3 )
pd( true, const localIndex, 4, 3 )
pd( true, const globalIndex, 4, 3 )
pd( true, const real32, 4, 3 )
pd( true, const real64, 4, 3 )
// 5d native types
pd( true, const int, 5, 0 )
pd( true, const localIndex, 5, 0 )
pd( true, const globalIndex, 5, 0 )
pd( true, const real32, 5, 0 )
pd( true, const real64, 5, 0 )

pd( true, const int, 5, 1 )
pd( true, const localIndex, 5, 1 )
pd( true, const globalIndex, 5, 1 )
pd( true, const real32, 5, 1 )
pd( true, const real64, 5, 1 )

pd( true, const int, 5, 2 )
pd( true, const localIndex, 5, 2 )
pd( true, const globalIndex, 5, 2 )
pd( true, const real32, 5, 2 )
pd( true, const real64, 5, 2 )

pd( true, const int, 5, 3 )
pd( true, const localIndex, 5, 3 )
pd( true, const globalIndex, 5, 3 )
pd( true, const real32, 5, 3 )
pd( true, const real64, 5, 3 )

pd( true, const int, 5, 4 )
pd( true, const localIndex, 5, 4 )
pd( true, const globalIndex, 5, 4 )
pd( true, const real32, 5, 4 )
pd( true, const real64, 5, 4 )
// 1d tensor-types
pd( true, const R1Tensor, 1, 0 )
pd( true, const R2Tensor, 1, 0 )
pd( true, const R2SymTensor, 1, 0 )
// 2d tensor-types
pd( true, const R1Tensor, 2, 0 )
pd( true, const R2Tensor, 2, 0 )
pd( true, const R2SymTensor, 2, 0 )

pd( true, const R1Tensor, 2, 1 )
pd( true, const R2Tensor, 2, 1 )
pd( true, const R2SymTensor, 2, 1 )

// 1d native types
upd( int, 1, 0 )
upd( localIndex, 1, 0 )
upd( globalIndex, 1, 0 )
upd( real32, 1, 0 )
upd( real64, 1, 0 )
// 2d native types
upd( int, 2, 0 )
upd( localIndex, 2, 0 )
upd( globalIndex, 2, 0 )
upd( real32, 2, 0 )
upd( real64, 2, 0 )

upd( int, 2, 1 )
upd( localIndex, 2, 1 )
upd( globalIndex, 2, 1 )
upd( real32, 2, 1 )
upd( real64, 2, 1 )

// 3d native types
upd( int, 3, 0 )
upd( localIndex, 3, 0 )
upd( globalIndex, 3, 0 )
upd( real32, 3, 0 )
upd( real64, 3, 0 )

upd( int, 3, 1 )
upd( localIndex, 3, 1 )
upd( globalIndex, 3, 1 )
upd( real32, 3, 1 )
upd( real64, 3, 1 )

upd( int, 3, 2 )
upd( localIndex, 3, 2 )
upd( globalIndex, 3, 2 )
upd( real32, 3, 2 )
upd( real64, 3, 2 )

// 4d native types
upd( int, 4, 0 )
upd( localIndex, 4, 0 )
upd( globalIndex, 4, 0 )
upd( real32, 4, 0 )
upd( real64, 4, 0 )

upd( int, 4, 1 )
upd( localIndex, 4, 1 )
upd( globalIndex, 4, 1 )
upd( real32, 4, 1 )
upd( real64, 4, 1 )

upd( int, 4, 2 )
upd( localIndex, 4, 2 )
upd( globalIndex, 4, 2 )
upd( real32, 4, 2 )
upd( real64, 4, 2 )

upd( int, 4, 3 )
upd( localIndex, 4, 3 )
upd( globalIndex, 4, 3 )
upd( real32, 4, 3 )
upd( real64, 4, 3 )

// 5d native types
upd( int, 5, 0 )
upd( localIndex, 5, 0 )
upd( globalIndex, 5, 0 )
upd( real32, 5, 0 )
upd( real64, 5, 0 )

upd( int, 5, 1 )
upd( localIndex, 5, 1 )
upd( globalIndex, 5, 1 )
upd( real32, 5, 1 )
upd( real64, 5, 1 )

upd( int, 5, 2 )
upd( localIndex, 5, 2 )
upd( globalIndex, 5, 2 )
upd( real32, 5, 2 )
upd( real64, 5, 2 )

upd( int, 5, 3 )
upd( localIndex, 5, 3 )
upd( globalIndex, 5, 3 )
upd( real32, 5, 3 )
upd( real64, 5, 3 )

upd( int, 5, 4 )
upd( localIndex, 5, 4 )
upd( globalIndex, 5, 4 )
upd( real32, 5, 4 )
upd( real64, 5, 4 )
// 1d tensor-types
upd( R1Tensor, 1, 0 )
upd( R2Tensor, 1, 0 )
upd( R2SymTensor, 1, 0 )
// 2d tensor-types
upd( R1Tensor, 2, 0 )
upd( R2Tensor, 2, 0 )
upd( R2SymTensor, 2, 0 )

upd( R1Tensor, 2, 1 )
upd( R2Tensor, 2, 1 )
upd( R2SymTensor, 2, 1 )
// 3d tensor-types
upd( R1Tensor, 3, 0 )
upd( R2Tensor, 3, 0 )
upd( R2SymTensor, 3, 0 )

upd( R1Tensor, 3, 1 )
upd( R2Tensor, 3, 1 )
upd( R2SymTensor, 3, 1 )

upd( R1Tensor, 3, 2 )
upd( R2Tensor, 3, 2 )
upd( R2SymTensor, 3, 2 )
// 1d native types
pbid( false, const int, 1 )
pbid( false, const localIndex, 1 )
pbid( false, const globalIndex, 1 )
pbid( false, const real32, 1 )
pbid( false, const real64, 1 )
// 2d native types
pbid( false, const int, 2 )
pbid( false, const localIndex, 2 )
pbid( false, const globalIndex, 2 )
pbid( false, const real32, 2 )
pbid( false, const real64, 2 )
// 3d native types
pbid( false, const int, 3 )
pbid( false, const localIndex, 3 )
pbid( false, const globalIndex, 3 )
pbid( false, const real32, 3 )
pbid( false, const real64, 3 )
// 1d tensor-types
pbid( false, const R1Tensor, 1 )
pbid( false, const R2Tensor, 1 )
pbid( false, const R2SymTensor, 1 )
// 2d tensor-types
pbid( false, const R1Tensor, 2 )
pbid( false, const R2Tensor, 2 )
pbid( false, const R2SymTensor, 2 )
// 3d tensor-types
pbid( false, const R1Tensor, 3 )
pbid( false, const R2Tensor, 3 )
pbid( false, const R2SymTensor, 3 )
// 4d native types
pbid( false, const int, 4 )
pbid( false, const localIndex, 4 )
pbid( false, const globalIndex, 4 )
pbid( false, const real32, 4 )
pbid( false, const real64, 4 )
// 5d native types
pbid( false, const int, 5 )
pbid( false, const localIndex, 5 )
pbid( false, const globalIndex, 5 )
pbid( false, const real32, 5 )
pbid( false, const real64, 5 )

// 1d native types
pbid( true, const int, 1 )
pbid( true, const localIndex, 1 )
pbid( true, const globalIndex, 1 )
pbid( true, const real32, 1 )
pbid( true, const real64, 1 )
// 2d native types
pbid( true, const int, 2 )
pbid( true, const localIndex, 2 )
pbid( true, const globalIndex, 2 )
pbid( true, const real32, 2 )
pbid( true, const real64, 2 )
// 3d native types
pbid( true, const int, 3 )
pbid( true, const localIndex, 3 )
pbid( true, const globalIndex, 3 )
pbid( true, const real32, 3 )
pbid( true, const real64, 3 )
// 1d tensor-types
pbid( true, const R1Tensor, 1 )
pbid( true, const R2Tensor, 1 )
pbid( true, const R2SymTensor, 1 )
// 2d tensor-types
pbid( true, const R1Tensor, 2 )
pbid( true, const R2Tensor, 2 )
pbid( true, const R2SymTensor, 2 )
// 3d tensor-types
pbid( true, const R1Tensor, 3 )
pbid( true, const R2Tensor, 3 )
pbid( true, const R2SymTensor, 3 )
// 4d native types
pbid( true, const int, 4 )
pbid( true, const localIndex, 4 )
pbid( true, const globalIndex, 4 )
pbid( true, const real32, 4 )
pbid( true, const real64, 4 )
// 5d native types
pbid( true, const int, 5 )
pbid( true, const localIndex, 5 )
pbid( true, const globalIndex, 5 )
pbid( true, const real32, 5 )
pbid( true, const real64, 5 )

// 1d native types
upbid( int, 1 )
upbid( localIndex, 1 )
upbid( globalIndex, 1 )
upbid( real32, 1 )
upbid( real64, 1 )
// 2d native types
upbid( int, 2 )
upbid( localIndex, 2 )
upbid( globalIndex, 2 )
upbid( real32, 2 )
upbid( real64, 2 )
// 3d native types
upbid( int, 3 )
upbid( localIndex, 3 )
upbid( globalIndex, 3 )
upbid( real32, 3 )
upbid( real64, 3 )
// 4d native types
upbid( int, 4 )
upbid( localIndex, 4 )
upbid( globalIndex, 4 )
upbid( real32, 4 )
upbid( real64, 4 )
// 5d native types
upbid( int, 5 )
upbid( localIndex, 5 )
upbid( globalIndex, 5 )
upbid( real32, 5 )
upbid( real64, 5 )
// 1d tensor-types
upbid( R1Tensor, 1 )
upbid( R2Tensor, 1 )
upbid( R2SymTensor, 1 )
// 2d tensor-types
upbid( R1Tensor, 2 )
upbid( R2Tensor, 2 )
upbid( R2SymTensor, 2 )

// 3d tensor-types
upbid( R1Tensor, 3 )
upbid( R2Tensor, 3 )
upbid( R2SymTensor, 3 )
}

}
