/*
 * packer_inline.hpp
 *
 *  Created on: Jan 5, 2018
 *      Author: settgast
 */

#ifndef SRC_COMPONENTS_CORE_SRC_MPI_COMMUNICATIONS_COMMBUFFEROPS_INLINE_HPP_
#define SRC_COMPONENTS_CORE_SRC_MPI_COMMUNICATIONS_COMMBUFFEROPS_INLINE_HPP_

namespace geosx
{

//template< typename T >
//localIndex CommBufferOps::PackSize( T const & var )
//{
//  return sizeof(T);
//}



template< bool DO_PACKING,
          typename T,
          typename INDEX_TYPE >
typename std::enable_if< std::is_arithmetic<T>::value, localIndex >::type
CommBufferOps::Pack( char*&  buffer,
                     T const * const var,
                     INDEX_TYPE const length )
{
  localIndex sizeOfPackedChars = 0;

  sizeOfPackedChars += Pack<DO_PACKING>( buffer, length );

  sizeOfPackedChars += length*sizeof(T);
  static_if( DO_PACKING )
  {
    memcpy( buffer, var, length*sizeof(T) );
    buffer += length*sizeof(T);
  });

  return sizeOfPackedChars;
}


template< bool DO_PACKING,
          typename T,
          typename INDEX_TYPE >
typename std::enable_if< !std::is_arithmetic<T>::value, localIndex >::type
CommBufferOps::Pack( char*&  buffer,
                     T const * const var,
                     INDEX_TYPE const length )
{
  localIndex sizeOfPackedChars = 0;

  sizeOfPackedChars += Pack<DO_PACKING>( buffer, length );

  for( INDEX_TYPE a=0 ; a<length ; ++a )
  {
    sizeOfPackedChars += Pack<DO_PACKING>( buffer, var[ a ] );
  }

  return sizeOfPackedChars;
}

template< typename T,
          typename INDEX_TYPE >
localIndex
CommBufferOps::Unpack( char const *& buffer,
                       T * const var,
                       INDEX_TYPE & length )
{
  localIndex sizeOfUnpackedChars = 0;

  sizeOfUnpackedChars += Unpack( buffer, length );

  char * const ptr_var = reinterpret_cast<char *>(var);
  memcpy( ptr_var, buffer, length * sizeof(T) );
  buffer += length * sizeof(T);

  return sizeOfUnpackedChars;
}


template< bool DO_PACKING,
          typename T,
          typename INDEX_TYPE >
localIndex
CommBufferOps::Pack( char*&  buffer,
                     T const * const var,
                     INDEX_TYPE const * const indices,
                     INDEX_TYPE const length )
{
  localIndex sizeOfPackedChars = 0;

  sizeOfPackedChars += Pack<DO_PACKING>( buffer, length );

  for( INDEX_TYPE a=0 ; a<length ; ++a )
  {
    sizeOfPackedChars += Pack<DO_PACKING>( buffer, var[ indices[a] ] );
  }

  return sizeOfPackedChars;
}

template< typename T,
          typename INDEX_TYPE >
localIndex
CommBufferOps::Unpack( char const *& buffer,
                       T * const var,
                       INDEX_TYPE const * const indices,
                       INDEX_TYPE & length )
{
  localIndex sizeOfUnpackedChars = 0;

  sizeOfUnpackedChars += Unpack( buffer, length );

  for( INDEX_TYPE a=0 ; a<length ; ++a )
  {
    sizeOfUnpackedChars += Unpack( buffer, var[ indices[a] ] );
  }

  return sizeOfUnpackedChars;
}



//******************************************************************************
/**
 * @author settgast
 * @tparam type of data to pack
 * @param var data to pack
 * @return size (in bytes) of packed data
 */
template< bool DO_PACKING,
          typename T >
localIndex CommBufferOps::Pack( array<char> & buffer, T const & var )
{
  localIndex const sizeOfPackedChars = sizeof(T);
  static_if( DO_PACKING )
  {
    buffer.resize(buffer.size() + sizeOfPackedChars);
    char* p_buffer = &(buffer.back()) - sizeOfPackedChars + 1;
    memcpy( p_buffer, &var, sizeOfPackedChars );
  });
  return sizeOfPackedChars;
}

/**
 *
 * @param buffer
 * @param var
 * @return
 */
template< bool DO_PACKING,
          typename T >
localIndex CommBufferOps::Pack( char*&  buffer, T const & var )
{
  localIndex const sizeOfPackedChars = sizeof(T);
  static_if( DO_PACKING )
  {
    memcpy( buffer, &var, sizeOfPackedChars );
    buffer += sizeOfPackedChars;
  });
  return sizeOfPackedChars;
}

/**
 * @author settgast
 * @param buffer
 * @param var
 * @return size (in bytes) of unpacked data
 */
template< typename T>
localIndex CommBufferOps::Unpack( char const *& buffer, T & var )
{
  localIndex const sizeOfUnpackedChars = sizeof(T);
  const T* const tbuffer = reinterpret_cast<const T*>(buffer);
  var = *tbuffer;
  buffer += sizeOfUnpackedChars;
  return sizeOfUnpackedChars;
}





template< bool DO_PACKING >
localIndex CommBufferOps::Pack( array<char> & buffer, string const & var )
{
  string::size_type sizeOfPackedChars = 0;

  Pack<DO_PACKING>( buffer, sizeOfPackedChars );

  const char* cvar = var.data();
  for( string::size_type i=0 ; i<var.size() ; ++i )
  {
    buffer.push_back(*(cvar++));
  }
  sizeOfPackedChars += var.size();

  return integer_conversion<localIndex>(sizeOfPackedChars);
}

template< bool DO_PACKING >
localIndex CommBufferOps::Pack( char*& buffer,  const std::string& var )
{
  string::size_type sizeOfPackedChars = var.size();

  Pack<DO_PACKING>( buffer, sizeOfPackedChars );

  static_if( DO_PACKING )
  {
    for( string::size_type i=0 ; i<sizeOfPackedChars ; ++i )
    {
      *buffer = var[i];
      buffer++;
    }
  });

  sizeOfPackedChars += sizeof( localIndex );
  return integer_conversion<localIndex>(sizeOfPackedChars);
}

template< bool DO_PACKING,
          typename T >
localIndex CommBufferOps::Pack( char *& buffer, set<T> const & var )
{

  localIndex sizeOfPackedChars = 0;

//  const localIndex length = integer_conversion<localIndex>(var.size());
//
//  sizeOfPackedChars += Pack<DO_PACKING>( buffer, length );
//
//  for( typename std::set<T>::const_iterator i=var.begin() ; i!=var.end() ; ++i )
//  {
//    sizeOfPackedChars += Pack<DO_PACKING>( buffer, *i);
//  }

  return sizeOfPackedChars;
}


template< typename T>
localIndex CommBufferOps::Unpack( char const *& buffer, set<T> & setToRead )
{
  setToRead.clear();

  localIndex sizeOfUnpackedChars = 0;
//
//  localIndex set_length;
//  sizeOfUnpackedChars += Unpack( buffer, set_length );
//
//
//  for( localIndex a=0 ; a<set_length ; ++a )
//  {
//    T temp;
//    sizeOfUnpackedChars += Unpack( buffer, temp );
//    setToRead.insert( temp );
//  }

  return sizeOfUnpackedChars;
}



//********************************************************************************************************************
template< bool DO_PACKING , typename T_KEY, typename T_VAL >
localIndex CommBufferOps::Pack( char*& buffer,
                                std::map<T_KEY,T_VAL> const & var )
{

  localIndex sizeOfPackedChars = 0;

  const typename std::map<T_KEY,T_VAL>::size_type length = var.size();

  sizeOfPackedChars += Pack<DO_PACKING>( buffer, length );

  for( typename std::map<T_KEY,T_VAL>::const_iterator i=var.begin() ; i!=var.end() ; ++i )
  {
    sizeOfPackedChars += Pack<DO_PACKING>( buffer, i->first );
    sizeOfPackedChars += Pack<DO_PACKING>( buffer, i->second );
  }

  return sizeOfPackedChars;
}

template< typename T_KEY, typename T_VAL >
localIndex CommBufferOps::Unpack( char const *& buffer, std::map<T_KEY,T_VAL>& map )
{
  map.clear();

  localIndex sizeOfUnpackedChars = 0;

  typename std::map<T_KEY,T_VAL>::size_type map_length;
  sizeOfUnpackedChars += Unpack( buffer, map_length );


  for( typename std::map<T_KEY,T_VAL>::size_type a=0 ; a<map_length ; ++a )
  {
    T_KEY key;
    T_VAL value;
    sizeOfUnpackedChars += Unpack( buffer, key );
    sizeOfUnpackedChars += Unpack( buffer, value );

    map[key] = value;
  }

  return sizeOfUnpackedChars;
}

template< bool DO_PACKING ,
          typename T_KEY,
          typename T_VAL,
          typename T_INDICES >
localIndex CommBufferOps::Pack( char*& buffer,
                                std::map<T_KEY,T_VAL> const & var,
                                T_INDICES const & packIndices )
{

  localIndex sizeOfPackedChars = 0;

  const typename std::map<T_KEY,T_VAL>::size_type length = var.size();

  sizeOfPackedChars += Pack<DO_PACKING>( buffer, length );

  for( typename std::map<T_KEY,T_VAL>::const_iterator i=var.begin() ; i!=var.end() ; ++i )
  {
    sizeOfPackedChars += Pack<DO_PACKING>( buffer, i->first );
    sizeOfPackedChars += Pack<DO_PACKING>( buffer, i->second, packIndices );
  }

  return sizeOfPackedChars;
}

template< typename T_KEY,
          typename T_VAL,
          typename T_INDICES >
localIndex CommBufferOps::Unpack( char const *& buffer,
                                  std::map<T_KEY,T_VAL>& map,
                                  T_INDICES const & unpackIndices )
{
  map.clear();

  localIndex sizeOfUnpackedChars = 0;

  typename std::map<T_KEY,T_VAL>::size_type map_length;
  sizeOfUnpackedChars += Unpack( buffer, map_length );


  for( typename std::map<T_KEY,T_VAL>::size_type a=0 ; a<map_length ; ++a )
  {
    T_KEY key;
    T_VAL value;
    sizeOfUnpackedChars += Unpack( buffer, key );
    sizeOfUnpackedChars += Unpack( buffer, value, unpackIndices );

    map[key] = value;
  }

  return sizeOfUnpackedChars;
}


template< bool DO_PACKING, typename T, int NDIM, typename INDEX_TYPE >
localIndex CommBufferOps::Pack( char*& buffer,
                                multidimensionalArray::ManagedArray<T,NDIM,INDEX_TYPE> const & var )
{
  localIndex sizeOfPackedChars = 0;

  sizeOfPackedChars += Pack<DO_PACKING>( buffer, var.dims(), NDIM );

  sizeOfPackedChars += Pack<DO_PACKING>( buffer, var.strides(), NDIM );

  const localIndex length = var.size();
  localIndex sizeOfPackedArrayChars = length*sizeof(T);

  sizeOfPackedChars += Pack<DO_PACKING>( buffer, var.data(), length );

  return sizeOfPackedChars;
}


template< typename T, int NDIM, typename INDEX_TYPE >
localIndex CommBufferOps::Unpack( char const *& buffer,
                                  multidimensionalArray::ManagedArray<T,NDIM,INDEX_TYPE> & var )
{
  localIndex sizeOfUnpackedChars = 0;

  int numDimsRead = 0;
  INDEX_TYPE dims[NDIM];
  sizeOfUnpackedChars += Unpack( buffer, dims, numDimsRead );

  if( numDimsRead != NDIM )
  {
    GEOS_ERROR( "error reading dims");
  }
  else
  {
    var.resize( NDIM, dims );
  }

  int numStrideRead = 0;
  INDEX_TYPE strides[NDIM];
  sizeOfUnpackedChars += Unpack( buffer, strides, numStrideRead );

  if( numStrideRead != NDIM )
  {
    GEOS_ERROR( "error reading strides");
  }
  else
  {
    var.resize( NDIM, strides );
  }

  localIndex numValuesRead;
  sizeOfUnpackedChars += Unpack( buffer, var.data(), numValuesRead );
  if( numValuesRead != var.size() )
  {
    GEOS_ERROR( "error reading data");
  }


  return sizeOfUnpackedChars;

}


template< bool DO_PACKING,
          typename T,
          int NDIM,
          typename T_indices,
          typename INDEX_TYPE >
localIndex CommBufferOps::Pack( char*& buffer,
                                multidimensionalArray::ManagedArray<T,NDIM,INDEX_TYPE> const & var,
                                const T_indices& indices )
{


  localIndex sizeOfPackedChars = 0;

  sizeOfPackedChars += Pack<DO_PACKING>( buffer, var.dims(), NDIM );

  sizeOfPackedChars += Pack<DO_PACKING>( buffer, var.strides(), NDIM );

  const localIndex length = indices.size();
  localIndex sizeOfPackedArrayChars = length*sizeof(T);

  sizeOfPackedChars += Pack<DO_PACKING>( buffer, var.data(), indices.data(), length );

  return sizeOfPackedChars;
}



template< typename T,
          int NDIM,
          typename T_indices,
          typename INDEX_TYPE >
localIndex CommBufferOps::Unpack( char const *& buffer,
                                  multidimensionalArray::ManagedArray<T,NDIM,INDEX_TYPE> & var,
                                  const T_indices& indices )
{
  localIndex sizeOfUnpackedChars = 0;

  int numDimsRead = 0;
  INDEX_TYPE dims[NDIM];
  sizeOfUnpackedChars += Unpack( buffer, dims, numDimsRead );

  if( numDimsRead != NDIM )
  {
    GEOS_ERROR( "error reading dims");
  }
  else
  {
    var.resize( NDIM, dims );
  }

  int numStrideRead = 0;
  INDEX_TYPE strides[NDIM];
  sizeOfUnpackedChars += Unpack( buffer, strides, numStrideRead );

  if( numStrideRead != NDIM )
  {
    GEOS_ERROR( "error reading strides");
  }
  else
  {
    var.resize( NDIM, strides );
  }

  localIndex numValuesRead;
  sizeOfUnpackedChars += Unpack( buffer, var.data(), indices.data(), numValuesRead );
  if( numValuesRead != var.size() )
  {
    GEOS_ERROR( "error reading data");
  }


  return sizeOfUnpackedChars;

}




}



#endif /* SRC_COMPONENTS_CORE_SRC_MPI_COMMUNICATIONS_COMMBUFFEROPS_INLINE_HPP_ */
