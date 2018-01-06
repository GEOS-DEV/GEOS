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



//******************************************************************************
/**
 * @author settgast
 * @tparam type of data to pack
 * @param var data to pack
 * @return size (in bytes) of packed data
 */
template< typename T >
localIndex CommBufferOps::Pack( array<char> & buffer, T const & var )
{
  localIndex const sizeOfPackedChars = sizeof(T);
  buffer.resize(buffer.size() + sizeOfPackedChars);
  char* p_buffer = &(buffer.back()) - sizeOfPackedChars + 1;
  memcpy( p_buffer, &var, sizeOfPackedChars );
  return sizeOfPackedChars;
}

/**
 *
 * @param buffer
 * @param var
 * @return
 */
template< typename T>
localIndex CommBufferOps::Pack( char*&  buffer, T const & var )
{
  localIndex const sizeOfPackedChars = sizeof(T);
  memcpy( buffer, &var, sizeOfPackedChars );
  buffer += sizeOfPackedChars;
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

//******************************************************************************
/**
 * @author settgast
 * @param container
 * @return
 */
template< typename T >
localIndex CommBufferOps::Pack( array<char> & buffer, array<T> const & var )
{
  localIndex sizeOfPackedChars = 0;

  const localIndex length = var.size();
  localIndex sizeOfPackedArrayChars = length*sizeof(T);


  sizeOfPackedChars += Pack( buffer, length );

  buffer.resize( buffer.size() + sizeOfPackedArrayChars );

  char* p_buffer = &(buffer.back()) - sizeOfPackedArrayChars + 1;
  memcpy( p_buffer, var.data(), sizeOfPackedArrayChars );
  sizeOfPackedChars += sizeOfPackedArrayChars;

  return sizeOfPackedChars;
}

template< typename T >
localIndex CommBufferOps::Pack( char *& buffer, array<T> const & var )
{
  localIndex sizeOfPackedChars = 0;

  const localIndex length = var.size();
  localIndex sizeOfPackedArrayChars = length*sizeof(T);

  sizeOfPackedChars += Pack( buffer, length );

  memcpy( buffer, var.data(), sizeOfPackedArrayChars );
  sizeOfPackedChars += sizeOfPackedArrayChars;
  buffer += sizeOfPackedArrayChars;

  return sizeOfPackedChars;
}




template< typename T>
localIndex CommBufferOps::Unpack( char const *& buffer, array<T>& array )
{

  localIndex sizeOfUnpackedChars = 0;

  localIndex array_length;
  sizeOfUnpackedChars += Unpack( buffer, array_length );
  array.resize(array_length);
  localIndex length = array_length * sizeof(T);

  memcpy( array.data(), buffer, length );
  buffer += length;
  sizeOfUnpackedChars += length;

  return sizeOfUnpackedChars;
}




template< typename T >
localIndex CommBufferOps::Pack( array<char> & buffer, const std::set<T>& var )
{

  localIndex sizeOfPackedChars = 0;

  const localIndex length = integer_conversion<localIndex>(var.size());

  sizeOfPackedChars += Pack( buffer, length );

  for( typename std::set<T>::const_iterator i=var.begin() ; i!=var.end() ; ++i )
  {
    sizeOfPackedChars += Pack( buffer, *i);
  }

  return sizeOfPackedChars;
}

template< typename T >
localIndex CommBufferOps::Pack( char *& buffer, std::set<T> const & var )
{

  localIndex sizeOfPackedChars = 0;

  const localIndex length = integer_conversion<localIndex>(var.size());

  sizeOfPackedChars += Pack( buffer, length );

  for( typename std::set<T>::const_iterator i=var.begin() ; i!=var.end() ; ++i )
  {
    sizeOfPackedChars += Pack( buffer, *i);
  }

  return sizeOfPackedChars;
}


template< typename T>
localIndex CommBufferOps::Unpack( char const *& buffer, std::set<T> & setToRead )
{
  setToRead.clear();

  localIndex sizeOfUnpackedChars = 0;

  localIndex set_length;
  sizeOfUnpackedChars += Unpack( buffer, set_length );


  for( localIndex a=0 ; a<set_length ; ++a )
  {
    T temp;
    sizeOfUnpackedChars += Unpack( buffer, temp );
    setToRead.insert( temp );
  }

  return sizeOfUnpackedChars;
}



  //********************************************************************************************************************
  template< typename T_KEY, typename T_VAL >
  localIndex Pack( array<char> & buffer,
                   std::map<T_KEY,T_VAL> const & var )
  {

    localIndex sizeOfPackedChars = 0;

    const typename std::map<T_KEY,T_VAL>::size_type length = var.size();

    sizeOfPackedChars += Pack( buffer, length );



    for( typename std::map<T_KEY,T_VAL>::const_iterator i=var.begin() ; i!=var.end() ; ++i )
    {
      sizeOfPackedChars += Pack( buffer, i->first );
      sizeOfPackedChars += Pack( buffer, i->second );
    }

    return sizeOfPackedChars;
  }



  template< typename T_KEY, typename T_VAL >
  static localIndex Unpack( char const *& buffer, std::map<T_KEY,T_VAL>& map )
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

//********************************************************************************************************************
template< typename T, typename T_indices >
localIndex CommBufferOps::Pack( array<char> & buffer,
                                   const array<T>& container,
                                   const T_indices& indices )
{
  localIndex sizeOfPackedChars = 0;

  const typename T_indices::size_type length = indices.size();
  localIndex sizeOfPackedArrayChars = length*sizeof(T);
  buffer.resize( buffer.size() + sizeOfPackedArrayChars );

  sizeOfPackedChars += Pack( buffer, length );

  char* p_buffer = &(buffer.back()) - sizeOfPackedArrayChars + 1;
  for( typename T_indices::const_iterator i=indices.begin() ; i!=indices.end() ; ++i )
  {
    Pack( p_buffer, container[*i] );
  }
  sizeOfPackedChars += sizeOfPackedArrayChars;

  return sizeOfPackedChars;
}


template< typename T, typename T_indices >
localIndex CommBufferOps::Pack( char *& buffer,
                                   const array<T>& container,
                                   const T_indices& indices )
{
  localIndex sizeOfPackedChars = 0;

  const typename T_indices::size_type length = indices.size();
  localIndex sizeOfPackedArrayChars = length*sizeof(T);

  sizeOfPackedChars += Pack( buffer, length );

  for( typename T_indices::const_iterator i=indices.begin() ; i!=indices.end() ; ++i )
  {
    Pack( buffer, container[*i] );
  }
  sizeOfPackedChars += sizeOfPackedArrayChars;

  return sizeOfPackedChars;
}

template< typename T, typename T_indices >
localIndex CommBufferOps::Unpack( char const *& buffer,
                                  array<T>& array,
                                  const T_indices& indices )
{

  localIndex sizeOfUnpackedChars = 0;

  typename T_indices::size_type array_length;

  sizeOfUnpackedChars += Unpack( buffer, array_length );

  if( array_length != indices.size() )
  {
    GEOS_ERROR("CommBuffer::PrivateUnpackArray(): incorrect number of data");
  }


  for( typename T_indices::const_iterator i=indices.begin() ; i!=indices.end() ; ++i )
  {
    const T* const tbuffer = static_cast<const T*>(buffer);
    array[*i] = *tbuffer;
    buffer += sizeof(T);
  }
  sizeOfUnpackedChars += array_length * sizeof(T);

  return sizeOfUnpackedChars;
}






//**********************************************************************************************************************

//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
////********************************************************************************************************************
//template< typename T>
//localIndex CommBuffer::PrivatePackGlobal( const T& container, globalIndex_array const& localToGlobal )
//{
//  const localIndex length = integer_conversion<localIndex>(container.size());
//  localIndex sizeOfPackedChars = 0;
//
////  std::cout<<"container.size() = "<<length<<std::endl;
//
//  sizeOfPackedChars += this->Pack(length);
//  for( typename T::const_iterator i=container.begin() ; i!=container.end() ; ++i )
//  {
//    const globalIndex temp = localToGlobal[*i];
//    sizeOfPackedChars += this->Pack( temp );
//  }
//  return sizeOfPackedChars;
//}
//
//
//inline localIndex CommBuffer::UnpackGlobal( char const *& buffer, const std::map<globalIndex,localIndex>& globalToLocal, localIndex_array& array )
//{
//  array.clear();
//  localIndex sizeOfUnpackedChars = 0;
//
//
//  localIndex array_length;
//  sizeOfUnpackedChars += Unpack( buffer, array_length );
//
//  array.resize(array_length);
//  globalIndex_array temp(array_length);
//  localIndex length = array_length * sizeof(localIndex);
//
//  memcpy( temp.data(), buffer, length );
//  buffer += length;
//  sizeOfUnpackedChars += length;
//
//
//  for( localIndex i=0 ; i<array_length ; ++i )
//  {
//    const localIndex li = stlMapLookup( globalToLocal, temp[i] );
//    array[i] = li;
//
//  }
//
//  return sizeOfUnpackedChars;
//}
//
//inline localIndex CommBuffer::UnpackGlobal( char const *& buffer, const std::map<globalIndex,localIndex>& globalToLocal, lSet& array )
//{
//  array.clear();
//  localIndex sizeOfUnpackedChars = 0;
//
//
//  localIndex array_length;
//  sizeOfUnpackedChars += Unpack( buffer, array_length );
//
////  std::cout<<"array_length = "<<array_length<<std::endl;
//
//  globalIndex_array temp(array_length);
//  localIndex length = array_length * sizeof(localIndex);
//
//  memcpy( temp.data(), buffer, length );
//  buffer += length;
//  sizeOfUnpackedChars += length;
//
//
//  for( auto i=0 ; i<array_length ; ++i )
//  {
//    const localIndex li = stlMapLookup( globalToLocal, temp[i] );
//    array.insert( li );
//  }
//
//  return sizeOfUnpackedChars;
//}
//
//
//
////**********************************************************************************************************************
//
//template< typename T, typename T_indices >
//inline localIndex CommBuffer::PrivatePackRelationT<T,T_indices>::Pack( CommBuffer& buffer, const T& relation, const T_indices& indices, const bool packGlobal )
//{
//  localIndex sizeOfPackedChars = 0;
//
//  if( packGlobal )
//  {
//    globalIndex_array const& localToGlobal = relation.RelatedObjectLocalToGlobal();
//    for( typename T_indices::const_iterator i = indices.begin() ; i != indices.end() ; ++i )
//    {
//      sizeOfPackedChars += buffer.PackGlobal( relation[*i], localToGlobal );
//    }
//  }
//  else
//  {
//    for( typename T_indices::const_iterator i = indices.begin() ; i != indices.end() ; ++i )
//    {
//      sizeOfPackedChars += buffer.Pack( relation[*i] );
//    }
//  }
//
//  return sizeOfPackedChars;
//}
//template<typename T_indices>
//inline localIndex CommBuffer::PrivatePackRelationT<OneToOneRelation,T_indices>::Pack( CommBuffer& buffer, const OneToOneRelation& relation,
//                                                                                     const T_indices& indices, const bool packGlobal )
//{
//  localIndex sizeOfPackedChars = 0;
//
//  if( packGlobal )
//  {
//    globalIndex_array const& localToGlobal = relation.RelatedObjectLocalToGlobal();
//    for( typename T_indices::const_iterator i = indices.begin() ; i != indices.end() ; ++i )
//    {
//      sizeOfPackedChars += buffer.Pack( localToGlobal[relation[*i]] );
//    }
//  }
//  else
//  {
//    for( typename T_indices::const_iterator i = indices.begin() ; i != indices.end() ; ++i )
//    {
//      sizeOfPackedChars += buffer.Pack( relation[*i] );
//    }
//  }
//  return sizeOfPackedChars;
//}
//template<typename T_indices>
//inline localIndex CommBuffer::PrivatePackRelationT<FixedOneToManyRelation,T_indices>::Pack( CommBuffer& buffer, const FixedOneToManyRelation& relation,
//                                                                                           const T_indices& indices, const bool packGlobal )
//{
//  localIndex sizeOfPackedChars = 0;
//
//  sizeOfPackedChars += buffer.Pack( static_cast<localIndex>(relation.size(1)) );
//
//  if( packGlobal )
//  {
//    globalIndex_array const& localToGlobal = relation.RelatedObjectLocalToGlobal();
//    for( typename T_indices::const_iterator i = indices.begin() ; i != indices.end() ; ++i )
//    {
//      for( auto j=0u ; j<relation.size(1) ; ++j )
//      {
//        sizeOfPackedChars += buffer.Pack( localToGlobal[relation[*i][j]] );
//      }
//    }
//  }
//  else
//  {
//    for( typename T_indices::const_iterator i = indices.begin() ; i != indices.end() ; ++i )
//    {
//      for( auto j=0u ; j<relation.size(1) ; ++j )
//      {
//        sizeOfPackedChars += buffer.Pack( relation[*i][j] );
//      }
//    }
//  }
//  return sizeOfPackedChars;
//}
//
//
//template< typename T >
//inline localIndex CommBuffer::PrivateUnpackRelation( char const *& buffer, T& relation, localIndex_array const& indices, const bool unpackGlobal )
//{
//  localIndex sizeOfUnpackedChars = 0;
//
//  if( unpackGlobal )
//  {
//    const std::map<globalIndex,localIndex>& globalToLocal = relation.RelatedObjectGlobalToLocal();
//    for( localIndex_array::const_iterator i = indices.begin() ; i != indices.end() ; ++i )
//    {
//      sizeOfUnpackedChars += CommBuffer::UnpackGlobal( buffer, globalToLocal, relation[*i] );
//    }
//  }
//  else
//  {
//    for( localIndex_array::const_iterator i = indices.begin() ; i != indices.end() ; ++i )
//    {
//      sizeOfUnpackedChars += CommBuffer::Unpack( buffer, relation[*i] );
//    }
//  }
//
//  return sizeOfUnpackedChars;
//}
//
//template<>
//inline localIndex CommBuffer::PrivateUnpackRelation( char const *& buffer, OneToOneRelation& relation, localIndex_array const& indices, const bool unpackGlobal )
//{
//  localIndex sizeOfUnpackedChars = 0;
//
//  if( unpackGlobal )
//  {
//    const std::map<globalIndex,localIndex>& globalToLocal = relation.RelatedObjectGlobalToLocal();
//    for( localIndex_array::const_iterator i = indices.begin() ; i != indices.end() ; ++i )
//    {
//      globalIndex temp;
//      sizeOfUnpackedChars += CommBuffer::Unpack( buffer, temp );
//      relation[*i] = stlMapLookup(globalToLocal,temp);
//    }
//  }
//  else
//  {
//    for( localIndex_array::const_iterator i = indices.begin() ; i != indices.end() ; ++i )
//    {
//      sizeOfUnpackedChars += CommBuffer::Unpack( buffer, relation[*i] );
//    }
//  }
//
//  return sizeOfUnpackedChars;
//}
//
//template<>
//inline localIndex CommBuffer::PrivateUnpackRelation( char const *& buffer, FixedOneToManyRelation& relation, localIndex_array const& indices,
//                                                    const bool unpackGlobal )
//{
//  localIndex sizeOfUnpackedChars = 0;
//
//  localIndex dimension;
//  sizeOfUnpackedChars += CommBuffer::Unpack( buffer, dimension );
//
//  if( dimension != relation.size(1))
//  {
//#ifdef USE_ATK
//    SLIC_ERROR("CommBuffer::PrivateUnpackRelation(): mismatched dimension");
//#endif
//  }
//
//  if( unpackGlobal )
//  {
//    const std::map<globalIndex,localIndex>& globalToLocal = relation.RelatedObjectGlobalToLocal();
//    for( localIndex_array::const_iterator i = indices.begin() ; i != indices.end() ; ++i )
//    {
//      for( localIndex j=0 ; j<dimension ; ++j )
//      {
//        globalIndex temp;
//        sizeOfUnpackedChars += CommBuffer::Unpack( buffer, temp );
//        relation[*i][j] = stlMapLookup(globalToLocal,temp);
//      }
//    }
//  }
//  else
//  {
//    for( localIndex_array::const_iterator i = indices.begin() ; i != indices.end() ; ++i )
//    {
//      for( localIndex j=0 ; j<dimension ; ++j )
//      {
//        sizeOfUnpackedChars += CommBuffer::Unpack( buffer, relation[*i][j] );
//      }
//    }
//  }
//
//  return sizeOfUnpackedChars;
//}
//
//
//inline localIndex CommBuffer::Unpack( char const *& buffer, OneToOneRelation& relation, localIndex_array const& indices, const bool unpackGlobal ) {
//  return PrivateUnpackRelation( buffer, relation, indices, unpackGlobal );
//}
//inline localIndex CommBuffer::Unpack( char const *& buffer, FixedOneToManyRelation& relation, localIndex_array const& indices, const bool unpackGlobal ) {
//  return PrivateUnpackRelation( buffer, relation, indices, unpackGlobal );
//}
//inline localIndex CommBuffer::Unpack( char const *& buffer, OrderedVariableOneToManyRelation& relation, localIndex_array const& indices,
//                                     const bool unpackGlobal ) {
//  return PrivateUnpackRelation( buffer, relation, indices, unpackGlobal );
//}
//inline localIndex CommBuffer::Unpack( char const *& buffer, UnorderedVariableOneToManyRelation& relation, localIndex_array const& indices,
//                                     const bool unpackGlobal ) { return PrivateUnpackRelation( buffer, relation, indices, unpackGlobal ); }
////inline localIndex CommBuffer::Unpack( char const *& buffer,
//// OrderedVariableOneToManyPairRelation& relation, localIndex_array const&
//// indices, const bool unpackGlobal ) { return PrivateUnpackRelation( buffer,
//// relation, indices, unpackGlobal ); }
////inline localIndex CommBuffer::Unpack( char const *& buffer,
//// UnorderedVariableOneToManyPairRelation& relation, localIndex_array const&
//// indices, const bool unpackGlobal ) { return PrivateUnpackRelation( buffer,
//// relation, indices, unpackGlobal ); }
//
//



}



#endif /* SRC_COMPONENTS_CORE_SRC_MPI_COMMUNICATIONS_COMMBUFFEROPS_INLINE_HPP_ */
