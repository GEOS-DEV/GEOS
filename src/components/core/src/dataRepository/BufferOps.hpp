// Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at
// the Lawrence Livermore National Laboratory. LLNL-CODE-746361. All Rights
// reserved. See file COPYRIGHT for details.
//
// This file is part of the GEOSX Simulation Framework.

//
// GEOSX is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License (as published by the Free
// Software Foundation) version 2.1 dated February 1999.


#ifndef DATAREPOSITORY_BUFFEROPS_H_
#define DATAREPOSITORY_BUFFEROPS_H_

#include "common/DataTypes.hpp"
//#include "codingUtilities/Utilities.hpp"
#include "common/Logger.hpp"
#include "codingUtilities/static_if.hpp"
#include "codingUtilities/GeosxTraits.hpp"
#include <type_traits>
namespace geosx
{



namespace bufferOps
{


  /** @name Packing/Unpacking functions for pointer arrays of arbitrary types
   */
  ///@{



  template< bool DO_PACKING,
            typename T,
            typename INDEX_TYPE >
  typename std::enable_if< std::is_trivial<T>::value, localIndex >::type
  Pack( char*&  buffer,
        T const * const var,
        INDEX_TYPE const length );


  template< bool DO_PACKING,
            typename T,
            typename INDEX_TYPE >
  typename std::enable_if< !std::is_trivial<T>::value, localIndex >::type
  Pack( char*&  buffer,
        T const * const var,
        INDEX_TYPE const length );


  template< typename T, typename INDEX_TYPE >
  typename std::enable_if< std::is_trivial<T>::value, localIndex >::type
  Unpack( char const *& buffer, T * const var, INDEX_TYPE const length );

  template< typename T, typename INDEX_TYPE >
  typename std::enable_if< !std::is_trivial<T>::value, localIndex >::type
  Unpack( char const *& buffer, T * const var, INDEX_TYPE const length );


  template< bool DO_PACKING, typename T, typename INDEX_TYPE >
  localIndex Pack( char*&  buffer,
                   T const * const var,
                   INDEX_TYPE const * const indices,
                   INDEX_TYPE const length  );

  template< typename T, typename INDEX_TYPE >
  localIndex Unpack( char const *& buffer,
                     T * const var,
                     INDEX_TYPE const * const indices,
                     INDEX_TYPE & length );


  ///@}








  /** @name Packing/Unpacking functions for arbitrary types
   */
  ///@{


  template< bool DO_PACKING, typename T >
  typename std::enable_if< std::is_trivial<T>::value, localIndex >::type
  Pack( char*&  buffer, T const & var );

  template< typename T >
  typename std::enable_if< std::is_trivial<T>::value, localIndex >::type
  Unpack( char const *& buffer, T & var );




  ///@}




  /** @name Packing/Unpacking functions for strings
   */
  ///@{

  template< bool DO_PACKING >
  localIndex Pack( char*& buffer, string const & var );

  localIndex Unpack( char const *& buffer, string& var );

  ///@}


  /** @name Packing/Unpacking functions for single point tensor objects
   */
  ///@{
  template< bool DO_PACKING, typename T >
  typename std::enable_if< traits::is_tensorT<T>::value, localIndex >::type
  Pack( char*& buffer, T const & var );

  template< typename T >
  typename std::enable_if< traits::is_tensorT<T>::value, localIndex >::type
  Unpack( char const *& buffer, T & var );
  ///@}


  /** @name Packing functions for entire collections of data
   */
  ///@{


  template< bool DO_PACKING,
            typename T,
            int NDIM,
            typename INDEX_TYPE=std::int_fast32_t >
  typename std::enable_if<
  is_packable_array< multidimensionalArray::ManagedArray<T,NDIM,INDEX_TYPE> >::value,
  localIndex >::type
  Pack( char*& buffer,
        multidimensionalArray::ManagedArray<T,NDIM,INDEX_TYPE> const & var );

  template< typename T,
            int NDIM,
            typename INDEX_TYPE=std::int_fast32_t >
  typename std::enable_if<
  is_packable_array< multidimensionalArray::ManagedArray<T,NDIM,INDEX_TYPE> >::value,
  localIndex >::type
  Unpack( char const *& buffer,
          multidimensionalArray::ManagedArray<T,NDIM,INDEX_TYPE> & var );



  template< bool DO_PACKING, typename T   >
  localIndex Pack( char*& buffer, set<T> const & var );


  template< typename T>
  localIndex Unpack( char const *& buffer, set<T> & setToRead );




  template< bool DO_PACKING, typename T_KEY, typename T_VAL >
  typename std::enable_if<
  is_packable_map< map<T_KEY,T_VAL> >::value,
  localIndex >::type
  Pack( char*& buffer,
        std::map<T_KEY,T_VAL> const & var );

  template< typename T_KEY, typename T_VAL >
  typename std::enable_if<
  is_packable_map< map<T_KEY,T_VAL> >::value,
  localIndex >::type
  Unpack( char const *& buffer,
          std::map<T_KEY,T_VAL>& map );


  ///@}

  /** @name Packing functions for filtered/masked collections of data
   */
  ///@{




  template< bool DO_PACKING,
            typename T,
            int NDIM,
            typename T_indices,
            typename INDEX_TYPE=std::int_fast32_t >
  typename std::enable_if<
  is_packable_array< multidimensionalArray::ManagedArray<T,NDIM,INDEX_TYPE> >::value,
  localIndex >::type
  Pack( char*& buffer,
        multidimensionalArray::ManagedArray<T,NDIM,INDEX_TYPE> const & var,
        const T_indices& indices );

  template< typename T,
            int NDIM,
            typename T_indices,
            typename INDEX_TYPE=std::int_fast32_t >
  localIndex Unpack( char const *& buffer,
                            multidimensionalArray::ManagedArray<T,NDIM,INDEX_TYPE> & var,
                            const T_indices & indices );


  template< bool DO_PACKING, typename T_KEY, typename T_VAL, typename T_INDICES >
  typename std::enable_if<
  is_packable_map< map<T_KEY,T_VAL> >::value && is_packable_by_index<T_VAL>::value,
  localIndex >::type
  Pack( char*& buffer,
        const std::map<T_KEY,T_VAL>& var,
        T_INDICES const & packIndices );

  template< typename T_KEY, typename T_VAL, typename T_INDICES >
  typename std::enable_if<
  is_packable_map< map<T_KEY,T_VAL> >::value && is_packable_by_index<T_VAL>::value,
  localIndex >::type
  Unpack( char const *& buffer,
          std::map<T_KEY,T_VAL>& map,
          T_INDICES const & unpackIndices );

  ///@}







//  template< bool DO_PACKING >
//  localIndex Pack( char*& buffer,
//            localIndex const * const var,
//            localIndex const length,
//            globalIndex_array const & localToGlobalMap );

//  localIndex Unpack( char const *& buffer,
//              localIndex_array & var,
//              map<globalIndex,localIndex> const & globalToLocalMap );


//  template< bool DO_PACKING >
//  localIndex Pack( char*& buffer,
//            multidimensionalArray::ManagedArray<localIndex,2,localIndex> const & var,
//            localIndex_array const & indices,
//            globalIndex_array const & localToGlobalMap );
//
//  localIndex Unpack( char const *& buffer,
//              multidimensionalArray::ManagedArray<localIndex,2,localIndex> & var,
//              localIndex_array const & indices,
//              globalIndex_array const & globalToLocalMap );

//
//  template< bool DO_PACKING >
//  int Pack( char*& buffer,
//            array< localIndex_array > const & var,
//            localIndex_array const & indices,
//            globalIndex_array const & localToGlobalMap );
//
//  int Unpack( char const *& buffer,
//              array< localIndex_array > & var,
//              localIndex_array const & indices,
//              map<globalIndex,localIndex> const & globalToLocalMap,
//              map<globalIndex,localIndex> const & relatedObjectGlobalToLocalMap );
//

















  template< typename... VARPACK >
  localIndex PackSize( VARPACK const &&...  pack )
  {
    char* junk = nullptr;
    return Pack<false>( junk, pack... );
  }

  template< typename... VARPACK >
  localIndex PackSize( VARPACK &&...  pack )
  {
    char* junk = nullptr;
    return Pack<false>( junk, pack... );
  }








  template< bool DO_PACKING, typename T >
  typename std::enable_if< !is_packable<T>::value, localIndex >::type
  Pack( char*&  buffer, T const & var )
  {
    GEOS_ERROR("Trying to unpack data type by index list, but type is not packable by index ");
    return 0;
  }


  template< typename T >
  typename std::enable_if< !is_packable<T>::value, localIndex >::type
  Unpack( char const *& buffer, T & var )
  {
    GEOS_ERROR("Trying to unpack data type ("<<typeid(T).name()<<"), but type does not have packable functionality ");
    return 0;
  }

  template< bool DO_PACKING, typename T, typename T_INDICES >
  typename std::enable_if< !is_packable_by_index<T>::value, localIndex >::type
  Pack( char*&  buffer, T const & var, T_INDICES const& )
  {
    GEOS_ERROR("Trying to pack data type by index list, but type is not packable by index ");
    return 0;
  }

  template< typename T, typename T_INDICES >
  typename std::enable_if< !is_packable_by_index<T>::value, localIndex >::type
  Unpack( char const *& buffer, T & var, T_INDICES const& )
  {
    GEOS_ERROR("Trying to unpack data type by index list, but type is not packable by index ");
    return 0;
  }


}



}

#include "BufferOps_inline.hpp"


#endif /* BUFVECTOR_H_ */
