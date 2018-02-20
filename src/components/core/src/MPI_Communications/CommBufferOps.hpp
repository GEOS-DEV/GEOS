

#ifndef COMMBUFFEROPS_H_
#define COMMBUFFEROPS_H_

#include "common/DataTypes.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/Logger.hpp"
#include "codingUtilities/static_if.hpp"

namespace geosx
{


class CommBufferOps
{
public:

  /** @name Packing functions for data objects
   */
  ///@{

  template< bool DO_PACKING,
            typename T,
            typename INDEX_TYPE >
  static typename std::enable_if< std::is_arithmetic<T>::value, localIndex >::type
  Pack( char*&  buffer,
                          T const * const var,
                          INDEX_TYPE const length );


  template< bool DO_PACKING,
            typename T,
            typename INDEX_TYPE >
  static typename std::enable_if< !std::is_arithmetic<T>::value, localIndex >::type
  Pack( char*&  buffer,
                          T const * const var,
                          INDEX_TYPE const length );


  template< typename T, typename INDEX_TYPE >
  static localIndex Unpack( char const *& buffer, T * const var, INDEX_TYPE & length );


  template< bool DO_PACKING, typename T, typename INDEX_TYPE >
  static localIndex Pack( char*&  buffer,
                          T const * const var,
                          INDEX_TYPE const * const indices,
                          INDEX_TYPE const length  );

  template< typename T, typename INDEX_TYPE >
  static localIndex Unpack( char const *& buffer,
                            T * const var,
                            INDEX_TYPE const * const indices,
                            INDEX_TYPE & length );




  template< bool DO_PACKING, typename T >
  static localIndex Pack( array<char> & buffer, T const & var );

  template< bool DO_PACKING, typename T >
  static localIndex Pack( char*&  buffer, T const & var );

  template< typename T >
  static localIndex Unpack( char const *& buffer, T & var );



  template< bool DO_PACKING >
  static localIndex Pack( array<char> & buffer, string const & var );

  template< bool DO_PACKING >
  static localIndex Pack( char*& buffer, string const & var );

  template< bool DO_PACKING, typename T_INDICES >
  static localIndex Pack( char*& buffer, string const & var, T_INDICES const& )
  {
    return Pack<DO_PACKING>(buffer,var);
  }

  static localIndex Unpack( char const *& buffer, string& var );

  ///@}



  /** @name Packing functions for entire collections of data
   */
  ///@{

  template< bool DO_PACKING, typename T   >
  static localIndex Pack( char*& buffer, set<T> const & var );

  template< bool DO_PACKING, typename T, typename T_INDICES >
  static localIndex Pack( char*& buffer, set<T> const & var, T_INDICES const & )
  { return 0; }

  template< typename T>
  static localIndex Unpack( char const *& buffer, std::set<T> & setToRead );


  template< bool DO_PACKING, typename T_KEY, typename T_VAL >
  static localIndex Pack( char*& buffer,
                          std::map<T_KEY,T_VAL> const & var );

  template< typename T_KEY, typename T_VAL >
  static localIndex Unpack( char const *& buffer,
                            std::map<T_KEY,T_VAL>& map );

  template< bool DO_PACKING, typename T_KEY, typename T_VAL, typename T_INDICES >
  static localIndex Pack( char*& buffer,
                          const std::map<T_KEY,T_VAL>& var,
                          T_INDICES const & packIndices );

  template< typename T_KEY, typename T_VAL, typename T_INDICES >
  static localIndex Unpack( char const *& buffer,
                            std::map<T_KEY,T_VAL>& map,
                            T_INDICES const & unpackIndices );


  ///@}

  /** @name Packing functions for filtered/masked collections of data
   */
  ///@{



  template< bool DO_PACKING,
            typename T,
            int NDIM,
            typename INDEX_TYPE=std::int_fast32_t >
  static localIndex Pack( char*& buffer,
                          multidimensionalArray::ManagedArray<T,NDIM,INDEX_TYPE> const & var );

  template< typename T,
            int NDIM,
            typename INDEX_TYPE=std::int_fast32_t >
  static localIndex Unpack( char const *& buffer,
                            multidimensionalArray::ManagedArray<T,NDIM,INDEX_TYPE> & var );

  template< bool DO_PACKING,
            typename T,
            int NDIM,
            typename T_indices,
            typename INDEX_TYPE=std::int_fast32_t >
  static localIndex Pack( char*& buffer,
                          multidimensionalArray::ManagedArray<T,NDIM,INDEX_TYPE> const & var,
                          const T_indices& indices );

  template< typename T,
            int NDIM,
            typename T_indices,
            typename INDEX_TYPE=std::int_fast32_t >
  static localIndex Unpack( char const *& buffer,
                            multidimensionalArray::ManagedArray<T,NDIM,INDEX_TYPE> & var,
                            const T_indices& indices );



  ///@}


  template< typename... VARPACK >
  static localIndex PackSize( VARPACK const &&...  pack )
  {
    char* junk = nullptr;
    return Pack<false>( junk, pack... );
  }

  template< typename... VARPACK >
  static localIndex PackSize( VARPACK &&...  pack )
  {
    char* junk = nullptr;
    return Pack<false>( junk, pack... );
  }


};



}

#include "CommBufferOps_inline.hpp"


#endif /* BUFVECTOR_H_ */
