

#ifndef COMMBUFFEROPS_H_
#define COMMBUFFEROPS_H_

#include "common/DataTypes.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/Logger.hpp"


namespace geosx
{


class CommBufferOps
{
public:

  /** @name Packing functions for data objects
   */
  ///@{

  template< typename T >
  static localIndex Pack( array<char> & buffer, T const & var );

  template< typename T >
  static localIndex Pack( char*&  buffer, T const & var );

  template< typename T >
  static localIndex Unpack( char const *& buffer, T & var );

  static localIndex Pack( array<char> & buffer, string const & var );

  static localIndex Pack( char*& buffer, string const & var );

  static localIndex Unpack( char const *& buffer, string& var );

  ///@}



  /** @name Packing functions for entire collections of data
   */
  ///@{

  template< typename T >
  static localIndex Pack( array<char> & buffer, array<T> const & var );

  template< typename T >
  static localIndex Pack( char*& buffer, array<T> const & var );

  template< typename T >
  static localIndex Unpack( char const *& buffer, array<T>& var );

  static localIndex Pack( array<char> & buffer, array<string> const & var );

  static localIndex Pack( char*& buffer, array<string> const & var );

  static localIndex Unpack( char const *& buffer, array<string> & var );

  template< typename T, int NDIM, typename INDEX_TYPE=std::int_fast32_t >
  static localIndex Pack( char*& buffer, ManagedArray<T> const & var );


  template< typename T >
  static localIndex Pack( array<char> & buffer, set<T> const & var );

  template< typename T >
  static localIndex Pack( char*& buffer, set<T> const & var );

  template< typename T>
  static localIndex Unpack( char const *& buffer, std::set<T> & setToRead );


  template< typename T_KEY, typename T_VAL >
  static localIndex Pack( array<char> & buffer, std::map<T_KEY,T_VAL> const & var );

  template< typename T_KEY, typename T_VAL >
  static localIndex Pack( char*& buffer, const std::map<T_KEY,T_VAL>& var );


  ///@}

  /** @name Packing functions for filtered/masked collections of data
   */
  ///@{
  template< typename T, typename T_INDICES >
  static localIndex Pack( array<char> & buffer,
                          array<T> const & var,
                          T_INDICES const & indices );

  template< typename T, typename T_INDICES >
  static localIndex Pack( char*& buffer,
                          array<T> const & var,
                          T_INDICES const & indices );

  template< typename T, typename T_indices >
  static localIndex Unpack( char const *& buffer,
                            array<T>& var,
                            const T_indices& indices );

  ///@}









//
//
//
//
//
//  template< typename T_INDICES >
//  static localIndex Pack( char*& buffer,
//                          array<string> const & var,
//                          T_INDICES const & indices  );
//
//
//  template< typename T_CONTAINER >
//  static localIndex PackGlobal( array<char> & buffer,
//                                T_CONTAINER const & container,
//                                globalIndex_array const& localToGlobal );
//
//  template< typename T_CONTAINER >
//  static localIndex PackGlobal( char*& buffer,
//                                T_CONTAINER const & container,
//                                globalIndex_array const& localToGlobal );
//

//
//  template< typename T_indices >
//  localIndex Pack( const OneToOneRelation& relation,
//                   const T_indices& indices,
//                   const bool packGlobal )
//  {
//    return PrivatePackRelation( relation, indices, packGlobal );
//  }
//  template< typename T_indices >
//  localIndex Pack( const FixedOneToManyRelation& relation,
//                   const T_indices& indices,
//                   const bool packGlobal )
//  {
//    return PrivatePackRelation( relation, indices, packGlobal );
//  }
//
//  template< typename T_indices >
//  localIndex Pack( const OrderedVariableOneToManyRelation& relation,
//                   const T_indices& indices,
//                   const bool packGlobal )
//  {
//    return PrivatePackRelation( relation, indices, packGlobal );
//  }
//  template< typename T_indices >
//  localIndex Pack( const UnorderedVariableOneToManyRelation& relation,
//                   const T_indices& indices,
//                   const bool packGlobal )
//  {
//    return PrivatePackRelation( relation, indices, packGlobal );
//  }
//  template< typename T_indices >
//  localIndex Pack( const OrderedVariableOneToManyPairRelation& relation,
//                   const T_indices& indices,
//                   const bool packGlobal )
//  {
//    return PrivatePackRelation( relation, indices, packGlobal );
//  }
//  template< typename T_indices >
//  localIndex Pack( const UnorderedVariableOneToManyPairRelation& relation,
//                   const T_indices& indices,
//                   const bool packGlobal )
//  {
//    return PrivatePackRelation( relation, indices, packGlobal );
//  }
//
//

//
//  static localIndex Unpack( const char*& buffer, int& var ) { return PrivateUnpack(buffer,var); }
//  static localIndex Unpack( const char*& buffer, long int& var ) { return PrivateUnpack(buffer,var); }
//  static localIndex Unpack( const char*& buffer, long long int& var ) { return PrivateUnpack(buffer,var); }
//  static localIndex Unpack( const char*& buffer, realT& var ) { return PrivateUnpack(buffer,var); }
////  static localIndex Unpack( const char*& buffer, localIndex& var ) { return
//// PrivateUnpack(buffer,var); }
////  static localIndex Unpack( const char*& buffer, globalIndex& var ) { return
//// PrivateUnpack(buffer,var); }
//  static localIndex Unpack( const char*& buffer, R1Tensor& var ) { return PrivateUnpack(buffer,var); }
//  static localIndex Unpack( const char*& buffer, R2Tensor& var ) { return PrivateUnpack(buffer,var); }
//  static localIndex Unpack( const char*& buffer, R2SymTensor& var ) { return PrivateUnpack(buffer,var); }
////  static localIndex Unpack( const char*& buffer, size_t& var ) { return
//// PrivateUnpack(buffer,var); }
//  static localIndex Unpack( const char*& buffer, std::string& var );
//
//  template< typename T >
//  static localIndex UnpackSerialObject( const char*& buffer, T& object ) { return PrivateUnpack(buffer,object); }
//
////  template< int NINT, int NR0, int NR1, int NR2, int NR2S >
////  static localIndex UnpackSerialObject( const char*& buffer,
//// EncapsulatedObjectBase< NINT, NR0, NR1, NR2, NR2S >& object ) { return
//// PrivateUnpack(buffer, object ); }
//
//
//
//
//
//  static localIndex Unpack( const char*& buffer, set<int>& var ) { return PrivateUnpackSet(buffer,var); }
//  static localIndex Unpack( const char*& buffer, set<long int>& var ) { return PrivateUnpackSet(buffer,var); }
//  static localIndex Unpack( const char*& buffer, set<long long int>& var ) { return PrivateUnpackSet(buffer,var); }
//
//  template< typename T_KEY, typename T_VAL >
//  static localIndex Unpack( const char*& buffer, std::map<T_KEY,T_VAL>& var ) { return PrivateUnpackMap( buffer, var );}
//
//
//  static localIndex UnpackGlobal( const char*& buffer, const std::map<globalIndex,localIndex>& globalToLocal, localIndex_array& array );
//  static localIndex UnpackGlobal( const char*& buffer, const std::map<globalIndex,localIndex>& globalToLocal, lSet& set );
//
//  static localIndex Unpack( const char*& buffer, OneToOneRelation& relation, localIndex_array const& indices, const bool unpackGlobal );
//  static localIndex Unpack( const char*& buffer, FixedOneToManyRelation& relation, localIndex_array const& indices, const bool unpackGlobal );
//  static localIndex Unpack( const char*& buffer, OrderedVariableOneToManyRelation& relation, localIndex_array const& indices, const bool unpackGlobal );
//  static localIndex Unpack( const char*& buffer, UnorderedVariableOneToManyRelation& relation, localIndex_array const& indices, const bool unpackGlobal );
////  static localIndex Unpack( const char*& buffer,
//// OrderedVariableOneToManyPairRelation& relation, localIndex_array const&
//// indices, const bool unpackGlobal );
////  static localIndex Unpack( const char*& buffer,
//// UnorderedVariableOneToManyPairRelation& relation, localIndex_array const&
//// indices, const bool unpackGlobal );
//
//  localIndex PackFieldname( std::string fieldname ); // nb not passed by
//                                                     // reference
//  static localIndex PackFieldname( char*& buffer,  std::string fieldname );
//  static bool FieldnameMatchesIdString(std::string fieldname, std::string id);
//  static const unsigned sizeOfPackedFieldString;
//private:
//
//  //********************************************************************************************************************
//
//
//  //********************************************************************************************************************

//
//
//
//  //********************************************************************************************************************
//
//
//
//
//  template< typename T>
//  localIndex PrivatePackGlobal( const T& container, globalIndex_array const& localToGlobal );
//
//
//
//  template< typename T, typename T_indices >
//  localIndex PrivatePackGlobal( const T& container, const T_indices& indices, globalIndex_array const& localToGlobal );
//
//
//
//  template< typename T, typename T_indices >
//  struct PrivatePackRelationT
//  {
//    static localIndex Pack( bufvector& buffer, const T& relation, const T_indices& indices, const bool packGlobal );
//  };
//  template< typename T_indices >
//  struct PrivatePackRelationT<OneToOneRelation,T_indices>
//  {
//    static localIndex Pack( bufvector& buffer, const OneToOneRelation& relation, const T_indices& indices, const bool packGlobal );
//  };
//  template< typename T_indices >
//  struct PrivatePackRelationT<FixedOneToManyRelation,T_indices>
//  {
//    static localIndex Pack( bufvector& buffer, const FixedOneToManyRelation& relation, const T_indices& indices, const bool packGlobal );
//  };
//
//  template< typename T, typename T_indices >
//  localIndex PrivatePackRelation( const T& relation, const T_indices& indices, const bool packGlobal )
//  {
//    return PrivatePackRelationT<T,T_indices>::Pack(*this,relation,indices,packGlobal);
//  }
//
//
//
//  template< typename T >
//  static localIndex PrivateUnpackRelation( const char*& buffer, T& relation, localIndex_array const& indices, const bool unpackGlobal );
//


};



}

#include "CommBufferOps_inline.hpp"


#endif /* BUFVECTOR_H_ */
