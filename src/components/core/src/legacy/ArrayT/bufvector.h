//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Copyright (c) 2015, Lawrence Livermore National Security, LLC.
//  Produced at the Lawrence Livermore National Laboratory
//
//  GEOS Computational Framework - Core Package, Version 3.0.0
//
//  Written by:
//  Randolph Settgast (settgast1@llnl.gov)
//  Stuart Walsh(walsh24@llnl.gov)
//  Pengcheng Fu (fu4@llnl.gov)
//  Joshua White (white230@llnl.gov)
//  Chandrasekhar Annavarapu Srinivas
//  Eric Herbold
//  Michael Homel
//
//
//  All rights reserved.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S. Department of Energy (DOE). This work was produced at Lawrence Livermore 
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National Security, LLC nor any of their employees, makes any warranty, express or 
//     implied, or assumes any liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or 
//     process disclosed, or represents that its use would not infringe privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or services by trade name, trademark, manufacturer or otherwise does not 
//     necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or Lawrence Livermore National Security, 
//     LLC. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or Lawrence 
//     Livermore National Security, LLC, and shall not be used for advertising or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @file bufvector.h
 * @author settgast1
 * @date Apr 4, 2011
 */

#ifndef BUFVECTOR_H_
#define BUFVECTOR_H_

#include <string>
#include <string.h>
//#include "legacy/Common/typedefs.h"
#include "common/DataTypes.hpp"
#include <map>
#include "codingUtilities/Utilities.hpp"
#include "common/InterObjectRelation.hpp"
//#include "legacy/DataStructures/EncapsulatedObjects/EncapsulatedObjectBase.h"
#include "VectorT_derived.h"

#ifdef USE_ATK
#include <slic/slic.hpp>
#endif

namespace geosx
{
class bufvector: public VectorT<char>
{
public:
  bufvector(){}
  virtual ~bufvector(){}

  bufvector( bufvector const & ) = default;

  bufvector & operator=( bufvector const & ) = default;


  localIndex Pack( const int& var )           { return PrivatePack(var); }
  localIndex Pack( const long int& var )           { return PrivatePack(var); }
  localIndex Pack( const long long int& var )           { return PrivatePack(var); }
  localIndex Pack( const realT& var )         { return PrivatePack(var); }
//  localIndex Pack( const localIndex& var )    { return PrivatePack(var); }
//  localIndex Pack( const globalIndex& var )   { return PrivatePack(var); }
  localIndex Pack( const R1Tensor& var )      { return PrivatePack(var); }
  localIndex Pack( const R2Tensor& var )      { return PrivatePack(var); }
  localIndex Pack( const R2SymTensor& var )   { return PrivatePack(var); }
//  localIndex Pack( const size_t& var )   { return PrivatePack(var); }


  localIndex Pack( const std::string& var );

  localIndex Pack( array<int> const & var )              { return 0;}//PrivatePackArray(var); }
  localIndex Pack( array<long int> const & var )              { return 0;}//PrivatePackArray(var); }
  localIndex Pack( array<long long int> const & var )              { return 0;}//PrivatePackArray(var); }
  localIndex Pack( real64_array const& var )              { return PrivatePackArray(var); }
//  localIndex Pack( localIndex_array const& var )              { return PrivatePackArray(var); }
//  localIndex Pack( globalIndex_array const& var )              { return PrivatePackArray(var); }
  localIndex Pack( r1_array const& var )    { return PrivatePackArray(var); }
  localIndex Pack( r2_array const& var )    { return PrivatePackArray(var); }
  localIndex Pack( r2Sym_array const& var ) { return PrivatePackArray(var); }
  localIndex Pack( const sArray1d& var );


  localIndex Pack( const set<int>& var )              { return PrivatePackSet(var); }
  localIndex Pack( const set<long int>& var )              { return PrivatePackSet(var); }
  localIndex Pack( const set<long long int>& var )              { return PrivatePackSet(var); }

  template< typename T_KEY, typename T_VAL >
  localIndex Pack( const std::map<T_KEY,T_VAL>& var ) { return PrivatePackMap(var);}

  template< typename T > localIndex PackSerialObject( const T& object ) { return PrivatePack(object); }
//  template< int NINT, int NR0, int NR1, int NR2, int NR2S >
//  localIndex PackSerialObject( const EncapsulatedObjectBase< NINT, NR0, NR1, NR2, NR2S >& object ) { return PrivatePack( object ); }
//  template< int NINT, int NR0, int NR1, int NR2, int NR2S >
//  static localIndex PackSerialObject( char*& buffer, const EncapsulatedObjectBase< NINT, NR0, NR1, NR2, NR2S >& object ) { return PrivatePack(buffer, object ); }

  template< typename T_indices > localIndex Pack( integer_array const & array, const T_indices& indices ) { return PrivatePackArray( array, indices ); }
  template< typename T_indices > localIndex Pack( real64_array const& array, const T_indices& indices ) { return PrivatePackArray( array, indices ); }
//  template< typename T_indices > localIndex Pack( localIndex_array const& array, const T_indices& indices ) { return PrivatePackArray( array, indices ); }
  template< typename T_indices > localIndex Pack( globalIndex_array const& array, const T_indices& indices ) { return PrivatePackArray( array, indices ); }
  template< typename T_indices > localIndex Pack( r1_array const& array, const T_indices& indices ) { return PrivatePackArray( array, indices ); }
  template< typename T_indices > localIndex Pack( r2_array const& array, const T_indices& indices ) { return PrivatePackArray( array, indices ); }
  template< typename T_indices > localIndex Pack( r2Sym_array const& array, const T_indices& indices ) { return PrivatePackArray( array, indices ); }



  static localIndex Pack( char*& buffer, const int& var ) { return PrivatePack(buffer,var); }
  static localIndex Pack( char*& buffer, const realT& var ) { return PrivatePack(buffer,var); }
//  static localIndex Pack( char*& buffer, const localIndex& var ) { return PrivatePack(buffer,var); }
  static localIndex Pack( char*& buffer, const globalIndex& var ) { return PrivatePack(buffer,var); }
  static localIndex Pack( char*& buffer, const R1Tensor& var ) { return PrivatePack(buffer,var); }
  static localIndex Pack( char*& buffer, const R2Tensor& var ) { return PrivatePack(buffer,var); }
  static localIndex Pack( char*& buffer, const R2SymTensor& var ) { return PrivatePack(buffer,var); }
  template< typename T > static localIndex PackSerialObject( char*& buffer, const T& object ) { return PrivatePack(buffer,object); }


  static localIndex Pack( char*& buffer, const std::string& var );

  localIndex PackGlobal( localIndex_array const& container, globalIndex_array const& localToGlobal ) { return PrivatePackGlobal(container,localToGlobal); }
  localIndex PackGlobal( const lSet& container, globalIndex_array const& localToGlobal ) { return PrivatePackGlobal(container,localToGlobal); }

  template< typename T_indices > localIndex Pack( const OneToOneRelation& relation, const T_indices& indices, const bool packGlobal ) { return PrivatePackRelation( relation, indices, packGlobal );}
  template< typename T_indices > localIndex Pack( const FixedOneToManyRelation& relation, const T_indices& indices, const bool packGlobal ) { return PrivatePackRelation( relation, indices, packGlobal );}
  template< typename T_indices > localIndex Pack( const OrderedVariableOneToManyRelation& relation, const T_indices& indices, const bool packGlobal ) { return PrivatePackRelation( relation, indices, packGlobal );}
  template< typename T_indices > localIndex Pack( const UnorderedVariableOneToManyRelation& relation, const T_indices& indices, const bool packGlobal ) { return PrivatePackRelation( relation, indices, packGlobal );}
  template< typename T_indices > localIndex Pack( const OrderedVariableOneToManyPairRelation& relation, const T_indices& indices, const bool packGlobal ) { return PrivatePackRelation( relation, indices, packGlobal );}
  template< typename T_indices > localIndex Pack( const UnorderedVariableOneToManyPairRelation& relation, const T_indices& indices, const bool packGlobal ) { return PrivatePackRelation( relation, indices, packGlobal );}


  static localIndex Unpack( const char*& buffer, int& var ) { return PrivateUnpack(buffer,var); }
  static localIndex Unpack( const char*& buffer, long int& var ) { return PrivateUnpack(buffer,var); }
  static localIndex Unpack( const char*& buffer, long long int& var ) { return PrivateUnpack(buffer,var); }
  static localIndex Unpack( const char*& buffer, realT& var ) { return PrivateUnpack(buffer,var); }
//  static localIndex Unpack( const char*& buffer, localIndex& var ) { return PrivateUnpack(buffer,var); }
//  static localIndex Unpack( const char*& buffer, globalIndex& var ) { return PrivateUnpack(buffer,var); }
  static localIndex Unpack( const char*& buffer, R1Tensor& var ) { return PrivateUnpack(buffer,var); }
  static localIndex Unpack( const char*& buffer, R2Tensor& var ) { return PrivateUnpack(buffer,var); }
  static localIndex Unpack( const char*& buffer, R2SymTensor& var ) { return PrivateUnpack(buffer,var); }
//  static localIndex Unpack( const char*& buffer, size_t& var ) { return PrivateUnpack(buffer,var); }
  static localIndex Unpack( const char*& buffer, std::string& var );

  template< typename T >
  static localIndex UnpackSerialObject( const char*& buffer, T& object ) { return PrivateUnpack(buffer,object); }

//  template< int NINT, int NR0, int NR1, int NR2, int NR2S >
//  static localIndex UnpackSerialObject( const char*& buffer, EncapsulatedObjectBase< NINT, NR0, NR1, NR2, NR2S >& object ) { return PrivateUnpack(buffer, object ); }



  static localIndex Unpack( const char*& buffer, integer_array& var ) { return 0;}//PrivateUnpackArray(buffer,var); }
  static localIndex Unpack( const char*& buffer, real64_array& var ) { return PrivateUnpackArray(buffer,var); }
//  static localIndex Unpack( const char*& buffer, localIndex_array& var ) { return PrivateUnpackArray(buffer,var); }
  static localIndex Unpack( const char*& buffer, globalIndex_array& var ) { return PrivateUnpackArray(buffer,var); }
  static localIndex Unpack( const char*& buffer, r1_array& var ) { return PrivateUnpackArray(buffer,var); }
  static localIndex Unpack( const char*& buffer, r2_array& var ) { return PrivateUnpackArray(buffer,var); }
  static localIndex Unpack( const char*& buffer, r2Sym_array& var ) { return PrivateUnpackArray(buffer,var); }
  static localIndex Unpack( const char*& buffer, sArray1d& var );

  static localIndex Unpack( const char*& buffer, set<int>& var ) { return PrivateUnpackSet(buffer,var); }
  static localIndex Unpack( const char*& buffer, set<long int>& var ) { return PrivateUnpackSet(buffer,var); }
  static localIndex Unpack( const char*& buffer, set<long long int>& var ) { return PrivateUnpackSet(buffer,var); }

  template< typename T_KEY, typename T_VAL >
  static localIndex Unpack( const char*& buffer, std::map<T_KEY,T_VAL>& var ) { return PrivateUnpackMap( buffer, var );}


  static localIndex UnpackGlobal( const char*& buffer, const std::map<globalIndex,localIndex>& globalToLocal, localIndex_array& array );
  static localIndex UnpackGlobal( const char*& buffer, const std::map<globalIndex,localIndex>& globalToLocal, lSet& set );

  static localIndex Unpack( const char*& buffer, OneToOneRelation& relation, localIndex_array const& indices, const bool unpackGlobal );
  static localIndex Unpack( const char*& buffer, FixedOneToManyRelation& relation, localIndex_array const& indices, const bool unpackGlobal );
  static localIndex Unpack( const char*& buffer, OrderedVariableOneToManyRelation& relation, localIndex_array const& indices, const bool unpackGlobal );
  static localIndex Unpack( const char*& buffer, UnorderedVariableOneToManyRelation& relation, localIndex_array const& indices, const bool unpackGlobal );
//  static localIndex Unpack( const char*& buffer, OrderedVariableOneToManyPairRelation& relation, localIndex_array const& indices, const bool unpackGlobal );
//  static localIndex Unpack( const char*& buffer, UnorderedVariableOneToManyPairRelation& relation, localIndex_array const& indices, const bool unpackGlobal );

  localIndex PackFieldname( std::string fieldname ); // nb not passed by reference
  static localIndex PackFieldname( char*& buffer,  std::string fieldname );
  static bool FieldnameMatchesIdString(std::string fieldname, std::string id);
  static const unsigned sizeOfPackedFieldString;
private:

  //********************************************************************************************************************
  /**
   * @author settgast
   * @tparam type of data to pack
   * @param var data to pack
   * @return size (in bytes) of packed data
   */
  template< typename T >
  localIndex PrivatePack( const T& var )
  {

    localIndex sizeOfPackedChars = sizeof(T);

  #if 1
    this->resize(this->size() + sizeOfPackedChars);

    char* p_buffer = &(this->back()) - sizeOfPackedChars + 1;

    memcpy( p_buffer, &var, sizeOfPackedChars );


  #else
    const char* cvar = (const char*)(&var);
    for( localIndex i=0 ; i<sizeOfPackedChars ; ++i )
    {
      this->push_back(*(cvar++));
    }
  #endif
    return sizeOfPackedChars;

  }

  /**
   * @author settgast
   * @param buffer
   * @param var
   * @return size (in bytes) of unpacked data
   */
  template< typename T>
  static localIndex PrivateUnpack( const char*& buffer, T& var )
  {
    localIndex sizeOfUnpackedChars = sizeof(T);

    const T* const tbuffer = reinterpret_cast<const T*>(buffer);
    var = *tbuffer;
    buffer += sizeOfUnpackedChars;

    return sizeOfUnpackedChars;
  }

  //********************************************************************************************************************
  /**
   * @author settgast
   * @param container
   * @return
   */
  template< typename T >
  localIndex PrivatePackArray( const array<T>& container )
  {
    localIndex sizeOfPackedChars = 0;

    const localIndex length = container.size();
    localIndex sizeOfPackedArrayChars = length*sizeof(T);


    sizeOfPackedChars += this->Pack( length );

    this->resize( this->size() + sizeOfPackedArrayChars );

    char* p_buffer = &(this->back()) - sizeOfPackedArrayChars + 1;
    memcpy( p_buffer, container.data(), sizeOfPackedArrayChars );
    sizeOfPackedChars += sizeOfPackedArrayChars;

    return sizeOfPackedChars;
  }




  template< typename T>
  static localIndex PrivateUnpackArray( const char*& buffer, array<T>& array )
  {

    localIndex sizeOfUnpackedChars = 0;

    localIndex array_length;
    sizeOfUnpackedChars += Unpack( buffer, array_length );
    array.resize(array_length);
    localIndex length = array_length * sizeof(T);

    memcpy( array.data() , buffer, length );
    buffer += length;
    sizeOfUnpackedChars += length;


    return sizeOfUnpackedChars;
  }

  //********************************************************************************************************************
  template< typename T, typename T_indices >
  localIndex PrivatePackArray( const array<T>& container, const T_indices& indices )
  {
    localIndex sizeOfPackedChars = 0;

    const typename T_indices::size_type length = indices.size();
    localIndex sizeOfPackedArrayChars = length*sizeof(T);
    this->resize( this->size() + sizeOfPackedArrayChars );

    sizeOfPackedChars += this->Pack( length );

    char* p_buffer = &(this->back()) - sizeOfPackedArrayChars + 1;
    for( typename T_indices::const_iterator i=indices.begin() ; i!=indices.end() ; ++i )
    {
      memcpy( p_buffer, &(container[*i]), sizeof(T) );
      p_buffer += sizeof(T);
    }
    sizeOfPackedChars += sizeOfPackedArrayChars;

    return sizeOfPackedChars;
  }


  template< typename T, typename T_indices >
  static localIndex PrivateUnpackArray( const char*& buffer, array<T>& array, const T_indices& indices )
  {

    localIndex sizeOfUnpackedChars = 0;

    typename T_indices::size_type array_length;

    sizeOfUnpackedChars += Unpack( buffer, array_length );

    if( array_length != indices.size() )
    {
#ifdef USE_ATK
      SLIC_ERROR("bufvector::PrivateUnpackArray(): incorrect number of data");
#endif
//      throw GPException("bufvector::PrivateUnpackArray(): incorrect number of data");
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

  //********************************************************************************************************************
  template< typename T >
  localIndex PrivatePackSet( const std::set<T>& var )
  {

    localIndex sizeOfPackedChars = 0;

    const localIndex length = var.size();

    sizeOfPackedChars += Pack( length );

    for( typename std::set<T>::const_iterator i=var.begin() ; i!=var.end() ; ++i )
    {
      sizeOfPackedChars += this->Pack(*i);
    }

    return sizeOfPackedChars;
  }



  template< typename T>
  static localIndex PrivateUnpackSet( const char*& buffer, std::set<T>& setToRead )
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
  localIndex PrivatePackMap( const std::map<T_KEY,T_VAL>& var )
  {

    localIndex sizeOfPackedChars = 0;

    const typename std::map<T_KEY,T_VAL>::size_type length = var.size();

    sizeOfPackedChars += Pack( length );



    for( typename std::map<T_KEY,T_VAL>::const_iterator i=var.begin() ; i!=var.end() ; ++i )
    {
      sizeOfPackedChars += this->Pack(i->first);
      sizeOfPackedChars += this->Pack(i->second);
    }

    return sizeOfPackedChars;
  }



  template< typename T_KEY, typename T_VAL >
  static localIndex PrivateUnpackMap( const char*& buffer, std::map<T_KEY,T_VAL>& map )
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

  template< typename T>
  static localIndex PrivatePack( char*& buffer, const T& var )
  {
    memcpy( buffer, &var, sizeof(T) );
    buffer += sizeof(T);
    return sizeof(T);
  }








  template< typename T>
  localIndex PrivatePackGlobal( const T& container, globalIndex_array const& localToGlobal );



  template< typename T, typename T_indices >
  localIndex PrivatePackGlobal( const T& container, const T_indices& indices, globalIndex_array const& localToGlobal );






  template< typename T, typename T_indices >
  struct PrivatePackRelationT
  {
    static localIndex Pack( bufvector& buffer, const T& relation, const T_indices& indices, const bool packGlobal );
  };
  template< typename T_indices >
  struct PrivatePackRelationT<OneToOneRelation,T_indices>
  {
    static localIndex Pack( bufvector& buffer, const OneToOneRelation& relation, const T_indices& indices, const bool packGlobal );
  };
  template< typename T_indices >
  struct PrivatePackRelationT<FixedOneToManyRelation,T_indices>
  {
    static localIndex Pack( bufvector& buffer, const FixedOneToManyRelation& relation, const T_indices& indices, const bool packGlobal );
  };

  template< typename T, typename T_indices >
  localIndex PrivatePackRelation( const T& relation, const T_indices& indices, const bool packGlobal )
  {
    return PrivatePackRelationT<T,T_indices>::Pack(*this,relation,indices,packGlobal);
  }



  template< typename T >
  static localIndex PrivateUnpackRelation( const char*& buffer, T& relation, localIndex_array const& indices, const bool unpackGlobal );



};


inline localIndex bufvector::Pack( const sArray1d& container )
{
  localIndex sizeOfPackedChars = 0;
  const localIndex arrayLength = container.size();

  sizeOfPackedChars += this->Pack( arrayLength );
  for( sArray1d::const_iterator i=container.begin() ; i!=container.end() ; ++i )
  {
    sizeOfPackedChars += this->Pack(*i);
  }

  return sizeOfPackedChars;
}



inline localIndex bufvector::Unpack( const char*& buffer, sArray1d& array )
{
  array.clear();
  localIndex sizeOfUnpackedChars = 0;

  localIndex arrayLength;
  sizeOfUnpackedChars += Unpack( buffer, arrayLength );
  array.resize(arrayLength);

  for( auto i=0 ; i<arrayLength ; ++i )
  {
    sizeOfUnpackedChars += Unpack( buffer, array[i] );
  }

  return sizeOfUnpackedChars;
}


//********************************************************************************************************************
template< typename T>
localIndex bufvector::PrivatePackGlobal( const T& container, globalIndex_array const& localToGlobal )
{
  const localIndex length = container.size();
  localIndex sizeOfPackedChars = 0;

//  std::cout<<"container.size() = "<<length<<std::endl;

  sizeOfPackedChars += this->Pack(length);
  for( typename T::const_iterator i=container.begin() ; i!=container.end() ; ++i )
  {
    const globalIndex temp = localToGlobal[*i];
    sizeOfPackedChars += this->Pack( temp );
  }
  return sizeOfPackedChars;
}


inline localIndex bufvector::UnpackGlobal( const char*& buffer, const std::map<globalIndex,localIndex>& globalToLocal, localIndex_array& array )
{
  array.clear();
  localIndex sizeOfUnpackedChars = 0;


  localIndex array_length;
  sizeOfUnpackedChars += Unpack( buffer, array_length );

  array.resize(array_length);
  globalIndex_array temp(array_length);
  localIndex length = array_length * sizeof(localIndex);

  memcpy( temp.data() , buffer, length );
  buffer += length;
  sizeOfUnpackedChars += length;


  for( localIndex i=0 ; i<array_length ; ++i )
  {
    const localIndex li = stlMapLookup( globalToLocal, temp[i] );
    array[i] = li;

  }

  return sizeOfUnpackedChars;
}

inline localIndex bufvector::UnpackGlobal( const char*& buffer, const std::map<globalIndex,localIndex>& globalToLocal, lSet& array )
{
  array.clear();
  localIndex sizeOfUnpackedChars = 0;


  localIndex array_length;
  sizeOfUnpackedChars += Unpack( buffer, array_length );

//  std::cout<<"array_length = "<<array_length<<std::endl;

  globalIndex_array temp(array_length);
  localIndex length = array_length * sizeof(localIndex);

  memcpy( temp.data() , buffer, length );
  buffer += length;
  sizeOfUnpackedChars += length;


  for( auto i=0 ; i<array_length ; ++i )
  {
    const localIndex li = stlMapLookup( globalToLocal, temp[i] );
    array.insert( li );
  }

  return sizeOfUnpackedChars;
}




//**********************************************************************************************************************

template< typename T, typename T_indices >
inline localIndex bufvector::PrivatePackRelationT<T,T_indices>::Pack( bufvector& buffer, const T& relation, const T_indices& indices, const bool packGlobal )
{
  localIndex sizeOfPackedChars = 0;

  if( packGlobal )
  {
    globalIndex_array const& localToGlobal = relation.RelatedObjectLocalToGlobal();
    for( typename T_indices::const_iterator i = indices.begin() ; i != indices.end() ; ++i )
    {
      sizeOfPackedChars += buffer.PackGlobal( relation[*i], localToGlobal );
    }
  }
  else
  {
    for( typename T_indices::const_iterator i = indices.begin() ; i != indices.end() ; ++i )
    {
      sizeOfPackedChars += buffer.Pack( relation[*i] );
    }
  }

  return sizeOfPackedChars;
}
template<typename T_indices>
inline localIndex bufvector::PrivatePackRelationT<OneToOneRelation,T_indices>::Pack( bufvector& buffer, const OneToOneRelation& relation, const T_indices& indices, const bool packGlobal )
{
  localIndex sizeOfPackedChars = 0;

  if( packGlobal )
  {
    globalIndex_array const& localToGlobal = relation.RelatedObjectLocalToGlobal();
    for( typename T_indices::const_iterator i = indices.begin() ; i != indices.end() ; ++i )
    {
      sizeOfPackedChars += buffer.Pack( localToGlobal[relation[*i]] );
    }
  }
  else
  {
    for( typename T_indices::const_iterator i = indices.begin() ; i != indices.end() ; ++i )
    {
      sizeOfPackedChars += buffer.Pack( relation[*i] );
    }
  }
  return sizeOfPackedChars;
}
template<typename T_indices>
inline localIndex bufvector::PrivatePackRelationT<FixedOneToManyRelation,T_indices>::Pack( bufvector& buffer, const FixedOneToManyRelation& relation, const T_indices& indices, const bool packGlobal )
{
  localIndex sizeOfPackedChars = 0;

  sizeOfPackedChars += buffer.Pack( static_cast<localIndex>(relation.Dimension(1)) );

  if( packGlobal )
  {
    globalIndex_array const& localToGlobal = relation.RelatedObjectLocalToGlobal();
    for( typename T_indices::const_iterator i = indices.begin() ; i != indices.end() ; ++i )
    {
      for( auto j=0u ; j<relation.Dimension(1) ; ++j )
      {
        sizeOfPackedChars += buffer.Pack( localToGlobal[relation[*i][j]] );
      }
    }
  }
  else
  {
    for( typename T_indices::const_iterator i = indices.begin() ; i != indices.end() ; ++i )
    {
      for( auto j=0u ; j<relation.Dimension(1) ; ++j )
      {
        sizeOfPackedChars += buffer.Pack( relation[*i][j] );
      }
    }
  }
  return sizeOfPackedChars;
}


template< typename T >
inline localIndex bufvector::PrivateUnpackRelation( const char*& buffer, T& relation, localIndex_array const& indices, const bool unpackGlobal )
{
  localIndex sizeOfUnpackedChars = 0;

  if( unpackGlobal )
  {
    const std::map<globalIndex,localIndex>& globalToLocal = relation.RelatedObjectGlobalToLocal();
    for( localIndex_array::const_iterator i = indices.begin() ; i != indices.end() ; ++i )
    {
      sizeOfUnpackedChars += bufvector::UnpackGlobal( buffer, globalToLocal, relation[*i] );
    }
  }
  else
  {
    for( localIndex_array::const_iterator i = indices.begin() ; i != indices.end() ; ++i )
    {
      sizeOfUnpackedChars += bufvector::Unpack( buffer, relation[*i] );
    }
  }

  return sizeOfUnpackedChars;
}

template<>
inline localIndex bufvector::PrivateUnpackRelation( const char*& buffer, OneToOneRelation& relation, localIndex_array const& indices, const bool unpackGlobal )
{
  localIndex sizeOfUnpackedChars = 0;

  if( unpackGlobal )
  {
    const std::map<globalIndex,localIndex>& globalToLocal = relation.RelatedObjectGlobalToLocal();
    for( localIndex_array::const_iterator i = indices.begin() ; i != indices.end() ; ++i )
    {
      globalIndex temp;
      sizeOfUnpackedChars += bufvector::Unpack( buffer, temp );
      relation[*i] = stlMapLookup(globalToLocal,temp);
    }
  }
  else
  {
    for( localIndex_array::const_iterator i = indices.begin() ; i != indices.end() ; ++i )
    {
      sizeOfUnpackedChars += bufvector::Unpack( buffer, relation[*i] );
    }
  }

  return sizeOfUnpackedChars;
}

template<>
inline localIndex bufvector::PrivateUnpackRelation( const char*& buffer, FixedOneToManyRelation& relation, localIndex_array const& indices, const bool unpackGlobal )
{
  localIndex sizeOfUnpackedChars = 0;

  localIndex dimension;
  sizeOfUnpackedChars += bufvector::Unpack( buffer, dimension );

  if( dimension != relation.Dimension(1)) {
#ifdef USE_ATK
    SLIC_ERROR("bufvector::PrivateUnpackRelation(): mismatched dimension");
#endif
  }

  if( unpackGlobal )
  {
    const std::map<globalIndex,localIndex>& globalToLocal = relation.RelatedObjectGlobalToLocal();
    for( localIndex_array::const_iterator i = indices.begin() ; i != indices.end() ; ++i )
    {
      for( localIndex j=0 ; j<dimension ; ++j )
      {
        globalIndex temp;
        sizeOfUnpackedChars += bufvector::Unpack( buffer, temp );
        relation[*i][j] = stlMapLookup(globalToLocal,temp);
      }
    }
  }
  else
  {
    for( localIndex_array::const_iterator i = indices.begin() ; i != indices.end() ; ++i )
    {
      for( localIndex j=0 ; j<dimension ; ++j )
      {
        sizeOfUnpackedChars += bufvector::Unpack( buffer, relation[*i][j] );
      }
    }
  }

  return sizeOfUnpackedChars;
}


inline localIndex bufvector::Unpack( const char*& buffer, OneToOneRelation& relation, localIndex_array const& indices, const bool unpackGlobal ) { return PrivateUnpackRelation( buffer, relation, indices, unpackGlobal ); }
inline localIndex bufvector::Unpack( const char*& buffer, FixedOneToManyRelation& relation, localIndex_array const& indices, const bool unpackGlobal ) { return PrivateUnpackRelation( buffer, relation, indices, unpackGlobal ); }
inline localIndex bufvector::Unpack( const char*& buffer, OrderedVariableOneToManyRelation& relation, localIndex_array const& indices, const bool unpackGlobal ) { return PrivateUnpackRelation( buffer, relation, indices, unpackGlobal ); }
inline localIndex bufvector::Unpack( const char*& buffer, UnorderedVariableOneToManyRelation& relation, localIndex_array const& indices, const bool unpackGlobal ) { return PrivateUnpackRelation( buffer, relation, indices, unpackGlobal ); }
//inline localIndex bufvector::Unpack( const char*& buffer, OrderedVariableOneToManyPairRelation& relation, localIndex_array const& indices, const bool unpackGlobal ) { return PrivateUnpackRelation( buffer, relation, indices, unpackGlobal ); }
//inline localIndex bufvector::Unpack( const char*& buffer, UnorderedVariableOneToManyPairRelation& relation, localIndex_array const& indices, const bool unpackGlobal ) { return PrivateUnpackRelation( buffer, relation, indices, unpackGlobal ); }





//**********************************************************************************************************************
inline localIndex bufvector::Pack( const std::string& var )
{
  localIndex sizeOfPackedChars = var.size();

  this->PrivatePack( sizeOfPackedChars );

  const char* cvar = var.data();
  for( localIndex i=0 ; i<sizeOfPackedChars ; ++i )
  {
    this->push_back(*(cvar++));
  }

  sizeOfPackedChars += sizeof( localIndex );
  return sizeOfPackedChars;
}

inline localIndex bufvector::Pack(char*& buffer,  const std::string& var )
{
  localIndex sizeOfPackedChars = var.size();

  PrivatePack(buffer,  sizeOfPackedChars );

  for( localIndex i=0 ; i<sizeOfPackedChars ; ++i )
  {
	*buffer = var[i];
	buffer++;
  }

  sizeOfPackedChars += sizeof( localIndex );
  return sizeOfPackedChars;
}

inline localIndex bufvector::Unpack( const char*& buffer, std::string& var )
{
  localIndex sizeOfUnpackedChars = 0;
  localIndex stringsize = 0;

  sizeOfUnpackedChars += PrivateUnpack( buffer, stringsize );

  var.assign( buffer, stringsize );
  buffer += stringsize;
  sizeOfUnpackedChars += stringsize;

  return sizeOfUnpackedChars;
}

}













#endif /* BUFVECTOR_H_ */
