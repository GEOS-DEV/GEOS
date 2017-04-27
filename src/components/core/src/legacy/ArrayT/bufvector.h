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
#include "legacy/Common/typedefs.h"
#include <map>
#include "codingUtilities/Utilities.hpp"
#include "legacy/DataStructures/InterObjectRelation.h"
#include "legacy/DataStructures/EncapsulatedObjects/EncapsulatedObjectBase.h"

class bufvector: public VectorT<char>
{
public:
  bufvector(){}
  virtual ~bufvector(){}

  bufvector( bufvector const & ) = default;

  bufvector & operator=( bufvector const & ) = default;


  unsigned int Pack( const int& var )           { return PrivatePack(var); }
  unsigned int Pack( const realT& var )         { return PrivatePack(var); }
//  unsigned int Pack( const localIndex& var )    { return PrivatePack(var); }
  unsigned int Pack( const globalIndex& var )   { return PrivatePack(var); }
  unsigned int Pack( const R1Tensor& var )      { return PrivatePack(var); }
  unsigned int Pack( const R2Tensor& var )      { return PrivatePack(var); }
  unsigned int Pack( const R2SymTensor& var )   { return PrivatePack(var); }
  unsigned int Pack( const size_t& var )   { return PrivatePack(var); }


  unsigned int Pack( const std::string& var );

  unsigned int Pack( const iArray1d& var )              { return PrivatePackArray(var); }
  unsigned int Pack( const rArray1d& var )              { return PrivatePackArray(var); }
//  unsigned int Pack( const lArray1d& var )              { return PrivatePackArray(var); }
  unsigned int Pack( const gArray1d& var )              { return PrivatePackArray(var); }
  unsigned int Pack( const Array1dT<R1Tensor>& var )    { return PrivatePackArray(var); }
  unsigned int Pack( const Array1dT<R2Tensor>& var )    { return PrivatePackArray(var); }
  unsigned int Pack( const Array1dT<R2SymTensor>& var ) { return PrivatePackArray(var); }
  unsigned int Pack( const sArray1d& var );


  unsigned int Pack( const lSet& var )              { return PrivatePackSet(var); }
  unsigned int Pack( const gSet& var )              { return PrivatePackSet(var); }

  template< typename T_KEY, typename T_VAL >
  unsigned int Pack( const std::map<T_KEY,T_VAL>& var ) { return PrivatePackMap(var);}

  template< typename T > unsigned int PackSerialObject( const T& object ) { return PrivatePack(object); }
  template< int NINT, int NR0, int NR1, int NR2, int NR2S >
  unsigned int PackSerialObject( const EncapsulatedObjectBase< NINT, NR0, NR1, NR2, NR2S >& object ) { return PrivatePack( object ); }
  template< int NINT, int NR0, int NR1, int NR2, int NR2S >
  static unsigned int PackSerialObject( char*& buffer, const EncapsulatedObjectBase< NINT, NR0, NR1, NR2, NR2S >& object ) { return PrivatePack(buffer, object ); }

  template< typename T_indices > unsigned int Pack( const iArray1d& array, const T_indices& indices ) { return PrivatePackArray( array, indices ); }
  template< typename T_indices > unsigned int Pack( const rArray1d& array, const T_indices& indices ) { return PrivatePackArray( array, indices ); }
//  template< typename T_indices > unsigned int Pack( const lArray1d& array, const T_indices& indices ) { return PrivatePackArray( array, indices ); }
  template< typename T_indices > unsigned int Pack( const gArray1d& array, const T_indices& indices ) { return PrivatePackArray( array, indices ); }
  template< typename T_indices > unsigned int Pack( const Array1dT<R1Tensor>& array, const T_indices& indices ) { return PrivatePackArray( array, indices ); }
  template< typename T_indices > unsigned int Pack( const Array1dT<R2Tensor>& array, const T_indices& indices ) { return PrivatePackArray( array, indices ); }
  template< typename T_indices > unsigned int Pack( const Array1dT<R2SymTensor>& array, const T_indices& indices ) { return PrivatePackArray( array, indices ); }



  static unsigned int Pack( char*& buffer, const int& var ) { return PrivatePack(buffer,var); }
  static unsigned int Pack( char*& buffer, const realT& var ) { return PrivatePack(buffer,var); }
//  static unsigned int Pack( char*& buffer, const localIndex& var ) { return PrivatePack(buffer,var); }
  static unsigned int Pack( char*& buffer, const globalIndex& var ) { return PrivatePack(buffer,var); }
  static unsigned int Pack( char*& buffer, const R1Tensor& var ) { return PrivatePack(buffer,var); }
  static unsigned int Pack( char*& buffer, const R2Tensor& var ) { return PrivatePack(buffer,var); }
  static unsigned int Pack( char*& buffer, const R2SymTensor& var ) { return PrivatePack(buffer,var); }
  template< typename T > static unsigned int PackSerialObject( char*& buffer, const T& object ) { return PrivatePack(buffer,object); }


  static unsigned int Pack( char*& buffer, const std::string& var );

  unsigned int PackGlobal( const lArray1d& container, const gArray1d& localToGlobal ) { return PrivatePackGlobal(container,localToGlobal); }
  unsigned int PackGlobal( const lSet& container, const gArray1d& localToGlobal ) { return PrivatePackGlobal(container,localToGlobal); }

  template< typename T_indices > unsigned int Pack( const OneToOneRelation& relation, const T_indices& indices, const bool packGlobal ) { return PrivatePackRelation( relation, indices, packGlobal );}
  template< typename T_indices > unsigned int Pack( const FixedOneToManyRelation& relation, const T_indices& indices, const bool packGlobal ) { return PrivatePackRelation( relation, indices, packGlobal );}
  template< typename T_indices > unsigned int Pack( const OrderedVariableOneToManyRelation& relation, const T_indices& indices, const bool packGlobal ) { return PrivatePackRelation( relation, indices, packGlobal );}
  template< typename T_indices > unsigned int Pack( const UnorderedVariableOneToManyRelation& relation, const T_indices& indices, const bool packGlobal ) { return PrivatePackRelation( relation, indices, packGlobal );}
  template< typename T_indices > unsigned int Pack( const OrderedVariableOneToManyPairRelation& relation, const T_indices& indices, const bool packGlobal ) { return PrivatePackRelation( relation, indices, packGlobal );}
  template< typename T_indices > unsigned int Pack( const UnorderedVariableOneToManyPairRelation& relation, const T_indices& indices, const bool packGlobal ) { return PrivatePackRelation( relation, indices, packGlobal );}


  static unsigned int Unpack( const char*& buffer, int& var ) { return PrivateUnpack(buffer,var); }
  static unsigned int Unpack( const char*& buffer, realT& var ) { return PrivateUnpack(buffer,var); }
//  static unsigned int Unpack( const char*& buffer, localIndex& var ) { return PrivateUnpack(buffer,var); }
  static unsigned int Unpack( const char*& buffer, globalIndex& var ) { return PrivateUnpack(buffer,var); }
  static unsigned int Unpack( const char*& buffer, R1Tensor& var ) { return PrivateUnpack(buffer,var); }
  static unsigned int Unpack( const char*& buffer, R2Tensor& var ) { return PrivateUnpack(buffer,var); }
  static unsigned int Unpack( const char*& buffer, R2SymTensor& var ) { return PrivateUnpack(buffer,var); }
  static unsigned int Unpack( const char*& buffer, size_t& var ) { return PrivateUnpack(buffer,var); }
  static unsigned int Unpack( const char*& buffer, std::string& var );

  template< typename T >
  static unsigned int UnpackSerialObject( const char*& buffer, T& object ) { return PrivateUnpack(buffer,object); }

  template< int NINT, int NR0, int NR1, int NR2, int NR2S >
  static unsigned int UnpackSerialObject( const char*& buffer, EncapsulatedObjectBase< NINT, NR0, NR1, NR2, NR2S >& object ) { return PrivateUnpack(buffer, object ); }



  static unsigned int Unpack( const char*& buffer, iArray1d& var ) { return PrivateUnpackArray(buffer,var); }
  static unsigned int Unpack( const char*& buffer, rArray1d& var ) { return PrivateUnpackArray(buffer,var); }
//  static unsigned int Unpack( const char*& buffer, lArray1d& var ) { return PrivateUnpackArray(buffer,var); }
  static unsigned int Unpack( const char*& buffer, gArray1d& var ) { return PrivateUnpackArray(buffer,var); }
  static unsigned int Unpack( const char*& buffer, Array1dT<R1Tensor>& var ) { return PrivateUnpackArray(buffer,var); }
  static unsigned int Unpack( const char*& buffer, Array1dT<R2Tensor>& var ) { return PrivateUnpackArray(buffer,var); }
  static unsigned int Unpack( const char*& buffer, Array1dT<R2SymTensor>& var ) { return PrivateUnpackArray(buffer,var); }
  static unsigned int Unpack( const char*& buffer, sArray1d& var );

  static unsigned int Unpack( const char*& buffer, lSet& var ) { return PrivateUnpackSet(buffer,var); }
  static unsigned int Unpack( const char*& buffer, gSet& var ) { return PrivateUnpackSet(buffer,var); }

  template< typename T_KEY, typename T_VAL >
  static unsigned int Unpack( const char*& buffer, std::map<T_KEY,T_VAL>& var ) { return PrivateUnpackMap( buffer, var );}


  static unsigned int UnpackGlobal( const char*& buffer, const std::map<globalIndex,localIndex>& globalToLocal, lArray1d& array );
  static unsigned int UnpackGlobal( const char*& buffer, const std::map<globalIndex,localIndex>& globalToLocal, lSet& set );

  static unsigned int Unpack( const char*& buffer, OneToOneRelation& relation, const lArray1d& indices, const bool unpackGlobal );
  static unsigned int Unpack( const char*& buffer, FixedOneToManyRelation& relation, const lArray1d& indices, const bool unpackGlobal );
  static unsigned int Unpack( const char*& buffer, OrderedVariableOneToManyRelation& relation, const lArray1d& indices, const bool unpackGlobal );
  static unsigned int Unpack( const char*& buffer, UnorderedVariableOneToManyRelation& relation, const lArray1d& indices, const bool unpackGlobal );
//  static unsigned int Unpack( const char*& buffer, OrderedVariableOneToManyPairRelation& relation, const lArray1d& indices, const bool unpackGlobal );
//  static unsigned int Unpack( const char*& buffer, UnorderedVariableOneToManyPairRelation& relation, const lArray1d& indices, const bool unpackGlobal );

  unsigned int PackFieldname( std::string fieldname ); // nb not passed by reference
  static unsigned int PackFieldname( char*& buffer,  std::string fieldname );
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
  unsigned int PrivatePack( const T& var )
  {

    unsigned int sizeOfPackedChars = sizeof(T);

  #if 1
    this->resize(this->size() + sizeOfPackedChars);

    char* p_buffer = &(this->back()) - sizeOfPackedChars + 1;

    memcpy( p_buffer, &var, sizeOfPackedChars );


  #else
    const char* cvar = (const char*)(&var);
    for( unsigned int i=0 ; i<sizeOfPackedChars ; ++i )
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
  static unsigned int PrivateUnpack( const char*& buffer, T& var )
  {
    unsigned int sizeOfUnpackedChars = sizeof(T);

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
  unsigned int PrivatePackArray( const Array1dT<T>& container )
  {
    unsigned int sizeOfPackedChars = 0;

    const typename Array1dT<T>::size_type length = container.size();
    unsigned int sizeOfPackedArrayChars = length*sizeof(T);


    sizeOfPackedChars += this->Pack( length );

    this->resize( this->size() + sizeOfPackedArrayChars );

    char* p_buffer = &(this->back()) - sizeOfPackedArrayChars + 1;
    memcpy( p_buffer, container.data(), sizeOfPackedArrayChars );
    sizeOfPackedChars += sizeOfPackedArrayChars;

    return sizeOfPackedChars;
  }




  template< typename T>
  static unsigned int PrivateUnpackArray( const char*& buffer, Array1dT<T>& array )
  {

    unsigned int sizeOfUnpackedChars = 0;

    typename Array1dT<T>::size_type array_length;
    sizeOfUnpackedChars += Unpack( buffer, array_length );
    array.resize(array_length);
    unsigned int length = array_length * sizeof(T);

    memcpy( array.data() , buffer, length );
    buffer += length;
    sizeOfUnpackedChars += length;


    return sizeOfUnpackedChars;
  }

  //********************************************************************************************************************
  template< typename T, typename T_indices >
  unsigned int PrivatePackArray( const Array1dT<T>& container, const T_indices& indices )
  {
    unsigned int sizeOfPackedChars = 0;

    const typename T_indices::size_type length = indices.size();
    unsigned int sizeOfPackedArrayChars = length*sizeof(T);
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
  static unsigned int PrivateUnpackArray( const char*& buffer, Array1dT<T>& array, const T_indices& indices )
  {

    unsigned int sizeOfUnpackedChars = 0;

    typename T_indices::size_type array_length;

    sizeOfUnpackedChars += Unpack( buffer, array_length );

    if( array_length != indices.size() )
    {
      SLIC_ERROR("bufvector::PrivateUnpackArray(): incorrect number of data");
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
  unsigned int PrivatePackSet( const std::set<T>& var )
  {

    unsigned int sizeOfPackedChars = 0;

    const typename std::set<T>::size_type length = var.size();

    sizeOfPackedChars += Pack( length );

    for( typename std::set<T>::const_iterator i=var.begin() ; i!=var.end() ; ++i )
    {
      sizeOfPackedChars += this->Pack(*i);
    }

    return sizeOfPackedChars;
  }



  template< typename T>
  static unsigned int PrivateUnpackSet( const char*& buffer, std::set<T>& setToRead )
  {
    setToRead.clear();

    unsigned int sizeOfUnpackedChars = 0;

    typename std::set<T>::size_type set_length;
    sizeOfUnpackedChars += Unpack( buffer, set_length );


    for( typename std::set<T>::size_type a=0 ; a<set_length ; ++a )
    {
      T temp;
      sizeOfUnpackedChars += Unpack( buffer, temp );
      setToRead.insert( temp );
    }

    return sizeOfUnpackedChars;
  }





  //********************************************************************************************************************
  template< typename T_KEY, typename T_VAL >
  unsigned int PrivatePackMap( const std::map<T_KEY,T_VAL>& var )
  {

    unsigned int sizeOfPackedChars = 0;

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
  static unsigned int PrivateUnpackMap( const char*& buffer, std::map<T_KEY,T_VAL>& map )
  {
    map.clear();

    unsigned int sizeOfUnpackedChars = 0;

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
  static unsigned int PrivatePack( char*& buffer, const T& var )
  {
    memcpy( buffer, &var, sizeof(T) );
    buffer += sizeof(T);
    return sizeof(T);
  }








  template< typename T>
  unsigned int PrivatePackGlobal( const T& container, const gArray1d& localToGlobal );



  template< typename T, typename T_indices >
  unsigned int PrivatePackGlobal( const T& container, const T_indices& indices, const gArray1d& localToGlobal );






  template< typename T, typename T_indices >
  struct PrivatePackRelationT
  {
    static unsigned int Pack( bufvector& buffer, const T& relation, const T_indices& indices, const bool packGlobal );
  };
  template< typename T_indices >
  struct PrivatePackRelationT<OneToOneRelation,T_indices>
  {
    static unsigned int Pack( bufvector& buffer, const OneToOneRelation& relation, const T_indices& indices, const bool packGlobal );
  };
  template< typename T_indices >
  struct PrivatePackRelationT<FixedOneToManyRelation,T_indices>
  {
    static unsigned int Pack( bufvector& buffer, const FixedOneToManyRelation& relation, const T_indices& indices, const bool packGlobal );
  };

  template< typename T, typename T_indices >
  unsigned int PrivatePackRelation( const T& relation, const T_indices& indices, const bool packGlobal )
  {
    return PrivatePackRelationT<T,T_indices>::Pack(*this,relation,indices,packGlobal);
  }



  template< typename T >
  static unsigned int PrivateUnpackRelation( const char*& buffer, T& relation, const lArray1d& indices, const bool unpackGlobal );



};


inline unsigned int bufvector::Pack( const sArray1d& container )
{
  unsigned int sizeOfPackedChars = 0;
  const localIndex arrayLength = container.size();

  sizeOfPackedChars += this->Pack( arrayLength );
  for( sArray1d::const_iterator i=container.begin() ; i!=container.end() ; ++i )
  {
    sizeOfPackedChars += this->Pack(*i);
  }

  return sizeOfPackedChars;
}



inline unsigned int bufvector::Unpack( const char*& buffer, sArray1d& array )
{
  array.clear();
  unsigned int sizeOfUnpackedChars = 0;

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
unsigned int bufvector::PrivatePackGlobal( const T& container, const gArray1d& localToGlobal )
{
  const typename T::size_type length = container.size();
  unsigned int sizeOfPackedChars = 0;

//  std::cout<<"container.size() = "<<length<<std::endl;

  sizeOfPackedChars += this->Pack(length);
  for( typename T::const_iterator i=container.begin() ; i!=container.end() ; ++i )
  {
    const globalIndex temp = localToGlobal[*i];
    sizeOfPackedChars += this->Pack( temp );
  }
  return sizeOfPackedChars;
}


inline unsigned int bufvector::UnpackGlobal( const char*& buffer, const std::map<globalIndex,localIndex>& globalToLocal, lArray1d& array )
{
  array.clear();
  unsigned int sizeOfUnpackedChars = 0;


  localIndex array_length;
  sizeOfUnpackedChars += Unpack( buffer, array_length );

  array.resize(array_length);
  gArray1d temp(array_length);
  unsigned int length = array_length * sizeof(localIndex);

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

inline unsigned int bufvector::UnpackGlobal( const char*& buffer, const std::map<globalIndex,localIndex>& globalToLocal, lSet& array )
{
  array.clear();
  unsigned int sizeOfUnpackedChars = 0;


  localIndex array_length;
  sizeOfUnpackedChars += Unpack( buffer, array_length );

//  std::cout<<"array_length = "<<array_length<<std::endl;

  gArray1d temp(array_length);
  unsigned int length = array_length * sizeof(localIndex);

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
inline unsigned int bufvector::PrivatePackRelationT<T,T_indices>::Pack( bufvector& buffer, const T& relation, const T_indices& indices, const bool packGlobal )
{
  unsigned int sizeOfPackedChars = 0;

  if( packGlobal )
  {
    const gArray1d& localToGlobal = relation.RelatedObjectLocalToGlobal();
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
inline unsigned int bufvector::PrivatePackRelationT<OneToOneRelation,T_indices>::Pack( bufvector& buffer, const OneToOneRelation& relation, const T_indices& indices, const bool packGlobal )
{
  unsigned int sizeOfPackedChars = 0;

  if( packGlobal )
  {
    const gArray1d& localToGlobal = relation.RelatedObjectLocalToGlobal();
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
inline unsigned int bufvector::PrivatePackRelationT<FixedOneToManyRelation,T_indices>::Pack( bufvector& buffer, const FixedOneToManyRelation& relation, const T_indices& indices, const bool packGlobal )
{
  unsigned int sizeOfPackedChars = 0;

  sizeOfPackedChars += buffer.Pack( static_cast<localIndex>(relation.Dimension(1)) );

  if( packGlobal )
  {
    const gArray1d& localToGlobal = relation.RelatedObjectLocalToGlobal();
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
inline unsigned int bufvector::PrivateUnpackRelation( const char*& buffer, T& relation, const lArray1d& indices, const bool unpackGlobal )
{
  unsigned int sizeOfUnpackedChars = 0;

  if( unpackGlobal )
  {
    const std::map<globalIndex,localIndex>& globalToLocal = relation.RelatedObjectGlobalToLocal();
    for( lArray1d::const_iterator i = indices.begin() ; i != indices.end() ; ++i )
    {
      sizeOfUnpackedChars += bufvector::UnpackGlobal( buffer, globalToLocal, relation[*i] );
    }
  }
  else
  {
    for( lArray1d::const_iterator i = indices.begin() ; i != indices.end() ; ++i )
    {
      sizeOfUnpackedChars += bufvector::Unpack( buffer, relation[*i] );
    }
  }

  return sizeOfUnpackedChars;
}

template<>
inline unsigned int bufvector::PrivateUnpackRelation( const char*& buffer, OneToOneRelation& relation, const lArray1d& indices, const bool unpackGlobal )
{
  unsigned int sizeOfUnpackedChars = 0;

  if( unpackGlobal )
  {
    const std::map<globalIndex,localIndex>& globalToLocal = relation.RelatedObjectGlobalToLocal();
    for( lArray1d::const_iterator i = indices.begin() ; i != indices.end() ; ++i )
    {
      globalIndex temp;
      sizeOfUnpackedChars += bufvector::Unpack( buffer, temp );
      relation[*i] = stlMapLookup(globalToLocal,temp);
    }
  }
  else
  {
    for( lArray1d::const_iterator i = indices.begin() ; i != indices.end() ; ++i )
    {
      sizeOfUnpackedChars += bufvector::Unpack( buffer, relation[*i] );
    }
  }

  return sizeOfUnpackedChars;
}

template<>
inline unsigned int bufvector::PrivateUnpackRelation( const char*& buffer, FixedOneToManyRelation& relation, const lArray1d& indices, const bool unpackGlobal )
{
  unsigned int sizeOfUnpackedChars = 0;

  localIndex dimension;
  sizeOfUnpackedChars += bufvector::Unpack( buffer, dimension );

  if( dimension != static_cast<int>(relation.Dimension(1)) )
    SLIC_ERROR("bufvector::PrivateUnpackRelation(): mismatched dimension");

  if( unpackGlobal )
  {
    const std::map<globalIndex,localIndex>& globalToLocal = relation.RelatedObjectGlobalToLocal();
    for( lArray1d::const_iterator i = indices.begin() ; i != indices.end() ; ++i )
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
    for( lArray1d::const_iterator i = indices.begin() ; i != indices.end() ; ++i )
    {
      for( localIndex j=0 ; j<dimension ; ++j )
      {
        sizeOfUnpackedChars += bufvector::Unpack( buffer, relation[*i][j] );
      }
    }
  }

  return sizeOfUnpackedChars;
}


inline unsigned int bufvector::Unpack( const char*& buffer, OneToOneRelation& relation, const lArray1d& indices, const bool unpackGlobal ) { return PrivateUnpackRelation( buffer, relation, indices, unpackGlobal ); }
inline unsigned int bufvector::Unpack( const char*& buffer, FixedOneToManyRelation& relation, const lArray1d& indices, const bool unpackGlobal ) { return PrivateUnpackRelation( buffer, relation, indices, unpackGlobal ); }
inline unsigned int bufvector::Unpack( const char*& buffer, OrderedVariableOneToManyRelation& relation, const lArray1d& indices, const bool unpackGlobal ) { return PrivateUnpackRelation( buffer, relation, indices, unpackGlobal ); }
inline unsigned int bufvector::Unpack( const char*& buffer, UnorderedVariableOneToManyRelation& relation, const lArray1d& indices, const bool unpackGlobal ) { return PrivateUnpackRelation( buffer, relation, indices, unpackGlobal ); }
//inline unsigned int bufvector::Unpack( const char*& buffer, OrderedVariableOneToManyPairRelation& relation, const lArray1d& indices, const bool unpackGlobal ) { return PrivateUnpackRelation( buffer, relation, indices, unpackGlobal ); }
//inline unsigned int bufvector::Unpack( const char*& buffer, UnorderedVariableOneToManyPairRelation& relation, const lArray1d& indices, const bool unpackGlobal ) { return PrivateUnpackRelation( buffer, relation, indices, unpackGlobal ); }





//**********************************************************************************************************************
inline unsigned int bufvector::Pack( const std::string& var )
{
  unsigned int sizeOfPackedChars = var.size();

  this->PrivatePack( sizeOfPackedChars );

  const char* cvar = var.data();
  for( unsigned int i=0 ; i<sizeOfPackedChars ; ++i )
  {
    this->push_back(*(cvar++));
  }

  sizeOfPackedChars += sizeof( unsigned int );
  return sizeOfPackedChars;
}

inline unsigned int bufvector::Pack(char*& buffer,  const std::string& var )
{
  unsigned int sizeOfPackedChars = var.size();

  PrivatePack(buffer,  sizeOfPackedChars );

  for( unsigned int i=0 ; i<sizeOfPackedChars ; ++i )
  {
	*buffer = var[i];
	buffer++;
  }

  sizeOfPackedChars += sizeof( unsigned int );
  return sizeOfPackedChars;
}

inline unsigned int bufvector::Unpack( const char*& buffer, std::string& var )
{
  unsigned int sizeOfUnpackedChars = 0;
  unsigned int stringsize = 0;

  sizeOfUnpackedChars += PrivateUnpack( buffer, stringsize );

  var.assign( buffer, stringsize );
  buffer += stringsize;
  sizeOfUnpackedChars += stringsize;

  return sizeOfUnpackedChars;
}















#endif /* BUFVECTOR_H_ */
