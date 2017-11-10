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
 * File: FileManagerT.h
 * Class provides file IO
 * created : RRS (10/11/2001)
 */
 
#ifndef BIN_STREAM_H
#define BIN_STREAM_H

// ***** Included Headers *****************************************************
#include "common/DataTypes.hpp"
#include <map>
#include <set>
#include <fstream>
#include <iostream>

#ifdef USE_ATK
#include "slic/slic.hpp"
#endif

namespace geosx
{
// ****************************************************************************
// ***** BINSTREAM CLASS DECLARATION ******************************************
// ****************************************************************************
class BinStream
{
public:
  BinStream(void);
  virtual ~BinStream(void);

  virtual void open( const char* filename, const bool truncate ) = 0;
  virtual void close(void) = 0;

protected:

};




// ****************************************************************************
// ***** oBINSTREAM CLASS DECLARATION *****************************************
// ****************************************************************************
class oBinStream : public BinStream
{
public:
  oBinStream(void);
  ~oBinStream(void);

  void open( const char* filename, const bool truncate = false );
  void close(void);

  template<class TYPE>
  void write( const TYPE* const p_var , const int var_length )
  {
    output.write( reinterpret_cast<const char*>(p_var) , sizeof(TYPE)*var_length );
  }
  
  template< typename TYPE >
  void write( const TYPE& val )
  {
    this->write( &val , 1 );
  }

  void write( const std::string& str )
  {
    std::string::size_type length = str.size();
    this->write( length );
    this->write( str.c_str(), str.size() );
  }

  template< typename TYPE >
  void write( const array<TYPE>& arr )
  {
    const typename array<TYPE>::size_type length = arr.size();
    this->write( length );
    this->write( arr.data(), arr.size() );
  }

  template< typename TYPE >
  void write( const std::set<TYPE>& set )
  {
    const typename std::set<TYPE>::size_type length = set.size();
    this->write( length );
    for( typename std::set<TYPE>::const_iterator i=set.begin() ; i!=set.end() ; ++i )
    {
      this->write( *i );
    }
  }


  template< typename TYPE >
  void write( const Array2dT<TYPE>& arr )
  {
    const typename Array2dT<TYPE>::size_type length = arr.size();
    this->write( length );

    const typename Array2dT<TYPE>::size_type dimension[2] = { arr.Dimension(0), arr.Dimension(1) };
    this->write( dimension, 2 );
    this->write( arr.data(), arr.size() );
  }


  template< typename TYPE >
  void write( const array<array<TYPE> >& arr )
  {
    const typename array<array<TYPE> >::size_type length0 = arr.size();
    this->write( length0 );

    for( typename array<array<TYPE> >::const_iterator i=arr.begin() ; i!=arr.end() ; ++i )
    {
      this->write(*i);
    }
  }


  template< typename TYPE >
  void write( const array<std::set<TYPE> >& arr )
  {

    const typename array<std::set<TYPE> >::size_type length0 = arr.size();
    this->write( reinterpret_cast<const char*>(&length0), sizeof(typename array<std::set<TYPE> >::size_type) );

    for( typename array<std::set<TYPE> >::const_iterator i=arr.begin() ; i!=arr.end() ; ++i )
    {
      this->write( *i );
    }
  }


  template< typename T1, typename T2 >
  void write( const std::map< T1, T2 >& datamap )
  {
    this->write( datamap.size() );
    for( typename std::map<T1,T2>::const_iterator i=datamap.begin() ; i!=datamap.end() ; ++i )
    {
      this->write( i->first );
      this->write( i->second );
    }
  }

  template< typename T >
  void write( const std::map< std::string, T>& member )
  {



    // first get the number of map entries that we are writing.
    typename std::map< std::string, T >::size_type numEntries = 0;
    // Iterate over all entries in the member map
    for( typename std::map< std::string, T >::const_iterator i = member.begin() ; i!=member.end() ; ++i )
    {
      // the field name is the key to the map
      const std::string fieldName = i->first;

      // check to see if the field should be written
      std::map<std::string, FieldBase*>::const_iterator fieldAttributes = FieldInfo::AttributesByName.find(fieldName);

      if( fieldAttributes == FieldInfo::AttributesByName.end() )
      {
        ++numEntries;
      }
      else if( fieldAttributes->second->m_WriteToRestart )
      {
        ++numEntries;
      }
    }
    // write the number of entries that we are writing out.
    this->write(numEntries);



    // iterate over all entries in the member map
    for( typename std::map< std::string, T >::const_iterator i = member.begin() ; i!=member.end() ; ++i )
    {
      // the field name is the key to the map
      const std::string fieldName = i->first;

      // check to see if the field should be written
      std::map<std::string, FieldBase*>::const_iterator fieldAttributes = FieldInfo::AttributesByName.find(fieldName);

      bool writeField = false;
      if( fieldAttributes == FieldInfo::AttributesByName.end() )
      {
        writeField = true;
      }
      else if( fieldAttributes->second->m_WriteToRestart )
      {
        writeField = true;
      }

      if( writeField )
      {
//        std::cout<<"    writing "<<fieldName<<std::endl;

        // the field data is mapped value
        const T& fieldData = i->second;

        // write the field name
        this->write( fieldName );

        // write the field data
        this->write( fieldData );
      }
    }
  }

  //***** Data Member Declarations **********************************************
private:
  std::ofstream output;

public:
  std::ofstream::pos_type tellp(void)
  { return output.tellp(); }

};








// ****************************************************************************
// ***** iBINSTREAM CLASS DECLARATION *****************************************
// ****************************************************************************
class iBinStream : public BinStream
{
public:
  iBinStream(void);
  ~iBinStream(void);

  void open( const char* filename, const bool truncate = false );
  void close(void);

  template<class TYPE>
  void read( TYPE* const p_var , const int var_length )
  {
    input.read( reinterpret_cast<char*>(p_var) , sizeof(TYPE)*var_length );
  }

  template<class TYPE>
  void read( TYPE& p_var )
  {
    this->read( &p_var , 1 );
  }

  void read( std::string& str )
  {
    std::string::size_type length;
    this->read( length );

    array<char> readstring;
    readstring.resize(length);
    this->read( readstring.data(), readstring.size() );

    str.assign( readstring.begin(),readstring.end() );

  }

  template< typename TYPE >
  void read( array<TYPE>& arr, const bool realloc = true )
  {
    const typename array<TYPE>::size_type length = arr.size();
    typename array<TYPE>::size_type readLength;
    this->read( readLength );

    if( readLength != length )
    {
     if( realloc )
     {
       arr.resize(readLength);
     }
     else
     {
#ifdef USE_ATK
      SLIC_ERROR( "BinStream::read(array<TYPE>& arr): length mismatch\n");
#endif
     }
    }

    this->read( arr.data(), arr.size() );
  }


  template< typename TYPE >
  void read( std::set<TYPE>& set, const bool realloc = true )
  {
    typename std::set<TYPE>::size_type readLength;
    this->read( readLength );

    if( readLength != set.size() && !realloc )
    {
#ifdef USE_ATK
      SLIC_ERROR( "BinStream::read(std::set<TYPE>& set): length mismatch\n");
#endif
    }

    for( typename std::set<TYPE>::size_type i=0 ; i<readLength ; ++i )
    {
      TYPE readVal;
      this->read( readVal );
      set.insert(readVal);
    }
  }

  template< typename TYPE >
  void read( Array2dT<TYPE>& arr, const bool realloc = true )
  {
    const typename Array2dT<TYPE>::size_type length = arr.size();
    typename Array2dT<TYPE>::size_type readLength;
    this->read( readLength );
    if( readLength != length ) {
#ifdef USE_ATK
      SLIC_ERROR( "BinStream::read(Array2dT<TYPE>& arr): length mismatch\n");
#endif
    }

    size_t dimension[2] = { arr.Dimension(0), arr.Dimension(1) };
    size_t readDimension[2];
    this->read( readDimension, 2 );

    if( dimension[0] != readDimension[0] || dimension[1] != readDimension[1] )
    {
      if( realloc )
      {
        arr.resize2( readDimension[0], readDimension[1] );
      }
      else
      {
#ifdef USE_ATK
        SLIC_ERROR( "BinStream::read(Array2dT<TYPE>& arr): dimension mismatch\n");
#endif
      }
    }

    this->read( arr.data(), arr.size() );
  }


  template< typename TYPE >
  void read( array<array<TYPE> >& arr, const bool realloc = true )
  {
    const typename array<array<TYPE> >::size_type length = arr.size();
    typename array<array<TYPE> >::size_type readLength;
    this->read( readLength );

    if( readLength != length )
    {
      if( realloc )
      {
        arr.resize(readLength);
      }
      else
      {
#ifdef USE_ATK
        SLIC_ERROR( "BinStream::read(array<array<TYPE> >& arr): length mismatch\n");
#endif
      }
    }

    for( typename array<array<TYPE> >::iterator i=arr.begin() ; i!=arr.end() ; ++i )
    {
      read( *i, realloc );
    }
  }


  template< typename TYPE >
  void read( array<std::set<TYPE> >& arr, const bool realloc = true )
  {
    const typename array<std::set<TYPE> >::size_type length = arr.size();
    typename array<std::set<TYPE> >::size_type readLength;
    this->read( readLength );

    if( readLength != length )
    {
      if( realloc )
      {
        arr.resize(readLength);
      }
      else
      {
#ifdef USE_ATK
        SLIC_ERROR( "BinStream::read(array<std::set<TYPE> >& arr): length mismatch\n");
#endif
      }
    }
    for( typename array<std::set<TYPE> >::iterator i=arr.begin() ; i!=arr.end() ; ++i )
    {
      this->read( *i, realloc );
    }
  }


  template< typename T1, typename T2 >
  void read( std::map< T1, T2 >& member )
  {
    typename std::map< T1, T2 >::size_type length;
    this->read( length );
    for( typename std::map< T1, T2 >::size_type i=0 ; i<length ; ++i  )
    {
      T1 key;
      T2 val;
      this->read( key );
      this->read( val );
      member[key] = val;
    }
  }


  template< typename T >
  void read( std::map< std::string, T>& member, const bool realloc = true )
  {

    typename std::map< std::string, T >::size_type numEntries = 0;
    this->read(numEntries);



    for( typename std::map< std::string, T >::size_type i=0 ; i<numEntries ; ++i )
    {
      // read the field name
      std::string readFieldName;
      this->read(readFieldName);

      std::map<std::string, FieldBase*>::iterator fieldAttributes = FieldInfo::AttributesByName.find(readFieldName);

      if( fieldAttributes == FieldInfo::AttributesByName.end()  )
      {
        if( realloc )
        {
          FieldInfo::AttributesByName[readFieldName]  = new FieldBase( FieldInfo::noKey, readFieldName, true, true );
          fieldAttributes = FieldInfo::AttributesByName.find(readFieldName);
        }
        else
        {
#ifdef USE_ATK
          SLIC_ERROR("ObjectDataStructureBaseT::ReadMapFromRestart: name not found\n");
#endif
        }
      }


//      std::cout<<"    reading "<<readFieldName<<std::endl;

      T& fieldData = member[readFieldName];

      // read the field data
      this->read( fieldData );


    }

    /*
    // iterate over all entries in the member map
    for( typename std::map< std::string, T >::iterator i = member.begin() ; i!=member.end() ; ++i )
    {
      // the field name is the key to the map
      const std::string fieldName = i->first;

      std::map<std::string, FieldBase*>::const_iterator fieldAttributes = FieldInfo::AttributesByName.find(fieldName);
      bool readField = false;
      if( fieldAttributes == FieldInfo::AttributesByName.end() )
      {
        readField = true;
      }
      else if( fieldAttributes->second->m_WriteToRestart )
      {
        readField = true;
      }

      if( readField )
      {
        std::cout<<"    reading "<<fieldName<<std::endl;
        // the field data is mapped value
        T& fieldData = i->second;

        // read the field name
        std::string readFieldName;
        this->read(readFieldName);

        if( fieldName != readFieldName )
#ifdef USE_ATK
          SLIC_ERROR("ObjectDataStructureBaseT::ReadMapFromRestart: field name mismatch\n");
#endif

        // write the field data
        this->read( fieldData );
      }
    }*/
  }

private:
  //***** Data Member Declarations **********************************************
  std::ifstream input;

public:
  std::ifstream::pos_type tellg(void)
  { return input.tellg(); }

};


}
#endif
