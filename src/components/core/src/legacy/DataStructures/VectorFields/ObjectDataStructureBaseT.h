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
 * @file ObjectDataStructureBaseT.h
 * @author Randolph Settgast
 * @date created on Sep 13, 2010
 */
 

#ifndef OBJECTDATASTUCTUREBASET_H_
#define OBJECTDATASTUCTUREBASET_H_

// keep this flag for time being - ultimately will replace with static flag in object manager
#define PACK_FIELD_NAMES_IN_MPI_BUFFER

#include "Common/Common.h"
#include "Common/typedefs.h"
#include "Common/intrinsic_typedefs.h"
#include <fstream>
#include <map>
#include <set>
#include <string>
#include "IO/silo/SiloFile.h"
//#include "BoundaryConditions/BoundaryConditions.h"
#include "Utilities/StringUtilities.h"
#include "Utilities/Utilities.h"
#include <limits.h>

//#include "DataStructures/InterObjectRelation.h"

#include "ArrayT/bufvector.h"

class oBinStream;
class iBinStream;
class BoundaryConditionBase;
class PhysicalDomainT;
class EncapsulatedObjectManagerBase;


// *********************************************************************************************************************
// *********************************************************************************************************************
/**
 * @author Randolph Settgast
 * @brief base class that manages the object data
 * The ObjectDataStructureBaseT class manages the object data by providing
 * members that are std::maps<key,array<TYPE> >, where the key is a std::string.
 * Because std::map is used there can only be unique combinations of key and TYPE.
 * It is possible to have identical keys with different TYPES, so be aware of this
 * when asking for data from the manager.
 *
 */
class ObjectDataStructureBaseT
{
public:


  //********************************************************************************************************************
  // Basic class functionality
  //********************************************************************************************************************

  /// this is to specify the type of derived object that can be used for checking base objects passed into a
  /// function
  enum ObjectType
  {
    ElementManager,
    ElementRegion,
    NodeManager,
    FaceManager,
    ExternalFaceManager,
    FaultElementManager,
    EdgeManager,
    DiscreteElementManager,
    EllipsoidalDiscreteElementManager,
    EllipsoidalContactManager,
    ContactBaseManager,
    CartesianGridManager,
    CrackFaceManager,
    CrackSurfaceManager,
    CrackSurfaceVertexManager,
    Temp,
    WellboreManager,
    MicroseismicManager
  };


  /// Constructor that sets object type
  ObjectDataStructureBaseT( const ObjectType objectType );

  /// Copy constructor
  ObjectDataStructureBaseT( const ObjectDataStructureBaseT& init );

  /// Destructor
  virtual ~ObjectDataStructureBaseT();

  /// to initialize the data object
  virtual void Initialize(  ) = 0 ;

  //********************************************************************************************************************
  // Allocation, resize, insertion, etc
  //********************************************************************************************************************

  void insert( const localIndex i, const globalIndex gIndex);


  virtual void erase( const localIndex i );

  /// function to allocate all data arrays to given size.
  virtual globalIndex resize( const localIndex size, const bool assignGlobals = false ) ;

protected:
  virtual globalIndex insert( const localIndex i, const bool assignGlobals = false );

  template< typename T >
  void AllocateTemporaryFields( const array<string>& names, array<array<T>* >& vars );

  template< typename T >
  void DeallocateTemporaryFields( const array<string>& names );

private:
  virtual void insertBase( const localIndex i );

  /// function to resize a single data array of the data structure
  template< typename T >
  void ResizeObjectField( std::map< std::string, T>& member );

  /// templated insert function for fields
  template< typename T >
  void insert( std::map< std::string, array<T> >& member, const localIndex i );

  /// insert in map to hold fixed size one-to-many maps
  void insert( std::map< std::string, FixedOneToManyRelation >& member, const localIndex i);

  template< typename T >
  void insert( std::map< std::string, InterObjectRelation<array<T> > >& member, const localIndex i );


  template< typename T >
  void erase( std::map< std::string, array<T> >& member, const localIndex i );

  template< typename T >
  void erase( std::map< std::string, InterObjectRelation<T> >& member, const localIndex i );

  /// insert in map to hold fixed size one-to-many maps
  template< typename T >
  void erase( std::map< std::string, Array2dT<T> >& member, const localIndex i);

public:




  //********************************************************************************************************************
  // registration, field/relationship additions
  //********************************************************************************************************************

  /// add a data field to a member
  template< FieldKey FIELDKEY >
  int AddKeyedDataField();

  /// add a data field to a member
  template< typename T >
  int AddKeylessDataField( const std::string& name, const bool restart = false, const bool plot = false );
  
  int AddKeylessDataField( FieldType type, const std::string& name, const bool restart = false, const bool plot = false );

  template< typename T >
  void AddMap( const std::string& name );

  template< typename T >
  void RemoveDataField( const std::string& name );


  //********************************************************************************************************************
  // access functions
  //********************************************************************************************************************

  /// returns reference to specified field
  template< FieldKey FIELDKEY>
  array<typename Field<FIELDKEY>::Type>& GetFieldData( )
  { return const_cast<array<typename Field<FIELDKEY>::Type>&>( static_cast<const ObjectDataStructureBaseT&>(*this).GetFieldData<FIELDKEY>()); }


  /// returns const reference to specified field
  template< FieldKey FIELDKEY>
  const array<typename Field<FIELDKEY>::Type>& GetFieldData( ) const;


  /// returns reference to specified field
  template< FieldKey FIELDKEY>
  array<typename Field<FIELDKEY>::Type>* GetFieldDataPointer( )
  { return const_cast<array<typename Field<FIELDKEY>::Type>*>( static_cast<const ObjectDataStructureBaseT&>(*this).GetFieldDataPointer<FIELDKEY>()); }


  /// returns const reference to specified field
  template< FieldKey FIELDKEY>
  const array<typename Field<FIELDKEY>::Type>* GetFieldDataPointer( ) const;


  /// returns reference to specified field
  template< typename TYPE >
  array<TYPE>& GetFieldData( const std::string& fieldName )
  { return const_cast<array<TYPE>&>( static_cast<const ObjectDataStructureBaseT&>(*this).GetFieldData<TYPE>(fieldName)); }

  /// returns const reference to specified field
  template< typename TYPE >
  const array<TYPE>& GetFieldData( const std::string& name ) const;

  /// returns reference to specified field
  template< typename TYPE >
  array<TYPE>* GetFieldDataPointer( const std::string& fieldName )
  { return const_cast<array<TYPE>*>( static_cast<const ObjectDataStructureBaseT&>(*this).GetFieldDataPointer<TYPE>(fieldName)); }

  /// returns const reference to specified field
  template< typename TYPE >
  const array<TYPE>* GetFieldDataPointer( const std::string& name ) const;

  
  /// returns reference to specified set
  lSet& GetSet( const std::string& name );

  /// returns const reference to specified set
  const lSet& GetSet( const std::string& name ) const;
  
  template<typename TYPE>
  void SetSubsetValue(const std::string& fieldName, const std::string& subsetname, const TYPE& value);

  /// returns true if field exists
  template< typename TYPE >
  bool HasField( const std::string& name ) const;

  /// returns the type of a named field 
  /// (nb. does not return multiple types if fields of different type have the same name)
  FieldType GetFieldType(const std::string& fieldName) const{
    FieldType rv = FieldInfo::numFieldTypes;
    //if( m_LocalIndexData.count(fieldName) == 1 ){ rv = FieldInfo::longField;}
    //else
    if( m_IntegerData.count(fieldName) == 1 ){ rv = FieldInfo::integerField;}
    else if( m_LocalIndexData.count(fieldName) == 1 ){ rv = FieldInfo::localIndexField; }
    else if( m_GlobalIndexData.count(fieldName) == 1 ){ rv = FieldInfo::globalIndexField; }
    else if( m_realData.count(fieldName) == 1 ){ rv = FieldInfo::realField;}
    else if( m_R1TensorData.count(fieldName) == 1 ){ rv = FieldInfo::R1TensorField;}
    else if( m_R2TensorData.count(fieldName) == 1 ){ rv = FieldInfo::R2TensorField;}
    else if( m_R2SymTensorData.count(fieldName) == 1 ){ rv = FieldInfo::R2SymTensorField;}

    return rv;
  };

  /// Checks if the field has only one type on the same object
  bool FieldHasSingleType(const std::string& fieldName) const{
    int count = 0;
    bool rv = false;

    if( m_LocalIndexData.count(fieldName) == 1 ){ ++count; }
    if( m_GlobalIndexData.count(fieldName) == 1 ){ ++count; }
    if( m_IntegerData.count(fieldName) == 1 ){ ++count; }
    if( m_realData.count(fieldName) == 1 ){ ++count;}
    if( m_R1TensorData.count(fieldName) == 1 ){ ++count;}
    if( m_R2TensorData.count(fieldName) == 1 ){ ++count;}
    if( m_R2SymTensorData.count(fieldName) == 1 ){ ++count;}
    if( m_ArrayR1TensorData.count(fieldName) == 1 ){ ++count;}

    if(count == 1) rv = true;

    return rv;
  };
  

  /**
   * @author Randolph Settgast
   * @tparam T the type of data in the map
   * @return a reference to the map
   *
   * This templated function returns a member map based on type. This is achieved through
   * explicit template specialization.
   */
  template< typename T >
  std::map< std::string, T >& GetDataMemberMap();



  /**
   * @author Randolph Settgast
   * @tparam T the type of data in the map
   * @return a reference to the map
   *
   * This templated function returns a const member map based on type. This is achieved through
   * explicit template specialization.
   */
  template< typename T >
  const std::map< std::string, T >& GetDataMemberMap() const;

  // sets the value of a const pointer to a keyed field location
  template< FieldKey FIELDKEY >
  void SetConstPointer( array<typename Field<FIELDKEY>::Type>* const& pointer );



  //********************************************************************************************************************
  // IO
  //********************************************************************************************************************

  void WriteSilo( SiloFile& siloFile,
                  const std::string& siloDirName,
                  const std::string& meshname,
                  const int centering,
                  const int cycleNum,
                  const realT problemTime,
                  const bool isRestart,
                  const std::string& regionName = "none",
                  const lArray1d& mask = lArray1d() );

  int ReadSilo( const SiloFile& siloFile,
                const std::string& siloDirName,
                const std::string& meshname,
                const int centering,
                const int cycleNum,
                const realT problemTime,
                const bool isRestart,
                const std::string& regionName = "none",
                const lArray1d& mask = lArray1d() );

  void UncheckAllFieldsForPlot();


protected:
  template< typename OUTPUTTYPE, typename T >
  void UncheckAllFieldsForPlot( const std::map< std::string, T>& member);

  void WriteSilo( SiloFile& siloFile,
                  const std::string& meshname,
                  const int centering,
                  const int cycleNum,
                  const realT problemTime,
                  const bool isRestart,
                  const std::string& multiRoot,
                  const std::string& regionName = "none",
                  const lArray1d& mask = lArray1d() ) const;


  void ReadSilo( const SiloFile& siloFile,
                 const std::string& meshname,
                 const int centering,
                 const int cycleNum,
                 const realT problemTime,
                 const bool isRestart,
                 const std::string& regionName = "none",
                 const lArray1d& mask = lArray1d() );


  void EncapsulatedObjectsToManagedData_PreRead( const EncapsulatedObjectManagerBase& eom );
  void EncapsulatedObjectsToManagedData_PostRead( const EncapsulatedObjectManagerBase& eom );

  void ManagedDataToEncapsulatedObjects_PreRead( const EncapsulatedObjectManagerBase& eom );

  void ManagedDataToEncapsulatedObjects_PostRead( EncapsulatedObjectManagerBase& eom );

public:
  template< typename T >
  void AllocateDummyFields( const array<string>& names, const bool plotFlag  )
  {
    for( array<string>::size_type i=0 ; i<names.size() ; ++i )
    {
      this->AddKeylessDataField<T>(names[i],true,plotFlag);
    }
  }


  template< typename T >
  void AllocateDummyFields( const array<string>& names, array<array<T>* >& vars, const bool plotFlag  )
  {
    vars.resize( names.size() );
    for( array<string>::size_type i=0 ; i<names.size() ; ++i )
    {
      this->AddKeylessDataField<T>(names[i],true,plotFlag);
      vars[i] = &this->GetFieldData<T>(names[i]);
    }
  }

  template< typename T >
  void SetDummyFieldPointers( const array<string>& names, array<array<T>* >& vars )
  {
    vars.resize( names.size() );
    for( array<string>::size_type i=0 ; i<names.size() ; ++i )
    {
      vars[i] = &this->GetFieldData<T>(names[i]);
    }
  }

  template< typename T >
  void DeallocateDummyFields( const array<string>& names )
  {
    for( array<string>::size_type i=0 ; i<names.size() ; ++i )
      this->RemoveDataField<T>(names[i]);
  }


protected:


  virtual void WriteNonManagedDataMembersToSilo( SiloFile& ,
                                                 const std::string& ,
                                                 const std::string& ,
                                                 const int ,
                                                 const int ,
                                                 const realT ,
                                                 const bool ,
                                                 const std::string& ,
                                                 const std::string& regionName = "none",
                                                 const lArray1d& mask = lArray1d() ) {}

  virtual void ReadNonManagedDataMembersFromSilo( const SiloFile& ,
                                                  const std::string& ,
                                                  const std::string& ,
                                                  const int ,
                                                  const int ,
                                                  const realT ,
                                                  const bool ,
                                                  const std::string& regionName = "none",
                                                  const lArray1d& mask = lArray1d() ){}

  virtual void WriteSiloSetup() {}
  virtual void WriteSiloCleanup() {}

  virtual void ReadSiloSetup() {}
  virtual void ReadSiloCleanup() {}



public:





  /// reads field from file
  template< typename T >
  void ReadAsciiFieldData(const std::string& fieldName, const std::string& fileName );
  
  void ReadAsciiFieldData( FieldType type, const std::string& fieldName, const std::string& fileName ){
  	switch( type ){
      case FieldInfo::integerField:     ReadAsciiFieldData<int>( fieldName, fileName);  break;
      case FieldInfo::localIndexField:  ReadAsciiFieldData<localIndex>( fieldName, fileName);  break;
      case FieldInfo::globalIndexField: ReadAsciiFieldData<globalIndex>( fieldName, fileName);  break;
      case FieldInfo::realField:        ReadAsciiFieldData<realT>( fieldName, fileName);       break;
      case FieldInfo::R1TensorField:    ReadAsciiFieldData<R1Tensor>( fieldName, fileName);    break;
      case FieldInfo::R2TensorField:    ReadAsciiFieldData<R2Tensor>( fieldName, fileName);    break;
      case FieldInfo::R2SymTensorField: ReadAsciiFieldData<R2SymTensor>( fieldName, fileName); break; 
      case FieldInfo::integerParameter:
      case FieldInfo::realParameter:
      case FieldInfo::numFieldTypes:
      default:
        throw GPException("ReadAsciiFieldData: Unrecognized field type "+ type);
    }
  };
  
  /// reads subset of field from file
  template< typename T >
  void ReadAsciiFieldData(const std::string& fieldName, const std::string& fileName, const lSet& subset);
  
  void ReadAsciiFieldData( FieldType type, const std::string& fieldName, const std::string& fileName, const lSet& subset){
  	switch( type ){
      case FieldInfo::integerField:     ReadAsciiFieldData<int>( fieldName, fileName, subset);  break;
      case FieldInfo::localIndexField:  ReadAsciiFieldData<localIndex>( fieldName, fileName, subset);  break;
      case FieldInfo::globalIndexField: ReadAsciiFieldData<globalIndex>( fieldName, fileName, subset);  break;
      case FieldInfo::realField:        ReadAsciiFieldData<realT>( fieldName, fileName, subset);       break;
      case FieldInfo::R1TensorField:    ReadAsciiFieldData<R1Tensor>( fieldName, fileName, subset);    break;
      case FieldInfo::R2TensorField:    ReadAsciiFieldData<R2Tensor>( fieldName, fileName, subset);    break;
      case FieldInfo::R2SymTensorField: ReadAsciiFieldData<R2SymTensor>( fieldName, fileName, subset); break; 
      case FieldInfo::integerParameter:
      case FieldInfo::realParameter:
      case FieldInfo::numFieldTypes:
      default:
        throw GPException("ReadAsciiFieldData: Unrecognized field type "+ type);
    }
  };

  /* reads field from file where data is referenced by global index
   eg.
   0 1.232
   10002 3.14
   20003 4.51
  */
  template< typename T >
  void ReadIndexedAsciiFieldData(const std::string& fieldName, const std::string& fileName );

  void ReadIndexedAsciiFieldData( FieldType type, const std::string& fieldName, const std::string& fileName ){
  	switch( type ){
      case FieldInfo::integerField:     ReadIndexedAsciiFieldData<int>( fieldName, fileName);  break;
      case FieldInfo::localIndexField:  ReadIndexedAsciiFieldData<localIndex>( fieldName, fileName);  break;
      case FieldInfo::globalIndexField: ReadIndexedAsciiFieldData<globalIndex>( fieldName, fileName);  break;
      case FieldInfo::realField:        ReadIndexedAsciiFieldData<realT>( fieldName, fileName);       break;
      case FieldInfo::R1TensorField:    ReadIndexedAsciiFieldData<R1Tensor>( fieldName, fileName);    break;
      case FieldInfo::R2TensorField:    ReadIndexedAsciiFieldData<R2Tensor>( fieldName, fileName);    break;
      case FieldInfo::R2SymTensorField: ReadIndexedAsciiFieldData<R2SymTensor>( fieldName, fileName); break;
      case FieldInfo::integerParameter:
      case FieldInfo::realParameter:
      case FieldInfo::numFieldTypes:
      default:
        throw GPException("ReadAsciiFieldData: Unrecognized field type "+ type);
    }
  };

/// write a single field to a file
  template< typename T >
  void WriteAsciiFieldData(const std::string& fieldName, const std::string& fileName, bool append = false);
  
  void WriteAsciiFieldData( FieldType type, const std::string& fieldName, const std::string& fileName, bool append ){
  	switch( type ){
      case FieldInfo::integerField:     WriteAsciiFieldData<int>( fieldName, fileName, append);  break;
      case FieldInfo::localIndexField:  WriteAsciiFieldData<localIndex>( fieldName, fileName, append);  break;
      case FieldInfo::globalIndexField: WriteAsciiFieldData<globalIndex>( fieldName, fileName, append);  break;
      case FieldInfo::realField:        WriteAsciiFieldData<realT>( fieldName, fileName, append);    break;
      case FieldInfo::R1TensorField:    WriteAsciiFieldData<R1Tensor>( fieldName, fileName, append);    break;
      case FieldInfo::R2TensorField:    WriteAsciiFieldData<R2Tensor>( fieldName, fileName, append);    break;
      case FieldInfo::R2SymTensorField: WriteAsciiFieldData<R2SymTensor>( fieldName, fileName, append); break; 
      case FieldInfo::integerParameter:
      case FieldInfo::realParameter:
      case FieldInfo::numFieldTypes:
      default:
        throw GPException("WriteAsciiFieldData: Unrecognized field type "+ type);
    }
  };
  
  /// reads subset of field from file
  template< typename T >
  void WriteAsciiFieldData(const std::string& fieldName, const std::string& fileName, const lSet& subset,bool append);
  
  void WriteAsciiFieldData( FieldType type, const std::string& fieldName, const std::string& fileName, const lSet& subset,bool append = false)
  {
    switch( type ){
      case FieldInfo::integerField:     WriteAsciiFieldData<int>( fieldName, fileName, subset, append);  break;
      case FieldInfo::localIndexField:  WriteAsciiFieldData<localIndex>( fieldName, fileName, subset, append);  break;
      case FieldInfo::globalIndexField: WriteAsciiFieldData<globalIndex>( fieldName, fileName, subset, append); break;
      case FieldInfo::realField:        WriteAsciiFieldData<realT>( fieldName, fileName, subset, append);       break;
      case FieldInfo::R1TensorField:    WriteAsciiFieldData<R1Tensor>( fieldName, fileName, subset, append);    break;
      case FieldInfo::R2TensorField:    WriteAsciiFieldData<R2Tensor>( fieldName, fileName, subset, append);    break;
      case FieldInfo::R2SymTensorField: WriteAsciiFieldData<R2SymTensor>( fieldName, fileName, subset, append); break; 
      case FieldInfo::integerParameter:
      case FieldInfo::realParameter:
      case FieldInfo::numFieldTypes:
      default:
        throw GPException("WriteAsciiFieldData: Unrecognized field type "+ type);
    }
  };


  void WriteAsciiFieldData( const std::vector<FieldType>& types, const array<string>& fieldNames, const std::string& fileName, bool append);

  void WriteAsciiFieldData( const std::vector<FieldType>& types, const array<string>& fieldNames, const std::string& fileName, const lSet& subset, bool append );
  
  /// sets field to constant value 
  template< typename T >
  void SetFieldToConstantFromString(const std::string& fieldName, const std::string& value, const bool& additive);
  
  void SetFieldToConstantFromString( FieldType type, const std::string& fieldName, const std::string& value, const bool& additive ){
  	switch(type){
      case FieldInfo::integerField:     SetFieldToConstantFromString<int>( fieldName, value, additive);  break;
      case FieldInfo::localIndexField:  SetFieldToConstantFromString<localIndex>( fieldName, value, additive);       break;
      case FieldInfo::globalIndexField: SetFieldToConstantFromString<globalIndex>( fieldName, value, additive);       break;
      case FieldInfo::realField:        SetFieldToConstantFromString<realT>( fieldName, value, additive);       break;
      case FieldInfo::R1TensorField:    SetFieldToConstantFromString<R1Tensor>( fieldName, value, additive);    break;
      case FieldInfo::R2TensorField:    SetFieldToConstantFromString<R2Tensor>( fieldName, value, additive);    break;
      case FieldInfo::R2SymTensorField: SetFieldToConstantFromString<R2SymTensor>( fieldName, value, additive); break;
      case FieldInfo::integerParameter:
      case FieldInfo::realParameter:
      case FieldInfo::numFieldTypes:
      default:
        throw GPException("SetConstantFieldFromString: Unrecognized field type "+ type);
    }
  };
  
  /// sets subset of field to constant value 
  template< typename T >
  void SetFieldToConstantFromString(const std::string& fieldName, const std::string& value, const lSet& subset, const bool& additive);
  
  void SetFieldToConstantFromString( FieldType type, const std::string& fieldName, const std::string& value, const lSet& subset, const bool& additive){
  	switch(type){
      case FieldInfo::integerField:     SetFieldToConstantFromString<int>( fieldName, value,subset, additive);  break;
      case FieldInfo::localIndexField:  SetFieldToConstantFromString<localIndex>( fieldName, value,subset, additive);       break;
      case FieldInfo::globalIndexField:  SetFieldToConstantFromString<globalIndex>( fieldName, value,subset, additive);       break;
      case FieldInfo::realField:        SetFieldToConstantFromString<realT>( fieldName, value,subset, additive);       break;
      case FieldInfo::R1TensorField:    SetFieldToConstantFromString<R1Tensor>( fieldName, value,subset, additive);    break;
      case FieldInfo::R2TensorField:    SetFieldToConstantFromString<R2Tensor>( fieldName, value,subset, additive);    break;
      case FieldInfo::R2SymTensorField: SetFieldToConstantFromString<R2SymTensor>( fieldName, value,subset, additive); break;
      case FieldInfo::integerParameter:
      case FieldInfo::realParameter:
      case FieldInfo::numFieldTypes:
      default:
        throw GPException("SetConstantFieldFromString: Unrecognized field type "+ type);
    }
  };
    
  /// Set the value of a field using a named function
  void SetFieldEqualToFunction(FieldType fieldType, const std::string& fieldName, 
                               const std::string& functionName,
                               const array<string>& variables, 
                               const array<FieldType>& variable_types,
                               int component=0, realT time=0.0, realT dt=0.0);
  /// Set the value of a subset of a Yuliyafield using a named function
  void SetFieldEqualToFunction(FieldType fieldType, const std::string& fieldName, 
                               const std::string& functionName,
                               const array<string>& variables, 
                               const array<FieldType>& variable_types,
                               const lSet& aSet,
                               int component=0, realT time=0.0, realT dt=0.0);


/// sets subset of field to constant value 
  template< typename T >
  void CopyFieldSubset(const std::string& fieldName, const lArray1d& source, const lArray1d& target);
  
  void CopyFieldSubset( FieldType type, const std::string& fieldName,  const lArray1d& source, const lArray1d& target){
  	switch(type){
      case FieldInfo::integerField:     CopyFieldSubset<int>( fieldName, source, target);  break;
      case FieldInfo::localIndexField:  CopyFieldSubset<localIndex>( fieldName,  source, target);       break;
      case FieldInfo::globalIndexField:  CopyFieldSubset<globalIndex>( fieldName,  source, target);       break;
      case FieldInfo::realField:        CopyFieldSubset<realT>( fieldName, source, target);       break;
      case FieldInfo::R1TensorField:    CopyFieldSubset<R1Tensor>( fieldName, source, target);    break;
      case FieldInfo::R2TensorField:    CopyFieldSubset<R2Tensor>( fieldName, source, target);    break;
      case FieldInfo::R2SymTensorField: CopyFieldSubset<R2SymTensor>( fieldName, source, target); break; 
      case FieldInfo::integerParameter:
      case FieldInfo::realParameter:
      case FieldInfo::numFieldTypes:
      default:
        throw GPException("SetConstantFieldFromString: Unrecognized field type "+ type);
    }
  };


  //********************************************************************************************************************
  // Utilities
  //********************************************************************************************************************

  virtual void ExtractMapFromObjectForAssignGlobalObjectNumbers( const ObjectDataStructureBaseT& compositionObjectManager,
                                                                 array<gArray1d>& objectToCompositionObject ) = 0;

  /// pure virtual function that sets what objects are on the boundary of the domain
  virtual void SetDomainBoundaryObjects( const ObjectDataStructureBaseT* const referenceObject = nullptr) = 0 ;

  /// pure virtual function that sets what objects are external
  virtual void SetIsExternal( const ObjectDataStructureBaseT* const referenceObject = nullptr) = 0 ;

  /// function to reset m_globalToLocalMap based on m_localToGlobalMap
  virtual void ResetGlobalToLocalMap();

  template< typename T1, typename T2 >
  void LocalToGlobal( const T1& locals, T2& globals );

  localIndex GetNumGhosts() const
  {

    const array<integer>& ghostRank = this->GetFieldData<FieldInfo::ghostRank>();


    return ( std::count_if( ghostRank.begin(), ghostRank.end(), isGTE0 ) );

  }

  localIndex GetParentIndex( const localIndex index ) const
  {
    localIndex parentIndex = index;
    //while( m_parentIndex[parentIndex] != LOCALINDEX_MAX )
    while( !IsParent( parentIndex ) )
    {
      localIndex newParentIndex = m_parentIndex[parentIndex];
      if(newParentIndex == parentIndex){
        throw GPException("ObjectDataStructureT::GetParentIndex: Self referential parent index \n ");
        break;
      } else{
        parentIndex = newParentIndex;
      }
    }
    return parentIndex;
  }

  // This should be called "DontHaveParent" instead
  inline bool IsParent( const localIndex index )const {
    return m_parentIndex[index] == LOCALINDEX_MAX;
  }

  // Yeah, let's make a "DontHaveParent" and use it hereafter
  inline bool DontHaveParent( const localIndex index )const {
    return m_parentIndex[index] == LOCALINDEX_MAX;
  }


  /// generate vector containing all field names
  void GetAllFieldNames( array<string>& fieldNames ) const;

  /// get list of objects on the boundary of a computational domain
  virtual void ConstructListOfBoundaryObjects( lArray1d& objectList ) const ;

  /// get list of objects on the boundary of a computational domain
  virtual void ConstructListOfBoundaryObjects( gArray1d& objectList ) const ;


  /// builds a new set on this object given another objects set and the map between them
  void ConstructSetFromSetAndMap( const lSet& inputSet,
                                  const lArray2d& map,
                                  const std::string& newSetName );

  /// builds a new set on this object given another objects set and the map between them
  void ConstructSetFromSetAndMap( const lSet& inputSet,
                                  const array<lArray1d>& map,
                                  const std::string& newSetName );

  bool SplitObject( const localIndex indexToSplit,
                    const int rank,
                    localIndex& newIndex);

  void SplitObjectExcludeSets( const localIndex indexToSplit,
                               const int rank,
                               localIndex& newIndex);

  virtual bool SplitObject( const localIndex indexToSplit,
                            const int rank,
                            localIndex newIndices[2],
                            const bool forceSplit );

  virtual void CreateObject( const int rank,
                             localIndex& newIndex);

  void CopyObject( const localIndex source, const localIndex destination );

  void CopyObject( const localIndex source, const localIndex destination0, const localIndex destination1 );

  void CopyObjectWithExcludedSets( const localIndex source, const localIndex destination, const array<string>& excludeSets );

  void CopyObjectFields( const localIndex source, const localIndex destination );

  template< typename T >
  void CopyMemberFields( const localIndex source, const localIndex destination, std::map< std::string, array<T> >& member );

  template< typename T >
  void CopyMemberFields( const localIndex source, const localIndex destination0, const localIndex destination1, std::map< std::string, array<T> >& member );

  template< typename T >
  void CopyMemberFields( const localIndex source, const lArray1d& destination, std::map< std::string, array<T> >& member );


  void UpdateMaximumGlobalNumber();

  virtual void DeserializeObjectField(const std::string& name, const array<real64>& field)
  {
    (void)name;
    (void)field;

    throw GPException("DeserializeObjectField - undefined in base");
  }

  virtual void DeserializeObjectFields(const array<string>& names, const array<array<real64>>& fields)
  {
    if(names.size() != fields.size())
      throw GPException("DeserializeObjectFields - input array sizes do not match");
    for(localIndex i = 0; i < names.size(); i++)
      DeserializeObjectField(names[i], fields[i]);
  }

  //********************************************************************************************************************
  // buffer packing/unpacking
  //********************************************************************************************************************


  /// packs a list of fields into a bufvector
  template< typename T_indices, typename T_buffer >
  unsigned int PackFieldsIntoBuffer( T_buffer& buffer,
                                     const array<string>& fieldNames,
                                     const T_indices& localIndices,
                                     const bool doBufferPacking=true ) const;

  /// unpacks a list of fields from a buffer pointer. Changes the pointer!
  template< typename T_indices >
  unsigned int UnpackFieldsFromBuffer( const char*& buffer,
                                       const array<string>& fieldNames,
                                       const T_indices& localIndices );

  /// a wrapper of PackFieldsIntoBuffer to pack ALL fields.
  template< typename T_indices, typename T_buffer >
  unsigned int PackAllFieldsIntoBuffer( T_buffer& buffer,
                                        const T_indices& localIndices ) const;

  /// a wrapper of UnpackFieldsFromBuffer to unpack ALL fields. Changes the pointer!
  template< typename T_indices >
  unsigned int UnpackAllFieldsFromBuffer( const char*& buffer,
                                          const T_indices& localIndices );



  template< typename T_indices >
  unsigned int PackRelationsIntoBuffer( bufvector& buffer,
                                        const T_indices& localIndices,
                                        const bool packRelationsToGlobal ) const;

  unsigned int UnpackRelationsFromBuffer( const char*& buffer,
                                          const lArray1d& localIndices,
                                          const bool unpackRelationsToLocal );

//  virtual void SetRelationPointers() = 0;


  template< typename T_indices >
  unsigned int PackSets( const T_indices& sendlist,
                         bufvector& buffer ) const;

  unsigned int UnpackSets( const char*& buffer ) ;

  ObjectType GetObjectType() const
  {
    return m_objectType;
  }

  /// function to check what ObjectType *this is.
  void CheckObjectType( const ObjectType type ) const
  {
    if( m_objectType != type )
    {
      throw GPException("ObjectDataStructureT::CheckObjectType: Failure of check \n ");
    }
  }






//  virtual unsigned int PackObjectData( bufvector& buffer, const lArray1d& indices ) = 0;
//  virtual unsigned int PackObjectData( bufvector& buffer, const lSet& indices ) = 0;
//  virtual unsigned int UnpackObjectData( const char*& buffer, lArray1d& unpackedLocalIndices ) = 0;
  template< typename T_indices >
  unsigned int PackBaseObjectData( bufvector& buffer,
                                   const T_indices& indices,
                                   const bool packFields,
                                   const bool packMaps,
                                   const bool packSets,
                                   const bool packRelationsToGlobal ) const;

  unsigned int UnpackBaseObjectData( const char*& buffer,
                                     lArray1d& unpackedLocalIndices,
                                     lArray1d& newLocalIndices,
                                     const bool unpackFields,
                                     const bool unpackMaps,
                                     const bool unpackSets,
                                     const bool unpackRelationsToLocal );


private:
  /// Default Constructor
  ObjectDataStructureBaseT();


  template< typename T, typename T_indices, typename T_buffer >
  unsigned int PackMemberFieldsIntoBuffer( T_buffer& buffer,
                                           const std::map< std::string, array<T> >& member,
                                           const T_indices& localIndices,
                                           const bool doBufferPacking ) const;

  /// pack the fields of member into the bufvector
  template< typename T, typename T_indices, typename T_buffer >
  unsigned int PackMemberFieldsIntoBuffer( T_buffer& buffer,
                                           const array<string>& fieldNames,
                                           const std::map< std::string, array<T> >& member,
                                           const T_indices& localIndices,
                                           const bool doBufferPacking ) const;


  /// unpack the fields of member from a buffer
  template< typename T, typename T_indices >
  unsigned int UnpackMemberFieldsFromBuffer( const char*& buffer,
                                             std::map< std::string, array<T> >& member,
                                             const T_indices& localIndices );

  /// unpack the fields of member from a buffer
  template< typename T, typename T_indices >
  unsigned int UnpackMemberFieldsFromBuffer( const char*& buffer,
                                             const array<string>& fieldNames,
                                             std::map< std::string, array<T> >& member,
                                             const T_indices& localIndices );

  /// pack a single field into a bufvector
  template< typename T, typename T_indices >
  unsigned int PackFieldIntoBuffer( bufvector& buffer,
                                    const std::string& fieldName,
                                    const array<T>& field,
                                    const T_indices& localIndices,
                                    const bool doBufferPacking ) const;

  /// pack a single field into a bufvector
  template< typename T, typename T_indices >
  unsigned int PackFieldIntoBuffer( char*& buffer,
                                    const std::string& fieldName,
                                    const array<T>& field,
                                    const T_indices& localIndices,
                                    const bool doBufferPacking ) const;

  /// unpack a single field from a buffer
  template< typename T, typename T_indices >
  unsigned int UnpackFieldFromBuffer( const char*& buffer,
                                      const std::string& fieldName,
                                      array<T>& field,
                                      const T_indices& localIndices );


  template< typename T, typename T_indices >
  unsigned int PackRelationIntoBuffer( bufvector& buffer,
                                       const std::map<std::string,T>& map,
                                       const T_indices& localIndices,
                                       const bool packRelationsToGlobal ) const;

  template< typename T >
  unsigned int UnpackRelationFromBuffer( const char*& buffer,
                                         std::map<std::string,T>& map,
                                         const lArray1d& localIndices,
                                         const bool unpackRelationsToLocal );



  //********************************************************************************************************************
  // data members
  //********************************************************************************************************************

  // these public members should be protected...or private...but I am lazy, and really no need yet.
public:

  /// lookup of global index from local indexes
  gArray1d m_localToGlobalMap;

  /// map of local index from global indices. Not a very fast operation.
  std::map<globalIndex,localIndex> m_globalToLocalMap;

  /// maximum global number
  globalIndex m_maxGlobalNumber;

  /// sets... i.e. nodeset, facesets, elementsets...
  std::map< std::string, lSet >             m_Sets;

  /// boundary conditions on this object
  array<BoundaryConditionBase*>  m_bcData;

  



  size_t DataLengths() const {return m_DataLengths;}
protected:

  /// type of object
  const ObjectType m_objectType;

  /// member to hold the length of all data arrays and maps
  size_t m_DataLengths;

  /// map to hold arrays of unsigned long
  std::map< std::string, lArray1d >             m_LocalIndexData;

  /// map to hold arrays of unsigned long long
  std::map< std::string, gArray1d >             m_GlobalIndexData;

  /// map to hold arrays of integers
  std::map< std::string, array<integer> >             m_IntegerData;

  /// map to hold arrays of reals
  std::map< std::string, array<real64> >             m_realData;

  /// map to hold arrays of R1Tensors
  std::map< std::string, array<R1Tensor> >    m_R1TensorData;

  /// map to hold arrays of R2Tensors
  std::map< std::string, array<R2Tensor> >    m_R2TensorData;

  /// map to hold arrays of R2SymTensors
  std::map< std::string, array<R2SymTensor> > m_R2SymTensorData;

  /// map to hold arrays of arrays of R1Tensors
  std::map< std::string, array<array<R1Tensor>>>  m_ArrayR1TensorData;

  /// map to hold one-to-one maps
  std::map< std::string, OneToOneRelation >             m_OneToOneMaps;

  /// map to hold fixed size one-to-many maps
  std::map< std::string, FixedOneToManyRelation >             m_FixedOneToManyMaps;

  /// map to hold variable sized one-to-many maps
  std::map< std::string, OrderedVariableOneToManyRelation >   m_VariableOneToManyMaps;

  /// map to hold variable sized one-to-many maps
  std::map< std::string, OrderedVariableOneToManyToManyRelation >   m_VariableOneToManyToManyMaps;

  /// map to hold unordered variable one-to-many maps
  std::map< std::string, UnorderedVariableOneToManyRelation > m_UnorderedVariableOneToManyMaps;


public:

  array<integer>& m_isExternal;
  array<integer>& m_processColor;
  localIndex m_processColorPad;

  OrderedVariableOneToManyRelation& m_childIndices;
  OneToOneRelation& m_parentIndex;
  int m_rank = -1;

  template< typename TYPE >
  void ConstructVariableOneToManyFromInverse( const std::string& mapname, const array<TYPE>& inverse );

  template< typename TYPE, typename TYPE2 >
  void AddToVariableOneToManyFromInverse( const std::string& mapname, const array<TYPE>& inverseMap, const TYPE2& newIndices );


  const OneToOneRelation& GetOneToOneMap( const std::string& name ) const { return stlMapLookup(m_OneToOneMaps,name); }
  OneToOneRelation& GetOneToOneMap( const std::string& name ) { return stlMapLookup(m_OneToOneMaps,name); }

  const FixedOneToManyRelation& GetFixedOneToManyRelation( const std::string& name ) const { return stlMapLookup(m_FixedOneToManyMaps,name); }
  FixedOneToManyRelation& GetFixedOneToManyRelation( const std::string& name ) { return stlMapLookup(m_FixedOneToManyMaps,name); }


  const OrderedVariableOneToManyRelation& GetVariableOneToManyMap( const std::string& name ) const { return stlMapLookup(m_VariableOneToManyMaps,name); }
  OrderedVariableOneToManyRelation& GetVariableOneToManyMap( const std::string& name ) { return stlMapLookup(m_VariableOneToManyMaps,name); }

  const OrderedVariableOneToManyToManyRelation& GetVariableOneToManyToManyMap( const std::string& name ) const { return stlMapLookup(m_VariableOneToManyToManyMaps,name); }
  OrderedVariableOneToManyToManyRelation& GetVariableOneToManyToManyMap( const std::string& name ) { return stlMapLookup(m_VariableOneToManyToManyMaps,name); }

  const UnorderedVariableOneToManyRelation& GetUnorderedVariableOneToManyMap( const std::string& name ) const { return stlMapLookup(m_UnorderedVariableOneToManyMaps,name); }
  UnorderedVariableOneToManyRelation& GetUnorderedVariableOneToManyMap( const std::string& name ) { return stlMapLookup(m_UnorderedVariableOneToManyMaps,name); }


  const array<lSet>* GetUnorderedVariableOneToManyMapPointer( const std::string& name ) const { return stlMapLookupPointer(m_UnorderedVariableOneToManyMaps,name); }
  array<lSet>* GetUnorderedVariableOneToManyMapPointer( const std::string& name ) { return stlMapLookupPointer(m_UnorderedVariableOneToManyMaps,name); }

  //const lSet* GetNodeMapPointer(localIndex indx){return 0;};

private:
  /// Assignment operator
  ObjectDataStructureBaseT& operator=( const ObjectDataStructureBaseT&);

};

typedef ObjectDataStructureBaseT::ObjectType ObjectManagerType;



// *********************************************************************************************************************
// *********************************************************************************************************************

template< FieldKey FIELDKEY >
int ObjectDataStructureBaseT::AddKeyedDataField()
{
  typedef typename Field<FIELDKEY>::Type T;


  // get data member map corresponding to the type
  std::map< std::string, array<T> >& DataArrayMember = GetDataMemberMap<array<T> >();

  std::string name = Field<FIELDKEY>::Name();

  // the [] will create the new map entry
  DataArrayMember[name];

  // resize the new field
  DataArrayMember[name].resize(m_DataLengths);

  return 0;
}
















// *********************************************************************************************************************
/**
 * @author R. Settgast
 * @tparam FIELDKEY the key of the field to return
 * @return const reference to the field
 *
 * returns a const reference to the field specified by the field key.
 *
 */
template< FieldKey FIELDKEY>
inline const array<typename Field<FIELDKEY>::Type>& ObjectDataStructureBaseT::GetFieldData( ) const
{
  // define a typedef of the "type" corresponding to FIELDKEY for easier reading
  typedef typename Field<FIELDKEY>::Type T;

  // get the data member map corresponding to the type
  const std::map< std::string, array<T> >& objectmap = GetDataMemberMap<array<T> >();

  // get an iterator to the map entry that corresponds with the field name and return it.
  const typename std::map< std::string, array<T> >::const_iterator iobjectmap = objectmap.find(Field<FIELDKEY>::Name());

  if( iobjectmap == objectmap.end() ){
  	std::string errStr = Field<FIELDKEY>::Name();
    throw GPException("didn't find() requested field '" + errStr + "' in ObjectDataStructureBaseT /n");
  }
  return(iobjectmap->second);

}


template< FieldKey FIELDKEY>
inline const array<typename Field<FIELDKEY>::Type>* ObjectDataStructureBaseT::GetFieldDataPointer( ) const
{
  // define a typedef of the "type" corresponding to FIELDKEY for easier reading
  typedef typename Field<FIELDKEY>::Type T;

  // get the data member map corresponding to the type
  const std::map< std::string, array<T> >& objectmap = GetDataMemberMap<array<T> >();

  // get an iterator to the map entry that corresponds with the field name and return it.
  const typename std::map< std::string, array<T> >::const_iterator iobjectmap = objectmap.find(Field<FIELDKEY>::Name());


  const array<typename Field<FIELDKEY>::Type>* rval = nullptr;
  if( iobjectmap != objectmap.end() )
  {
    rval = &(iobjectmap->second);
  }

  return(rval);

}


// *********************************************************************************************************************
/**
 * @author R. Settgast
 * @tparam FIELDKEY the key of the field to return
 * @return const reference to the field
 *
 * returns a const reference to the field specified by the field key.
 *
 */
template< typename T >
inline const array<T>& ObjectDataStructureBaseT::GetFieldData( const std::string& fieldName ) const
{
  // get the data member map corresponding to the type
  const std::map< std::string, array<T> >& objectmap = GetDataMemberMap<array<T> >();

  // get an iterator to the map entry that corresponds with the field name and return it.
  const typename std::map< std::string, array<T> >::const_iterator iobjectmap = objectmap.find(fieldName);
  if( iobjectmap == objectmap.end() )
    throw GPException("didn't find() requested field '" + fieldName + "' in ObjectDataStructureBaseT /n");


  return(iobjectmap->second);
}

template< typename T >
inline const array<T>* ObjectDataStructureBaseT::GetFieldDataPointer( const std::string& fieldName ) const
{
  // get the data member map corresponding to the type
  const std::map< std::string, array<T> >& objectmap = GetDataMemberMap<array<T> >();

  // get an iterator to the map entry that corresponds with the field name and return it.
  const typename std::map< std::string, array<T> >::const_iterator iobjectmap = objectmap.find(fieldName);

  const array<T>* rval = nullptr;
  if( iobjectmap != objectmap.end() )
  {
    rval = &(iobjectmap->second);
  }


  return(rval);
}

// *********************************************************************************************************************
/**
 * @author S. Walsh
 * @return reference to the set
 *
 * returns a reference to the set specified by the setName.
 *
 */
inline
lSet& ObjectDataStructureBaseT::GetSet( const std::string& setName )
{
  // get an iterator to the map entry that corresponds with the field name and return it.
  std::map< std::string, lSet >::iterator itr = m_Sets.find(setName);
  if( itr == m_Sets.end() )
    throw GPException("GetSet: didn't find() requested set '" + setName + "' in ObjectDataStructureBaseT /n");
  return(itr->second);
}

// *********************************************************************************************************************
/**
 * @author S. Walsh
 * @return const reference to the set
 *
 * returns a const reference to the set specified by the setName.
 *
 */
inline
const lSet& ObjectDataStructureBaseT::GetSet( const std::string& setName ) const
{
  // get an iterator to the map entry that corresponds with the field name and return it.
  const std::map< std::string, lSet >::const_iterator itr = m_Sets.find(setName);
  if( itr == m_Sets.end() )
    throw GPException("GetSet: didn't find() requested set '" + setName + "' in ObjectDataStructureBaseT /n");
  return(itr->second);
}


template< typename T >
void ObjectDataStructureBaseT::SetSubsetValue(const std::string& fieldName,
                                              const std::string& setName,
                                              const T& value){
    lSet& theSet = GetSet(setName);
   std::set<localIndex>::iterator itr = theSet.begin();
   std::set<localIndex>::iterator iend = theSet.end();

   array<T>& field = GetFieldData<T>(fieldName);
   for(; itr != iend;++itr){
     field[*itr] = value;
   }
}

/**
 * @author R. Settgast
 * @param[in] fieldName the name of the desired field
 * @return true if the field exists
 */
template< typename T >
bool ObjectDataStructureBaseT::HasField( const std::string& fieldName ) const
{
  // get the data member map corresponding to the type
  const std::map< std::string, array<T> >& objectmap = GetDataMemberMap<array<T> >();

  // get an iterator to the map entry that corresponds with the field name
  const typename std::map< std::string, array<T> >::const_iterator iobjectmap = objectmap.find(fieldName);

  return(iobjectmap != objectmap.end());
}

/**
 * @author R. Settgast
 * @param pointer[in] pointer to set
 * this function sets the value of a const* to the location of a field array
 */
template< FieldKey FIELDKEY >
void ObjectDataStructureBaseT::SetConstPointer( array<typename Field<FIELDKEY>::Type>* const& pointer )
{
  // define a typedef of the "type" corresponding to FIELDKEY for easier reading
  typedef typename Field<FIELDKEY>::Type T;

  // get the data member map corresponding to the type
  std::map< std::string, array<T> >& objectmap = GetDataMemberMap<array<T> >();

  const typename std::map< std::string, array<T> >::iterator iobjectmap = objectmap.find(Field<FIELDKEY>::Name());

  array<T>** temp = const_cast< array<T>**>(&pointer);

  if( iobjectmap == objectmap.end() )
    *temp = nullptr;
  else
    *temp = &(iobjectmap->second);

  return;
}




//**********************************************************************************************************************
/**
 * @author R. Settgast
 *
 * @param[out] buffer the bufvector that is used as the packing buffer
 * @param[in] fieldNames the list of field names that are to be packed into buffer
 * @param[in] localIndices the list of local indices to pack into the buffer
 * @return the total number of characters packed into buffer by this function
 */
template< typename T_indices, typename T_buffer >
unsigned int ObjectDataStructureBaseT::PackFieldsIntoBuffer( T_buffer& buffer,
                                                             const array<string>& fieldNames,
                                                             const T_indices& localIndices,
                                                             const bool doBufferPacking ) const
{
  unsigned int packedSize = 0;

  packedSize += PackMemberFieldsIntoBuffer( buffer, fieldNames, m_IntegerData, localIndices, doBufferPacking );
  packedSize += PackMemberFieldsIntoBuffer( buffer, fieldNames, m_LocalIndexData, localIndices, doBufferPacking );
  packedSize += PackMemberFieldsIntoBuffer( buffer, fieldNames, m_GlobalIndexData, localIndices, doBufferPacking );
  packedSize += PackMemberFieldsIntoBuffer( buffer, fieldNames, m_realData, localIndices, doBufferPacking );
  packedSize += PackMemberFieldsIntoBuffer( buffer, fieldNames, m_R1TensorData, localIndices, doBufferPacking );
  packedSize += PackMemberFieldsIntoBuffer( buffer, fieldNames, m_R2TensorData, localIndices, doBufferPacking );
  packedSize += PackMemberFieldsIntoBuffer( buffer, fieldNames, m_R2SymTensorData, localIndices, doBufferPacking );

  return packedSize;
}




/**
 * @author R. Settgast
 * @param[out] buffer the bufvector that is used as the packing buffer
 * @param[in] localIndices the list of local indices to pack into the buffer
 * @return the total number of characters packed into buffer by this function
 */
template< typename T_indices, typename T_buffer >
unsigned int ObjectDataStructureBaseT::PackAllFieldsIntoBuffer( T_buffer& buffer,
                                                                const T_indices& localIndices ) const
{
  unsigned int packedSize = 0;

  packedSize += PackMemberFieldsIntoBuffer( buffer, m_IntegerData,localIndices, true );
  packedSize += PackMemberFieldsIntoBuffer( buffer, m_LocalIndexData,localIndices, true );
  packedSize += PackMemberFieldsIntoBuffer( buffer, m_GlobalIndexData,localIndices, true );
  packedSize += PackMemberFieldsIntoBuffer( buffer, m_realData,localIndices, true );
  packedSize += PackMemberFieldsIntoBuffer( buffer, m_R1TensorData,localIndices, true );
  packedSize += PackMemberFieldsIntoBuffer( buffer, m_R2TensorData,localIndices, true );
  packedSize += PackMemberFieldsIntoBuffer( buffer, m_R2SymTensorData,localIndices, true );

  return( packedSize );
}

/**
 * @author R. Settgast
 * @param[in,out] buffer character pointer that holds the data to be unpacked. It is incremented
 * @param[in] fieldNames the list of field names that are to be packed into buffer
 * @param[in] localIndices the list of local indices unpack the data to
 * @return total number of characters unpacked by the function
 */
template< typename T_indices >
unsigned int ObjectDataStructureBaseT::UnpackFieldsFromBuffer( const char*& buffer,
                                                              const array<string>& fieldNames,
                                                              const T_indices& localIndices )
{
  unsigned int sizeOfUnpackedChars = 0;

  sizeOfUnpackedChars += UnpackMemberFieldsFromBuffer( buffer, fieldNames, m_IntegerData, localIndices );
  sizeOfUnpackedChars += UnpackMemberFieldsFromBuffer( buffer, fieldNames, m_LocalIndexData, localIndices );
  sizeOfUnpackedChars += UnpackMemberFieldsFromBuffer( buffer, fieldNames, m_GlobalIndexData, localIndices );
  sizeOfUnpackedChars += UnpackMemberFieldsFromBuffer( buffer, fieldNames, m_realData, localIndices ) ;
  sizeOfUnpackedChars += UnpackMemberFieldsFromBuffer( buffer, fieldNames, m_R1TensorData, localIndices );
  sizeOfUnpackedChars += UnpackMemberFieldsFromBuffer( buffer, fieldNames, m_R2TensorData, localIndices );
  sizeOfUnpackedChars += UnpackMemberFieldsFromBuffer( buffer, fieldNames, m_R2SymTensorData, localIndices );

  return sizeOfUnpackedChars;
}

/**
 * @author R. Settgast
 * @param[in,out] buffer character pointer that holds the data to be unpacked. It is incremented
 * @param[in] localIndices the list of local indices unpack the data to
 * @return total number of characters unpacked by the function
 */
template< typename T_indices >
unsigned int ObjectDataStructureBaseT::UnpackAllFieldsFromBuffer( const char*& buffer,
                                                                  const T_indices& localIndices )
{
  unsigned int sizeOfUnpackedChars = 0;

  sizeOfUnpackedChars += UnpackMemberFieldsFromBuffer( buffer, m_IntegerData, localIndices );
  sizeOfUnpackedChars += UnpackMemberFieldsFromBuffer( buffer, m_LocalIndexData, localIndices );
  sizeOfUnpackedChars += UnpackMemberFieldsFromBuffer( buffer, m_GlobalIndexData, localIndices );
  sizeOfUnpackedChars += UnpackMemberFieldsFromBuffer( buffer, m_realData, localIndices ) ;
  sizeOfUnpackedChars += UnpackMemberFieldsFromBuffer( buffer, m_R1TensorData, localIndices );
  sizeOfUnpackedChars += UnpackMemberFieldsFromBuffer( buffer, m_R2TensorData, localIndices );
  sizeOfUnpackedChars += UnpackMemberFieldsFromBuffer( buffer, m_R2SymTensorData, localIndices );

  return( sizeOfUnpackedChars );
}



template< typename T, typename T_indices, typename T_buffer >
unsigned int ObjectDataStructureBaseT::PackMemberFieldsIntoBuffer( T_buffer& buffer,
                                                                   const std::map< std::string, array<T> >& member,
                                                                   const T_indices& localIndices,
                                                                   const bool doBufferPacking ) const
{
  unsigned int sizeOfPackedChars = 0;

  for( typename std::map< std::string, array<T> >::const_iterator i = member.begin() ; i!=member.end() ; ++i )
  {
    const std::string fieldName = i->first;
    const array<T>& field = i->second;
    sizeOfPackedChars += PackFieldIntoBuffer( buffer, fieldName, field, localIndices, doBufferPacking );
  }
  return sizeOfPackedChars ;
}


// *********************************************************************************************************************
/**
 * @author R. Settgast
 * @tparam T the type of data that the member holds.
 * @param[out] buffer the bufvector that is used as the packing buffer
 * @param[in] fieldNames the list of field names that are to be packed into buffer
 * @param[in] member the data member that should hold fields that are keyed on fieldNames
 * @param[in] localIndices the list of local indices to pack into the buffer
 * @return the total number of characters packed into buffer by this function
 */
template< typename T, typename T_indices, typename T_buffer  >
unsigned int ObjectDataStructureBaseT::PackMemberFieldsIntoBuffer( T_buffer& buffer,
                                                                   const array<string>& fieldNames,
                                                                   const std::map< std::string, array<T> >& member,
                                                                   const T_indices& localIndices,
                                                                   const bool doBufferPacking ) const
{
  unsigned int sizeOfPackedChars = 0;

  array<string>::const_iterator fieldName=fieldNames.begin() ;

  for(  ; fieldName!=fieldNames.end(); ++fieldName )
  {
    const typename std::map< std::string, array<T> >::const_iterator i = member.find(*fieldName);
    if( i != member.end() )
    {
      const array<T>& field = i->second;
      sizeOfPackedChars += PackFieldIntoBuffer( buffer, *fieldName, field, localIndices, doBufferPacking );
    }
  }
  return sizeOfPackedChars ;
}

template< typename T, typename T_indices >
unsigned int ObjectDataStructureBaseT::UnpackMemberFieldsFromBuffer( const char*& buffer,
                                                                     std::map< std::string, array<T> >& member,
                                                                     const T_indices& localIndices )
{
  unsigned int sizeOfUnpackedChars = 0;

  for( typename std::map< std::string, array<T> >::iterator i = member.begin() ; i!=member.end() ; ++i )
  {
    const std::string fieldName = i->first;
    array<T>& field = i->second;

    sizeOfUnpackedChars += UnpackFieldFromBuffer(  buffer, fieldName, field, localIndices );
  }
  return sizeOfUnpackedChars;
}


// *********************************************************************************************************************
/**
 * @author R. Settgast
 * @tparam T the type of data that the member holds.
 * @param[in,out] buffer character pointer that holds the data to be unpacked. It is incremented
 * @param[in] fieldNames the list of field names that are to be packed into buffer
 * @param[in] member the data member that should hold fields that are keyed on fieldNames
 * @param[in] localIndices the list of local indices unpack the data to
 * @return total number of characters unpacked by the function.
 */
template< typename T, typename T_indices >
unsigned int ObjectDataStructureBaseT::UnpackMemberFieldsFromBuffer( const char*& buffer,
                                                            const array<string>& fieldNames,
                                                            std::map< std::string, array<T> >& member,
                                                            const T_indices& localIndices )
{
  unsigned int sizeOfUnpackedChars = 0;

  array<string>::const_iterator fieldName=fieldNames.begin() ;

  for(  ; fieldName!=fieldNames.end(); ++fieldName )
  {
    const typename std::map< std::string, array<T> >::iterator i = member.find(*fieldName);
    if( i != member.end() )
    {
      array<T>& field = i->second;
      sizeOfUnpackedChars += UnpackFieldFromBuffer(  buffer, *fieldName, field, localIndices );
    }
  }
  return sizeOfUnpackedChars;
}

// *********************************************************************************************************************
/**
 * @author R. Settgast
 * @tparam T the type of data that the member holds.
 * @param[out] buffer the bufvector that is used as the packing buffer
 * @param[in] fieldName the name of the field to be packed
 * @param field[out] the individual field data to be packed
 * @param[in] localIndices the list of local indices to pack into the buffer
 * @return the total number of characters packed into buffer by this function
 */
template< typename T, typename T_indices >
unsigned int ObjectDataStructureBaseT::PackFieldIntoBuffer( bufvector& buffer,
                                                            const std::string& fieldName ,
                                                            const array<T>& field,
                                                            const T_indices& localIndices,
                                                            const bool doBufferPacking ) const
{
  unsigned int sizeOfPackedChars = 0;

  if( doBufferPacking )
  {

#ifdef PACK_FIELD_NAMES_IN_MPI_BUFFER
  sizeOfPackedChars += buffer.PackFieldname(fieldName);
#endif

    for( typename T_indices::const_iterator i = localIndices.begin() ;
        i != localIndices.end() ; ++i )
    {
      sizeOfPackedChars += buffer.Pack( field[*i] );
    }
  }
  else
  {

#ifdef PACK_FIELD_NAMES_IN_MPI_BUFFER
    sizeOfPackedChars += sizeof(unsigned int);
    /*sizeOfPackedChars += fieldName.size();*/
    sizeOfPackedChars += bufvector::sizeOfPackedFieldString;

#endif
    sizeOfPackedChars += localIndices.size()*sizeof(T);
  }
  return sizeOfPackedChars;
}

template< typename T, typename T_indices >
unsigned int ObjectDataStructureBaseT::PackFieldIntoBuffer( char*& buffer,
                                                            const std::string& fieldName ,
                                                            const array<T>& field,
                                                            const T_indices& localIndices,
                                                            const bool doBufferPacking ) const
{
  unsigned int sizeOfPackedChars = 0;

  if( doBufferPacking )
  {
#ifdef PACK_FIELD_NAMES_IN_MPI_BUFFER
    sizeOfPackedChars += bufvector::PackFieldname(buffer, fieldName);
#endif


    for( typename T_indices::const_iterator i = localIndices.begin() ;
         i != localIndices.end() ; ++i )
    {
      sizeOfPackedChars += bufvector::Pack( buffer, field[*i] );
    }
  }
  else
  {

#ifdef PACK_FIELD_NAMES_IN_MPI_BUFFER
    sizeOfPackedChars += sizeof(unsigned int);
    /*sizeOfPackedChars += fieldName.size();*/
    sizeOfPackedChars += bufvector::sizeOfPackedFieldString;
#endif
    sizeOfPackedChars += localIndices.size()*sizeof(T);
  }


  return sizeOfPackedChars;

}


// *********************************************************************************************************************
/**
 * @author R. Settgast
 * @tparam T the type of data that the member holds.
 * @param[in,out] buffer character pointer that holds the data to be unpacked. It is incremented
 * @param[in] fieldName the name of the field to be unpacked
 * @param field[in] the individual field data to be unpacked
 * @param[in] localIndices the list of local indices to pack into the buffer
 * @return total number of characters unpacked by the function.
 */
template< typename T, typename T_indices >
unsigned int ObjectDataStructureBaseT::UnpackFieldFromBuffer( const char*& buffer,
                                                              const std::string& fieldName ,
                                                              array<T>& field,
                                                              const T_indices& localIndices )
{
  unsigned int sizeOfUnpackedChars = 0;
  std::string readFieldName;

#ifdef PACK_FIELD_NAMES_IN_MPI_BUFFER
  sizeOfUnpackedChars += bufvector::Unpack( buffer, readFieldName );

  if( !bufvector::FieldnameMatchesIdString(fieldName,readFieldName) )
    throw GPException("ObjectDataStuctureBaseT::UnpackFieldFromBuffer: field names ("+ fieldName + "," +readFieldName +") do not match up!\n");
#endif

#if 0
  for( lArray1d::const_iterator localIndex = localIndices.begin() ;
      localIndex != localIndices.end() ; ++localIndex )
  {
    sizeOfUnpackedChars += bufvector::Unpack( buffer, field[*localIndex] );
  }
#else

  const T* tbuffer = reinterpret_cast<const T*>(buffer);
  for( typename T_indices::const_iterator i = localIndices.begin() ;
      i != localIndices.end() ; ++i )
  {
    sizeOfUnpackedChars += sizeof(T);
    field[*i] = *tbuffer;
    ++tbuffer;
  }
  buffer = reinterpret_cast<const char*>(tbuffer);

#endif

  return sizeOfUnpackedChars;
}






// *********************************************************************************************************************
/**
 * @author S. Walsh
 * @tparam fileName  the file containing the field
 * @tparam fieldName the key of the field to read
 *
 */
template< typename T >
void ObjectDataStructureBaseT::ReadAsciiFieldData(const std::string& fieldName, const std::string& fileName)
{
  // get the data member map corresponding to the type
  std::map< std::string, array<T> >& objectmap = GetDataMemberMap<array<T> >();

  // get an iterator to the map entry that corresponds with the field name
  typename std::map< std::string, array<T> >::iterator iobjectmap = objectmap.find(fieldName);
  if( iobjectmap == objectmap.end() )
    throw GPException("didn't find() requested field "+ fieldName +" in ObjectDataStructureBaseT /n");

  std::ifstream fStream(fileName.c_str());
  array<T>& field = iobjectmap->second;

  // read data
  for( size_t i =0; i < field.size(); ++i){
  	if(fStream.good()){ fStream >> field[i];}
  	else {
      throw GPException("Error ReadAsciiFieldData: Data vector in "+ fileName +" (size " +toString(i+1) + ") is shorter than field " + fieldName + "(size " +toString(field.size()+1) + ") /n");
    }
  }

  fStream.close();
}


// ascii field data referenced by global index
// format: globalId  data
template< typename T >
void ObjectDataStructureBaseT::ReadIndexedAsciiFieldData(const std::string& fieldName, const std::string& fileName)
{
  // get the data member map corresponding to the type
  std::map< std::string, array<T> >& objectmap = GetDataMemberMap<array<T> >();

  // get an iterator to the map entry that corresponds with the field name
  typename std::map< std::string, array<T> >::iterator iobjectmap = objectmap.find(fieldName);
  if( iobjectmap == objectmap.end() )
    throw GPException("didn't find() requested field "+ fieldName +" in ObjectDataStructureBaseT /n");

  std::ifstream fStream(fileName.c_str());
  array<T>& field = iobjectmap->second;

  // read data
  std::string lineStr;
  localIndex aLocalIndex;
  globalIndex aGlobalIndex;
  T datum;
  while( std::getline(fStream,lineStr) ){ //fStream.good()
	std::istringstream iss(lineStr);
	iss >> aGlobalIndex >> datum;
    if( !iss.fail() ){ // not comment line
  		if( isMember(aGlobalIndex,m_globalToLocalMap) ){  // could be on another processor
  		  aLocalIndex = m_globalToLocalMap[aGlobalIndex];
  		  field[aLocalIndex] = datum;
  		}
  	}
  }

  fStream.close();
}

template< typename T >
void ObjectDataStructureBaseT::ReadAsciiFieldData(const std::string& fieldName, const std::string& fileName,
                                                  const lSet& subset)
{
  // get the data member map corresponding to the type
  std::map< std::string, array<T> >& objectmap = GetDataMemberMap<array<T> >();

  // get an iterator to the map entry that corresponds with the field name
  typename std::map< std::string, array<T> >::iterator iobjectmap = objectmap.find(fieldName);
  if( iobjectmap == objectmap.end() )
    throw GPException("didn't find() requested field "+ fieldName +" in ObjectDataStructureBaseT /n");

  std::ifstream fStream(fileName.c_str());
  array<T>& field = iobjectmap->second;

  // read data

  lSet::const_iterator isubset=subset.begin();
  for( size_t i =0; i < subset.size(); ++i, ++isubset){
  	if(fStream.good()){ fStream >> field[*isubset];}
  	else {
      throw GPException("Error ReadAsciiFieldData: Data vector in "+ fileName +" (size " +toString(i+1) + ") is shorter than field " + fieldName + " set (size " +toString(subset.size()+1) + ") /n");
    }
  }

  fStream.close();
}

// *********************************************************************************************************************
/**
 * @author S. Walsh
 * @tparam fileName  the file containing the field
 * @tparam fieldName the key of the field to read
 *
 */
template< typename T >
void ObjectDataStructureBaseT::WriteAsciiFieldData(const std::string& fieldName, const std::string& fileName, bool append)
{
  // get the data member map corresponding to the type
  std::map< std::string, array<T> >& objectmap = GetDataMemberMap<array<T> >();

  // get an iterator to the map entry that corresponds with the field name
  typename std::map< std::string, array<T> >::iterator iobjectmap = objectmap.find(fieldName);
  if( iobjectmap == objectmap.end() )
    throw GPException("didn't find() requested field "+ fieldName +" in ObjectDataStructureBaseT /n");

  std::ofstream fStream;

  if(append){
    fStream.open(fileName.c_str(), std::ios::out | std::ios::app );
  } else {
    fStream.open(fileName.c_str(), std::ios::out );
  }

  fStream.precision(std::numeric_limits<realT>::digits10 + 2);

  array<T>& field = iobjectmap->second;

  // write data
  for( size_t i =0; i < field.size(); ++i){
  	 fStream << field[i];
  	 fStream << "\n";
  }

  fStream.close();
}

template< typename T >
void ObjectDataStructureBaseT::WriteAsciiFieldData(const std::string& fieldName, const std::string& fileName,
                                                  const lSet& subset, bool append)
{
  // get the data member map corresponding to the type
  std::map< std::string, array<T> >& objectmap = GetDataMemberMap<array<T> >();

  // get an iterator to the map entry that corresponds with the field name
  typename std::map< std::string, array<T> >::iterator iobjectmap = objectmap.find(fieldName);
  if( iobjectmap == objectmap.end() )
    throw GPException("didn't find() requested field "+ fieldName +" in ObjectDataStructureBaseT /n");

  std::ofstream fStream;

  if(append){
    fStream.open(fileName.c_str(), std::ios::out | std::ios::app );
  } else {
    fStream.open(fileName.c_str(), std::ios::out );
  }

  fStream.precision(std::numeric_limits<realT>::digits10 + 2);

  array<T>& field = iobjectmap->second;

  // write data
  lSet::const_iterator isubset=subset.begin();
  for( size_t i =0; i < subset.size(); ++i, ++isubset){
  	 fStream << field[*isubset];;
  	 fStream << "\n";
  }

  fStream.close();
}





// *********************************************************************************************************************/
/**
 * @author S. Walsh
 * @tparam valueStr  the string value
 * @tparam fieldName the key of the field
 *
 */
template< typename T >
void ObjectDataStructureBaseT::SetFieldToConstantFromString(const std::string& fieldName,
                                                            const std::string& valueStr,
                                                            const bool& additive)
{
  // get the data member map corresponding to the type
  std::map< std::string, array<T> >& objectmap = GetDataMemberMap<array<T> >();

  // check if field exisits
  typename std::map< std::string, array<T> >::iterator iobjectmap = objectmap.find(fieldName);
  if( iobjectmap == objectmap.end() )
    throw GPException("didn't find() requested field "+ fieldName +" in ObjectDataStructureBaseT /n");

  T value = fromString<T>(valueStr);
  if (additive)
  {
    objectmap[fieldName] += value;
  }
  else
  {
    objectmap[fieldName] = value;
  }

}

template< typename T >
void ObjectDataStructureBaseT::SetFieldToConstantFromString(const std::string& fieldName,
                                                            const std::string& valueStr,
                                                            const lSet& subset,
                                                            const bool& additive)
{
  // get the data member map corresponding to the type
  std::map< std::string, array<T> >& objectmap = GetDataMemberMap<array<T> >();

  // check if field exisits
  typename std::map< std::string, array<T> >::iterator iobjectmap = objectmap.find(fieldName);
  if( iobjectmap == objectmap.end() )
    throw GPException("didn't find() requested field "+ fieldName +" in ObjectDataStructureBaseT /n");

  T value = fromString<T>(valueStr);

  if (additive)
  {
    for(lSet::const_iterator si=subset.begin() ; si!=subset.end() ; ++si )
      iobjectmap->second[ *si ] += value;
  }
  else
  {
    for(lSet::const_iterator si=subset.begin() ; si!=subset.end() ; ++si )
      iobjectmap->second[ *si ] = value;
  }

}


// *********************************************************************************************************************

template< typename T >
void ObjectDataStructureBaseT::CopyFieldSubset(const std::string& fieldName,
                                                            const lArray1d& source,
                                                            const lArray1d& target)
{
  // get the data member map corresponding to the type
  std::map< std::string, array<T> >& objectmap = GetDataMemberMap<array<T> >();

  // check if field exisits
  typename std::map< std::string, array<T> >::iterator iobjectmap = objectmap.find(fieldName);
  if( iobjectmap == objectmap.end() )
    throw GPException("didn't find() requested field "+ fieldName +" in ObjectDataStructureBaseT /n");


  if( source.size() != target.size() )
    throw GPException("ObjectDataStructureBaseT::CopyFieldSubset source and target should be the same size. /n");
  
  for(lArray1d::const_iterator si=source.begin(), ti=target.begin() ; si!=source.end() ; ++si,++ti )
    iobjectmap->second[ *ti ] = iobjectmap->second[ *si ];
  
}

// *********************************************************************************************************************



template< typename TYPE >
void ObjectDataStructureBaseT::ConstructVariableOneToManyFromInverse( const std::string& mapname, const array<TYPE>& inverseMap )
{
  array<lSet>& map = this->m_UnorderedVariableOneToManyMaps[mapname];

  map.resize(this->m_DataLengths);
  for( typename array<TYPE>::const_iterator i=inverseMap.begin() ; i!=inverseMap.end() ; ++i )
  {
    for( typename TYPE::const_iterator j=i->begin() ; j!=i->end() ; ++j  )
    {
      map[*j].insert(*i);
    }
  }
}


template< typename TYPE, typename TYPE2 >
void ObjectDataStructureBaseT::AddToVariableOneToManyFromInverse( const std::string& mapname, const array<TYPE>& inverseMap, const TYPE2& newIndices )
{
  array<lSet>& map = this->m_UnorderedVariableOneToManyMaps[mapname];

  for( typename TYPE2::const_iterator i=newIndices.begin() ; i!=newIndices.end() ; ++i )
  {
    for( typename TYPE::const_iterator j=inverseMap[*i].begin() ; j!=inverseMap[*i].end() ; ++j  )
    {
      map[*j].insert(*i);
    }
  }
}



template< typename T >
void ObjectDataStructureBaseT::CopyMemberFields( const localIndex source, const localIndex destination, std::map< std::string, array<T> >& member )
{
  for( typename std::map< std::string, array<T> >::iterator i=member.begin() ; i!=member.end() ; ++i )
  {
    array<T>& field = i->second;
    field[destination] = field[source];

  }
}

template< typename T >
void ObjectDataStructureBaseT::CopyMemberFields( const localIndex source,
                                                 const localIndex destination0,
                                                 const localIndex destination1,
                                                 std::map< std::string, array<T> >& member )
{
  for( typename std::map< std::string, array<T> >::iterator i=member.begin() ; i!=member.end() ; ++i )
  {
    array<T>& field = i->second;
    field[destination0] = field[source];
    field[destination1] = field[source];

  }
}

template< typename T >
void ObjectDataStructureBaseT::CopyMemberFields( const localIndex source, const lArray1d& destination, std::map< std::string, array<T> >& member )
{

  for( typename std::map< std::string, array<T> >::iterator i=member.begin(); i!=member.end() ; ++i )
  {
    array<T>& field = i->second;
    for( lArray1d::const_iterator iterDest=destination.begin() ; iterDest!=destination.end() ; ++iterDest )
    {
      field[*iterDest] = field[source];
    }

  }
}

template< typename T >
void ObjectDataStructureBaseT::AllocateTemporaryFields( const array<string>& names, array<array<T>* >& vars )
{
  vars.resize( names.size() );
  for( array<string>::size_type i=0 ; i<names.size() ; ++i )
  {
    this->AddKeylessDataField<T>(names[i],true,true);
    vars[i] = &(this->GetFieldData<T>(names[i]));
  }
}

template< typename T >
void ObjectDataStructureBaseT::DeallocateTemporaryFields( const array<string>& names )
{
  for( array<string>::size_type i=0 ; i<names.size() ; ++i )
  {
      this->RemoveDataField<T>(names[i]);
  }
}

// *********************************************************************************************************************
/** @name Template Specializations
 */
///@{

template<>
inline std::map< std::string, lArray1d>& ObjectDataStructureBaseT::GetDataMemberMap<lArray1d>()
{ return m_LocalIndexData; }
template<>
inline const std::map< std::string, lArray1d>& ObjectDataStructureBaseT::GetDataMemberMap<lArray1d>() const
{ return m_LocalIndexData; }

template<>
inline std::map< std::string, gArray1d>& ObjectDataStructureBaseT::GetDataMemberMap<gArray1d>()
{ return m_GlobalIndexData; }
template<>
inline const std::map< std::string, gArray1d>& ObjectDataStructureBaseT::GetDataMemberMap<gArray1d>() const
{ return m_GlobalIndexData; }

template<>
inline std::map< std::string, array<int> >& ObjectDataStructureBaseT::GetDataMemberMap<array<integer>>()
{ return m_IntegerData; }
template<>
inline const std::map< std::string, array<int> >& ObjectDataStructureBaseT::GetDataMemberMap<array<integer>>() const
{ return m_IntegerData; }

template<>
inline std::map< std::string, array<realT> >& ObjectDataStructureBaseT::GetDataMemberMap<array<real64>>()
{ return m_realData; }
template<>
inline const std::map< std::string, array<realT> >& ObjectDataStructureBaseT::GetDataMemberMap<array<real64>>() const
{ return m_realData; }

template<>
inline std::map< std::string, array<R1Tensor> >& ObjectDataStructureBaseT::GetDataMemberMap<array<R1Tensor> >()
{ return m_R1TensorData; }
template<>
inline const std::map< std::string, array<R1Tensor> >& ObjectDataStructureBaseT::GetDataMemberMap<array<R1Tensor> >() const
{ return m_R1TensorData; }

template<>
inline std::map< std::string, array<R2Tensor> >& ObjectDataStructureBaseT::GetDataMemberMap<array<R2Tensor> >()
{ return m_R2TensorData; }
template<>
inline const std::map< std::string, array<R2Tensor> >& ObjectDataStructureBaseT::GetDataMemberMap<array<R2Tensor> >() const
{ return m_R2TensorData; }

template<>
inline std::map< std::string, array<R2SymTensor> >& ObjectDataStructureBaseT::GetDataMemberMap<array<R2SymTensor> >()
{ return m_R2SymTensorData; }
template<>
inline const std::map< std::string, array<R2SymTensor> >& ObjectDataStructureBaseT::GetDataMemberMap<array<R2SymTensor> >() const
{ return m_R2SymTensorData; }

template<>
inline std::map< std::string, array<array<R1Tensor> > >& ObjectDataStructureBaseT::GetDataMemberMap<array<array<R1Tensor> > >()
{ return m_ArrayR1TensorData; }
template<>
inline const std::map< std::string, array<array<R1Tensor> > >& ObjectDataStructureBaseT::GetDataMemberMap<array<array<R1Tensor> > >() const
{ return m_ArrayR1TensorData; }



///@}


#endif /* OBJECTDATASTUCTUREBASET_H_ */
