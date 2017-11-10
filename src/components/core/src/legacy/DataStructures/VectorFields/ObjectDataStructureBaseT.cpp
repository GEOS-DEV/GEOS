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
 * @file ObjectDataStructureBaseT.cpp
 * @author Randolph Settgast
 * @date created on Sep 13, 2010
 */


#include "ObjectDataStructureBaseT.h"

#include "../../../codingUtilities/Functions.hpp"
#include "ObjectManagers/FunctionManager.h"
#include "Utilities/Utilities.h"
#include "Utilities/FieldTypeMultiPtr.h"
#include "BoundaryConditions/BoundaryConditions.h"
#include "IO/BinStream.h"

#include "DataStructures/EncapsulatedObjects/EncapsulatedObjectManagerBase.h"

//**********************************************************************************************************************
/**
 * @author Scott Johnson
 * Make sure m_DataLengths is initialized (see Myers Item 4)
 * */
/*
ObjectDataStructureBaseT::ObjectDataStructureBaseT():
m_objectType(),
m_DataLengths(0),
m_isExternal(m_IntegerData[Field<FieldInfo::isExternal>::Name()]),
m_childIndices( m_VariableOneToManyMaps["childIndices"] ),
m_parentIndex( m_OneToOneMaps["parentIndex"] )
{
  this->AddKeyedDataField<FieldInfo::isDomainBoundary>(  );
  this->AddKeyedDataField<FieldInfo::ghostRank>( );
  this->AddKeyedDataField<FieldInfo::ownedByRank>( );


  m_parentIndex = INT_MAX;
}
*/

//**********************************************************************************************************************
/**
 * @author R. Settgast
 * @param objectType the type of object this is as defined in the enum ObjectDataStructureBaseT::ObjectType
 * */
ObjectDataStructureBaseT::ObjectDataStructureBaseT( const ObjectType objectType ):
m_localToGlobalMap(),
m_globalToLocalMap(),
m_maxGlobalNumber(std::numeric_limits<globalIndex>::max()),
m_Sets(),
m_bcData(),
m_objectType(objectType),
m_DataLengths(0),
m_LocalIndexData(),
m_GlobalIndexData(),
m_IntegerData(),
m_realData(),
m_R1TensorData(),
m_R2TensorData(),
m_R2SymTensorData(),
m_OneToOneMaps(),
m_FixedOneToManyMaps(),
m_VariableOneToManyMaps(),
m_UnorderedVariableOneToManyMaps(),
m_isExternal(m_IntegerData[Field<FieldInfo::isExternal>::Name()]),
m_processColor(m_IntegerData[Field<FieldInfo::processColor>::Name()]),
m_processColorPad(),
m_childIndices( m_VariableOneToManyMaps["childIndices"] ),
m_parentIndex( m_OneToOneMaps["parentIndex"] )
{
  this->AddKeyedDataField<FieldInfo::isDomainBoundary>(  );
  this->AddKeyedDataField<FieldInfo::processColor>(  );
  this->AddKeyedDataField<FieldInfo::ghostRank>( );
  this->AddKeyedDataField<FieldInfo::ownedByRank>( );

  array<Field<FieldInfo::isDomainBoundary>::Type>& isDomainBoundary = this->GetFieldData<FieldInfo::isDomainBoundary>();
  array<Field<FieldInfo::ghostRank>::Type>& ghostRank = this->GetFieldData<FieldInfo::ghostRank>();
  array<Field<FieldInfo::ownedByRank>::Type>& isLocallyOwned = this->GetFieldData<FieldInfo::ownedByRank>();

  isDomainBoundary = 0;
  m_isExternal = 0;
  ghostRank = -1;
  isLocallyOwned = 0;
  m_processColor = -1;
  m_processColorPad = 100;

  m_parentIndex = LOCALINDEX_MAX;

  m_parentIndex.SetRelatedObject( this );

  MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);

}


//**********************************************************************************************************************
/**
 *
 * @author R. Settgast
 * @param init object used for initialization
 * @return none
 *
 * Copy constructor
 */
ObjectDataStructureBaseT::ObjectDataStructureBaseT( const ObjectDataStructureBaseT& init ):
m_localToGlobalMap( init.m_localToGlobalMap ),
m_globalToLocalMap( init.m_globalToLocalMap ),
m_maxGlobalNumber( init.m_maxGlobalNumber ),
m_Sets( init.m_Sets ),
m_bcData( init.m_bcData ),
m_objectType(init.m_objectType),
m_DataLengths(init.m_DataLengths ),
m_LocalIndexData( init.m_LocalIndexData ),
m_GlobalIndexData( init.m_GlobalIndexData ),
m_IntegerData( init.m_IntegerData ),
m_realData( init.m_realData ),
m_R1TensorData( init.m_R1TensorData ),
m_R2TensorData( init.m_R2TensorData ),
m_R2SymTensorData( init.m_R2SymTensorData ),
m_OneToOneMaps( init.m_OneToOneMaps ),
m_FixedOneToManyMaps( init.m_FixedOneToManyMaps ),
m_VariableOneToManyMaps( init.m_VariableOneToManyMaps ),
m_UnorderedVariableOneToManyMaps( init.m_UnorderedVariableOneToManyMaps ),
m_isExternal(m_IntegerData[Field<FieldInfo::isExternal>::Name()]),
m_processColor(m_IntegerData[Field<FieldInfo::processColor>::Name()]),
m_processColorPad( init.m_processColorPad ),
m_childIndices( m_VariableOneToManyMaps["childIndices"] ),
m_parentIndex( m_OneToOneMaps["parentIndex"] ),
m_rank(init.m_rank)
{
  m_parentIndex.SetRelatedObject( this );
}


/**
 * @author R. Settgast
 * @return
 */
ObjectDataStructureBaseT::~ObjectDataStructureBaseT()
{
  for( array<BoundaryConditionBase*>::iterator i=m_bcData.begin() ; i!=m_bcData.end() ; ++i )
  {
    delete *i;
  }
}



//********************************************************************************************************************
// Allocation, resize, insertion, etc
//********************************************************************************************************************

/**
 *
 * @author R. Settgast
 * @param size new length of data arrays
 * @return
 *
 * Sets the size of each data member array by calling appropriate resize function on each object.
 */
globalIndex ObjectDataStructureBaseT::resize( const localIndex size, const bool assignGlobals )
{
  const localIndex oldSize = m_DataLengths;
  m_DataLengths = size;

  m_localToGlobalMap.resize(m_DataLengths);
  ResizeObjectField(m_IntegerData);
  ResizeObjectField(m_LocalIndexData);
  ResizeObjectField(m_GlobalIndexData);
  ResizeObjectField(m_realData);
  ResizeObjectField(m_R1TensorData);
  ResizeObjectField(m_ArrayR1TensorData);
  ResizeObjectField(m_R2TensorData);
  ResizeObjectField(m_R2SymTensorData);
  ResizeObjectField(m_OneToOneMaps);
  ResizeObjectField(m_FixedOneToManyMaps);
  ResizeObjectField(m_VariableOneToManyMaps);
  ResizeObjectField(m_VariableOneToManyToManyMaps);
  ResizeObjectField(m_UnorderedVariableOneToManyMaps);



  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  for( localIndex a=oldSize ; a<size ; ++a )
  {
    m_parentIndex[a] = LOCALINDEX_MAX;
    if( assignGlobals )
    {
      m_localToGlobalMap[a] = GlobalIndexManager::Index(rank,a);
      m_globalToLocalMap[ m_localToGlobalMap[a] ] = a;
    }
    else
    {
      m_localToGlobalMap[a] = GLOBALINDEX_MAX;
    }
  }

  globalIndex rval = 0;
  if( size>oldSize )
  {
    rval = m_localToGlobalMap[oldSize];
  }
  return rval;
}

// *********************************************************************************************************************
/**
 * @author R. Settgast
 * @param member data member that will be changing size
 *
 * This function will change the size of a single data array to be of length m_DataLengths.
 * Obviously this is only to be called by the public resize() function.
 */
template< typename T >
void ObjectDataStructureBaseT::ResizeObjectField( std::map< std::string, T>& member )
{
//  const typename T::size_type oldSize = this->m_DataLengths;

  for( typename std::map< std::string, T >::iterator i = member.begin() ;
       i!=member.end() ;
       ++i )
  {
    T& container = (i->second);
    container.resize(m_DataLengths);

//    std::fill( container.begin()+oldSize-1, container.end(), std::numeric_limits<typename T::value_type>::max() );
  }
}
/*
template void ObjectDataStructureBaseT::ResizeObjectField( std::map< std::string, array<integer>>& );
template void ObjectDataStructureBaseT::ResizeObjectField( std::map< std::string, array<real64>>& );
template void ObjectDataStructureBaseT::ResizeObjectField( std::map< std::string, array<R1Tensor> >& );
template void ObjectDataStructureBaseT::ResizeObjectField( std::map< std::string, array<R2Tensor> >& );
template void ObjectDataStructureBaseT::ResizeObjectField( std::map< std::string, array<R2SymTensor> >& );
template void ObjectDataStructureBaseT::ResizeObjectField( std::map< std::string, OneToOneRelationship >& );
template void ObjectDataStructureBaseT::ResizeObjectField( std::map< std::string, FixedOneToManyRelationship >& );
template void ObjectDataStructureBaseT::ResizeObjectField( std::map< std::string, array<lArray1d> >& );
template void ObjectDataStructureBaseT::ResizeObjectField( std::map< std::string, array< lSet > >& );
*/



void ObjectDataStructureBaseT::insert( const localIndex i, const globalIndex gIndex)
{
  // insert global mappings
  {
    gArray1d::iterator iter0 = m_localToGlobalMap.begin() + i;
    iter0 = m_localToGlobalMap.insert(iter0, gIndex);
    localIndex j = i;

    //move all successive local indices ahead by one
    while(iter0 != m_localToGlobalMap.end())
    {
      m_globalToLocalMap[*iter0] = j++;
      ++iter0;
    }
  }

  // call normal insert command
  insertBase(i);
}

globalIndex ObjectDataStructureBaseT::insert( const localIndex i, const bool assignGlobals )
{
  insertBase(i);

  m_localToGlobalMap.resize(m_DataLengths);

  if( assignGlobals )
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    m_localToGlobalMap[i] = GlobalIndexManager::Index(rank,i);
    m_globalToLocalMap[ m_localToGlobalMap[i] ] = i;
  }
  else
  {
    m_localToGlobalMap[i] = GLOBALINDEX_MAX;
  }

  return m_localToGlobalMap[i];
}

void ObjectDataStructureBaseT::insertBase( const localIndex i)
{
  ++m_DataLengths;

  insert( this->m_IntegerData, i );
  insert( this->m_LocalIndexData, i );
  insert( this->m_GlobalIndexData, i );
  insert( this->m_realData, i );
  insert( this->m_R1TensorData, i );
  insert( this->m_R2TensorData, i );
  insert( this->m_R2SymTensorData, i );
  insert( this->m_OneToOneMaps, i );

  insert( this->m_VariableOneToManyMaps, i );
  insert( this->m_FixedOneToManyMaps, i );
  insert( this->m_UnorderedVariableOneToManyMaps, i );
}



template< typename T >
void ObjectDataStructureBaseT::insert( std::map< std::string, array<T> >& member, const localIndex index )
{
  for( typename std::map< std::string, array<T> >::iterator i=member.begin(); i!=member.end() ; ++i )
  {
    array<T>& field = i->second;
    typename array<T>::iterator iter0 = field.begin() + index;
    T zero ;
    zero = 0;
    field.insert(iter0,zero);
  }
}

void ObjectDataStructureBaseT::insert( std::map< std::string, FixedOneToManyRelation >& member, const localIndex index )
{
  for( std::map< std::string, FixedOneToManyRelation >::iterator i=member.begin(); i!=member.end() ; ++i )
  {
    FixedOneToManyRelation& field = i->second;

    //exception handling
    if(index > field.Dimension(0))
      throw GPException("Index out of range for ObjectDataStructureBaseT::insert( std::map< std::string, Array2dT<T> >!");

    //insert dim1 entries at iter0, where iter0 is the entry at index*dim1
    localIndex dim1 = field.Dimension(1);
    FixedOneToManyRelation::iterator iter0 = field.begin() + index * dim1;
    localIndex zero = 0;
    for(localIndex j = 0; j < dim1; ++j)
      field.insert(iter0, zero);
  }
}

template< typename T >
void ObjectDataStructureBaseT::insert( std::map< std::string, InterObjectRelation<array<T> > >& member, const localIndex index )
{
  for( typename std::map< std::string, InterObjectRelation<array<T> > >::iterator i=member.begin(); i!=member.end() ; ++i )
  {
    InterObjectRelation<array<T> >& field = i->second;
    typename InterObjectRelation<array<T> >::iterator iter0 = field.begin() + index;
    field.insert(iter0, T());
  }
}


void ObjectDataStructureBaseT::erase( const localIndex i )
{
  if( i >= m_DataLengths )
    throw GPException("ObjectDataStructureBaseT::erase - index is greater than size");
  else
    --m_DataLengths;

  // delete global mappings
  {
    gArray1d::iterator iter0 = m_localToGlobalMap.begin() + i;
    m_globalToLocalMap.erase(*iter0);
    iter0 = m_localToGlobalMap.erase(iter0);

    //move all successive local indices back by one
    localIndex j = i;
    while(iter0 != m_localToGlobalMap.end())
    {
      m_globalToLocalMap[*iter0] = j++;
      ++iter0;
    }
  }

  // erase, adjust sets.
  for( std::map< std::string, lSet >::iterator set=m_Sets.begin() ; set!=m_Sets.end() ; ++set )
  {
    // TODO
    throw GPException("ObjectDataStructureBaseT::erase not implemented for Sets!\n");
  }

  erase( this->m_IntegerData, i );
  erase( this->m_LocalIndexData, i );
  erase( this->m_GlobalIndexData, i );
  erase( this->m_realData, i );
  erase( this->m_R1TensorData, i );
  erase( this->m_R2TensorData, i );
  erase( this->m_R2SymTensorData, i );
  erase( this->m_OneToOneMaps, i );
  erase( this->m_VariableOneToManyMaps, i );
  erase( this->m_FixedOneToManyMaps, i );
  erase( this->m_UnorderedVariableOneToManyMaps, i );
}


template< typename T >
void ObjectDataStructureBaseT::erase( std::map< std::string, array<T> >& member, const localIndex index )
{
  for( typename std::map< std::string, array<T> >::iterator i=member.begin(); i!=member.end() ; ++i )
  {
    array<T>& field = i->second;
    typename array<T>::iterator iter0 = field.begin() + index;
    field.erase(iter0);
  }

}

template< typename T >
void ObjectDataStructureBaseT::erase( std::map< std::string, InterObjectRelation<T> >& member, const localIndex index )
{
  for( typename std::map< std::string, InterObjectRelation<T> >::iterator i=member.begin(); i!=member.end() ; ++i )
  {
    InterObjectRelation<T>& field = i->second;
    typename InterObjectRelation<T>::iterator iter0 = field.begin() + index;
    field.erase(iter0);
  }
}

template< typename T >
void ObjectDataStructureBaseT::erase( std::map< std::string, Array2dT<T> >& member, const localIndex index )
{
  //note: for each 2D array, the assumption is that you want to eliminate all elements in a "row" (dimension 0) denoted by index
  for( typename std::map< std::string, Array2dT<T> >::iterator i=member.begin(); i!=member.end() ; ++i )
  {
    Array2dT<T>& field = i->second;

    //exception handling
    if(index >= field.Dimension(0))
      throw GPException("Index out of range for ObjectDataStructureBaseT::erase( std::map< std::string, Array2dT<T> >!");

    //insert dim1 entries at iter0, where iter0 is the entry at index*dim1
    localIndex dim1 = field.Dimension(1);
    typename Array2dT<T>::iterator iter0 = field.begin() + index * dim1;
    for(localIndex j = 0; j < dim1; ++j)
      iter0 = field.erase(iter0);
  }
}


//********************************************************************************************************************
// registration, and field access
//********************************************************************************************************************


/**
 * @author R. Settgast
 * @tparam T the type of data being held in the field
 * @param[in]  fieldname the name of the new field
 * @return none
 *
 * This function adds a field to one of the member arrays. Which member is determined
 * by the type specification given by the caller. Thus:
 *
 *  AddDataField<int>("newIntField"); will add a field to m_IntergerData with the key of "newIntField" and
 *
 *  AddDataField<R1Tensor>("newR1TensorField"); will add a field to m_R1TensorData with the key of "newR1TensorField" and
 */
template< typename T >
int ObjectDataStructureBaseT::AddKeylessDataField( const std::string& name, const bool restart, const bool plot )
{
  // get data member map corresponding to the type
  std::map< std::string, array<T> >& DataArrayMember = GetDataMemberMap<array<T> >();

  // the [] will create the new map entry
  DataArrayMember[name];

  // resize the new field
  DataArrayMember[name].resize(m_DataLengths);


  std::map<std::string, FieldBase*>::iterator i = FieldInfo::AttributesByName.find(name);
  if( i == FieldInfo::AttributesByName.end() )
  {
    FieldInfo::AttributesByName[name]  = new FieldBase( FieldInfo::noKey, name, restart, plot );
  }
  else
  {
    FieldInfo::AttributesByName[name]->m_name = name;
    FieldInfo::AttributesByName[name]->m_WriteToRestart = restart;
    FieldInfo::AttributesByName[name]->m_WriteToPlot = plot;
//    throw GPException("ObjectDataStructureBaseT::AddKeylessDataField(): field already exists");
  }

  if( !FieldHasSingleType(name) ){
    throw GPException("ObjectDataStructureBaseT::AddKeylessDataField(): field " + name + " defined with multiple types on the same object.");
  }

  return 0;
}
template int ObjectDataStructureBaseT::AddKeylessDataField<localIndex>( const std::string&, const bool, const bool );
template int ObjectDataStructureBaseT::AddKeylessDataField<globalIndex>( const std::string&, const bool, const bool );
template int ObjectDataStructureBaseT::AddKeylessDataField<int>( const std::string&, const bool, const bool );
template int ObjectDataStructureBaseT::AddKeylessDataField<realT>( const std::string&, const bool, const bool );
template int ObjectDataStructureBaseT::AddKeylessDataField<R1Tensor>( const std::string&, const bool, const bool );
template int ObjectDataStructureBaseT::AddKeylessDataField<R2Tensor>( const std::string&, const bool, const bool );
template int ObjectDataStructureBaseT::AddKeylessDataField<R2SymTensor>( const std::string&, const bool, const bool );
template int ObjectDataStructureBaseT::AddKeylessDataField<array<R1Tensor> >( const std::string&, const bool, const bool );


int ObjectDataStructureBaseT::AddKeylessDataField( FieldType type, const std::string& name, const bool restart, const bool plot )
{
  int rvalue;
  switch(type)
  {
    case FieldInfo::integerField:     rvalue = AddKeylessDataField<int>( name, restart, plot );  break;
    case FieldInfo::localIndexField:  rvalue = AddKeylessDataField<localIndex>( name, restart, plot );  break;
    case FieldInfo::globalIndexField: rvalue = AddKeylessDataField<globalIndex>( name, restart, plot );  break;
    case FieldInfo::realField:        rvalue = AddKeylessDataField<realT>( name, restart, plot );       break;
    case FieldInfo::R1TensorField:    rvalue = AddKeylessDataField<R1Tensor>(name, restart, plot );     break;
    case FieldInfo::R2TensorField:    rvalue = AddKeylessDataField<R2Tensor>( name, restart, plot );    break;
    case FieldInfo::R2SymTensorField: rvalue = AddKeylessDataField<R2SymTensor>( name, restart, plot ); break;
    case FieldInfo::integerParameter:
    case FieldInfo::realParameter:
    case FieldInfo::numFieldTypes:
    default:
      throw GPException("AddKeylessDataField: Unrecognized field type "+ type);
  }
  return rvalue;
}

template<  >
void ObjectDataStructureBaseT::AddMap<OneToOneRelation>( const std::string& name )
{
  OneToOneRelation& map = m_OneToOneMaps[name];
  map.resize( this->m_DataLengths );
}

template<  >
void ObjectDataStructureBaseT::AddMap< OrderedVariableOneToManyRelation >( const std::string& name )
{
  OrderedVariableOneToManyRelation& map = m_VariableOneToManyMaps[name];
  map.resize( this->m_DataLengths );
}

template<  >
void ObjectDataStructureBaseT::AddMap< OrderedVariableOneToManyToManyRelation >( const std::string& name )
{
  OrderedVariableOneToManyToManyRelation& map = m_VariableOneToManyToManyMaps[name];
  map.resize( this->m_DataLengths );
}



template<  >
void ObjectDataStructureBaseT::AddMap< FixedOneToManyRelation >( const std::string& name )
{
  FixedOneToManyRelation& map = m_FixedOneToManyMaps[name];
  map.resize( this->m_DataLengths );
}

template<  >
void ObjectDataStructureBaseT::AddMap< UnorderedVariableOneToManyRelation >( const std::string& name )
{
  UnorderedVariableOneToManyRelation& map = m_UnorderedVariableOneToManyMaps[name];
  map.resize( this->m_DataLengths );
}


template< typename T >
void ObjectDataStructureBaseT::RemoveDataField( const std::string& name )
{
  std::map< std::string, array<T> >& DataArrayMember = GetDataMemberMap<array<T> >();

  // the [] will create the new map entry
  DataArrayMember.erase(name);
}
template void ObjectDataStructureBaseT::RemoveDataField<int>( const std::string& );
template void ObjectDataStructureBaseT::RemoveDataField<realT>( const std::string& );
template void ObjectDataStructureBaseT::RemoveDataField<R1Tensor>( const std::string& );
template void ObjectDataStructureBaseT::RemoveDataField<R2Tensor>( const std::string& );
template void ObjectDataStructureBaseT::RemoveDataField<R2SymTensor>( const std::string& );




//********************************************************************************************************************
// IO
//********************************************************************************************************************

void ObjectDataStructureBaseT::WriteSilo( SiloFile& siloFile,
                                          const std::string& siloDirName,
                                          const std::string& meshname,
                                          const int centering,
                                          const int cycleNum,
                                          const realT problemTime,
                                          const bool isRestart,
                                          const std::string& regionName,
                                          const lArray1d& mask )
{

  std::string subDirectory = siloDirName;
  std::string rootDirectory = "/" + siloDirName;

  {
  std::string shortsubdir(siloDirName);
  std::string::size_type pos = siloDirName.find_last_of("//");

  if( pos != shortsubdir.npos )
  {
    shortsubdir.erase(0,pos+1);
  }

  siloFile.MakeSubDirectory( shortsubdir, rootDirectory );
  DBSetDir(siloFile.m_dbFilePtr, shortsubdir.c_str());
  }



  WriteSilo( siloFile, meshname, centering, cycleNum, problemTime, isRestart, rootDirectory, regionName, mask);

  WriteNonManagedDataMembersToSilo( siloFile, siloDirName, meshname, centering, cycleNum, problemTime, isRestart, rootDirectory, regionName, mask);

  DBSetDir(siloFile.m_dbFilePtr, "..");

}


int ObjectDataStructureBaseT::ReadSilo( const SiloFile& siloFile,
                                        const std::string& siloDirName,
                                        const std::string& meshname,
                                        const int centering,
                                        const int cycleNum,
                                        const realT problemTime,
                                        const bool isRestart,
                                        const std::string& regionName,
                                        const lArray1d& mask )
{
  if( DBSetDir(siloFile.m_dbFilePtr, siloDirName.c_str()) != -1 )
  {

    ObjectDataStructureBaseT::ReadSilo( siloFile, meshname, centering, cycleNum, problemTime, isRestart, regionName, mask);

    ReadNonManagedDataMembersFromSilo( siloFile, siloDirName, meshname, centering, cycleNum, problemTime, isRestart, regionName, mask );

    DBSetDir(siloFile.m_dbFilePtr, "..");

    return 0;
  }
  else
  {
    return 1;
  }
}


void ObjectDataStructureBaseT::WriteSilo( SiloFile& siloFile,
                                          const std::string& meshname,
                                          const int centering,
                                          const int cycleNum,
                                          const realT problemTime,
                                          const bool isRestart,
                                          const std::string& multiRoot,
                                          const std::string& regionName,
                                          const lArray1d& mask ) const
{

  if( isRestart )
  {
	// This data will not be visualized
    siloFile.DBWriteWrapper("m_objectType",static_cast<int>(m_objectType));
    siloFile.DBWriteWrapper("m_DataLengths",m_DataLengths);
    siloFile.DBWriteWrapper("m_localToGlobalMap",m_localToGlobalMap);
    siloFile.DBWriteWrapper("m_globalToLocalMap",m_globalToLocalMap);
    siloFile.DBWriteWrapper("m_maxGlobalNumber",m_maxGlobalNumber);


    siloFile.DBWriteWrapper("m_LocalIndexData",m_LocalIndexData);
    //siloFile.DBWriteWrapper("m_GlobalIndexData",m_GlobalIndexData);
    siloFile.DBWriteWrapper("m_OneToOneMaps",m_OneToOneMaps);
    siloFile.DBWriteWrapper("m_FixedOneToManyMaps",m_FixedOneToManyMaps);
    siloFile.DBWriteWrapper("m_VariableOneToManyMaps",m_VariableOneToManyMaps);
    siloFile.DBWriteWrapper("m_UnorderedVariableOneToManyMaps",m_UnorderedVariableOneToManyMaps);

    array<string> setNames;
    for( std::map<std::string,lSet>::const_iterator i=m_Sets.begin() ; i!=m_Sets.end() ; ++i )
    {
      setNames.push_back( i->first );
    }
    siloFile.DBWriteWrapper("setNames", setNames );
    siloFile.DBWriteWrapper("m_Sets", m_Sets);
  }

  // Data for visualization
  siloFile.WriteFieldMapToSilo<int>(   meshname, m_IntegerData,     centering, cycleNum, problemTime, isRestart, multiRoot, regionName, mask);
  siloFile.WriteFieldMapToSilo<localIndex>(   meshname, m_LocalIndexData,     centering, cycleNum, problemTime, isRestart, multiRoot, regionName, mask);
  siloFile.WriteFieldMapToSilo<globalIndex>(  meshname, m_GlobalIndexData,     centering, cycleNum, problemTime, isRestart, multiRoot, regionName, mask);
  siloFile.WriteFieldMapToSilo<realT>( meshname, m_realData,        centering, cycleNum, problemTime, isRestart,  multiRoot, regionName, mask );
  siloFile.WriteFieldMapToSilo<realT>( meshname, m_R1TensorData,    centering, cycleNum, problemTime, isRestart,  multiRoot, regionName, mask );
  siloFile.WriteFieldMapToSilo<realT>( meshname, m_R2TensorData,    centering, cycleNum, problemTime, isRestart,  multiRoot, regionName, mask );
  siloFile.WriteFieldMapToSilo<realT>( meshname, m_R2SymTensorData, centering, cycleNum, problemTime, isRestart,  multiRoot, regionName, mask );

}



void ObjectDataStructureBaseT::ReadSilo( const SiloFile& siloFile,
                                         const std::string& meshname,
                                         const int centering,
                                         const int cycleNum,
                                         const realT problemTime,
                                         const bool isRestart,
                                         const std::string& regionName,
                                         const lArray1d& mask )
{

  if( isRestart )
  {
	// This data will not be visualized when reading the silo file
    int temp;
    siloFile.DBReadWrapper("m_objectType",temp);
    if( m_objectType != temp )
    {
      throw GPException("ObjectDataStructureBaseT::ReadSilo: m_objectType is wrong");
    }

    unsigned long long dataLengths;
    siloFile.DBReadWrapper("m_DataLengths",dataLengths);
    this->resize( dataLengths );

    siloFile.DBReadWrapper("m_localToGlobalMap",m_localToGlobalMap);
    siloFile.DBReadWrapper("m_globalToLocalMap",m_globalToLocalMap);
    siloFile.DBReadWrapper("m_maxGlobalNumber",m_maxGlobalNumber);


    siloFile.DBReadWrapper("m_LocalIndexData",m_LocalIndexData);
   // siloFile.DBReadWrapper("m_GlobalIndexData",m_GlobalIndexData);

    siloFile.DBReadWrapper("m_OneToOneMaps",m_OneToOneMaps);
    siloFile.DBReadWrapper("m_FixedOneToManyMaps",m_FixedOneToManyMaps);

    std::string junk("m_VariableOneToManyMaps");
    siloFile.DBReadWrapper(junk,m_VariableOneToManyMaps);
//    siloFile.DBReadWrapper("m_VariableOneToManyMaps",m_VariableOneToManyMaps);
    siloFile.DBReadWrapper("m_UnorderedVariableOneToManyMaps",m_UnorderedVariableOneToManyMaps);


    array<string> setNames;
    siloFile.DBReadWrapper("setNames", setNames );

    for( array<string>::const_iterator i=setNames.begin() ; i!=setNames.end() ; ++i )
    {
      m_Sets[*i];
    }


    siloFile.DBReadWrapper("m_Sets", m_Sets);
  }


  if( m_DataLengths>0 )
  {
    // This data can be visualized
    siloFile.ReadFieldMapFromSilo<int>(    m_IntegerData,     meshname, centering, cycleNum, problemTime, isRestart, regionName, mask);
    siloFile.ReadFieldMapFromSilo<localIndex>(    m_LocalIndexData,     meshname, centering, cycleNum, problemTime, isRestart, regionName, mask);
    siloFile.ReadFieldMapFromSilo<globalIndex>(   m_GlobalIndexData,     meshname, centering, cycleNum, problemTime, isRestart, regionName, mask);
    siloFile.ReadFieldMapFromSilo<realT>(  m_realData,        meshname, centering, cycleNum, problemTime, isRestart, regionName, mask );
    siloFile.ReadFieldMapFromSilo<realT>(  m_R1TensorData,    meshname, centering, cycleNum, problemTime, isRestart, regionName, mask );
    siloFile.ReadFieldMapFromSilo<realT>(  m_R2TensorData,    meshname, centering, cycleNum, problemTime, isRestart, regionName, mask );
    siloFile.ReadFieldMapFromSilo<realT>(  m_R2SymTensorData, meshname, centering, cycleNum, problemTime, isRestart, regionName, mask );
  }
}


void ObjectDataStructureBaseT::UncheckAllFieldsForPlot()
{
  UncheckAllFieldsForPlot<int>(m_IntegerData);
  UncheckAllFieldsForPlot<localIndex>(m_LocalIndexData);
  UncheckAllFieldsForPlot<globalIndex>(m_GlobalIndexData);
  UncheckAllFieldsForPlot<realT>(m_realData);
  UncheckAllFieldsForPlot<realT>(m_R1TensorData  );
  UncheckAllFieldsForPlot<realT>(m_R2TensorData  );
  UncheckAllFieldsForPlot<realT>(m_R2SymTensorData);
}


template< typename OUTPUTTYPE, typename T >
void ObjectDataStructureBaseT::UncheckAllFieldsForPlot( const std::map< std::string, T>& member)
{
  for( typename std::map< std::string, T >::const_iterator iter = member.begin(); iter!=member.end(); ++iter )
  {
    const std::string fieldName = iter->first;
    FieldBase*& fieldAttributesMap = stlMapLookup( FieldInfo::AttributesByName, fieldName, "" );
    if (fieldAttributesMap == NULL)
    {
      std::cout << fieldName << " not found in unchecking all fields for plot!";
    }
    else
    {
    fieldAttributesMap-> m_WriteToPlot = false;
    }
  }

}




void ObjectDataStructureBaseT::EncapsulatedObjectsToManagedData_PreRead( const EncapsulatedObjectManagerBase& eom )
{

  array<string> intVarNames;
  array<string> realVarNames;
  array<string> R1TensorVarNames;
  array<string> R2TensorVarNames;
  array<string> R2SymTensorVarNames;

  array<array<integer>*> intVars;
  array<array<real64>*> realVars;
  array<array<R1Tensor>*> R1Vars;
  array<array<R2Tensor>*> R2Vars;
  array<array<R2SymTensor>*> R2SymVars;

  eom.GetVariableNames( intVarNames, realVarNames, R1TensorVarNames, R2TensorVarNames, R2SymTensorVarNames );

  AllocateDummyFields( intVarNames, intVars, false );
  AllocateDummyFields( realVarNames, realVars, false );
  AllocateDummyFields( R1TensorVarNames, R1Vars, false );
  AllocateDummyFields( R2TensorVarNames, R2Vars, false );
  AllocateDummyFields( R2SymTensorVarNames, R2SymVars, false );

  eom.Serialize( intVars, realVars, R1Vars, R2Vars, R2SymVars);

}

void ObjectDataStructureBaseT::EncapsulatedObjectsToManagedData_PostRead( const EncapsulatedObjectManagerBase& eom )
{

  array<string> intVarNames;
  array<string> realVarNames;
  array<string> R1TensorVarNames;
  array<string> R2TensorVarNames;
  array<string> R2SymTensorVarNames;

  eom.GetVariableNames( intVarNames, realVarNames, R1TensorVarNames, R2TensorVarNames, R2SymTensorVarNames );

  DeallocateDummyFields<int>(intVarNames);
  DeallocateDummyFields<realT>(realVarNames);
  DeallocateDummyFields<R1Tensor>(R1TensorVarNames);
  DeallocateDummyFields<R2Tensor>(R2TensorVarNames);
  DeallocateDummyFields<R2SymTensor>(R2SymTensorVarNames);

}

void ObjectDataStructureBaseT::ManagedDataToEncapsulatedObjects_PreRead( const EncapsulatedObjectManagerBase& eom )
{
  array<string> intVarNames;
  array<string> realVarNames;
  array<string> R1TensorVarNames;
  array<string> R2TensorVarNames;
  array<string> R2SymTensorVarNames;

  eom.GetVariableNames( intVarNames, realVarNames, R1TensorVarNames, R2TensorVarNames, R2SymTensorVarNames );

  AllocateDummyFields<int>( intVarNames, false );
  AllocateDummyFields<realT>( realVarNames, false );
  AllocateDummyFields<R1Tensor>( R1TensorVarNames, false );
  AllocateDummyFields<R2Tensor>( R2TensorVarNames, false );
  AllocateDummyFields<R2SymTensor>( R2SymTensorVarNames, false );

}

void ObjectDataStructureBaseT::ManagedDataToEncapsulatedObjects_PostRead( EncapsulatedObjectManagerBase& eom )
{

  array<array<integer>*> intVars;
  array<array<real64>*> realVars;
  array<array<R1Tensor>*> R1Vars;
  array<array<R2Tensor>*> R2Vars;
  array<array<R2SymTensor>*> R2SymVars;

  array<string> intVarNames;
  array<string> realVarNames;
  array<string> R1TensorVarNames;
  array<string> R2TensorVarNames;
  array<string> R2SymTensorVarNames;

  eom.GetVariableNames( intVarNames, realVarNames, R1TensorVarNames, R2TensorVarNames, R2SymTensorVarNames );

  SetDummyFieldPointers( intVarNames, intVars );
  SetDummyFieldPointers( realVarNames, realVars );
  SetDummyFieldPointers( R1TensorVarNames, R1Vars );
  SetDummyFieldPointers( R2TensorVarNames, R2Vars );
  SetDummyFieldPointers( R2SymTensorVarNames, R2SymVars );

  eom.Deserialize( intVars, realVars, R1Vars, R2Vars, R2SymVars );


  DeallocateDummyFields<int>(intVarNames);
  DeallocateDummyFields<realT>(realVarNames);
  DeallocateDummyFields<R1Tensor>(R1TensorVarNames);
  DeallocateDummyFields<R2Tensor>(R2TensorVarNames);
  DeallocateDummyFields<R2SymTensor>(R2SymTensorVarNames);

}


/******************************************************/

void ObjectDataStructureBaseT::WriteAsciiFieldData( const std::vector<FieldType>& types, const array<string>& fieldNames, const std::string& fileName, bool append ){

  std::ofstream fStream;

  if(append){
    fStream.open(fileName.c_str(), std::ios::out | std::ios::app );
  } else {
    fStream.open(fileName.c_str(), std::ios::out );
  }


  fStream.precision(std::numeric_limits<realT>::digits10 + 2);

  unsigned nVars = types.size();
  std::vector<FieldTypeMultiPtr> fieldPtr(nVars );

 if(nVars == 1) {
   // faster if only one field
   WriteAsciiFieldData( types[0], fieldNames[0], fileName,append);
   return;
  }

  int nComponents = 0;   
  for(unsigned j =0; j < nVars; ++j){
    nComponents += FieldSize( types[j] );
    fieldPtr[j].SetFieldPtr(*this, types[j], fieldNames[j]);	
  }
    
  array<real64> x(nComponents );

  for(localIndex i = 0; i < m_DataLengths; ++i){
    std::vector<realT>::iterator xItr = x.begin(); 
    for(unsigned j =0; j < nVars; ++j) xItr = fieldPtr[j].CopyValues(i, xItr);
    xItr = x.begin(); 
    fStream << *xItr ; ++xItr;
    for(; xItr != x.end(); ++xItr) fStream << " " << *xItr; 
    fStream  << "\n";
  }

  fStream.close();

}

void ObjectDataStructureBaseT::WriteAsciiFieldData( const std::vector<FieldType>& types, const array<string>& fieldNames, const std::string& fileName, const lSet& subset, bool append ){

  std::ofstream fStream;

  if(append){
    fStream.open(fileName.c_str(), std::ios::out | std::ios::app );
  } else {
    fStream.open(fileName.c_str(), std::ios::out );
  }

  fStream.precision(std::numeric_limits<realT>::digits10 + 2);

  unsigned nVars = types.size();
  std::vector<FieldTypeMultiPtr> fieldPtr(nVars );

  if(nVars == 1) {
    // faster to use single field version if possible
    WriteAsciiFieldData( types[0], fieldNames[0], fileName, subset, append);
    return;
  }

  int nComponents = 0;   
  for(unsigned j =0; j < nVars; ++j){
    nComponents += FieldSize( types[j] );
    fieldPtr[j].SetFieldPtr(*this, types[j], fieldNames[j]);	
  }
    
  array<real64> x(nComponents );

  for(lSet::const_iterator itr = subset.begin() ; itr != subset.end(); ++itr){
    
    std::vector<realT>::iterator xItr = x.begin(); 
    for(unsigned j =0; j < nVars; ++j) xItr = fieldPtr[j].CopyValues(*itr, xItr);

    xItr = x.begin(); 
    fStream << *xItr ; ++xItr;
    for(; xItr != x.end(); ++xItr) fStream << " " << *xItr; 
    fStream  << "\n";
  }

  fStream.close();

}





//********************************************************************************************************************
// Utilities
//********************************************************************************************************************



//**********************************************************************************************************************
/**
 * @author R. Settgast
 * @param[out] fieldNames vector of field names
 *
 * This function generates a vector of all field names by iterating through each data member and doing a
 * push_back operation on the argument list vector.
**/
void ObjectDataStructureBaseT::GetAllFieldNames( array<string>& fieldNames ) const
{
  fieldNames.clear();
  for( std::map<std::string,array<integer>>::const_iterator i=m_IntegerData.begin() ; i!=m_IntegerData.end() ; ++i )
    fieldNames.push_back(i->first);

  for( std::map<std::string,lArray1d>::const_iterator i=m_LocalIndexData.begin() ; i!=m_LocalIndexData.end() ; ++i )
    fieldNames.push_back(i->first);

  for( std::map<std::string,gArray1d>::const_iterator i=m_GlobalIndexData.begin() ; i!=m_GlobalIndexData.end() ; ++i )
    fieldNames.push_back(i->first);

  for( std::map<std::string,array<real64>>::const_iterator i=m_realData.begin() ; i!=m_realData.end() ; ++i )
    fieldNames.push_back(i->first);

  for( std::map<std::string,array<R1Tensor> >::const_iterator i=m_R1TensorData.begin() ; i!=m_R1TensorData.end() ; ++i )
    fieldNames.push_back(i->first);

  for( std::map<std::string,array<R2Tensor> >::const_iterator i=m_R2TensorData.begin() ; i!=m_R2TensorData.end() ; ++i )
    fieldNames.push_back(i->first);

  for( std::map<std::string,array<R2SymTensor> >::const_iterator i=m_R2SymTensorData.begin() ; i!=m_R2SymTensorData.end() ; ++i )
    fieldNames.push_back(i->first);
}

/**
 * @author R. Settgast
 * @return none
 *
 * This function clears and sets the values of m_globalToLocalMap
 *
 */
void ObjectDataStructureBaseT::ResetGlobalToLocalMap(  )
{
  m_globalToLocalMap.clear();
  for( size_t i=0 ; i<m_localToGlobalMap.size() ; ++i )
  {
    // will not overwrite in the event of duplicate global indices
    m_globalToLocalMap.insert(std::pair<globalIndex,localIndex>( m_localToGlobalMap[i],i) );
  }
}

template< typename T1, typename T2 >
void ObjectDataStructureBaseT::LocalToGlobal( const T1& locals, T2& globals )
{
  for( typename T1::const_iterator a=locals.begin() ; a!=locals.end() ; ++a )
  {
    const globalIndex gIndex = m_localToGlobalMap[*a];
    globals.push_back( gIndex );
  }
}
template void ObjectDataStructureBaseT::LocalToGlobal( const lSet&, gArray1d& );


/**
 * @author R. Settgast
 * @param objectList array to put indexes of boundary objects
 *
 * This function fills the objectList with the localIndex of all boundary objects as
 * specified by the field "isDomainBoundary".
 */
void ObjectDataStructureBaseT::ConstructListOfBoundaryObjects( lArray1d& objectList ) const
{
  const array<integer>& isDomainBoundary = this->GetFieldData<int>("isDomainBoundary");
  for( localIndex k=0 ; k<this->m_DataLengths ; ++k )
  {
    if( isDomainBoundary[k] == 1 )
    {
      objectList.push_back( k );
    }
  }
}

void ObjectDataStructureBaseT::ConstructListOfBoundaryObjects( gArray1d& objectList ) const
{
  const array<integer>& isDomainBoundary = this->GetFieldData<int>("isDomainBoundary");
  for( localIndex k=0 ; k<this->m_DataLengths ; ++k )
  {
    if( isDomainBoundary[k] == 1 )
    {
      objectList.push_back( this->m_localToGlobalMap[k] );
    }
  }
  std::sort( objectList.begin(), objectList.end() );
}

/**
 * @author R. Settgast
 * @param inputSet
 * @param map
 * @param newSetName
 */
void ObjectDataStructureBaseT::ConstructSetFromSetAndMap( const lSet& inputSet,
                                                          const lArray2d& map,
                                                          const std::string& newSetName )
{


  std::map< std::string, lSet >::iterator i = m_Sets.find( newSetName );
  if( i != m_Sets.end() )
  {
    throw GPException("ObjectDataStructureBaseT::ConstructSetFromSetAndMap: Set " + newSetName + " already exists /n");
  }

  lSet& newset = m_Sets[newSetName];
  
  int size = map.Dimension(1);
  for( localIndex ka=0 ; ka<m_DataLengths ; ++ka )
  {
    const localIndex* const sublist = map[ka];
    int addToSet = 0;
    for( int a=0 ; a<size ; ++a )
    {
      if( inputSet.count( sublist[a] ) == 1 )
      {
        ++addToSet;
      }
    }
    if( addToSet == size )
    {
      newset.insert( ka );
    }
  }
}

void ObjectDataStructureBaseT::ConstructSetFromSetAndMap( const lSet& inputSet,
                                                          const array<lArray1d>& map,
                                                          const std::string& newSetName )
{

  std::map< std::string, lSet >::iterator i = m_Sets.find( newSetName );
  if( i != m_Sets.end() )
  {
    throw GPException("ObjectDataStructureBaseT::ConstructSetFromSetAndMap: Set " + newSetName + " already exists /n");
  }

  lSet& newset = m_Sets[newSetName];
  for( localIndex ka=0 ; ka<m_DataLengths ; ++ka )
  {
    int addToSet = 0;
    int size = map[ka].size();
    for( int a=0 ; a<size ; ++a )
    {
      if( inputSet.count( map[ka][a] ) == 1 )
      {
        ++addToSet;
      }
    }
    if( addToSet == size )
    {
      newset.insert( ka );
    }
  }
}

void ObjectDataStructureBaseT::UpdateMaximumGlobalNumber( )
{
//  const globalIndex localMax = m_globalToLocalMap.end()->first;

//  MPI
}


/**
 *
 * @author R. Settgast
 * @param source
 * @param destination
 */
void ObjectDataStructureBaseT::CopyObject( const localIndex source, const localIndex destination )
{

  CopyMemberFields( source, destination, this->m_IntegerData );
  CopyMemberFields( source, destination, this->m_LocalIndexData );
  CopyMemberFields( source, destination, this->m_GlobalIndexData );
  CopyMemberFields( source, destination, this->m_realData );
  CopyMemberFields( source, destination, this->m_R1TensorData );
  CopyMemberFields( source, destination, this->m_R2TensorData );
  CopyMemberFields( source, destination, this->m_R2SymTensorData );

  m_localToGlobalMap[destination] = m_localToGlobalMap[source];
  // TODO m_globalToLocalMap
  // TODO m_maxGlobalNumber

  for( std::map< std::string, lSet >::iterator i=m_Sets.begin() ; i!=m_Sets.end() ; ++i )
  {
    lSet& set = i->second;
    if(set.count(source) > 0) set.insert(destination);
      /*
    for( lSet::iterator j=set.begin() ; j!=set.end() ; ++j )
    {
      if( *j == source )
      {
        set.insert(destination);
        break;
      }
    }*/
  }
}


void ObjectDataStructureBaseT::CopyObjectWithExcludedSets( const localIndex source, const localIndex destination,
                                           const array<string>& excludedSets )
{

  CopyMemberFields( source, destination, this->m_IntegerData );
  CopyMemberFields( source, destination, this->m_LocalIndexData );
  CopyMemberFields( source, destination, this->m_GlobalIndexData );
  CopyMemberFields( source, destination, this->m_realData );
  CopyMemberFields( source, destination, this->m_R1TensorData );
  CopyMemberFields( source, destination, this->m_R2TensorData );
  CopyMemberFields( source, destination, this->m_R2SymTensorData );

  m_localToGlobalMap[destination] = m_localToGlobalMap[source];
  // TODO m_globalToLocalMap
  // TODO m_maxGlobalNumber

  for( std::map< std::string, lSet >::iterator i=m_Sets.begin() ; i!=m_Sets.end() ; ++i )
  {
    const std::string& setName = i->first;
    if( find(excludedSets.begin(),excludedSets.end(),setName) == excludedSets.end() ){ // ie not an excluded set
      lSet& set = i->second;
      if(set.count(source) > 0) set.insert(destination);
    }
  }
}

void ObjectDataStructureBaseT::CopyObjectFields( const localIndex source, const localIndex destination )
{

  CopyMemberFields( source, destination, this->m_IntegerData );
  CopyMemberFields( source, destination, this->m_LocalIndexData );
  CopyMemberFields( source, destination, this->m_GlobalIndexData );
  CopyMemberFields( source, destination, this->m_realData );
  CopyMemberFields( source, destination, this->m_R1TensorData );
  CopyMemberFields( source, destination, this->m_R2TensorData );
  CopyMemberFields( source, destination, this->m_R2SymTensorData );

}

void ObjectDataStructureBaseT::CopyObject( const localIndex source,
                                           const localIndex destination0,
                                           const localIndex destination1 )
{

  CopyMemberFields( source, destination0, destination1, this->m_IntegerData );
  CopyMemberFields( source, destination0, destination1, this->m_LocalIndexData );
  CopyMemberFields( source, destination0, destination1, this->m_GlobalIndexData );
  CopyMemberFields( source, destination0, destination1, this->m_realData );
  CopyMemberFields( source, destination0, destination1, this->m_R1TensorData );
  CopyMemberFields( source, destination0, destination1, this->m_R2TensorData );
  CopyMemberFields( source, destination0, destination1, this->m_R2SymTensorData );
  CopyMemberFields( source, destination0, destination1, this->m_ArrayR1TensorData );

  //m_localToGlobalMap[destination] = m_localToGlobalMap[source];
  // TODO m_globalToLocalMap
  // TODO m_maxGlobalNumber

  for( std::map< std::string, lSet >::iterator i=m_Sets.begin() ; i!=m_Sets.end() ; ++i )
  {
    lSet& set = i->second;
    for( lSet::iterator j=set.begin() ; j!=set.end() ; ++j )
    {
      if( *j == source )
      {
        set.insert(destination0);
        set.insert(destination1);
        break;
      }
    }
  }
}



bool ObjectDataStructureBaseT::SplitObject( const localIndex indexToSplit,
                                            const int rank,
                                            localIndex newIndices[2],
                                            const bool forceSplit )
{
  bool rval=false;

  // if the object index has a zero sized childIndices entry or it is forced to split, then this object can be split into two
  // new objects.
  if( m_childIndices[indexToSplit].size() == 0 || forceSplit)
  {
    // the new indices are tacked on to the end of the arrays
    newIndices[0] = m_DataLengths ;
    newIndices[1] = m_DataLengths + 1;
    this->resize( m_DataLengths + 2 );

    // copy the fields
    CopyObject( indexToSplit, newIndices[0], newIndices[1] );
    m_parentIndex[newIndices[0]] = indexToSplit;
    m_parentIndex[newIndices[1]] = indexToSplit;
    m_childIndices[indexToSplit].push_back( newIndices[0] );
    m_childIndices[indexToSplit].push_back( newIndices[1] );

    m_isExternal[newIndices[0]] = 1;
    m_isExternal[newIndices[1]] = 1;

    rval = true;
  }
  // otherwise this object has already been split, and can't be split again. In this case just set the "new" indices
  // as the existing child indices and return a value of false.
  else
  {
    newIndices[0] = m_childIndices[indexToSplit][0];
    newIndices[1] = m_childIndices[indexToSplit][1];
  }


  const int parentRank = GlobalIndexManager::OwningRank(m_localToGlobalMap[indexToSplit]);

  if( parentRank == rank )
  {
    m_localToGlobalMap[newIndices[0]] = GlobalIndexManager::Index( rank, newIndices[0] );
    m_globalToLocalMap[m_localToGlobalMap[newIndices[0]]] = newIndices[0];
    m_localToGlobalMap[newIndices[1]] = GlobalIndexManager::Index( rank, newIndices[1] );
    m_globalToLocalMap[m_localToGlobalMap[newIndices[1]]] = newIndices[1];
  }
  else
  {
    m_localToGlobalMap[newIndices[0]] = GLOBALINDEX_MAX;
    m_localToGlobalMap[newIndices[1]] = GLOBALINDEX_MAX;
  }

#if 0
  // otherwise this object has already been split, and will now be split into between the original object, and
  // a single new object.
  else
  {
    // the new indices are tacked on to the end of the arrays
    newIndices[0] = m_DataLengths ;
    newIndices[1] = -1;
    this->resize( m_DataLengths + 1 );

    // copy the fields
    CopyObject( indexToSplit, newIndices[0] );
    parentIndex[newIndices[0]] = indexToSplit;
    childIndices[indexToSplit].push_back( newIndices[0] );

    m_isExternal[newIndices[0]] = 1;

    rval = true;

  }
#endif


  return rval;

}


bool ObjectDataStructureBaseT::SplitObject( const localIndex indexToSplit,
                                            const int rank,
                                            localIndex& newIndex)
{
  bool rval=false;

  // if the object index has a zero sized childIndices entry, then this object can be split into two
  // new objects

  // the new indices are tacked on to the end of the arrays
  newIndex = m_DataLengths ;
  this->resize( m_DataLengths + 1 );

  // copy the fields
  CopyObject( indexToSplit, newIndex );
  m_parentIndex[newIndex] = indexToSplit;

  m_childIndices[indexToSplit].push_back( newIndex );

  const int parentRank = GlobalIndexManager::OwningRank(m_localToGlobalMap[indexToSplit]);

  if( parentRank == rank )
  {
    m_localToGlobalMap[newIndex] = GlobalIndexManager::Index( rank, newIndex );
    m_globalToLocalMap[m_localToGlobalMap[newIndex]] = newIndex;
  }
  else
  {
    m_localToGlobalMap[newIndex] = GLOBALINDEX_MAX;
  }


  for( std::map< std::string, lSet >::iterator i=m_Sets.begin() ; i!=m_Sets.end() ; ++i )
  {
    lSet& set = i->second;
    if( set.count( indexToSplit ) > 0 )
    {
      set.insert(newIndex);
    }
  }

  m_isExternal[indexToSplit] = 1;
  m_isExternal[newIndex]     = 1;

  rval = true;


  return rval;

}

void ObjectDataStructureBaseT::SplitObjectExcludeSets( const localIndex indexToSplit,
                                                       const int rank,
                                                       localIndex& newIndex)
{

  // if the object index has a zero sized childIndices entry, then this object can be split into two
  // new objects

  // the new indices are tacked on to the end of the arrays
  newIndex = m_DataLengths ;
  this->resize( m_DataLengths + 1 );

  // copy the fields
  CopyObjectFields( indexToSplit, newIndex );
  m_parentIndex[newIndex] = indexToSplit;

  m_childIndices[indexToSplit].push_back( newIndex );

  const int parentRank = GlobalIndexManager::OwningRank(m_localToGlobalMap[indexToSplit]);

  if( parentRank == rank )
  {
    m_localToGlobalMap[newIndex] = GlobalIndexManager::Index( rank, newIndex );
    m_globalToLocalMap[m_localToGlobalMap[newIndex]] = newIndex;
  }
  else
  {
    m_localToGlobalMap[newIndex] = GLOBALINDEX_MAX;
  }

//  for( std::map< std::string, lSet >::iterator i=m_Sets.begin() ; i!=m_Sets.end() ; ++i )
//  {
//    lSet& set = i->second;
//    if( set.count( indexToSplit ) > 0 )
//    {
//      set.insert(newIndex);
//    }
//  }

  m_isExternal[indexToSplit] = 1;
  m_isExternal[newIndex]     = 1;

}

void ObjectDataStructureBaseT::CreateObject( const int rank,
                                            localIndex& newIndex)
{

  // the new indices are tacked on to the end of the arrays
  newIndex = m_DataLengths ;
  this->resize( m_DataLengths + 1 );
  m_isExternal[newIndex] = 1;

  // copy the fields: More involved than a simple copying of all the fields now
  // For vertices, we will have to interpolate the fields. For now, leave these empty
  // CopyObject( indexToSplit, newIndices[0], newIndices[1] );

  // New guy so has no parent
  //  m_parentIndex[newIndices[0]] = indexToSplit;
  //  m_parentIndex[newIndices[1]] = indexToSplit;
  //  m_childIndices[indexToSplit].push_back( newIndices[0] );
  //  m_childIndices[indexToSplit].push_back( newIndices[1] );

  /* Speak with Randy about this:
  const int parentRank = GlobalIndexManager::OwningRank(m_localToGlobalMap[indexToSplit]);

  if( parentRank == rank )
  {
    m_localToGlobalMap[newIndices[0]] = GlobalIndexManager::Index( rank, newIndices[0] );
    m_globalToLocalMap[m_localToGlobalMap[newIndices[0]]] = newIndices[0];
    m_localToGlobalMap[newIndices[1]] = GlobalIndexManager::Index( rank, newIndices[1] );
    m_globalToLocalMap[m_localToGlobalMap[newIndices[1]]] = newIndices[1];
  }
  else
  {
    m_localToGlobalMap[newIndices[0]] = GLOBALINDEX_MAX;
    m_localToGlobalMap[newIndices[1]] = GLOBALINDEX_MAX;
  }
  */

}


//**********************************************************************************************************************

 /*
   * @author S. Walsh
   * 
   * Use a function to set the values of a field on the object. The function 
   * must be defined in terms of variables given in fields on the same object. 
   * 
   * @param fieldType, the type of the field
   * @param fieldName, the name of the field
   * @param variables, the function variable names 
   * @param variable_types, the function variable types
   * @param component, the tensor component to set (ignored for scalar values)
   * 
   * 
   */

void ObjectDataStructureBaseT::SetFieldEqualToFunction(FieldType fieldType, const std::string& fieldName, 
                             const std::string& functionName,
                             const array<string>& variables, const array<FieldType>& variable_types,
                             int component, realT time, realT dt)
{
    FunctionManager& functionManager = FunctionManager::Instance();   
    Function& func = functionManager.GetFunction(functionName);
    int nVars = variables.size();
    std::vector<FieldTypeMultiPtr> fieldPtr(nVars );
    
    FieldTypeMultiPtr theFieldPtr;
    theFieldPtr.SetFieldPtr(*this, fieldType, fieldName);
  
    int xLength = 0;   
    for(int j =0; j < nVars; ++j){
      if(streq(variables[j],"time")){
        xLength += 1;
        fieldPtr[j].SetFieldPtr(&time);
      } else if(streq(variables[j],"dt")){
        xLength += 1;
        fieldPtr[j].SetFieldPtr(&dt);
      } else {
        xLength += FieldSize( variable_types[j] );
        fieldPtr[j].SetFieldPtr(*this, variable_types[j], variables[j]);
      }
    }
    
    std::vector<realT> x(xLength );
	
    for(localIndex i = 0; i < this->m_DataLengths; ++i){
  	  // pack field values into x
      std::vector<realT>::iterator xItr = x.begin(); 
      for(int j =0; j < nVars; ++j) xItr = fieldPtr[j].CopyValues(i, xItr);
      
      // calculate value
      theFieldPtr.SetValue(i,component, func(x[0]) );
    }
}

void ObjectDataStructureBaseT::SetFieldEqualToFunction(FieldType fieldType, const std::string& fieldName, 
                                                       const std::string& functionName,
                                                       const array<string>& variables, const array<FieldType>& variable_types,
                                                       const lSet& subset,
                                                       int component, realT time , realT dt )
{
  FunctionManager& functionManager = FunctionManager::Instance();
  Function& func = functionManager.GetFunction(functionName);
  int nVars = variables.size();
  std::vector<FieldTypeMultiPtr> fieldPtr(nVars );
    
  FieldTypeMultiPtr theFieldPtr;
  theFieldPtr.SetFieldPtr(*this, fieldType, fieldName);
  
  int xLength = 0;
  for(int j =0; j < nVars; ++j){
    if(streq(variables[j],"time")){
      xLength += 1;
      fieldPtr[j].SetFieldPtr(&time);
    } else if(streq(variables[j],"dt")){
      xLength += 1;
      fieldPtr[j].SetFieldPtr(&dt);
    } else {
      xLength += FieldSize( variable_types[j] );
      fieldPtr[j].SetFieldPtr(*this, variable_types[j], variables[j]);
    }
  }

  std::vector<realT> x(xLength );

	for( lSet::const_iterator si=subset.begin() ; si!=subset.end() ; ++si )
	{
    // pack field values into x
    std::vector<realT>::iterator xItr = x.begin();
    for(int j =0; j < nVars; ++j) xItr = fieldPtr[j].CopyValues(*si, xItr);
    
    // calculate value
    theFieldPtr.SetValue(*si,component, func(x[0]) );
  }
}



template< typename T_indices >
unsigned int ObjectDataStructureBaseT::PackSets( const T_indices& sendlist,
                                                 bufvector& buffer ) const
{
  unsigned int sizeOfPackedChars = 0;


  sizeOfPackedChars += buffer.Pack(m_Sets.size());

  // loop over all sets
  for( std::map<std::string,lSet>::const_iterator iset=m_Sets.begin() ; iset!=m_Sets.end() ; ++iset )
  {
    const std::string& setname = iset->first;
    const lSet& set = iset->second;

    // find how many set members are in the sendlist
    gSet sendSet;
    for( typename T_indices::const_iterator a=sendlist.begin() ; a!=sendlist.end() ; ++a )
    {
      // if the entry from the sendlist is found in the set
//      if( std::find(set.begin(), set.end(), *a ) != set.end() )
      if( set.count(*a) != 0 )
      {
        sendSet.insert( this->m_localToGlobalMap(*a) );
      }
    }

    sizeOfPackedChars += buffer.Pack(setname);
    sizeOfPackedChars += buffer.Pack(sendSet.size());
    for( gSet::const_iterator b=sendSet.begin() ; b!=sendSet.end() ; ++b )
    {
      sizeOfPackedChars += buffer.Pack(*b);
    }
  }

  return sizeOfPackedChars;
}


unsigned int ObjectDataStructureBaseT::UnpackSets( const char*& buffer )
{
  unsigned int sizeOfUnpackedChars = 0;

  lSet::size_type numSets;
  sizeOfUnpackedChars += bufvector::Unpack( buffer, numSets );

  for( lSet::size_type k=0 ; k<numSets ; ++k )
  {
    std::string setname;
    gSet::size_type numEntries;


    sizeOfUnpackedChars += bufvector::Unpack(buffer, setname);
    sizeOfUnpackedChars += bufvector::Unpack(buffer, numEntries);
    for( gSet::size_type a=0 ; a<numEntries ; ++a )
    {
      globalIndex globalNumber;
      sizeOfUnpackedChars += bufvector::Unpack(buffer, globalNumber);

      std::map<globalIndex,localIndex>::iterator iterG2L = m_globalToLocalMap.find(globalNumber);
      if( iterG2L != m_globalToLocalMap.end() )
      {
        localIndex localID = iterG2L->second;
        m_Sets[setname].insert(localID);
      }
    }


  }

  return sizeOfUnpackedChars;

}


template< typename T_indices >
unsigned int ObjectDataStructureBaseT::PackBaseObjectData( bufvector& buffer,
                                                           const T_indices& indices,
                                                           const bool packFields,
                                                           const bool packMaps,
                                                           const bool packSets,
                                                           const bool packRelationsToGlobal ) const
{
  unsigned int sizeOfPackedChars = 0;

  int size = indices.size();
  sizeOfPackedChars += buffer.Pack( size );

  for( typename T_indices::const_iterator index=indices.begin() ; index!=indices.end() ; ++index )
  {

    sizeOfPackedChars += buffer.Pack( m_localToGlobalMap[*index] );

    const localIndex parentIndex = GetParentIndex( *index );
    sizeOfPackedChars += buffer.Pack(m_localToGlobalMap[parentIndex]);

    sizeOfPackedChars += buffer.Pack( m_childIndices[*index].size() );
    for( lArray1d::const_iterator i=m_childIndices[*index].begin() ; i!=m_childIndices[*index].end() ; ++i )
    {
      sizeOfPackedChars += buffer.Pack( m_localToGlobalMap[*i] );
    }
  }

  if( packFields )
  {
    sizeOfPackedChars += this->PackAllFieldsIntoBuffer( buffer, indices );
  }

  if( packMaps )
  {
    sizeOfPackedChars += this->PackRelationsIntoBuffer( buffer, indices, packRelationsToGlobal );
  }

  if( packSets )
  {
    sizeOfPackedChars += this->PackSets( indices, buffer );
  }

return sizeOfPackedChars;
}
template unsigned int ObjectDataStructureBaseT::PackBaseObjectData( bufvector&, const lSet&, const bool, const bool, const bool, const bool ) const;
template unsigned int ObjectDataStructureBaseT::PackBaseObjectData( bufvector&, const lArray1d&, const bool, const bool, const bool, const bool ) const;


unsigned int ObjectDataStructureBaseT::UnpackBaseObjectData( const char*& buffer,
                                                             lArray1d& unpackedLocalIndices,
                                                             lArray1d& newLocalIndices,
                                                             const bool unpackFields,
                                                             const bool unpackMaps,
                                                             const bool unpackSets,
                                                             const bool unpackRelationsToLocal  )
{
  unsigned int sizeOfUnpacked = 0;
  int rank=0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank );


  // unpack the number of received objects
  int numReceivedObjects;
  sizeOfUnpacked += bufvector::Unpack( buffer, numReceivedObjects );
  const localIndex oldSize = this->m_DataLengths;

  // resize the new local indices array
  unpackedLocalIndices.resize( numReceivedObjects );
  // arrays to provide temporary storage to the unpacked data
  gArray1d            globalIndices(numReceivedObjects);
  gArray1d            parentIndices(numReceivedObjects);
  array<gArray1d>  childIndices(numReceivedObjects);

  gArray1d newGlobalIndices;


  // unpack the object information
  int numNewObjects = 0;
  for( int a=0 ; a<numReceivedObjects ; ++a )
  {
    // upack the global index, and the parent index
    sizeOfUnpacked += bufvector::Unpack( buffer, globalIndices[a] );
    sizeOfUnpacked += bufvector::Unpack( buffer, parentIndices[a] );

    // unpack the child indices
    array<lArray1d>::size_type childSize;
    sizeOfUnpacked += bufvector::Unpack( buffer, childSize );
    childIndices[a].resize(childSize);
    for( array<lArray1d>::size_type i=0 ; i<childSize ; ++i )
    {
      sizeOfUnpacked += bufvector::Unpack( buffer, childIndices[a][i] );
    }


    // check to see if the object already exists by checking for the global index in m_globalToLocalMap. If it doesn't,
    // then add the object

    std::map<globalIndex,localIndex>::iterator iterG2L = m_globalToLocalMap.find(globalIndices[a]);
    if( iterG2L == m_globalToLocalMap.end() )
    {
      // object does not exist on this domain
      const localIndex newLocalIndex = oldSize + numNewObjects;

      if( globalIndices[a]==GLOBALINDEX_MAX )
      {
        globalIndices[a] = GlobalIndexManager::Index( rank, newLocalIndex );
      }

      // add the global index of the new object to the globalToLocal map
      m_globalToLocalMap[globalIndices[a]] = newLocalIndex;

      unpackedLocalIndices(a) = newLocalIndex;

      newGlobalIndices.push_back( globalIndices[a] );

      ++numNewObjects;
    }
    else
    {
      // object already exists on this domain

      // get the local index of the node
      localIndex b = iterG2L->second;
      unpackedLocalIndices(a) = b;
    }
  }

  // figure out new size of object container, and resize object
  const localIndex newSize = oldSize + numNewObjects;
  this->resize( newSize );

  // now distribute unpacked data to appropriate locations

  newLocalIndices.resize( numNewObjects );

  // add the new objects to the localToGlobalMap.
  for( int a=0 ; a<numNewObjects ; ++a )
  {
    localIndex b = oldSize + a;
    m_localToGlobalMap[b] = newGlobalIndices(a);
    newLocalIndices[a] = b;
  }


  // now take care of the parent/child data
  {
  lArray1d localChildIndices;
  localChildIndices.reserve(4);
  for( int a=0 ; a<numReceivedObjects ; ++a )
  {
    this->m_parentIndex[unpackedLocalIndices[a]] = stlMapLookup( m_globalToLocalMap, parentIndices[a] );

    if( this->m_parentIndex[unpackedLocalIndices[a]] == unpackedLocalIndices[a] )
    {
      this->m_parentIndex[unpackedLocalIndices[a]] = LOCALINDEX_MAX;
    }

    localChildIndices.resize( childIndices[a].size() );
    for( lArray1d::size_type i=0 ; i<childIndices[a].size() ; ++i )
    {
      if( childIndices[a][i] != LOCALINDEX_MAX )
      {
        localChildIndices[i] = stlMapLookup( m_globalToLocalMap, childIndices[a][i] );
      }
      else
      {
        localChildIndices[i] = LOCALINDEX_MAX;
      }
    }
    this->m_childIndices[unpackedLocalIndices[a]] = localChildIndices;
  }
  }



  if( unpackFields )
  {
    sizeOfUnpacked += UnpackAllFieldsFromBuffer( buffer, unpackedLocalIndices );
  }
  if( unpackMaps )
  {
    sizeOfUnpacked += UnpackRelationsFromBuffer( buffer, unpackedLocalIndices, unpackRelationsToLocal );
  }
  if( unpackSets )
  {
    sizeOfUnpacked += UnpackSets( buffer );
  }


  // have to modify the ghostRank field manually, as the data comes across in reference to the sending
  // rank.
  array<Field<FieldInfo::ghostRank>::Type>& ghostRank = this->GetFieldData<FieldInfo::ghostRank>();

  for( int a=0 ; a<numReceivedObjects ; ++a )
  {
//    ghostRank[ unpackedLocalIndices[a] ] = ghostRank[ GetParentIndex(unpackedLocalIndices[a]) ];
    const int owningRank = GlobalIndexManager::OwningRank( m_localToGlobalMap[unpackedLocalIndices[a]]);

    if( owningRank == rank )
    {
      ghostRank[ unpackedLocalIndices[a] ] = -1;
    }
    else
    {
      ghostRank[ unpackedLocalIndices[a] ] = owningRank;

    }
  }


  return sizeOfUnpacked;
}





template< typename T, typename T_indices >
unsigned int ObjectDataStructureBaseT::PackRelationIntoBuffer( bufvector& buffer,
                                                               const std::map<std::string,T>& map,
                                                               const T_indices& localIndices,
                                                               const bool packRelationsToGlobal ) const
{
  unsigned int sizeOfPackedChars = 0;


  for( typename std::map<std::string,T>::const_iterator i=map.begin() ; i!=map.end() ; ++i )
  {
    const T& relation = i->second;
    sizeOfPackedChars += buffer.Pack( relation, localIndices, packRelationsToGlobal );
  }

  return sizeOfPackedChars;
}





template< typename T >
unsigned int ObjectDataStructureBaseT::UnpackRelationFromBuffer( const char*& buffer,
                                                                 std::map<std::string,T>& map,
                                                                 const lArray1d& localIndices,
                                                                 const bool unpackRelationsToLocal )
{
  unsigned int sizeOfUnpackedChars = 0;

  for( typename std::map<std::string,T>::iterator i=map.begin() ; i!=map.end() ; ++i )
  {
    T& relation = i->second;

    sizeOfUnpackedChars += bufvector::Unpack( buffer, relation, localIndices, unpackRelationsToLocal );
  }
  return sizeOfUnpackedChars;
}

template< typename T_indices >
unsigned int ObjectDataStructureBaseT::PackRelationsIntoBuffer( bufvector& ,
                                                                const T_indices& ,
                                                                const bool  ) const
{
  unsigned int sizeOfPackedChars = 0;

//  sizeOfPackedChars += PackRelationIntoBuffer( buffer, m_OneToOneMaps, localIndices, packRelationsToGlobal );
//  sizeOfPackedChars += PackRelationIntoBuffer( buffer, m_FixedOneToManyMaps, localIndices, packRelationsToGlobal );
//  sizeOfPackedChars += PackRelationIntoBuffer( buffer, m_VariableOneToManyMaps, localIndices, packRelationsToGlobal );
//  sizeOfPackedChars += PackRelationIntoBuffer( buffer, m_UnorderedVariableOneToManyMaps, localIndices, packRelationsToGlobal );

  return sizeOfPackedChars;
}
template unsigned int ObjectDataStructureBaseT::PackRelationsIntoBuffer( bufvector&, const lArray1d&, bool ) const;
template unsigned int ObjectDataStructureBaseT::PackRelationsIntoBuffer( bufvector&, const lSet&, bool ) const;

unsigned int ObjectDataStructureBaseT::UnpackRelationsFromBuffer( const char*& ,
                                                                  const lArray1d& ,
                                                                  const bool  )
{
  unsigned int sizeOfUnpackedChars = 0;

//  sizeOfUnpackedChars += UnpackRelationFromBuffer( buffer, m_OneToOneMaps, localIndices, unpackRelationsToLocal );
//  sizeOfUnpackedChars += UnpackRelationFromBuffer( buffer, m_FixedOneToManyMaps, localIndices, unpackRelationsToLocal );
//  sizeOfUnpackedChars += UnpackRelationFromBuffer( buffer, m_VariableOneToManyMaps, localIndices, unpackRelationsToLocal );
//  sizeOfUnpackedChars += UnpackRelationFromBuffer( buffer, m_UnorderedVariableOneToManyMaps, localIndices, unpackRelationsToLocal );

  return sizeOfUnpackedChars;
}

