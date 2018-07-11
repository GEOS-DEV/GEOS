/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file SiloFile.hpp
 */

#ifndef SILOFILE_HPP_
#define SILOFILE_HPP_

#include "common/DataTypes.hpp"
#include "silo.h"
#include <vector>

#if USE_MPI
#include <mpi.h>
#endif

#include "mpi.h"
#include "pmpio.h"
#include "common/Logger.hpp"

/////////////////////////////////////////////////


#include "mesh/InterObjectRelation.hpp"



/// typedef my callback function for use with WriteMultiXXXX. This function is
// essentially a placeholder
/// for DB,PutMultimesh, DBPutMultivar, etc.
typedef struct _PMPIO_baton_t PMPIO_baton_t;
//typedef int (*DBPutMultimeshType)(DBfile *, char DB_CONSTARR1, int, char
// DB_CONSTARR2, int DB_CONSTARR1, DBoptlist const *);
//typedef int (*DBPutMultivarType)(DBfile *, const char *, int, char **, int *,
// DBoptlist *);


namespace geosx
{

class DomainPartition;
class MeshLevel;
class ElementRegionManager;
namespace constitutive
{
class ConstitutiveManager;
}

// *********************************************************************************************************************
// *********************************************************************************************************************
/**
 * @author settgast
 * This class serves as a wrapper to isolate the code from the specifics of SILO
 * output/input. Its members contain all
 * the necessary information for reading/writing a group of SILO files.
 */
class SiloFile
{

public:

  /// Default Constructor
  SiloFile();

  /// Destructor
  virtual ~SiloFile();

  void MakeSiloDirectories();

  /// Initializes the silo library
  void Initialize( const PMPIO_iomode_t readwrite );

  /// finishes up the silo library usage
  void Finish();

  /// Wait for the Baton when doing PMPIO
  void WaitForBatonWrite( int const domainNumber, int const cycleNum, bool const isRestart );

  void WaitForBaton( int const domainNumber, string const & restartFileName );

  /// Hand off the Baton when doing PMPIO
  void HandOffBaton();


  void MakeSubDirectory( string const & subdir, string const & rootdir )
  {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if( rank == 0 )
    {
      DBMkDir(m_dbBaseFilePtr, rootdir.c_str());
    }
    DBMkDir(m_dbFilePtr, subdir.c_str());
  }

  /// Write out a single mesh object
  void WriteMeshObject(string const & meshName,
                       const localIndex nnodes,
                       real64* coords[3],
                       const globalIndex*,
                       int const numRegions,
                       int const * shapecnt,
                       const localIndex* const * const meshConnectivity,
                       const globalIndex* const * const globalElementNum,
                       int const * const * const,
                       int const * const shapetype,
                       int const * const shapesize,
                       int const cycleNumber,
                       real64 const problemTime);

  /// Write out a polygon mesh object
  void WritePolygonMeshObject(string const & meshName,
                              const localIndex nnodes,
                              real64* coords[3],
                              const globalIndex*,
                              int const numRegions,
                              int const * shapecnt,
                              const localIndex* const * const meshConnectivity,
                              const globalIndex* const * const globalElementNum,
                              int const * const * const,
                              int const * const shapetype,
                              int const * const shapesize,
                              int const cycleNumber,
                              real64 const problemTime,
                              int const lnodelist);

  void WriteDomainPartition( DomainPartition const & domain,
                             int const cycleNum,
                             real64 const problemTime,
                             bool const isRestart );

  void WriteMeshLevel( MeshLevel const * const meshLevel,
                       constitutive::ConstitutiveManager const * const constitutiveManager,
                       int const cycleNum,
                       real64 const problemTime,
                       bool const isRestart );

  void TestPolyhedralCells();
  void XFEMMesh(std::vector<double> xcoords,
                std::vector<double> ycoords,
                std::vector<double> zcoords,
                std::vector<int> nodelist,
                int lnodelist,
                std::vector<int> shapesize,
                std::vector<int> shapecounts,
                std::vector<int> shapetype,
                int nshapetypes,
                int nnodes,
                int nzones,
                int ndims,
                int const cycleNumber,
                real64 const problemTime);

  void WriteDiscreteElementCSGObject(string const & meshName,
                                     const array<R1Tensor>& x,
                                     const array<R1Tensor>& r,
                                     const array<R2Tensor>& R,
                                     int const cycleNumber,
                                     real64 const problemTime);

  void WritePointMesh( string const & meshName,
                       const localIndex numPoints,
                       real64* coords[3],
                       int const cycleNumber,
                       real64 const problemTime );

  void WriteBeamMesh(string const & meshName,
                     const localIndex nnodes,
                     real64* coords[3],
                     const localIndex_array& node1,
                     const localIndex_array& node2,
                     int const cycleNumber,
                     real64 const problemTime);

  void WriteBeamMesh(string const & meshName,
                     const localIndex nnodes,
                     real64* coords[3],
                     const std::map<int, int>& connectivity,
                     int const cycleNumber,
                     real64 const problemTime);

  void WriteBeamMesh(string const & meshName,
                     const localIndex nnodes,
                     real64* coords[3],
                     integer_array& nodelist,
                     int const cycleNumber,
                     real64 const problemTime);

  /**
   * @brief Write out a single "quad" (i.e. hex) mesh object, and record the
   * mesh dimensions for subsequent write data field calls.
   * @author walsh24
   */
  int  WriteQuadMeshObject(string const & meshName,
                           const localIndex nX,
                           const localIndex nY,
                           const localIndex nZ,
                           real64* coords[3],
                           int const cycleNumber,
                           real64 const problemTime);

//  void TestWriteDiscreteElementMeshObject();

  void WriteRegionSpecifications( ElementRegionManager const * const elementManager,
                                  constitutive::ConstitutiveManager const * const constitutiveManager,
                                  string const & meshName,
                                  int const cycleNumber,
                                  real64 const problemTime);


  void WriteManagedGroupSilo( dataRepository::ManagedGroup const * group,
                              string const & siloDirName,
                              string const & meshname,
                              int const centering,
                              int const cycleNum,
                              real64 const problemTime,
                              bool const isRestart,
                              string const & materialName,
                              array<localIndex> const & zoneToMatMap,
                              const localIndex_array& mask );


  /// writes out fields in a data member map
  template< typename OUTPUTTYPE >
  void WriteViewWrappersToSilo( string const & meshname,
                                const dataRepository::ManagedGroup::viewWrapperMap & wrappers,
                                int const centering,
                                int const cycleNum,
                                real64 const problemTime,
                                bool const isRestart,
                                string const & multiRoot,
                                string const & materialName,
                                array<localIndex> const & zoneToMatMap,
                                const localIndex_array& mask );

  template< typename INPUTTYPE, typename TYPE >
  void ReadFieldMapFromSilo( std::map< string, array<TYPE> >& member,
                             string const & meshname,
                             int const centering,
                             int const cycleNum,
                             real64 const problemTime,
                             bool const isRestart,
                             string const & regionName,
                             const localIndex_array& mask ) const;

  /// Write out a data field
  template<typename OUTTYPE, typename TYPE>
  void WriteDataField( string const & meshName,
                       string const & fieldName,
                       const TYPE& field,
                       int const centering,
                       int const cycleNumber,
                       real64 const problemTime,
                       string const & multiRoot,
                       string const & materialName );

  template<typename OUTTYPE, typename TYPE>
  void WriteDataField( string const & meshName,
                       string const & fieldName,
                       const array<TYPE>& field,
                       int const centering,
                       int const cycleNumber,
                       real64 const problemTime,
                       string const & multiRoot,
                       string const & materialName,
                       array<localIndex> const & zoneToMatMap );

  /// Write out a data field
  template< typename INPUTTYPE, typename TYPE>
  void ReadDataField( array<TYPE>& field,
                      string const & meshName,
                      string const & fieldName,
                      int const centering,
                      int const cycleNumber,
                      real64 const problemTime,
                      string const & regionName ) const;

  int GetMeshType( string const & meshName ) const
  {
    int meshType = -1;
    {
      // in order to get mesh type, we might have to go up a few levels in the
      // silo directory structure
      // before we can find the mesh.
      char pwd[256];
      DBGetDir(m_dbFilePtr, pwd);
      for( int i=0 ; i<3 ; ++i )
      {
        meshType = DBInqMeshtype(m_dbFilePtr,meshName.c_str());
        if( meshType != -1 && meshType != 610 )
          break;
        else
          DBSetDir(m_dbFilePtr,"..");
      }
      DBSetDir(m_dbFilePtr,pwd);
    }
    return meshType;
  }

  template <typename TYPE>
  void** GetDataVar( string const & fieldName,
                     string const & meshName,
                     const typename array<TYPE>::size_type nels,
                     int const centering,
                     int const cycleNumber,
                     real64 const problemTime,
                     string const & ) const;

  /// Write out a multi-mesh object
  template< typename CBF >
  void WriteMultiXXXX(const DBObjectType type, CBF DBPutMultiCB,
                      int const centering, string const name, int const cycleNumber,
                      string const & multiRoot, const DBoptlist* optlist = nullptr);


  /// dummy function to get rid of compiler warnings about unused functions from PMPIO's use of static function
  /// definitions in the header. This should go away as we fully utilize PMPIO.
  void StopSiloCompilerWarnings();

  void ClearEmptiesFromMultiObjects(int const cycleNum);

  //private:

  /// pointer to the DBfile that this class is working on
  DBfile* m_dbFilePtr;

  /// pointer to the Master DBfile
  DBfile* m_dbBaseFilePtr;

  /// total number of "parallel" files to write out
  int m_numGroups;

  /// the pmpio baton. A processor needs this to write to the file.
  PMPIO_baton_t *m_baton;

  /// which database to use. DB_PDB or DB_HDF5
  int m_driver;

  /// root of the filename that we will be reading/writing
  string m_plotFileRoot;

  string m_restartFileRoot;

  string const m_siloDirectory = "siloFiles";

  string const m_siloDataSubDirectory = "data";

  string m_fileName;

  string m_baseFileName;

  bool m_markGhosts;

  array<string> m_emptyMeshes;
  array<string> m_emptyVariables;

  integer_array SiloNodeOrdering();


private:
  int m_quadMeshDims[3];
  int m_quadMeshNDims;

};

// *********************************************************************************************************************
// *********************************************************************************************************************
/**
 * Namespace to hold some utilities needed by the functions in SiloFile.
 */
namespace SiloFileUtilities
{
/**
 * @author settgast
 * @tparam OUTTYPE the type of data to write out (int,float,real64)
 * @return the integer identifier associated with the enum DBdatatype defined in
 * the silo.h file.
 *
 * This templated function is a "specialization only" definition. There is no
 * general
 * definition, only specializations for predetermined data types.
 */
template<typename OUTTYPE>
int DB_TYPE();

/**
 * @author settgast
 * @tparam TYPE the data type in question
 * @return the number of "variables" in a data field. For instance, 1 for a int,
 * float or real64.
 * 3 for a R1Tensor...etc.
 */
template<typename TYPE>
int GetNumberOfVariablesInField();

template<typename TYPE>
int GetTensorRank();

template<typename OUTTYPE, typename TYPE>
OUTTYPE CastField(const TYPE& field, int const i = 0);     // avoids compiler
                                                           // warning

template<typename OUTTYPE, typename TYPE>
OUTTYPE CastField(const TYPE& field, int const i)
{
  return static_cast<OUTTYPE>(field.Data()[i]);
//    return field;
}

template<> inline int CastField<int, int> (const int& field, int const )
{
  return field;
}


//  template<> inline int CastField<int, localIndex> (const localIndex& field,
// int const )
//  {
//    return static_cast<int>(field);
//  }

//  template<> inline localIndex CastField<localIndex, localIndex> (const
// localIndex& field, int const )
//  {
//    return field;
//  }


template<> inline globalIndex CastField<globalIndex, globalIndex> (const globalIndex& field, int const )
{
  return field;
}

template<> inline int CastField<int, long long unsigned int > (const long long unsigned int& field, int const )
{
  return static_cast<int>(field);
}


template<> inline real64 CastField<real64, real64> (const real64& field, int const )
{
  return field;
}
template<> inline float CastField<float, real64> (const real64& field, int const )
{
  return static_cast<float> (field);
}


template< typename TYPE, typename PTR_TYPE> inline PTR_TYPE* DataPtr( TYPE& data)
{
  return static_cast<PTR_TYPE*>(&data);
}

template<> inline real64* DataPtr( R1Tensor& data)        {    return data.begin();  }
template<> inline real64* DataPtr( R2Tensor& data)        {    return data.begin();  }
template<> inline real64* DataPtr( R2SymTensor& data)  {    return data.begin();  }



template<typename TYPE>
void SetVariableNames(string const & fieldName, array<string>& varnamestring, char* varnames[]);

template<typename OBJECT_TYPE>
int FieldCentering();

void SetCenteringSubdir(int const centering, string& subdir);

}


template< typename OUTPUTTYPE >
void SiloFile::WriteViewWrappersToSilo( string const & meshname,
                                        const dataRepository::ManagedGroup::viewWrapperMap & wrappers,
                                        int const centering,
                                        int const cycleNum,
                                        real64 const problemTime,
                                        bool const isRestart,
                                        string const & multiRoot,
                                        string const & materialName,
                                        array<localIndex> const & zoneToMatMap,
                                        const localIndex_array& mask )
{

  // iterate over all entries in the member map
  for( auto const & wrapperIter : wrappers )
  {
    auto const & wrapper = wrapperIter.second;
    // the field name is the key to the map
    string const fieldName = wrapper->getName();

    std::type_info const & typeID = wrapper->get_typeid();

    // TODO This is wrong. problem with uniqueness
    if( typeID==typeid(real64_array) )
    {
      auto const & viewWrapperT = dynamic_cast< dataRepository::ViewWrapper<real64_array> const & >( *wrapper );
      this->WriteDataField<real64>(meshname.c_str(), fieldName,
                                   viewWrapperT.reference(), centering, cycleNum, problemTime, multiRoot, materialName, zoneToMatMap );
    }
    if( typeID==typeid(r1_array) )
    {
      auto const & viewWrapperT = dynamic_cast< dataRepository::ViewWrapper<r1_array> const & >( *wrapper );
      this->WriteDataField<real64>(meshname.c_str(), fieldName,
                                   viewWrapperT.reference(), centering, cycleNum, problemTime, multiRoot, materialName, zoneToMatMap );
    }
    if( typeID==typeid(integer_array) )
    {
      auto const & viewWrapperT = dynamic_cast< dataRepository::ViewWrapper<integer_array> const & >( *wrapper );
      this->WriteDataField<integer>(meshname.c_str(), fieldName,
                                    viewWrapperT.reference(), centering, cycleNum, problemTime, multiRoot, materialName, zoneToMatMap );
    }
    if( typeID==typeid(localIndex_array) )
    {
      auto const & viewWrapperT = dynamic_cast< dataRepository::ViewWrapper<localIndex_array> const & >( *wrapper );
      this->WriteDataField<localIndex>(meshname.c_str(), fieldName,
                                       viewWrapperT.reference(), centering, cycleNum, problemTime, multiRoot, materialName, zoneToMatMap );
    }

  }
}


/**
 *
 * @param meshName
 * @param fieldName
 * @param field
 * @param centering
 * @param cycleNumber
 * @param problemTime
 */

template<typename OUTTYPE, typename TYPE>
void SiloFile::WriteDataField( string const & meshName,
                               string const & fieldName,
                               const TYPE& field,
                               int const centering,
                               int const cycleNumber,
                               real64 const problemTime,
                               string const & multiRoot,
                               string const & materialName )
{}

template<typename OUTTYPE, typename TYPE>
void SiloFile::WriteDataField( string const & meshName,
                               string const & fieldName,
                               const array<TYPE>& field,
                               int const centering,
                               int const cycleNumber,
                               real64 const problemTime,
                               string const & multiRoot,
                               string const & materialName,
                               array<localIndex> const & zoneToMatMap )
{
  int const nvars = SiloFileUtilities::GetNumberOfVariablesInField<TYPE>();
  int nels = field.size();

  int const meshType = GetMeshType( meshName );


  DBoptlist *optlist = DBMakeOptlist(5);
  DBAddOption(optlist, DBOPT_CYCLE, const_cast<int*> (&cycleNumber));
  DBAddOption(optlist, DBOPT_DTIME, const_cast<real64*> (&problemTime));

  char *regionpnames[2] =
  { nullptr, nullptr };


  array<char*> varnames(nvars);
  array<void*> vars(nvars);


  array<string> varnamestring(nvars);
  std::vector<std::vector<OUTTYPE> > castedField(nvars);

  if( materialName != "none" )
  {
    regionpnames[0] = const_cast<char*> (materialName.c_str());
    regionpnames[1] = nullptr;
    DBAddOption(optlist, DBOPT_REGION_PNAMES, &regionpnames);

    nels=0;
    for( localIndex i=0 ; i<zoneToMatMap.size() ; ++i )
    {
      if( zoneToMatMap[i] > -1 )
      {
        ++nels;
      }
    }

    for( localIndex i = 0 ; i < nvars ; ++i )
    {
      castedField[i].resize(nels);
      vars[i] = static_cast<void*> (&(castedField[i][0]));
      for( int k = 0 ; k < nels ; ++k )
      {
        castedField[i][k] = SiloFileUtilities::CastField<OUTTYPE>(field[zoneToMatMap[k]], i);
      }
    }

  }
  else
  {
    for( int i = 0 ; i < nvars ; ++i )
    {
      castedField[i].resize(nels);
      vars[i] = static_cast<void*> (&(castedField[i][0]));
      for( int k = 0 ; k < nels ; ++k )
      {
        castedField[i][k] = SiloFileUtilities::CastField<OUTTYPE>(field[k], i);
      }
    }

  }

  // if the number of elements is zero, then record the path to the var. This
  // will be used later to delete the entry
  // from the multivar.
  if( nels == 0 )
  {
    char pwd[256];
    DBGetDir(m_dbFilePtr, pwd);
    string emptyObject = pwd;
    emptyObject += "/" + fieldName;
    m_emptyVariables.push_back(emptyObject);
  }
  else
  {

    SiloFileUtilities::SetVariableNames<TYPE>(fieldName, varnamestring, varnames.data() );


    int err = -2;
    if( meshType == DB_UCDMESH )
    {
      err = DBPutUcdvar( m_dbFilePtr, fieldName.c_str(), meshName.c_str(), nvars, varnames.data(), reinterpret_cast<float**>(vars.data()),
                         nels, nullptr, 0, SiloFileUtilities::DB_TYPE<OUTTYPE>(), centering, optlist);
    }
    else if( meshType == DB_POINTMESH )
    {
      err = DBPutPointvar( m_dbFilePtr, fieldName.c_str(), meshName.c_str(), nvars, reinterpret_cast<float**>(vars.data()),
                           nels, SiloFileUtilities::DB_TYPE<OUTTYPE>(), optlist);
    }
    else if( meshType == DB_QUADCURV )
    {
      err = DBPutQuadvar( m_dbFilePtr, fieldName.c_str(), meshName.c_str(), nvars, varnames.data(), reinterpret_cast<float**>(vars.data()),
                          m_quadMeshDims,m_quadMeshNDims,nullptr, 0,  SiloFileUtilities::DB_TYPE<OUTTYPE>(),centering, optlist);
    }
    if( err < 0 )
    {
      if( err < -1 )
      {
        GEOS_ERROR("unhandled case in SiloFile::WriteDataField A\n");
      }
      else
      {
        GEOS_ERROR("unhandled failure in adding variable during SiloFile::WriteDataField\n");
      }
    }
  }

  // write multimesh object
  int rank = 0;
#if USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  if( rank == 0 )
  {
    int tensorRank = SiloFileUtilities::GetTensorRank<TYPE>();
    DBAddOption(optlist, DBOPT_TENSOR_RANK, const_cast<int*> (&tensorRank));
    DBAddOption(optlist, DBOPT_MMESH_NAME, const_cast<char*> (meshName.c_str()));

    DBObjectType vartype = DB_INVALID_OBJECT;

    if( meshType == DB_UCDMESH )
    {
      vartype = DB_UCDVAR;
    }
    else if( meshType == DB_POINTMESH )
    {
      vartype = DB_POINTVAR;
    }
    else if( meshType == DB_QUADCURV )
    {
      vartype = DB_QUADVAR;
    }
    else
    {
      GEOS_ERROR("unhandled case in SiloFile::WriteDataField B\n");
    }


    WriteMultiXXXX(vartype, DBPutMultivar, centering, fieldName.c_str(), cycleNumber, multiRoot,
                   optlist);
  }

  DBFreeOptlist(optlist);

}

template< typename INPUTTYPE, typename TYPE >
void SiloFile::ReadFieldMapFromSilo( std::map< string, array<TYPE> >& member,
                                     string const & meshname,
                                     int const centering,
                                     int const cycleNum,
                                     real64 const problemTime,
                                     bool const isRestart,
                                     string const & regionName,
                                     const localIndex_array& mask ) const
{
  // iterate over all entries in the member map
  for( typename std::map< string, array<TYPE> >::iterator iter = member.begin() ; iter!=member.end() ; ++iter )
  {
    // the field name is the key to the map
    string const fieldName = iter->first;

    // check to see if the field should have been written
    if( FieldInfo::AttributesByName.find(fieldName) != FieldInfo::AttributesByName.end() )
    {
      if( (  isRestart && FieldInfo::AttributesByName[fieldName]->m_WriteToRestart) ||
          ( !isRestart && FieldInfo::AttributesByName[fieldName]->m_WriteToPlot ) )
      {
        // the field data is mapped value
        array<TYPE> & fieldData = iter->second;

        if( !(mask.empty()) && !isRestart )
        {
          array<TYPE> dataToRead( mask.size() );
          // write the data field
          ReadDataField<INPUTTYPE>( dataToRead, meshname.c_str(), fieldName, centering, cycleNum, problemTime, regionName );

          for( localIndex_array::size_type i = 0 ; i < mask.size() ; ++i )
          {
            fieldData[mask[i]] = dataToRead[i];
          }

        }
        else
        {
          ReadDataField<INPUTTYPE>( fieldData, meshname.c_str(), fieldName, centering, cycleNum, problemTime, regionName );
        }
      }
    }
  }
}

template< typename INPUTTYPE, typename TYPE>
void SiloFile::ReadDataField( array<TYPE>& field,
                              string const & meshName,
                              string const & fieldName,
                              int const centering,
                              int const cycleNumber,
                              real64 const problemTime,
                              string const & regionName ) const
{

  INPUTTYPE** var = static_cast<INPUTTYPE**>( GetDataVar<TYPE>( fieldName, meshName, field.size(), centering, cycleNumber, problemTime, regionName ) );

  for( typename array<TYPE>::size_type a=0 ; a<field.size() ; ++a )
  {
    TYPE temp;
    INPUTTYPE* ptemp = SiloFileUtilities::DataPtr<TYPE,INPUTTYPE>( temp );

    for( int i=0 ; i<SiloFileUtilities::GetNumberOfVariablesInField<TYPE>() ; ++i )
    {
      ptemp[i] = var[i][a];
    }
    field[a] = temp;
  }

}

/**
 *
 * @param type
 * @param DBPutMultiCB
 * @param centering
 * @param name
 * @param cycleNumber
 * @param optlist
 */
template< typename CBF >
void SiloFile::WriteMultiXXXX( const DBObjectType type,
                               CBF DBPutMultiCB,
                               int const centering,
                               string const name,
                               const int,
                               string const & multiRoot,
                               const DBoptlist* optlist)
{
  (void)centering;

  int size = 1;
#if USE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

  array<string> vBlockNames(size);
  std::vector<char*> BlockNames(size);
  std::vector<int> blockTypes(size);
  char tempBuffer[1024];
  char currentDirectory[256];

  DBGetDir(m_dbBaseFilePtr, currentDirectory);
  DBSetDir(m_dbBaseFilePtr, multiRoot.c_str());


  string multiRootString(multiRoot);
  if( !(multiRootString.compare("/")) )
  {
    multiRootString.clear();
  }

  for( int i = 0 ; i < size ; ++i )
  {
    int groupRank = PMPIO_GroupRank(m_baton, i);


    sprintf( tempBuffer,
             "%s%s%s.%03d:/domain_%05d%s/%s",
             m_siloDataSubDirectory.c_str(),
             "/",
             m_baseFileName.c_str(),
             groupRank,
             i,
             multiRootString.c_str(),
             name.c_str());

    vBlockNames[i] = tempBuffer;
    BlockNames[i] = const_cast<char*>( vBlockNames[i].c_str() );
    blockTypes[i] = type;
  }

  string multiName = name;
  DBPutMultiCB(m_dbBaseFilePtr, multiName.c_str(), size, BlockNames.data(), blockTypes.data(),
               const_cast<DBoptlist*> (optlist));

  DBSetDir(m_dbBaseFilePtr, currentDirectory);

}



}
#endif /* SILOFILE_H_ */
