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

#ifdef GEOSX_USE_MPI
#include <mpi.h>
#endif

#include "mpi.h"
#include "pmpio.h"

#include "mesh/ElementRegionManager.hpp"
#include "mesh/InterObjectRelation.hpp"


typedef struct _PMPIO_baton_t PMPIO_baton_t;


namespace geosx
{

class DomainPartition;
class MeshLevel;
namespace constitutive
{
class ConstitutiveManager;
}

// *********************************************************************************************************************
// *********************************************************************************************************************
/**
 * This class serves as a wrapper to isolate the code from the specifics of SILO
 * output/input. Its members contain all the necessary information for reading/writing a group of SILO files.
 */
class SiloFile
{

public:

  /// Default Constructor
  SiloFile();

  /// Destructor
  virtual ~SiloFile();

  /**
   * @brief function to setup directories where silo files will be written
   */
  void MakeSiloDirectories();

  /**
   * @brief Initializes silo for input/output
   * @param readwrite input/output specifier
   */
  void Initialize( const PMPIO_iomode_t readwrite, int const numGroups=1 );

  /**
   * @brief finishes/closes up the silo interface
   */
  void Finish();

  /**
   * @brief Wait for the Baton when writing using PMPIO
   * @param domainNumber domain partition number
   * @param cycleNum  cycle number of simulation
   * @param isRestart whether or not we are writing a restart file
   *
   * This function requests the write baton from silo PMPIO. The involves determining
   * the file names, and opening the file for write.
   */
  void WaitForBatonWrite( int const domainNumber,
                          int const cycleNum,
                          integer const eventCounter,
                          bool const isRestart );

  /**
   * @brief Wait for the Baton when reading using PMPIO
   * @param domainNumber domain partition number
   * @param restartFileName base of the restart file to open
   */
  void WaitForBaton( int const domainNumber, string const & restartFileName );

  /**
   * @brief Hand off the Baton when done writing to file
   */
  void HandOffBaton();


  /**
   * @brief Make a subdirectory within the silo file.
   * @param subdir the new directory name
   * @param rootdir the root directory path
   */
  void MakeSubDirectory( string const & subdir, string const & rootdir )
  {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_GEOSX, &rank);

    char dirname[100];
    if( rank == 0 )
    {
//      DBGetDir(m_dbBaseFilePtr, dirname );
      DBMkDir(m_dbBaseFilePtr, rootdir.c_str());
    }

//    DBGetDir (m_dbFilePtr, dirname );
    DBMkDir(m_dbFilePtr, subdir.c_str());
  }

  /**
   * @brief Write out a single silo mesh object
   * @param meshName name of the mesh in the silo db
   * @param nnodes number of nodes
   * @param coords array[3] of pointers to x, y, and z.
   * @param globalNodeNum array to the global node numbers. This might be redundant as there is a field for this.
   * @param numShapes number of element zone type (i.e. number of zone types with different topology)
   * @param shapecnt pointer to array that contains the number of zones per shape type
   * @param meshConnectivity pointer to array that contains the zone to element map for each  zone type
   * @param globalElementNum pointer to array of global zone numbers for each shape type
   * @param
   * @param shapetype pointer to array containing the shape types
   * @param shapesize pointer to array containing the number of nodes in each zone in the shape types
   * @param cycleNumber the current cycle number
   * @param problemTime the current problem time
   *
   * This function takes in the the required data to call a silo::DBPutUcdMesh() and a
   * DBPutZonelist2, and calls those functions to create a silo mesh object. In addition
   * the MultiVar is written in the root file.
   */
  void WriteMeshObject(string const & meshName,
                       const localIndex nnodes,
                       real64* coords[3],
                       const globalIndex* globalNodeNum,
                       char const * const ghostNodeName,
                       char const * const ghostZoneName,
                       int const numShapes,
                       int const * shapecnt,
                       const localIndex* const * const meshConnectivity,
                       const globalIndex* const * const globalElementNum,
                       int const * const shapetype,
                       int const * const shapesize,
                       int const cycleNumber,
                       real64 const problemTime);


  void WritePolygonMeshObject(const std::string& meshName,
                                        const localIndex nnodes,
                                        realT* coords[3],
                                        const globalIndex*,
                                        const int numRegions,
                                        const int* shapecnt,
                                        const localIndex* const * const meshConnectivity,
                                        const globalIndex* const * const globalElementNum,
                                        const int* const * const,
                                        const int* const shapetype,
                                        const int* const shapesize,
                                        const int cycleNumber,
                                        const realT problemTime,
                                        const int lnodelist);
/**
 * @brief write a domain parititon out to silo file
 * @param domain the domain partition to write
 * @param cycleNum the current cycle number
 * @param problemTime the current problem time
 * @param isRestart whether or not we want to write restart only data
 */
  void WriteDomainPartition( DomainPartition const & domain,
                             int const cycleNum,
                             real64 const problemTime,
                             bool const isRestart );

  /**
   * @brief write a mesh level out to the silo file
   * @param meshLevel the meshLevel to write out
   * @param constitutiveManager the constitutive manager object that holds the constitutive data
   * @param cycleNum the current cycle number
   * @param problemTime the current problem time
   * @param isRestart whether or not we want to write restart only data
   */
  void WriteMeshLevel( MeshLevel const * const meshLevel,
                       constitutive::ConstitutiveManager const * const constitutiveManager,
                       int const cycleNum,
                       real64 const problemTime,
                       bool const isRestart );

  /**
   * @param meshName name of the mesh to write
   * @param numPoints number of points to write
   * @param coords array[3] of pointers to x, y, and z.
   * @param cycleNumber current cycle number
   * @param problemTime current problem time
   */
  void WritePointMesh( string const & meshName,
                       const localIndex numPoints,
                       real64* coords[3],
                       int const cycleNumber,
                       real64 const problemTime );

  /**
   * @brief write a beam mesh to silo file
   * @param meshName name of the mesh to write
   * @param nnodes number of nodes
   * @param coords array[3] of pointers to x, y, and z.
   * @param node1 array of the first node for each element
   * @param node2 array of the second node for each element
   * @param cycleNumber current cycle number
   * @param problemTime current problem time
   */
  void WriteBeamMesh(string const & meshName,
                     const localIndex nnodes,
                     real64* coords[3],
                     const localIndex_array& node1,
                     const localIndex_array& node2,
                     int const cycleNumber,
                     real64 const problemTime);

  /**
   *
   * @param meshName name of the mesh to write
   * @param nnodes number of nodes
   * @param coords array[3] of pointers to x, y, and z.
   * @param nodelist nodal connectivity array
   * @param cycleNumber current cycle number
   * @param problemTime current problem time
   */
  void WriteBeamMesh(string const & meshName,
                     const localIndex nnodes,
                     real64* coords[3],
                     integer_array& nodelist,
                     int const cycleNumber,
                     real64 const problemTime);


  /**
   *
   * @param elementManager the element region manager
   * @param constitutiveManager the constitutive manager
   * @param meshName the name of the mesh that this write applies to
   * @param cycleNumber current cycle number
   * @param problemTime current problem time
   */
  void WriteMaterialMapsCompactStorage( ElementRegionManager const * const elementManager,
                          constitutive::ConstitutiveManager const * const constitutiveManager,
                          string const & meshName,
                          int const cycleNumber,
                          real64 const problemTime);


  void WriteMaterialMapsFullStorage( ElementRegionManager const * const elementManager,
                          constitutive::ConstitutiveManager const * const constitutiveManager,
                          string const & meshName,
                          int const cycleNumber,
                          real64 const problemTime);
  /**
   *
   * @param group the group that holds the data to be written to the silo file
   * @param siloDirName the name of the silo directory to put this data into (i.e. nodalFields, elemFields)
   * @param meshname the name of the mesh attach this write to
   * @param centering the centering of the data (e.g. DB_ZONECENT)
   * @param cycleNumber current cycle number
   * @param problemTime current problem time
   * @param isRestart write restart only data
   * @param mask indices to write out to the silo file
   */
  void WriteManagedGroupSilo( dataRepository::ManagedGroup const * group,
                              string const & siloDirName,
                              string const & meshname,
                              int const centering,
                              int const cycleNum,
                              real64 const problemTime,
                              bool const isRestart,
                              const localIndex_array& mask );


  /**
   * Writes the contents of a group of ViewWrapper objects
   * @tparam the output varaible type
   * @param meshname the name of the mesh attach this write to
   * @param wrappers a group of wrappers
   * @param centering the silo centering to use for this operation (DB_NODECENT, DB_ZONECENT)
   * @param cycleNum the current cycle number
   * @param problemTime the current problem time
   * @param isRestart write restart only data
   * @param multiRoot location to write the multivar entries
   * @param mask indices to write out to the silo file
   */
  template< typename OUTPUTTYPE >
  void WriteViewWrappersToSilo( string const & meshname,
                                const dataRepository::ManagedGroup::viewWrapperMap & wrappers,
                                int const centering,
                                int const cycleNum,
                                real64 const problemTime,
                                bool const isRestart,
                                string const & multiRoot,
                                const localIndex_array& mask );

  /**
   *
   * @param meshname the name of the mesh attach this write to
   * @param fieldName name of the field to write
   * @param field field data
   * @param centering the silo centering to use for this operation (DB_NODECENT, DB_ZONECENT)
   * @param cycleNum the current cycle number
   * @param problemTime the current problem time
   * @param multiRoot location to write the multivar entries
   */
  template<typename OUTTYPE, typename TYPE>
  void WriteDataField( string const & meshName,
                       string const & fieldName,
                       const array1d<TYPE>& field,
                       int const centering,
                       int const cycleNumber,
                       real64 const problemTime,
                       string const & multiRoot );

  template<typename OUTTYPE, typename TYPE>
  void WriteMaterialDataField( string const & meshName,
                               string const & fieldName,
                               ElementRegionManager::MaterialViewAccessor< arrayView2d<TYPE> > const & field,
                               ElementRegionManager const * const elementManager,
                               constitutive::ConstitutiveManager const * const constitutiveManager,
                               int const centering,
                               int const cycleNumber,
                               real64 const problemTime,
                               string const & multiRoot,
                               string_array const & materialNames );

  /**
   * find the silo mesh type that we are attempting to reference
   * @param meshName the name of the mesh object we are attaching data to
   * @return integer that results from a call to DBInqMeshtype()
   */
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

  /**
   * write out a multi mesh or multivar object assocaited with a mesh or var object.
   * @param type the object type (DB_UCDVAR, ...)
   * @param DBPutMultiCB a callback function for the placement of the multi object. ( DBPutMultivar, DBPutMultimesh...)
   * @param centering the centering of the data (DB_NODECENT, DB_ZONECENT, ...)
   * @param name the name of the multivar
   * @param cycleNumber the current cycle number
   * @param multiRoot the location in the silo file to put hte multiXXXX
   * @param optlist the option list assocaited with the multiXXXX
   */
  template< typename CBF >
  void WriteMultiXXXX(const DBObjectType type, CBF DBPutMultiCB,
                      int const centering, string const name, int const cycleNumber,
                      string const & multiRoot, const DBoptlist* optlist = nullptr);


  /**
   * fucntion to clear any empty multi-objects
   * @param cycleNum
   *
   * When we write our multimesh and multivar objects, we do it assuming that there is a non-empty multimesh
   * or multivar object on each domain. This is incorrect, so we must modify the rootfile to remove empty references
   * to the multivar or multimesh objects and replace their path with "EMPTY"
   */
  void ClearEmptiesFromMultiObjects(int const cycleNum);


  void setNumGroups( int const numGroups )
  {
    m_numGroups = numGroups;
  }

  void setPlotLevel( int const plotLevel )
  {
    m_plotLevel = dataRepository::IntToPlotLevel(plotLevel);
  }
private:

  /// pointer to the DBfile that this class is working on
  DBfile* m_dbFilePtr;

  /// pointer to the Master DBfile
  DBfile* m_dbBaseFilePtr;

  /// total number of "parallel" files to write out
  int m_numGroups;

  /// the pmpio baton. A processor needs this to write to the file.
  PMPIO_baton_t *m_baton;

  /// which database to use. DB_PDB or DB_HDF5
  int const m_driver;

  /// root of the filename that we will be reading/writing
  string m_plotFileRoot;

  string m_restartFileRoot;

  string const m_siloDirectory = "siloFiles";

  string const m_siloDataSubDirectory = "data";

  string m_fileName;

  string m_baseFileName;

  string_array m_emptyMeshes;
  string_array m_emptyVariables;

  dataRepository::PlotLevel m_plotLevel;

  /**
   *
   * @return returns the ordering of nodes for a silo zone type.
   */
  integer_array SiloNodeOrdering(const string & elementType);







};

/**
 * Namespace to hold some utilities needed by the functions in SiloFile.
 */
namespace SiloFileUtilities
{
/**
 * @tparam OUTTYPE the type of data to write out (int,float,real64)
 * @return the integer identifier associated with the enum DBdatatype defined in
 * the silo.h file.
 *
 * This templated function is a "specialization only" definition. There is no
 * general definition, only specializations for predetermined data types.
 */
template<typename OUTTYPE>
int DB_TYPE();

/**
 * @tparam TYPE the data type in question
 * @return the number of "variables" in a data field. For instance, 1 for a int,
 * float or real64.
 * 1 for a scalar
 * 3 for a R1Tensor...etc.
 */
template<typename TYPE>
int GetNumberOfVariablesInField();

/**
 * @tparam TYPE the data type in question
 * @return the silo DB_VARTYPE specifier that describes the rank of the variable.
 */
template<typename TYPE>
int GetTensorRank();

/**
 * @tparam OUTTYPE the type to cast to
 * @tparam TYPE the type to cast from
 * @param field the value to cast
 * @param i the component of the varaible to cast, assuming there is dimensionaliy to the variable
 * @return the casted value
 */
template<typename OUTTYPE, typename TYPE>
OUTTYPE CastField(const TYPE& field, int const i = 0);     // avoids compiler
                                                           // warning

template<typename OUTTYPE, typename TYPE>
OUTTYPE CastField(const TYPE& field, int const i)
{
  return field.Data()[i];
}

template<> inline int CastField<int, int> (const int& field, int const )
{
  return field;
}


template<> inline long int CastField<long int, long int> (const long int& field, int const )
{
  return field;
}

template<> inline int CastField<int, long int> (const long int& field, int const )
{
  return integer_conversion<int>(field);
}

template<> inline long long int CastField<long long int, long long int> (const long long int& field, int const )
{
  return field;
}

template<> inline int CastField<int, long long int> (const long long int& field, int const )
{
  return integer_conversion<int>(field);
}

template<> inline real64 CastField<real64, real64> (const real64& field, int const )
{
  return field;
}
template<> inline float CastField<float, real64> (const real64& field, int const )
{
  return static_cast<float>(field);
}



/**
 * @tparam the type of the field
 * @param fieldName the name of the field
 * @param varnamestring an array of strings that hold the generated names
 * @param varnames a pointer to each of the strings that hold the generated names
 *
 * This function sets variable names for a field. In some cases, this will just fill varnamestring with
 * the input variable, and in other cases such as in the case of a tensor, the names will have a suffix
 * with the component of the tensor.
 */
template<typename TYPE>
void SetVariableNames(string const & fieldName, string_array& varnamestring, char const* varnames[]);


}


template< typename OUTPUTTYPE >
void SiloFile::WriteViewWrappersToSilo( string const & meshname,
                                        const dataRepository::ManagedGroup::viewWrapperMap & wrappers,
                                        int const centering,
                                        int const cycleNum,
                                        real64 const problemTime,
                                        bool const isRestart,
                                        string const & multiRoot,
                                        const localIndex_array& mask )
{

  // iterate over all entries in the member map
  for( auto const & wrapperIter : wrappers )
  {
    auto const & wrapper = wrapperIter.second;

    if( wrapper->getPlotLevel() < m_plotLevel )
    {
      // the field name is the key to the map
      string const fieldName = wrapper->getName();

      std::type_info const & typeID = wrapper->get_typeid();

      // TODO This is wrong. problem with uniqueness
      if( typeID==typeid(real64_array) )
      {
        auto const & viewWrapperT = dynamic_cast< dataRepository::ViewWrapper<real64_array> const & >( *wrapper );
        this->WriteDataField<real64>(meshname.c_str(), fieldName,
                                     viewWrapperT.reference(), centering, cycleNum, problemTime, multiRoot );
      }
      if( typeID==typeid(r1_array) )
      {
        auto const & viewWrapperT = dynamic_cast< dataRepository::ViewWrapper<r1_array> const & >( *wrapper );
        this->WriteDataField<real64>(meshname.c_str(), fieldName,
                                     viewWrapperT.reference(), centering, cycleNum, problemTime, multiRoot );
      }
      if( typeID==typeid(integer_array) )
      {
        auto const & viewWrapperT = dynamic_cast< dataRepository::ViewWrapper<integer_array> const & >( *wrapper );
        this->WriteDataField<integer>(meshname.c_str(), fieldName,
                                      viewWrapperT.reference(), centering, cycleNum, problemTime, multiRoot );
      }
      if( typeID==typeid(localIndex_array) )
      {
        auto const & viewWrapperT = dynamic_cast< dataRepository::ViewWrapper<localIndex_array> const & >( *wrapper );
        this->WriteDataField<localIndex>(meshname.c_str(), fieldName,
                                         viewWrapperT.reference(), centering, cycleNum, problemTime, multiRoot );
      }
      if( typeID==typeid(globalIndex_array) )
      {
        auto const & viewWrapperT = dynamic_cast< dataRepository::ViewWrapper<globalIndex_array> const & >( *wrapper );
        this->WriteDataField<globalIndex>(meshname.c_str(), fieldName,
                                         viewWrapperT.reference(), centering, cycleNum, problemTime, multiRoot );
      }
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
                               const array1d<TYPE>& field,
                               int const centering,
                               int const cycleNumber,
                               real64 const problemTime,
                               string const & multiRoot )
{
  int const nvars = SiloFileUtilities::GetNumberOfVariablesInField<TYPE>();
  int nels = field.size();

  int const meshType = GetMeshType( meshName );


  DBoptlist *optlist = DBMakeOptlist(5);
  DBAddOption(optlist, DBOPT_CYCLE, const_cast<int*> (&cycleNumber));
  DBAddOption(optlist, DBOPT_DTIME, const_cast<real64*> (&problemTime));

  array1d<char const*> varnames(nvars);
  array1d<void*> vars(nvars);


  string_array varnamestring(nvars);
  std::vector<std::vector<OUTTYPE> > castedField(nvars);


  for( int i = 0 ; i < nvars ; ++i )
  {
    if( std::is_same<OUTTYPE,TYPE>::value )
    {
      vars[i] = const_cast<void*>(static_cast<void const*> (&(field[0])+i));
    }
    else
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
      err = DBPutUcdvar( m_dbFilePtr,
                         fieldName.c_str(),
                         meshName.c_str(),
                         nvars,
                         varnames.data(),
                         reinterpret_cast<float**>(vars.data()),
                         nels,
                         nullptr,
                         0,
                         SiloFileUtilities::DB_TYPE<OUTTYPE>(),
                         centering,
                         optlist);
    }
    else if( meshType == DB_POINTMESH )
    {
      err = DBPutPointvar( m_dbFilePtr,
                           fieldName.c_str(),
                           meshName.c_str(),
                           nvars,
                           reinterpret_cast<float**>(vars.data()),
                           nels,
                           SiloFileUtilities::DB_TYPE<OUTTYPE>(),
                           optlist);
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
#ifdef GEOSX_USE_MPI
  MPI_Comm_rank(MPI_COMM_GEOSX, &rank);
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

template<typename OUTTYPE, typename TYPE>
void SiloFile::WriteMaterialDataField( string const & meshName,
                                       string const & fieldName,
                                       ElementRegionManager::MaterialViewAccessor< arrayView2d<TYPE> > const & field,
                                       ElementRegionManager const * const elementManager,
                                       constitutive::ConstitutiveManager const * const constitutiveManager,
                                       int const centering,
                                       int const cycleNumber,
                                       real64 const problemTime,
                                       string const & multiRoot,
                                       string_array const & materialNames )
{
  int const nvars = SiloFileUtilities::GetNumberOfVariablesInField<TYPE>();

  int const meshType = GetMeshType( meshName );

  string_array activeMaterialNames;

//  double missingValue = 0.0;
  DBoptlist *optlist = DBMakeOptlist(5);
  DBAddOption(optlist, DBOPT_CYCLE, const_cast<int*> (&cycleNumber));
  DBAddOption(optlist, DBOPT_DTIME, const_cast<real64*> (&problemTime));
//  DBAddOption(optlist, DBOPT_MISSING_VALUE, &missingValue);

  char * regionpnames[ 100 ];

  // if the number of elements is zero, then record the path to the var. This
  // will be used later to delete the entry
  // from the multivar.
  if( field.size() == 0 )
  {
    char pwd[256];
    DBGetDir(m_dbFilePtr, pwd);
    string emptyObject = pwd;
    emptyObject += "/" + fieldName;
    m_emptyVariables.push_back(emptyObject);
  }
  else
  {

    string_array varnamestring(nvars);
    array1d<char const*> varnames(nvars);

    SiloFileUtilities::SetVariableNames<TYPE>(fieldName, varnamestring, varnames.data() );

    int nels = 0;
    localIndex mixlen = 0;
    for( localIndex er=0 ; er<elementManager->numRegions() ; ++er )
    {
      ElementRegion const * const elemRegion = elementManager->GetRegion(er);
      localIndex const numMatInRegion = integer_conversion<localIndex>(elemRegion->getMaterialList().size());

      array1d<localIndex> matIndices(numMatInRegion);
      for( localIndex a=0 ; a<numMatInRegion ; ++a )
      {
        matIndices[a] = constitutiveManager->
                        GetConstitituveRelation( elemRegion->getMaterialList()[a] )->
                        getIndexInParent();
      }

      for( localIndex esr=0 ; esr<elemRegion->numSubRegions() ; ++esr )
      {
        CellBlockSubRegion const * const subRegion = elemRegion->GetSubRegion(esr);

        nels += subRegion->size();

        for( localIndex matIndex=0 ; matIndex<numMatInRegion ; ++matIndex )
        {
          if( field[er][esr][matIndices[matIndex]].size() > 0 )
          {
            activeMaterialNames.push_back( constitutiveManager->GetConstitituveRelation( matIndices[matIndex] )->getName() );
            mixlen += subRegion->size();
          }
        }
      }
    }

    std::vector<void*> vars(nvars);
    std::vector<std::vector<OUTTYPE> > varsData(nvars);

    std::vector<void*> mixvars(nvars);
    std::vector<std::vector<OUTTYPE> > mixvarsData(nvars);

    for( int a=0 ; a<nvars ; ++a )
    {
      varsData[a].resize(nels);
      mixvarsData[a].resize(mixlen);

      vars[a] = static_cast<void*> (&(varsData[a][0]));
      mixvars[a] = static_cast<void*> (&(mixvarsData[a][0]));

    }

    localIndex mixlen2 = 0;
    localIndex nels2 = 0;

    for( localIndex er=0 ; er<elementManager->numRegions() ; ++er )
    {
      ElementRegion const * const elemRegion = elementManager->GetRegion(er);
      localIndex const numMatInRegion = integer_conversion<localIndex>(elemRegion->getMaterialList().size());

      array1d<localIndex> matIndices(numMatInRegion);

      for( localIndex a=0 ; a<numMatInRegion ; ++a )
      {
        matIndices[a] = constitutiveManager->
                        GetConstitituveRelation( elemRegion->getMaterialList()[a] )->
                        getIndexInParent();
      }

      for( localIndex esr=0 ; esr<elemRegion->numSubRegions() ; ++esr )
      {
        CellBlockSubRegion const * const subRegion = elemRegion->GetSubRegion(esr);

        if( numMatInRegion == 1 )
        {
          for( int i = 0 ; i < nvars ; ++i )
          {
            for( localIndex k = 0 ; k < subRegion->size() ; ++k )
            {
              varsData[i][nels2++] = SiloFileUtilities::CastField<OUTTYPE>(field[er][esr][matIndices[0]][k][0], i);
              mixvarsData[i][mixlen2++] = SiloFileUtilities::CastField<OUTTYPE>(field[er][esr][matIndices[0]][k][0], i);
            }
          }
        }
        else if( numMatInRegion > 1 )
        {
          for( int i = 0 ; i < nvars ; ++i )
          {
            for( localIndex k = 0 ; k < subRegion->size() ; ++k )
            {
              varsData[i][nels2++] = SiloFileUtilities::CastField<OUTTYPE>(field[er][esr][matIndices[0]][k][0], i);
              for( localIndex a=0 ; a<numMatInRegion ; ++a )
              {
                if( field[er][esr][matIndices[a]].size() > 0 )
                {
                  mixvarsData[i][mixlen2++] = SiloFileUtilities::CastField<OUTTYPE>(field[er][esr][matIndices[a]][k][0], i);
                }
                else
                {
                  mixvarsData[i][mixlen2++] = 0.0;
                }
              }
            }
          }
        }
      }
    }
    
    for( string_array::size_type a=0 ; a<activeMaterialNames.size() ; ++a )
    {
      regionpnames[a] = const_cast<char*> (activeMaterialNames[a].c_str());
    }
    regionpnames[activeMaterialNames.size()] = nullptr;
    DBAddOption(optlist, DBOPT_REGION_PNAMES, &regionpnames );


    int err = -2;
    if( meshType == DB_UCDMESH )
    {
      err = DBPutUcdvar( m_dbFilePtr,
                         fieldName.c_str(),
                         meshName.c_str(),
                         nvars,
                         varnames.data(),
                         reinterpret_cast<void**>(vars.data()),
                         nels,
                         reinterpret_cast<void**>(mixvars.data()),
                         mixlen,
                         SiloFileUtilities::DB_TYPE<OUTTYPE>(),
                         centering,
                         optlist);
    }
    else if( meshType == DB_POINTMESH )
    {
      err = DBPutPointvar( m_dbFilePtr,
                           fieldName.c_str(),
                           meshName.c_str(),
                           nvars,
                           reinterpret_cast<float**>(vars.data()),
                           nels,
                           SiloFileUtilities::DB_TYPE<OUTTYPE>(),
                           optlist);
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
#ifdef GEOSX_USE_MPI
  MPI_Comm_rank(MPI_COMM_GEOSX, &rank);
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
#ifdef GEOSX_USE_MPI
  MPI_Comm_size(MPI_COMM_GEOSX, &size);
#endif

  string_array vBlockNames(size);
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
