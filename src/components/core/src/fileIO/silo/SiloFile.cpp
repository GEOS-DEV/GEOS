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
 * @file SiloFile.cpp
 */

#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>

#include "managers/DomainPartition.hpp"
#include "mesh/MeshBody.hpp"


#pragma GCC diagnostic push

#ifdef __clang__
#pragma GCC diagnostic ignored "-Wshorten-64-to-32"
#endif

#pragma GCC diagnostic ignored "-Wold-style-cast"

#include "SiloFile.hpp"
//#include "ObjectManagers/ElementManagerT.h"

#include "common/Logger.hpp"

#include "codingUtilities/Utilities.hpp"

#include "constitutive/ConstitutiveManager.hpp"

#include "mesh/MeshBody.hpp"

/// forward declaration of NodeManagerT for use as a template argument
class NodeManager;

namespace geosx
{


/**
 *
 * @return
 */
namespace SiloFileUtilities
{

template<> int DB_TYPE<int> ()
{
  return DB_INT;
}
template<> int DB_TYPE<unsigned int> ()
{
  return DB_INT;
}
template<> int DB_TYPE<float> ()
{
  return DB_FLOAT;
}
template<> int DB_TYPE<realT> ()
{
  return DB_DOUBLE;
}
template<> int DB_TYPE<R1Tensor> ()
{
  return DB_DOUBLE;
}
template<> int DB_TYPE<R2Tensor> ()
{
  return DB_DOUBLE;
}
template<> int DB_TYPE<R2SymTensor> ()
{
  return DB_DOUBLE;
}
template<> int DB_TYPE<unsigned long> ()
{
  return DB_LONG;
}
template<> int DB_TYPE<long> ()
{
  return DB_LONG;
}
template<> int DB_TYPE<long long> ()
{
  return DB_LONG_LONG;
}
template<> int DB_TYPE<string> ()
{
  return DB_CHAR;
}

template<> int GetNumberOfVariablesInField<int> ()
{
  return 1;
}
template<> int GetNumberOfVariablesInField<unsigned int> ()
{
  return 1;
}
template<> int GetNumberOfVariablesInField<unsigned long> ()
{
  return 1;
}
template<> int GetNumberOfVariablesInField<long> ()
{
  return 1;
}
template<> int GetNumberOfVariablesInField<float> ()
{
  return 1;
}
template<> int GetNumberOfVariablesInField<realT> ()
{
  return 1;
}
template<> int GetNumberOfVariablesInField<long long unsigned int> ()
{
  return 1;
}
template<> int GetNumberOfVariablesInField<long long int> ()
{
  return 1;
}
template<> int GetNumberOfVariablesInField<R1Tensor> ()
{
  return R1Tensor::Length();
}
template<> int GetNumberOfVariablesInField<R2Tensor> ()
{
  return R2Tensor::Length();
}
template<> int GetNumberOfVariablesInField<R2SymTensor> ()
{
  return R2SymTensor::Length();
}
template<> int GetNumberOfVariablesInField<string> ()
{
  return 1;
}

template<typename TYPE>
void SetVariableNames(const std::string& fieldName, array<string>& varnamestring, char* varnames[])
{
  varnamestring.resize(GetNumberOfVariablesInField<TYPE> ());
  int count = 0;
  for (array<string>::iterator i = varnamestring.begin() ; i != varnamestring.end() ; ++i)
  {
    *i = fieldName;
    varnames[count++] =  const_cast<char*>((*i).c_str());
  }
}
template void SetVariableNames<int> (const std::string& fieldName, array<string>& varnamestring,
                                     char* varnames[]);
template void SetVariableNames<unsigned long> (const std::string& fieldName, array<string>& varnamestring,
                                               char* varnames[]);
template void SetVariableNames<realT> (const std::string& fieldName, array<string>& varnamestring,
                                       char* varnames[]);
template void SetVariableNames<long long unsigned int> (const std::string& fieldName, array<string>& varnamestring,
                                                        char* varnames[]);



template<>
void SetVariableNames<R1Tensor> (const std::string& fieldName, array<string>& varnamestring,
                                 char* varnames[])
{
  varnamestring.resize(GetNumberOfVariablesInField<R1Tensor> ());
  varnamestring[0] = fieldName + "_1";
  varnamestring[1] = fieldName + "_2";
  varnamestring[2] = fieldName + "_3";
  varnames[0] = const_cast<char*>( varnamestring[0].c_str() );
  varnames[1] = const_cast<char*>( varnamestring[1].c_str() );
  varnames[2] = const_cast<char*>( varnamestring[2].c_str() );
}

template<>
void SetVariableNames<R2Tensor> (const std::string& fieldName, array<string>& varnamestring,
                                 char* varnames[])
{
  varnamestring.resize(GetNumberOfVariablesInField<R2Tensor> ());
  varnamestring[0] = fieldName + "_11";
  varnamestring[1] = fieldName + "_12";
  varnamestring[2] = fieldName + "_13";
  varnamestring[3] = fieldName + "_21";
  varnamestring[4] = fieldName + "_22";
  varnamestring[5] = fieldName + "_23";
  varnamestring[6] = fieldName + "_31";
  varnamestring[7] = fieldName + "_32";
  varnamestring[8] = fieldName + "_33";
  varnames[0] = const_cast<char*>( varnamestring[0].c_str() );
  varnames[1] = const_cast<char*>( varnamestring[1].c_str() );
  varnames[2] = const_cast<char*>( varnamestring[2].c_str() );
  varnames[3] = const_cast<char*>( varnamestring[3].c_str() );
  varnames[4] = const_cast<char*>( varnamestring[4].c_str() );
  varnames[5] = const_cast<char*>( varnamestring[5].c_str() );
  varnames[6] = const_cast<char*>( varnamestring[6].c_str() );
  varnames[7] = const_cast<char*>( varnamestring[7].c_str() );
  varnames[8] = const_cast<char*>( varnamestring[8].c_str() );
}

template<>
void SetVariableNames<R2SymTensor> (const std::string& fieldName, array<string>& varnamestring,
                                    char* varnames[])
{
  varnamestring.resize(GetNumberOfVariablesInField<R2Tensor> ());
  varnamestring[0] = fieldName + "_11";
  varnamestring[1] = fieldName + "_21";
  varnamestring[2] = fieldName + "_22";
  varnamestring[3] = fieldName + "_31";
  varnamestring[4] = fieldName + "_32";
  varnamestring[5] = fieldName + "_33";
  varnames[0] = const_cast<char*>( varnamestring[0].c_str() );
  varnames[1] = const_cast<char*>( varnamestring[1].c_str() );
  varnames[2] = const_cast<char*>( varnamestring[2].c_str() );
  varnames[3] = const_cast<char*>( varnamestring[3].c_str() );
  varnames[4] = const_cast<char*>( varnamestring[4].c_str() );
  varnames[5] = const_cast<char*>( varnamestring[5].c_str() );
}

template<> int FieldCentering<NodeManager> ()
{
  return DB_NODECENT;
}

template<> int GetTensorRank<int> ()
{
  return DB_VARTYPE_SCALAR;
}
template<> int GetTensorRank<unsigned int> ()
{
  return DB_VARTYPE_SCALAR;
}
template<> int GetTensorRank<long> ()
{
  return DB_VARTYPE_SCALAR;
}
template<> int GetTensorRank<unsigned long> ()
{
  return DB_VARTYPE_SCALAR;
}
template<> int GetTensorRank<real32> ()
{
  return DB_VARTYPE_SCALAR;
}
template<> int GetTensorRank<realT> ()
{
  return DB_VARTYPE_SCALAR;
}
template<> int GetTensorRank<R1Tensor> ()
{
  return DB_VARTYPE_VECTOR;
}
template<> int GetTensorRank<R2Tensor> ()
{
  return DB_VARTYPE_TENSOR;
}
template<> int GetTensorRank<R2SymTensor> ()
{
  return DB_VARTYPE_SYMTENSOR;
}
template<> int GetTensorRank<long long unsigned int> ()
{
  return DB_VARTYPE_SCALAR;
}
template<> int GetTensorRank<string> ()
{
  return DB_VARTYPE_SCALAR;
}



void SetCenteringSubdir(const int centering, std::string& subdir)
{

  if (centering == DB_NODECENT)
    subdir = "node_fields";
  else if (centering == DB_ZONECENT)
    subdir = "zone_fields";
  else if (centering == DB_FACECENT)
    subdir = "face_fields";
  else if (centering == DB_EDGECENT)
    subdir = "edge_fields";
}
}



using namespace constitutive;
using namespace dataRepository;

// *********************************************************************************************************************
/// Default Constructor
SiloFile::SiloFile():
  m_dbFilePtr(nullptr),
  m_numGroups(1),
  m_baton(nullptr),
  m_driver(DB_HDF5),
  m_fileRoot("silo"),
  m_restartFileRoot("restart"),
  m_fileName(),
  m_baseFileName(),
  m_markGhosts(0)
{}

// *********************************************************************************************************************
/// Destructor
SiloFile::~SiloFile()
{}

// *********************************************************************************************************************
/**
 * @author settgast
 *
 */
void SiloFile::Initialize( const PMPIO_iomode_t readwrite )
{
#if USE_MPI
  // Ensure all procs agree on numGroups, driver and file_ext
  m_numGroups = 2;

  MPI_Bcast(&m_numGroups, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&m_driver, 1, MPI_INT, 0, MPI_COMM_WORLD);
  // Initialize PMPIO, pass a pointer to the driver type as the user data.
  m_baton = PMPIO_Init(m_numGroups, readwrite, MPI_COMM_WORLD, 1, PMPIO_DefaultCreate, PMPIO_DefaultOpen, PMPIO_DefaultClose, &m_driver);
#else
  m_numGroups = 1;
  // Initialize PMPIO, pass a pointer to the driver type as the user data.
  m_baton = PMPIO_Init(m_numGroups, PMPIO_WRITE, nullptr, 1, PMPIO_DefaultCreate, PMPIO_DefaultOpen,
                       PMPIO_DefaultClose, &m_driver);
#endif

}

// *********************************************************************************************************************
/**
 * @author settgast
 *
 */
void SiloFile::Finish()
{
  m_baseFileName.clear();
  PMPIO_Finish(m_baton);
}

// *********************************************************************************************************************
/**
 *
 * @param domainNumber
 * @param cycleNum
 * @param fileName
 * @param nsName
 */
void SiloFile::WaitForBaton( const int domainNumber, const int cycleNum, const bool isRestart )
{

  int rank = 0;
#if USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  const int groupRank = PMPIO_GroupRank(m_baton, rank);
  char fileName[200] = { 0 };
  char baseFileName[200] = { 0 };
  char dirName[200] = { 0 };

  if( isRestart )
  {

    sprintf(baseFileName, "%s_%06d", m_restartFileRoot.c_str(), cycleNum);
    if (groupRank == 0)
    {
      sprintf(fileName, "%s_%06d", m_restartFileRoot.c_str(), cycleNum);
    }
    else
    {
      if (m_slaveDirectory.empty())
      {
        sprintf(fileName, "%s_%06d.%03d", m_restartFileRoot.c_str(), cycleNum, groupRank);
      }
      else
      {
        sprintf(fileName, "%s%s%s_%06d.%03d", m_slaveDirectory.c_str(), "/", m_restartFileRoot.c_str(), cycleNum, groupRank);
      }
    }

  }
  else
  {
    sprintf(baseFileName, "%s_%06d", m_fileRoot.c_str(), cycleNum);
    if (groupRank == 0)
      sprintf(fileName, "%s_%06d", m_fileRoot.c_str(), cycleNum);
    else if (m_slaveDirectory.empty())
    {
      sprintf(fileName, "%s_%06d.%03d", m_fileRoot.c_str(), cycleNum, groupRank);
    }
    else
    {
      sprintf(fileName, "%s%s%s_%06d.%03d", m_slaveDirectory.c_str(), "/", m_fileRoot.c_str(), cycleNum, groupRank);
    }
  }
  sprintf(dirName, "domain_%04d", domainNumber);

  m_dbFilePtr = static_cast<DBfile *>( PMPIO_WaitForBaton(m_baton, fileName, dirName) );

  m_fileName = fileName;
  m_baseFileName = baseFileName;
}


void SiloFile::WaitForBaton( const int domainNumber, const std::string& restartFileName )
{

  int rank = 0;
#if USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  const int groupRank = PMPIO_GroupRank(m_baton, rank);
  char fileName[200] = { 0 };
  char baseFileName[200] = { 0 };
  char dirName[200] = { 0 };


  sprintf(baseFileName, "%s", restartFileName.c_str());
  if (groupRank == 0)
    sprintf(fileName, "%s", restartFileName.c_str());
  else
  {
    if (m_slaveDirectory.empty())
    {
      sprintf(fileName, "%s.%03d", restartFileName.c_str(), groupRank);
    }
    else
    {
      sprintf(fileName, "%s%s%s.%03d", m_slaveDirectory.c_str(), "/", restartFileName.c_str(), groupRank);
    }

  }

  sprintf(dirName, "domain_%04d", domainNumber);

  m_dbFilePtr = (DBfile *) PMPIO_WaitForBaton(m_baton, fileName, dirName);

  m_fileName = fileName;
  m_baseFileName = baseFileName;
}
/**
 *
 */
void SiloFile::HandOffBaton()
{
  PMPIO_HandOffBaton(m_baton, m_dbFilePtr);
}

/**
 *
 * @param meshName
 * @param nnodes
 * @param coords
 * @param globalNodeNum
 * @param numRegions
 * @param shapecnt
 * @param meshConnectivity
 * @param globalElementNum
 * @param ghostFlag
 * @param shapetype
 * @param shapesize
 * @param cycleNumber
 * @param problemTime
 */
void SiloFile::WriteMeshObject(const std::string& meshName,
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
                               const realT problemTime)
{

  const DBdatatype datatype = DB_DOUBLE;


//  DBfacelist* facelist;
//  std::string facelistName;
//  facelistName = meshName + "_facelist";

  DBoptlist* optlist = DBMakeOptlist(4);
//  DBAddOption(optlist, DBOPT_NODENUM, const_cast<globalIndex*>
// (globalNodeNum));
  DBAddOption(optlist, DBOPT_CYCLE, const_cast<int*> (&cycleNumber));
  DBAddOption(optlist, DBOPT_DTIME, const_cast<realT*> (&problemTime));

  int numTotZones = 0;
  int lnodelist = 0;
  for (int i = 0 ; i < numRegions ; ++i)
  {
    numTotZones += shapecnt[i];
    //  if shapesize <= 0, that signals that we are using arbitrary polygons.
    if( shapesize[i] > 0 )
      lnodelist += shapecnt[i] * shapesize[i];
    else
      lnodelist += -shapesize[i];
  }


  if (numTotZones == 0)
  {
    char pwd[256];
    DBGetDir(m_dbFilePtr, pwd);
    std::string emptyObject = pwd;
    emptyObject += "/" + meshName;
    m_emptyMeshes.push_back(emptyObject);
  }
  else
  {

    std::string zonelistName;
    zonelistName = meshName + "_zonelist";

//  if (type_name<globalIndex>::name() == type_name<long long>::name())
//    DBAddOption(optlist, DBOPT_LLONGNZNUM, const_cast<int*> (&one));


    DBPutUcdmesh(m_dbFilePtr, meshName.c_str(), 3, nullptr, (float**) coords, nnodes, numTotZones,
                 zonelistName.c_str(), nullptr, datatype, optlist);

    DBClearOptlist(optlist);

    array<integer> nodelist(lnodelist);
    globalIndex_array globalZoneNumber(lnodelist);

    int count = 0;
    int elemCount = 0;
    for (int i = 0 ; i < numRegions ; ++i)
    {
      int n;
      if( shapesize[i] > 0 )
        n = shapecnt[i] * shapesize[i];
      else
        n = -shapesize[i];
      for (int j = 0 ; j < n ; ++j)
      {
        nodelist[count++] = meshConnectivity[i][j];
      }

      if( globalElementNum != nullptr )
      {
        for (int j = 0 ; j < shapecnt[i] ; ++j)
        {
          globalZoneNumber[elemCount++] = globalElementNum[i][j];
        }
        // write zonelist
        //      DBAddOption(optlist, DBOPT_ZONENUM, const_cast<globalIndex*>
        // (globalZoneNumber.data()));
        if ( std::is_same<globalIndex,long long>::value )
          DBAddOption(optlist, DBOPT_LLONGNZNUM, const_cast<int*> (&one));
      }
    }

    integer_array shapesize2(numRegions);
    for (int i = 0 ; i < numRegions ; ++i)
    {
      if( shapesize[i] < 0 )
      {
        shapesize2[i] = 0;
      }
      else
        shapesize2[i] = shapesize[i];
    }

    int hi_offset = 0;

    DBPutZonelist2( m_dbFilePtr, zonelistName.c_str(), numTotZones, 3, nodelist.data(), lnodelist, 0, 0,
                    hi_offset, const_cast<int*>(shapetype), const_cast<int*>(shapesize2.data()),
                    const_cast<int*>(shapecnt), numRegions,
                    optlist);

    DBClearOptlist(optlist);

    /*
       facelist = DBCalcExternalFacelist2(nodelist.data(), lnodelist,
                                         0, hi_offset,
                                         0, (int*) shapetype,
                                         (int*) shapesize2.data(), (int*)
                                            shapecnt, numRegions,
                                         nullptr, 0) ;

       DBPutFacelist(m_dbFilePtr, facelistName.c_str(), facelist->nfaces,
          facelist->ndims,
                              facelist->nodelist, facelist->lnodelist,
                              facelist->origin, facelist->zoneno,
                              facelist->shapesize, facelist->shapecnt,
                              1, nullptr, nullptr, 0);
     */
  }

  // write multimesh object
  int rank = 0;
#if USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  if (rank == 0)
  {
    DBAddOption(optlist, DBOPT_CYCLE, const_cast<int*> (&cycleNumber));
    DBAddOption(optlist, DBOPT_DTIME, const_cast<realT*> (&problemTime));

    WriteMultiXXXX(DB_UCDMESH, DBPutMultimesh, 0, meshName, cycleNumber, "/", optlist);
  }

  DBFreeOptlist(optlist);
}

// Arbitrary polygon. Have to deal with this separately
void SiloFile::WritePolygonMeshObject(const std::string& meshName,
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
                                      const int lnodelist)
{

  const DBdatatype datatype = DB_DOUBLE;


//  DBfacelist* facelist;
//  std::string facelistName;
//  facelistName = meshName + "_facelist";

  DBoptlist* optlist = DBMakeOptlist(4);
//  DBAddOption(optlist, DBOPT_NODENUM, const_cast<globalIndex*>
// (globalNodeNum));
  DBAddOption(optlist, DBOPT_CYCLE, const_cast<int*> (&cycleNumber));
  DBAddOption(optlist, DBOPT_DTIME, const_cast<realT*> (&problemTime));

  int numTotZones = shapecnt[0];
  if (numTotZones == 0)
  {
    char pwd[256];
    DBGetDir(m_dbFilePtr, pwd);
    std::string emptyObject = pwd;
    emptyObject += "/" + meshName;
    m_emptyMeshes.push_back(emptyObject);
  }
  else
  {
    std::string zonelistName;
    zonelistName = meshName + "_zonelist";


    DBPutUcdmesh(m_dbFilePtr, meshName.c_str(), 3, nullptr, (float**) coords, nnodes, numTotZones,
                 zonelistName.c_str(), nullptr, datatype, optlist);

    DBClearOptlist(optlist);

    array<integer> nodelist(lnodelist);
    globalIndex_array globalZoneNumber(lnodelist);

    int elemCount = 0;
    for (int j = 0 ; j < lnodelist ; ++j)
    {
      nodelist[j] = meshConnectivity[0][j];
    }

    {
      if( globalElementNum != nullptr )
      {
        for (int j = 0 ; j < shapecnt[0] ; ++j)
        {
          globalZoneNumber[elemCount++] = globalElementNum[0][j];
        }
        // write zonelist
        //      DBAddOption(optlist, DBOPT_ZONENUM, const_cast<globalIndex*>
        // (globalZoneNumber.data()));
        if( std::is_same<globalIndex,long long>::value )
          DBAddOption(optlist, DBOPT_LLONGNZNUM, const_cast<int*> (&one));
      }
    }



    int hi_offset = 0;

    DBPutZonelist2( m_dbFilePtr, zonelistName.c_str(), numTotZones, 3, nodelist.data(), lnodelist, 0, 0,
                    hi_offset, const_cast<int*>(shapetype), const_cast<int*>(shapesize), const_cast<int*>(shapecnt), numRegions,
                    optlist);

    DBClearOptlist(optlist);

    /*
       facelist = DBCalcExternalFacelist2(nodelist.data(), lnodelist,
                                         0, hi_offset,
                                         0, (int*) shapetype,
                                         (int*) shapesize2.data(), (int*)
                                            shapecnt, numRegions,
                                         nullptr, 0) ;

       DBPutFacelist(m_dbFilePtr, facelistName.c_str(), facelist->nfaces,
          facelist->ndims,
                              facelist->nodelist, facelist->lnodelist,
                              facelist->origin, facelist->zoneno,
                              facelist->shapesize, facelist->shapecnt,
                              1, nullptr, nullptr, 0);
     */
  }

  // write multimesh object
  int rank = 0;
#if USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  if (rank == 0)
  {
    DBAddOption(optlist, DBOPT_CYCLE, const_cast<int*> (&cycleNumber));
    DBAddOption(optlist, DBOPT_DTIME, const_cast<realT*> (&problemTime));

    WriteMultiXXXX(DB_UCDMESH, DBPutMultimesh, 0, meshName, cycleNumber, "/", optlist);
  }

  DBFreeOptlist(optlist);
}


/**
 * @author walsh24
 *
 * @param meshName
 * @param nX
 * @param nY
 * @param nZ
 * @param coords
 * @param cycleNumber
 * @param problemTime
 */
int SiloFile::WriteQuadMeshObject(const std::string& meshName,
                                  const localIndex nX,
                                  const localIndex nY,
                                  const localIndex nZ,
                                  realT* coords[3],
                                  const int cycleNumber,
                                  const realT problemTime)
{
  const DBdatatype datatype = DB_DOUBLE;

  DBoptlist* optlist = DBMakeOptlist(4);
//  DBAddOption(optlist, DBOPT_NODENUM, const_cast<globalIndex*>
// (globalNodeNum));
  DBAddOption(optlist, DBOPT_CYCLE, const_cast<int*> (&cycleNumber));
  DBAddOption(optlist, DBOPT_DTIME, const_cast<realT*> (&problemTime));

//  realT *coord[3];

  // record quad mesh dims - ugly but required to later write data to quad mesh
  // data file
  m_quadMeshDims[0] = nX;
  m_quadMeshDims[1] = nY;
  m_quadMeshDims[2] = nZ;
  m_quadMeshNDims = 3;


  int ret = DBPutQuadmesh(m_dbFilePtr, meshName.c_str(), nullptr,(float**) coords, m_quadMeshDims, 3,
                          datatype, DB_NONCOLLINEAR, optlist);
  if(ret == -1)
    return ret;

  DBClearOptlist(optlist);


  //----write multimesh object
  {
    int rank = 0;
  #if USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #endif
    if (rank == 0)
    {
      DBAddOption(optlist, DBOPT_CYCLE, const_cast<int*> (&cycleNumber));
      DBAddOption(optlist, DBOPT_DTIME, const_cast<realT*> (&problemTime));
      WriteMultiXXXX(DB_QUADMESH, DBPutMultimesh, 0, meshName.c_str(), cycleNumber, "/", optlist);
    }
  }

  //---free the option list
  DBFreeOptlist(optlist);

  return ret;
}

void SiloFile::WriteBeamMesh(const std::string& meshName,
                             const localIndex nnodes,
                             realT* coords[3],
                             const localIndex_array& node1,
                             const localIndex_array& node2,
                             const int cycleNumber,
                             const realT problemTime)
{
  // Connectivity.
  integer_array nodelist;
  {
    nodelist.reserve(2*node1.size());
    localIndex_array::const_iterator it2 = node2.begin();
    for (localIndex_array::const_iterator it = node1.begin() ;
         it != node1.end() ; ++it, ++it2)
    {
      nodelist.push_back(static_cast<int>(*it));
      nodelist.push_back(static_cast<int>(*it2));
    }
  }

  WriteBeamMesh( meshName, nnodes, coords, nodelist, cycleNumber, problemTime);
}


void SiloFile::WriteBeamMesh(const std::string& meshName,
                             const localIndex nnodes,
                             realT* coords[3],
                             const std::map<int, int>& connectivity,
                             const int cycleNumber,
                             const realT problemTime)
{
  // Connectivity.
  integer_array nodelist;
  {
    nodelist.reserve(2*connectivity.size());
    for (std::map<int,int>::const_iterator it = connectivity.begin() ;
         it != connectivity.end() ; ++it)
    {
      nodelist.push_back(it->first);
      nodelist.push_back(it->second);
    }
  }

  WriteBeamMesh( meshName, nnodes, coords, nodelist, cycleNumber, problemTime);
}

void SiloFile::WriteBeamMesh(const std::string& meshName,
                             const localIndex nnodes,
                             realT* coords[3],
                             integer_array& nodelist,
                             const int cycleNumber,
                             const realT problemTime)
{
  const DBdatatype datatype = DB_DOUBLE;

  DBoptlist* optlist = DBMakeOptlist(4);
  DBAddOption(optlist, DBOPT_CYCLE, const_cast<int*> (&cycleNumber));
  DBAddOption(optlist, DBOPT_DTIME, const_cast<realT*> (&problemTime));

  int lnodelist = nodelist.size();
  int nshapetypes = 1;
  int nzones = lnodelist/2;
  int shapecnt[] = {nzones};
  int shapesize[] = {2};
  int shapetype[] = {DB_ZONETYPE_BEAM};

  // Write out connectivity information.
  const int origin = 0;
  const int lo_offset = 0;
  const int hi_offset = 0;
  const int ndims = 3;
  const char* coordnames[3]  = { "xcoords", "ycoords", "zcoords" };

  // Write out connectivity information.
  DBPutZonelist2(m_dbFilePtr, "zonelist", nzones, ndims, lnodelist > 0 ? &nodelist[0] : nullptr,
                 lnodelist, origin, lo_offset, hi_offset,
                 shapetype, shapesize, shapecnt, nshapetypes, optlist);

  // Write an unstructured mesh.
  DBPutUcdmesh(m_dbFilePtr, meshName.c_str(), ndims, const_cast<char**>(coordnames), &coords[0], nnodes, nzones,
               "zonelist", nullptr, datatype, nullptr);
  DBClearOptlist(optlist);

  //----write multimesh object
  {
    int rank = 0;
  #if USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #endif
    if (rank == 0)
    {
      DBAddOption(optlist, DBOPT_CYCLE, const_cast<int*> (&cycleNumber));
      DBAddOption(optlist, DBOPT_DTIME, const_cast<realT*> (&problemTime));
      WriteMultiXXXX(DB_UCDMESH, DBPutMultimesh, 0, meshName.c_str(), cycleNumber, "/", optlist);
    }
  }

  //---free the option list
  DBFreeOptlist(optlist);
}

#if 0
void SiloFile::WriteArbitratryPolyhedralMeshObject( const std::string& meshName,
                                                    const localIndex nnodes,
                                                    realT* coords[3],
                                                    const globalIndex* globalNodeNum,
                                                    const int nDiscreteElements,
                                                    const int nfaces,
                                                    int* nodecnts,
                                                    const int sumnodecnts,
                                                    int* nodelist,
                                                    int* facecnts,
                                                    const int sumfacecnts,
                                                    int* facelist,
                                                    const globalIndex* const * const globalElementNum,
                                                    const int* ghostFlag,
                                                    const int cycleNumber,
                                                    const realT problemTime)
{
  //----create option list (done)
  DBoptlist* optlist = DBMakeOptlist(5);//allocate option list with 4 options

  //----create unstructured mesh (done)
  {
    {
//      DBAddOption(optlist, DBOPT_NODENUM, const_cast<globalIndex*>
// (globalNodeNum));//map local to global node index
      DBAddOption(optlist, DBOPT_CYCLE, const_cast<int*> (&cycleNumber));//cycle
                                                                         // number
                                                                         // in
                                                                         // the
                                                                         // problem
      DBAddOption(optlist, DBOPT_DTIME, const_cast<realT*> (&problemTime));//problem
                                                                           // time
                                                                           // value
//      if (type_name<globalIndex>::name() == type_name<long
// long>::name())//flag that the array for NODENUM is long long
//        DBAddOption(optlist, DBOPT_LLONGNZNUM, const_cast<int*> (&one));

      char temp[] = "dezonelist";
      DBAddOption(optlist, DBOPT_PHZONELIST, temp);//cycle number in the problem

    }

    //write a UCD mesh to the Silo file
    //UCD = unstructured cell data ... used to define any unstructured mesh
    DBPutUcdmesh(m_dbFilePtr, meshName.c_str(), nsdof, nullptr, (float**) coords, nnodes, nDiscreteElements,
                 nullptr, nullptr, DB_DOUBLE, optlist);
    DBClearOptlist(optlist);
  }

  //----create arbitrary polyhedral zone list
  {
    /*
       //get nodelist and globalZoneNumber
       ivector globalZoneNumber(nnodes);
       {
       int elemCount = 0;
       for (int i = 0; i < numRegions; ++i)
       {
        for (int j = 0; j < shapecnt[i] * shapesize[i]; ++j)
        {
          nodelist[count++] = meshConnectivity[i][j];
        }

        for (int j = 0; j < shapecnt[i]; ++j)
        {
          globalZoneNumber[elemCount++] = globalElementNum[i][j];
        }
       }
       }
     */

    //create option list
    /*
       DBAddOption(optlist, DBOPT_ZONENUM, const_cast<globalIndex*>
          (globalZoneNumber.data()));
       if (type_name<globalIndex>::name() == type_name<long long>::name())
       DBAddOption(optlist, DBOPT_LLONGNZNUM, const_cast<int*> (&one));
     */
    {
      //define all faces in the DE's as external faces
      char* extface = new char[nfaces];
      for(int i=0 ; i < nfaces ; i++)
        extface[i] = '1';

      DBPutPHZonelist(m_dbFilePtr, "dezonelist", nfaces, nodecnts, sumnodecnts, nodelist, extface, nDiscreteElements,
                      facecnts, sumfacecnts, facelist, 0, 0, nDiscreteElements-1, optlist);
      delete[] extface;
    }
    DBClearOptlist( optlist);
  }

  //----write multimesh object
  {
    int rank = 0;
  #if USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #endif
    if (rank == 0)
    {
      DBAddOption(optlist, DBOPT_CYCLE, const_cast<int*> (&cycleNumber));
      DBAddOption(optlist, DBOPT_DTIME, const_cast<realT*> (&problemTime));
      WriteMultiXXXX(DB_UCDMESH, DBPutMultimesh, 0, "demesh", cycleNumber, "/", optlist);
    }
  }

  //---free the option list
  DBFreeOptlist(optlist);
}



#endif



void SiloFile::WritePointMesh( const std::string& meshName,
                               const localIndex numPoints,
                               realT* coords[3],
                               const int cycleNumber,
                               const realT problemTime)
{
  const DBdatatype datatype = DB_DOUBLE;
  DBoptlist* optlist = DBMakeOptlist(2);
  DBAddOption(optlist, DBOPT_CYCLE, const_cast<int*> (&cycleNumber));
  DBAddOption(optlist, DBOPT_DTIME, const_cast<realT*> (&problemTime));
  DBPutPointmesh (m_dbFilePtr, meshName.c_str(), 3, coords, numPoints, datatype, optlist);


  //----write multimesh object
  {
    int rank = 0;
  #if USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  #endif
    if (rank == 0)
    {
      WriteMultiXXXX(DB_POINTMESH, DBPutMultimesh, 0, meshName.c_str(), cycleNumber, "/", optlist);
    }
  }

}

//void SiloFile::TestWriteDiscreteElementMeshObject()
//{
//  //----create option list (done)
//  DBoptlist* optlist = DBMakeOptlist(5);//allocate option list with 4 options
//  const char temp[] = "de_zonelist";
//
//  // define the problem
//  const std::string meshName = "de_mesh";
//  const localIndex nnodes = 4;
//  realT coords[nsdof][4];
//  const globalIndex globalNodeNum[] = {1,2,3,4};
//  const int nDiscreteElements = 1;
//  const int nfaces = 4;
//  const int nodecnts[] = {3,3,3,3};
//  const int sumnodecnts = 12;
//  const int nodelist[] = {0,2,3,1,3,2,0,1,2,1,0,3};
//  const int facecnts[] = {4};
//  const int sumfacecnts = 4;
//  const int facelist[] = {0,1,2,3};
//  const int cycleNumber = 0;
//  const realT problemTime = 0.;
//  {
//    coords[0][0] = -0.5;
//    coords[1][0] = -0.5;
//    coords[2][0] = -0.5;
//    coords[0][1] = 0.5;
//    coords[1][1] = -0.5;
//    coords[2][1] = -0.5;
//    coords[0][2] = 0;
//    coords[1][2] = -0.5;
//    coords[2][2] = 0.5;
//    coords[0][3] = 0;
//    coords[1][3] = 0;
//    coords[2][3] = 0.5;
//  }
//
//  //----create unstructured mesh (done)
//  {
//    {
//      DBAddOption(optlist, DBOPT_CYCLE, const_cast<int*>
// (&cycleNumber));//cycle number in the problem
//      DBAddOption(optlist, DBOPT_DTIME, const_cast<realT*>
// (&problemTime));//problem time value
//      DBAddOption(optlist, DBOPT_PHZONELIST, temp);//cycle number in the
// problem
//
//    }
//
//    //write a UCD mesh to the Silo file
//    //UCD = unstructured cell data ... used to define any unstructured mesh
//    DBPutUcdmesh(m_dbFilePtr, meshName.c_str(), nsdof, nullptr, (float**) coords,
// nnodes, nDiscreteElements,
//                 nullptr, nullptr, DB_DOUBLE, optlist);
//    DBClearOptlist(optlist);
//  }
//
//  //----create arbitrary polyhedral zone list
//  {
//    {
//      //define all faces in the DE's as external faces
//      char* extface = new char[nfaces];
//      for(int i=0; i < nfaces; i++)
//        extface[i] = '1';
//
//      DBPutPHZonelist(m_dbFilePtr, temp, nfaces, nodecnts, sumnodecnts,
// nodelist, extface, nDiscreteElements,
//                      facecnts, sumfacecnts, facelist, 0, 0,
// nDiscreteElements-1, optlist);
//      delete[] extface;
//    }
//    DBClearOptlist( optlist);
//  }
//
//  //----write multimesh object
//  {
//    int rank = 0;
//  #if USE_MPI
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//  #endif
//    if (rank == 0)
//    {
//      DBAddOption(optlist, DBOPT_CYCLE, const_cast<int*> (&cycleNumber));
//      DBAddOption(optlist, DBOPT_DTIME, const_cast<realT*> (&problemTime));
//      WriteMultiXXXX(DB_UCDMESH, DBPutMultimesh, 0, meshName.c_str(),
// cycleNumber, "/", optlist);
//    }
//  }
//
//  //---free the option list
//  DBFreeOptlist(optlist);
//}



void SiloFile::StopSiloCompilerWarnings()
{
  PMPIO_RankInGroup(m_baton, 0);

}

/**
 *
 * @param elementManager
 * @param cycleNumber
 * @param problemTime
 */
void SiloFile::WriteRegionSpecifications( ElementRegionManager const * const elementManager,
                                          ConstitutiveManager const * const constitutiveManager,
                                          const std::string& meshName,
                                          const int cycleNumber,
                                          const realT problemTime)
{

  auto const
  constitutiveMap = elementManager->
                    ConstructViewAccessor<
                    std::pair< Array2dT<localIndex>,Array2dT<localIndex> >
                    >( CellBlockSubRegion::viewKeyStruct::constitutiveMapString );

  std::string name = "Regions";
  int const nmat = constitutiveManager->GetSubGroups().size();
  array<int> matnos(nmat);
  std::vector<string> materialNameStrings(nmat);
  array<char const*> materialNames(nmat+1);
  materialNames.back() = nullptr;

  for( int matIndex=0 ; matIndex<nmat ; ++matIndex)
  {
    matnos[matIndex] = matIndex;
    materialNameStrings[matIndex] = constitutiveManager->GetGroup(matIndex)->getName();
    materialNames[matIndex] = materialNameStrings[matIndex].c_str();
  }

  int ndims = 1;
  int dims = elementManager->getNumberOfElements();

  array<integer> matlist(dims * nmat);


  int elemCount = 0;
  int regionCount = 0;
  for( localIndex er=0 ; er<elementManager->numRegions() ; ++er )
  {
    ElementRegion const * const elemRegion = elementManager->GetRegion(er);
    for( localIndex esr=0 ; esr<elemRegion->numSubRegions() ; ++esr )
    {
      CellBlockSubRegion const * const subRegion = elemRegion->GetSubRegion(esr);
      for (localIndex k = 0 ; k < subRegion->size() ; ++k)
      {
        // matIndex1 is the index of the material contained in the element
        localIndex const matIndex1 = constitutiveMap[er][esr].get().first[k][0];
        // matIndex2 is the index of the point within material specified in matIndex1
        localIndex const matIndex2 = constitutiveMap[er][esr].get().second[k][0];

        matlist[elemCount++] = matIndex1;
      }
    }
  }

  {
    DBoptlist* optlist = DBMakeOptlist(3);
    DBAddOption(optlist, DBOPT_MATNAMES, materialNames.data());
    DBAddOption(optlist, DBOPT_CYCLE, const_cast<int*> (&cycleNumber));
    DBAddOption(optlist, DBOPT_DTIME, const_cast<realT*> (&problemTime));

    DBPutMaterial( m_dbFilePtr,
                   name.c_str(),
                   meshName.c_str(),
                   nmat,
                   matnos.data(),
                   matlist.data(),
                   &dims,
                   ndims,
                   nullptr,
                   nullptr,
                   nullptr,
                   nullptr,
                   0,
                   DB_DOUBLE,
                   optlist);

    DBFreeOptlist(optlist);
  }
  // write multimesh object
  int rank = 0;
#if USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
  if (rank == 0)
  {

    int size = 1;
#if USE_MPI
    MPI_Comm_size(MPI_COMM_WORLD, &size);
#endif

    array<string> vBlockNames(size);
    std::vector<char*> BlockNames(size);
    char tempBuffer[1024];
    char currentDirectory[256];

    DBGetDir(m_dbFilePtr, currentDirectory);
    DBSetDir(m_dbFilePtr, "/");

    for (int i = 0 ; i < size ; ++i)
    {
      int groupRank = PMPIO_GroupRank(m_baton, i);

      if (groupRank == 0)
      {
        /* this mesh block is in the file 'root' owns */
        sprintf(tempBuffer, "/domain_%04d/%s", i, name.c_str());

      }
      else
      {
        /* this mesh block is another file */
        if (m_slaveDirectory.empty())
        {
          sprintf(tempBuffer, "%s.%03d:/domain_%04d/%s",
                  m_baseFileName.c_str(),
                  groupRank, i, name.c_str());
        }
        else
        {
          sprintf(tempBuffer, "%s%s%s.%03d:/domain_%04d/%s", m_slaveDirectory.c_str(), "/",
                  m_baseFileName.c_str(),
                  groupRank, i, name.c_str());
        }

      }
      vBlockNames[i] = tempBuffer;
      BlockNames[i] = const_cast<char*>( vBlockNames[i].c_str() );
    }

    {
      DBoptlist* optlist = DBMakeOptlist(5);
      DBAddOption(optlist, DBOPT_MATNAMES, materialNames.data());
      DBAddOption(optlist, DBOPT_CYCLE, const_cast<int*> (&cycleNumber));
      DBAddOption(optlist, DBOPT_DTIME, const_cast<realT*> (&problemTime));
      DBAddOption(optlist, DBOPT_NMATNOS, const_cast<int*>(&nmat) );
      DBAddOption(optlist, DBOPT_MATNOS, matnos.data() );

      DBPutMultimat(m_dbFilePtr, name.c_str(), size, BlockNames.data(),
                    const_cast<DBoptlist*> (optlist));
      DBFreeOptlist(optlist);

    }

    DBSetDir(m_dbFilePtr, currentDirectory);

  }


}

void SiloFile::ClearEmptiesFromMultiObjects(const int cycleNum)
{

  int size = 1;
  int rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::string sendbufferVars;
  std::string sendbufferMesh;

  if( rank != 0 )
  {
    for( array<string>::const_iterator emptyObject=m_emptyVariables.begin() ;
         emptyObject!=m_emptyVariables.end() ; ++emptyObject )
    {
      sendbufferVars += *emptyObject + ' ';
    }

    for( array<string>::const_iterator emptyObject=m_emptyMeshes.begin() ;
         emptyObject!=m_emptyMeshes.end() ; ++emptyObject )
    {
      sendbufferMesh += *emptyObject + ' ';
    }

  }

  int sizeOfSendBufferVars = sendbufferVars.size();
  int sizeOfSendBufferMesh = sendbufferMesh.size();

  integer_array rcounts(size);
  integer_array displs(size);
  MPI_Gather( &sizeOfSendBufferVars, 1, MPI_INT, rcounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

  int sizeOfReceiveBuffer = 0;
  displs[0] = 0;
  for ( int i=1 ; i<size ; ++i)
  {
    displs[i] = displs[i-1]+rcounts[i-1];
    sizeOfReceiveBuffer += rcounts[i];
  }
  std::string receiveBufferVars(sizeOfReceiveBuffer,'\0');

  MPI_Gatherv ( &sendbufferVars[0], sizeOfSendBufferVars, MPI_CHAR,
                &receiveBufferVars[0], rcounts.data(), displs.data(),
                MPI_CHAR, 0, MPI_COMM_WORLD );


  MPI_Gather( &sizeOfSendBufferMesh, 1, MPI_INT, rcounts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

  int sizeOfReceiveBufferMesh = 0;
  displs[0] = 0;
  for ( int i=1 ; i<size ; ++i)
  {
    displs[i] = displs[i-1]+rcounts[i-1];
    sizeOfReceiveBufferMesh += rcounts[i];
  }
  std::string receiveBufferMesh(sizeOfReceiveBufferMesh,'\0');

  MPI_Gatherv ( &sendbufferMesh[0], sizeOfSendBufferMesh, MPI_CHAR,
                &receiveBufferMesh[0], rcounts.data(), displs.data(),
                MPI_CHAR, 0, MPI_COMM_WORLD );



  if( rank== 0 )
  {
    std::istringstream iss(receiveBufferVars);
    copy(std::istream_iterator<std::string>(iss),
         std::istream_iterator<std::string>(),
         std::back_inserter< array<string> >(m_emptyVariables));

    std::istringstream issm(receiveBufferMesh);
    copy(std::istream_iterator<std::string>(issm),
         std::istream_iterator<std::string>(),
         std::back_inserter< array<string> >(m_emptyMeshes));
  }

  if (rank == 0)
  {
    DBfile *siloFile = DBOpen(m_baseFileName.c_str(), DB_UNKNOWN, DB_APPEND);
    std::string empty("EMPTY");

    for (array<string>::iterator emptyObject = m_emptyVariables.begin() ; emptyObject
         != m_emptyVariables.end() ; ++emptyObject)
    {
      size_t pathBegin = emptyObject->find_first_of('/', 1);
      size_t pathEnd = emptyObject->find_last_of('/');
      std::string domainString(*emptyObject, 1, pathBegin);
      std::string pathToMultiObj(*emptyObject, pathBegin, pathEnd - pathBegin);
      std::string varName(*emptyObject, pathEnd + 1);

      DBSetDir(siloFile, pathToMultiObj.c_str());

      DBmultivar* multiVar = DBGetMultivar(siloFile, varName.c_str());

      if( multiVar != nullptr )
      {
        array<const char*> newvarnames(multiVar->nvars);

        for (int i = 0 ; i < multiVar->nvars ; ++i)
        {

          std::string path(multiVar->varnames[i]);
          if (path.find(domainString) != std::string::npos)
          {
            newvarnames(i) = empty.c_str();
          }
          else
          {
            newvarnames(i) = multiVar->varnames[i];
          }
        }

        DBoptlist *optlist = DBMakeOptlist(5);
        DBAddOption(optlist, DBOPT_CYCLE, const_cast<int*> (&cycleNum));
        //      DBAddOption( optlist, DBOPT_DTIME,
        // const_cast<realT*>(&problemTime) );
        DBAddOption(optlist, DBOPT_REGION_PNAMES, multiVar->region_pnames);
        DBAddOption(optlist, DBOPT_TENSOR_RANK, &multiVar->tensor_rank);
        DBAddOption(optlist, DBOPT_MMESH_NAME, multiVar->mmesh_name);

        DBPutMultivar(siloFile, varName.c_str(), multiVar->nvars,
                      const_cast<char**> (newvarnames.data()), multiVar->vartypes, optlist);
        DBFreeOptlist(optlist);
      }
    }


    for (array<string>::iterator emptyObject = m_emptyMeshes.begin() ; emptyObject
         != m_emptyMeshes.end() ; ++emptyObject)
    {
      size_t pathBegin = emptyObject->find_first_of('/', 1);
      size_t pathEnd = emptyObject->find_last_of('/');
      std::string domainString(*emptyObject, 1, pathBegin);
      std::string pathToMultiObj(*emptyObject, pathBegin, pathEnd - pathBegin);
      std::string varName(*emptyObject, pathEnd + 1);

      if( !(pathToMultiObj.compare("")) )
      {
        pathToMultiObj = "/";
      }

      DBSetDir(siloFile, pathToMultiObj.c_str());

      DBmultimesh* multiMesh = DBGetMultimesh(siloFile, varName.c_str());

      if( multiMesh != nullptr )
      {
        array<const char*> newmeshnames(multiMesh->nblocks);

        for (int i = 0 ; i < multiMesh->nblocks ; ++i)
        {

          std::string path(multiMesh->meshnames[i]);
          if (path.find(domainString) != std::string::npos)
          {
            newmeshnames(i) = empty.c_str();
          }
          else
          {
            newmeshnames(i) = multiMesh->meshnames[i];
          }
        }

        DBoptlist *optlist = DBMakeOptlist(2);
        DBAddOption( optlist, DBOPT_CYCLE, const_cast<int*> (&cycleNum));
//        DBAddOption( optlist, DBOPT_DTIME, const_cast<realT*>(&problemTime) );

        DBPutMultimesh( siloFile, varName.c_str(), multiMesh->nblocks,
                        const_cast<char**> (newmeshnames.data()), multiMesh->meshtypes, optlist);
        DBFreeOptlist(optlist);
      }
    }
    DBClose(siloFile);
  }
  m_emptyVariables.clear();
  m_emptyMeshes.clear();

}



template<typename TYPE>
void SiloFile::DBWriteWrapper( const std::string& name, const TYPE& data )
{
  int dims = SiloFileUtilities::GetNumberOfVariablesInField<TYPE>();
  TYPE temp = data;
  DBWrite( m_dbFilePtr, name.c_str(), &temp, &dims, 1, SiloFileUtilities::DB_TYPE<TYPE>() );
}
template void SiloFile::DBWriteWrapper( const std::string&, const int& );
template void SiloFile::DBWriteWrapper( const std::string&, const unsigned int& );
template void SiloFile::DBWriteWrapper( const std::string&, const double& );
//template void SiloFile::DBWriteWrapper( const std::string&, const localIndex&
// );
template void SiloFile::DBWriteWrapper( const std::string&, const globalIndex& );
template void SiloFile::DBWriteWrapper( const std::string&, const R1Tensor& );
template<>
void SiloFile::DBWriteWrapper( const std::string& name, const std::string& data )
{
  int dims = data.size();
  std::string temp = data;
  DBWrite( m_dbFilePtr, name.c_str(), const_cast<char*>(temp.data()), &dims, 1, DB_CHAR );
}


void SiloFile::DBWriteWrapper( const std::string& name, const array<std::string>& data )
{

  if( !data.empty() )
  {
    std::string temp;

    for( array<std::string>::const_iterator i=data.begin() ; i!=data.end() ; ++i )
    {
      temp += *i + ' ';

    }

    int dims = temp.size();
    DBWrite( m_dbFilePtr, name.c_str(), const_cast<char*>(temp.data()), &dims, 1, DB_CHAR );
  }
}


template<typename TYPE>
void SiloFile::DBWriteWrapper( const std::string& name, const set<TYPE>& data )
{
  if( !data.empty() )
  {
    array<TYPE> dataArray;
    dataArray.resize( integer_conversion<localIndex>(data.size()) );
    std::copy( data.begin(), data.end(), dataArray.begin() );
    DBWriteWrapper(name,dataArray);
  }
}
template void SiloFile::DBWriteWrapper( const std::string&, const lSet& );
//template void SiloFile::DBWriteWrapper( const std::string&, const iSet& );


template<typename TYPE>
void SiloFile::DBWriteWrapper( const std::string& name, const Array2dT<TYPE>& data )
{
  if( !data.empty() )
  {
    int dims[3];
    dims[0] = SiloFileUtilities::GetNumberOfVariablesInField<TYPE>();
    dims[1] = static_cast<int>(data.size(1));
    dims[2] = static_cast<int>(data.size(0));

    DBWrite( m_dbFilePtr, name.c_str(), const_cast<TYPE*>(data.data()), dims,
             3, SiloFileUtilities::DB_TYPE<TYPE>() );
  }
}
template void SiloFile::DBWriteWrapper( const std::string&, const Array2dT<real64>& );
template void SiloFile::DBWriteWrapper( const std::string&, const Array2dT<R1Tensor>& );
template void SiloFile::DBWriteWrapper( const std::string&, const Array2dT<R2Tensor>& );

template<typename TYPE>
void SiloFile::DBWriteWrapper( const std::string& name, const array<array<TYPE> >& data )
{

  const std::string sizeName = name+"_sizes";
  const std::string dataName = name+"_data";


  array<typename array<TYPE>::size_type> dataSizes;
  array<TYPE> dataSerial;

  for( typename array<array<TYPE> >::const_iterator i=data.begin() ; i!=data.end() ; ++i )
  {
    dataSizes.push_back(i->size());
    dataSerial.insert( dataSerial.end(), i->begin(), i->end() );
  }

  if( !(dataSerial.empty()) )
  {
    DBWriteWrapper( sizeName, dataSizes );
    DBWriteWrapper( dataName, dataSerial );
  }

}
template void SiloFile::DBWriteWrapper( const std::string&, const array<array<realT> >& );
template void SiloFile::DBWriteWrapper( const std::string&, const array<array<R2Tensor> >& );
template void SiloFile::DBWriteWrapper( const std::string&, const array<array<int> >& );


template<typename TYPE>
void SiloFile::DBWriteWrapper( const std::string& name, const array<Array2dT<TYPE> >& data )
{

  const std::string sizeName0 = name+"_sizes0";
  const std::string sizeName1 = name+"_sizes1";
  const std::string dataName = name+"_data";


  array<typename array<Array2dT<TYPE> >::size_type> dataSizes0;
  array<typename Array2dT<TYPE>::size_type> dataSizes1;
  array<TYPE> dataSerial;

  for( typename array<Array2dT<TYPE> >::const_iterator i=data.begin() ; i!=data.end() ; ++i )
  {
    dataSizes0.push_back(i->size(0));
    dataSizes1.push_back(i->size(1));
    dataSerial.insert( dataSerial.end(), i->begin(), i->end() );
  }

  if( !(dataSerial.empty()) )
  {
    DBWriteWrapper( sizeName0, dataSizes0 );
    DBWriteWrapper( sizeName1, dataSizes1 );
    DBWriteWrapper( dataName, dataSerial );
  }

}
template void SiloFile::DBWriteWrapper( const std::string&, const array<Array2dT<R1Tensor> >& );


template<typename TYPE>
void SiloFile::DBWriteWrapper( const std::string& name, const array<array<array<TYPE> > >& data )
{
  const std::string sizeName0 = name+"_sizes0";
  const std::string sizeName1 = name+"_sizes1";
  const std::string dataName = name+"_data";


  array<typename array<array<TYPE> >::size_type> dataSizes0;
  array<typename array<TYPE>::size_type> dataSizes1;
  array<TYPE> dataSerial;

  for( typename array<array<array<TYPE> > >::const_iterator i=data.begin() ; i!=data.end() ; ++i )
  {
    dataSizes0.push_back(i->size());

    for( typename array<array<TYPE> >::const_iterator j=i->begin() ; j!=i->end() ; ++j )
    {
      dataSizes1.push_back(j->size());
      dataSerial.insert( dataSerial.end(), j->begin(), j->end() );
    }
  }

  if( !(dataSerial.empty()) )
  {
    DBWriteWrapper( sizeName0, dataSizes0 );
    DBWriteWrapper( sizeName1, dataSizes1 );
    DBWriteWrapper( dataName, dataSerial );
  }

}
template void SiloFile::DBWriteWrapper( const std::string& name, const array<array<array<R1Tensor> > >& data );



template<typename TYPE>
void SiloFile::DBWriteWrapper( const std::string& name, const array<set<TYPE> >& data )
{
  const std::string sizeName = name+"_sizes";
  const std::string dataName = name+"_data";


  array<typename array<TYPE>::size_type> dataSizes;
  array<TYPE> dataSerial;

  for( typename array<set<TYPE> >::const_iterator i=data.begin() ; i!=data.end() ; ++i )
  {
    dataSizes.push_back(i->size());
    dataSerial.insert( dataSerial.end(), i->begin(), i->end() );
  }
  if( !(dataSerial.empty()) )
  {
    DBWriteWrapper( sizeName, dataSizes );
    DBWriteWrapper( dataName, dataSerial );
  }
}

template< typename T1, typename T2 >
void SiloFile::DBWriteWrapper( const std::string& name, const std::map< T1, T2 >& datamap )
{
  if( !(datamap.empty()) )
  {
    const std::string keyName = name+"_keyval";
    const std::string mapName = name+"_mapval";

    array<T1> keyArray;
    keyArray.reserve(datamap.size());

    array<T2> mapArray;
    mapArray.reserve(datamap.size());


    for( typename std::map< T1, T2 >::const_iterator i=datamap.begin() ; i!=datamap.end() ; ++i )
    {
      keyArray.push_back(i->first);
      mapArray.push_back(i->second);
    }

    DBWriteWrapper( keyName, keyArray);
    DBWriteWrapper( mapName, mapArray);
  }
}
template void SiloFile::DBWriteWrapper( const std::string&, const std::map< globalIndex, localIndex >& );


template<typename TYPE>
void SiloFile::DBWriteWrapper( const std::string& subdir, const std::map< std::string, TYPE>& member )
{

  DBMkDir( m_dbFilePtr, subdir.c_str() );
  DBSetDir( m_dbFilePtr, subdir.c_str() );

  array<string> memberNames;


  // iterate over all entries in the member map
  for( typename std::map< std::string, TYPE >::const_iterator i = member.begin() ; i!=member.end() ; ++i )
  {
    // the field name is the key to the map
    const std::string fieldName = i->first;

    memberNames.push_back( fieldName );

    // check to see if the field should be written
//    std::map<std::string, FieldBase*>::const_iterator fieldAttributes =
// FieldInfo::AttributesByName.find(fieldName);

    bool writeField = true;
//    if( fieldAttributes == FieldInfo::AttributesByName.end() )
//    {
//      writeField = true;
//    }
//    else if( fieldAttributes->second->m_WriteToRestart )
//    {
//      writeField = true;
//    }

    if( writeField )
    {
      // the field data is mapped value
      const TYPE& fieldData = i->second;

      DBWriteWrapper( fieldName, fieldData );
    }
  }

  DBWriteWrapper( "memberNames", memberNames );

  DBSetDir( m_dbFilePtr, ".." );

}
template void SiloFile::DBWriteWrapper( const std::string&, const std::map< std::string, array<int> >& );
template void SiloFile::DBWriteWrapper( const std::string&, const std::map< std::string, array<long int> >& );
template void SiloFile::DBWriteWrapper( const std::string&, const std::map< std::string, array<long long int> >& );
template void SiloFile::DBWriteWrapper( const std::string&, const std::map< std::string, lArray2d>& );
template void SiloFile::DBWriteWrapper( const std::string&, const std::map< std::string, array<localIndex_array> >& );
template void SiloFile::DBWriteWrapper( const std::string&, const std::map< std::string, lSet>& );
template void SiloFile::DBWriteWrapper( const std::string&, const std::map< std::string, array<string> >& );



template<typename TYPE>
void SiloFile::DBWriteWrapper( const std::string& subdir, const std::map< std::string, InterObjectRelation<TYPE> >& member )
{
  DBMkDir( m_dbFilePtr, subdir.c_str() );
  DBSetDir( m_dbFilePtr, subdir.c_str() );

  // iterate over all entries in the member map
  for( typename std::map< std::string, InterObjectRelation<TYPE> >::const_iterator i = member.begin() ; i!=member.end() ; ++i )
  {
    // the field name is the key to the map
    const std::string fieldName = i->first;

    // the field data is mapped value
    const TYPE& fieldData = i->second;

    DBWriteWrapper( fieldName, fieldData );
  }

  DBSetDir( m_dbFilePtr, ".." );

}
template void SiloFile::DBWriteWrapper( const std::string&, const std::map< std::string, OneToOneRelation>& );
template void SiloFile::DBWriteWrapper( const std::string&, const std::map< std::string, FixedOneToManyRelation>& );
template void SiloFile::DBWriteWrapper( const std::string&, const std::map< std::string, OrderedVariableOneToManyRelation>& );
template void SiloFile::DBWriteWrapper( const std::string&, const std::map< std::string, UnorderedVariableOneToManyRelation>& );



template<typename TYPE>
void SiloFile::DBWriteWrapper( const std::string& name, const TYPE* const data, const int size )
{
  if( size!=0 )
  {
    int dims[2];
    dims[0] = SiloFileUtilities::GetNumberOfVariablesInField<TYPE>();
    dims[1] = size;

    DBWrite( m_dbFilePtr, name.c_str(), const_cast<TYPE*>(data), dims,
             2, SiloFileUtilities::DB_TYPE<TYPE>() );
  }
}
template void SiloFile::DBWriteWrapper( const std::string&, const int* const, const int );



template<typename TYPE>
void SiloFile::DBReadWrapper( const std::string& name, TYPE& data ) const
{
  if( DBInqVarExists( m_dbFilePtr, name.c_str() ) )
  {

    if( SiloFileUtilities::DB_TYPE<TYPE>() != DBGetVarType( m_dbFilePtr, name.c_str() ) )
    {
      GEOS_ERROR("SiloFile::DBReadWrapper: variable "+ name +" is incorrect type" );
    }
    if( sizeof(TYPE) != DBGetVarByteLength( m_dbFilePtr, name.c_str() ) )
    {
      GEOS_ERROR("SiloFile::DBReadWrapper: variable "+ name +" is incorrect size" );
    }

    DBReadVar( m_dbFilePtr, name.c_str(), &data );
  }
  else
  {
    //GEOS_ERROR("SiloFile::DBReadWrapper: variable "+ name +" does not exist in
    // silo file" );
  }
}
template void SiloFile::DBReadWrapper( const std::string&, int& ) const;
template void SiloFile::DBReadWrapper( const std::string&, unsigned int& ) const;
template void SiloFile::DBReadWrapper( const std::string&, double& ) const;
//template void SiloFile::DBReadWrapper( const std::string&, localIndex& )
// const;
template void SiloFile::DBReadWrapper( const std::string&, globalIndex& ) const;
template void SiloFile::DBReadWrapper( const std::string&, R1Tensor& ) const;
template<>
void SiloFile::DBReadWrapper( const std::string& name, std::string& data ) const
{
  if( DBInqVarExists( m_dbFilePtr, name.c_str() ) )
  {

    if( DB_CHAR != DBGetVarType( m_dbFilePtr, name.c_str() ) )
    {
      GEOS_ERROR("SiloFile::DBReadWrapper: variable "+ name +" is incorrect type" );
    }

    data.resize( DBGetVarByteLength( m_dbFilePtr, name.c_str() ) );
    DBReadVar( m_dbFilePtr, name.c_str(), const_cast<char*>( data.data() ) );
  }
  else
  {
    //GEOS_ERROR("SiloFile::DBReadWrapper: variable "+ name +" does not exist in
    // silo file" );
  }
}

template<typename TYPE>
void SiloFile::DBReadWrapper( const std::string& name, array<TYPE>& data ) const
{


  if( DBInqVarExists( m_dbFilePtr, name.c_str() ) )
  {
    if( SiloFileUtilities::DB_TYPE<TYPE>() != DBGetVarType( m_dbFilePtr, name.c_str() ) )
    {
      GEOS_ERROR("SiloFile::DBReadWrapper: variable "+ name +" is incorrect type" );
    }

    if( 1 )
    {
      int dims[2];
      DBGetVarDims( m_dbFilePtr, name.c_str(), 2, dims );

      data.resize( dims[1] );
    }
    else if( static_cast<int>(sizeof(TYPE)*data.size()) != DBGetVarByteLength( m_dbFilePtr, name.c_str() ) )
    {
      GEOS_ERROR("SiloFile::DBReadWrapper: variable "+ name +" is incorrect size" );
    }

    DBReadVar( m_dbFilePtr, name.c_str(), data.data() );
  }
  else
  {
    //GEOS_ERROR("SiloFile::DBReadWrapper: variable "+ name +" does not exist in
    // silo file" );
  }
}
template void SiloFile::DBReadWrapper( const std::string&, array<int>& ) const;
template void SiloFile::DBReadWrapper( const std::string&, array<long int>& ) const;
template void SiloFile::DBReadWrapper( const std::string&, array<long long int>& ) const;
//template void SiloFile::DBReadWrapper( const std::string&, integer_array& )
// const;
template void SiloFile::DBReadWrapper( const std::string&, array<real64>& ) const;
template void SiloFile::DBReadWrapper( const std::string&, array<R1Tensor>& ) const;
template void SiloFile::DBReadWrapper( const std::string&, array<R2Tensor>& ) const;
template void SiloFile::DBReadWrapper( const std::string&, array<R2SymTensor>& ) const;

template<>
void SiloFile::DBReadWrapper( const std::string& name, array<std::string>& data ) const
{

  if( DBInqVarExists( m_dbFilePtr, name.c_str() ) )
  {

    array<char> temp( DBGetVarByteLength( m_dbFilePtr, name.c_str() ) );

    DBReadVar( m_dbFilePtr, name.c_str(), temp.data() );

    std::string tmpStr;
    for( array<char>::const_iterator i=temp.begin() ; i!=temp.end() ; ++i )
    {
      if( *i != ' ')
      {
        tmpStr.push_back(*i);
      }
      else
      {
        data.push_back(tmpStr);
        tmpStr.clear();
      }
    }

  }
  else
  {
    //GEOS_ERROR("SiloFile::DBReadWrapper: variable "+ name +" does not exist in
    // silo file" );
  }
}



template<typename TYPE>
void SiloFile::DBReadWrapper( const std::string& name, set<TYPE>& data ) const
{

  if( DBInqVarExists( m_dbFilePtr, name.c_str() ) )
  {
    array<TYPE> dataArray;

    int dims[2];
    DBGetVarDims( m_dbFilePtr, name.c_str(), 2, dims );
    dataArray.resize( dims[1] );

    DBReadWrapper(name,dataArray);

    data.clear();
    data.insert( dataArray.begin(), dataArray.end() );
  }

}
template void SiloFile::DBReadWrapper( const std::string&, lSet& ) const;
//template void SiloFile::DBReadWrapper( const std::string&, iSet& ) const;


template<typename TYPE>
void SiloFile::DBReadWrapper( const std::string& name, Array2dT<TYPE>& data ) const
{

  if( DBInqVarExists( m_dbFilePtr, name.c_str() ) )
  {
    if( SiloFileUtilities::DB_TYPE<TYPE>() != DBGetVarType( m_dbFilePtr, name.c_str() ) )
    {
      GEOS_ERROR("SiloFile::DBReadWrapper: variable "+ name +" is incorrect type" );
    }
    if( static_cast<int>(sizeof(TYPE)*data.size()) != DBGetVarByteLength( m_dbFilePtr, name.c_str() ) )
    {
      GEOS_ERROR("SiloFile::DBReadWrapper: variable "+ name +" is incorrect size" );
    }

    int dims[3];
    DBGetVarDims( m_dbFilePtr, name.c_str(), 3, dims);
    if( dims[0] != SiloFileUtilities::GetNumberOfVariablesInField<TYPE>() ||
        dims[1] != integer_conversion<int>(data.size(1)) ||
        dims[2] != integer_conversion<int>(data.size(0)) )
    {
      GEOS_ERROR("SiloFile::DBReadWrapper: variable "+ name +" dimensions are incorrect" );
    }

    DBReadVar( m_dbFilePtr, name.c_str(), data.data() );
  }
  else
  {
    //GEOS_ERROR("SiloFile::DBReadWrapper: variable "+ name +" does not exist in
    // silo file" );
  }
}
template void SiloFile::DBReadWrapper( const std::string& name, Array2dT<real64>& data ) const;
template void SiloFile::DBReadWrapper( const std::string& name, Array2dT<R1Tensor>& data ) const;
template void SiloFile::DBReadWrapper( const std::string& name, Array2dT<R2Tensor>& data ) const;


template<typename TYPE>
void SiloFile::DBReadWrapper( const std::string& name, array<array<TYPE> >& data ) const
{

  const std::string sizeName = name+"_sizes";
  const std::string dataName = name+"_data";


  array<typename array<TYPE>::size_type> dataSizes(data.size());

  DBReadWrapper( sizeName, dataSizes );

  typename array<TYPE>::size_type serialSize = 0;
  for( integer i=0 ; i<data.size() ; ++i )
  {
    data[i].resize( dataSizes[i]);
    serialSize += dataSizes[i];
  }

  if( serialSize > 0 )
  {
    array<TYPE> dataSerial;
    dataSerial.resize(serialSize);
    DBReadWrapper( dataName, dataSerial );

    typename array<TYPE>::const_iterator idataSerial = dataSerial.begin();

    for( integer i=0 ; i<data.size() ; ++i )
    {
      for( integer j=0 ; j<data[i].size() ; ++j, ++idataSerial )
      {
        data[i][j] = *idataSerial;
      }
    }
  }

}
template void SiloFile::DBReadWrapper( const std::string&, array<array<int> >& ) const;
template void SiloFile::DBReadWrapper( const std::string&, array<array<realT> >& ) const;
template void SiloFile::DBReadWrapper( const std::string&, array<array<R2Tensor> >& ) const;


template<typename TYPE>
void SiloFile::DBReadWrapper( const std::string& name, array<Array2dT<TYPE> >& data ) const
{
  const std::string sizeName0 = name+"_sizes0";
  const std::string sizeName1 = name+"_sizes1";
  const std::string dataName = name+"_data";


  array<typename array<Array2dT<TYPE> >::size_type> dataSizes0(data.size());
  array<typename Array2dT<TYPE>::size_type> dataSizes1(data.size());

  DBReadWrapper( sizeName0, dataSizes0 );
  DBReadWrapper( sizeName1, dataSizes1 );


  typename array<TYPE>::size_type serialSize = 0;
  for( integer i=0 ; i<data.size() ; ++i )
  {
    data[i].resize(dataSizes0[i],dataSizes1[i]);
    serialSize += dataSizes0[i] * dataSizes1[i];
  }

  if( serialSize > 0 )
  {
    array<TYPE> dataSerial(serialSize);
    DBReadWrapper( dataName, dataSerial );

    typename array<TYPE>::const_iterator idataSerial = dataSerial.begin();

    for( integer i=0 ; i< static_cast<integer>(data.size()) ; ++i )
    {
      for( integer j=0 ; j< static_cast<integer>(data[i].size(0)) ; ++j )
      {
        for( integer k=0 ; k< static_cast<integer>(data[i].size(1)); ++k )
        {
          data[i][j][k] = *idataSerial;
          ++idataSerial;
        }
      }


    }
  }
}
template void SiloFile::DBReadWrapper( const std::string&, array<Array2dT<R1Tensor> >& ) const;



template<typename TYPE>
void SiloFile::DBReadWrapper( const std::string& name, array<array<array<TYPE> > >& data ) const
{
  const std::string sizeName0 = name+"_sizes0";
  const std::string sizeName1 = name+"_sizes1";
  const std::string dataName = name+"_data";


  array<typename array<array<TYPE> >::size_type> dataSizes0(data.size());
  array<typename array<TYPE>::size_type> dataSizes1;

  DBReadWrapper( sizeName0, dataSizes0 );

  typename array<TYPE>::size_type sizeSize1 = 0;
  for( integer i=0 ; i<data.size() ; ++i )
  {
    data[i].resize( dataSizes0[i]);
    sizeSize1 += dataSizes0[i];
  }

  dataSizes1.resize( sizeSize1 );
  DBReadWrapper( sizeName1, dataSizes1 );

  typename array<TYPE>::size_type serialSize = 0;
  typename array<TYPE>::size_type size1count = 0;
  for( integer i=0 ; i<data.size() ; ++i )
  {
    for( integer j=0 ; j<data[i].size() ; ++j )
    {
      data[i][j].resize( dataSizes1[size1count]);
      serialSize += dataSizes1[size1count];
      ++size1count;
    }
  }

  if( serialSize > 0 )
  {
    array<TYPE> dataSerial(serialSize);
    DBReadWrapper( dataName, dataSerial );

    typename array<TYPE>::const_iterator idataSerial = dataSerial.begin();

    for( integer i=0 ; i<data.size() ; ++i )
    {
      for( integer j=0 ; j<data[i].size() ; ++j )
      {
        for( integer k=0 ; k<data[i][j].size() ; ++k, ++idataSerial )
        {
          data[i][j][k] = *idataSerial;
        }
      }
    }
  }
}
template void SiloFile::DBReadWrapper( const std::string& name, array<array<array<R1Tensor> > >& data ) const;

template<typename TYPE>
void SiloFile::DBReadWrapper( const std::string& name, array<set<TYPE> >& data ) const
{

  array<array<TYPE> > vdata(data.size());

  DBReadWrapper( name, vdata );

  for( integer i=0 ; i<data.size() ; ++i )
  {
    data[i].clear();
    data[i].insert( vdata[i].begin(), vdata[i].end() );
  }
}

template< typename T1, typename T2 >
void SiloFile::DBReadWrapper( const std::string& name, std::map< T1, T2 >& datamap ) const
{
  const std::string keyName = name+"_keyval";
  const std::string mapName = name+"_mapval";

  if( DBInqVarExists( m_dbFilePtr, keyName.c_str() ) )
  {

    unsigned int size = DBGetVarLength ( m_dbFilePtr, keyName.c_str() );

    array<T1> keyArray;
    keyArray.resize(size);

    array<T2> mapArray;
    mapArray.resize(size);


    DBReadWrapper( keyName, keyArray);
    DBReadWrapper( mapName, mapArray);

    datamap.clear();


    for( unsigned int i=0 ; i<size ; ++i  )
    {
      datamap[keyArray[i]] = mapArray[i];
    }

  }
}
template void SiloFile::DBReadWrapper( const std::string&, std::map< globalIndex, localIndex >& ) const;



template<typename TYPE>
void SiloFile::DBReadWrapper( const std::string& subdir, std::map< std::string, TYPE>& member ) const
{

  DBSetDir( m_dbFilePtr, subdir.c_str() );

  array<string> memberNames;

  DBReadWrapper( "memberNames", memberNames );

  for( array<string>::const_iterator i=memberNames.begin() ; i!=memberNames.end() ; ++i )
  {
    const std::string fieldName = *i;
    // check to see if the field should be written
//    std::map<std::string, FieldBase*>::const_iterator fieldAttributes =
// FieldInfo::AttributesByName.find(fieldName);

    bool writeField = true;
//    if( fieldAttributes == FieldInfo::AttributesByName.end() )
//    {
//      writeField = true;
//    }
//    else if( fieldAttributes->second->m_WriteToRestart )
//    {
//      writeField = true;
//    }

    if( writeField )
    {
      // the field data is mapped value
      TYPE& fieldData = member[fieldName];

      DBReadWrapper( fieldName, fieldData );
    }
  }

  /*
     // iterate over all entries in the member map
     for( typename std::map< std::string, TYPE >::iterator i = member.begin() ;
        i!=member.end() ; ++i )
     {
     // the field name is the key to the map
     const std::string fieldName = i->first;

     // check to see if the field should be written
     std::map<std::string, FieldBase*>::const_iterator fieldAttributes =
        FieldInfo::AttributesByName.find(fieldName);

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
      // the field data is mapped value
      TYPE& fieldData = i->second;

      DBReadWrapper( fieldName, fieldData );
     }
     }
   */

  DBSetDir( m_dbFilePtr, ".." );

}
template void SiloFile::DBReadWrapper( const std::string&, std::map< std::string, array<integer> >& ) const;
template void SiloFile::DBReadWrapper( const std::string&, std::map< std::string, array<long int> >& ) const;
template void SiloFile::DBReadWrapper( const std::string&, std::map< std::string, array<long long int> >& ) const;
template void SiloFile::DBReadWrapper( const std::string&, std::map< std::string, lArray2d>& ) const;
template void SiloFile::DBReadWrapper( const std::string&, std::map< std::string, array<localIndex_array> >& ) const;
template void SiloFile::DBReadWrapper( const std::string&, std::map< std::string, lSet>& ) const;
template void SiloFile::DBReadWrapper( const std::string&, std::map< std::string, array<string> >& ) const;



template<typename TYPE>
void SiloFile::DBReadWrapper( const std::string& subdir, std::map< std::string, InterObjectRelation<TYPE> >& member ) const
{
  DBSetDir( m_dbFilePtr, subdir.c_str() );

  // iterate over all entries in the member map
  for( typename std::map< std::string, InterObjectRelation<TYPE> >::iterator i = member.begin() ; i!=member.end() ; ++i )
  {
    // the field name is the key to the map
    const std::string fieldName = i->first;

    // the field data is mapped value
    TYPE& fieldData = i->second;

    DBReadWrapper( fieldName, fieldData );
  }

  DBSetDir( m_dbFilePtr, ".." );

}
template void SiloFile::DBReadWrapper( const std::string&, std::map< std::string, OneToOneRelation>& ) const;
template void SiloFile::DBReadWrapper( const std::string&, std::map< std::string, FixedOneToManyRelation>& ) const;
template void SiloFile::DBReadWrapper( const std::string&, std::map< std::string, OrderedVariableOneToManyRelation>& ) const;
template void SiloFile::DBReadWrapper( const std::string&, std::map< std::string, UnorderedVariableOneToManyRelation>& ) const;


template<typename TYPE>
void SiloFile::DBReadWrapper( const std::string& name, TYPE* const data, const int size ) const
{
  if( DBInqVarExists( m_dbFilePtr, name.c_str() ) )
  {
    if( SiloFileUtilities::DB_TYPE<TYPE>() != DBGetVarType( m_dbFilePtr, name.c_str() ) )
    {
      GEOS_ERROR("SiloFile::DBReadWrapper: variable "+ name +" is incorrect type" );
    }
    if( static_cast<int>(sizeof(TYPE)*size) != DBGetVarByteLength( m_dbFilePtr, name.c_str() ) )
    {
      GEOS_ERROR("SiloFile::DBReadWrapper: variable "+ name +" is incorrect size" );
    }

    DBReadVar( m_dbFilePtr, name.c_str(), data );
  }
  else
  {
    //GEOS_ERROR("SiloFile::DBReadWrapper: variable "+ name +" does not exist in
    // silo file" );
  }
}
template void SiloFile::DBReadWrapper( const std::string&, int* const, const int ) const;



template <typename TYPE>
void** SiloFile::GetDataVar( const std::string& fieldName,
                             const std::string& meshName,
                             const typename array<TYPE>::size_type nels,
                             const int centering,
                             const int cycleNumber,
                             const realT problemTime,
                             const std::string& ) const
{

  const int nvars = SiloFileUtilities::GetNumberOfVariablesInField<TYPE>();

  void** rval = nullptr;

  const int meshType = GetMeshType( meshName );


  if( meshType == DB_UCDMESH )
  {
    const DBucdvar* const ucdVar = DBGetUcdvar( m_dbFilePtr, fieldName.data() );

    if( meshName.compare( ucdVar->meshname ) )
    {
      GEOS_ERROR("SiloFile::GetDataVar: meshname is not consistent for " + meshName + ":" + fieldName );
    }
    if( (int)nels != ucdVar->nels )
    {
      GEOS_ERROR("SiloFile::GetDataVar: size is not consistent for " + meshName + ":" + fieldName);
    }
    if( centering != ucdVar->centering )
    {
      GEOS_ERROR("SiloFile::GetDataVar: centering is not consistent for " + meshName + ":" + fieldName);
    }
    if( cycleNumber != ucdVar->cycle )
    {
      GEOS_ERROR("SiloFile::GetDataVar: cycleNumber is not consistent for " + meshName + ":" + fieldName);
    }
    if( !isEqual(problemTime,ucdVar->dtime) )
    {
      GEOS_ERROR("SiloFile::GetDataVar: problemTime is not consistent for " + meshName + ":" + fieldName);
    }
    if( nvars != ucdVar->nvals )
    {
      GEOS_ERROR("SiloFile::GetDataVar: size is not consistent for " + meshName + ":" + fieldName);
    }
    if( SiloFileUtilities::DB_TYPE<TYPE>() != ucdVar->datatype )
    {
      GEOS_ERROR("SiloFile::GetDataVar: datatype is not consistent for " + meshName + ":" + fieldName);
    }

    rval = ucdVar->vals;
  }
  else if( meshType == DB_POINTMESH )
  {
    const DBmeshvar* const meshVar = DBGetPointvar( m_dbFilePtr, fieldName.data() );
    rval = meshVar->vals;

  }
  else if( meshType == DB_QUADCURV )
  {
    const DBquadvar* const quadVar = DBGetQuadvar( m_dbFilePtr, fieldName.data() );
    rval = quadVar->vals;

  }
  else
  {
    GEOS_ERROR("SiloFile::GetDataVar: invalid meshtype");
  }


  return rval;
}
template void** SiloFile::GetDataVar<int>( const std::string&, const std::string&, const array<int>::size_type, const int, const int, const realT,
                                           const std::string& ) const;
template void** SiloFile::GetDataVar<long int>( const std::string&, const std::string&, const array<long int>::size_type, const int, const int, const realT,
                                                const std::string& ) const;
template void** SiloFile::GetDataVar<long long int>( const std::string&, const std::string&, const array<long long int>::size_type, const int, const int,
                                                     const realT, const std::string& ) const;

//template void** SiloFile::GetDataVar<int>( const std::string&, const
// std::string&, const array<int>::size_type , const int, const int, const
// realT, const std::string& ) const;
template void** SiloFile::GetDataVar<realT>( const std::string&, const std::string&, const array<realT>::size_type, const int, const int, const realT,
                                             const std::string& ) const;
template void** SiloFile::GetDataVar<R1Tensor>( const std::string&, const std::string&, const array<R1Tensor>::size_type, const int, const int, const realT,
                                                const std::string& ) const;
template void** SiloFile::GetDataVar<R2Tensor>( const std::string&, const std::string&, const array<R2Tensor>::size_type, const int, const int, const realT,
                                                const std::string& ) const;
template void** SiloFile::GetDataVar<R2SymTensor>( const std::string&, const std::string&, const array<R2SymTensor>::size_type, const int, const int,
                                                   const realT, const std::string& ) const;



integer_array SiloFile::SiloNodeOrdering()
{

  integer_array nodeOrdering;

//  if( !m_elementGeometryID.compare(0, 4, "CPE2") )
//  {
//    nodeOrdering.resize(2);
//    nodeOrdering[0] = 0;
//    nodeOrdering[1] = 1;
//  }
//  else if( !m_elementGeometryID.compare(0, 4, "CPE3") )
//  {
//    nodeOrdering.resize(3);
//    nodeOrdering[0] = 0;
//    nodeOrdering[1] = 1;
//    nodeOrdering[2] = 2;
//    //    throw GPException("ElementRegionT::AllocateElementLibrary(): CPE3
// unimplemented");
//  }
//  else if (!m_elementGeometryID.compare(0, 4, "CPE4"))
//  {
//    nodeOrdering.resize(4);
//    nodeOrdering[0] = 0;
//    nodeOrdering[1] = 1;
//    nodeOrdering[2] = 3;
//    nodeOrdering[3] = 2;
//  }
//  else if (!m_elementGeometryID.compare(0, 4, "C3D4"))
//  {
//    nodeOrdering.resize(4);
//    nodeOrdering[0] = 1;
//    nodeOrdering[1] = 0;
//    nodeOrdering[2] = 2;
//    nodeOrdering[3] = 3;
//  }
//  else if (!m_elementGeometryID.compare(0, 4, "C3D8") ||
// !m_elementGeometryID.compare(0, 4, "C3D6"))
//  {
  nodeOrdering.resize(8);
  nodeOrdering[0] = 0;
  nodeOrdering[1] = 1;
  nodeOrdering[2] = 3;
  nodeOrdering[3] = 2;
  nodeOrdering[4] = 4;
  nodeOrdering[5] = 5;
  nodeOrdering[6] = 7;
  nodeOrdering[7] = 6;
//  }
//  else if (!m_elementGeometryID.compare(0, 4, "STRI"))
//  {
//    nodeOrdering.resize(3);
//    nodeOrdering[0] = 0;
//    nodeOrdering[1] = 1;
//    nodeOrdering[2] = 2;
//  }
//  else if (!m_elementGeometryID.compare(0, 3, "S4R"))
//  {
//    nodeOrdering.resize(4);
//    nodeOrdering[0] = 0;
//    nodeOrdering[1] = 1;
//    nodeOrdering[2] = 2;
//    nodeOrdering[3] = 3;
//  }
//  else if (!m_elementGeometryID.compare(0, 4, "TRSH"))
//  {
//    nodeOrdering.resize(4);
//    nodeOrdering[0] = 0;
//    nodeOrdering[1] = 1;
//    nodeOrdering[2] = 2;
//  }
  return nodeOrdering;
}



void SiloFile::WriteManagedGroupSilo( ManagedGroup const * group,
                                      const std::string& siloDirName,
                                      const std::string& meshname,
                                      const int centering,
                                      const int cycleNum,
                                      const realT problemTime,
                                      const bool isRestart,
                                      string const & regionName,
                                      const localIndex_array& mask )
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

    MakeSubDirectory( shortsubdir, rootDirectory );
    DBSetDir(m_dbFilePtr, shortsubdir.c_str());
  }

  WriteViewWrappersToSilo<real64>( meshname,
                                   group->wrappers(),
                                   centering,
                                   cycleNum,
                                   problemTime,
                                   isRestart,
                                   rootDirectory,
                                   regionName,
                                   mask);


  DBSetDir(m_dbFilePtr, "..");

}



void SiloFile::WriteDomainPartition( DomainPartition const & domain,
                                     const int cycleNum,
                                     const realT problemTime,
                                     const bool isRestart )
{

  MeshLevel const * const mesh = domain.getMeshBody(0)->getMeshLevel(0);
  ConstitutiveManager const * const
  constitutiveManager = domain.GetGroup<ConstitutiveManager>(keys::ConstitutiveManager);

  WriteMeshLevel( mesh, constitutiveManager, cycleNum, problemTime, isRestart );


//  WriteCommonPlanes( siloFile, cycleNum, problemTime, isRestart, writeCP );

//  WriteCartesianGrid( siloFile, cycleNum, problemTime, isRestart, writeCG );
//  m_wellboreManager.WriteWellboreSilo( siloFile, cycleNum, problemTime,
// isRestart );

  if( isRestart )
  {
//    siloFile.DBWriteWrapper("m_globalDomainNumber",m_globalDomainNumber);
  }

}

void SiloFile::WriteMeshLevel( MeshLevel const * const meshLevel,
                               ConstitutiveManager const * const constitutiveManager,
                               const int cycleNum,
                               const realT problemTime,
                               const bool isRestart )
{
  int rank = 0;
#ifdef USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  //--------------WRITE FE DATA-----------------
//  if (m_feElementManager->m_numElems > 0)
  {

    NodeManager const * const nodeManager = meshLevel->getNodeManager();
    localIndex const numNodes = nodeManager->size();

    r1_array const & referencePosition = nodeManager->getReference<r1_array>(keys::referencePositionString);

    bool writeArbitraryPolygon(false);
    const std::string meshName("MeshLevel");

    //set the nodal coordinate data structure
    realT* coords[3];
    array<real64> xcoords(numNodes);
    array<real64> ycoords(numNodes);
    array<real64> zcoords(numNodes);
    for (localIndex a = 0 ; a < numNodes ; ++a)
    {
      R1Tensor nodePosition;
      nodePosition = referencePosition[a];

      xcoords[a] = nodePosition(0);
      ycoords[a] = nodePosition(1);
      zcoords[a] = nodePosition(2);
    }

    coords[0] = xcoords.data();
    coords[1] = ycoords.data();
    coords[2] = zcoords.data();

    ElementRegionManager const * const elementManager = meshLevel->getElemManager();
    const localIndex numElementRegions = elementManager->GetGroup(keys::elementRegions)->GetSubGroups().size();
    array<localIndex*> meshConnectivity(numElementRegions);
    array<int*> isGhostElement(numElementRegions);
    array<globalIndex*> globalElementNumbers(numElementRegions);
    array<integer> shapecnt(numElementRegions);
    array<integer> shapetype(numElementRegions);
    array<integer> shapesize(numElementRegions);

    array<FixedOneToManyRelation> elementToNodeMap;
    elementToNodeMap.resize( numElementRegions );

    int count = 0;
    ManagedGroup const * elementRegions = elementManager->GetGroup(dataRepository::keys::elementRegions);

    for( auto const & region : elementRegions->GetSubGroups() )
    {
      ManagedGroup const * cellBlockSubRegions = region.second->GetGroup(dataRepository::keys::cellBlockSubRegions);
      for( auto const & iterCellBlocks : cellBlockSubRegions->GetSubGroups() )
      {
        CellBlockSubRegion const * cellBlock = cellBlockSubRegions->GetGroup<CellBlockSubRegion>(iterCellBlocks.first);

        lArray2d const & elemsToNodes = cellBlock->getWrapper<FixedOneToManyRelation>(cellBlock->viewKeys().nodeList)->reference();// getData<lArray2d>(keys::nodeList);

        // The following line seems to be redundant. It's actual function is to
        // size this temp array.(pfu)
        elementToNodeMap[count].resize(elemsToNodes.size(0),elemsToNodes.size(1));

        for (localIndex k = 0 ; k < cellBlock->size() ; ++k)
        {
          arrayView1d<localIndex const> const elemToNodeMap = elemsToNodes[k];

          const integer_array nodeOrdering = SiloNodeOrdering();
          integer numNodesPerElement = integer_conversion<int>(elemsToNodes.size(1));
          for (localIndex a = 0 ; a < numNodesPerElement ; ++a)
          {
            elementToNodeMap[count](k, a) = elemToNodeMap[nodeOrdering[a]];
          }
        }


        meshConnectivity[count] = elementToNodeMap[count].data();

//      isGhostElement[count] =
// (cellBlock->GetFieldData<FieldInfo::ghostRank>()).data();

//      globalElementNumbers[count] = elementRegion.m_localToGlobalMap.data();
        shapecnt[count] = static_cast<int>(cellBlock->size());

//      if ( !elementRegion.m_elementGeometryID.compare(0, 4, "C3D8") )
//      {
        shapetype[count] = DB_ZONETYPE_HEX;
//      }
//      else if ( !elementRegion.m_elementGeometryID.compare(0, 4, "C3D6") )
//      {
//        shapetype[count] = DB_ZONETYPE_HEX;
//        writeArbitraryPolygon = true;
//      }
//      else if ( !elementRegion.m_elementGeometryID.compare(0, 4, "C3D4") )
//      {
//        shapetype[count] = DB_ZONETYPE_TET;
//      }
//      else if ( !elementRegion.m_elementGeometryID.compare(0, 4, "CPE4") ||
// !elementRegion.m_elementGeometryID.compare(0, 3, "S4R") )
//      {
//        shapetype[count] = DB_ZONETYPE_QUAD;
//      }
//      else if ( !elementRegion.m_elementGeometryID.compare(0, 4, "STRI") ||
// !elementRegion.m_elementGeometryID.compare(0, 4, "TRSH") ||
// !elementRegion.m_elementGeometryID.compare(0, 4, "CPE3"))
//      {
//        shapetype[count] = DB_ZONETYPE_TRIANGLE;
//      }
//      else if ( !elementRegion.m_elementGeometryID.compare(0, 4, "CPE2") )
//      {
//        shapetype[count] = DB_ZONETYPE_TRIANGLE;
//      }
//      else
//      {
//        GEOS_ERROR("PhysicalDomainT::WriteFiniteElementMesh: Do not recognize
// geometry type " + elementRegion.m_elementGeometryID + " \n");
//      }

        shapesize[count] = integer_conversion<int>(elemsToNodes.size(1));
        count++;
      }
    }

    WriteMeshObject(meshName,
                    numNodes,
                    coords,
                    nullptr,
                    integer_conversion<int>(numElementRegions),
                    shapecnt.data(),
                    meshConnectivity.data(),
                    nullptr /*globalElementNumbers.data()*/,
                    isGhostElement.data(),
                    shapetype.data(),
                    shapesize.data(),
                    cycleNum, problemTime);


    // write node fields in silo mesh, and all restart data as unassociated
    // variables.



    WriteManagedGroupSilo( nodeManager,
                           "NodalFields",
                           meshName,
                           DB_NODECENT,
                           cycleNum,
                           problemTime,
                           isRestart,
                           "none",
                           localIndex_array());








    {
    ElementRegionManager const * const elemManager = meshLevel->getElemManager();


    this->WriteRegionSpecifications( elemManager,
                                     constitutiveManager,
                                     meshName,
                                     cycleNum,
                                     problemTime );

    for( localIndex er=0 ; er<elemManager->numRegions() ; ++er )
    {
      ElementRegion const * const elemRegion = elemManager->GetRegion(er);
      for( localIndex esr=0 ; esr<elemRegion->numSubRegions() ; ++esr )
      {
        CellBlockSubRegion const * const subRegion = elemRegion->GetSubRegion(esr);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        string regionName = elemRegion->getName() + "_" + subRegion->getName();

        WriteManagedGroupSilo( subRegion,
                               regionName,
                               meshName,
                               DB_ZONECENT,
                               cycleNum,
                               problemTime,
                               isRestart,
                               "none",
                               localIndex_array());
      }
    }
    for( localIndex matIndex=0 ; matIndex<constitutiveManager->GetSubGroups().size() ; ++matIndex )
    {
      ConstitutiveBase const * const
      constitutiveModel = constitutiveManager->GetGroup<ConstitutiveBase>(matIndex);

      WriteManagedGroupSilo( constitutiveModel->GetStateData(),
                             constitutiveModel->getName(),
                             meshName,
                             DB_ZONECENT,
                             cycleNum,
                             problemTime,
                             isRestart,
                             constitutiveModel->getName(),
                             localIndex_array());
    }

    }
//    m_feElementManager->WriteSilo( siloFile, meshName, cycleNum, problemTime,
// isRestart );


  }//end FE write
}

}
#pragma GCC diagnostic pop
