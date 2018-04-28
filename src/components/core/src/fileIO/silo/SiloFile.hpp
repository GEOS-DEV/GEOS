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
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
//  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL
// SECURITY,
//  LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
//  AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
// TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
//  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
//  1. This notice is required to be provided under our contract with the U.S.
// Department of Energy (DOE). This work was produced at Lawrence Livermore
//     National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
//  2. Neither the United States Government nor Lawrence Livermore National
// Security, LLC nor any of their employees, makes any warranty, express or
//     implied, or assumes any liability or responsibility for the accuracy,
// completeness, or usefulness of any information, apparatus, product, or
//     process disclosed, or represents that its use would not infringe
// privately-owned rights.
//  3. Also, reference herein to any specific commercial products, process, or
// services by trade name, trademark, manufacturer or otherwise does not
//     necessarily constitute or imply its endorsement, recommendation, or
// favoring by the United States Government or Lawrence Livermore National
// Security,
//     LLC. The views and opinions of authors expressed herein do not
// necessarily state or reflect those of the United States Government or
// Lawrence
//     Livermore National Security, LLC, and shall not be used for advertising
// or product endorsement purposes.
//
//  This Software derives from a BSD open source release LLNL-CODE-656616. The
// BSD  License statment is included in this distribution in src/bsd_notice.txt.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @file SiloFile.h
 * @author settgast1
 * @date Jan 24, 2011
 */

#ifndef SILOFILE_H_
#define SILOFILE_H_

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


#include "common/InterObjectRelation.hpp"



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

  /// Initializes the silo library
  void Initialize( const PMPIO_iomode_t readwrite );

  /// finishes up the silo library usage
  void Finish();

  /// Wait for the Baton when doing PMPIO
  void WaitForBaton( const int domainNumber, const int cycleNum, const bool isRestart );

  void WaitForBaton( const int domainNumber, const std::string& restartFileName );

  /// Hand off the Baton when doing PMPIO
  void HandOffBaton();


  void MakeSubDirectory( const std::string& subdir, const std::string& rootdir )
  {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0)
    {
      DBMkDir(m_dbFilePtr, rootdir.c_str());
    }
    DBMkDir(m_dbFilePtr, subdir.c_str());
  }

  /// Write out a single mesh object
  void WriteMeshObject(const std::string& meshName,
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
                       const realT problemTime);

  /// Write out a polygon mesh object
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
   * @brief Write out a single discrete element mesh object
   * @author Scott Johnson
   */
  void WriteDiscreteElementMeshObject(const std::string& meshName,
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
                                      const realT problemTime);

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
                const int cycleNumber,
                const realT problemTime);

  /**
   * @brief Write out a single ellipsoidal discrete element mesh object
   * @author Scott Johnson
   */
  void WriteDiscreteElementCSGObject(const std::string& meshName,
                                     const array<R1Tensor>& x,
                                     const array<R1Tensor>& r,
                                     const array<R2Tensor>& R,
                                     const int cycleNumber,
                                     const realT problemTime);

  void WritePointMesh( const std::string& meshName,
                       const localIndex numPoints,
                       realT* coords[3],
                       const int cycleNumber,
                       const realT problemTime );

  void WriteBeamMesh(const std::string& meshName,
                     const localIndex nnodes,
                     realT* coords[3],
                     const localIndex_array& node1,
                     const localIndex_array& node2,
                     const int cycleNumber,
                     const realT problemTime);

  void WriteBeamMesh(const std::string& meshName,
                     const localIndex nnodes,
                     realT* coords[3],
                     const std::map<int, int>& connectivity,
                     const int cycleNumber,
                     const realT problemTime);

  void WriteBeamMesh(const std::string& meshName,
                     const localIndex nnodes,
                     realT* coords[3],
                     integer_array& nodelist,
                     const int cycleNumber,
                     const realT problemTime);

  /**
   * @brief Write out a single "quad" (i.e. hex) mesh object, and record the
   * mesh dimensions for subsequent write data field calls.
   * @author walsh24
   */
  int  WriteQuadMeshObject(const std::string& meshName,
                           const localIndex nX,
                           const localIndex nY,
                           const localIndex nZ,
                           realT* coords[3],
                           const int cycleNumber,
                           const realT problemTime);

//  void TestWriteDiscreteElementMeshObject();

  void WriteRegionSpecifications(const ElementRegionManager*  elementManager,
                                 constitutive::ConstitutiveManager const * constitutiveManager,
                                 const std::string& meshName,
                                 const int cycleNumber,
                                 const realT problemTime);


  void WriteManagedGroupSilo( dataRepository::ManagedGroup const * group,
                              const std::string& siloDirName,
                              const std::string& meshname,
                              const int centering,
                              const int cycleNum,
                              const realT problemTime,
                              const bool isRestart,
                              const localIndex_array& mask );



  void WriteManagedGroupSilo( dataRepository::ManagedGroup const * group,
                              const std::string& meshname,
                              const int centering,
                              const int cycleNum,
                              const realT problemTime,
                              const bool isRestart,
                              const std::string& multiRoot,
                              const localIndex_array& mask );

  /// writes out fields in a data member map
  template< typename OUTPUTTYPE >
  void WriteViewWrappersToSilo( const std::string& meshname,
                                const dataRepository::ManagedGroup::viewWrapperMap & wrappers,
                                const int centering,
                                const int cycleNum,
                                const realT problemTime,
                                const bool isRestart,
                                const std::string& multiRoot,
                                const std::string& regionName,
                                const localIndex_array& mask );

  template< typename INPUTTYPE, typename TYPE >
  void ReadFieldMapFromSilo( std::map< std::string, array<TYPE> >& member,
                             const std::string& meshname,
                             const int centering,
                             const int cycleNum,
                             const realT problemTime,
                             const bool isRestart,
                             const std::string& regionName,
                             const localIndex_array& mask ) const;

  /// Write out a data field
  template<typename OUTTYPE, typename TYPE>
  void WriteDataField( const std::string& meshName,
                       const std::string& fieldName,
                       const TYPE& field,
                       const int centering,
                       const int cycleNumber,
                       const realT problemTime,
                       const std::string& multiRoot,
                       const std::string& regionName );

  template<typename OUTTYPE, typename TYPE>
  void WriteDataField( const std::string& meshName,
                       const std::string& fieldName,
                       const array<TYPE>& field,
                       const int centering,
                       const int cycleNumber,
                       const realT problemTime,
                       const std::string& multiRoot,
                       const std::string& regionName );

  /// Write out a data field
  template< typename INPUTTYPE, typename TYPE>
  void ReadDataField( array<TYPE>& field,
                      const std::string& meshName,
                      const std::string& fieldName,
                      const int centering,
                      const int cycleNumber,
                      const realT problemTime,
                      const std::string& regionName ) const;

  int GetMeshType( const std::string& meshName ) const
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
  void** GetDataVar( const std::string& fieldName,
                     const std::string& meshName,
                     const typename array<TYPE>::size_type nels,
                     const int centering,
                     const int cycleNumber,
                     const realT problemTime,
                     const std::string& ) const;

  /// Write out a multi-mesh object
  template< typename CBF >
  void WriteMultiXXXX(const DBObjectType type, CBF DBPutMultiCB,
                      const int centering, const std::string name, const int cycleNumber,
                      const std::string& multiRoot, const DBoptlist* optlist = NULL);


  template<typename TYPE>
  void DBWriteWrapper( const std::string& name, const TYPE& data );

  void DBWriteWrapper( const std::string& name, const string_array& data );

  template<typename TYPE>
  void DBWriteWrapper( const std::string& name, const array<TYPE>& data );


  template<typename TYPE>
  void DBWriteWrapper( const std::string& name, const set<TYPE>& data );

  template<typename TYPE>
  void DBWriteWrapper( const std::string& name, const Array2dT<TYPE>& data );

  template<typename TYPE>
  void DBWriteWrapper( const std::string& name, const array<array<TYPE> >& data );

  template<typename TYPE>
  void DBWriteWrapper( const std::string& name, const array<Array2dT<TYPE> >& data );

  template<typename TYPE>
  void DBWriteWrapper( const std::string& name, const array<array<array<TYPE> > >& data );

  template<typename TYPE>
  void DBWriteWrapper( const std::string& name, const array<set<TYPE> >& data );

  template< typename T1, typename T2 >
  void DBWriteWrapper( const std::string& name, const std::map< T1, T2 >& datamap );

  template<typename TYPE>
  void DBWriteWrapper( const std::string& subdir, const std::map< std::string, TYPE>& member );

  template<typename TYPE>
  void DBWriteWrapper( const std::string& subdir, const std::map< std::string, InterObjectRelation<TYPE> >& member );

  template<typename TYPE>
  void DBWriteWrapper( const std::string& name, const TYPE* const data, const int size );

  template<typename TYPE>
  void DBReadWrapper( const std::string& name, TYPE& data ) const;

  template<typename TYPE>
  void DBReadWrapper( const std::string& name, array<TYPE>& data ) const;

  template<typename TYPE>
  void DBReadWrapper( const std::string& name, set<TYPE>& data ) const;

  template<typename TYPE>
  void DBReadWrapper( const std::string& name, Array2dT<TYPE>& data ) const;

  template<typename TYPE>
  void DBReadWrapper( const std::string& name, array<array<TYPE> >& data ) const;

  template<typename TYPE>
  void DBReadWrapper( const std::string& name, array<Array2dT<TYPE> >& data ) const;

  template<typename TYPE>
  void DBReadWrapper( const std::string& name, array<array<array<TYPE> > >& data ) const;

  template<typename TYPE>
  void DBReadWrapper( const std::string& name, array<set<TYPE> >& data ) const;

  template< typename T1, typename T2 >
  void DBReadWrapper( const std::string& name, std::map< T1, T2 >& datamap ) const;

  template<typename TYPE>
  void DBReadWrapper( const std::string& subdir, std::map< std::string, TYPE>& member ) const;

  template<typename TYPE>
  void DBReadWrapper( const std::string& subdir, std::map< std::string, InterObjectRelation<TYPE> >& member ) const;

  template<typename TYPE>
  void DBReadWrapper( const std::string& name, TYPE* const data, const int size ) const;

  /// dummy function to get rid of compiler warnings about unused functions from PMPIO's use of static function
  /// definitions in the header. This should go away as we fully utilize PMPIO.
  void StopSiloCompilerWarnings();

  void ClearEmptiesFromMultiObjects(const int cycleNum);

  //private:

  /// pointer to the DBfile that this class is working on
  DBfile* m_dbFilePtr;

  /// total number of "parallel" files to write out
  int m_numGroups;

  /// the pmpio baton. A processor needs this to write to the file.
  PMPIO_baton_t *m_baton;

  /// which database to use. DB_PDB or DB_HDF5
  int m_driver;

  /// root of the filename that we will be reading/writing
  std::string m_fileRoot;

  std::string m_restartFileRoot;

  std::string m_slaveDirectory = "plotFiles";

  std::string m_fileName;

  std::string m_baseFileName;

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
 * @tparam OUTTYPE the type of data to write out (int,float,realT)
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
 * float or realT.
 * 3 for a R1Tensor...etc.
 */
template<typename TYPE>
int GetNumberOfVariablesInField();

template<typename TYPE>
int GetTensorRank();

template<typename OUTTYPE, typename TYPE>
OUTTYPE CastField(const TYPE& field, const int i = 0);     // avoids compiler
                                                           // warning

template<typename OUTTYPE, typename TYPE>
OUTTYPE CastField(const TYPE& field, const int i)
{
  return static_cast<OUTTYPE>(field.Data()[i]);
//    return field;
}

template<> inline int CastField<int, int> (const int& field, const int )
{
  return field;
}


//  template<> inline int CastField<int, localIndex> (const localIndex& field,
// const int )
//  {
//    return static_cast<int>(field);
//  }

//  template<> inline localIndex CastField<localIndex, localIndex> (const
// localIndex& field, const int )
//  {
//    return field;
//  }


template<> inline globalIndex CastField<globalIndex, globalIndex> (const globalIndex& field, const int )
{
  return field;
}

template<> inline int CastField<int, long long unsigned int > (const long long unsigned int& field, const int )
{
  return static_cast<int>(field);
}


template<> inline realT CastField<realT, realT> (const realT& field, const int )
{
  return field;
}
template<> inline float CastField<float, realT> (const realT& field, const int )
{
  return static_cast<float> (field);
}


template< typename TYPE, typename PTR_TYPE> inline PTR_TYPE* DataPtr( TYPE& data)
{
  return static_cast<PTR_TYPE*>(&data);
}

template<> inline realT* DataPtr( R1Tensor& data)        {    return data.begin();  }
template<> inline realT* DataPtr( R2Tensor& data)        {    return data.begin();  }
template<> inline realT* DataPtr( R2SymTensor& data)  {    return data.begin();  }



template<typename TYPE>
void SetVariableNames(const std::string& fieldName, array<string>& varnamestring, char* varnames[]);

template<typename OBJECT_TYPE>
int FieldCentering();

void SetCenteringSubdir(const int centering, std::string& subdir);

}


template< typename OUTPUTTYPE >
void SiloFile::WriteViewWrappersToSilo( const std::string& meshname,
                                        const dataRepository::ManagedGroup::viewWrapperMap & wrappers,
                                        const int centering,
                                        const int cycleNum,
                                        const realT problemTime,
                                        const bool isRestart,
                                        const std::string& multiRoot,
                                        const std::string& regionName,
                                        const localIndex_array& mask )
{

  // iterate over all entries in the member map
  for( auto const & wrapperIter : wrappers )
  {
    auto const & wrapper = wrapperIter.second;
    // the field name is the key to the map
    const std::string fieldName = wrapper->getName();

    std::type_info const & typeID = wrapper->get_typeid();

    // TODO This is wrong. problem with uniqueness
    if( typeID==typeid(real64_array) )
    {
      auto const & viewWrapperT = dynamic_cast< dataRepository::ViewWrapper<real64_array> const & >( *wrapper );
      this->WriteDataField<real64>(meshname.c_str(), fieldName, viewWrapperT.reference(), centering, cycleNum, problemTime, multiRoot, regionName );
    }
    if( typeID==typeid(r1_array) )
    {
      auto const & viewWrapperT = dynamic_cast< dataRepository::ViewWrapper<r1_array> const & >( *wrapper );
      this->WriteDataField<real64>(meshname.c_str(), fieldName, viewWrapperT.reference(), centering, cycleNum, problemTime, multiRoot, regionName );
    }
    if( typeID==typeid(integer_array) )
    {
      auto const & viewWrapperT = dynamic_cast< dataRepository::ViewWrapper<integer_array> const & >( *wrapper );
      this->WriteDataField<integer>(meshname.c_str(), fieldName, viewWrapperT.reference(), centering, cycleNum, problemTime, multiRoot, regionName );
    }
    if( typeID==typeid(localIndex_array) )
    {
      auto const & viewWrapperT = dynamic_cast< dataRepository::ViewWrapper<localIndex_array> const & >( *wrapper );
      this->WriteDataField<localIndex>(meshname.c_str(), fieldName, viewWrapperT.reference(), centering, cycleNum, problemTime, multiRoot, regionName );
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
void SiloFile::WriteDataField( const std::string& meshName,
                               const std::string& fieldName,
                               const TYPE& field,
                               const int centering,
                               const int cycleNumber,
                               const realT problemTime,
                               const std::string& multiRoot,
                               const std::string& regionName )
{}

template<typename OUTTYPE, typename TYPE>
void SiloFile::WriteDataField( const std::string& meshName,
                               const std::string& fieldName,
                               const array<TYPE>& field,
                               const int centering,
                               const int cycleNumber,
                               const realT problemTime,
                               const std::string& multiRoot,
                               const std::string& regionName )
{
  const int nvars = SiloFileUtilities::GetNumberOfVariablesInField<TYPE>();
  const int nels = field.size();

  const int meshType = GetMeshType( meshName );


  DBoptlist *optlist = DBMakeOptlist(5);
  DBAddOption(optlist, DBOPT_CYCLE, const_cast<int*> (&cycleNumber));
  DBAddOption(optlist, DBOPT_DTIME, const_cast<realT*> (&problemTime));

  char *regionpnames[2] =
  { NULL, NULL };



  if (regionName != "none")
  {
    regionpnames[0] = const_cast<char*> (regionName.c_str());
    regionpnames[1] = NULL;
    DBAddOption(optlist, DBOPT_REGION_PNAMES, &regionpnames);
  }

  // if the number of elements is zero, then record the path to the var. This
  // will be used later to delete the entry
  // from the multivar.
  if (nels == 0)
  {
    char pwd[256];
    DBGetDir(m_dbFilePtr, pwd);
    std::string emptyObject = pwd;
    emptyObject += "/" + fieldName;
    m_emptyVariables.push_back(emptyObject);
  }
  else
  {
    array<char*> varnames(nvars);
    array<void*> vars(nvars);


    array<string> varnamestring(nvars);
    std::vector<std::vector<OUTTYPE> > castedField(nvars);

    SiloFileUtilities::SetVariableNames<TYPE>(fieldName, varnamestring, varnames.data() );

    for (int i = 0 ; i < nvars ; ++i)
    {
      castedField[i].resize(nels);
      vars[i] = static_cast<void*> (&(castedField[i][0]));
      for (int k = 0 ; k < nels ; ++k)
      {
        castedField[i][k] = SiloFileUtilities::CastField<OUTTYPE>(field[k], i);
      }
    }

    int err = -2;
    if( meshType == DB_UCDMESH )
    {
      err = DBPutUcdvar( m_dbFilePtr, fieldName.c_str(), meshName.c_str(), nvars, varnames.data(), reinterpret_cast<float**>(vars.data()),
                         nels, NULL, 0, SiloFileUtilities::DB_TYPE<OUTTYPE>(), centering, optlist);
    }
    else if( meshType == DB_POINTMESH )
    {
      err = DBPutPointvar( m_dbFilePtr, fieldName.c_str(), meshName.c_str(), nvars, reinterpret_cast<float**>(vars.data()),
                           nels, SiloFileUtilities::DB_TYPE<OUTTYPE>(), optlist);
    }
    else if( meshType == DB_QUADCURV )
    {
      err = DBPutQuadvar( m_dbFilePtr, fieldName.c_str(), meshName.c_str(), nvars, varnames.data(), reinterpret_cast<float**>(vars.data()),
                          m_quadMeshDims,m_quadMeshNDims,NULL, 0,  SiloFileUtilities::DB_TYPE<OUTTYPE>(),centering, optlist);
    }
    if(err < 0)
    {
      if(err < -1)
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
  if (rank == 0)
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
void SiloFile::ReadFieldMapFromSilo( std::map< std::string, array<TYPE> >& member,
                                     const std::string& meshname,
                                     const int centering,
                                     const int cycleNum,
                                     const realT problemTime,
                                     const bool isRestart,
                                     const std::string& regionName,
                                     const localIndex_array& mask ) const
{
  // iterate over all entries in the member map
  for( typename std::map< std::string, array<TYPE> >::iterator iter = member.begin() ; iter!=member.end() ; ++iter )
  {
    // the field name is the key to the map
    const std::string fieldName = iter->first;

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

          for( localIndex_array::size_type i = 0 ; i < mask.size() ; ++i)
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
                              const std::string& meshName,
                              const std::string& fieldName,
                              const int centering,
                              const int cycleNumber,
                              const realT problemTime,
                              const std::string& regionName ) const
{

  INPUTTYPE** var = static_cast<INPUTTYPE**>( GetDataVar<TYPE>( fieldName, meshName, field.size(), centering, cycleNumber, problemTime, regionName ) );

  for( typename array<TYPE>::size_type a=0 ; a<field.size() ; ++a )
  {
    TYPE temp;
    INPUTTYPE* ptemp = SiloFileUtilities::DataPtr<TYPE,INPUTTYPE>( temp );

    for( int i=0 ; i<SiloFileUtilities::GetNumberOfVariablesInField<TYPE>() ; ++i  )
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
                               const int centering,
                               const std::string name,
                               const int,
                               const std::string& multiRoot,
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

  DBGetDir(m_dbFilePtr, currentDirectory);
  DBSetDir(m_dbFilePtr, multiRoot.c_str());


  std::string multiRootString(multiRoot);
  if( !(multiRootString.compare("/")) )
  {
    multiRootString.clear();
  }

  for (int i = 0 ; i < size ; ++i)
  {
    int groupRank = PMPIO_GroupRank(m_baton, i);

    if (groupRank == 0)
    {

      sprintf(tempBuffer, "/domain_%04d%s/%s", i, multiRootString.c_str(), name.c_str());
      //      sprintf(tempBuffer, "%s_%04d:/domain_%04d%s/%s",
      // m_fileRoot.c_str(),
      //              cycleNumber, i, multiRoot.c_str(), name.c_str() );

    }
    else
    {
      if (m_slaveDirectory.empty())
      {
        sprintf(tempBuffer, "%s.%03d:/domain_%04d%s/%s", m_baseFileName.c_str(),
                groupRank, i, multiRootString.c_str(), name.c_str());
      }
      else
      {
        sprintf(tempBuffer, "%s%s%s.%03d:/domain_%04d%s/%s", m_slaveDirectory.c_str(), "/", m_baseFileName.c_str(),
                groupRank, i, multiRootString.c_str(), name.c_str());

      }

    }
    vBlockNames[i] = tempBuffer;
    BlockNames[i] = const_cast<char*>( vBlockNames[i].c_str() );
    blockTypes[i] = type;
  }

  std::string multiName = name;
  DBPutMultiCB(m_dbFilePtr, multiName.c_str(), size, BlockNames.data(), blockTypes.data(),
               const_cast<DBoptlist*> (optlist));

  DBSetDir(m_dbFilePtr, currentDirectory);

}


template<typename TYPE>
void SiloFile::DBWriteWrapper( const std::string& name, const array<TYPE>& data )
{
  if( !data.empty() )
  {
    int dims[2];
    dims[0] = SiloFileUtilities::GetNumberOfVariablesInField<TYPE>();
    dims[1] = data.size();

    DBWrite( m_dbFilePtr, name.c_str(), const_cast<TYPE*>(data.data()), dims,
             2, SiloFileUtilities::DB_TYPE<TYPE>() );
  }

}


}
#endif /* SILOFILE_H_ */
