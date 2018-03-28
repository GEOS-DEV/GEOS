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
 * @file TableManager.h
 * @author settgast1
 * @date Mar 4, 2011
 */

#ifndef TABLEMANAGER_H_
#define TABLEMANAGER_H_

#include "Common/typedefs.h"
#include "DataStructures/Tables/Table.h"
#include "Utilities/StringUtilities.h"
#include "Utilities/IOUtilities.h"
#include "Utilities/Utilities.h"
#include "ObjectManagers/FunctionManager.h"
#include "DataStructures/Tables/TableTypes.h"

#include <map>

#include "../../codingUtilities/Functions.hpp"
#include "../IO/ticpp/HierarchicalDataNode.h.old"

class TableManager
{
public:

  static TableManager& Instance()
  {
    static TableManager theTableManager;
    return theTableManager;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Tables

  //CONST
  template < unsigned int dim>
  inline const std::map<std::string,  Table<dim, realT> >& Tables() const { throw GPException("Cannot call base specialization"); }

  //NON-CONST
  template < unsigned int dim>
  inline std::map<std::string,  Table<dim, realT> >& Tables() { throw GPException("Cannot call base specialization"); }

  //ADD
  template < unsigned int dim, class ARRAY, class ARRAY2>
  inline void NewTable(const std::string& name, const ARRAY2& x, const ARRAY& values, TableInterpolation::Order interp) {
    throw GPException("Cannot call base specialization");
  }

  template <unsigned int dim, class ARRAY>
  inline Table<dim,realT> * GetTable( const std::string& tableName )
  {
    typename std::map<std::string, Table<dim, realT> >::iterator table = Tables<dim>().find(tableName);
    if (table == Tables<dim>().end())
      throw GPException("Table name " + tableName + " not found.\n");
    return &(table->second);

  }

  //LOOKUP
  template <unsigned int dim, class ARRAY>
  inline realT LookupTable(const std::string& tableName,
                           const ARRAY& key,
                           TableInterpolation::Order interpolate = TableInterpolation::linear) const
  {
    typename std::map<std::string, Table<dim, realT> >::const_iterator table =
      Tables<dim>().find(tableName);
    if (table == Tables<dim>().end())
      throw GPException("Table name " + tableName + " not found.\n");
    return table->second.Lookup(key, interpolate);
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // VectorFields

  //CONST
  template < unsigned int dim>
  inline const std::map<std::string,  Table<dim, R1Tensor> >& VectorFields() const { throw GPException("Cannot call base specialization"); }

  //NON-CONST
  template < unsigned int dim>
  inline std::map<std::string,  Table<dim, R1Tensor> >& VectorFields() { throw GPException("Cannot call base specialization"); }

  //ADD
  template < unsigned int dim, class ARRAY, class ARRAY2>
  inline void NewVectorField(const std::string& name, const ARRAY2& x, const ARRAY& values) { throw GPException("Cannot call base specialization"); }

  //LOOKUP
  template <unsigned int dim, class ARRAY>
  inline R1Tensor LookupVectorField(const std::string& tableName,
                                    const ARRAY& key,
                                    const bool interpolate = true) const
  {
    typename std::map<std::string, Table<dim, R1Tensor> >::const_iterator table =
      VectorFields<dim>().find(tableName);
    if (table == Tables<dim>().end())
      throw GPException("VectorField name " + tableName + " not found.\n");
    return table->second.Lookup(key, interpolate);
  }


  void ReadXML(xmlWrapper::xmlNode* hdn);

private:
  std::map<std::string, Table1D > m_tables1;
  std::map<std::string, Table2D > m_tables2;
  std::map<std::string, Table3D > m_tables3;
  std::map<std::string, Table4D > m_tables4;

  std::map<std::string, VectorField1D > m_vectorFields1;
  std::map<std::string, VectorField2D > m_vectorFields2;
  std::map<std::string, VectorField3D > m_vectorFields3;
  std::map<std::string, VectorField4D > m_vectorFields4;

  TableManager():
    m_tables1(),
    m_tables2(),
    m_tables3(),
    m_tables4(),
    m_vectorFields1(),
    m_vectorFields2(),
    m_vectorFields3(),
    m_vectorFields4()
  {}

  virtual ~TableManager() {}

  TableManager(const TableManager&);
  TableManager& operator=(const TableManager&);

  template < class ARRAY >
  void ReadVoxelFile(const std::string& filename, localIndex nComponents, ARRAY& values);

  template < class ARRAY >
  void ReadTimeVoxelFile(const std::string& filename, localIndex nComponents, ARRAY& values);

  template < class ARRAY, class ARRAY2 >
  void ReadNUFTFile(const std::string& filename,
                    ARRAY2& x,
                    ARRAY& values);

};



template <>
inline const std::map<std::string,  Table1D >& TableManager::Tables<1>() const { return m_tables1; }

template <>
inline const std::map<std::string,  Table2D >& TableManager::Tables<2>() const { return m_tables2; }

template <>
inline const std::map<std::string,  Table3D >& TableManager::Tables<3>() const { return m_tables3; }

template <>
inline const std::map<std::string,  Table4D >& TableManager::Tables<4>() const { return m_tables4; }



template <>
inline std::map<std::string,  Table1D >& TableManager::Tables<1>() { return m_tables1; }

template <>
inline std::map<std::string,  Table2D >& TableManager::Tables<2>() { return m_tables2; }

template <>
inline std::map<std::string,  Table3D >& TableManager::Tables<3>() { return m_tables3; }

template <>
inline std::map<std::string,  Table4D >& TableManager::Tables<4>() { return m_tables4; }



template <>
inline void TableManager::NewTable<1>(const std::string& name, const array<array<real64> >& x, const array<real64>& values, TableInterpolation::Order interp)
{
  m_tables1.insert(std::make_pair(name, Table1D()));
  Table1D& table = m_tables1[name];
  table.SetGrid(x);
  table.SetValues(values);
  table.SetInterpolation(interp);
  Function* tableFuncPtr = new ::Lookup1DTable(name, &table);
  FunctionManager::Instance().AddFunction(name, tableFuncPtr);
}

template <>
inline void TableManager::NewTable<2>(const std::string& name, const array<array<real64> >& x, const array<real64>& values,TableInterpolation::Order interp)
{
  m_tables2.insert(std::make_pair(name, Table2D()));
  Table2D& table = m_tables2[name];
  table.SetGrid(x);
  table.SetValues(values);
  table.SetInterpolation(interp);
  Function* tableFuncPtr = new ::Lookup2DTable(name, &table);
  FunctionManager::Instance().AddFunction(name, tableFuncPtr);
}

template <>
inline void TableManager::NewTable<3>(const std::string& name, const array<array<real64> >& x, const array<real64>& values,TableInterpolation::Order interp)
{
  m_tables3.insert(std::make_pair(name, Table3D()));
  Table3D& table = m_tables3[name];
  table.SetGrid(x);
  table.SetValues(values);
  table.SetInterpolation(interp);
  Function* tableFuncPtr = new ::Lookup3DTable(name, &table);
  FunctionManager::Instance().AddFunction(name, tableFuncPtr);
}

template <>
inline void TableManager::NewTable<4>(const std::string& name, const array<array<real64> >& x, const array<real64>& values,TableInterpolation::Order interp)
{
  m_tables4.insert(std::make_pair(name, Table4D()));
  Table4D& table = m_tables4[name];
  table.SetGrid(x);
  table.SetValues(values);
  table.SetInterpolation(interp);
  Function* tableFuncPtr = new ::Lookup4DTable(name, &table);
  FunctionManager::Instance().AddFunction(name, tableFuncPtr);
}



template <>
inline const std::map<std::string,  VectorField1D >& TableManager::VectorFields<1>() const { return m_vectorFields1; }

template <>
inline const std::map<std::string,  VectorField2D >& TableManager::VectorFields<2>() const { return m_vectorFields2; }

template <>
inline const std::map<std::string,  VectorField3D >& TableManager::VectorFields<3>() const { return m_vectorFields3; }

template <>
inline const std::map<std::string,  VectorField4D >& TableManager::VectorFields<4>() const { return m_vectorFields4; }



template <>
inline std::map<std::string,  VectorField1D >& TableManager::VectorFields<1>() { return m_vectorFields1; }

template <>
inline std::map<std::string,  VectorField2D >& TableManager::VectorFields<2>() { return m_vectorFields2; }

template <>
inline std::map<std::string,  VectorField3D >& TableManager::VectorFields<3>() { return m_vectorFields3; }

template <>
inline std::map<std::string,  VectorField4D >& TableManager::VectorFields<4>() { return m_vectorFields4; }



template <>
inline void TableManager::NewVectorField<1>(const std::string& name, const array<array<real64> >& x, const array<R1Tensor>& values)
{
  m_vectorFields1.insert(std::make_pair(name, VectorField1D()));
  VectorField1D& table = m_vectorFields1[name];
  table.SetGrid(x);
  table.SetValues(values);
  //TODO: should add automatic function addition once functions can return
  // non-scalars
}

template <>
inline void TableManager::NewVectorField<2>(const std::string& name, const array<array<real64> >& x, const array<R1Tensor>& values)
{
  m_vectorFields2.insert(std::make_pair(name, VectorField2D()));
  VectorField2D& table = m_vectorFields2[name];
  table.SetGrid(x);
  table.SetValues(values);
  //TODO: should add automatic function addition once functions can return
  // non-scalars
}

template <>
inline void TableManager::NewVectorField<3>(const std::string& name, const array<array<real64> >& x, const array<R1Tensor>& values)
{
  m_vectorFields3.insert(std::make_pair(name, VectorField3D()));
  VectorField3D& table = m_vectorFields3[name];
  table.SetGrid(x);
  table.SetValues(values);
  //TODO: should add automatic function addition once functions can return
  // non-scalars
}

template <>
inline void TableManager::NewVectorField<4>(const std::string& name, const array<array<real64> >& x, const array<R1Tensor>& values)
{
  m_vectorFields4.insert(std::make_pair(name, VectorField4D()));
  VectorField4D& table = m_vectorFields4[name];
  table.SetGrid(x);
  table.SetValues(values);
  //TODO: should add automatic function addition once functions can return
  // non-scalars
}



/// Read from a space deliminated file into a vector of vectors
/// nX nY nZ
/// v(0,0,0) v(1,0,0) v(2,0,0) .. v(nX,0,0) etc.
template < class ARRAY >
void TableManager::ReadVoxelFile(const std::string& filename, localIndex nComponents, ARRAY& values)
{
  std::ifstream inputStream(filename.c_str());
  if (inputStream)
  {
    // get dimensions of voxels
    localIndex nX, nY, nZ, ii=0;
    inputStream >> nX >> nY >> nZ;
    values.resize(nX * nY * nZ * nComponents);

    for (localIndex k = 0 ; k < nZ ; ++k)
      for (localIndex j = 0 ; j < nY ; ++j)
        for (localIndex i = 0 ; i < nX ; ++i)
          for(localIndex a = 0 ; a < nComponents ; ++a, ++ii)
            inputStream >> values[ii];
    inputStream.close();
  }
  else
  {
    throw GPException("ReadVoxelFile: Failed to load file:" + filename + " \n");
  }
}

/// Read from a space deliminated file into a vector of vectors
/// nX nY nZ nT
/// v(0,0,0,0) v(1,0,0,0) v(2,0,0,0) .. v(nX,0,0,0) etc.
template < class ARRAY >
void TableManager::ReadTimeVoxelFile(const std::string& filename, localIndex nComponents, ARRAY& values)
{
  std::ifstream inputStream(filename.c_str());
  if (inputStream)
  {
    // get dimensions of voxels
    localIndex nX, nY, nZ, nT, ii=0;
    inputStream >> nX >> nY >> nZ >> nT;
    values.resize(nX * nY * nZ * nT * nComponents);

    for (localIndex l = 0 ; l < nT ; ++l)
      for (localIndex k = 0 ; k < nZ ; ++k)
        for (localIndex j = 0 ; j < nY ; ++j)
          for (localIndex i = 0 ; i < nX ; ++i)
            for(localIndex a = 0 ; a < nComponents ; ++a, ++ii)
              inputStream >> values[ii];
    inputStream.close();
  }
  else
  {
    throw GPException("ReadTimeVoxelFile: Failed to load file:" + filename + " \n");
  }
}

/// Read from NUFT file
template < class ARRAY, class ARRAY2 >
void TableManager::ReadNUFTFile(const std::string& filename,
                                ARRAY2& x,
                                ARRAY& values)
{
  x.clear();
  x.resize(4);

  std::string ss, ss1, ss2, ss3, ss4, ss5, ss6;
  realT rr;
  int nlines;
  realT xcurr[] =
  { 0, 0, 0 };
  unsigned int icurr[] =
  { 0, 0, 0 };
  std::vector<unsigned int> dims(4);

  //read in grid and set values size
  unsigned int imin[] =
  { std::numeric_limits<unsigned int>::max(),
    std::numeric_limits<unsigned int>::max(),
    std::numeric_limits<unsigned int>::max() };
  {
    std::ifstream inputStream(filename.c_str());
    if (inputStream)
    {
      //GET YEAR AND NUMBER OF ELEMENTS TO READ PER TIME SLICE
      inputStream >> rr >> ss >> nlines >> ss1;
      x[3].push_back(rr);

      //GET HEADERS ... IGNORE
      inputStream >> ss >> ss1 >> ss2 >> ss3 >> ss4 >> ss5 >> ss6;

      set<realT> xs, ys, zs;
      for (int i = 0; i < nlines; i++)
      {
        inputStream >> icurr[0] >> icurr[1] >> icurr[2] >> xcurr[0] >> xcurr[1] >> xcurr[2] >> rr;
        for (int j = 0 ; j < 3 ; j++)
        {
          if (imin[j] > icurr[j])
            imin[j] = icurr[j];
          xs.insert(xcurr[0]);
          ys.insert(xcurr[1]);
          zs.insert(xcurr[2]);
        }
      }
      x[0].resize(xs.size());
      std::copy(xs.begin(), xs.end(), x[1].begin());
      x[1].resize(ys.size());
      std::copy(ys.begin(), ys.end(), x[2].begin());
      x[2].resize(zs.size());
      std::copy(zs.begin(), zs.end(), x[3].begin());

      //GET THE REST OF THE YEARS
      while (!inputStream.eof())
      {
        inputStream >> rr >> ss >> nlines >> ss1;
        x[3].push_back(rr);

        inputStream >> ss >> ss1 >> ss2 >> ss3 >> ss4 >> ss5 >> ss6;
        for (int i = 0 ; i < nlines ; i++)
        {
          inputStream >> icurr[0] >> icurr[1] >> icurr[2] >> xcurr[0] >> xcurr[1] >> xcurr[2] >> rr;
        }
      }
      inputStream.close();

      dims[0] = x[0].size();
      dims[1] = x[1].size();
      dims[2] = x[2].size();
      dims[3] = x[3].size();

      values.resize(dims[0] * dims[1] * dims[2] * dims[3]);
    }
    else
    {
      throw GPException("readTimeVoxelFile: Failed to load file:" + filename + " \n");
    }
  }

  //read in values
  {
    std::ifstream inputStream(filename.c_str());
    std::vector<unsigned int> ii(4);

    //GET THE VALUES
    unsigned int it = 0;
    while (!inputStream.eof())
    {
      inputStream >> rr >> ss >> nlines >> ss1;
      inputStream >> ss >> ss1 >> ss2 >> ss3 >> ss4 >> ss5 >> ss6;
      for (int i = 0 ; i < nlines ; i++)
      {
        inputStream >> icurr[0] >> icurr[1] >> icurr[2] >> xcurr[0] >> xcurr[1] >> xcurr[2] >> rr;
        ii[0] = icurr[0] - imin[0];
        ii[1] = icurr[1] - imin[1];
        ii[2] = icurr[2] - imin[2];
        ii[3] = it;
        unsigned int itable = 0;
        for (unsigned int idim = 0 ; idim < dims.size() ; idim++)
        {
          unsigned int ictable = ii[idim];
          for (unsigned int jdim = 0 ; jdim < idim ; jdim++)
            ictable *= dims[jdim];
          itable += ictable;
        }
        values[itable] = rr;
      }
      ++it;
    }
    inputStream.close();
  }
}


#endif /* TABLEMANAGER_H_ */
