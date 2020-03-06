/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOSX_FILEIO_VTK_CELLVTKDATAARRAY_HPP_
#define GEOSX_FILEIO_VTK_CELLVTKDATAARRAY_HPP_

#include "DummyVTKDataArray.hpp"

#include "common/DataTypes.hpp"
#include "mesh/CellElementSubRegion.hpp"


namespace geosx
{
namespace vtk
{

class CellConnectivityVTKDataArray : public vtkAOSDataArrayTemplate<long long>
{
  public:
  CellConnectivityVTKDataArray( CellElementSubRegion const * const esr ) :
    m_esr( esr )
  {
    auto size = esr->size();
    auto dimension = esr->nodeList()[0].size();
    this->MaxId = size * dimension;
  }

  virtual void GetTuple(vtkIdType tupleIdx, double* tuple)
  {
    tuple[0] = static_cast< double >( m_esr->nodeList()[tupleIdx][0]);
  }

  virtual double * GetTuple(vtkIdType tupleIdx )
  {
  //  return m_geosxArray2d[ tupleIdx].data();
   tupleIdx++;
   double * dumy = nullptr;
   return dumy;
  }
    
  private:

  /// Encapuslated CellElementSubRegion
  CellElementSubRegion const * const m_esr;


};

class CellOffsetVTKDataArray : public vtkAOSDataArrayTemplate<long long>
{
  public:
  CellOffsetVTKDataArray( CellElementSubRegion const * const esr ) :
    m_esr( esr )
  {
    auto size = esr->size();
    this->MaxId = size;
  }

  virtual void GetTuple(vtkIdType tupleIdx, double* tuple)
  {
    tuple[0] = (tupleIdx+ 1) * 8;
  }

  virtual double * GetTuple(vtkIdType tupleIdx )
  {
  //  return m_geosxArray2d[ tupleIdx].data();
   tupleIdx++;
   double * dumy = nullptr;
   std::cout << m_esr->size() << std::endl;
   return dumy;
  }
    
  private:

  /// Encapuslated CellElementSubRegion
  CellElementSubRegion const * const m_esr;


};

}
}

#endif
