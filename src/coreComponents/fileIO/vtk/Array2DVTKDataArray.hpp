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

#ifndef GEOSX_FILEIO_VTK_ARRAY2DVTKDATAARRAY_HPP_
#define GEOSX_FILEIO_VTK_ARRAY2DVTKDATAARRAY_HPP_

#include "DummyVTKDataArray.hpp"

#include "common/DataTypes.hpp"

#include<vtkArrayIterator.h>
#include <vtkArrayIteratorTemplate.h>


namespace geosx
{
namespace vtk
{

template< typename DATA_TYPE, int USD >
class Array2DVTKDataArray : public DummyVTKDataArray< DATA_TYPE >
{
  public:
  Array2DVTKDataArray(arrayView2d< DATA_TYPE, USD > const & geosxArray2D ) :
    m_geosxArray2d( geosxArray2D )
  {
    auto size = geosxArray2D.dims()[0];
    auto dimension = geosxArray2D.dims()[1];
    std::cout << "size "<< size << std::endl;
    std::cout << "dim "<< dimension << std::endl;
    std::cout << "total ?? "<< geosxArray2D.size() << std::endl;
    this->SetNumberOfComponents( dimension );
    //this->SetNumberOfValues( size * dimension);
    this->MaxId = size * dimension;
  }
  virtual void GetTuple(vtkIdType tupleIdx, double* tuple)
  {
    std::cout << tupleIdx << std::endl;
    auto dimension = m_geosxArray2d.dims()[1];
    for( auto i = 0; i < dimension; ++i )
    {
      tuple[i] = m_geosxArray2d[tupleIdx][i];
    }
  }
  virtual double * GetTuple(vtkIdType tupleIdx )
  {
  //  return m_geosxArray2d[ tupleIdx].data();
   tupleIdx++;
   double * dumy = nullptr;
   return dumy;
  }
  virtual vtkArrayIterator* NewIterator()
  {
    vtkArrayIterator* iter = vtkArrayIteratorTemplate<double>::New();
    iter->Initialize( this );
    return iter;
  };
  private:

  /// A view to the encapsulated array2d
  arrayView2d< DATA_TYPE, USD > const & m_geosxArray2d;


};

}
}

#endif
