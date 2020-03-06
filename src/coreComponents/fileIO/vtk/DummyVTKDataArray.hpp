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

#ifndef GEOSX_FILEIO_VTK_DUMMYVTKDATAARRAY_HPP_
#define GEOSX_FILEIO_VTK_DUMMYVTKDATAARRAY_HPP_

#include "common/DataTypes.hpp"

#include "VTKTypes.hpp"

#include <vtkAOSDataArrayTemplate.h>

namespace geosx
{
namespace vtk
{

template< typename DATA_TYPE >
class DummyVTKDataArray : public vtkDataArray
{
  public:
  DummyVTKDataArray() = default;

  virtual int GetDataType() const
  {
    return GetVTKType( std::type_index( typeid( DATA_TYPE ) ) );
  }

  virtual int GetDataTypeSize() const
  {
    return sizeof( DATA_TYPE );
  }

  virtual void Initialize()
  {
    GEOSX_ERROR("Initialize()");
  };

  virtual void SetNumberOfTuples( vtkIdType GEOSX_UNUSED_PARAM( numTuples ) )
  {
    GEOSX_ERROR("SetNumberOfTuples()");
  }

  virtual void *GetVoidPointer( vtkIdType GEOSX_UNUSED_PARAM( valueIdx ) )
  {
    GEOSX_ERROR("GetVoidPointer()");
    void * dummy = nullptr;
    return dummy;
  }

  virtual vtkTypeBool Resize( vtkIdType GEOSX_UNUSED_PARAM( numTuples) )
  {
    GEOSX_ERROR("Resize()");
    vtkTypeBool dummy = false;
    return dummy;
  }

  virtual void SetArrayFreeFunction( void ( *GEOSX_UNUSED_PARAM( callback ) )( void * ) )
  {
    GEOSX_ERROR("SetArrayFreeFunction()");
  }

  virtual vtkArrayIterator* NewIterator()
  {
    GEOSX_ERROR("NewIterator()");
    vtkArrayIterator * dummy = nullptr;
    return dummy;
  };

  virtual vtkIdType LookupValue( vtkVariant GEOSX_UNUSED_PARAM( value) )
  {
    GEOSX_ERROR("LookupValue()");
    vtkIdType dummy = 0;
    return dummy;
  }

  virtual void LookupValue( vtkVariant GEOSX_UNUSED_PARAM( value ), vtkIdList* GEOSX_UNUSED_PARAM( valueIds ) )
  {
    GEOSX_ERROR("LookupValue()");
  }

  virtual void InsertVariantValue( vtkIdType GEOSX_UNUSED_PARAM( valueIdx ),
                                   vtkVariant GEOSX_UNUSED_PARAM( value ) )
  {
    GEOSX_ERROR("InsertVariantValue()");
  }

  virtual void SetVariantValue( vtkIdType GEOSX_UNUSED_PARAM( valueIdx ),
                                vtkVariant GEOSX_UNUSED_PARAM( value) )
  {
    GEOSX_ERROR("SetVariantValue()");
  }

  virtual void DataChanged()
  {
    GEOSX_ERROR("DataChanged()");
  }

  virtual void ClearLookup()
  {
    GEOSX_ERROR("ClearLookup()");
  }

  virtual double *GetTuple( vtkIdType GEOSX_UNUSED_PARAM ( tupleIdx) )
  {
    GEOSX_ERROR("GetTuple()");
    double * dummy = nullptr;
    return dummy;
  }

  virtual void InsertTuple( vtkIdType GEOSX_UNUSED_PARAM( tupleIdx ),
                            const float * GEOSX_UNUSED_PARAM( tuple ) )
  {
    GEOSX_ERROR("InsertTuple()");
  }

  virtual void InsertTuple( vtkIdType GEOSX_UNUSED_PARAM ( tupleIdx ),
                            const double * GEOSX_UNUSED_PARAM( tuple ) )
  {
    GEOSX_ERROR("InsertTuple()");
  }

  virtual vtkIdType InsertNextTuple( const float * GEOSX_UNUSED_PARAM( tuple ) )
  {
    GEOSX_ERROR("InsertNextTuple()");
    vtkIdType dummy = 0;
    return dummy;
  }

  virtual vtkIdType InsertNextTuple( const double * GEOSX_UNUSED_PARAM( tuple ) )
  {
    GEOSX_ERROR("InsertNextTuple()");
    vtkIdType dummy = 0;
    return dummy;
  }

  virtual void RemoveTuple( vtkIdType GEOSX_UNUSED_PARAM( tupleIdx) )
  {
    GEOSX_ERROR("RemoveTuple()");
  }

  virtual void* WriteVoidPointer( vtkIdType GEOSX_UNUSED_PARAM( valueIdx ),
                                  vtkIdType GEOSX_UNUSED_PARAM( numValues ) )
  {
    GEOSX_ERROR("WriteVoidPointer()");
    void * dummy = nullptr;
    return dummy;
  }

  virtual vtkTypeBool Allocate( vtkIdType GEOSX_UNUSED_PARAM( numValues ),
                                vtkIdType GEOSX_UNUSED_PARAM( ext ) )
  {
    GEOSX_ERROR("Allocate()");
    vtkTypeBool dummy = false;
    return dummy;
  }

  virtual void Squeeze()
  {
    GEOSX_ERROR("Squeeze()");
  }

  virtual void SetVoidArray( void *vtkNotUsed(array), vtkIdType vtkNotUsed(size), int vtkNotUsed(save))
  {
    GEOSX_ERROR("SetVoidArray()");
  }
};
}
}

#endif
