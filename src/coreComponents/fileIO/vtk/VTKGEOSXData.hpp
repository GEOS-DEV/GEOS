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

#ifndef GEOSX_FILEIO_VTK_VTKGEOSXDATA_HPP_
#define GEOSX_FILEIO_VTK_VTKGEOSXDATA_HPP_

#include "common/DataTypes.hpp"

#include <vtkAOSDataArrayTemplate.h>

namespace geosx
{
namespace vtk
{

class VTKGEOSXData : public vtkAOSDataArrayTemplate< real64 >
{
  public:
    static VTKGEOSXData *New()
    {
      VTK_STANDARD_NEW_BODY(VTKGEOSXData);
    }
    template< typename T >
    void InsertValue2( localIndex index,  T const & val )
    {
      this->InsertValue( index, val );
    }
};

template<>
void VTKGEOSXData::InsertValue2< R1Tensor>( localIndex index, R1Tensor const & val );

template<>
void VTKGEOSXData::InsertValue2< R2Tensor>( localIndex index, R2Tensor const & val );

template<>
void VTKGEOSXData::InsertValue2< R2SymTensor >( localIndex index, R2SymTensor const & val );
}
}

#endif
