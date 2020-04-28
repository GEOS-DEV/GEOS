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

#include "VTKGEOSXData.hpp"


namespace geosx
{
namespace vtk
{
template<>
void VTKGEOSXData::CustomInsertValue< R1Tensor >( localIndex index, R1Tensor const & val )
{
  for( localIndex j = 0; j < 3; j++ )
  {
    this->InsertValue( 3 * index + j, val[j] );
  }
}

template<>
void VTKGEOSXData::CustomInsertValue< R2Tensor >( localIndex GEOSX_UNUSED_PARAM( index ), R2Tensor const & GEOSX_UNUSED_PARAM( val ) )
{
  GEOSX_WARNING( "Can't output R2Tensor with VTK export" );
}


template<>
void VTKGEOSXData::CustomInsertValue< R2SymTensor >( localIndex GEOSX_UNUSED_PARAM( index ), R2SymTensor const & GEOSX_UNUSED_PARAM ( val ) )
{
  GEOSX_WARNING( "Can't output R2Tensor with VTK export" );
}
}
}
