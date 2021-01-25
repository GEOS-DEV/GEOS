/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
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
void VTKGEOSXData::customInsertValue< R1Tensor >( localIndex index, R1Tensor const & val )
{
  for( localIndex j = 0; j < 3; j++ )
  {
    this->InsertValue( 3 * index + j, val[j] );
  }
}
}
}
