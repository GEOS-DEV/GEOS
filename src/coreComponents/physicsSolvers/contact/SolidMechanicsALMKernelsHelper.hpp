/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


/**
 * @file SolidMechanicsALMKernelsHelper.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSALMKERNELSHELPER_HPP_
#define GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSALMKERNELSHELPER_HPP_

#include "common/GeosxMacros.hpp"

namespace geos
{

namespace solidMechanicsALMKernelsHelper
{

template< int I_SIZE,
          int J_SIZE,
          int NUM_NODES >
GEOS_HOST_DEVICE
inline
void accumulateAtuLocalOperator( real64 ( & matrix )[I_SIZE][J_SIZE],
                                 real64 ( & N )[NUM_NODES],
                                 int const ( & perm )[NUM_NODES],
                                 real64 const detJ )
{
  //LvArray::tensorOps::fill< I_SIZE, J_SIZE >( matrix, 0 );  //make 0
  for( int a=0; a < NUM_NODES; ++a )
  {
    for( int i=0; i < I_SIZE; ++i )
    {
      matrix[i][ a*3 + i ] += -N[ perm[ a ] ] * detJ;
      matrix[i][ 3*NUM_NODES + a*3 + i ] += N[ perm[ a ] ] * detJ;
    }
  }
}

}

} // geosx

#endif /* GEOS_PHYSICSSOLVERS_CONTACT_SOLIDMECHANICSALMKERNELSHELPER_HPP_ */