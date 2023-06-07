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
 * @file SolutionCheckKernel.hpp
 */

#ifndef GEOSX_SOLUTIONCHECKKERNELH_HPP
#define GEOSX_SOLUTIONCHECKKERNELH_HPP

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"


namespace geos
{

namespace compositionalMultiphaseHybridFVMKernels
{

/******************************** SolutionCheckKernel ********************************/

struct SolutionCheckKernel
{

  template< typename POLICY >
  static localIndex
  launch( arrayView1d< real64 const > const & localSolution,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & dofNumber,
          arrayView1d< integer const > const & ghostRank,
          arrayView1d< real64 const > const & facePres,
          real64 const scalingFactor )
  {
    RAJA::ReduceMin< ReducePolicy< POLICY >, integer > check( 1 );

    forAll< POLICY >( dofNumber.size(), [=] GEOS_HOST_DEVICE ( localIndex const iface )
    {
      if( ghostRank[iface] < 0 && dofNumber[iface] >= 0 )
      {
        localIndex const localRow = LvArray::integerConversion< localIndex >( dofNumber[iface] - rankOffset );
        {
          real64 const newFacePres = facePres[iface] + scalingFactor * localSolution[localRow];
          check.min( newFacePres >= 0.0 );
        }
      }
    } );
    return check.get();
  }

};

}

}

#endif //GEOSX_SOLUTIONCHECKKERNELH_HPP
