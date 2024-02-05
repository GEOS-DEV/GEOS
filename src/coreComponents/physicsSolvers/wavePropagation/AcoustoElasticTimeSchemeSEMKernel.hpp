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
 * @file AcoustoElasticTimeSchemeSEMKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTOELASTICTIMESCHEMESEMKERNEL_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTOELASTICTIMESCHEMESEMKERNEL_HPP_

namespace geos
{

struct AcoustoElasticTimeSchemeSEM
{
  
  using EXEC_POLICY = parallelDevicePolicy< >;
  using ATOMIC_POLICY = parallelDeviceAtomic;

  static void LeapFrog(real64 const dt,
                       arrayView1d< real32 > const ux_np1,
                       arrayView1d< real32 > const uy_np1,
                       arrayView1d< real32 > const uz_np1,
                       arrayView1d< real32 const > const p_n,
                       arrayView1d< real32 const > const elasticMass,
                       arrayView1d< real32 const > const atoex,
                       arrayView1d< real32 const > const atoey,
                       arrayView1d< real32 const > const atoez,
                       arrayView1d< localIndex const > const elasticFSNodeIndicator,
                       SortedArrayView< localIndex const > const interfaceNodesSet)
  {
    real64 const dt2 = pow( dt, 2 );
    forAll< EXEC_POLICY >( interfaceNodesSet.size(), [=] GEOS_HOST_DEVICE ( localIndex const n )
    {
      localIndex const a = interfaceNodesSet[n];
      if( elasticFSNodeIndicator[a] == 1 )
        return;

      real32 const aux = -p_n[a] / elasticMass[a];
      real32 const localIncrementx = dt2 * atoex[a] * aux;
      real32 const localIncrementy = dt2 * atoey[a] * aux;
      real32 const localIncrementz = dt2 * atoez[a] * aux;

      RAJA::atomicAdd< ATOMIC_POLICY >( &ux_np1[a], localIncrementx );
      RAJA::atomicAdd< ATOMIC_POLICY >( &uy_np1[a], localIncrementy );
      RAJA::atomicAdd< ATOMIC_POLICY >( &uz_np1[a], localIncrementz );
    } );

  };

};

}

#endif //GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTOELASTICTIMESCHEMESEMKERNEL_HPP_
