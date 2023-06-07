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

#ifndef GEOSX_SOLUTIONCHECKKERNEL_HPP
#define GEOSX_SOLUTIONCHECKKERNEL_HPP

#include "physicsSolvers/SolverBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/fields/FlowSolverBaseFields.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "constitutive/solid/CoupledSolidBase.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "common/DataTypes.hpp"

namespace geos
{

namespace singlePhaseBaseKernels
{

/******************************** SolutionCheckKernel ********************************/

struct SolutionCheckKernel
{
  template< typename POLICY >
  static geos::localIndex launch( geos::arrayView1d< geos::real64 const > const & localSolution,
                                  geos::globalIndex const rankOffset,
                                  geos::arrayView1d< geos::globalIndex const > const & dofNumber,
                                  geos::arrayView1d< geos::integer const > const & ghostRank,
                                  geos::arrayView1d< geos::real64 const > const & pres,
                                  geos::real64 const scalingFactor )
  {
    RAJA::ReduceMin< geos::ReducePolicy< POLICY >, geos::localIndex > minVal( 1 );

    geos::forAll< POLICY >( dofNumber.size(), [=] GEOS_HOST_DEVICE ( geos::localIndex const ei )
    {
      if( ghostRank[ei] < 0 && dofNumber[ei] >= 0 )
      {
        geos::localIndex const lid = dofNumber[ei] - rankOffset;
        geos::real64 const newPres = pres[ei] + scalingFactor * localSolution[lid];

        if( newPres < 0.0 )
        {
          minVal.min( 0 );
        }
      }

    } );

    return minVal.get();
  }

};

}

}

#endif //GEOSX_SOLUTIONCHECKKERNEL_HPP
