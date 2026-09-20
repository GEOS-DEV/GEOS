/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SolutionScalingKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_REACTIVE_SOLUTIONSCALINGKERNEL_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_REACTIVE_SOLUTIONSCALINGKERNEL_HPP

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geos
{

namespace singlePhaseReactiveBaseKernels
{

/**
 * @brief Bound the change in log primary species concentration over a Newton iteration.
 */
struct SolutionScalingKernel
{
  /**
   * @brief Compute the factor that keeps the largest log concentration change within the bound.
   * @tparam POLICY the execution policy
   * @param localSolution the solution vector
   * @param rankOffset the offset of this rank in the solution vector
   * @param dofNumber the first degree of freedom of each element
   * @param ghostRank the ghost rank of each element
   * @param numDofPerCell the number of degrees of freedom per cell
   * @param speciesOffset the offset of the first species within a cell, 1 or 2 when thermal
   * @param numPrimarySpecies the number of primary species
   * @param maxAbsoluteLogConcChange the bound, in ln units. Zero or less disables the scaling.
   * @return the scaling factor and the largest log concentration change seen, before scaling
   */
  template< typename POLICY >
  static std::pair< real64, real64 > launch( arrayView1d< real64 const > const & localSolution,
                                             globalIndex const rankOffset,
                                             arrayView1d< globalIndex const > const & dofNumber,
                                             arrayView1d< integer const > const & ghostRank,
                                             integer const numDofPerCell,
                                             integer const speciesOffset,
                                             integer const numPrimarySpecies,
                                             real64 const maxAbsoluteLogConcChange )
  {
    GEOS_UNUSED_VAR( numDofPerCell );

    RAJA::ReduceMin< ReducePolicy< POLICY >, real64 > scalingFactor( 1.0 );
    RAJA::ReduceMax< ReducePolicy< POLICY >, real64 > maxDeltaLogConc( 0.0 );

    forAll< POLICY >( dofNumber.size(), [=] GEOS_HOST_DEVICE ( localIndex const ei )
    {
      if( ghostRank[ei] < 0 && dofNumber[ei] >= 0 )
      {
        localIndex const lid = dofNumber[ei] - rankOffset;

        for( integer ic = 0; ic < numPrimarySpecies; ++ic )
        {
          real64 const absLogConcChange = LvArray::math::abs( localSolution[lid + speciesOffset + ic] );
          maxDeltaLogConc.max( absLogConcChange );

          if( maxAbsoluteLogConcChange > 0.0 && absLogConcChange > maxAbsoluteLogConcChange )
          {
            scalingFactor.min( maxAbsoluteLogConcChange / absLogConcChange );
          }
        }
      }

    } );

    return { scalingFactor.get(), maxDeltaLogConc.get() };
  }

};

} // namespace singlePhaseReactiveBaseKernels

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_REACTIVE_SOLUTIONSCALINGKERNEL_HPP
