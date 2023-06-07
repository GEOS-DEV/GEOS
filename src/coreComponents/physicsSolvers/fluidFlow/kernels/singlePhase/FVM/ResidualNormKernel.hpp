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
 * @file ResidualNormKernel.hpp
 */

#ifndef GEOSX_RESIDUALNORMKERNEL_HPP
#define GEOSX_RESIDUALNORMKERNEL_HPP

#include "physicsSolvers/SolverBaseKernels.hpp"


namespace geos
{

namespace singlePhaseBaseKernels
{

/******************************** ResidualNormKernel ********************************/

/**
 * @class ResidualNormKernel
 */
class ResidualNormKernel : public geos::solverBaseKernels::ResidualNormKernelBase< 1 >
{
public:

  using Base = ResidualNormKernelBase< 1 >;
  using Base::minNormalizer;
  using Base::m_rankOffset;
  using Base::m_localResidual;
  using Base::m_dofNumber;

  ResidualNormKernel( geos::globalIndex const rankOffset,
                      geos::arrayView1d< geos::real64 const > const & localResidual,
                      geos::arrayView1d< geos::globalIndex const > const & dofNumber,
                      geos::arrayView1d< geos::localIndex const > const & ghostRank,
                      geos::ElementSubRegionBase const & subRegion,
                      geos::constitutive::SingleFluidBase const & fluid,
                      geos::constitutive::CoupledSolidBase const & solid )
    : Base( rankOffset,
            localResidual,
            dofNumber,
            ghostRank ),
    m_volume( subRegion.getElementVolume() ),
    m_porosity_n( solid.getPorosity_n() ),
    m_density_n( fluid.density_n() )
  { }

  GEOS_HOST_DEVICE
  virtual void computeLinf( geos::localIndex const ei,
                            LinfStackVariables & stack ) const override
  {
    geos::real64 const massNormalizer = LvArray::math::max( minNormalizer, m_density_n[ei][0] * m_porosity_n[ei][0] * m_volume[ei] );
    geos::real64 const valMass = LvArray::math::abs( m_localResidual[stack.localRow] ) / massNormalizer;
    if( valMass > stack.localValue[0] )
    {
      stack.localValue[0] = valMass;
    }
  }

  GEOS_HOST_DEVICE
  virtual void computeL2( geos::localIndex const ei,
                          L2StackVariables & stack ) const override
  {
    geos::real64 const massNormalizer = LvArray::math::max( minNormalizer, m_density_n[ei][0] * m_porosity_n[ei][0] * m_volume[ei] );
    stack.localValue[0] += m_localResidual[stack.localRow] * m_localResidual[stack.localRow];
    stack.localNormalizer[0] += massNormalizer;
  }


protected:

  /// View on the volume
  geos::arrayView1d< geos::real64 const > const m_volume;

  /// View on porosity at the previous converged time step
  geos::arrayView2d< geos::real64 const > const m_porosity_n;

  /// View on total mass/molar density at the previous converged time step
  geos::arrayView2d< geos::real64 const > const m_density_n;

};

/**
 * @class ResidualNormKernelFactory
 */
class ResidualNormKernelFactory
{
public:

  /**
   * @brief Create a new kernel and launch
   * @tparam POLICY the policy used in the RAJA kernel
   * @param[in] normType the type of norm used (Linf or L2)
   * @param[in] rankOffset the offset of my MPI rank
   * @param[in] dofKey the string key to retrieve the degress of freedom numbers
   * @param[in] localResidual the residual vector on my MPI rank
   * @param[in] subRegion the element subregion
   * @param[in] fluid the fluid model
   * @param[in] solid the solid model
   * @param[out] residualNorm the residual norm on the subRegion
   * @param[out] residualNormalizer the residual normalizer on the subRegion
   */
  template< typename POLICY >
  static void
  createAndLaunch( geos::solverBaseKernels::NormType const normType,
                   geos::globalIndex const rankOffset,
                   geos::string const dofKey,
                   geos::arrayView1d< geos::real64 const > const & localResidual,
                   geos::ElementSubRegionBase const & subRegion,
                   geos::constitutive::SingleFluidBase const & fluid,
                   geos::constitutive::CoupledSolidBase const & solid,
                   geos::real64 (& residualNorm)[1],
                   geos::real64 (& residualNormalizer)[1] )
  {
    geos::arrayView1d< geos::globalIndex const > const dofNumber = subRegion.getReference< geos::array1d< geos::globalIndex > >( dofKey );
    geos::arrayView1d< geos::integer const > const ghostRank = subRegion.ghostRank();

    ResidualNormKernel kernel( rankOffset, localResidual, dofNumber, ghostRank, subRegion, fluid, solid );
    if( normType == geos::solverBaseKernels::NormType::Linf )
    {
      ResidualNormKernel::launchLinf< POLICY >( subRegion.size(), kernel, residualNorm );
    }
    else // L2 norm
    {
      ResidualNormKernel::launchL2< POLICY >( subRegion.size(), kernel, residualNorm, residualNormalizer );
    }
  }

};

}

}

#endif //GEOSX_RESIDUALNORMKERNEL_HPP
