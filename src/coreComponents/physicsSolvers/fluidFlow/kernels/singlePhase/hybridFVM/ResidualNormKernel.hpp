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

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_HYBRIDFVM_RESIDUALNORMKERNEL_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_HYBRIDFVM_RESIDUALNORMKERNEL_HPP

#include "physicsSolvers/SolverBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/StencilAccessors.hpp"
#include "constitutive/solid/porosity/PorosityFields.hpp"
#include "constitutive/solid/porosity/PorosityBase.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidFields.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidBase.hpp"


namespace geos
{

namespace singlePhaseHybridFVMKernels
{

/******************************** ResidualNormKernel ********************************/

/**
 * @class ResidualNormKernel
 */
class ResidualNormKernel : public solverBaseKernels::ResidualNormKernelBase< 1 >
{
public:

  using Base = ResidualNormKernelBase< 1 >;
  using Base::minNormalizer;
  using Base::m_rankOffset;
  using Base::m_localResidual;
  using Base::m_dofNumber;

  /**
   * @brief The type for element-based non-constitutive data parameters.
   * Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  using SinglePhaseFlowAccessors =
    StencilAccessors< fields::elementVolume >;

  using SinglePhaseFluidAccessors =
    StencilMaterialAccessors< constitutive::SingleFluidBase,
                              fields::singlefluid::density_n >;
  using PorosityAccessors =
    StencilMaterialAccessors< constitutive::PorosityBase,
                              fields::porosity::porosity_n >;


  ResidualNormKernel( globalIndex const rankOffset,
                      arrayView1d< real64 const > const & localResidual,
                      arrayView1d< globalIndex const > const & dofNumber,
                      arrayView1d< localIndex const > const & ghostRank,
                      SortedArrayView< localIndex const > const & regionFilter,
                      FaceManager const & faceManager,
                      SinglePhaseFlowAccessors const & singlePhaseFlowAccessors,
                      SinglePhaseFluidAccessors const & singlePhaseFluidAccessors,
                      PorosityAccessors const & porosityAccessors,
                      real64 const & defaultViscosity,
                      real64 const & dt )
    : Base( rankOffset,
            localResidual,
            dofNumber,
            ghostRank ),
    m_dt( dt ),
    m_regionFilter( regionFilter ),
    m_defaultViscosity( defaultViscosity ),
    m_elemRegionList( faceManager.elementRegionList() ),
    m_elemSubRegionList( faceManager.elementSubRegionList() ),
    m_elemList( faceManager.elementList() ),
    m_volume( singlePhaseFlowAccessors.get( fields::elementVolume {} ) ),
    m_porosity_n( porosityAccessors.get( fields::porosity::porosity_n {} ) ),
    m_density_n( singlePhaseFluidAccessors.get( fields::singlefluid::density_n {} ) )
  {}

  GEOS_HOST_DEVICE
  void computeMassNormalizer( localIndex const kf,
                              real64 & massNormalizer,
                              real64 & multiplier ) const
  {
    integer elemCounter = 0;

    for( integer k = 0; k < m_elemRegionList.size( 1 ); ++k )
    {
      localIndex const er  = m_elemRegionList[kf][k];
      localIndex const esr = m_elemSubRegionList[kf][k];
      localIndex const ei  = m_elemList[kf][k];
      bool const onBoundary = (er == -1 || esr == -1 || ei == -1);
      bool const isInTarget = m_regionFilter.contains( er );

      // if not on boundary, increment the normalizer
      if( !onBoundary && isInTarget )
      {
        massNormalizer += m_density_n[er][esr][ei][0] * m_porosity_n[er][esr][ei][0] * m_volume[er][esr][ei];
        multiplier += m_density_n[er][esr][ei][0];
        elemCounter++;
      }
    }
    massNormalizer /= elemCounter; // average mass in the adjacent cells at the previous converged time step
    multiplier *= m_dt / elemCounter / m_defaultViscosity;  // average dt * mobility at the previous converged time step
  }

  GEOS_HOST_DEVICE
  virtual void computeLinf( localIndex const kf,
                            LinfStackVariables & stack ) const override
  {
    // if the face is adjacent to target region, compute the local values
    if( m_dofNumber[kf] >= 0 )
    {
      real64 multiplier = 0, massNormalizer = 0;
      computeMassNormalizer( kf, massNormalizer, multiplier );

      // scaled residual to be in mass units (needed because element and face residuals are blended in a single norm)
      stack.localValue[0] += LvArray::math::abs( m_localResidual[stack.localRow] * multiplier ) / LvArray::math::max( minNormalizer, massNormalizer );
    }
  }

  GEOS_HOST_DEVICE
  virtual void computeL2( localIndex const kf,
                          L2StackVariables & stack ) const override
  {
    // if the face is adjacent to target region, compute the local values
    if( m_dofNumber[kf] >= 0 )
    {
      real64 multiplier = 0, massNormalizer = 0;
      computeMassNormalizer( kf, massNormalizer, multiplier );

      // scaled residual to be in mass units (needed because element and face residuals are blended in a single norm)
      real64 const valMass = m_localResidual[stack.localRow] * multiplier;
      stack.localValue[0] += valMass * valMass;
      stack.localNormalizer[0] += massNormalizer;
    }
  }


protected:

  /// Time step size
  real64 const m_dt;

  /// Filter to identify the target regions of the solver
  SortedArrayView< localIndex const > const m_regionFilter;

  /// Value of the default viscosity
  real64 const m_defaultViscosity;

  /// Views on the maps face to elements
  arrayView2d< localIndex const > const m_elemRegionList;
  arrayView2d< localIndex const > const m_elemSubRegionList;
  arrayView2d< localIndex const > const m_elemList;

  /// View on the volume
  ElementViewConst< arrayView1d< real64 const > > const m_volume;

  /// View on porosity at the previous converged time step
  ElementViewConst< arrayView2d< real64 const > > const m_porosity_n;

  /// View on total mass density at the previous converged time step
  ElementViewConst< arrayView2d< real64 const > > const m_density_n;

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
   * @param[in] regionFilter filter to identify the target regions of the solver
   * @param[in] localResidual the residual vector on my MPI rank
   * @param[in] solverName the name of the solver
   * @param[in] elemManager reference to the element region manager
   * @param[in] faceManager reference to the face manager
   * @param[in] defaultViscosity the viscosity used for normalization
   * @param[in] dt time step size
   * @param[out] residualNorm the residual norm on the subRegion
   * @param[out] residualNormalizer the residual normalizer on the subRegion
   */
  template< typename POLICY >
  static void
  createAndLaunch( solverBaseKernels::NormType const normType,
                   globalIndex const rankOffset,
                   string const dofKey,
                   arrayView1d< real64 const > const & localResidual,
                   SortedArrayView< localIndex const > const & regionFilter,
                   string const & solverName,
                   ElementRegionManager const & elemManager,
                   FaceManager const & faceManager,
                   real64 const & defaultViscosity,
                   real64 const & dt,
                   real64 (& residualNorm)[1],
                   real64 (& residualNormalizer)[1] )
  {
    arrayView1d< globalIndex const > const dofNumber = faceManager.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const ghostRank = faceManager.ghostRank();

    using kernelType = ResidualNormKernel;
    typename kernelType::SinglePhaseFlowAccessors flowAccessors( elemManager, solverName );
    typename kernelType::SinglePhaseFluidAccessors fluidAccessors( elemManager, solverName );
    typename kernelType::PorosityAccessors poroAccessors( elemManager, solverName );

    ResidualNormKernel kernel( rankOffset, localResidual, dofNumber, ghostRank,
                               regionFilter, faceManager, flowAccessors, fluidAccessors, poroAccessors, defaultViscosity, dt );
    if( normType == solverBaseKernels::NormType::Linf )
    {
      ResidualNormKernel::launchLinf< POLICY >( faceManager.size(), kernel, residualNorm );
    }
    else // L2 norm
    {
      ResidualNormKernel::launchL2< POLICY >( faceManager.size(), kernel, residualNorm, residualNormalizer );
    }
  }

};

} // namespace singlePhaseHybridFVMKernels

} // namespace geos

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_HYBRIDFVM_RESIDUALNORMKERNEL_HPP
