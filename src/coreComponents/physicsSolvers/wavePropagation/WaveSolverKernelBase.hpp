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
 * @file WaveSolverKernelBase.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_WAVESOLVERKERNELBASE_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_WAVESOLVERKERNELBASE_HPP_
#include "finiteElement/kernelInterface/KernelBase.hpp"
#include "WaveSolverUtils.hpp"
#if !defined( GEOS_USE_HIP )
#include "finiteElement/elementFormulations/Qk_Hexahedron_Lagrange_GaussLobatto.hpp"
#endif

namespace geos
{
namespace finiteElement
{
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class WaveSolverKernelBase: public finiteElement::KernelBase< SUBREGION_TYPE,
                                                         CONSTITUTIVE_TYPE,
                                                         FE_TYPE,
                                                         1,
                                                         1 >
{
public:
  /// Alias for the base class;
  using Base = finiteElement::KernelBase< SUBREGION_TYPE,
                                          CONSTITUTIVE_TYPE,
                                          FE_TYPE,
                                          1,
                                          1 >;

  /// Maximum number of nodes per element, which is equal to the maxNumTestSupportPointPerElem and
  /// maxNumTrialSupportPointPerElem by definition. When the FE_TYPE is not a Virtual Element, this
  /// will be the actual number of nodes per element.
  static constexpr int numNodesPerElem = Base::maxNumTestSupportPointsPerElem;

  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;
  using Base::m_elemsToNodes;
  using Base::m_elemGhostRank;
  using Base::m_constitutiveUpdate;
  using Base::m_finiteElementSpace;

  /**
   * @brief Constructor
   * @copydoc geos::finiteElement::KernelBase::KernelBase
   */
  WaveSolverKernelBase( SUBREGION_TYPE const & elementSubRegion,
                        FE_TYPE const & finiteElementSpace,
                        CONSTITUTIVE_TYPE & inputConstitutiveType ):
    Base( elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType ) { }

  template< typename POLICY,
            typename KERNEL_TYPE,
            typename K = std::enable_if_t< std::is_base_of_v< WaveSolverKernelBase, KERNEL_TYPE >, KERNEL_TYPE > >
  static
  real64
  kernelLaunch( localIndex const numElems,
                K const & kernelComponent )
  {
    GEOS_MARK_FUNCTION;
    forAll< POLICY >( numElems, [=] GEOS_HOST_DEVICE ( localIndex const k )
    {
      typename K::StackVariables stack;

      kernelComponent.setup( k, stack );
      
      FE_TYPE::constForOnCell( [&] ( auto q )
      {
        kernelComponent.template quadraturePointKernel< q >( k, stack );
      } );
      kernelComponent.complete( k, stack );
    } );
    return 0;
  }
  //END_kernelLauncher
};

} // namespace finiteElement
} // namespace geos
#endif
