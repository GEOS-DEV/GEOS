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
 * @file AquiferBCKernel.hpp
 */

#ifndef GEOSX_AQUIFERBCKERNEL_HPP
#define GEOSX_AQUIFERBCKERNEL_HPP


namespace geos
{

namespace singlePhaseFVMKernels
{

/******************************** AquiferBCKernel ********************************/

/**
 * @brief Functions to assemble aquifer boundary condition contributions to residual and Jacobian
 */
struct AquiferBCKernel
{

  /**
   * @brief The type for element-based data. Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewConstAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = geos::ElementRegionManager::ElementViewConst< VIEWTYPE >;

  GEOS_HOST_DEVICE
  static void
  compute( geos::real64 const & aquiferVolFlux,
           geos::real64 const & dAquiferVolFlux_dPres,
           geos::real64 const & aquiferDens,
           geos::real64 const & dens,
           geos::real64 const & dDens_dPres,
           geos::real64 const & dt,
           geos::real64 & localFlux,
           geos::real64 & localFluxJacobian )
  {
    if( aquiferVolFlux > 0 ) // aquifer is upstream
    {
      localFlux -= dt * aquiferVolFlux * aquiferDens;
      localFluxJacobian -= dt * dAquiferVolFlux_dPres * aquiferDens;
    }
    else // reservoir is upstream
    {
      localFlux -= dt * aquiferVolFlux * dens;
      localFluxJacobian -= dt * ( dAquiferVolFlux_dPres * dens + aquiferVolFlux * dDens_dPres );
    }
  }

  static void
  launch( geos::BoundaryStencil const & stencil,
          geos::globalIndex const rankOffset,
          ElementViewConst< geos::arrayView1d< geos::globalIndex const > > const & dofNumber,
          ElementViewConst< geos::arrayView1d< geos::integer const > > const & ghostRank,
          geos::AquiferBoundaryCondition::KernelWrapper const & aquiferBCWrapper,
          geos::real64 const & aquiferDens,
          ElementViewConst< geos::arrayView1d< geos::real64 const > > const & pres,
          ElementViewConst< geos::arrayView1d< geos::real64 const > > const & pres_n,
          ElementViewConst< geos::arrayView1d< geos::real64 const > > const & gravCoef,
          ElementViewConst< geos::arrayView2d< geos::real64 const > > const & dens,
          ElementViewConst< geos::arrayView2d< geos::real64 const > > const & dDens_dPres,
          geos::real64 const & timeAtBeginningOfStep,
          geos::real64 const & dt,
          geos::CRSMatrixView< geos::real64, geos::globalIndex const > const & localMatrix,
          geos::arrayView1d< geos::real64 > const & localRhs )
  {
    using Order = geos::BoundaryStencil::Order;

    geos::BoundaryStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
    geos::BoundaryStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
    geos::BoundaryStencil::IndexContainerViewConstType const & sefi = stencil.getElementIndices();
    geos::BoundaryStencil::WeightContainerViewConstType const & weight = stencil.getWeights();

    geos::forAll< geos::parallelDevicePolicy<> >( stencil.size(), [=] GEOS_HOST_DEVICE ( geos::localIndex const iconn )
    {
      // working variables
      geos::real64 localFlux = 0.0;
      geos::real64 localFluxJacobian = 0.0;

      geos::localIndex const er = seri( iconn, Order::ELEM );
      geos::localIndex const esr = sesri( iconn, Order::ELEM );
      geos::localIndex const ei = sefi( iconn, Order::ELEM );
      geos::real64 const areaFraction = weight( iconn, Order::ELEM );

      // compute the aquifer influx rate using the pressure influence function and the aquifer props
      geos::real64 dAquiferVolFlux_dPres = 0.0;
      geos::real64 const aquiferVolFlux = aquiferBCWrapper.compute( timeAtBeginningOfStep,
                                                                    dt,
                                                                    pres[er][esr][ei],
                                                                    pres_n[er][esr][ei],
                                                                    gravCoef[er][esr][ei],
                                                                    areaFraction,
                                                                    dAquiferVolFlux_dPres );

      // compute the phase/component aquifer flux
      AquiferBCKernel::compute( aquiferVolFlux,
                                dAquiferVolFlux_dPres,
                                aquiferDens,
                                dens[er][esr][ei][0],
                                dDens_dPres[er][esr][ei][0],
                                dt,
                                localFlux,
                                localFluxJacobian );

      // Add to residual/jacobian
      if( ghostRank[er][esr][ei] < 0 )
      {
        geos::globalIndex const globalRow = dofNumber[er][esr][ei];
        geos::localIndex const localRow = LvArray::integerConversion< geos::localIndex >( globalRow - rankOffset );
        GEOS_ASSERT_GE( localRow, 0 );
        GEOS_ASSERT_GT( localMatrix.numRows(), localRow );

        RAJA::atomicAdd( geos::parallelDeviceAtomic{}, &localRhs[localRow], localFlux );
        localMatrix.addToRow< geos::parallelDeviceAtomic >( localRow,
                                                            &dofNumber[er][esr][ei],
                                                            &localFluxJacobian,
                                                            1 );
      }
    } );
  }

};

}

}

#endif //GEOSX_AQUIFERBCKERNEL_HPP
