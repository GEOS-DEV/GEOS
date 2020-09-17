/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file SinglePhaseBaseKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEBASEKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEBASEKERNELS_HPP

#include "common/DataTypes.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"

namespace geosx
{

namespace SinglePhaseBaseKernels
{

/******************************** MobilityKernel ********************************/

struct MobilityKernel
{
  static void
  Compute( real64 const & dens,
           real64 const & dDens_dPres,
           real64 const & visc,
           real64 const & dVisc_dPres,
           real64 & mob,
           real64 & dMob_dPres );

  static void
  Compute( real64 const & dens,
           real64 const & visc,
           real64 & mob );

  static void Launch( localIndex const size,
                      arrayView2d< real64 const > const & dens,
                      arrayView2d< real64 const > const & dDens_dPres,
                      arrayView2d< real64 const > const & visc,
                      arrayView2d< real64 const > const & dVisc_dPres,
                      arrayView1d< real64 > const & mob,
                      arrayView1d< real64 > const & dMob_dPres );

  static void Launch( SortedArrayView< localIndex const > targetSet,
                      arrayView2d< real64 const > const & dens,
                      arrayView2d< real64 const > const & dDens_dPres,
                      arrayView2d< real64 const > const & visc,
                      arrayView2d< real64 const > const & dVisc_dPres,
                      arrayView1d< real64 > const & mob,
                      arrayView1d< real64 > const & dMob_dPres );

  static void Launch( localIndex const size,
                      arrayView2d< real64 const > const & dens,
                      arrayView2d< real64 const > const & visc,
                      arrayView1d< real64 > const & mob );

  static void Launch( SortedArrayView< localIndex const > targetSet,
                      arrayView2d< real64 const > const & dens,
                      arrayView2d< real64 const > const & visc,
                      arrayView1d< real64 > const & mob );
};

/******************************** AccumulationKernel ********************************/

template< bool ISPORO >
struct AssembleAccumulationTermsHelper;

template<>
struct AssembleAccumulationTermsHelper< true >
{
  inline static constexpr void
  porosityUpdate( real64 & poro,
                  real64 & dPoro_dPres,
                  real64 const biotCoefficient,
                  real64 const poroOld,
                  real64 const bulkModulus,
                  real64 const totalMeanStress,
                  real64 const oldTotalMeanStress,
                  real64 const dPres,
                  real64 const GEOSX_UNUSED_PARAM( poroRef ),
                  real64 const GEOSX_UNUSED_PARAM( pvmult ),
                  real64 const GEOSX_UNUSED_PARAM( dPVMult_dPres ) )
  {
    dPoro_dPres = (biotCoefficient - poroOld) / bulkModulus;
    poro = poroOld + dPoro_dPres * (totalMeanStress - oldTotalMeanStress + dPres);
  }
};

template<>
struct AssembleAccumulationTermsHelper< false >
{
  inline static constexpr void
  porosityUpdate( real64 & poro,
                  real64 & dPoro_dPres,
                  real64 const GEOSX_UNUSED_PARAM( biotCoefficient ),
                  real64 const GEOSX_UNUSED_PARAM( poroOld ),
                  real64 const GEOSX_UNUSED_PARAM( bulkModulus ),
                  real64 const GEOSX_UNUSED_PARAM( totalMeanStress ),
                  real64 const GEOSX_UNUSED_PARAM( oldTotalMeanStress ),
                  real64 const GEOSX_UNUSED_PARAM( dPres ),
                  real64 const poroRef,
                  real64 const pvmult,
                  real64 const dPVMult_dPres )
  {
    poro = poroRef * pvmult;
    dPoro_dPres = dPVMult_dPres * poroRef;
  }
};


template< typename REGIONTYPE >
struct AccumulationKernel
{};

template<>
struct AccumulationKernel< CellElementSubRegion >
{


  template< bool COUPLED >
  inline static void
  Compute( real64 const & dPres,
           real64 const & densNew,
           real64 const & densOld,
           real64 const & dDens_dPres,
           real64 const & volume,
           real64 const & dVol,
           real64 const & poroRef,
           real64 const & poroOld,
           real64 const & pvMult,
           real64 const & dPVMult_dPres,
           real64 const & biotCoefficient,
           real64 const & bulkModulus,
           real64 const & totalMeanStress,
           real64 const & oldTotalMeanStress,
           real64 & poroNew,
           real64 & localAccum,
           real64 & localAccumJacobian )
  {
    real64 const volNew = volume + dVol;

    // TODO porosity update needs to be elsewhere...
    real64 dPoro_dPres;
    AssembleAccumulationTermsHelper< COUPLED >::porosityUpdate( poroNew,
                                                                dPoro_dPres,
                                                                biotCoefficient,
                                                                poroOld,
                                                                bulkModulus,
                                                                totalMeanStress,
                                                                oldTotalMeanStress,
                                                                dPres,
                                                                poroRef,
                                                                pvMult,
                                                                dPVMult_dPres );


    // Residual contribution is mass conservation in the cell
    localAccum = poroNew * densNew * volNew - poroOld * densOld * volume;

    // Derivative of residual wrt to pressure in the cell
    localAccumJacobian = (dPoro_dPres * densNew + dDens_dPres * poroNew) * volNew;
  }
};


template<>
struct AccumulationKernel< FaceElementSubRegion >
{

  template< bool COUPLED >
  inline static void
  Compute( real64 const & densNew,
           real64 const & densOld,
           real64 const & dDens_dPres,
           real64 const & volume,
           real64 const & dVol,
           real64 & localAccum,
           real64 & localAccumJacobian )
  {
    real64 const volNew = volume + dVol;

    // Residual contribution is mass conservation in the cell
    localAccum = densNew * volNew - densOld * volume;

    // Derivative of residual wrt to pressure in the cell
    localAccumJacobian =  dDens_dPres * volNew;
  }
};

} // namespace SinglePhaseBaseKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEBASEKERNELS_HPP
