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
 * @file SinglePhaseFVMKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEFVMKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEFVMKERNELS_HPP

#include "common/DataTypes.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"

namespace geosx
{

namespace SinglePhaseFVMKernels
{

/******************************** FluxKernel ********************************/

/**
 * @struct Structure to contain helper functions for the FluxKernel struct.
 */
struct FluxKernelHelper
{

  /**
   * @tparam INTEGRATION_OPTION This specifies the choice of integration rule for the aperture term
   *         in a lubrication permeability.
   * @param[in] aper0 The beginning of step aperture
   * @param[in] aper The current approximation to the end of step aperture
   * @param[out] aperTerm The resulting
   * @param[out] dAperTerm_dAper
   *
   * Typically in lubrication theory, the permeabilty involves a \f$ aperture^3 \f$ term. The
   * flow residual equation assumes a constant value for all parameters over \f$ dt \f$, which
   * may introduce significant errors given the highly nonlinear nature of the cubic aperture term.
   * The template parameter provides options:
   *  - (0) Forward Euler. This results in no non-linearity since the beginning of step aperture
   *    does not change.
   *  - (1) Exact/Simpson's Rule. This is the result of taking
   *    \f$ \int_{0}^{1} (aperture0 + (aperture-aperture0)x)x^3 dx \f$. This results
   *    in a cubic non-linearity in the resulting set of equations.
   *  .
   *  @note The use of option (1) does not imply that the time integration of the residual
   *        equation is exact, or applying Simpson's Rule. It only means that the integral of
   *        the cubic aperture term in the permeablity is exact. All other components of the
   *        residual equation are assumed constant over the integral, or use a backward
   *        Euler as the case may be. Also, we omit the use of a backward Euler option as
   *        it offers no benefit over the exact integration.
   */
  template< int INTEGRATION_OPTION >
  void static apertureForPermeablityCalculation( real64 const aper0,
                                                 real64 const aper,
                                                 real64 & aperTerm,
                                                 real64 & dAperTerm_dAper );


};

template<>
inline void
FluxKernelHelper::apertureForPermeablityCalculation< 0 >( real64 const aper0,
                                                          real64 const,
                                                          real64 & aperTerm,
                                                          real64 & dAperTerm_dAper )
{
  aperTerm = aper0*aper0*aper0;
  dAperTerm_dAper = 0.0;
}

template<>
inline void
FluxKernelHelper::apertureForPermeablityCalculation< 1 >( real64 const aper0,
                                                          real64 const aper,
                                                          real64 & aperTerm,
                                                          real64 & dAperTerm_dAper )
{
  aperTerm = 0.25 * ( aper0*aper0*aper0 +
                      aper0*aper0*aper +
                      aper0*aper*aper +
                      aper*aper*aper );

  dAperTerm_dAper = 0.25 * ( aper0*aper0 +
                             2*aper0*aper +
                             3*aper*aper );


  //printf( "aper0, aper, Kf = %4.2e, %4.2e, %4.2e \n", aper0, aper, aperTerm );
}


template<>
inline void
FluxKernelHelper::apertureForPermeablityCalculation< 2 >( real64 const,
                                                          real64 const aper,
                                                          real64 & aperTerm,
                                                          real64 & dAperTerm_dAper )
{
  aperTerm = aper*aper*aper;
  dAperTerm_dAper = 3.0*aper*aper;
}


struct FluxKernel
{
  /**
   * @brief The type for element-based non-constitutive data parameters.
   * Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementView = typename ElementRegionManager::ElementViewAccessor< VIEWTYPE >::ViewTypeConst;

  /**
   * @brief The type for element-based constitutive data parameters.
   * Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::MaterialViewAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using MaterialView = typename ElementRegionManager::MaterialViewAccessor< VIEWTYPE >::ViewTypeConst;

  /**
   * @brief launches the kernel to assemble the flux contributions to the linear system.
   * @tparam STENCIL_TYPE The type of the stencil that is being used.
   * @param[in] stencil The stencil object.
   * @param[in] dt The timestep for the integration step.
   * @param[in] dofNumber The dofNumbers for each element
   * @param[in] pres The pressures in each element
   * @param[in] dPres The change in pressure for each element
   * @param[in] gravCoef The factor for gravity calculations (g*H)
   * @param[in] dens The material density in each element
   * @param[in] dDens_dPres The change in material density for each element
   * @param[in] mob The fluid mobility in each element
   * @param[in] dMob_dPres The derivative of mobility wrt pressure in each element
   * @param[out] jacobian The linear system matrix
   * @param[out] residual The linear system residual
   */
  template< typename STENCIL_TYPE >
  static void
    Launch( STENCIL_TYPE const & stencil,
            real64 const dt,
            ElementView< arrayView1d< globalIndex const > > const & dofNumber,
            ElementView< arrayView1d< real64 const > > const & pres,
            ElementView< arrayView1d< real64 const > > const & dPres,
            ElementView< arrayView1d< real64 const > > const & gravCoef,
            ElementView< arrayView2d< real64 const > > const & dens,
            ElementView< arrayView2d< real64 const > > const & dDens_dPres,
            ElementView< arrayView1d< real64 const > > const & mob,
            ElementView< arrayView1d< real64 const > > const & dMob_dPres,
            ElementView< arrayView1d< real64 const > > const & aperture0,
            ElementView< arrayView1d< real64 const > > const & aperture,
            ElementView< arrayView1d< R1Tensor const > > const & transTMultiplier,
            R1Tensor const gravityVector,
            real64 const meanPermCoeff,
#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
            ElementView< arrayView1d< real64 const > > const & s,
            ElementView< arrayView1d< real64 const > > const & dSdAper,
#endif
            ParallelMatrix * const jacobian,
            ParallelVector * const residual,
            CRSMatrixView< real64, localIndex > const & dR_dAper );


  /**
   * @brief Compute flux and its derivatives for a given connection
   *
   * This is a general version that assumes different element regions.
   * See below for a specialized version for fluxes within a region.
   */
  inline static void
  Compute( localIndex const stencilSize,
           arraySlice1d< localIndex const > const & seri,
           arraySlice1d< localIndex const > const & sesri,
           arraySlice1d< localIndex const > const & sei,
           arraySlice1d< real64 const > const & stencilWeights,
           ElementView< arrayView1d< real64 const > > const & pres,
           ElementView< arrayView1d< real64 const > > const & dPres,
           ElementView< arrayView1d< real64 const > > const & gravCoef,
           ElementView< arrayView2d< real64 const > > const & dens,
           ElementView< arrayView2d< real64 const > > const & dDens_dPres,
           ElementView< arrayView1d< real64 const > > const & mob,
           ElementView< arrayView1d< real64 const > > const & dMob_dPres,
           real64 const dt,
           arraySlice1d< real64 > const & flux,
           arraySlice2d< real64 > const & fluxJacobian )
  {
    localIndex constexpr numElems = CellElementStencilTPFA::NUM_POINT_IN_FLUX;
    localIndex constexpr maxStencil = CellElementStencilTPFA::MAX_STENCIL_SIZE;

    stackArray1d< real64, numElems >   densWeight( numElems );
    stackArray1d< real64, maxStencil > dDensMean_dP( stencilSize );
    stackArray1d< real64, maxStencil > dFlux_dP( stencilSize );

    // clear working arrays
    dDensMean_dP = 0.0;

    // density averaging weights
    densWeight = 0.5;

    // calculate quantities on primary connected cells
    real64 densMean = 0.0;
    for( localIndex ke = 0; ke < numElems; ++ke )
    {
      // density
      real64 const density = dens[seri[ke]][sesri[ke]][sei[ke]][0];
      real64 const dDens_dP = dDens_dPres[seri[ke]][sesri[ke]][sei[ke]][0];

      // average density
      densMean        += densWeight[ke] * density;
      dDensMean_dP[ke] = densWeight[ke] * dDens_dP;
    }

    // compute potential difference MPFA-style
    real64 potDif = 0.0;
    real64 sumWeightGrav = 0.0;
    for( localIndex ke = 0; ke < stencilSize; ++ke )
    {
      localIndex const er  = seri[ke];
      localIndex const esr = sesri[ke];
      localIndex const ei  = sei[ke];

      real64 weight = stencilWeights[ke];

      real64 const gravD = gravCoef[er][esr][ei];
      real64 const gravTerm = densMean * gravD;
      sumWeightGrav += weight * gravD;

      potDif += weight * (pres[er][esr][ei] + dPres[er][esr][ei] - gravTerm);
    }

    // upwinding of fluid properties (make this an option?)
    localIndex const k_up = (potDif >= 0) ? 0 : 1;

    localIndex er_up  = seri[k_up];
    localIndex esr_up = sesri[k_up];
    localIndex ei_up  = sei[k_up];

    real64 const mobility     = mob[er_up][esr_up][ei_up];
    real64 const dMobility_dP = dMob_dPres[er_up][esr_up][ei_up];

    // compute the final flux and derivatives
    real64 const fluxVal = mobility * potDif;
    for( localIndex ke = 0; ke < stencilSize; ++ke )
    {
      real64 const weight = stencilWeights[ke];
      dFlux_dP[ke] = mobility * ( weight - dDensMean_dP[ke] * sumWeightGrav);
    }

    dFlux_dP[k_up] += dMobility_dP * potDif;

    // populate local flux vector and derivatives
    flux[0] = dt * fluxVal;
    flux[1] = -flux[0];

    for( localIndex ke = 0; ke < stencilSize; ++ke )
    {
      fluxJacobian[0][ke] = dt * dFlux_dP[ke];
      fluxJacobian[1][ke] = -fluxJacobian[0][ke];
    }
  }


  /**
   * @brief Compute flux and its derivatives for a given connection
   *.
   * This is a specialized version for fluxes within the same region.
   * See above for a general version.
   */
  inline static void
  Compute( localIndex const stencilSize,
           arraySlice1d< localIndex const > const &,
           arraySlice1d< localIndex const > const &,
           arraySlice1d< localIndex const > const & stencilElementIndices,
           arraySlice1d< real64 const > const & stencilWeights,
           arrayView1d< real64 const > const & pres,
           arrayView1d< real64 const > const & dPres,
           arrayView1d< real64 const > const & gravCoef,
           arrayView2d< real64 const > const & dens,
           arrayView2d< real64 const > const & dDens_dPres,
           arrayView1d< real64 const > const & mob,
           arrayView1d< real64 const > const & dMob_dPres,
           real64 const dt,
           arraySlice1d< real64 > const & flux,
           arraySlice2d< real64 > const & fluxJacobian )
  {
    localIndex constexpr numElems = CellElementStencilTPFA::NUM_POINT_IN_FLUX;
    localIndex constexpr maxStencil = CellElementStencilTPFA::MAX_STENCIL_SIZE;

    stackArray1d< real64, numElems > densWeight( numElems );
    stackArray1d< real64, maxStencil > dDensMean_dP( stencilSize );
    stackArray1d< real64, maxStencil > dFlux_dP( stencilSize );

    // clear working arrays
    dDensMean_dP = 0.0;

    // density averaging weights
    densWeight = 1.0 / numElems;

    // calculate quantities on primary connected cells
    real64 densMean = 0.0;
    for( localIndex i = 0; i < numElems; ++i )
    {
      // density
      real64 const density = dens[stencilElementIndices[i]][0];
      real64 const dDens_dP = dDens_dPres[stencilElementIndices[i]][0];

      // average density
      densMean += densWeight[i] * density;
      dDensMean_dP[i] = densWeight[i] * dDens_dP;
    }

    // compute potential difference MPFA-style
    real64 potDif = 0.0;
    real64 sumWeightGrav = 0.0;
    for( localIndex ke = 0; ke < stencilSize; ++ke )
    {
      localIndex const ei = stencilElementIndices[ke];
      real64 const weight = stencilWeights[ke];

      real64 const gravD = gravCoef[ei];
      real64 const gravTerm = densMean * gravD;
      sumWeightGrav += weight * gravD;
      potDif += weight * (pres[ei] + dPres[ei] - gravTerm);
    }



    // upwinding of fluid properties (make this an option?)
    localIndex const k_up = (potDif >= 0) ? 0 : 1;

    localIndex ei_up  = stencilElementIndices[k_up];

    real64 const mobility     = mob[ei_up];
    real64 const dMobility_dP = dMob_dPres[ei_up];

    // compute the final flux and derivatives
    real64 const fluxVal = mobility * potDif;
    for( localIndex ke = 0; ke < stencilSize; ++ke )
    {
      real64 const weight = stencilWeights[ke];
      dFlux_dP[ke] = mobility * ( weight - dDensMean_dP[ke] * sumWeightGrav);
    }

    dFlux_dP[k_up] += dMobility_dP * potDif;

    // populate local flux vector and derivatives
    flux[0] = dt * fluxVal;
    flux[1] = -flux[0];

    for( localIndex ke = 0; ke < stencilSize; ++ke )
    {
      fluxJacobian[0][ke] = dt * dFlux_dP[ke];
      fluxJacobian[1][ke] = -fluxJacobian[0][ke];
    }
  }



  /**
   * @brief Compute flux and its derivatives for a given multi-element connector.
   *
   * This is a specialized version that flux in a single region, and uses
   * element pairing instead of a proper junction.
   */
  inline static void
  ComputeJunction( localIndex const numFluxElems,
                   arraySlice1d< localIndex const > const & stencilElementIndices,
                   arraySlice1d< real64 const > const & stencilWeights,
                   arrayView1d< real64 const > const & pres,
                   arrayView1d< real64 const > const & dPres,
                   arrayView1d< real64 const > const & gravCoef,
                   arrayView2d< real64 const > const & dens,
                   arrayView2d< real64 const > const & dDens_dPres,
                   arrayView1d< real64 const > const & mob,
                   arrayView1d< real64 const > const & dMob_dPres,
                   arrayView1d< real64 const > const & aperture0,
                   arrayView1d< real64 const > const & aperture,
                   real64 const meanPermCoeff,
#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
                   arrayView1d< real64 const > const &,//s,
                   arrayView1d< real64 const > const &,//dSdAper,
#endif
                   real64 const dt,
                   arraySlice1d< real64 > const & flux,
                   arraySlice2d< real64 > const & fluxJacobian,
                   arraySlice2d< real64 > const & dFlux_dAperture )
  {
    real64 sumOfWeights = 0;
    real64 aperTerm[10];
    real64 dAperTerm_dAper[10];

    for( localIndex k=0; k<numFluxElems; ++k )
    {

#define PERM_CALC 1
//      real64 const aperAdd = aperture0[stencilElementIndices[k]] < 0.09e-3 ? ( 0.09e-3 -
// aperture0[stencilElementIndices[k]] ) : 0.0;
#if PERM_CALC==1
      FluxKernelHelper::
        apertureForPermeablityCalculation< 1 >( aperture0[stencilElementIndices[k]],
                                                aperture[stencilElementIndices[k]],
                                                aperTerm[k],
                                                dAperTerm_dAper[k] );
#elif PERM_CALC==2

      if( s[k] >= 1.0 )
      {
        aperTerm[k] = aperture[stencilElementIndices[k]] * aperture[stencilElementIndices[k]] * aperture[stencilElementIndices[k]];
        dAperTerm_dAper[k] = 3*aperture[stencilElementIndices[k]]*aperture[stencilElementIndices[k]];
      }
      else
      {
        aperTerm[k] = aperture[stencilElementIndices[k]] * aperture[stencilElementIndices[k]] * aperture[stencilElementIndices[k]]/s[k];
        dAperTerm_dAper[k] = 3*aperture[stencilElementIndices[k]]*aperture[stencilElementIndices[k]]/s[k]
                             - aperture[stencilElementIndices[k]] * aperture[stencilElementIndices[k]] * aperture[stencilElementIndices[k]]/(s[k]*s[k]) *
                             dSdAper[k];
      }
#endif
//      aperTerm[k] += aperAdd*aperAdd*aperAdd;


      sumOfWeights += aperTerm[k] * stencilWeights[k];
    }

    localIndex k[2];
    for( k[0]=0; k[0]<numFluxElems; ++k[0] )
    {
      for( k[1]=k[0]+1; k[1]<numFluxElems; ++k[1] )
      {
        real64 dFlux_dP[2] = {0, 0};

        localIndex const ei[2] = { stencilElementIndices[k[0]],
                                   stencilElementIndices[k[1]] };
#if 0
        real64 const weight = ( stencilWeights[k[0]]*aperTerm[k[0]] ) *
                              ( stencilWeights[k[1]]*aperTerm[k[1]] ) / sumOfWeights;

        real64 const
        dWeight_dAper[2] =
        { ( 1 / aperTerm[k[0]]  - stencilWeights[k[0]] / sumOfWeights ) * weight * dAperTerm_dAper[k[0]],
          ( 1 / aperTerm[k[1]]  - stencilWeights[k[1]] / sumOfWeights ) * weight * dAperTerm_dAper[k[1]]};
#else
        real64 const c = meanPermCoeff;

        real64 const harmonicWeight = ( stencilWeights[k[0]]*aperTerm[k[0]] ) *
                                      ( stencilWeights[k[1]]*aperTerm[k[1]] ) / sumOfWeights;

        real64 const weight = c * harmonicWeight
                              + (1.0 - c) * 0.25 * ( stencilWeights[k[0]]*aperTerm[k[0]] + stencilWeights[k[1]]*aperTerm[k[1]] );

        real64 const
        dHarmonicWeight_dAper[2] =
        { ( 1 / aperTerm[k[0]]  - stencilWeights[k[0]] / sumOfWeights ) * harmonicWeight * dAperTerm_dAper[k[0]],
          ( 1 / aperTerm[k[1]]  - stencilWeights[k[1]] / sumOfWeights ) * harmonicWeight * dAperTerm_dAper[k[1]]};

        real64 const
        dWeight_dAper[2] =
        { c * dHarmonicWeight_dAper[0] + 0.25 * ( 1.0 - c )*stencilWeights[k[0]]*dAperTerm_dAper[k[0]],
          c * dHarmonicWeight_dAper[1] + 0.25 * ( 1.0 - c )*stencilWeights[k[1]]*dAperTerm_dAper[k[1]] };

#endif
        // average density
        real64 const densMean = 0.5 * ( dens[ei[0]][0] + dens[ei[1]][0] );

        real64 const dDensMean_dP[2] = { 0.5 * dDens_dPres[ei[0]][0],
                                         0.5 * dDens_dPres[ei[1]][0] };

        real64 const potDif =  ( ( pres[ei[0]] + dPres[ei[0]] ) - ( pres[ei[1]] + dPres[ei[1]] ) -
                                 densMean * ( gravCoef[ei[0]] - gravCoef[ei[1]] ) );

        // upwinding of fluid properties (make this an option?)
        localIndex const k_up = (potDif >= 0) ? 0 : 1;

        localIndex ei_up  = stencilElementIndices[k[k_up]];

        real64 const mobility     = mob[ei_up];
        real64 const dMobility_dP = dMob_dPres[ei_up];

        // Compute flux and fill flux rval
        real64 const fluxVal = mobility * weight * potDif * dt;
        flux[k[0]] += fluxVal;
        flux[k[1]] -= fluxVal;

        // compute and fill dFlux_dP
        dFlux_dP[0] = mobility * weight * (  1 - dDensMean_dP[0] * ( gravCoef[ei[0]] - gravCoef[ei[1]] ) ) * dt;
        dFlux_dP[1] = mobility * weight * ( -1 - dDensMean_dP[1] * ( gravCoef[ei[0]] - gravCoef[ei[1]] ) ) * dt;
        dFlux_dP[k_up] += dMobility_dP * weight * potDif * dt;

        fluxJacobian[k[0]][k[0]] += dFlux_dP[0];
        fluxJacobian[k[0]][k[1]] += dFlux_dP[1];
        fluxJacobian[k[1]][k[0]] -= dFlux_dP[0];
        fluxJacobian[k[1]][k[1]] -= dFlux_dP[1];

        real64 const dFlux_dAper[2] = { mobility * dWeight_dAper[0] * potDif * dt,
                     mobility * dWeight_dAper[1] * potDif * dt };
        dFlux_dAperture[k[0]][k[0]] += dFlux_dAper[0];
        dFlux_dAperture[k[0]][k[1]] += dFlux_dAper[1];
        dFlux_dAperture[k[1]][k[0]] -= dFlux_dAper[0];
        dFlux_dAperture[k[1]][k[1]] -= dFlux_dAper[1];
      }
    }
  }

};


} // namespace SinglePhaseFVMKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEFVMKERNELS_HPP
