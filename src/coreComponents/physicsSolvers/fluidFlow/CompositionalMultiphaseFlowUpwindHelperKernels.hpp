/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file CompositionalMultiphaseFlowUpwindHelperKernels.hpp
 */
#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEFLOWUPWINDHELPERKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEFLOWUPWINDHELPERKERNELS_HPP

namespace geosx
{

namespace CompositionalMultiphaseFlowUpwindHelperKernels
{

/******************************** UpwindingHelpers ********************************/
struct UpwindHelper
  {

  /**
   * @brief The type for element-based non-constitutive data parameters.
   * Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  /**
   * @brief Form gravitational head for phase from gravity and massDensities
   * @tparam NC number of components
   * @tparam NUM_ELEMS number of elements neighbors of considered face
   * @param phaseIndex phase concerned
   * @param stencilSize number of points in the stencil
   * @param seri arraySlice of the stencil implied element region index
   * @param sesri arraySlice of the stencil implied element subregion index
   * @param sei arraySlice of the stencil implied element index
   * @param stencilWeights weights associated with elements in the stencil
   */
  template< localIndex NC, localIndex NUM_ELEMS>
  GEOSX_HOST_DEVICE
  static void formGravHead( localIndex const phaseIndex,
                            localIndex const stencilSize,
                            arraySlice1d< localIndex const > const seri,
                            arraySlice1d< localIndex const > const sesri,
                            arraySlice1d< localIndex const > const sei,
                            arraySlice1d< real64 const > const stencilWeights,
                            ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                            ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
                            ElementViewConst< arrayView3d< real64 const > > const & phaseMassDens,
                            ElementViewConst< arrayView3d< real64 const > > const & dPhaseMassDens_dPres,
                            ElementViewConst< arrayView4d< real64 const > > const & dPhaseMassDens_dComp,
                            real64& gravHead,
                            real64 (& dGravHead_dPres) [NUM_ELEMS],
                            real64 (& dGravHead_dComp) [NUM_ELEMS][NC],
                            real64 (& dProp_dComp) [NC])
  {
    //working arrays
    real64 densMean {};
    real64 dDensMean_dPres [NUM_ELEMS]{};
    real64 dDensMean_dComp [NUM_ELEMS][NC]{};

    //init
    gravHead = 0.0;
    for( localIndex i = 0; i < NUM_ELEMS; ++i )
    {
      dGravHead_dPres[i] = 0.0;
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dGravHead_dComp[i][jc] = 0.0;
        dProp_dComp[jc] = 0.0;
      }
    }

    for( localIndex i = 0; i < NUM_ELEMS; ++i )
    {
      localIndex const er  = seri[i];
      localIndex const esr = sesri[i];
      localIndex const ei  = sei[i];

      // density
      real64 const density  = phaseMassDens[er][esr][ei][0][phaseIndex];
      real64 const dDens_dPres = dPhaseMassDens_dPres[er][esr][ei][0][phaseIndex];

      applyChainRule( NC,
                      dCompFrac_dCompDens[er][esr][ei],
                      dPhaseMassDens_dComp[er][esr][ei][0][phaseIndex],
                      dProp_dComp );

      // average density and derivatives
      densMean += 0.5 * density;
      dDensMean_dPres[i] = 0.5 * dDens_dPres;
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dDensMean_dComp[i][jc] = 0.5 * dProp_dComp[jc];
      }
    }

    // compute potential difference MPFA-style
    for( localIndex i = 0; i < stencilSize; ++i )
    {
      localIndex const er  = seri[i];
      localIndex const esr = sesri[i];
      localIndex const ei  = sei[i];
      real64 const weight  = stencilWeights[i];

      real64 const gravD = weight * gravCoef[er][esr][ei];
      gravHead += densMean * gravD;

      // need to add contributions from both cells the mean density depends on
      for( localIndex j = 0; j < NUM_ELEMS; ++j )
      {
        dGravHead_dPres[j] += dDensMean_dPres[j] * gravD;
        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dGravHead_dComp[j][jc] += dDensMean_dComp[j][jc] * gravD;
        }
      }
    }

  }

  /**
   * @brief Distribute phaseFlux onto component fluxes
   * @tparam NC number of components
   * @tparam MAX_STENCIL maximum number of points in the stencil
   * @param phaseIndex concerned phase index
   * @param upwindDir upwind direction of the phase
   * @param stencilSize number of points in the stencil
   * @param seri arraySlice of the stencil implied element region index
   * @param sesri arraySlice of the stencil implied element subregion index
   * @param sei arraySlice of the stencil implied element index
   */
  template< localIndex NC, localIndex MAX_STENCIL >
  GEOSX_HOST_DEVICE
  static void
  formPhaseComp( localIndex const phaseIndex,
                 localIndex const upwindDir,
                 localIndex const stencilSize,
                 arraySlice1d< localIndex const > const seri,
                 arraySlice1d< localIndex const > const sesri,
                 arraySlice1d< localIndex const > const sei,
                 ElementViewConst< arrayView4d< real64 const > > const & phaseCompFrac,
                 ElementViewConst< arrayView4d< real64 const > > const & dPhaseCompFrac_dPres,
                 ElementViewConst< arrayView5d< real64 const > > const & dPhaseCompFrac_dComp,
                 ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
                 real64 const& phaseFlux,
                 real64 const (& dPhaseFlux_dPres) [MAX_STENCIL],
                 real64 const (& dPhaseFlux_dComp) [MAX_STENCIL][NC],
                 real64 (& compFlux) [NC],
                 real64 (& dCompFlux_dPres) [MAX_STENCIL][NC],
                 real64 (& dCompFlux_dComp) [MAX_STENCIL][NC][NC] )
  {
    /*update phaseComp from grav part*/
    localIndex const er_up = seri[upwindDir];
    localIndex const esr_up = sesri[upwindDir];
    localIndex const ei_up = sei[upwindDir];

    arraySlice1d< real64 const > phaseCompFracSub = phaseCompFrac[er_up][esr_up][ei_up][0][phaseIndex];
    arraySlice1d< real64 const > dPhaseCompFrac_dPresSub = dPhaseCompFrac_dPres[er_up][esr_up][ei_up][0][phaseIndex];
    arraySlice2d< real64 const > dPhaseCompFrac_dCompSub = dPhaseCompFrac_dComp[er_up][esr_up][ei_up][0][phaseIndex];

    real64 dProp_dC[NC] {};

    // compute component fluxes and derivatives using upstream cell composition
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      real64 const ycp = phaseCompFracSub[ic];
      compFlux[ic] += phaseFlux * ycp;

      // derivatives stemming from phase flux
      for( localIndex ke = 0; ke < stencilSize; ++ke )
      {
        dCompFlux_dPres[ke][ic] += dPhaseFlux_dPres[ke] * ycp;
        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dCompFlux_dComp[ke][ic][jc] += dPhaseFlux_dComp[ke][jc] * ycp;
        }
      }

      // additional derivatives stemming from upstream cell phase composition
      dCompFlux_dPres[upwindDir][ic] += phaseFlux * dPhaseCompFrac_dPresSub[ic];

      // convert derivatives of component fraction w.r.t. component fractions to derivatives w.r.t. component
      // densities
      applyChainRule( NC, dCompFrac_dCompDens[er_up][esr_up][ei_up], dPhaseCompFrac_dCompSub[ic], dProp_dC );
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dCompFlux_dComp[upwindDir][ic][jc] += phaseFlux * dProp_dC[jc];
      }
    }
  }


    /************************* UPWIND ******************/

    /// This enum to select the proper physics in Upwind class specialization
    enum term { Viscous, Gravity, Capillary };

  /**
   * @brief Template base class for different upwind Scheme
   * @tparam NC number of components
   * @tparam NUM_ELEMS number of elements neighbors of considered face
   * @tparam T physics concerned by the scheme if specialized
   */
  template< localIndex NC, localIndex NUM_ELEMS, term T >
  class UpwindScheme{
    protected :
        real64 pot_;
        localIndex source;

    protected:
        UpwindScheme( localIndex const GEOSX_UNUSED_PARAM(numPhase),
                      localIndex const GEOSX_UNUSED_PARAM(phaseIndex),
                      localIndex const GEOSX_UNUSED_PARAM(stencilSize),
                      arraySlice1d< localIndex const > const GEOSX_UNUSED_PARAM(seri),
                      arraySlice1d< localIndex const > const GEOSX_UNUSED_PARAM(sesri),
                      arraySlice1d< localIndex const > const GEOSX_UNUSED_PARAM(sei),
                      arraySlice1d< real64 const > const GEOSX_UNUSED_PARAM(stencilWeights),
                      real64 const GEOSX_UNUSED_PARAM(totFlux), //in fine should be a ElemnetViewConst once seq form are in place
                      ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM(pres),
                      ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM(dPres),
                      ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM(gravCoef),
                      ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM(phaseMob),
                      ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM(dCompFrac_dCompDens),
                      ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM(phaseMassDens),
                      ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM(dPhaseMassDens_dPres),
                      ElementViewConst< arrayView4d< real64 const > > const & GEOSX_UNUSED_PARAM(dPhaseMassDens_dComp) ):
                      pot_(0),source(0)
                      {};

    public:
        localIndex operator()() { return (pot_>0)? source : ((source ==0) ? 1 : 0); };
  };

/**
 * @brief Class describing the classical Phase Potential Upwind Scheme as studied in Sammon. "An analysis of upstream
 *        differencing." SPE reservoir engineering (1988)
 * @tparam NC number of components
 * @tparam NUM_ELEMS number of elements neighbors of considered face
 * @tparam T physics concerned by the scheme if specialized
 */
  template< localIndex NC, localIndex NUM_ELEMS, term T >
  class PPU : public UpwindScheme<NC,NUM_ELEMS,T>
  {
    public:
    PPU() = delete;
    PPU(localIndex const numPhase,
        localIndex const phaseIndex,
        localIndex const stencilSize,
        arraySlice1d< localIndex const > const seri,
        arraySlice1d< localIndex const > const sesri,
        arraySlice1d< localIndex const > const sei,
        arraySlice1d< real64 const > const stencilWeights,
        real64 const totFlux, //in fine should be a ElemnetViewConst once seq form are in place
        ElementViewConst< arrayView1d< real64 const > > const & pres,
        ElementViewConst< arrayView1d< real64 const > > const & dPres,
        ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
        ElementViewConst< arrayView2d< real64 const > > const & phaseMob,
        ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
        ElementViewConst< arrayView3d< real64 const > > const & phaseMassDens,
        ElementViewConst< arrayView3d< real64 const > > const & dPhaseMassDens_dPres,
        ElementViewConst< arrayView4d< real64 const > > const & dPhaseMassDens_dComp )
        :
        UpwindScheme<NC,NUM_ELEMS,T>(numPhase,phaseIndex,stencilSize,seri,sesri,sei,stencilWeights,totFlux,pres,dPres,gravCoef,
                  phaseMob,dCompFrac_dCompDens,phaseMassDens,dPhaseMassDens_dPres,dPhaseMassDens_dComp)
      {
        //compute presGrad
        real64 presGrad{};

      for( localIndex i = 0; i < stencilSize; ++i )
      {
        localIndex const er = seri[i];
        localIndex const esr = sesri[i];
        localIndex const ei = sei[i];
        real64 const weight = stencilWeights[i];

        //TODO add capillary
        presGrad += weight * (pres[er][esr][ei] + dPres[er][esr][ei]);// - capPressure);

      }

        // then form GravHead
      real64 gravHead{};
      real64 dGravHead_dP[NUM_ELEMS]{};
      real64 dGravHead_dC[NUM_ELEMS][NC]{};
      real64 dProp_dC[NC]{};

      formGravHead( phaseIndex,
                    stencilSize,
                    seri,
                    sesri,
                    sei,
                    stencilWeights,
                    gravCoef,
                    dCompFrac_dCompDens,
                    phaseMassDens,
                    dPhaseMassDens_dPres,
                    dPhaseMassDens_dComp,
                    gravHead,
                    dGravHead_dP,
                    dGravHead_dC,
                    dProp_dC);

      UpwindScheme<NC,NUM_ELEMS,T>::pot_ = presGrad - gravHead;
    }

  };

  /**
   * @brief Class describing the Hybrid Upwinding as in Lee, Efendiev, and Tchelepi. "Hybrid upwind
   *        discretization of nonlinear two-phase flow with gravity." Advances in Water Resources (2015).
   * @tparam NC number of components
   * @tparam NUM_ELEMS number of elements neighbors of considered face
   * @tparam T physics concerned enum
   */
  template< localIndex NC, localIndex NUM_ELEMS, term T>
  class HU : public UpwindScheme<NC, NUM_ELEMS, T>
    {
    public:
      HU() = delete;
      HU( localIndex const GEOSX_UNUSED_PARAM(numPhase),
          localIndex const GEOSX_UNUSED_PARAM(phaseIndex),
          localIndex const GEOSX_UNUSED_PARAM(stencilSize),
          arraySlice1d< localIndex const > const GEOSX_UNUSED_PARAM(seri),
          arraySlice1d< localIndex const > const GEOSX_UNUSED_PARAM(sesri),
          arraySlice1d< localIndex const > const GEOSX_UNUSED_PARAM(sei),
          arraySlice1d< real64 const > const GEOSX_UNUSED_PARAM(stencilWeights),
          real64 const GEOSX_UNUSED_PARAM(totFlux), //in fine should be a ElemnetViewConst once seq form are in place
          ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM(pres),
          ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM(dPres),
          ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM(gravCoef),
          ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM(phaseMob),
          ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM(dCompFrac_dCompDens),
          ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM(phaseMassDens),
          ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM(dPhaseMassDens_dPres),
          ElementViewConst< arrayView4d< real64 const > > const & GEOSX_UNUSED_PARAM(dPhaseMassDens_dComp) ) = delete;

  };


  /**
   * @brief Class describing the Potential Upwind as in Brenier and Jaffré. "Upstream differencing for
   *        multiphase flow in reservoir simulation." SIAM journal on numerical analysis (1991)
   * @tparam NC number of components
   * @tparam NUM_ELEMS number of elements neighbors of considered face
   * @tparam T physics concerned enum
   */
  template< localIndex NC, localIndex NUM_ELEMS, term T >
  class PU : public UpwindScheme<NC,NUM_ELEMS,T>
  {
    PU() = delete;
    PU( localIndex const GEOSX_UNUSED_PARAM(numPhase),
        localIndex const GEOSX_UNUSED_PARAM(phaseIndex),
        localIndex const GEOSX_UNUSED_PARAM(stencilSize),
        arraySlice1d< localIndex const > const GEOSX_UNUSED_PARAM(seri),
        arraySlice1d< localIndex const > const GEOSX_UNUSED_PARAM(sesri),
        arraySlice1d< localIndex const > const GEOSX_UNUSED_PARAM(sei),
        arraySlice1d< real64 const > const GEOSX_UNUSED_PARAM(stencilWeights),
        real64 const GEOSX_UNUSED_PARAM(totFlux), //in fine should be a ElemnetViewConst once seq form are in place
        ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM(pres),
        ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM(dPres),
        ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM(gravCoef),
        ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM(phaseMob),
        ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM(dCompFrac_dCompDens),
        ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM(phaseMassDens),
        ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM(dPhaseMassDens_dPres),
        ElementViewConst< arrayView4d< real64 const > > const & GEOSX_UNUSED_PARAM(dPhaseMassDens_dComp) ) = delete;


  };


  /**
   * @brief Specialization of PU for viscous terms
   * @tparam NC number of components
   * @tparam NUM_ELEMS number of elements neighbors of considered face
   */
template< localIndex NC, localIndex NUM_ELEMS>
class PU<NC,NUM_ELEMS,term::Viscous> : public UpwindScheme<NC,NUM_ELEMS,term::Viscous>
    {

    public:
    PU( localIndex const numPhase,
        localIndex const phaseIndex,
        localIndex const stencilSize,
        arraySlice1d< localIndex const > const seri,
        arraySlice1d< localIndex const > const sesri,
        arraySlice1d< localIndex const > const sei,
        arraySlice1d< real64 const > const stencilWeights,
        real64 const totFlux, //in fine should be a ElemnetViewConst once seq form are in place
        ElementViewConst< arrayView1d< real64 const > > const & pres,
        ElementViewConst< arrayView1d< real64 const > > const & dPres,
        ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
        ElementViewConst< arrayView2d< real64 const > > const & phaseMob,
        ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
        ElementViewConst< arrayView3d< real64 const > > const & phaseMassDens,
        ElementViewConst< arrayView3d< real64 const > > const & dPhaseMassDens_dPres,
        ElementViewConst< arrayView4d< real64 const > > const & dPhaseMassDens_dComp )
        :
        UpwindScheme<NC,NUM_ELEMS,term::Viscous>(numPhase,phaseIndex,stencilSize,seri,sesri,sei,stencilWeights,totFlux,pres,dPres,gravCoef,
                                  phaseMob,dCompFrac_dCompDens,phaseMassDens,dPhaseMassDens_dPres,dPhaseMassDens_dComp)
    {
      //Form total velocity
      UpwindScheme < NC, NUM_ELEMS, term::Viscous > ::pot_ = totFlux;
      // reinit source and target depending on gravity orientation
      // form gravHead for the currentPhase
      // then form GravHead
      real64 gravHead{};
      real64 dGravHead_dP[NUM_ELEMS]{};
      real64 dGravHead_dC[NUM_ELEMS][NC]{};
      real64 dProp_dC[NC]{};

      formGravHead( phaseIndex,
                    stencilSize,
                    seri,
                    sesri,
                    sei,
                    stencilWeights,
                    gravCoef,
                    dCompFrac_dCompDens,
                    phaseMassDens,
                    dPhaseMassDens_dPres,
                    dPhaseMassDens_dComp,
                    gravHead,
                    dGravHead_dP,
                    dGravHead_dC,
                    dProp_dC );


      localIndex const k_up = 0;
      localIndex const k_dw = 1;

      //loop other other phases to form
      for( localIndex ip = 0; ip < numPhase; ++ip )
      {

        if( ip != phaseIndex )
        {
          // then form GravHead
          real64 gravHeadOther{};
          real64 dGravHeadOther_dP[NUM_ELEMS]{};
          real64 dGravHeadOther_dC[NUM_ELEMS][NC]{};
          real64 dPropOther_dC[NC]{};

          localIndex const er_up = seri[k_up];
          localIndex const esr_up = sesri[k_up];
          localIndex const ei_up = sei[k_up];

          localIndex const er_dw = seri[k_dw];
          localIndex const esr_dw = sesri[k_dw];
          localIndex const ei_dw = sei[k_dw];

          formGravHead( ip,
                        stencilSize,
                        seri,
                        sesri,
                        sei,
                        stencilWeights,
                        gravCoef,
                        dCompFrac_dCompDens,
                        phaseMassDens,
                        dPhaseMassDens_dPres,
                        dPhaseMassDens_dComp,
                        gravHeadOther,
                        dGravHeadOther_dP,
                        dGravHeadOther_dC,
                        dPropOther_dC );

          real64 const mob_up = phaseMob[er_up][esr_up][ei_up][ip];
          real64 const mob_dw = phaseMob[er_dw][esr_dw][ei_dw][ip];

          UpwindScheme < NC, NUM_ELEMS, term::Viscous > ::pot_ +=
                                          ( gravHead - gravHeadOther >= 0 ) ? mob_dw * ( gravHead - gravHeadOther ) :
                                          mob_up * ( gravHead - gravHeadOther );

        }
      }
    }
  };

/**
 * @brief Specialization of PU  for gravitational terms
 * @tparam NC number of components
 * @tparam NUM_ELEMS number of elements neighbors of considered face
 */
  template<localIndex NC, localIndex NUM_ELEMS>
  class PU< NC,NUM_ELEMS,term::Gravity > : public UpwindScheme<NC,NUM_ELEMS,term::Gravity>
  {
  public :
    PU( localIndex const numPhase,
        localIndex const phaseIndex,
        localIndex const stencilSize,
        arraySlice1d< localIndex const > const seri,
        arraySlice1d< localIndex const > const sesri,
        arraySlice1d< localIndex const > const sei,
        arraySlice1d< real64 const > const stencilWeights,
        real64 const totFlux, //in fine should be a ElemnetViewConst once seq form are in place
        ElementViewConst< arrayView1d< real64 const > > const & pres,
        ElementViewConst< arrayView1d< real64 const > > const & dPres,
        ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
        ElementViewConst< arrayView2d< real64 const > > const & phaseMob,
        ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
        ElementViewConst< arrayView3d< real64 const > > const & phaseMassDens,
        ElementViewConst< arrayView3d< real64 const > > const & dPhaseMassDens_dPres,
        ElementViewConst< arrayView4d< real64 const > > const & dPhaseMassDens_dComp )
        :
        UpwindScheme<NC,NUM_ELEMS,term::Gravity>(numPhase,phaseIndex,stencilSize,seri,sesri,sei,stencilWeights,totFlux,pres,dPres,gravCoef,
         phaseMob,dCompFrac_dCompDens,phaseMassDens,dPhaseMassDens_dPres,dPhaseMassDens_dComp)
    {

      //Form total velocity
      UpwindScheme<NC,NUM_ELEMS,term::Gravity>::pot_ = totFlux;
      // reinit source and target depending on gravity orientation
      // form gravHead for the currentPhase
      // then form GravHead
      real64 gravHead{};
      real64 dGravHead_dP[NUM_ELEMS]{};
      real64 dGravHead_dC[NUM_ELEMS][NC]{};
      real64 dProp_dC[NC]{};

      formGravHead( phaseIndex,
                    stencilSize,
                    seri,
                    sesri,
                    sei,
                    stencilWeights,
                    gravCoef,
                    dCompFrac_dCompDens,
                    phaseMassDens,
                    dPhaseMassDens_dPres,
                    dPhaseMassDens_dComp,
                    gravHead,
                    dGravHead_dP,
                    dGravHead_dC,
                    dProp_dC);


      localIndex const k_up = 0;
      localIndex const k_dw = 1;

      //loop other other phases to form
      for(localIndex ip=0; ip<numPhase; ++ip)
      {
        if( ip!=phaseIndex ){
          // then form GravHead
          real64 gravHeadOther{};
          real64 dGravHeadOther_dP[NUM_ELEMS]{};
          real64 dGravHeadOther_dC[NUM_ELEMS][NC]{};
          real64 dPropOther_dC[NC]{};

          localIndex const er_up = seri[k_up];
          localIndex const esr_up = sesri[k_up];
          localIndex const ei_up = sei[k_up];

          localIndex const er_dw = seri[k_dw];
          localIndex const esr_dw = sesri[k_dw];
          localIndex const ei_dw = sei[k_dw];

          formGravHead( ip,
                        stencilSize,
                        seri,
                        sesri,
                        sei,
                        stencilWeights,
                        gravCoef,
                        dCompFrac_dCompDens,
                        phaseMassDens,
                        dPhaseMassDens_dPres,
                        dPhaseMassDens_dComp,
                        gravHeadOther,
                        dGravHeadOther_dP,
                        dGravHeadOther_dC,
                        dPropOther_dC);

          UpwindScheme<NC,NUM_ELEMS,term::Gravity>::source = ( gravCoef[er_up][esr_up][ei_up] > 0 ) ? 0 : 1;//heavy is upwinded with above value (resp. the lighter with below)
          real64 const mob_up = phaseMob[er_up][esr_up][ei_up][ip];
          real64 const mob_dw = phaseMob[er_dw][esr_dw][ei_dw][ip];

          UpwindScheme<NC,NUM_ELEMS,term::Gravity>::pot_ += ( gravHead - gravHeadOther >= 0 ) ? mob_dw * ( gravHead - gravHeadOther ) : mob_up * ( gravHead - gravHeadOther) ;
        }
      }
    }

  };

/**
 * @brief Specialization of HU for viscous terms
 * @tparam NC number of components
 * @tparam NUM_ELEMS number of elements neighbors of considered face
 */
   template< localIndex NC, localIndex NUM_ELEMS >
   class HU<NC,NUM_ELEMS,term::Viscous> : public UpwindScheme<NC, NUM_ELEMS, term::Viscous>
  {
    public:
    HU( localIndex const numPhase,
        localIndex const phaseIndex,
        localIndex const stencilSize,
        arraySlice1d< localIndex const > const seri,
        arraySlice1d< localIndex const > const sesri,
        arraySlice1d< localIndex const > const sei,
        arraySlice1d< real64 const > const stencilWeights,
        real64 const totFlux, //in fine should be a ElemnetViewConst once seq form are in place
        ElementViewConst< arrayView1d< real64 const > > const & pres,
        ElementViewConst< arrayView1d< real64 const > > const & dPres,
        ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
        ElementViewConst< arrayView2d< real64 const > > const & phaseMob,
        ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
        ElementViewConst< arrayView3d< real64 const > > const & phaseMassDens,
        ElementViewConst< arrayView3d< real64 const > > const & dPhaseMassDens_dPres,
        ElementViewConst< arrayView4d< real64 const > > const & dPhaseMassDens_dComp )
     :
     UpwindScheme<NC,NUM_ELEMS,term::Viscous>(numPhase,phaseIndex,stencilSize,seri,sesri,sei,stencilWeights,totFlux,pres,dPres,gravCoef,
                                               phaseMob,dCompFrac_dCompDens,phaseMassDens,dPhaseMassDens_dPres,dPhaseMassDens_dComp)
    {
        //Form total velocity
        UpwindScheme<NC,NUM_ELEMS,term::Viscous>::pot_ = totFlux;
    }

  };

/**
 * @brief Specialization of HU for gravitational terms
 * @tparam NC number of components
 * @tparam NUM_ELEMS number of elements neighbors of considered face
 */
  template< localIndex NC, localIndex NUM_ELEMS>
  class HU<NC,NUM_ELEMS,term::Gravity> : public UpwindScheme<NC, NUM_ELEMS, term::Gravity>
  {
  public:
    HU( localIndex const numPhase,
        localIndex const phaseIndex,
        localIndex const stencilSize,
        arraySlice1d< localIndex const > const seri,
        arraySlice1d< localIndex const > const sesri,
        arraySlice1d< localIndex const > const sei,
        arraySlice1d< real64 const > const stencilWeights,
        real64 const totFlux, //in fine should be a ElemnetViewConst once seq form are in place
        ElementViewConst< arrayView1d< real64 const > > const & pres,
        ElementViewConst< arrayView1d< real64 const > > const & dPres,
        ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
        ElementViewConst< arrayView2d< real64 const > > const & phaseMob,
        ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
        ElementViewConst< arrayView3d< real64 const > > const & phaseMassDens,
        ElementViewConst< arrayView3d< real64 const > > const & dPhaseMassDens_dPres,
        ElementViewConst< arrayView4d< real64 const > > const & dPhaseMassDens_dComp )
        :
        UpwindScheme<NC,NUM_ELEMS,term::Gravity>(numPhase,phaseIndex,stencilSize,seri,sesri,sei,stencilWeights,totFlux,pres,dPres,gravCoef,
          phaseMob,dCompFrac_dCompDens,phaseMassDens,dPhaseMassDens_dPres,dPhaseMassDens_dComp)
    {
      //Form total velocity
      UpwindScheme<NC,NUM_ELEMS,term::Gravity>::pot_ = 0;
      // reinit source and target depending on gravity orientation

      //form gravHead for the currentPhase
      // then form GravHead
      real64 gravHead{};
      real64 dGravHead_dP[NUM_ELEMS]{};
      real64 dGravHead_dC[NUM_ELEMS][NC]{};
      real64 dProp_dC[NC]{};

      formGravHead( phaseIndex,
                    stencilSize,
                    seri,
                    sesri,
                    sei,
                    stencilWeights,
                    gravCoef,
                    dCompFrac_dCompDens,
                    phaseMassDens,
                    dPhaseMassDens_dPres,
                    dPhaseMassDens_dComp,
                    gravHead,
                    dGravHead_dP,
                    dGravHead_dC,
                    dProp_dC);

  localIndex const k_up = 0;
  localIndex const k_dw = 1;

      //loop other other phases to form
      for(localIndex ip=0; ip<numPhase; ++ip)
      {

        if(ip!=phaseIndex){
            // then form GravHead
            real64 gravHeadOther{};
            real64 dGravHeadOther_dP[NUM_ELEMS]{};
            real64 dGravHeadOther_dC[NUM_ELEMS][NC]{};
            real64 dPropOther_dC[NC]{};

            localIndex const er_up = seri[k_up];
            localIndex const esr_up = sesri[k_up];
            localIndex const ei_up = sei[k_up];

            localIndex const er_dw = seri[k_dw];
            localIndex const esr_dw = sesri[k_dw];
            localIndex const ei_dw = sei[k_dw];

            formGravHead( ip,
                          stencilSize,
                          seri,
                          sesri,
                          sei,
                          stencilWeights,
                          gravCoef,
                          dCompFrac_dCompDens,
                          phaseMassDens,
                          dPhaseMassDens_dPres,
                          dPhaseMassDens_dComp,
                          gravHeadOther,
                          dGravHeadOther_dP,
                          dGravHeadOther_dC,
                          dPropOther_dC);

              UpwindScheme<NC,NUM_ELEMS,term::Gravity>::source = ( gravCoef[er_up][esr_up][ei_up] > 0 ) ? 0 : 1;//heavy is upwinded with above value (resp. the lighter with below)
              real64 const mob_up = phaseMob[er_up][esr_up][ei_up][ip];
              real64 const mob_dw = phaseMob[er_dw][esr_dw][ei_dw][ip];
              UpwindScheme<NC,NUM_ELEMS,term::Gravity>::pot_ += ( gravHead - gravHeadOther >= 0 ) ? mob_dw * ( gravHead - gravHeadOther ) : mob_up * ( gravHead - gravHeadOther) ;

        }
      }

    }
  };


  /*********** upwinded mobilities and fractional flow ************/

  /**
   * @brief Function returning upwinded mobility (and derivatives)  of specified phase as well as uppwind direction
   *        according to specified upwind scheme
   * @tparam NC number of components
   * @tparam UpwindScheme Desciption of how to construct potential used to decide upwind direction
   * @param numPhase total number of phases
   * @param phaseIndex concerned phase index
   * @param stencilSizei number of points in the stencil
   * @param seri arraySlice of the stencil implied element region index
   * @param sesri arraySlice of the stencil implied element subregion index
   * @param sei arraySlice of the stencil implied element index
   * @param stencilWeights weights associated with elements in the stencil
   * @param totFlux total flux signed value
   */
  template< localIndex NC, localIndex NUM_ELEMS, term T, template<localIndex,localIndex,term> class UPWIND_SCHEME >
  GEOSX_HOST_DEVICE
  static void
  upwindMob(    localIndex const numPhase,
                localIndex const phaseIndex,
                localIndex const stencilSize,
                arraySlice1d< localIndex const > const seri,
                arraySlice1d< localIndex const > const sesri,
                arraySlice1d< localIndex const > const sei,
                arraySlice1d< real64 const > const stencilWeights,
                real64 const totFlux, //in fine should be a ElemnetViewConst once seq form are in place
                ElementViewConst< arrayView1d< real64 const > > const & pres,
                ElementViewConst< arrayView1d< real64 const > > const & dPres,
                ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
                ElementViewConst< arrayView3d< real64 const > > const & phaseMassDens,
                ElementViewConst< arrayView3d< real64 const > > const & dPhaseMassDens_dPres,
                ElementViewConst< arrayView4d< real64 const > > const & dPhaseMassDens_dComp,
                ElementViewConst< arrayView2d< real64 const > > const & phaseMob,
                ElementViewConst< arrayView2d< real64 const > > const & dPhaseMob_dPres,
                ElementViewConst< arrayView3d< real64 const > > const & dPhaseMob_dComp,
                localIndex & upwindDir,
                real64  & mob,
                real64 (& dMob_dP),
                real64 (& dMob_dC)[NC] )
  {

    //reinit
    mob = 0.0;
    dMob_dP = 0.0;
    for( localIndex ic=0;ic<NC; ++ic )
      dMob_dC[ic] = 0.0;


    UPWIND_SCHEME<NC, NUM_ELEMS, T> scheme( numPhase,
                         phaseIndex,
                         stencilSize,
                         seri,
                         sesri,
                         sei,
                         stencilWeights,
                         totFlux, //in fine should be a ElementViewConst once seq form are in place
                         pres,
                         dPres,
                         gravCoef,
                         phaseMob,
                         dCompFrac_dCompDens,
                         phaseMassDens,
                         dPhaseMassDens_dPres,
                         dPhaseMassDens_dComp );

    upwindDir = scheme();

    localIndex const er_up = seri[upwindDir];
    localIndex const esr_up = sesri[upwindDir];
    localIndex const ei_up = sei[upwindDir];

    if( std::fabs(phaseMob[er_up][esr_up][ei_up][phaseIndex])>1e-20 )
    {
      mob = phaseMob[er_up][esr_up][ei_up][phaseIndex];
      dMob_dP = dPhaseMob_dPres[er_up][esr_up][ei_up][phaseIndex];
      for( localIndex ic = 0; ic < NC; ++ic )
        dMob_dC[ic] = dPhaseMob_dComp[er_up][esr_up][ei_up][phaseIndex][ic];
    }
  }

  /**
    * @brief Function returning upwinded fractional flow (and derivatives) as well as upwind direction of specified phase
    *        according to specified upwind scheme
    * @tparam NC number of components
    * @tparam UpwindScheme Desciption of how to construct potential used to decide upwind direction
    * @param numPhase total number of phases
    * @param phaseIndex concerned phase index
    * @param stencilSizei number of points in the stencil
    * @param seri arraySlice of the stencil implied element region index
    * @param sesri arraySlice of the stencil implied element subregion index
    * @param sei arraySlice of the stencil implied element index
    * @param stencilWeights weights associated with elements in the stencil
    * @param totFlux total flux signed value
    */
  template< localIndex NC, localIndex NUM_ELEMS, localIndex MAX_STENCIL, term T, template<localIndex,localIndex,term> class UPWIND_SCHEME >
  GEOSX_HOST_DEVICE
  static void
  formFracFlow( localIndex const numPhase,
                localIndex const phaseIndex,
                localIndex const stencilSize,
                arraySlice1d< localIndex const > const seri,
                arraySlice1d< localIndex const > const sesri,
                arraySlice1d< localIndex const > const sei,
                arraySlice1d< real64 const > const stencilWeights,
                real64 const totFlux, //in fine should be a ElemnetViewConst once seq form are in place
                ElementViewConst< arrayView1d< real64 const > > const & pres,
                ElementViewConst< arrayView1d< real64 const > > const & dPres,
                ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
                ElementViewConst< arrayView3d< real64 const > > const & phaseMassDens,
                ElementViewConst< arrayView3d< real64 const > > const & dPhaseMassDens_dPres,
                ElementViewConst< arrayView4d< real64 const > > const & dPhaseMassDens_dComp,
                ElementViewConst< arrayView2d< real64 const > > const & phaseMob,
                ElementViewConst< arrayView2d< real64 const > > const & dPhaseMob_dPres,
                ElementViewConst< arrayView3d< real64 const > > const & dPhaseMob_dComp,
                localIndex& k_up_main,
                real64  & fflow,
                real64 (& dFflow_dP) [MAX_STENCIL],
                real64 (& dFflow_dC) [MAX_STENCIL][NC] )
  {
    // get var to memorized the numerator mobility properly upwinded
    real64 mainMob {};
    real64 dMMob_dP{};
    real64 dMMob_dC[NC] {};

    real64 totMob {};
    real64 dTotMob_dP[MAX_STENCIL] {};
    real64 dTotMob_dC[MAX_STENCIL][NC] {};

    //reinit
    //fractional flow too low to let the upstream phase flow
    k_up_main = -1; //to throw error if unmodified
    fflow = 0;
    for( localIndex ke = 0; ke < stencilSize; ++ke )
    {
      dFflow_dP[ke] = 0;
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dFflow_dC[ke][jc] = 0;
      }
    }

    //Form totMob
    for(localIndex ip=0; ip<numPhase ; ++ip)
    {

      localIndex k_up;
      real64 mob {};
      real64 dMob_dP {};
      real64 dMob_dC[NC] {};

      upwindMob<NC, NUM_ELEMS, T, UPWIND_SCHEME>( numPhase,
                                   ip,
                                   stencilSize,
                                   seri,
                                   sesri,
                                   sei,
                                   stencilWeights,
                                   totFlux, //in fine should be a ElemnetViewConst once seq form are in place
                                   pres,
                                   dPres,
                                   gravCoef,
                                   dCompFrac_dCompDens,
                                   phaseMassDens,
                                   dPhaseMassDens_dPres,
                                   dPhaseMassDens_dComp,
                                   phaseMob,
                                   dPhaseMob_dPres,
                                   dPhaseMob_dComp,
                                   k_up,
                                   mob,
                                   dMob_dP,
                                   dMob_dC);


      totMob += mob;
      dTotMob_dP[k_up] += dMob_dP;
      for(localIndex ic=0; ic< NC; ++ic)
        dTotMob_dC[k_up][ic] += dMob_dC[ic];

      if(ip == phaseIndex)
      {
        k_up_main = k_up;
        mainMob = mob;
        dMMob_dP = dMob_dP;
        for(localIndex ic=0; ic< NC; ++ic)
          dMMob_dC[ic] = dMob_dC[ic];
      }
    }

    //guard against no flow region
    if( std::fabs(mainMob) > 1e-20 )
    {
      fflow = mainMob / totMob;
      dFflow_dP[k_up_main] = dMMob_dP/ totMob;
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dFflow_dC[k_up_main][jc] = dMMob_dC[jc] / totMob;

      }

      for( localIndex ke = 0; ke < stencilSize; ++ke )
      {
        dFflow_dP[ke] -= fflow * dTotMob_dP[ke] / totMob;

        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dFflow_dC[ke][jc] -= fflow * dTotMob_dC[ke][jc] / totMob;
        }
      }
    }
  }


  };//end of struct UpwindHelper



}//end namespace CompositionalMultiphaseFlowUpwindHelperKernels

}//end namespace geosx
#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEFLOWUPWINDHELPERKERNELS_HP
