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

#ifndef GEOSX_KILLOUGHHYSTERESIS_HPP
#define GEOSX_KILLOUGHHYSTERESIS_HPP

#include "constitutive/ConstitutiveBase.hpp"
#include "functions/TableFunction.hpp"

#include "relativePermeability/Layouts.hpp"
#include "capillaryPressure/Layouts.hpp"



namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

/***
 * @brief KilloughHysteresis is designed to hold Killough hystereis model parameters and
 *        be in charge of all compuration related to this model (trapped Saturation,Land Coefficient?)
 */



//should be up to constitutiveBase or some new SCALConstitutiveBase but for now let's POC on relativePermeabilityBase
class KilloughHysteresis : public ConstitutiveBase
{
public:

  struct PhaseWettability
  {
    enum : integer
    {
      WETTING = 0,
      NONWETTING = 1
    };
  };

  struct HysteresisCurve_t
  {
    real64 oppositeBoundSat = -1.;
    real64 imbibitionExtremaSat = -1.;
    real64 drainageExtremaSat = -1.;

    real64 oppositeBoundValue = -1.;
    real64 imbibitionExtremaValue = -1.;
    real64 drainageExtremaValue = -1.;

    HysteresisCurve_t() = default;

    HysteresisCurve_t( std::pair< real64, real64 > const & opp, std::pair< real64, real64 > const & imbE, std::pair< real64, real64 > const & drainE )
    {
      setPoints( opp, imbE, drainE );
    }

    void setPoints( std::pair< real64, real64 > const & opp, std::pair< real64, real64 > const & imbE, std::pair< real64, real64 > const & drainE )
    {
      oppositeBoundSat = opp.first; imbibitionExtremaSat = imbE.first; drainageExtremaSat = drainE.first;
      oppositeBoundValue = opp.second; imbibitionExtremaValue = imbE.second; drainageExtremaValue = drainE.second;
    }

    bool isWetting() const
    {
      return ((drainageExtremaSat > oppositeBoundSat) ? PhaseWettability::WETTING : PhaseWettability::NONWETTING) == PhaseWettability::WETTING;
    }
    bool isZero() const
    {
      return (oppositeBoundSat<=0.0) && (imbibitionExtremaSat <= 0.0) && (drainageExtremaSat <= 0.0);
    }

  };

  KilloughHysteresis( std::string const & name, Group * const parent );

  void setRelPermParameters( real64 const & jerauldA, real64 const & jerauldB, real64 const & relpermCurv );

  static std::string catalogName() { return "KilloughHysteresis"; }

  virtual string getCatalogName() const override { return catalogName(); }

  virtual void postProcessInput() override;

//  virtual void resizeFields( localIndex const size, localIndex const numPts ) override;

  GEOSX_HOST_DEVICE
  void computeLandCoefficient( HysteresisCurve_t const & hcruve, real64 & landParam );


  class KernelKilloughHysteresisBase
  {

public:
    /// To avoid division by zero, this is the min Scrd-Scri used in the computation of the Land constant
    static constexpr real64 minScriMinusScrd = 1e-12;
    /// To avoid frequent changes from drainage to imbibition and vice versa, we use this buffer
    static constexpr real64 flowReversalBuffer = 1e-12;


    /**
     * @brief Constructor for the kernel wrapper updating Killough model parameters from relperm
     * @param[in] jerauldParam_a first (modification) parameter proposed by Jerauld
     * @param[in] jerauldParam_b second (exponent) parameter proposed by Jerauld
     * @param[in] landParam Land trapping parameter
     * @param[in] drainageMinPhaseVolFraction drainage minimum volume fraction for each phase
     * @param[in] imbibitionMinPhaseVolFraction imbibition minimum volume fraction for the wetting and non-wetting phase
     * @param[in] drainageRelPermEndPoint drainage end-point relperm for each phase
     * @param[in] imbibitionRelPermEndPoint imbibition end-point relperm for the wetting and non-wetting phase
     * @param[in] drainageMaxPhaseVolFraction drainage maximum volume fraction for each phase
     * @param[in] imbibitionMaxPhaseVolFraction imbibition maximum volume fraction for the wetting and non-wetting phase
     * @param[out] phaseTrappedVolFrac trapped saturation for each phase
     */


    KernelKilloughHysteresisBase( arrayView1d< real64 const > const & landParam,
                                  real64 const & jerauldParamA,
                                  real64 const & jerauldParamB,
                                  real64 const & killoughCurvatureParam,
                                  arrayView3d< real64, cappres::USD_CAPPRES > const & phaseTrappedVolFrac );


    /// @cond DO_NOT_DOCUMENT
    /// We need these SMFs to enable array1d< KernelWrapper > and avoid
    /// host-device errors with CUDA. Otherwise rule of 0 would be fine.

    KernelKilloughHysteresisBase() = default;
    KernelKilloughHysteresisBase( KernelKilloughHysteresisBase const & ) = default;
    KernelKilloughHysteresisBase( KernelKilloughHysteresisBase && ) = default;
    KernelKilloughHysteresisBase & operator=( KernelKilloughHysteresisBase const & ) = default;

    /**
     * @brief Function computing the trapped critical phase volume fraction (Sgcrt)
     * @param[in] Scrd the drainage critical phase volume fraction
     * @param[in] Shy the max historical phase volume fraction
     * @param[in] Smx the max phase volume fraction (= end-point phase volume fraction)
       //  * @param[in] jerauldParam_a first (modification) parameter proposed by Jerauld
       //  * @param[in] jerauldParam_b second (exponent) parameter proposed by Jerauld
     * @param[in] landParam Land trapping parameter
     * @param[out] Scrt the trapped critical phase volume fraction
     */
    GEOSX_HOST_DEVICE
    void computeTrappedCriticalPhaseVolFraction( HysteresisCurve_t const & hcurve,
                                                 real64 const & Shy,
                                                 real64 const & landParam,
                                                 real64 & Scrt ) const;


    real64 getJerauldParamA() const;
    real64 getJerauldParamB() const;
    real64 getCurvatureParam() const;

private:
//from overnested
    /// Parameter a introduced by Jerauld in the Land model
    real64 m_jerauldParam_a;
    /// Parameter b introduced by Jerauld in the Land model
    real64 m_jerauldParam_b;

    real64 m_killoughCurvatureParamRelPerm;

    //from main Relperm Class
    // Trapping parameter from the Land model (typically called C)
    arrayView1d< real64 const > m_landParam;

  };

  KernelKilloughHysteresisBase createKernelWrapper( arrayView1d< real64 const > const & landParam,
                                                    arrayView3d< real64, cappres::USD_CAPPRES > const & phaseTrappedVolFrac ) const;



  struct viewKeyStruct
  {
    static constexpr char const * jerauldParameterAString() { return "jerauldParameterA"; }
    static constexpr char const * jerauldParameterBString() { return "jerauldParameterB"; }

    static constexpr char const * killoughCurvatureParameterString() { return "killoughCurvatureParameter"; }
    static constexpr char const * killoughCurvatureCapPresParameterString() { return "killoughCurvatureCapPresParameter";}

  };

private:

  /// Parameter a introduced by Jerauld in the Land model
  real64 m_jerauldParam_a;

  /// Parameter b introduced by Jerauld in the Land model
  real64 m_jerauldParam_b;

  /// Curvature parameter introduced for wetting phase hysteresis in Killough
  real64 m_killoughCurvatureParamRelPerm;

};

//inlines

GEOSX_HOST_DEVICE
inline void
KilloughHysteresis::KernelKilloughHysteresisBase::
  computeTrappedCriticalPhaseVolFraction( HysteresisCurve_t const & hcurve,
                                          real64 const & Shy,
                                          real64 const & landParam,
                                          real64 & Scrt ) const
{

  if( hcurve.isWetting())
  {
    //unpack values
    real64 const Smxd = hcurve.drainageExtremaSat;
    real64 const Swc = hcurve.oppositeBoundSat;

    real64 const A = 1 + m_jerauldParam_a * ( Shy - Swc );
    real64 const numerator = Shy - Smxd;
    real64 const denom = A + landParam * pow( ( Smxd - Shy ) / ( Smxd - Swc ), 1 + m_jerauldParam_b/landParam );
    Scrt = Smxd + numerator / denom;
  }
  else
  {
    //unpack values
    real64 const Scrd = hcurve.drainageExtremaSat;
    real64 const Smx = hcurve.oppositeBoundSat;

    real64 const A = 1 + m_jerauldParam_a * ( Smx - Shy );
    real64 const numerator = Shy - Scrd;
    real64 const denom = A + landParam * pow( ( Shy - Scrd ) / ( Smx - Scrd ), 1 + m_jerauldParam_b / landParam );
    Scrt = LvArray::math::max( 0.0, Scrd + numerator / denom ); // trapped critical saturation from equation 2.162
  }

}


}

}

#endif //GEOSX_KILLOUGHHYSTERESIS_HPP
