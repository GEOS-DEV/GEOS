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
 * @file TableRelativePermeabilityHysteresis.hpp
 */

#ifndef GEOSX_KILLOUGHHYSTERESIS_HPP
#define GEOSX_KILLOUGHHYSTERESIS_HPP

#include "relativePermeability/Layouts.hpp"
#include "capillaryPressure/Layouts.hpp"

#include "constitutive/ConstitutiveBase.hpp"
#include "functions/TableFunction.hpp"


namespace geosx
{

namespace constitutive
{

/***
 * @brief KilloughHysteresis is designed to hold Killough hystereis model parameters and
 *        be in charge of all compuration related to this model (trapped Saturation,Land Coefficient?)
 */



//should be up to constitutiveBase or some new SCALConstitutiveBase but for now let's POC on relativePermeabilityBase
class KilloughHysteresis
{
public:

  static constexpr real64 minScriMinusScrd = 1e-12;
  /// To avoid frequent changes from drainage to imbibition and vice versa, we use this buffer
  static constexpr real64 flowReversalBuffer = 1e-12;

  struct PhaseWettability
  {
    enum : integer
    {
      WETTING = 0,
      NONWETTING = 1
    };
  };

  /**
   * @brief  struct to represent hysteresis curves (relperm or capillary pressure)
   * whether they are wetting or non wetting, storing key point as pairs of
   * saturations and value, being either the relperm value (S,kr) or capillary pressure value (S,pc).
   * this way we can tell wetting curve from non wetting by the ordering of drainage/imbibition key values.
   * Indeed if imbibition comes at lower saturation than drainage then it is wetting curve, on the opposite
   * if drainage comes before imbibition this is a non-wetting hysteresis. This is completed by an opposite
   * point that is either the connate wetting saturation Swc or the maximum non wetting saturation Sgmax.
   * @param  oppositeBoundPhaseVolFraction represents either Swc or Sgmax depending if a wetting curve or nonwetting is described
   * @param imbibitionExtremaPhaseVolFraction represents in wetting case the imibibition max and in non-wetting the imbibition residual
   * @param drainageExtremaPhaseVolFraction represents in wetting case the drainage max and in non-wetting the drainage residual
   * @param oppositeBoundSCALValue represents the associate relperm or capillary pressure value
   * @param imbibitionExtremaSCALValue represents the associate relperm or capillary pressure value
   * @param drainageExtremaSCALValue represents the associate relperm or capillary pressure value
   */
  struct HysteresisCurve
  {
    real64 oppositeBoundPhaseVolFraction = -1.;
    real64 imbibitionExtremaPhaseVolFraction = -1.;
    real64 drainageExtremaPhaseVolFraction = -1.;

    real64 oppositeBoundSCALValue = -1.;
    real64 imbibitionExtremaSCALValue = -1.;
    real64 drainageExtremaSCALValue = -1.;

    HysteresisCurve() = default;

    HysteresisCurve( std::pair< real64, real64 > const & opp, std::pair< real64, real64 > const & imbE, std::pair< real64, real64 > const & drainE )
    {
      setPoints( opp, imbE, drainE );
    }

    void setPoints( std::pair< real64, real64 > const & opp, std::pair< real64, real64 > const & imbE, std::pair< real64, real64 > const & drainE )
    {
      oppositeBoundPhaseVolFraction = opp.first;
      imbibitionExtremaPhaseVolFraction = imbE.first;
      drainageExtremaPhaseVolFraction = drainE.first;

      oppositeBoundSCALValue = opp.second;
      imbibitionExtremaSCALValue = imbE.second;
      drainageExtremaSCALValue = drainE.second;
    }

    GEOSX_HOST_DEVICE
    inline bool isWetting() const
    {
      return ((drainageExtremaPhaseVolFraction > oppositeBoundPhaseVolFraction) ? PhaseWettability::WETTING : PhaseWettability::NONWETTING) == PhaseWettability::WETTING;
    }


    bool isZero() const
    {
      return (oppositeBoundPhaseVolFraction <= 0.0) && (imbibitionExtremaPhaseVolFraction <= 0.0) && (drainageExtremaPhaseVolFraction <= 0.0);
    }

  };

//  void setRelPermParameters( real64 const & jerauldA, real64 const & jerauldB, real64 const & relpermCurv );

  static std::string catalogName() { return "KilloughHysteresis"; }

  static void postProcessInput( real64 const & jerauldParam_a, real64 const & jerauldParam_b,
                                real64 const & killoughCurvatureParamRelPerm );

  GEOSX_HOST_DEVICE
  static void computeLandCoefficient( HysteresisCurve const & hcruve, real64 & landParam );
  /**
   * @brief Function computing the trapped critical phase volume fraction
   * @param[in] hcurve the hysteresis curve to be used and dispatched on
   * @param[in] Shy the max historical phase volume fraction
   * @param[in] landParam Land trapping parameter
   * @param[in] jerauldParam_a jerauld expononent
   * @param[in] jerauldParam_b jerauld expononent
   * @param[out] Scrt the trapped critical phase volume fraction
   */
  GEOSX_HOST_DEVICE
  static void computeTrappedCriticalPhaseVolFraction( HysteresisCurve const & hcurve,
                                                      real64 const & Shy,
                                                      real64 const & landParam,
                                                      real64 const & jerauldParam_a,
                                                      real64 const & jerauldParam_b,
                                                      real64 & Scrt );

  struct viewKeyStruct
  {
    static constexpr char const * jerauldParameterAString() { return "jerauldParameterA"; }
    static constexpr char const * jerauldParameterBString() { return "jerauldParameterB"; }

    static constexpr char const * killoughCurvatureParameterString() { return "killoughCurvatureParameter"; }
  };

};

        GEOSX_HOST_DEVICE
        inline void KilloughHysteresis::computeLandCoefficient( KilloughHysteresis::HysteresisCurve const & hcurve,
                                                         real64 & landParam )
        {

            // Note: for simplicity, the notations are taken from IX documentation (although this breaks our phaseVolFrac naming convention)

            // Step 1: Land parameter for the wetting phase
            if( hcurve.isWetting() )
            {
                real64 const Scrd = hcurve.oppositeBoundPhaseVolFraction;
                real64 const Smxd = hcurve.drainageExtremaPhaseVolFraction;
                real64 const Smxi = hcurve.imbibitionExtremaPhaseVolFraction;
                real64 const Swc = Scrd;
                GEOSX_ERROR_IF(  (Smxi - Smxd) > 0,
                                 GEOSX_FMT( "{}: For wetting phase hysteresis, imbibition end-point saturation Smxi( {} ) must be smaller than the drainage saturation end-point Smxd( {} ).\n"
                                            "Crossing relative permeability curves.\n",
                                            catalogName(),
                                            Smxi,
                                            Smxd ));

                landParam = ( Smxd - Swc ) / LvArray::math::max( KilloughHysteresis::minScriMinusScrd, ( Smxd - Smxi ) ) - 1.0;
            }
            else
                // Step 2: Land parameter for the non-wetting phase

            {
                real64 const Smx =  hcurve.oppositeBoundPhaseVolFraction;
                real64 const Scrd = hcurve.drainageExtremaPhaseVolFraction;
                real64 const Scri = hcurve.imbibitionExtremaPhaseVolFraction;
                GEOSX_ERROR_IF( (Scrd - Scri) > 0,
                                GEOSX_FMT( "{}: For non-wetting phase hysteresis, drainage trapped saturation Scrd( {} ) must be smaller than the imbibition saturation Scri( {} ).\n"
                                           "Crossing relative permeability curves.\n",
                                           catalogName(),
                                           Scrd,
                                           Scri ));

                landParam = ( Smx - Scrd ) / LvArray::math::max( KilloughHysteresis::minScriMinusScrd, ( Scri - Scrd ) ) - 1.0;
            }
        }

        GEOSX_HOST_DEVICE
        inline void
        KilloughHysteresis::
        computeTrappedCriticalPhaseVolFraction( HysteresisCurve const & hcurve,
                                                real64 const & Shy,
                                                real64 const & landParam,
                                                real64 const & jerauldParam_a,
                                                real64 const & jerauldParam_b,
                                                real64 & Scrt )
        {

            if( hcurve.isWetting())
            {
                //unpack values
                real64 const Smxd = hcurve.drainageExtremaPhaseVolFraction;
                real64 const Swc = hcurve.oppositeBoundPhaseVolFraction;

                real64 const A = 1 + jerauldParam_a * (Shy - Swc);
                real64 const numerator = Shy - Smxd;
                real64 const denom = A + landParam * pow((Smxd - Shy) / (Smxd - Swc), 1 + jerauldParam_b / landParam );
                Scrt = Smxd + numerator / denom;
            }
            else
            {
                //unpack values
                real64 const Scrd = hcurve.drainageExtremaPhaseVolFraction;
                real64 const Smx = hcurve.oppositeBoundPhaseVolFraction;

                real64 const A = 1 + jerauldParam_a * (Smx - Shy);
                real64 const numerator = Shy - Scrd;
                real64 const denom = A + landParam * pow((Shy - Scrd) / (Smx - Scrd), 1 + jerauldParam_b / landParam );
                Scrt = LvArray::math::max( 0.0,
                                           Scrd + numerator / denom );           // trapped critical saturation from equation 2.162
            }

        }

}

}

#endif //GEOSX_KILLOUGHHYSTERESIS_HPP
