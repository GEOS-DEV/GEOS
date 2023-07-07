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
 * @file KilloughHysteresis.hpp
 */

#ifndef GEOS_CONSTITUTIVE_RELATIVEPERMEABILITY_KILLOUGHHYSTERESIS_HPP
#define GEOS_CONSTITUTIVE_RELATIVEPERMEABILITY_KILLOUGHHYSTERESIS_HPP

#include "Layouts.hpp"
#include "constitutive/capillaryPressure/Layouts.hpp"

#include "constitutive/ConstitutiveBase.hpp"
#include "functions/TableFunction.hpp"


namespace geos
{

namespace constitutive
{

/***
 * @brief KilloughHysteresis is designed to hold Killough hysteresis model parameters and
 *        be in charge of all computation related to this model (trapped saturation, Land coefficient)
 */

class KilloughHysteresis
{
public:

  /// To avoid division by zero, this is the min Scrd-Scri used in the computation of the Land constant
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
   * The struct is used to identify wetting or non wetting curves, and store key points as pairs of
   * saturation and value, where this value is either the relperm value (S,kr) or capillary pressure value (S,pc).
   * This way we can tell wetting curve from non wetting by the ordering of drainage/imbibition key values.
   * Indeed if critcal saturation for (Scri) imbibition comes at lower saturation than (Scrd) drainage then
   * it is wetting curve, on the opposite if critical saturation (Scrd) drainage comes before imbibition
   * this is a non-wetting hysteresis. This is completed by an extremum
   * point (Sextr) that is either the connate wetting saturation Swc or the maximum non wetting saturation Sgmax.
   *
   *     1 +-------------------------------+    1 +-------------------------------+
   *       |                            *##|      |                  ##        ***|
   *       |                          **## |      |                 ##       ***  |
   *       |                        **###  |      |               ##       ***    |
   *       |                      ***##    |      |              ##      ***      |
   *       |                   ****###     |      |            ###    ****        |
   *       |                **** ###       |      |          ###   ****           |
   *       |            *****  ###         |      |        ### *****              |
   *       |     + ****** +####            |      |  + ####*****  +               |
   *     0 +-------------------------------+    0 +-------------------------------+
   *                   Scrd     Scri             Sextr     Sextr        Scri            Scrd
   *        Fig. Left non wetting hysteresis for quadratic relperm / Right wetting hysteresis for quadratic relperm
   *
   * @param  m_extremumPhaseVolFraction represents either Swc or Sgmax depending if a wetting curve or nonwetting is described
   * @param m_criticalImbibitionPhaseVolFraction represents in wetting case the imibibition max and in non-wetting the imbibition residual
   * @param m_criticalDrainagePhaseVolFraction represents in wetting case the drainage max and in non-wetting the drainage residual
   * @param m_extremumValue represents the associate relperm or capillary pressure value
   * @param m_criticalImbibitionValue represents the associate relperm or capillary pressure value
   * @param m_criticalDrainageValue represents the associate relperm or capillary pressure value
   */
  struct HysteresisCurve
  {
    real64 m_extremumPhaseVolFraction = -1.;
    real64 m_criticalImbibitionPhaseVolFraction = -1.;
    real64 m_criticalDrainagePhaseVolFraction = -1.;

    real64 m_extremumValue = -1.;
    real64 m_criticalImbibitionValue = -1.;
    real64 m_criticalDrainageValue = -1.;

    bool m_isWetting = false;

    HysteresisCurve() = default;

    HysteresisCurve( real64 const & extremumPhaseVoFraction,
                     real64 const & extermumValue,
                     real64 const & criticalImbibitionPhaseVolFraction,
                     real64 const & criticalImbibitionValue,
                     real64 const & criticalDrainagePhaseVolFraction,
                     real64 const & criticalDrainageValue )
    {
      setPoints( extremumPhaseVoFraction,
                 extermumValue,
                 criticalImbibitionPhaseVolFraction,
                 criticalImbibitionValue,
                 criticalDrainagePhaseVolFraction,
                 criticalDrainageValue );

      m_isWetting =  ((m_criticalDrainagePhaseVolFraction > m_extremumPhaseVolFraction) ?
                      PhaseWettability::WETTING :
                      PhaseWettability::NONWETTING) == PhaseWettability::WETTING;

    }

    void setPoints( real64 const & extremumPhaseVoFraction,
                    real64 const & extermumValue,
                    real64 const & criticalImbibitionPhaseVolFraction,
                    real64 const & criticalImbibitionValue,
                    real64 const & criticalDrainagePhaseVolFraction,
                    real64 const & criticalDrainageValue )
    {
      m_extremumPhaseVolFraction = extremumPhaseVoFraction;
      m_criticalImbibitionPhaseVolFraction = criticalImbibitionPhaseVolFraction;
      m_criticalDrainagePhaseVolFraction = criticalDrainagePhaseVolFraction;

      m_extremumValue = extermumValue;
      m_criticalImbibitionValue = criticalImbibitionValue;
      m_criticalDrainageValue = criticalDrainageValue;
    }

    GEOS_HOST_DEVICE
    inline bool isWetting() const
    {
      return m_isWetting;
    }


    bool isInitialized() const
    {
      return (m_extremumPhaseVolFraction <= 0.0) && (m_criticalImbibitionPhaseVolFraction <= 0.0) && (m_criticalDrainagePhaseVolFraction <= 0.0);
    }

  };

  static std::string catalogName() { return "KilloughHysteresis"; }

  static void postProcessInput( real64 const & jerauldParam_a,
                                real64 const & jerauldParam_b,
                                real64 const & killoughCurvatureParamRelPerm );

  GEOS_HOST_DEVICE
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
  GEOS_HOST_DEVICE
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


/**
 * @brief Compute the Land coefficient for the wetting and non-wetting phases
 */
GEOS_HOST_DEVICE
inline void KilloughHysteresis::computeLandCoefficient( KilloughHysteresis::HysteresisCurve const & hcurve,
                                                        real64 & landParam )
{

  // Note: for simplicity, the notations are taken from IX documentation (although this breaks our phaseVolFrac naming convention)

  // Step 1: Land parameter for the wetting phase
  if( hcurve.isWetting() )
  {
    real64 const Scrd = hcurve.m_extremumPhaseVolFraction;
    real64 const Smxd = hcurve.m_criticalDrainagePhaseVolFraction;
    real64 const Smxi = hcurve.m_criticalImbibitionPhaseVolFraction;
    real64 const Swc = Scrd;
    GEOS_ERROR_IF(  (Smxi - Smxd) > 0,
                    GEOS_FMT( "{}: For wetting phase hysteresis, imbibition end-point saturation Smxi( {} ) must be smaller than the drainage saturation end-point Smxd( {} ).\n"
                              "Crossing relative permeability curves.\n",
                              catalogName(),
                              Smxi,
                              Smxd ));

    landParam = ( Smxd - Swc ) / LvArray::math::max( KilloughHysteresis::minScriMinusScrd, ( Smxd - Smxi ) ) - 1.0;
  }
  else
  // Step 2: Land parameter for the non-wetting phase

  {
    real64 const Smx =  hcurve.m_extremumPhaseVolFraction;
    real64 const Scrd = hcurve.m_criticalDrainagePhaseVolFraction;
    real64 const Scri = hcurve.m_criticalImbibitionPhaseVolFraction;
    GEOS_ERROR_IF( (Scrd - Scri) > 0,
                   GEOS_FMT( "{}: For non-wetting phase hysteresis, drainage trapped saturation Scrd( {} ) must be smaller than the imbibition saturation Scri( {} ).\n"
                             "Crossing relative permeability curves.\n",
                             catalogName(),
                             Scrd,
                             Scri ));

    landParam = ( Smx - Scrd ) / LvArray::math::max( KilloughHysteresis::minScriMinusScrd, ( Scri - Scrd ) ) - 1.0;
  }
}

/**
 * @brief Compute critical saturation (i.e. trapped) of the wetting and non-wetting phases
 */

GEOS_HOST_DEVICE
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
    real64 const Smxd = hcurve.m_criticalDrainagePhaseVolFraction;
    real64 const Swc = hcurve.m_extremumPhaseVolFraction;

    real64 const A = 1 + jerauldParam_a * (Shy - Swc);
    real64 const numerator = Shy - Smxd;
    real64 const denom = A + landParam * pow((Smxd - Shy) / (Smxd - Swc), 1 + jerauldParam_b / landParam );
    Scrt = Smxd + numerator / denom;
  }
  else
  {
    //unpack values
    real64 const Scrd = hcurve.m_criticalDrainagePhaseVolFraction;
    real64 const Smx = hcurve.m_extremumPhaseVolFraction;

    real64 const A = 1 + jerauldParam_a * (Smx - Shy);
    real64 const numerator = Shy - Scrd;
    real64 const denom = A + landParam * pow((Shy - Scrd) / (Smx - Scrd), 1 + jerauldParam_b / landParam );
    Scrt = LvArray::math::max( 0.0,
                               Scrd + numerator / denom );                       // trapped critical saturation from equation 2.162
  }

}

}

}

#endif //GEOS_CONSTITUTIVE_RELATIVEPERMEABILITY_KILLOUGHHYSTERESIS_HPP
