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

#include "relativePermeability/layouts.hpp"
#include "capillaryPressure/layouts.hpp"



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


  //move to private ?
  struct HysteresisCurve_t
  {
    real64 oppositeBound;
    real64 imbibitionExtrema;
    real64 drainageExtrema;
    HysteresisCurve_t( real64 const & opp, real64 const & imbE, real64 const & drainE )
    {
      oppositeBound = opp; imbibitionExtrema = imbE; drainageExtrema = drainE;
    }
    bool isWetting() const
    {
      return ((drainageExtrema > oppositeBound) ?  PhaseWettability::WETTING : PhaseWettability::NONWETTING) == PhaseWettability::WETTING;
    }
    bool isZero() const
    {
      return (oppositeBound<=0.0) && (imbibitionExtrema<=0.0) && (drainageExtrema<=0.0);
    }



  };

  static array1d< real64 > toDrainagePhaseMinVolFraction( HysteresisCurve_t const & c1, HysteresisCurve_t const & c2 )
  {
    array1d< real64 > ret;
    ret.resize( 2 );

    if( c1.isWetting() && !c2.isWetting())
    {
      ret[0] = c1.oppositeBound;
      ret[1] = c2.drainageExtrema;
    }
    else if( !c1.isWetting() && c2.isWetting() )
    {
      return toDrainagePhaseMinVolFraction( c2, c1 );
    }
    else
    {
      GEOSX_THROW( "not the good type", InputError );
    }

    return ret;
  }

  static array1d< real64 > toDrainagePhaseMaxVolFraction( HysteresisCurve_t const & c1, HysteresisCurve_t const & c2 )
  {
    array1d< real64 > ret;
    ret.resize( 2 );

    if( c1.isWetting() && !c2.isWetting())
    {
      ret[0] = c1.drainageExtrema;
      ret[1] = c2.oppositeBound;
    }
    else if( !c1.isWetting() && c2.isWetting() )
    {
      return toDrainagePhaseMinVolFraction( c2, c1 );
    }
    else
    {
      GEOSX_THROW( "not the good type", InputError );
    }

    return ret;
  }

  static array1d< real64 > toImbibitionPhaseMinVolFraction( HysteresisCurve_t const & c1, HysteresisCurve_t const & c2 )
  {
    array1d< real64 > ret;
    ret.resize( 2 );

    if( c1.isWetting() && !c2.isWetting())
    {
      ret[0] = c1.oppositeBound;
      ret[1] = c2.imbibitionExtrema;
    }
    else if( !c1.isWetting() && c2.isWetting() )
    {
      return toDrainagePhaseMinVolFraction( c2, c1 );
    }
    else
    {
      GEOSX_THROW( "not the good type", InputError );
    }

    return ret;
  }

  static array1d< real64 > toImbibitionPhaseMaxVolFraction( HysteresisCurve_t const & c1, HysteresisCurve_t const & c2 )
  {
    array1d< real64 > ret;
    ret.resize( 2 );

    if( c1.isWetting() && !c2.isWetting())
    {
      ret[0] = c1.imbibitionExtrema;
      ret[1] = c2.oppositeBound;
    }
    else if( !c1.isWetting() && c2.isWetting() )
    {
      return toDrainagePhaseMinVolFraction( c2, c1 );
    }
    else
    {
      GEOSX_THROW( "not the good type", InputError );
    }

    return ret;
  }

  KilloughHysteresis( std::string const & name, Group * const parent );

  static std::string catalogName() { return "KilloughHysteresis"; }

  virtual string getCatalogName() const override { return catalogName(); }

  virtual void postProcessInput() override;

//  virtual void resizeFields( localIndex const size, localIndex const numPts ) override;

//  GEOSX_HOST_DEVICE
//  void computeLandCoefficient();


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
    //TODO refactor end point and min as a struct _RelpermData_
    KernelKilloughHysteresisBase( real64 const & jerauldParam_a,
                                  real64 const & jerauldParam_b,
                                  real64 const & killoughCruvParam,
                                  real64 const & killoughCruvPcParam,
                                  arrayView1d< real64 const > const & landParam,
                                  arrayView1d< real64 const > const & drainageMinPhaseVolFraction,
                                  arrayView1d< real64 const > const & imbibitionMinPhaseVolFraction,
                                  arrayView1d< real64 const > const & drainageRelPermEndPoint,
                                  arrayView1d< real64 const > const & imbibitionRelPermEndPoint,
                                  arrayView1d< real64 const > const & drainageMaxPhaseVolFraction,
                                  arrayView1d< real64 const > const & imbibitionMaxPhaseVolFraction,
                                  arrayView3d< real64, relperm::USD_RELPERM > const & phaseTrapppedVolFrac );

    KernelKilloughHysteresisBase( arrayView1d< real64 const > const & landParam,
                                  real64 const & killoughCurvaturePCParam,
                                  HysteresisCurve_t const & wettingCurve,
                                  HysteresisCurve_t const & nonWettingCurve,
                                  arrayView3d< real64, cappres::USD_CAPPRES > const & phaseTrappedVolFrac );


    KernelKilloughHysteresisBase() = default;

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
    real64 getCurvatureParamPc() const;

private:
//from overnested
    /// Parameter a introduced by Jerauld in the Land model
    real64 m_jerauldParam_a;
    /// Parameter b introduced by Jerauld in the Land model
    real64 m_jerauldParam_b;

    real64 m_killoughCurvatureParam;

    real64 m_killoughCurvatureParamCapPres;

    //from main Relperm Class
    // Trapping parameter from the Land model (typically called C)
    arrayView1d< real64 const > m_landParam;
    /// Minimum volume fraction for each phase in drainage (deduced from the drainage table)
    arrayView1d< real64 const > m_drainagePhaseMinVolFraction;

    /// Minimum volume fraction for each phase in imbibition (deduced from the imbibition table)
    arrayView1d< real64 const > m_imbibitionPhaseMinVolFraction;

    /// Relperm endpoint for each phase in drainage (deduced from the drainage table)
    arrayView1d< real64 const > m_drainagePhaseRelPermEndPoint;

    /// Relperm endpoint for each phase in imbibition (deduced from the imbibition table)
    arrayView1d< real64 const > m_imbibitionPhaseRelPermEndPoint;

    /// Maximum volume fraction for each phase
    arrayView1d< real64 const > m_drainagePhaseMaxVolFraction;

    /// Maximum volume fraction for each phase
    arrayView1d< real64 const > m_imbibitionPhaseMaxVolFraction;


  };

  KernelKilloughHysteresisBase createKernelWrapper( arrayView1d< real64 const > const & landParam,
                                                    arrayView1d< real64 const > const & drainageMinPhaseVolFraction,
                                                    arrayView1d< real64 const > const & imbibitionMinPhaseVolFraction,
                                                    arrayView1d< real64 const > const & drainageRelPermEndPoint,
                                                    arrayView1d< real64 const > const & imbibitionRelPermEndPoint,
                                                    arrayView1d< real64 const > const & drainageMaxPhaseVolFraction,
                                                    arrayView1d< real64 const > const & imbibitionMaxPhaseVolFraction,
                                                    arrayView3d< real64, relperm::USD_RELPERM > const & phaseTrapppedVolFrac ) const;

  KernelKilloughHysteresisBase createKernelWrapper( arrayView1d< real64 const > const & landParam,
                                                    HysteresisCurve_t const & wettingCurve,
                                                    HysteresisCurve_t const & nonWettingCurve,
                                                    arrayView3d< real64, cappres::USD_CAPPRES > const & phaseTrappedVolFrac ) const;



  struct viewKeyStruct
  {
    static constexpr char const * jerauldParameterAString() { return "jerauldParameterA"; }
    static constexpr char const * jerauldParameterBString() { return "jerauldParameterB"; }

    static constexpr char const * killoughCurvatureParameterString() { return "killoughCurvatureParameter"; }
    static constexpr char const * killoughCurvatureCapPresParameterString() { return "killoughCurvatureCapPresParameter"; }
  };

private:

  /// Parameter a introduced by Jerauld in the Land model
  real64 m_jerauldParam_a;

  /// Parameter b introduced by Jerauld in the Land model
  real64 m_jerauldParam_b;

  /// Curvature parameter introduced for wetting phase hysteresis in Killough
  real64 m_killoughCurvatureParam;
  real64 m_killoughCurvatureParamCapPres;


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
    real64 const Smxd = hcurve.drainageExtrema;
    real64 const Swc = hcurve.oppositeBound;

    real64 const A = 1 + m_jerauldParam_a * ( Shy - Swc );
    real64 const numerator = Shy - Smxd;
    real64 const denom = A + landParam * pow( ( Smxd - Shy ) / ( Smxd - Swc ), 1 + m_jerauldParam_b/landParam );
    Scrt = Smxd + numerator / denom;
  }
  else
  {
    //unpack values
    real64 const Scrd = hcurve.drainageExtrema;
    real64 const Smx = hcurve.oppositeBound;

    real64 const A = 1 + m_jerauldParam_a * ( Smx - Shy );
    real64 const numerator = Shy - Scrd;
    real64 const denom = A + landParam * pow( ( Shy - Scrd ) / ( Smx - Scrd ), 1 + m_jerauldParam_b / landParam );
    Scrt = LvArray::math::max( 0.0, Scrd + numerator / denom ); // trapped critical saturation from equation 2.162
  }

}


}

}

#endif //GEOSX_KILLOUGHHYSTERESIS_HPP
