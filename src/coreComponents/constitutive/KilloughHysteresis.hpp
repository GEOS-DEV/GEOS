//
// Created by root on 11/9/22.
//

#ifndef GEOSX_KILLOUGHHYSTERESIS_HPP
#define GEOSX_KILLOUGHHYSTERESIS_HPP

#include "constitutive/ConstitutiveBase.hpp"
//#include "common/GEOS_RAJA_Interface.hpp"
#include "functions/TableFunction.hpp"

#include "constitutive/relativePermeability/layouts.hpp"
//??
#include "constitutive/capillaryPressure/layouts.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

class KilloughHysteresis : public ConstitutiveBase
{
public:

  KilloughHysteresis(std::string const& name , Group * const  parent);

  virtual string getCatalogName() const = 0;

  virtual void postProcessInput();

//TODO  resize


//  /// by phase compute (inline ?)
//  GEOSX_HOST_DEVICE
//  void computeLandCoefficient(real64 const & Scrd,
//                              real64 const & Shy,
//                              real64 const & Smx,
//                              real64 & landParam );


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
                                   arrayView1d< real64 const > const & landParam,
                                   arrayView1d< real64 const > const & drainageMinPhaseVolFraction,
                                   arrayView1d< real64 const > const & imbibitionMinPhaseVolFraction,
                                   arrayView1d< real64 const > const & drainageRelPermEndPoint,
                                   arrayView1d< real64 const > const & imbibitionRelPermEndPoint,
                                   arrayView1d< real64 const > const & drainageMaxPhaseVolFraction,
                                   arrayView1d< real64 const > const & imbibitionMaxPhaseVolFraction,
                                   arrayView3d<real64, relperm::USD_RELPERM> const & phaseTrapppedVolFrac );




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
    void computeTrappedCriticalPhaseVolFraction( real64 const & Scrd,
                                                 real64 const & Shy,
                                                 real64 const & Smx,
                                                 real64 const & landParam,
                                                 real64 & Scrt ) const;



    /// Parameter a introduced by Jerauld in the Land model
    real64 const m_jerauldParam_a;

    /// Parameter b introduced by Jerauld in the Land model
    real64 const m_jerauldParam_b;

    /// Trapping parameter from the Land model (typically called C)
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

  struct viewKeyStruct
  {
    static constexpr char const * jerauldParameterAString() { return "jerauldParameterA"; }
    static constexpr char const * jerauldParameterBString() { return "jerauldParameterB"; }
    static constexpr char const * landParameterString() { return "landParameter"; }

    static constexpr char const * drainagePhaseRelPermEndPointString() { return "drainagePhaseRelPermEndPoint"; }
    static constexpr char const * imbibitionPhaseRelPermEndPointString() { return "imbibitionPhaseRelPermEndPoint"; }
    static constexpr char const * drainagePhaseMinVolumeFractionString() { return "drainagePhaseMinVolumeFraction"; }
    static constexpr char const * imbibitionPhaseMinVolumeFractionString() { return "imbibitionPhaseMinVolumeFraction"; }
    static constexpr char const * drainagePhaseMaxVolumeFractionString() { return "drainagePhaseMaxVolumeFraction"; }
    static constexpr char const * imbibitionPhaseMaxVolumeFractionString() { return "imbibitionPhaseMaxVolumeFraction"; }
  };

private:

  /// Parameter a introduced by Jerauld in the Land model
  real64 m_jerauldParam_a;

  /// Parameter b introduced by Jerauld in the Land model
  real64 m_jerauldParam_b;

  /// Trapping parameter from the Land model (typically called C)
  array1d< real64 > m_landParam;

  /// Minimum volume fraction for each phase in drainage (deduced from the drainage table)
  array1d< real64 > m_drainagePhaseMinVolFraction;

  /// Minimum volume fraction for each phase in imbibition (deduced from the imbibition table)
  array1d< real64 > m_imbibitionPhaseMinVolFraction;

  /// Relperm endpoint for each phase in drainage (deduced from the drainage table)
  array1d< real64 > m_drainagePhaseRelPermEndPoint;

  /// Relperm endpoint for each phase in imbibition (deduced from the imbibition table)
  array1d< real64 > m_imbibitionPhaseRelPermEndPoint;

  /// Maximum volume fraction for each phase
  array1d< real64 > m_drainagePhaseMaxVolFraction;

  /// Maximum volume fraction for each phase
  array1d< real64 > m_imbibitionPhaseMaxVolFraction;


};

//inlines

GEOSX_HOST_DEVICE
inline void
KilloughHysteresis::KernelKilloughHysteresisBase::
computeTrappedCriticalPhaseVolFraction( real64 const & Scrd,
                                        real64 const & Shy,
                                        real64 const & Smx,
                                        real64 const & landParam,
                                        real64 & Scrt ) const
{
  real64 const A = 1 + m_jerauldParam_a * ( Smx - Shy );
  real64 const numerator = Shy - Scrd;
  real64 const denom = A + landParam * pow( ( Shy - Scrd ) / ( Smx - Scrd ), 1 + m_jerauldParam_b / landParam );
  Scrt = LvArray::math::max( 0.0, Scrd + numerator / denom ); // trapped critical saturation from equation 2.162
}




}

}

#endif //GEOSX_KILLOUGHHYSTERESIS_HPP
