//
// Created by root on 11/9/22.
//

#ifndef GEOSX_KILLOUGHHYSTERESISCAPILLARYPRESSURE_HPP
#define GEOSX_KILLOUGHHYSTERESISCAPILLARYPRESSURE_HPP

#include "constitutive/KilloughHysteresis.hpp"
#include "functions/FunctionManager.hpp"

namespace geosx
{
namespace constitutive
{


class KilloughHysteresisCapillaryPressure : public KilloughHysteresis
{

public:

  KilloughHysteresisCapillaryPressure( string const & name, Group * const parent );

  void postProcessInput() override;
//  KernelWrapper createKernelWrapper();
  static std::string catalogName() { return "KilloughHysteresisCapillaryPressure"; }

  virtual string getCatalogName() const override { return catalogName(); }

  class KilloughHysteresisCapillaryPressureKernel : public KernelKilloughHysteresisBase
  {

    KilloughHysteresisCapillaryPressureKernel( real64 const & jerauldParam_a,
                                               real64 const & jerauldParam_b,
                                               real64 const & killoughCurvatureParam,
                                               arrayView1d< real64 const > const & landParam,
                                               const arrayView1d< const geosx::real64 > & drainageMinPhaseVolFraction,
                                               const arrayView1d< const geosx::real64 > & imbibitionMinPhaseVolFraction,
                                               const arrayView1d< const geosx::real64 > & drainageRelPermEndPoint,
                                               const arrayView1d< const geosx::real64 > & imbibitionRelPermEndPoint,
                                               const arrayView1d< const geosx::real64 > & drainageMaxPhaseVolFraction,
                                               const arrayView1d< const geosx::real64 > & imbibitionMaxPhaseVolFraction,
                                               const arrayView3d< geosx::real64, relperm::USD_RELPERM > & phaseTrapppedVolFrac );


    GEOSX_HOST_DEVICE
    inline void
    computeImbibitionWettingCapPres( TableFunction::KernelWrapper const & drainageRelPermKernelWrapper,
                                     TableFunction::KernelWrapper const & imbibitionRelPermKernelWrapper,
//                                     real64 const & jerauldParam_a,
//                                     real64 const & jerauldParam_b,
                                     real64 const & landParam,
                                     real64 const & phaseVolFraction,
                                     real64 const & phaseMinHistoricalVolFraction,
                                     real64 const & imbibitionPhaseMinWettingVolFraction,
                                     real64 const & drainagePhaseMaxVolFraction,
                                     real64 const & imbibitionPhaseMaxVolFraction,
                                     real64 const & drainageRelPermEndPoint,
                                     real64 const & imbibitionRelPermEndPoint,
                                     real64 & phaseTrappedVolFrac,
                                     real64 & phaseRelPerm,
                                     real64 & dPhaseRelPerm_dPhaseVolFrac ) const;
    GEOSX_HOST_DEVICE
    inline void
    computeImbibitionNonWettingCapPres( TableFunction::KernelWrapper const & drainageRelPermKernelWrapper,
                                        TableFunction::KernelWrapper const & imbibitionRelPermKernelWrapper,
//                                        real64 const & jerauldParam_a,
//                                        real64 const & jerauldParam_b,
                                        real64 const & landParam,
                                        real64 const & phaseVolFraction,
                                        real64 const & phaseMaxHistoricalVolFraction,
                                        real64 const & drainagePhaseMinVolFraction,
                                        real64 const & imbibitionPhaseMinVolFraction,
                                        real64 const & drainagePhaseMaxVolFraction,
                                        real64 const & drainageRelPermEndPoint,
                                        real64 & phaseTrappedVolFrac,
                                        real64 & phaseRelPerm,
                                        real64 & dPhaseRelPerm_dPhaseVolFrac ) const;



private:
    /// Curvature parameter introduced for wetting phase hysteresis in Killough (for PC)
    real64 const m_killoughCurvatureParamCapPres;

  };

  struct viewKeyStruct : public KilloughHysteresis::viewKeyStruct
  {
    static constexpr char const * killoughCurvatureCapPresParameterString() { return "killoughCurvatureCapPresParameter"; }
  };

private:

  real64 m_killoughCurvatureParamCapPres;


};


GEOSX_HOST_DEVICE
inline void
KilloughHysteresisCapillaryPressure::KilloughHysteresisCapillaryPressureKernel::
  computeImbibitionWettingCapPres( TableFunction::KernelWrapper const & drainageRelPermKernelWrapper,
                                   TableFunction::KernelWrapper const & imbibitionRelPermKernelWrapper,
//                                     real64 const & jerauldParam_a,
//                                     real64 const & jerauldParam_b,
                                   real64 const & landParam,
                                   real64 const & phaseVolFraction,
                                   real64 const & phaseMinHistoricalVolFraction,
                                   real64 const & imbibitionPhaseMinWettingVolFraction,
                                   real64 const & drainagePhaseMaxVolFraction,
                                   real64 const & imbibitionPhaseMaxVolFraction,
                                   real64 const & drainageRelPermEndPoint,
                                   real64 const & imbibitionRelPermEndPoint,
                                   real64 & phaseTrappedVolFrac,
                                   real64 & phaseRelPerm,
                                   real64 & dPhaseRelPerm_dPhaseVolFrac ) const
{
  // use curvature param
}

GEOSX_HOST_DEVICE
inline void
KilloughHysteresisCapillaryPressure::KilloughHysteresisCapillaryPressureKernel::
  computeImbibitionNonWettingCapPres( TableFunction::KernelWrapper const & drainageRelPermKernelWrapper,
                                      TableFunction::KernelWrapper const & imbibitionRelPermKernelWrapper,
//                                        real64 const & jerauldParam_a,
//                                        real64 const & jerauldParam_b,
                                      real64 const & landParam,
                                      real64 const & phaseVolFraction,
                                      real64 const & phaseMaxHistoricalVolFraction,
                                      real64 const & drainagePhaseMinVolFraction,
                                      real64 const & imbibitionPhaseMinVolFraction,
                                      real64 const & drainagePhaseMaxVolFraction,
                                      real64 const & drainageRelPermEndPoint,
                                      real64 & phaseTrappedVolFrac,
                                      real64 & phaseRelPerm,
                                      real64 & dPhaseRelPerm_dPhaseVolFrac ) const
{
  //use curvature param
}


}
}

#endif //GEOSX_KILLOUGHHYSTERESISCAPILLARYPRESSURE_HPP
