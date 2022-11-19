//
// Created by root on 11/9/22.
//

#include "KilloughHysteresisCapillaryPressure.hpp"

namespace geosx
{
namespace constitutive
{

KilloughHysteresisCapillaryPressure::KilloughHysteresisCapillaryPressureKernel::KilloughHysteresisCapillaryPressureKernel( const geosx::real64 & jerauldParam_a,
                                                                                                                           const geosx::real64 & jerauldParam_b,
                                                                                                                           const geosx::real64 & killoughCurvatureParam,
                                                                                                                           const arrayView1d< const geosx::real64 > & landParam,
                                                                                                                           const arrayView1d< const geosx::real64 > & drainageMinPhaseVolFraction,
                                                                                                                           const arrayView1d< const geosx::real64 > & imbibitionMinPhaseVolFraction,
                                                                                                                           const arrayView1d< const geosx::real64 > & drainageRelPermEndPoint,
                                                                                                                           const arrayView1d< const geosx::real64 > & imbibitionRelPermEndPoint,
                                                                                                                           const arrayView1d< const geosx::real64 > & drainageMaxPhaseVolFraction,
                                                                                                                           const arrayView1d< const geosx::real64 > & imbibitionMaxPhaseVolFraction,
                                                                                                                           const arrayView3d< geosx::real64,
                                                                                                                                              relperm::USD_RELPERM > & phaseTrapppedVolFrac )
  :
  KernelKilloughHysteresisBase( jerauldParam_a,
                                jerauldParam_b,
                                landParam,
                                drainageMinPhaseVolFraction,
                                imbibitionMinPhaseVolFraction,
                                drainageRelPermEndPoint,
                                imbibitionRelPermEndPoint,
                                drainageMaxPhaseVolFraction,
                                imbibitionMaxPhaseVolFraction,
                                phaseTrapppedVolFrac
                                ),
  m_killoughCurvatureParamCapPres( killoughCurvatureParam )
{ };

KilloughHysteresisCapillaryPressure::KilloughHysteresisCapillaryPressure( const geosx::string & name,
                                                                          geosx::dataRepository::Group * const parent )
  : KilloughHysteresis( name, parent )
{

  registerWrapper( viewKeyStruct::killoughCurvatureCapPresParameterString(), &m_killoughCurvatureParamCapPres ).
    setInputFlag(
    InputFlags::OPTIONAL ).
    setApplyDefaultValue(
    .1 ).
    setDescription(
    "Curvature parameter introduced by Killough for capillary pressure (see RTD documentation)." );
}


void KilloughHysteresisCapillaryPressure::postProcessInput()
{

  KilloughHysteresis::postProcessInput();

  GEOSX_THROW_IF( m_killoughCurvatureParamCapPres < 0,
                  GEOSX_FMT( "{}: the paramater {} must be postitive",
                             getFullName(),
                             viewKeyStruct::killoughCurvatureCapPresParameterString() ),
                  InputError );

}

}//end namespace constitutive
}//end namespace geosx
