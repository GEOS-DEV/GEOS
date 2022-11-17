//
// Created by root on 11/9/22.
//

#include "KilloughHysteresis.hpp"


namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

KilloughHysteresis::KernelKilloughHysteresisBase::KernelKilloughHysteresisBase( real64 const & jerauldParam_a,
                                                                                real64 const & jerauldParam_b,
                                                                                real64 const & killoughCurvatureParam,
                                                                                real64 const & killoughCurvaturePCParam,
                                                                                arrayView1d< real64 const > const & landParam,
                                                                                arrayView1d< const real64 > const & drainageMinPhaseVolFraction,
                                                                                arrayView1d< const real64 > const& imbibitionMinPhaseVolFraction,
                                                                                arrayView1d< const real64 > const & drainageRelPermEndPoint,
                                                                                arrayView1d< const real64 > const & imbibitionRelPermEndPoint,
                                                                                arrayView1d< const real64 > const & drainageMaxPhaseVolFraction,
                                                                                arrayView1d< const real64 > const & imbibitionMaxPhaseVolFraction,
                                                                                arrayView3d< real64, relperm::USD_RELPERM > const & phaseTrapppedVolFrac )
  :
  m_jerauldParam_a( jerauldParam_a ),
  m_jerauldParam_b( jerauldParam_b ),
  m_killoughCurvatureParam( killoughCurvatureParam ),
  m_killoughCurvatureParamCapPres( killoughCurvaturePCParam ),
  m_landParam( landParam ),
  m_drainagePhaseMinVolFraction( drainageMinPhaseVolFraction ),
  m_imbibitionPhaseMinVolFraction( imbibitionMinPhaseVolFraction ),
  m_drainagePhaseRelPermEndPoint( drainageRelPermEndPoint ),
  m_imbibitionPhaseRelPermEndPoint( imbibitionRelPermEndPoint ),
  m_drainagePhaseMaxVolFraction( drainageMaxPhaseVolFraction ),
  m_imbibitionPhaseMaxVolFraction( imbibitionMaxPhaseVolFraction )
{ }

KilloughHysteresis::KernelKilloughHysteresisBase::KernelKilloughHysteresisBase( const arrayView1d< const geosx::real64 > & landParam,
                                                                                real64 const & killoughCurvaturePCParam,
                                                                                const geosx::constitutive::KilloughHysteresis::HysteresisCurve_t & wettingCurve,
                                                                                const geosx::constitutive::KilloughHysteresis::HysteresisCurve_t & nonWettingCurve,
                                                                                const arrayView3d< geosx::real64, cappres::USD_CAPPRES > & phaseTrappedVolFrac ):
                                                                                m_killoughCurvatureParamCapPres(killoughCurvaturePCParam),
                                                                                m_landParam( landParam ),
                                                                                m_drainagePhaseMinVolFraction(KilloughHysteresis::toDrainagePhaseMinVolFraction(wettingCurve,nonWettingCurve)),
                                                                                m_imbibitionPhaseMinVolFraction(KilloughHysteresis::toDrainagePhaseMinVolFraction(wettingCurve,nonWettingCurve)),
                                                                                m_drainagePhaseMaxVolFraction(KilloughHysteresis::toDrainagePhaseMinVolFraction(wettingCurve,nonWettingCurve)),
                                                                                m_imbibitionPhaseMaxVolFraction(KilloughHysteresis::toDrainagePhaseMinVolFraction(wettingCurve,nonWettingCurve))
{}



real64 KilloughHysteresis::KernelKilloughHysteresisBase::getJerauldParamA() const
{
  return m_jerauldParam_a;
}

real64 KilloughHysteresis::KernelKilloughHysteresisBase::getJerauldParamB() const
{
  return m_jerauldParam_b;
}

real64 KilloughHysteresis::KernelKilloughHysteresisBase::getCurvatureParam() const
{
  return m_killoughCurvatureParam;
}

real64 KilloughHysteresis::KernelKilloughHysteresisBase::getCurvatureParamPc() const
{
  return m_killoughCurvatureParamCapPres;
}


KilloughHysteresis::KilloughHysteresis( const std::string & name,
                                        geosx::dataRepository::Group * const parent )
  :
  ConstitutiveBase( name, parent )
{

  registerWrapper( viewKeyStruct::jerauldParameterAString(), &m_jerauldParam_a ).
                                                                                  setInputFlag( InputFlags::OPTIONAL).
                                                                                  setApplyDefaultValue( 0.1 ).
                                                                                  setDescription(
    "First parameter (modification parameter) introduced by Jerauld in the Land trapping model (see RTD documentation)." );

  registerWrapper( viewKeyStruct::jerauldParameterBString(), &m_jerauldParam_b ).
                                                                                  setInputFlag( InputFlags::OPTIONAL ).
                                                                                  setApplyDefaultValue( 0.0 ).
                                                                                  setDescription(
    "Second parameter (modification parameter) introduced by Jerauld in the Land trapping model (see RTD documentation)." );

  registerWrapper( viewKeyStruct::killoughCurvatureParameterString(), &m_killoughCurvatureParam ).
                                                                                                   setInputFlag(
    InputFlags::OPTIONAL ).
                                                                                                   setApplyDefaultValue(
    1.0 ).
                                                                                                   setDescription(
    "Curvature parameter introduced by Killough for wetting-phase hysteresis (see RTD documentation)." );

  registerWrapper( viewKeyStruct::killoughCurvatureCapPresParameterString(), &m_killoughCurvatureParamCapPres ).
                                                                                                   setInputFlag(
    InputFlags::OPTIONAL ).
                                                                                                   setApplyDefaultValue(
    .1 ).
                                                                                                   setDescription(
    "Curvature parameter introduced by Killough for capillary pressure hysteresis (see RTD documentation)." );



}


void KilloughHysteresis::postProcessInput()
{
  GEOSX_THROW_IF( m_jerauldParam_a < 0,
                  GEOSX_FMT( "{}: the parameter {} must be positive",
                             getFullName(),
                             viewKeyStruct::jerauldParameterAString() ),
                  InputError );

  GEOSX_THROW_IF( m_jerauldParam_b < 0,
                  GEOSX_FMT( "{}: the paramater {} must be postitive",
                             getFullName(),
                             viewKeyStruct::jerauldParameterBString() ),
                  InputError );

  GEOSX_THROW_IF( m_killoughCurvatureParam < 0,
                  GEOSX_FMT( "{}: the paramater {} must be postitive",
                             getFullName(),
                             viewKeyStruct::killoughCurvatureParameterString() ),
                  InputError );


  GEOSX_THROW_IF( m_killoughCurvatureParamCapPres < 0,
                  GEOSX_FMT( "{}: the paramater {} must be postitive",
                             getFullName(),
                             viewKeyStruct::killoughCurvatureCapPresParameterString() ),
                  InputError );

}


KilloughHysteresis::KernelKilloughHysteresisBase
KilloughHysteresis::createKernelWrapper( arrayView1d< const geosx::real64 > const & landParam,
                                         arrayView1d< const geosx::real64 > const & drainageMinPhaseVolFraction,
                                         arrayView1d< const geosx::real64 > const & imbibitionMinPhaseVolFraction,
                                         arrayView1d< const geosx::real64 > const & drainageRelPermEndPoint,
                                         arrayView1d< const geosx::real64 > const & imbibitionRelPermEndPoint,
                                         arrayView1d< const geosx::real64 > const & drainageMaxPhaseVolFraction,
                                         arrayView1d< const geosx::real64 > const & imbibitionMaxPhaseVolFraction,
                                         arrayView3d< real64, relperm::USD_RELPERM > const & phaseTrapppedVolFrac )
{
  return KilloughHysteresis::KernelKilloughHysteresisBase( m_jerauldParam_a,
                                                           m_jerauldParam_b,
                                                           m_killoughCurvatureParam,
                                                           m_killoughCurvatureParamCapPres,
                                                           landParam,
                                                           drainageMinPhaseVolFraction,
                                                           imbibitionMinPhaseVolFraction,
                                                           drainageRelPermEndPoint,
                                                           imbibitionRelPermEndPoint,
                                                           drainageMaxPhaseVolFraction,
                                                           imbibitionMaxPhaseVolFraction,
                                                           phaseTrapppedVolFrac );
}


KilloughHysteresis::KernelKilloughHysteresisBase
KilloughHysteresis::createKernelWrapper( const arrayView1d< const geosx::real64 > & landParam,
                                         const geosx::constitutive::KilloughHysteresis::HysteresisCurve_t & wettingCurve,
                                         const geosx::constitutive::KilloughHysteresis::HysteresisCurve_t & nonWettingCurve,
                                         const arrayView3d< geosx::real64, cappres::USD_CAPPRES > & phaseTrappedVolFrac )
{
  return KilloughHysteresis::KernelKilloughHysteresisBase( landParam,
                                                           m_killoughCurvatureParamCapPres,
                                                           wettingCurve,
                                                           nonWettingCurve,
                                                           phaseTrappedVolFrac );
}


//void KilloughHysteresis::computeLandCoefficient()



REGISTER_CATALOG_ENTRY( ConstitutiveBase, KilloughHysteresis, std::string const &, Group * const )
}//end namespace
}//end namespace

