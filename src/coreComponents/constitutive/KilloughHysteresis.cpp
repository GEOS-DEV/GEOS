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
  m_landParam( landParam ),
  m_drainagePhaseMinVolFraction( drainageMinPhaseVolFraction ),
  m_imbibitionPhaseMinVolFraction( imbibitionMinPhaseVolFraction ),
  m_drainagePhaseRelPermEndPoint( drainageRelPermEndPoint ),
  m_imbibitionPhaseRelPermEndPoint( imbibitionRelPermEndPoint ),
  m_drainagePhaseMaxVolFraction( drainageMaxPhaseVolFraction ),
  m_imbibitionPhaseMaxVolFraction( imbibitionMaxPhaseVolFraction )
{ }

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
};


KilloughHysteresis::KilloughHysteresis( const std::string & name,
                                        geosx::dataRepository::Group * const parent )
  :
  RelativePermeabilityBase( name, parent )
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
                                                           landParam,
                                                           drainageMinPhaseVolFraction,
                                                           imbibitionMinPhaseVolFraction,
                                                           drainageRelPermEndPoint,
                                                           imbibitionRelPermEndPoint,
                                                           drainageMaxPhaseVolFraction,
                                                           imbibitionMaxPhaseVolFraction,
                                                           phaseTrapppedVolFrac );
}

//void KilloughHysteresis::computeLandCoefficient()



REGISTER_CATALOG_ENTRY( ConstitutiveBase, KilloughHysteresis, std::string const &, Group * const )
}//end namespace
}//end namespace

