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
                                                                                arrayView1d< real64 const > const & landParam,
                                                                                const arrayView1d< const geosx::real64 > & drainageMinPhaseVolFraction,
                                                                                const arrayView1d< const geosx::real64 > & imbibitionMinPhaseVolFraction,
                                                                                const arrayView1d< const geosx::real64 > & drainageRelPermEndPoint,
                                                                                const arrayView1d< const geosx::real64 > & imbibitionRelPermEndPoint,
                                                                                const arrayView1d< const geosx::real64 > & drainageMaxPhaseVolFraction,
                                                                                const arrayView1d< const geosx::real64 > & imbibitionMaxPhaseVolFraction,
                                                                                const arrayView3d< geosx::real64, relperm::USD_RELPERM > & phaseTrapppedVolFrac )
  :
  m_jerauldParam_a( jerauldParam_a ),
  m_jerauldParam_b( jerauldParam_b ),
  m_landParam( landParam ),
  m_drainagePhaseMinVolFraction( drainageMinPhaseVolFraction ),
  m_imbibitionPhaseMinVolFraction( imbibitionMinPhaseVolFraction ),
  m_drainagePhaseRelPermEndPoint( drainageRelPermEndPoint ),
  m_imbibitionPhaseRelPermEndPoint( imbibitionRelPermEndPoint ),
  m_drainagePhaseMaxVolFraction( drainageMaxPhaseVolFraction ),
  m_imbibitionPhaseMaxVolFraction( imbibitionMaxPhaseVolFraction )
{};


KilloughHysteresis::KilloughHysteresis( const std::string & name,
                                        geosx::dataRepository::Group * const parent )
  :
  ConstitutiveBase( name, parent )
{
  registerWrapper( viewKeyStruct::jerauldParameterAString(), &m_jerauldParam_a ).
                                                                                  setInputFlag( InputFlags::OPTIONAL ).
                                                                                  setApplyDefaultValue( 0.1 ).
                                                                                  setDescription(
    "First parameter (modification parameter) introduced by Jerauld in the Land trapping model (see RTD documentation)." );

  registerWrapper( viewKeyStruct::jerauldParameterBString(), &m_jerauldParam_b ).
                                                                                  setInputFlag( InputFlags::OPTIONAL ).
                                                                                  setApplyDefaultValue( 0.0 ).
                                                                                  setDescription(
    "Second parameter (modification parameter) introduced by Jerauld in the Land trapping model (see RTD documentation)." );

  // internal class data

  registerWrapper( viewKeyStruct::drainagePhaseMinVolumeFractionString(), &m_drainagePhaseMinVolFraction ).
                                                                                                            setInputFlag(
    InputFlags::FALSE ). // will be deduced from tables
                                                                                                            setSizedFromParent(
    0 );

  registerWrapper( viewKeyStruct::imbibitionPhaseMinVolumeFractionString(), &m_imbibitionPhaseMinVolFraction ).
                                                                                                                setInputFlag(
    InputFlags::FALSE ). // will be deduced from tables
                                                                                                                setSizedFromParent(
    0 );

  registerWrapper( viewKeyStruct::drainagePhaseRelPermEndPointString(), &m_drainagePhaseRelPermEndPoint ).
                                                                                                           setInputFlag(
    InputFlags::FALSE ). // will be deduced from tables
                                                                                                           setSizedFromParent(
    0 );

  registerWrapper( viewKeyStruct::imbibitionPhaseRelPermEndPointString(), &m_imbibitionPhaseRelPermEndPoint ).
                                                                                                               setInputFlag(
    InputFlags::FALSE ). // will be deduced from tables
                                                                                                               setSizedFromParent(
    0 );

  registerWrapper( viewKeyStruct::drainagePhaseMaxVolumeFractionString(), &m_drainagePhaseMaxVolFraction ).
                                                                                                            setInputFlag(
    InputFlags::FALSE ). // will be deduced from tables
                                                                                                            setSizedFromParent(
    0 );

  registerWrapper( viewKeyStruct::imbibitionPhaseMaxVolumeFractionString(), &m_imbibitionPhaseMaxVolFraction ).
                                                                                                                setInputFlag(
    InputFlags::FALSE ). // will be deduced from tables
                                                                                                                setSizedFromParent(
    0 );

  registerWrapper( viewKeyStruct::landParameterString(), &m_landParam ).
                                                                         setInputFlag( InputFlags::FALSE )
                                                                       . // will be deduced from tables
                                                                         setSizedFromParent( 0 );

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
}

//void KilloughHysteresis::computeLandCoefficient(real64 const & Scrd,
//                            real64 const & Shy,
//                            real64 const & Smx,
//                            real64 & landParam )
//{
//
//
//
//}



}//end namespace
}//end namespace

