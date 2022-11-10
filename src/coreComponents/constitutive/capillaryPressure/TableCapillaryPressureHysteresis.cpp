//
// Created by root on 11/8/22.
//

#include "TableCapillaryPressureHysteresis.hpp"

#include "constitutive/capillaryPressure/CapillaryPressureExtrinsicData.hpp"
#include "constitutive/capillaryPressure/TableCapillaryPressureHelpers.hpp"
#include "functions/FunctionManager.hpp"


namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

TableCapillaryPressureHysteresis::TableCapillaryPressureHysteresis( const std::string & name,
                                                                    dataRepository::Group * const parent )
  : CapillaryPressureBase( name, parent )
{

  registerWrapper( viewKeyStruct::drainagePhaseMinVolumeFractionString(), &m_drainagePhaseMinVolFraction ).
  setInputFlag( InputFlags::OPTIONAL).
  setDescription("");
  registerWrapper( viewKeyStruct::imbibitionPhaseMinVolumeFractionString(), &m_imbibitionPhaseMinVolFraction ).
  setInputFlag( InputFlags::OPTIONAL).
  setDescription("");
  //2phase
  registerWrapper( viewKeyStruct::drainageWettingNonWettingCapPresTableNameString(), &m_drainageWettingNonWettingCapPresTableName).
  setInputFlag( InputFlags::OPTIONAL ).
  setDescription("");
  registerWrapper( viewKeyStruct::imbibitionWettingNonWettingCapPresTableNameString(), &m_imbibitionWettingNonWettingCapPresTableName ).
  setInputFlag( InputFlags::OPTIONAL ).
                                        setDescription("");
  //3phase
  registerWrapper( viewKeyStruct::drainageWettingIntermediateCapPresTableNameString(),& m_drainageWettingIntermediateCapPresTableName ).
  setInputFlag( InputFlags::OPTIONAL ).
                                        setDescription("");
  registerWrapper( viewKeyStruct::drainageNonWettingIntermediateCapPresTableNameString(), &m_drainageNonWettingIntermediateCapPresTableName ).
  setInputFlag( InputFlags::OPTIONAL ).
                                        setDescription("");
  registerWrapper( viewKeyStruct::imbibitionWettingIntermediateCapPresTableNameString(), &m_imbibitionWettingIntermediateCapPresTableName ).
  setInputFlag( InputFlags::OPTIONAL ).
                                        setDescription("");
  registerWrapper( viewKeyStruct::imbibitionNonWettingIntermediateCapPresTableNameString(), &m_imbibitionNonWettingIntermediateCapPresTableName ).
  setInputFlag( InputFlags::OPTIONAL ).
                                        setDescription("");

  // ?
  registerWrapper( viewKeyStruct::drainageCapPresWrappersString(), & m_drainageCapillaryPressureKernelWrappers ).
  setInputFlag( InputFlags::OPTIONAL ).
                                        setDescription("");
  registerWrapper( viewKeyStruct::imbibitionCapPresWrappersString(), &m_imbibitionCapillaryPressureKernelWrappers ).
  setInputFlag( InputFlags::OPTIONAL ).
                                        setDescription("");

  // defaulted to 0.1
  registerWrapper( viewKeyStruct::killoughCurvatureParameterString(), &m_killoughCurvaturePc ).
  setInputFlag( InputFlags::OPTIONAL ).
  setDefaultValue(0.1).
                                        setDescription("");
  //used to compute different re-traversal path going drainage-imbibition-drainage
  registerWrapper( viewKeyStruct::tCurveOptionString(), &m_tCurveOption).
  setInputFlag( InputFlags::OPTIONAL ).
  setDefaultValue( 1 ).
                                        setDescription("");


};





}
} // geosx