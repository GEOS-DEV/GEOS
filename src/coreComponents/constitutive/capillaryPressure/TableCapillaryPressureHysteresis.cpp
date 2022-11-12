//
// Created by root on 11/8/22.
//

#include "TableCapillaryPressureHysteresis.hpp"

#include "constitutive/capillaryPressure/CapillaryPressureExtrinsicData.hpp"
#include "constitutive/capillaryPressure/TableCapillaryPressureHelpers.hpp"
#include "functions/FunctionManager.hpp"
#include "constitutive/ConstitutiveManager.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

TableCapillaryPressureHysteresis::TableCapillaryPressureHysteresis( const std::string & name,
                                                                    dataRepository::Group * const parent )
  : CapillaryPressureBase( name, parent )
{

  registerWrapper( viewKeyStruct::wettingPhaseMinVolumeFractionString(), &m_wettingPhaseMinVolumeFraction ).
  setInputFlag(InputFlags::OPTIONAL).
  setDescription("");

  registerWrapper( viewKeyStruct::nonWettingPhaseMaxVolumeFractionString(), &m_nonWettingPhaseMaxVolumeFraction ).
                                         setInputFlag(InputFlags::OPTIONAL).
                                         setDescription("");

  registerWrapper( viewKeyStruct::drainageWettingPhaseMaxVolumeFractionString(), &m_drainageWettingPhaseMaxVolumeFraction ).
                                         setInputFlag(InputFlags::OPTIONAL).
                                         setDescription("");
  registerWrapper( viewKeyStruct::imbibitionWettingPhaseMaxVolumeFractionString(), &m_imbibitionWettingPhaseMaxVolumeFraction ).
                                         setInputFlag(InputFlags::OPTIONAL).
                                         setDescription("");
  registerWrapper( viewKeyStruct::drainageNonWettingPhaseMinVolumeFractionString(), &m_drainageNonWettingPhaseMinVolumeFraction ).
                                         setInputFlag(InputFlags::OPTIONAL).
                                         setDescription("");

  registerWrapper( viewKeyStruct::imbibitionNonWettingPhaseMinVolumeFractionString(), &m_imbibitionNonWettingPhaseMinVolumeFraction ).
                                         setInputFlag(InputFlags::OPTIONAL).
                                         setDescription("");
  //2phase
  registerWrapper( viewKeyStruct::drainageWettingNonWettingCapPresTableNameString(), &m_drainageWettingNonWettingCapPresTableName).
  setInputFlag(InputFlags::OPTIONAL).
  setDescription("");
  registerWrapper( viewKeyStruct::imbibitionWettingNonWettingCapPresTableNameString(), &m_imbibitionWettingNonWettingCapPresTableName ).
  setInputFlag(InputFlags::OPTIONAL).
                                        setDescription("");
  //3phase
  registerWrapper( viewKeyStruct::drainageWettingIntermediateCapPresTableNameString(),& m_drainageWettingIntermediateCapPresTableName ).
  setInputFlag(InputFlags::OPTIONAL).
                                        setDescription("");
  registerWrapper( viewKeyStruct::drainageNonWettingIntermediateCapPresTableNameString(), &m_drainageNonWettingIntermediateCapPresTableName ).
  setInputFlag(InputFlags::OPTIONAL).
                                        setDescription("");
  registerWrapper( viewKeyStruct::imbibitionWettingIntermediateCapPresTableNameString(), &m_imbibitionWettingIntermediateCapPresTableName ).
  setInputFlag(InputFlags::OPTIONAL).
                                        setDescription("");
  registerWrapper( viewKeyStruct::imbibitionNonWettingIntermediateCapPresTableNameString(), &m_imbibitionNonWettingIntermediateCapPresTableName ).
  setInputFlag(InputFlags::OPTIONAL).
                                        setDescription("");

  // ?
  registerWrapper( viewKeyStruct::drainageCapPresWrappersString(), & m_drainageWettingNonWettingCapillaryPressureKernelWrappers ).
  setInputFlag(InputFlags::OPTIONAL).
                                        setDescription("");
  registerWrapper( viewKeyStruct::imbibitionCapPresWrappersString(), &m_imbibitionWettingNonWettingCapillaryPressureKernelWrappers ).
  setInputFlag(InputFlags::OPTIONAL).
                                        setDescription("");

};

/// usual utils

void TableCapillaryPressureHysteresis::postProcessInput()
{
  CapillaryPressureBase::postProcessInput();

  using TPP = ThreePhasePairPhaseType;

  integer const numPhases = m_phaseNames.size();
  GEOSX_THROW_IF( numPhases != 2 && numPhases != 3,
                  GEOSX_FMT( "{}: the expected number of fluid phases is either two, or three",
                             getFullName() ),
                  InputError );

  m_phaseHasHysteresis.resize( 2 );

  if( numPhases == 2 )
  {
    GEOSX_THROW_IF( m_drainageWettingNonWettingCapPresTableName.empty(),
                   GEOSX_FMT( "{}: for a two-phase flow simulation, we must use {} to specify the capillary pressure table for the drainage pair (wetting phase, non-wetting phase)",
                               getFullName(),
                               viewKeyStruct::drainageWettingNonWettingCapPresTableNameString() ),
                    InputError );



    m_phaseHasHysteresis[TPP::INTERMEDIATE_WETTING] = ( m_imbibitionWettingNonWettingCapPresTableName.empty() ||
                                                        m_imbibitionWettingNonWettingCapPresTableName == m_drainageWettingNonWettingCapPresTableName )
                                                      ? 0 : 1;
    m_phaseHasHysteresis[TPP::INTERMEDIATE_NONWETTING] = m_phaseHasHysteresis[TPP::INTERMEDIATE_WETTING];
  }
  else if( numPhases == 3 )
  {



   GEOSX_THROW_IF( m_drainageWettingIntermediateCapPresTableName.empty() || m_drainageNonWettingIntermediateCapPresTableName.empty(),
    GEOSX_FMT( "{}: for a three-phase flow simulation, we must use {} to specify the capillary pressure table "
               "for the pair (wetting phase, intermediate phase), and {} to specify the capillary pressure table "
               "for the pair (non-wetting phase, intermediate phase)",
                               getFullName(),
                               viewKeyStruct::drainageWettingIntermediateCapPresTableNameString(),
                               viewKeyStruct::drainageNonWettingIntermediateCapPresTableNameString()  ),
                    InputError );

    m_phaseHasHysteresis[TPP::INTERMEDIATE_WETTING] = ( m_imbibitionWettingIntermediateCapPresTableName.empty() ||
                                           m_imbibitionWettingIntermediateCapPresTableName == m_drainageWettingIntermediateCapPresTableName )
                                         ? 0 : 1;

    m_phaseHasHysteresis[TPP::INTERMEDIATE_NONWETTING] = ( m_imbibitionNonWettingIntermediateCapPresTableName.empty() ||
                             m_imbibitionNonWettingIntermediateCapPresTableName == m_drainageNonWettingIntermediateCapPresTableName )
                           ? 0 : 1;
  }

  GEOSX_THROW_IF( m_phaseHasHysteresis[TPP::INTERMEDIATE_WETTING] == 0 && m_phaseHasHysteresis[TPP::INTERMEDIATE_NONWETTING] == 0,
                  GEOSX_FMT( "{}: we must use {} (2-phase) / {} or {} (3-phase) to specify at least one imbibition relative permeability table",
                             getFullName(),
                             viewKeyStruct::imbibitionWettingNonWettingCapPresTableNameString(),
                             viewKeyStruct::imbibitionWettingIntermediateCapPresTableNameString(),
                             viewKeyStruct::imbibitionNonWettingIntermediateCapPresTableNameString() ),
                  InputError );

}


void TableCapillaryPressureHysteresis::initializePreSubGroups()
{
  CapillaryPressureBase::initializePreSubGroups();

  integer const numPhases = m_phaseNames.size();
  FunctionManager const & functionManager = FunctionManager::getInstance();

  // Step 1: check sanity of drainage tables
  if( numPhases == 2 )
  {
    {
      GEOSX_THROW_IF( !functionManager.hasGroup( m_drainageWettingNonWettingCapPresTableName ),
                      GEOSX_FMT( "{}: the table function named {} could not be found",
                                 getFullName(),
                                 m_drainageWettingNonWettingCapPresTableName ),
                      InputError );
      TableFunction const
        & capPresTable = functionManager.getGroup< TableFunction >( m_drainageWettingNonWettingCapPresTableName );
      bool const capPresMustBeIncreasing = ( m_phaseOrder[PhaseType::WATER] < 0 )
                                           ? true   // pc on the gas phase, function must be increasing
                                           : false; // pc on the water phase, function must be decreasing
      TableCapillaryPressureHelpers::validateCapillaryPressureTable( capPresTable, getFullName(),
                                                                     capPresMustBeIncreasing );
    }
//should be an imibibition as we use this model
//define scope to avoid differentiate temp var (lazy)
    {
      GEOSX_THROW_IF( !functionManager.hasGroup( m_imbibitionWettingNonWettingCapPresTableName ),
                      GEOSX_FMT( "{}: the table function named {} could not be found",
                                 getFullName(),
                                 m_imbibitionWettingNonWettingCapPresTableName ),
                      InputError );
      TableFunction const
        & capPresTable = functionManager.getGroup< TableFunction >( m_imbibitionWettingNonWettingCapPresTableName );
      bool const capPresMustBeIncreasing = ( m_phaseOrder[PhaseType::WATER] < 0 )
                                           ? true   // pc on the gas phase, function must be increasing
                                           : false; // pc on the water phase, function must be decreasing
      TableCapillaryPressureHelpers::validateCapillaryPressureTable( capPresTable, getFullName(),
                                                                     capPresMustBeIncreasing );
    }

  }
  else if( numPhases == 3 )
  {
//define scope to avoid differentiate temp var (lazy)
    {
      GEOSX_THROW_IF( !functionManager.hasGroup( m_drainageWettingIntermediateCapPresTableName ),
                      GEOSX_FMT( "{}: the table function named {} could not be found",
                                 getFullName(),
                                 m_drainageWettingIntermediateCapPresTableName ),
                      InputError );
      TableFunction const
        & capPresTableWI = functionManager.getGroup< TableFunction >( m_drainageWettingIntermediateCapPresTableName );
      TableCapillaryPressureHelpers::validateCapillaryPressureTable( capPresTableWI, getFullName(), false );

      GEOSX_THROW_IF( !functionManager.hasGroup( m_drainageNonWettingIntermediateCapPresTableName ),
                      GEOSX_FMT( "{}: the table function named {} could not be found",
                                 getFullName(),
                                 m_drainageNonWettingIntermediateCapPresTableName ),
                      InputError );
      TableFunction const & capPresTableNWI =
        functionManager.getGroup< TableFunction >( m_drainageNonWettingIntermediateCapPresTableName );
      TableCapillaryPressureHelpers::validateCapillaryPressureTable( capPresTableNWI, getFullName(), true );
    }

    if( !m_imbibitionWettingIntermediateCapPresTableName.empty() )
    {

      GEOSX_THROW_IF( !functionManager.hasGroup( m_imbibitionWettingIntermediateCapPresTableName ),
                      GEOSX_FMT( "{}: the table function named {} could not be found",
                                 getFullName(),
                                 m_imbibitionWettingIntermediateCapPresTableName ),
                      InputError );
      TableFunction const
        & capPresTableWI = functionManager.getGroup< TableFunction >( m_imbibitionWettingIntermediateCapPresTableName );
      TableCapillaryPressureHelpers::validateCapillaryPressureTable( capPresTableWI, getFullName(), false );
    }

    if( !m_imbibitionNonWettingIntermediateCapPresTableName.empty() )
    {

      GEOSX_THROW_IF( !functionManager.hasGroup( m_imbibitionNonWettingIntermediateCapPresTableName ),
                      GEOSX_FMT( "{}: the table function named {} could not be found",
                                 getFullName(),
                                 m_imbibitionNonWettingIntermediateCapPresTableName ),
                      InputError );
      TableFunction const & capPresTableNWI =
        functionManager.getGroup< TableFunction >( m_imbibitionNonWettingIntermediateCapPresTableName );
      TableCapillaryPressureHelpers::validateCapillaryPressureTable( capPresTableNWI, getFullName(), true );
    }
  }

  // Step 2: check the sanity btw drainage and imbibition
  //TODO

  // Step 3: compute the Land coefficient
  computeLandCoefficient();
}

/// Land coeff (tb refactored out in KilloughHysteresis) and saved cvgd

void TableCapillaryPressureHysteresis::computeLandCoefficient()
{
  // For now, we keep two separate Land parameters for the wetting and non-wetting phases
  // For two-phase flow, we make sure that they are equal
  m_landParam.resize( 2 );

  integer ipWetting = 0;
  integer ipNonWetting = 0;

  integer const numPhases = m_phaseNames.size();
  if( numPhases == 2 )
  {
    ipWetting = ( m_phaseOrder[PhaseType::WATER] >= 0 ) ? m_phaseOrder[PhaseType::WATER] : m_phaseOrder[PhaseType::OIL];
    ipNonWetting = ( m_phaseOrder[PhaseType::GAS] >= 0 ) ? m_phaseOrder[PhaseType::GAS] : m_phaseOrder[PhaseType::OIL];
  }
  else if( numPhases == 3 )
  {
    ipWetting = m_phaseOrder[PhaseType::WATER];
    ipNonWetting = m_phaseOrder[PhaseType::GAS];
  }

  // Note: for simplicity, the notations are taken from IX documentation (although this breaks our phaseVolFrac naming convention)

  // Step 1: Land parameter for the wetting phase

  using TPP = ThreePhasePairPhaseType;

  {
    real64 const Scrd = m_wettingPhaseMinVolumeFraction;
    real64 const Smxd = m_drainageWettingPhaseMaxVolumeFraction;
    real64 const Smxi = m_imbibitionWettingPhaseMaxVolumeFraction;
    real64 const Swc = Scrd;
    GEOSX_THROW_IF(  (Smxi - Smxd) > 0,
                     GEOSX_FMT( "{}: For wetting phase hysteresis, imbibition end-point saturation Smxi( {} ) must be smaller than the drainage saturation end-point Smxd( {} ).\n"
                                "Crossing relative permeability curves.\n",
                                getFullName(),
                                Smxi,
                                Smxd ),
                     InputError );

    m_landParam[TPP::INTERMEDIATE_WETTING] = ( Smxd - Swc ) / LvArray::math::max( KilloughHysteresis::KernelKilloughHysteresisBase::minScriMinusScrd, ( Smxd - Smxi ) ) - 1.0;
  }

  // Step 2: Land parameter for the non-wetting phase

  {
    real64 const Scrd = m_drainageNonWettingPhaseMinVolumeFraction;
    real64 const Scri = m_imbibitionNonWettingPhaseMinVolumeFraction;
    real64 const Smx = m_nonWettingPhaseMaxVolumeFraction;
    GEOSX_THROW_IF( (Scrd - Scri) > 0,
                    GEOSX_FMT( "{}: For non-wetting phase hysteresis, drainage trapped saturation Scrd( {} ) must be smaller than the imbibition saturation Scri( {} ).\n"
                               "Crossing relative permeability curves.\n",
                               getFullName(),
                               Scrd,
                               Scri ),
                    InputError );

    m_landParam[TPP::INTERMEDIATE_NONWETTING] = ( Smx - Scrd ) / LvArray::math::max( KilloughHysteresis::KernelKilloughHysteresisBase::minScriMinusScrd, ( Scri - Scrd ) ) - 1.0;
  }
}


void TableCapillaryPressureHysteresis::saveConvergedPhaseVolFractionState( arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFraction ) const
{
  CapillaryPressureBase::saveConvergedState();

  arrayView2d< real64, compflow::USD_PHASE > phaseMaxHistoricalVolFraction = m_phaseMaxHistoricalVolFraction.toView();
  arrayView2d< real64, compflow::USD_PHASE > phaseMinHistoricalVolFraction = m_phaseMinHistoricalVolFraction.toView();

  localIndex const numElems = phaseVolFraction.size( 0 );
  integer const numPhases = numFluidPhases();

  forAll< parallelDevicePolicy<> >( numElems, [=] GEOSX_HOST_DEVICE ( localIndex const ei )
  {
    for( integer ip = 0; ip < numPhases; ++ip )
    {
      phaseMaxHistoricalVolFraction[ei][ip] = LvArray::math::max( phaseVolFraction[ei][ip], phaseMaxHistoricalVolFraction[ei][ip] );
      phaseMinHistoricalVolFraction[ei][ip] = LvArray::math::min( phaseVolFraction[ei][ip], phaseMinHistoricalVolFraction[ei][ip] );
    }
  } );

}

/// kernel creation
KilloughHysteresis::KernelKilloughHysteresisBase TableCapillaryPressureHysteresis::createKilloughKernelWrapper()
{
  ConstitutiveManager
    & constitutiveManager = this->getGroupByPath< ConstitutiveManager >( "/Problem/domain/Constitutive" );
  KilloughHysteresis & KilloughModel = constitutiveManager.getGroup< KilloughHysteresis >( m_KilloughModelName );
  //can use move semantic
  //TODO change Interface
  return KilloughModel.createKernelWrapper( m_landParam,
                                            m_drainagePhaseMinVolFraction,
                                            m_imbibitionPhaseMinVolFraction,
                                            m_drainagePhaseRelPermEndPoint,
                                            m_imbibitionPhaseRelPermEndPoint,
                                            m_drainagePhaseMaxVolFraction,
                                            m_imbibitionPhaseMaxVolFraction,
                                            m_phaseTrappedVolFrac );

}

TableCapillaryPressureHysteresis::KernelWrapper
TableCapillaryPressureHysteresis::createKernelWrapper()
{

  // we want to make sure that the wrappers are always up-to-date, so we recreate them everytime
  createAllTableKernelWrappers();

  m_KilloughKernel = createKilloughKernelWrapper();

  // then we create the actual TableRelativePermeabilityHysteresis::KernelWrapper
  return KernelWrapper( m_drainageRelPermKernelWrappers,
                        m_imbibitionRelPermKernelWrappers,
                        m_KilloughKernel,
                        m_phaseHasHysteresis,
                        m_landParam,
                        m_drainagePhaseMinVolFraction,
                        m_imbibitionPhaseMinVolFraction,
                        m_drainagePhaseMaxVolFraction,
                        m_imbibitionPhaseMaxVolFraction,
                        m_drainagePhaseRelPermEndPoint,
                        m_imbibitionPhaseRelPermEndPoint,
                        m_phaseTypes,
                        m_phaseOrder,
                        m_phaseMinHistoricalVolFraction,
                        m_phaseMaxHistoricalVolFraction,
                        m_phaseTrappedVolFrac,
                        m_phaseRelPerm,
                        m_dPhaseRelPerm_dPhaseVolFrac );
}

void TableCapillaryPressureHysteresis::createAllTableKernelWrappers()
{
  using TPP = ThreePhasePairPhaseType;

  FunctionManager const & functionManager = FunctionManager::getInstance();

  integer const numPhases = m_phaseNames.size();

  // we want to make sure that the wrappers are always up-to-date, so we recreate them everytime

  m_wettingNonWettingCapillaryPressureKernelWrappers.clear();
  m_wettingIntermediateCapillaryPressureKernelWrappers.clear();
  m_nonWettingIntermediateCapillaryPressureKernelWrappers.clear();

  if( numPhases == 2 )
  {

    TableFunction const & drainageCapPresTable = functionManager.getGroup< TableFunction >( m_drainageWettingNonWettingCapPresTableName );
    m_wettingNonWettingCapillaryPressureKernelWrappers.emplace_back( drainageCapPresTable.createKernelWrapper() );

    TableFunction const & imbibitionWettingCapPresTable = m_phaseHasHysteresis[TPP::INTERMEDIATE_WETTING]
                                                          ? functionManager.getGroup< TableFunction >( m_imbibitionWettingNonWettingCapPresTableName )
                                                          : functionManager.getGroup< TableFunction >( m_drainageWettingNonWettingCapPresTableName );
    m_wettingNonWettingCapillaryPressureKernelWrappers.emplace_back( imbibitionWettingCapPresTable.createKernelWrapper() );

  }
  else if( numPhases == 3 )
  {
      TableFunction const & drainageWICapPres = functionManager.getGroup< TableFunction >( m_drainageWettingIntermediateCapPresTableName );
      m_wettingIntermediateCapillaryPressureKernelWrappers.emplace_back( drainageWICapPres.createKernelWrapper() );

      TableFunction const & drainageNWICapPres = functionManager.getGroup< TableFunction >( m_drainageNonWettingIntermediateCapPresTableName );
      m_nonWettingIntermediateCapillaryPressureKernelWrappers.emplace_back( drainageNWICapPres.createKernelWrapper() );

    TableFunction const &  imbibitionWICapPres = m_phaseHasHysteresis[TPP::INTERMEDIATE_WETTING]
                                                          ? functionManager.getGroup< TableFunction >( m_imbibitionWettingIntermediateCapPresTableName )
                                                          : functionManager.getGroup< TableFunction >( m_drainageWettingIntermediateCapPresTableName );
    m_wettingIntermediateCapillaryPressureKernelWrappers.emplace_back( imbibitionWICapPres.createKernelWrapper() );

    TableFunction const & imbibitionNWICapPres = m_phaseHasHysteresis[TPP::INTERMEDIATE_NONWETTING]
                                                             ? functionManager.getGroup< TableFunction >( m_imbibitionNonWettingIntermediateCapPresTableName )
                                                             : functionManager.getGroup< TableFunction >( m_drainageNonWettingIntermediateCapPresTableName );
    m_nonWettingIntermediateCapillaryPressureKernelWrappers.emplace_back( imbibitionNWICapPres.createKernelWrapper() );
  }

}

/// resize Fields
void TableCapillaryPressureHysteresis::resizeFields( localIndex const size, localIndex const numPts )
{
  CapillaryPressureBase::resizeFields( size, numPts );

  integer const numPhases = numFluidPhases();

  m_phaseMaxHistoricalVolFraction.resize( size, numPhases );
  m_phaseMinHistoricalVolFraction.resize( size, numPhases );
  m_phaseMaxHistoricalVolFraction.setValues< parallelDevicePolicy<> >( 0.0 );
  m_phaseMinHistoricalVolFraction.setValues< parallelDevicePolicy<> >( 1.0 );
}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, TableCapillaryPressureHysteresis, std::string const &, Group * const )

}
} // geosx