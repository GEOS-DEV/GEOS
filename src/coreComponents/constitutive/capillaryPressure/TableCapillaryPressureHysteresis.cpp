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

  //2p
  registerWrapper( viewKeyStruct::wettingNonWettingCapillaryPressureKernelWrappersString(), &m_wettingNonWettingCapillaryPressureKernelWrappers).
  setInputFlag(InputFlags::OPTIONAL).setDescription("");
  //3p
  registerWrapper( viewKeyStruct::wettingIntermediateCapillaryPressureKernelWrappersString(), &m_wettingIntermediateCapillaryPressureKernelWrappers).
  setInputFlag(InputFlags::OPTIONAL).setDescription("");
  registerWrapper( viewKeyStruct::nonWettingIntermediateCapillaryPressureKernelWrappersString(), &m_nonWettingIntermediateCapillaryPressureKernelWrappers).
  setInputFlag(InputFlags::OPTIONAL).setDescription("");

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

void
TableCapillaryPressureHysteresis::KernelWrapper::computeImbibitionWettingCapillaryPressure( arrayView1d< TableFunction::KernelWrapper const > const & wettingKernelWapper,
                                                                                            const KilloughHysteresis::HysteresisCurve_t & wettingCurve,
                                                                                            const geosx::real64 & landParam,
                                                                                            const geosx::real64 & phaseVolFraction,
                                                                                            const geosx::real64 & phaseMinHistoricalVolFraction,
                                                                                            const real64 & phaseIntermediateMinVolFraction,
                                                                                            geosx::real64 & phaseTrappedVolFrac,
                                                                                            geosx::real64 & phaseCapPressure,
                                                                                            geosx::real64 & dPhaseCapPressure_dPhaseVolFrac ) const
{
  GEOSX_ASSERT(wettingCurve.isWetting());
  real64 const S = phaseVolFraction;
  real64 const Smxi = wettingCurve.imbibitionExtrema;
  real64 const Smxd = wettingCurve.drainageExtrema;
  real64 const Smin = wettingCurve.oppositeBound;

  if( S<=Smin)
  {
    //below accessible range
    phaseCapPressure = CAP_INF;
    dPhaseCapPressure_dPhaseVolFrac = - CAP_INF_DERIV;
  }
  else if( S >= Smxd)
  {
    //above accessible range
    phaseCapPressure = - CAP_INF;
    dPhaseCapPressure_dPhaseVolFrac = - CAP_INF_DERIV;
  }
  else
  {
    //drainage to imbibition
    real64 dpci_dS, dpcd_dS;
    real64 const pci = wettingKernelWapper[ModeIndexType::IMBIBITION].compute(&S, &dpci_dS);
    real64 const pcd = wettingKernelWapper[ModeIndexType::DRAINAGE].compute(&S,&dpcd_dS);
    real64 const Somin = phaseIntermediateMinVolFraction;

    // Step 1: get the trapped from wetting data
    real64 const Shy = ( phaseMinHistoricalVolFraction > Smin ) ? phaseMinHistoricalVolFraction : Smin;

    real64 Scrt = 0.0;
    m_KilloughKernel.computeTrappedCriticalPhaseVolFraction( wettingCurve, Shy, landParam,
                                                                                 Scrt );
    real64 const Swma = 1 - Scrt - Somin;

    real64 const E = m_KilloughKernel.getCurvatureParamPc();

    //Step 2. compute F as in (EQ 34.15) F = (1/(Sw-Shy+E)-1/E) / (1/(Swma-Shy+E)-1/E)
    real64 const F = (1./(S-Shy+E) - 1./E) / (1./(Swma-Shy+E) - 1./E);

    //Step 3. compute dF_dS
    real64 dF_dS = (-1./(S*S)) / (1./(Swma-Shy+E) - 1./E);

    //Step 4. Eventually assemble everything following (EQ. 34.14)
    phaseCapPressure = pcd + F * ( pci - pcd );
    dPhaseCapPressure_dPhaseVolFrac = dpci_dS + F * ( dpci_dS - dpcd_dS);
    dPhaseCapPressure_dPhaseVolFrac += dF_dS * ( pci -pcd );

    //TODO reverse flow

  }

}

void TableCapillaryPressureHysteresis::KernelWrapper::computeTwoPhase( const geosx::integer ipWetting,
                                                                       const geosx::integer ipNonWetting,
                                                                       const arraySlice1d< const geosx::real64,
                                                                         compflow::USD_PHASE - 1 > & phaseVolFraction,
                                                                       const arraySlice1d< const geosx::real64,
                                                                         compflow::USD_PHASE
                                                                         - 1 > & phaseMaxHistoricalVolFraction,
                                                                       const arraySlice1d< const geosx::real64,
                                                                         compflow::USD_PHASE
                                                                         - 1 > & phaseMinHistoricalVolFraction,
                                                                       const arraySlice1d< geosx::real64,
                                                                         relperm::USD_RELPERM
                                                                         - 2 > & phaseTrappedVolFrac,
                                                                       arraySlice1d< geosx::real64,
                                                                         relperm::USD_RELPERM
                                                                         - 2 > const & phaseCapPressure,
                                                                       arraySlice2d< geosx::real64,
                                                                         relperm::USD_RELPERM_DS
                                                                         - 2 > const & dPhaseCapPressure_dPhaseVolFrac ) const
{
  using TTP = ThreePhasePairPhaseType;

  //--- wetting  cap pressure
  if( !m_phaseHasHysteresis[TTP::INTERMEDIATE_WETTING] || phaseVolFraction[ipWetting] <= phaseMinHistoricalVolFraction[ipWetting] + flowReversalBuffer )
  {
    phaseTrappedVolFrac[ipWetting] = LvArray::math::min( phaseVolFraction[ipWetting], m_wettingCurve.oppositeBound );
    computeDrainageCapillaryPressure( m_wettingNonWettingCapillaryPressureKernelWrappers[ModeIndexType::DRAINAGE],
                                      phaseVolFraction[ipWetting],
                                      phaseCapPressure[ipWetting],
                                      dPhaseCapPressure_dPhaseVolFrac[ipWetting][ipWetting] );

  }
  else
  {
      computeImbibitionWettingCapillaryPressure( m_wettingNonWettingCapillaryPressureKernelWrappers,
                                                 m_wettingCurve,
                                                 m_landParam[ipWetting],
                                                 phaseVolFraction[ipWetting],
                                                 phaseMinHistoricalVolFraction[ipWetting],
                                                 m_phaseIntermediateMinVolFraction,
                                                 phaseTrappedVolFrac[ipWetting],
                                                 phaseCapPressure[ipWetting],
                                                 dPhaseCapPressure_dPhaseVolFrac[ipWetting][ipWetting]
                                                 );

    }

  //-- non wetting

  if( !m_phaseHasHysteresis[TTP::INTERMEDIATE_NONWETTING] || phaseVolFraction[ipNonWetting] >= phaseMaxHistoricalVolFraction[ipNonWetting] + flowReversalBuffer )
    {
      phaseTrappedVolFrac[ipNonWetting] = LvArray::math::min( phaseVolFraction[ipNonWetting], m_nonWettingCurve.drainageExtrema );
      computeDrainageCapillaryPressure( m_wettingNonWettingCapillaryPressureKernelWrappers[ModeIndexType::DRAINAGE],
                                        phaseVolFraction[ipNonWetting],
                                        phaseCapPressure[ipNonWetting],
                                        dPhaseCapPressure_dPhaseVolFrac[ipNonWetting][ipNonWetting] );
      // when pc is on the gas phase, we need to multiply user input by -1
      // because CompositionalMultiphaseFVM does: pres_gas = pres_oil - pc_og, so we need a negative pc_og
      phaseCapPressure[ipNonWetting] *=-1;
      dPhaseCapPressure_dPhaseVolFrac[ipNonWetting][ipNonWetting] *= -1;

    }

  {
      computeImbibitionNonWettingCapillaryPressure( m_wettingNonWettingCapillaryPressureKernelWrappers,
                                                    m_nonWettingCurve,
                                                    m_landParam[ipNonWetting],
                                                    phaseVolFraction[ipNonWetting],
                                                    phaseMaxHistoricalVolFraction[ipNonWetting],
                                                    phaseTrappedVolFrac[ipNonWetting],
                                                    phaseCapPressure[ipNonWetting],
                                                    dPhaseCapPressure_dPhaseVolFrac[ipNonWetting][ipNonWetting]
                                                    );

    // when pc is on the gas phase, we need to multiply user input by -1
    // because CompositionalMultiphaseFVM does: pres_gas = pres_oil - pc_og, so we need a negative pc_og
    phaseCapPressure[ipNonWetting] *= -1;
    dPhaseCapPressure_dPhaseVolFrac[ipNonWetting][ipNonWetting] *= -1;
    }

}


void TableCapillaryPressureHysteresis::KernelWrapper::computeThreePhase( const geosx::integer ipWetting,
                                                                         const geosx::integer ipInter,
                                                                         const geosx::integer ipNonWetting,
                                                                         const arraySlice1d< const geosx::real64,
                                                                           compflow::USD_PHASE - 1 > & phaseVolFraction,
                                                                         const arraySlice1d< const geosx::real64,
                                                                           compflow::USD_PHASE
                                                                           - 1 > & phaseMaxHistoricalVolFraction,
                                                                         const arraySlice1d< const geosx::real64,
                                                                           compflow::USD_PHASE
                                                                           - 1 > & phaseMinHistoricalVolFraction,
                                                                         const arraySlice1d< geosx::real64,
                                                                           relperm::USD_RELPERM
                                                                           - 2 > & phaseTrappedVolFrac,
                                                                         const arraySlice1d< geosx::real64,
                                                                           relperm::USD_RELPERM - 2 > & phaseCapPressure,
                                                                         const arraySlice2d< geosx::real64,
                                                                           relperm::USD_RELPERM_DS
                                                                           - 2 > & dPhaseCapPressure_dPhaseVolFrac ) const
{


  LvArray::forValuesInSlice( dPhaseCapPressure_dPhaseVolFrac, []( real64 & val ){ val = 0.0; } );
  using TTP = ThreePhasePairPhaseType;

  // -- wetting curve if drainage only
     if( !m_phaseHasHysteresis[TTP::INTERMEDIATE_WETTING] || phaseVolFraction[ipWetting] <= phaseMinHistoricalVolFraction[ipWetting] + flowReversalBuffer )
  {
    //TODO trapped
       // water-oil capillary pressure
    phaseTrappedVolFrac[ipWetting] = LvArray::math::min( phaseVolFraction[ipWetting], m_wettingCurve.oppositeBound );
    phaseCapPressure[ipWetting] =
    m_wettingIntermediateCapillaryPressureKernelWrappers[ModeIndexType::DRAINAGE].compute( &(phaseVolFraction)[ipWetting],
                                                                &(dPhaseCapPressure_dPhaseVolFrac)[ipWetting][ipWetting] );
  }
     else
     {
       computeImbibitionWettingCapillaryPressure( m_wettingIntermediateCapillaryPressureKernelWrappers,
                                                  m_wettingCurve,
                                                  m_landParam[ipWetting],
                                                  phaseVolFraction[ipWetting],
                                                  phaseMinHistoricalVolFraction[ipWetting],
                                                  m_phaseIntermediateMinVolFraction,
                                                  phaseTrappedVolFrac[ipWetting],
                                                  phaseCapPressure[ipWetting],
                                                  dPhaseCapPressure_dPhaseVolFrac[ipWetting][ipWetting]);


     }


  // -- non wetting cure if drainage only
  // gas-oil capillary pressure
  if( !m_phaseHasHysteresis[TTP::INTERMEDIATE_NONWETTING] || phaseVolFraction[ipNonWetting] >= phaseMaxHistoricalVolFraction[ipNonWetting] + flowReversalBuffer )
  {
    //TODO trapped
    phaseTrappedVolFrac[ipNonWetting] = LvArray::math::min( phaseVolFraction[ipNonWetting], m_nonWettingCurve.drainageExtrema );
    phaseCapPressure[ipNonWetting] =
    m_nonWettingIntermediateCapillaryPressureKernelWrappers[ModeIndexType::IMBIBITION].compute( &(phaseVolFraction)[ipNonWetting],
                                                                   &(dPhaseCapPressure_dPhaseVolFrac)[ipNonWetting][ipNonWetting] );


    // when pc is on the gas phase, we need to multiply user input by -1
    // because CompositionalMultiphaseFVM does: pres_gas = pres_oil - pc_og, so we need a negative pc_og
    phaseCapPressure[ipNonWetting] *=-1;
    dPhaseCapPressure_dPhaseVolFrac[ipNonWetting][ipNonWetting] *= -1;
  }
  else
  {

    computeImbibitionNonWettingCapillaryPressure( m_nonWettingIntermediateCapillaryPressureKernelWrappers,
                                                                                                m_nonWettingCurve,
                                                                                                m_landParam[ipNonWetting],
                                                                                                phaseVolFraction[ipNonWetting],
                                                                                                phaseMinHistoricalVolFraction[ipNonWetting],
                                                                                                phaseTrappedVolFrac[ipNonWetting],
                                                                                                phaseCapPressure[ipNonWetting],
                                                                                                dPhaseCapPressure_dPhaseVolFrac[ipNonWetting][ipNonWetting] );
    // when pc is on the gas phase, we need to multiply user input by -1
    // because CompositionalMultiphaseFVM does: pres_gas = pres_oil - pc_og, so we need a negative pc_og
    phaseCapPressure[ipNonWetting] *=-1;
    dPhaseCapPressure_dPhaseVolFrac[ipNonWetting][ipNonWetting] *= -1;
  }




}


void
TableCapillaryPressureHysteresis::KernelWrapper::computeImbibitionNonWettingCapillaryPressure( arrayView1d< TableFunction::KernelWrapper const > const & nonWettingKernelWrapper,
                                                                                               const KilloughHysteresis::HysteresisCurve_t & nonWettingCurve,
                                                                                               const geosx::real64 & landParam,
                                                                                               const geosx::real64 & phaseVolFraction,
                                                                                               const geosx::real64 & phaseMaxHistoricalVolFraction,
                                                                                               geosx::real64 & phaseTrappedVolFrac,
                                                                                               geosx::real64 & phaseCapPressure,
                                                                                               geosx::real64 & dPhaseCapPressure_dPhaseVolFrac ) const
{
  GEOSX_ASSERT(!nonWettingCurve.isWetting());
  real64 const S = phaseVolFraction;
  real64 const Smii = nonWettingCurve.imbibitionExtrema;
  real64 const Smid = nonWettingCurve.drainageExtrema;
  real64 const Smax = nonWettingCurve.oppositeBound;

  if( S >= Smax )
  {
    //above accessible range
    phaseCapPressure = CAP_INF;
    dPhaseCapPressure_dPhaseVolFrac = CAP_INF_DERIV;
  }
  else if( S <= Smid )
  {
    //below accessible range
    phaseCapPressure = - CAP_INF;
    dPhaseCapPressure_dPhaseVolFrac = CAP_INF_DERIV;
  }
  else
  {
    //drainage to imbibition
    real64 dpci_dS, dpcd_dS;
    real64 const pci = nonWettingKernelWrapper[ModeIndexType::IMBIBITION].compute(&S,&dpci_dS);
    real64 const pcd = nonWettingKernelWrapper[ModeIndexType::DRAINAGE].compute(&S,&dpcd_dS);

    // Step 1: get the trapped from wetting data
    real64 const Shy = ( phaseMaxHistoricalVolFraction < Smax ) ? phaseMaxHistoricalVolFraction : Smax;

    real64 Scrt = 0.0;
    m_KilloughKernel.computeTrappedCriticalPhaseVolFraction( nonWettingCurve, Shy, landParam,
                                                             Scrt );

    real64 const E = m_KilloughKernel.getCurvatureParamPc();

    //Set 2. compute F as in (EQ 34.21) F = (1/(Shy-S+E)-1/E) / (1/(Shy - Sgcr +E)-1/E)
    real64 const F = (1./(Shy-S+E) - 1./E) / (1./(Shy - Scrt +E) - 1./E);

    //Step 3. compute dF_dS
    real64 dF_dS = (1./(S*S)) / (1./(Shy - Scrt +E) - 1./E);

    //Step 4. Eventually assemble everything following (EQ. 34.20)
    phaseCapPressure = pcd + F * ( pci - pcd );
    dPhaseCapPressure_dPhaseVolFrac = dpci_dS + F * ( dpci_dS - dpcd_dS);
    dPhaseCapPressure_dPhaseVolFrac += dF_dS * ( pci -pcd );
    //TODO reverse flow

  }
}

void TableCapillaryPressureHysteresis::KernelWrapper::computeDrainageCapillaryPressure( const TableFunction::KernelWrapper & drainageRelpermWrapper,
                                                                                        const geosx::real64 & phaseVolFraction,
                                                                                        geosx::real64 & phaseCapPressure,
                                                                                        geosx::real64 & dPhaseCapPressure_dPhaseVolFrac ) const
{
  phaseCapPressure = drainageRelpermWrapper.compute( &phaseVolFraction,
                                                     &dPhaseCapPressure_dPhaseVolFrac );

}



/// Land coeff (tb refactored out in KilloughHysteresis) and saved cvgd

void TableCapillaryPressureHysteresis::computeLandCoefficient()
{
  // For now, we keep two separate Land parameters for the wetting and non-wetting phases
  // For two-phase flow, we make sure that they are equal
  m_landParam.resize( 2 );

/*  integer ipWetting = 0;
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
  }*/

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
  //extract phase min and max from data structures
  return KilloughModel.createKernelWrapper( m_landParam,
                                            KilloughHysteresis::HysteresisCurve_t(m_wettingPhaseMinVolumeFraction,m_imbibitionWettingPhaseMaxVolumeFraction,m_drainageWettingPhaseMaxVolumeFraction),
                                            KilloughHysteresis::HysteresisCurve_t(m_nonWettingPhaseMaxVolumeFraction,m_imbibitionNonWettingPhaseMinVolumeFraction,m_drainageNonWettingPhaseMinVolumeFraction),
                                            m_phaseTrappedVolFrac );

}

TableCapillaryPressureHysteresis::KernelWrapper
TableCapillaryPressureHysteresis::createKernelWrapper()
{

  // we want to make sure that the wrappers are always up-to-date, so we recreate them everytime
  createAllTableKernelWrappers();

  m_KilloughKernel = createKilloughKernelWrapper();

  // then we create the actual TableRelativePermeabilityHysteresis::KernelWrapper
  return KernelWrapper( m_wettingNonWettingCapillaryPressureKernelWrappers,
                        m_wettingIntermediateCapillaryPressureKernelWrappers,
                        m_nonWettingIntermediateCapillaryPressureKernelWrappers,
                        m_KilloughKernel,
                        m_phaseHasHysteresis,
                        m_landParam,
                        m_phaseIntermediateMinVolFraction,
                        KilloughHysteresis::HysteresisCurve_t(m_wettingPhaseMinVolumeFraction,m_imbibitionWettingPhaseMaxVolumeFraction,m_drainageWettingPhaseMaxVolumeFraction),
                        KilloughHysteresis::HysteresisCurve_t(m_nonWettingPhaseMaxVolumeFraction,m_imbibitionNonWettingPhaseMinVolumeFraction,m_drainageNonWettingPhaseMinVolumeFraction),
                        m_phaseMinHistoricalVolFraction,
                        m_phaseMaxHistoricalVolFraction,
                        m_phaseTypes,
                        m_phaseOrder,
                        m_phaseTrappedVolFrac,
                        m_phaseCapPressure,
                        m_dPhaseCapPressure_dPhaseVolFrac );
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


///kernel ctor
TableCapillaryPressureHysteresis::KernelWrapper::KernelWrapper( arrayView1d< const TableFunction::KernelWrapper > const & wettingNonWettingCapillaryPressureKernelWrappers,
                                                                arrayView1d< const TableFunction::KernelWrapper > const & wettingIntermediateCapillaryPressureKernelWrappers,
                                                                arrayView1d< const TableFunction::KernelWrapper > const & nonWettingIntermediateCapillaryPressureKernelWrappers,
                                                                const KilloughHysteresis::KernelKilloughHysteresisBase & KilloughKernel,
                                                                const arrayView1d< const geosx::integer > & phaseHasHysteresis,
                                                                const arrayView1d< const geosx::real64 > & landParam,
                                                                const geosx::real64 & phaseIntermediateMinVolFraction,
                                                                const KilloughHysteresis::HysteresisCurve_t & wettingCurve,
                                                                const KilloughHysteresis::HysteresisCurve_t & nonWettingCurve,
                                                                const arrayView2d< const geosx::real64, compflow::USD_PHASE > & phaseMinHistoricalVolFraction,
                                                                const arrayView2d< const geosx::real64, compflow::USD_PHASE > & phaseMaxHistoricalVolFraction,
                                                                arrayView1d< integer const > const & phaseTypes,
                                                                arrayView1d< integer const > const & phaseOrder,
                                                                arrayView3d< real64, cappres::USD_CAPPRES > const & phaseTrapped,
                                                                const arrayView3d< geosx::real64, relperm::USD_RELPERM > & phaseCapPressure,
                                                                const arrayView4d< geosx::real64, relperm::USD_RELPERM_DS > & dPhaseCapPressure_dPhaseVolFrac )
                                                                :
  CapillaryPressureBaseUpdate( phaseTypes,
                               phaseOrder,
                               phaseTrapped,
                               phaseCapPressure,
                               dPhaseCapPressure_dPhaseVolFrac ),
                                                                m_wettingNonWettingCapillaryPressureKernelWrappers(wettingNonWettingCapillaryPressureKernelWrappers),
                                                                m_wettingIntermediateCapillaryPressureKernelWrappers(wettingIntermediateCapillaryPressureKernelWrappers),
                                                                m_nonWettingIntermediateCapillaryPressureKernelWrappers(nonWettingIntermediateCapillaryPressureKernelWrappers),
                                                                m_KilloughKernel(KilloughKernel),
                                                                m_phaseHasHysteresis(phaseHasHysteresis),
                                                                m_landParam(landParam),
                                                                m_phaseIntermediateMinVolFraction(phaseIntermediateMinVolFraction),
                                                                m_wettingCurve(wettingCurve),
                                                                m_nonWettingCurve(nonWettingCurve),
                                                                m_phaseMinHistoricalVolFraction(phaseMinHistoricalVolFraction),
                                                                m_phaseMaxHistoricalVolFraction(phaseMaxHistoricalVolFraction)
{}




REGISTER_CATALOG_ENTRY( ConstitutiveBase, TableCapillaryPressureHysteresis, std::string const &, Group * const )

}
} // geosx