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


#ifndef GEOSX_CONSTITUTIVE_TABLECAPILLARYPRESSUREHYSTERESIS_HPP
#define GEOSX_CONSTITUTIVE_TABLECAPILLARYPRESSUREHYSTERESIS_HPP

#include "constitutive/capillaryPressure/CapillaryPressureBase.hpp"
#include "functions/TableFunction.hpp"

#include "constitutive/KilloughHysteresis.hpp"

namespace geosx
{

namespace constitutive
{

class TableCapillaryPressureHysteresis : public CapillaryPressureBase
{
  /// useful constant
  static constexpr real64 CAP_INF = 1e9;
//          std::numeric_limits< real64 >::max();
  static constexpr real64 CAP_INF_DERIV = 1e9;
//          std::numeric_limits< real64 >::max();

public:

  /// order of the phase properties for three-phase flow
  struct ThreePhasePairPhaseType
  {
    enum : integer
    {
      INTERMEDIATE_WETTING = 0,   ///< index for intermediate-wetting
      INTERMEDIATE_NONWETTING = 1 ///< index for intermediate-non-wetting
    };
  };

  enum ModeIndexType : integer
  {
//    enum : integer
//    {
      DRAINAGE = 0,//to be used in array of Kernels
      IMBIBITION = 1,
      DRAINAGE_TO_IMBIBITION = 2,
      IMBIBITION_TO_DRAINAGE = 3
//    };
  };


  TableCapillaryPressureHysteresis( std::string const & name,
                                    dataRepository::Group * const parent );

  static std::string catalogName(){ return "TableCapillaryPressureHysteresis"; }
  virtual string getCatalogName() const override { return catalogName(); }

  ///Kernel
  class KernelWrapper final : public CapillaryPressureBaseUpdate
  {
public:

    KernelWrapper(
      arrayView1d< TableFunction::KernelWrapper const > const & wettingNonWettingCapillaryPressureKernelWrappers,
      arrayView1d< TableFunction::KernelWrapper const > const & wettingIntermediateCapillaryPressureKernelWrappers,
      arrayView1d< TableFunction::KernelWrapper const > const & nonWettingIntermediateCapillaryPressureKernelWrappers,
      KilloughHysteresis::KernelKilloughHysteresisBase const & KilloughKernel,
      arrayView1d< integer const > const & phaseHasHysteresis,
      arrayView1d< real64 const > const & landParam,
      real64 const & phaseIntermediateMinVolFraction,
      KilloughHysteresis::HysteresisCurve_t const & wettingCurve,
      KilloughHysteresis::HysteresisCurve_t const & nonWettingCurve,
      arrayView2d< real64 const, compflow::USD_PHASE > const & phaseMinHistoricalVolFraction,
      arrayView2d< real64 const, compflow::USD_PHASE > const & phaseMaxHistoricalVolFraction,
      arrayView1d< integer const > const & phaseTypes,
      arrayView1d< integer const > const & phaseOrder,
      ModeIndexType & mode,
      arrayView3d< real64, relperm::USD_RELPERM > const & phaseTrappedVolFrac,
      arrayView3d< real64, relperm::USD_RELPERM > const & phaseCapPressure,
      arrayView4d< real64, relperm::USD_RELPERM_DS > const & dPhaseCapPressure_dPhaseVolFrac );

    //actual workers
    GEOSX_HOST_DEVICE
    void computeBoundCapillaryPressure(TableFunction::KernelWrapper const & drainageCapPressureWrapper,
                                       real64 const & phaseVolFraction,
                                       real64 & phaseCapPressure,
                                       real64 & dPhaseCapPressure_dPhaseVolFrac ) const;

    GEOSX_HOST_DEVICE
    void
    computeImbibitionWettingCapillaryPressure(
      const arrayView1d< const TableFunction::KernelWrapper > & wettingKernelWapper,
      const KilloughHysteresis::HysteresisCurve_t & wettingCurve,
      const KilloughHysteresis::HysteresisCurve_t & nonWettingCurve,
      const geosx::real64 & landParam,
      const geosx::real64 & phaseVolFraction, const geosx::real64 & phaseMinHistoricalVolFraction,
      const real64 & phaseIntermediateMinVolFraction, geosx::real64 & phaseTrappedVolFrac,
      geosx::real64 & phaseCapPressure, geosx::real64 & dPhaseCapPressure_dPhaseVolFrac ) const;

    GEOSX_HOST_DEVICE
    void
    computeImbibitionNonWettingCapillaryPressure(
            const arrayView1d<const TableFunction::KernelWrapper> &nonWettingKernelWrapper,
            const KilloughHysteresis::HysteresisCurve_t &nonWettingCurve,
            const KilloughHysteresis::HysteresisCurve_t &wettingCurve, const geosx::real64 &landParam,
            const geosx::real64 &phaseVolFraction, const geosx::real64 &phaseMaxHistoricalVolFraction,
            geosx::real64 &phaseTrappedVolFrac, geosx::real64 &phaseCapPressure,
            geosx::real64 &dPhaseCapPressure_dPhaseVolFrac) const;



    //wrapper call wrt number of phase
    GEOSX_HOST_DEVICE
    void computeTwoPhaseWetting( integer const ipWetting,
                          integer const ipNonWetting,
                          arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
                          arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseMaxHistoricalVolFraction,
                          arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseMinHistoricalVolFraction,
                          arraySlice1d< real64, relperm::USD_RELPERM - 2 > const & phaseTrappedVolFrac,
                          arraySlice1d< real64, relperm::USD_RELPERM - 2 > const & phaseCapPressure,
                          arraySlice2d< real64,
                                        relperm::USD_RELPERM_DS - 2 > const & dPhaseCapPressure_dPhaseVolFrac ) const;


      GEOSX_HOST_DEVICE
      void computeTwoPhaseNonWetting( integer const ipWetting,
                            integer const ipNonWetting,
                            arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
                            arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseMaxHistoricalVolFraction,
                            arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseMinHistoricalVolFraction,
                            arraySlice1d< real64, relperm::USD_RELPERM - 2 > const & phaseTrappedVolFrac,
                            arraySlice1d< real64, relperm::USD_RELPERM - 2 > const & phaseCapPressure,
                            arraySlice2d< real64,
                                    relperm::USD_RELPERM_DS - 2 > const & dPhaseCapPressure_dPhaseVolFrac ) const;

    GEOSX_HOST_DEVICE
    void computeThreePhase( integer const ipWetting,
                            integer const ipInter,
                            integer const ipNonWetting,
                            arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
                            arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseMaxHistoricalVolFraction,
                            arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseMinHistoricalVolFraction,
                            arraySlice1d< real64, relperm::USD_RELPERM - 2 > const & phaseTrappedVolFrac,
                            arraySlice1d< real64, relperm::USD_RELPERM - 2 > const & phaseCapPressure,
                            arraySlice2d< real64, relperm::USD_RELPERM_DS - 2 > const & dPhaseCapPressure_dPhaseVolFrac ) const;

    //uppermost call-wrappers
    GEOSX_HOST_DEVICE
    virtual void compute( arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
                          arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseMaxHistoricalVolFraction,
                          arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseMinHistoricalVolFraction,
                          arraySlice1d< real64, cappres::USD_CAPPRES - 2 > const & phaseTrappedVolFrac,
                          arraySlice1d< real64, cappres::USD_CAPPRES - 2 > const & phaseCapPressure,
                          arraySlice2d< real64, cappres::USD_CAPPRES_DS - 2 > const & dPhaseCapPressure_dPhaseVolFrac ) const;

    GEOSX_HOST_DEVICE
    virtual void update( localIndex const k,
                         localIndex const q,
                         arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction ) const override;

private:


    static constexpr real64 flowReversalBuffer = KilloughHysteresis::KernelKilloughHysteresisBase::flowReversalBuffer;
    ModeIndexType& m_mode;

    //2p
    arrayView1d< TableFunction::KernelWrapper const > const m_wettingNonWettingCapillaryPressureKernelWrappers;
    //3p
    arrayView1d< TableFunction::KernelWrapper const > const m_wettingIntermediateCapillaryPressureKernelWrappers;
    arrayView1d< TableFunction::KernelWrapper const > const m_nonWettingIntermediateCapillaryPressureKernelWrappers;

    KilloughHysteresis::KernelKilloughHysteresisBase const & m_KilloughKernel;

    ///Land Coeff
    arrayView1d< integer const > m_phaseHasHysteresis;
    arrayView1d< real64 const > m_landParam;

    /// needed in 3p-wetting hysteresis as we need to get the max accessible pore space
    real64 const m_phaseIntermediateMinVolFraction;

    KilloughHysteresis::HysteresisCurve_t const m_wettingCurve;
    KilloughHysteresis::HysteresisCurve_t const m_nonWettingCurve;

    /// Minimum historical phase volume fraction for each phase
    arrayView2d< real64 const, compflow::USD_PHASE > m_phaseMinHistoricalVolFraction;

    /// Maximum historical phase volume fraction for each phase
    arrayView2d< real64 const, compflow::USD_PHASE > m_phaseMaxHistoricalVolFraction;


  };

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

  //might need it to be virtual one level higher --> from Killough/Hysteresis common class
  virtual void saveConvergedPhaseVolFractionState( arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFraction ) const override;


  struct viewKeyStruct : CapillaryPressureBase::viewKeyStruct
  {

    ///Killough Kernel
    static constexpr char const * KilloughModelNameString()
    { return "KilloughModelName"; }

    static constexpr char const * KilloughModelWrapperString()
    { return "KilloughWrappers"; }

    // defaulted to 0.1
//    static constexpr char const * killoughCurvatureParameterString() { return "killoughCurvatureParameter"; }
    //used to compute different re-traversal path going drainage-imbibition-drainage
//    static constexpr char const * tCurveOptionString() { return "tCurveOption";};
    ///Land Coeff
    static constexpr char const * landParameterString()
    { return "landParameter"; }

///pivot points
    static constexpr char const * wettingPhaseMinVolumeFractionString()
    { return "wettingPhaseMinVolumeFraction"; }
    static constexpr char const * drainageWettingPhaseMaxVolumeFractionString()
    { return "drainageWettingPhaseMaxVolumeFraction"; }
    static constexpr char const * imbibitionWettingPhaseMaxVolumeFractionString()
    { return "imbibitionWettingPhaseMaxVolumeFraction"; }
    static constexpr char const * nonWettingPhaseMaxVolumeFractionString()
    { return "nonWettingPhaseMaxVolumeFraction"; }
    static constexpr char const * drainageNonWettingPhaseMinVolumeFractionString()
    { return "drainageNonWettingPhaseMinVolumeFraction"; }
    static constexpr char const * imbibitionNonWettingPhaseMinVolumeFractionString()
    { return "imbibitionNonWettingPhaseMinVolumeFraction"; }

    ///flag
    static constexpr char const * phaseHasHysteresisString()
    { return "phaseHasHysteresis"; }

    ///tables and assoc. wrappers
    //2phase
    static constexpr char const * drainageWettingNonWettingCapPresTableNameString()
    { return "drainageWettingNonWettingCapPressureTableName"; }
    static constexpr char const * imbibitionWettingNonWettingCapPresTableNameString()
    { return "imbibitionWettingNonWettingCapPressureTableName"; }
    //3phase
    static constexpr char const * drainageWettingIntermediateCapPresTableNameString()
    { return "drainageWettingIntermediateCapPressureTableName"; }
    static constexpr char const * drainageNonWettingIntermediateCapPresTableNameString()
    { return "drainageNonWettingIntermediateCapPressureTableName"; }
    static constexpr char const * imbibitionWettingIntermediateCapPresTableNameString()
    { return "imbibitionWettingIntermediateCapPressureTableName"; }
    static constexpr char const * imbibitionNonWettingIntermediateCapPresTableNameString()
    { return "imbibitionNonWettingIntermediateCapPressureTableName"; }
    static constexpr char const * wettingNonWettingCapillaryPressureKernelWrappersString()
    { return "wettingNonWettingCapillaryPressureKernelWrappers"; }
    static constexpr char const * wettingIntermediateCapillaryPressureKernelWrappersString()
    { return "wettingIntermediateCapillaryPressureKernelWrappers"; }
    static constexpr char const * nonWettingIntermediateCapillaryPressureKernelWrappersString()
    { return "nonWettingIntermediateCapillaryPressureKernelWrappers"; }

    //misc
    static constexpr char const * phaseIntermediateMinVolFractionString()
    { return "phaseIntermediateMinVolFraction";}
    //to decide wheter drainage/drainage to imbibition or imbibition/imbibition to drainage
      static constexpr char const * modeTypeString()
      { return "modeType";}

  };


private:
  virtual void postProcessInput() override;

  virtual void initializePreSubGroups() override;

  void resizeFields( localIndex const size,
                     localIndex const numPts ) override;


  /**
   * @brief Create all the table kernel wrappers needed for the simulation (for all the phases present)
   */
  void createAllTableKernelWrappers();

  KilloughHysteresis::KernelKilloughHysteresisBase
  createKilloughKernelWrapper( const KilloughHysteresis::HysteresisCurve_t & wettingCurve,
                               const KilloughHysteresis::HysteresisCurve_t & nonWettingCurve );

  /**
   * @brief Compute the Land coefficient for the wetting and non-wetting phases
   */
  void computeLandCoefficient();

  ///data members

  // Hysteresis parameters
  KilloughHysteresis::KernelKilloughHysteresisBase m_KilloughKernel;
  string m_KilloughModelName;

  ModeIndexType m_mode;

  //TODO impl
//  array1d< integer >  m_tCurveOption;

  //might be further packed in CapillaryCurve_t
  real64 m_wettingPhaseMinVolumeFraction;
  real64 m_drainageWettingPhaseMaxVolumeFraction;
  real64 m_imbibitionWettingPhaseMaxVolumeFraction;

  real64 m_nonWettingPhaseMaxVolumeFraction;
  real64 m_drainageNonWettingPhaseMinVolumeFraction;
  real64 m_imbibitionNonWettingPhaseMinVolumeFraction;

  ///tables
  //2p
  string m_drainageWettingNonWettingCapPresTableName;
  string m_imbibitionWettingNonWettingCapPresTableName;
  //3p
  string m_drainageWettingIntermediateCapPresTableName;
  string m_drainageNonWettingIntermediateCapPresTableName;
  string m_imbibitionWettingIntermediateCapPresTableName;
  string m_imbibitionNonWettingIntermediateCapPresTableName;
  // kernel wrappers
  /// Imbibition kernel wrappers for relative permeabilities in the following order:
  /// 0- drainage
  /// 1- imbibition (cf. struct ModeIndexType)
  //2p
  array1d< TableFunction::KernelWrapper > m_wettingNonWettingCapillaryPressureKernelWrappers;
  //3p
  array1d< TableFunction::KernelWrapper > m_wettingIntermediateCapillaryPressureKernelWrappers;
  array1d< TableFunction::KernelWrapper > m_nonWettingIntermediateCapillaryPressureKernelWrappers;


  /// Flag to specify whether the phase has hysteresis or not (deduced from table input)
  array1d< integer > m_phaseHasHysteresis;

  /// Trapping parameter from the Land model (typically called C)
  array1d< real64 > m_landParam;

  // Max historical saturations
  /// Minimum historical phase volume fraction for each phase
  array2d< real64, compflow::LAYOUT_PHASE > m_phaseMinHistoricalVolFraction;

  /// Maximum historical phase volume fraction for each phase
  array2d< real64, compflow::LAYOUT_PHASE > m_phaseMaxHistoricalVolFraction;

  //needed in hysteresis of wetting phase
  real64 m_phaseIntermediateMinVolFraction;

};

GEOSX_HOST_DEVICE
inline void TableCapillaryPressureHysteresis::KernelWrapper::compute( arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
                                                                      arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseMaxHistoricalVolFraction,
                                                                      arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseMinHistoricalVolFraction,
                                                                      arraySlice1d< real64, cappres::USD_CAPPRES - 2 > const & phaseTrappedVolFrac,
                                                                      arraySlice1d< real64, cappres::USD_CAPPRES - 2 > const & phaseCapPressure,
                                                                      arraySlice2d< real64, cappres::USD_CAPPRES_DS - 2 > const & dPhaseCapPressure_dPhaseVolFrac

                                                                      ) const
{
  LvArray::forValuesInSlice( dPhaseCapPressure_dPhaseVolFrac, []( real64 & val ){ val = 0.0; } );

  using PT = CapillaryPressureBase::PhaseType;
  integer const ipWater = m_phaseOrder[PT::WATER];
  integer const ipOil   = m_phaseOrder[PT::OIL];
  integer const ipGas   = m_phaseOrder[PT::GAS];

  if( ipWater >= 0 && ipOil >= 0 && ipGas >= 0 )
  {
    computeThreePhase( ipWater, // wetting
                       ipOil,   // intermediate
                       ipGas,   // non-wetting
                       phaseVolFraction,
                       phaseMaxHistoricalVolFraction,
                       phaseMinHistoricalVolFraction,
                       phaseTrappedVolFrac,
                       phaseCapPressure,
                       dPhaseCapPressure_dPhaseVolFrac );

  }
  else if( ipWater < 0 )
  {
      computeTwoPhaseNonWetting(ipOil, // leading
                                ipGas, // deduced
                                phaseVolFraction,
                                phaseMaxHistoricalVolFraction,
                                phaseMinHistoricalVolFraction,
                                phaseTrappedVolFrac,
                                phaseCapPressure,
                                dPhaseCapPressure_dPhaseVolFrac);
  }
  else if( ipOil < 0 )
  {
    computeTwoPhaseWetting( ipWater, // leading
                     ipGas,   // deduced
                     phaseVolFraction,
                     phaseMaxHistoricalVolFraction,
                     phaseMinHistoricalVolFraction,
                     phaseTrappedVolFrac,
                     phaseCapPressure,
                     dPhaseCapPressure_dPhaseVolFrac );
  }
  else if( ipGas < 0 )
  {
      computeTwoPhaseWetting(ipWater, //leading
                             ipOil,   //deduced
                             phaseVolFraction,
                             phaseMaxHistoricalVolFraction,
                             phaseMinHistoricalVolFraction,
                             phaseTrappedVolFrac,
                             phaseCapPressure,
                             dPhaseCapPressure_dPhaseVolFrac);
  }


}

GEOSX_HOST_DEVICE
inline void TableCapillaryPressureHysteresis::KernelWrapper::update( const geosx::localIndex k,
                                                                     const geosx::localIndex q,
                                                                     const arraySlice1d< const geosx::real64,
                                                                                         compflow::USD_PHASE
                                                                                         - 1 > & phaseVolFraction ) const
{
  compute( phaseVolFraction,
           m_phaseMaxHistoricalVolFraction[k],
           m_phaseMinHistoricalVolFraction[k],
           m_phaseTrappedVolFrac[k][q],
           m_phaseCapPressure[k][q],
           m_dPhaseCapPressure_dPhaseVolFrac[k][q] );
}


} //constitutive
} // geosx

#endif //GEOSX_CONSTITUTIVE_TABLECAPILLARYPRESSUREHYSTERESIS_HPP
