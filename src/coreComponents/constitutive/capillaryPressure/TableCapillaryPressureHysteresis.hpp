//
// Created by root on 11/8/22.
//

#ifndef GEOSX_TABLECAPILLARYPRESSUREHYSTERESIS_HPP
#define GEOSX_TABLECAPILLARYPRESSUREHYSTERESIS_HPP

#include "constitutive/capillaryPressure/CapillaryPressureBase.hpp"
#include "functions/TableFunction.hpp"

#include "constitutive/KilloughHysteresis.hpp"

namespace geosx
{

namespace constitutive
{

class TableCapillaryPressureHysteresis : public CapillaryPressureBase
{

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

  struct ModeIndexType
  {
    enum: integer
    {
      DRAINAGE = 0,//to be used in array of Kernels
      IMBIBITION = 1
    };
  };



  TableCapillaryPressureHysteresis(std::string const & name, dataRepository::Group * const parent);
  static std::string catalogName() { return "TableCapillaryPressureHysteresis"; }
  virtual string getCatalogName() const override { return catalogName(); }

  ///Kernel
  class KernelWrapper final : public CapillaryPressureBaseUpdate
  {
  public:

    KernelWrapper();

    //actual workers
    GEOSX_HOST_DEVICE
    void computeImbibitionNonWettingCapillaryPressure();

    //wrapper call wrt number of phase
    GEOSX_HOST_DEVICE
    void computeTwoPhase();

    GEOSX_HOST_DEVICE
    void computeThreePhase();

    //uppermost call-warppers
    GEOSX_HOST_DEVICE
    virtual void compute() const;

    GEOSX_HOST_DEVICE
    virtual void update( localIndex const k,
                         localIndex const q,
                         arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction ) const override;


  };

  /**
 * @brief Create an update kernel wrapper.
 * @return the wrapper
 */
  KernelWrapper createKernelWrapper();

  KilloughHysteresis::KernelKilloughHysteresisBase createKilloughKernelWrapper();

  //might need it to be virtual one level higher --> from Killough/Hysteresis common class
  virtual void saveConvergedPhaseVolFractionState( arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFraction ) const;

  
  
  struct viewKeyStruct : CapillaryPressureBase::viewKeyStruct
  {

    ///Killough Kernel
    static constexpr char const * KilloughModelNameString() { return "KilloughModelName"; }
    static constexpr char const * KilloughModelWrapperString() { return "KilloughWrappers"; }
    // defaulted to 0.1
//    static constexpr char const * killoughCurvatureParameterString() { return "killoughCurvatureParameter"; }
    //used to compute different re-traversal path going drainage-imbibition-drainage
//    static constexpr char const * tCurveOptionString() { return "tCurveOption";};
    ///Land Coeff
    static constexpr char const * landParameterString() { return "landParameter"; }


///pivot points
    static constexpr char const * wettingPhaseMinVolumeFractionString() { return "wettingPhaseMinVolumeFraction"; }
    static constexpr char const * drainageWettingPhaseMaxVolumeFractionString() { return "drainageWettingPhaseMaxVolumeFraction"; }
    static constexpr char const * imbibitionWettingPhaseMaxVolumeFractionString() { return "imbibitionWettingPhaseMaxVolumeFraction"; }

    static constexpr char const * nonWettingPhaseMaxVolumeFractionString() { return "nonWettingPhaseMaxVolumeFraction"; }
    static constexpr char const * drainageNonWettingPhaseMinVolumeFractionString() { return "drainageNonWettingPhaseMinVolumeFraction"; }
    static constexpr char const * imbibitionNonWettingPhaseMinVolumeFractionString() { return "imbibitionNonWettingPhaseMinVolumeFraction"; }

    ///flag
    static constexpr char const * phaseHasHysteresisString() { return "phaseHasHysteresis"; }

    ///tables and assoc. wrappers
    //2phase
    static constexpr char const * drainageWettingNonWettingCapPresTableNameString() { return "drainageWettingNonWettingCapPressureTableName"; }
    static constexpr char const * imbibitionWettingNonWettingCapPresTableNameString() { return "imbibitionWettingNonWettingCapPressureTableName"; }
    //3phase
    static constexpr char const * drainageWettingIntermediateCapPresTableNameString() { return "drainageWettingIntermediateCapPressureTableName"; }
    static constexpr char const * drainageNonWettingIntermediateCapPresTableNameString() { return "drainageNonWettingIntermediateCapPressureTableName"; }
    static constexpr char const * imbibitionWettingIntermediateCapPresTableNameString() { return "imbibitionWettingIntermediateCapPressureTableName"; }
    static constexpr char const * imbibitionNonWettingIntermediateCapPresTableNameString() { return "imbibitionNonWettingIntermediateCapPressureTableName"; }

    static constexpr char const * drainageCapPresWrappersString() { return "drainageCapPresWrappers"; }
    static constexpr char const * imbibitionCapPresWrappersString() { return "imbibitionCapPresWrappers"; }

  };


private:


  void resizeFields( localIndex const size, localIndex const numPts ) override;

  virtual void postProcessInput() override;

  virtual void initializePreSubGroups() override;

  /**
   * @brief Create all the table kernel wrappers needed for the simulation (for all the phases present)
   */
  void createAllTableKernelWrappers();

  KilloughHysteresis::KernelKilloughHysteresisBase createKilloughKernelWrapper();

  /**
 * @brief Compute the Land coefficient for the wetting and non-wetting phases
 */
  void computeLandCoefficient();

  ///data members

  /// Minimum volume fraction for each phase in drainage (deduced from the drainage table)
  arrayView1d< real64 const > m_drainagePhaseMinVolFraction;
  /// Minimum volume fraction for each phase in imbibition (deduced from the imbibition table)
  arrayView1d< real64 const > m_imbibitionPhaseMinVolFraction;


  // Hysteresis parameters
  KilloughHysteresis::KernelKilloughHysteresisBase m_KilloughKernel;
  string m_KilloughModelName;
  //TODO impl
//  array1d< integer >  m_tCurveOption;

  // kernel wrappers
  /// Imbibition kernel wrappers for relative permeabilities in the following order:
  /// 0- drainage
  /// 1- imbibition (cf. struct ModeIndexType)
  //2p
  array1d< TableFunction::KernelWrapper > m_wettingNonWettingCapillaryPressureKernelWrappers;
  //3p
  array1d< TableFunction::KernelWrapper > m_wettingIntermediateCapillaryPressureKernelWrappers;
  array1d< TableFunction::KernelWrapper > m_nonWettingIntermediateCapillaryPressureKernelWrappers;



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
  /// Flag to specify whether the phase has hysteresis or not (deduced from table input)
  array1d< integer > m_phaseHasHysteresis;

  /// Trapping parameter from the Land model (typically called C)
  array1d< real64 > m_landParam;

  // Max historical saturations
  /// Minimum historical phase volume fraction for each phase
  array2d< real64, compflow::LAYOUT_PHASE > m_phaseMinHistoricalVolFraction;

  /// Maximum historical phase volume fraction for each phase
  array2d< real64, compflow::LAYOUT_PHASE > m_phaseMaxHistoricalVolFraction;

};

GEOSX_HOST_DEVICE
inline void TableCapillaryPressureHysteresis::KernelWrapper::compute() const
{}

GEOSX_HOST_DEVICE
inline void TableCapillaryPressureHysteresis::KernelWrapper::update( const geosx::localIndex k,
                                                                     const geosx::localIndex q,
                                                                     const arraySlice1d< const geosx::real64,
                                                                       compflow::USD_PHASE
                                                                       - 1 > & phaseVolFraction ) const
{
  compute();
}




} //constitutive
} // geosx

#endif //GEOSX_TABLECAPILLARYPRESSUREHYSTERESIS_HPP
