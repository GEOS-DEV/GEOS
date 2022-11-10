//
// Created by root on 11/8/22.
//

#ifndef GEOSX_TABLECAPILLARYPRESSUREHYSTERESIS_HPP
#define GEOSX_TABLECAPILLARYPRESSUREHYSTERESIS_HPP

#include "constitutive/capillaryPressure/CapillaryPressureBase.hpp"
#include "constitutive/capillaryPressure/KilloughHysteresisCapillaryPressure.hpp"
#include "functions/TableFunction.hpp"

namespace geosx
{

namespace constitutive
{

class TableCapillaryPressureHysteresis : public CapillaryPressureBase, virtual public KilloughHysteresisCapillaryPressure
{

public:

  TableCapillaryPressureHysteresis(std::string const & name, dataRepository::Group * const parent);
  static std::string catalogName() { return "TableCapillaryPressureHysteresis"; }

  virtual string getCatalogName() const override { return catalogName(); }

  ///Kernel
  class KernelWrapper final : public CapillaryPressureBaseUpdate, virtual public KilloughHysteresisCapillaryPressure
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

  //might need it to be virtual one level higher --> from Killough/Hysteresis common class
  virtual void saveConvergedPhaseVolFractionState( arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFraction ) const;

  
  
  struct viewKeyStruct : CapillaryPressureBase::viewKeyStruct
  {
    static constexpr char const * drainagePhaseMinVolumeFractionString() { return "drainagePhaseMinVolumeFraction"; }
    static constexpr char const * imbibitionPhaseMinVolumeFractionString() { return "imbibitionPhaseMinVolumeFraction"; }
    //2phase
    static constexpr char const * drainageWettingNonWettingCapPresTableNameString() { return "drainageWettingNonWettingCapPressureTableName"; }
    static constexpr char const * imbibitionWettingNonWettingCapPresTableNameString() { return "imbibitionWettingNonWettingCapPressureTableName"; }
    //3phase
    static constexpr char const * drainageWettingIntermediateCapPresTableNameString() { return "drainageWettingIntermediateCapPressureTableName"; }
    static constexpr char const * drainageNonWettingIntermediateCapPresTableNameString() { return "drainageNonWettingIntermediateCapPressureTableName"; }
    static constexpr char const * imbibitionWettingIntermediateCapPresTableNameString() { return "imbibitionWettingIntermediateCapPressureTableName"; }
    static constexpr char const * imbibitionNonWettingIntermediateCapPresTableNameString() { return "imbibitionNonWettingIntermediateCapPressureTableName"; }

    // ?
    static constexpr char const * drainageCapPresWrappersString() { return "drainageCapPresWrappers"; }
    static constexpr char const * imbibitionCapPresWrappersString() { return "imbibitionCapPresWrappers"; }

    // defaulted to 0.1
    static constexpr char const * killoughCurvatureParameterString() { return "killoughCurvatureParameter"; }
    //used to compute different re-traversal path going drainage-imbibition-drainage
    static constexpr char const * tCurveOptionString() { return "tCurveOption";};
  };


private:

  virtual void postProcessInput() override;

  virtual void initializePreSubGroups() override;
  
  /**
   * @brief Create all the table kernel wrappers needed for the simulation (for all the phases present)
   */
  void createAllTableKernelWrappers();

  /// TODO data members

  /// Minimum volume fraction for each phase in drainage (deduced from the drainage table)
  arrayView1d< real64 const > m_drainagePhaseMinVolFraction;
  /// Minimum volume fraction for each phase in imbibition (deduced from the imbibition table)
  arrayView1d< real64 const > m_imbibitionPhaseMinVolFraction;

  // kernel wrappers
  // 2p :
  // ... etc
  array1d< TableFunction::KernelWrapper > m_drainageCapillaryPressureKernelWrappers;

  /// Imbibition kernel wrappers for relative permeabilities in the following order:
  ///  0- wetting-phase
  ///  1- non-wetting-phase
  array1d< TableFunction::KernelWrapper > m_imbibitionCapillaryPressureKernelWrappers;

  ///tables
 //2p
  array1d< string> m_drainageWettingNonWettingCapPresTableName;
 array1d< string > m_imbibitionWettingNonWettingCapPresTableName;
 //3p
 array1d< string > m_drainageWettingIntermediateCapPresTableName;
 array1d< string > m_drainageNonWettingIntermediateCapPresTableName;
 array1d< string > m_imbibitionWettingIntermediateCapPresTableName;
 array1d< string > m_imbibitionNonWettingIntermediateCapPresTableName;
  /// Flag to specify whether the phase has hysteresis or not (deduced from table input)
  array1d< integer > m_phaseHasHysteresis;

  // Max historical saturations
  /// Minimum historical phase volume fraction for each phase
  array2d< real64, compflow::LAYOUT_PHASE > m_phaseMinHistoricalVolFraction;

  /// Maximum historical phase volume fraction for each phase
  array2d< real64, compflow::LAYOUT_PHASE > m_phaseMaxHistoricalVolFraction;

  //misc
  real64 m_killoughCurvaturePc;
  array1d< integer >  m_tCurveOption;



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
