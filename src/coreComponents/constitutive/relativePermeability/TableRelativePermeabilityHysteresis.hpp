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

/**
 * @file TableRelativePermeabilityHysteresis.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_TABLERELATIVEPERMEABILITYHYSTERESIS_HPP
#define GEOSX_CONSTITUTIVE_TABLERELATIVEPERMEABILITYHYSTERESIS_HPP

#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityInterpolators.hpp"
#include "functions/TableFunction.hpp"

namespace geosx
{
namespace constitutive
{

class TableRelativePermeabilityHysteresis : public RelativePermeabilityBase
{
public:

  /// order of the phase properties for three-phase flow
  struct ThreePhasePairPhaseType
  {
    enum : integer
    {
      WETTING = 0,                ///< wetting phase property
      INTERMEDIATE_WETTING = 1,   ///< intermediate phase property
      NONWETTING = 2,             ///< non-wetting phase property
      INTERMEDIATE_NONWETTING = 3 ///< intermediate phase property
    };
  };

  /// order of the phase properties in the wetting-non-wetting data
  struct TwoPhasePairPhaseType
  {
    enum : integer
    {
      WETTING = 0,   ///< wetting phase property
      NONWETTING = 1 ///< non-wetting phase property
    };
  };

  /// order of the phase properties in the imbibition data
  struct ImbibitionPhasePairPhaseType
  {
    enum : integer
    {
      WETTING = 0,   ///< wetting phase property
      NONWETTING = 1 ///< non-wetting phase property
    };
  };


  TableRelativePermeabilityHysteresis( std::string const & name, dataRepository::Group * const parent );

  static std::string catalogName() { return "TableRelativePermeabilityHysteresis"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  class KernelWrapper final : public RelativePermeabilityBaseUpdate
  {
public:

    /// To avoid division by zero, this is the min Scrd-Scri used in the computation of the Land constant
    static constexpr real64 minScriMinusScrd = 1e-12;

    /// To avoid frequent changes from drainage to imbibition and vice versa, we use this buffer
    static constexpr real64 flowReversalBuffer = 1e-12;


    /**
     * @brief Constructor for the kernel wrapper updating the relative permeabilities
     * @param[in] drainageRelPermKernelWrappers kernel wrappers storing the drainage relperms (see below for the distinction between 2 and 3
     * phase flow)
     * @param[in] imbibitionRelPermKernelWrappers kernel wrappers storing the imbibition relperms (see below for the distinction between 2
     * and 3 phase flow)
     * @param[in] jerauldParam_a first (modification) parameter proposed by Jerauld
     * @param[in] jerauldParam_b second (exponent) parameter proposed by Jerauld
     * @param[in] killoughCurvatureParam curvature parameter proposed by Killough
     * @param[in] phaseHasHysteresis flag indicating whether a phase has hysteresis or not
     * @param[in] landParam Land trapping parameter
     * @param[in] drainageMinPhaseVolFraction drainage minimum volume fraction for each phase
     * @param[in] imbibitionMinPhaseVolFraction imbibition minimum volume fraction for the wetting and non-wetting phase
     * @param[in] drainageMaxPhaseVolFraction drainage maximum volume fraction for each phase
     * @param[in] imbibitionMaxPhaseVolFraction imbibition maximum volume fraction for the wetting and non-wetting phase
     * @param[in] drainageRelPermEndPoint drainage end-point relperm for each phase
     * @param[in] imbibitionRelPermEndPoint imbibition end-point relperm for the wetting and non-wetting phase
     * @param[in] phaseTypes the phase types
     * @param[in] phaseOrder the phase order
     * @param[in] phaseMinHistoricalPhaseVolFraction minimum historical saturation for each phase
     * @param[in] phaseMaxHistoricalPhaseVolFraction maximum historical saturation for each phase
     * @param[out] phaseRelPerm relative permeability for each phase
     * @param[out] dPhaseRelPerm_dPhaseVolFrac derivative of relative permeability wrt phase volume fraction for each phase
     */
    KernelWrapper( arrayView1d< TableFunction::KernelWrapper const > const & drainageRelPermKernelWrappers,
                   arrayView1d< TableFunction::KernelWrapper const > const & imbibitionRelPermKernelWrappers,
                   real64 const & jerauldParam_a,
                   real64 const & jerauldParam_b,
                   real64 const & killoughCurvatureParam,
                   arrayView1d< integer const > const & phaseHasHysteresis,
                   arrayView1d< real64 const > const & landParam,
                   arrayView1d< real64 const > const & drainageMinPhaseVolFraction,
                   arrayView1d< real64 const > const & imbibitionMinPhaseVolFraction,
                   arrayView1d< real64 const > const & drainageMaxPhaseVolFraction,
                   arrayView1d< real64 const > const & imbibitionMaxPhaseVolFraction,
                   arrayView1d< real64 const > const & drainageRelPermEndPoint,
                   arrayView1d< real64 const > const & imbibitionRelPermEndPoint,
                   arrayView1d< integer const > const & phaseTypes,
                   arrayView1d< integer const > const & phaseOrder,
                   arrayView2d< real64 const, compflow::USD_PHASE > const & phaseMinHistoricalVolFraction,
                   arrayView2d< real64 const, compflow::USD_PHASE > const & phaseMaxHistoricalVolFraction,
                   arrayView3d< real64, relperm::USD_RELPERM > const & phaseRelPerm,
                   arrayView4d< real64, relperm::USD_RELPERM_DS > const & dPhaseRelPerm_dPhaseVolFrac );

    /**
     * @brief Function updating the relperm (and derivative) for a phase using the drainage table
     * @param[in] drainageRelPermKernelWrapper kernel wrapper storing the drainage relperm table for the phase we want to update here
     * @param[in] phaseVolFraction volume fraction of the phase we want to update here
     * @param[out] phaseRelPerm relative permeability of the phase we want to update here
     * @param[out] dPhaseRelPerm_dPhaseVolFrac derivative of the relative permeability wrt phase volume fraction for the phase we want to
     * update here
     */
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void computeDrainageRelPerm( TableFunction::KernelWrapper const & drainageRelPermKernelWrapper,
                                 real64 const & phaseVolFraction,
                                 real64 & phaseRelPerm,
                                 real64 & dPhaseRelPerm_dPhaseVolFrac ) const;

    /**
     * @brief Function updating the relperm (and derivative) for the wetting phase in imbibition using Killough's method
     * @param[in] drainageRelPermKernelWrapper kernel wrapper storing the drainage relperm table for the wetting phase
     * @param[in] imbibitionRelPermKernelWrapper kernel wrapper storing the imbibition relperm table for the wetting phase
     * @param[in] jerauldParam_a first (modification) parameter proposed by Jerauld
     * @param[in] jerauldParam_b second (exponent) parameter proposed by Jerauld
     * @param[in] landParam Land trapping parameter
     * @param[in] phaseVolFraction volume fraction for this phase
     * @param[in] phaseMinHistoricalVolFraction min historical volume fraction for this phase
     * @param[in] imbibitionPhaseMinWettingVolFraction imbibition minimum volume fraction for this phase
     * @param[in] drainagePhaseMaxVolFraction drainage maximum volume fraction for this phase
     * @param[in] imbibitionPhaseMaxVolFraction imbibition maximum volume fraction for this phase
     * @param[in] drainageRelPermEndPoint drainage end-point relperm for this phase
     * @param[in] imbibitionRelPermEndPoint imbibition end-point relperm for this phase
     * @param[out] phaseRelPerm relative permeability of the wetting phase
     * @param[out] dPhaseRelPerm_dPhaseVolFrac derivative of the relative permeability wrt phase volume fraction for the wetting phase
     */
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void computeImbibitionWettingRelPerm( TableFunction::KernelWrapper const & drainageRelPermKernelWrapper,
                                          TableFunction::KernelWrapper const & imbibitionRelPermKernelWrapper,
                                          real64 const & jerauldParam_a,
                                          real64 const & jerauldParam_b,
                                          real64 const & landParam,
                                          real64 const & phaseVolFraction,
                                          real64 const & phaseMinHistoricalVolFraction,
                                          real64 const & imbibitionPhaseMinWettingVolFraction,
                                          real64 const & drainagePhaseMaxVolFraction,
                                          real64 const & imbibitionPhaseMaxVolFraction,
                                          real64 const & drainageRelPermEndPoint,
                                          real64 const & imbibitionRelPermEndPoint,
                                          real64 & phaseRelPerm,
                                          real64 & dPhaseRelPerm_dPhaseVolFrac ) const;

    /**
     * @brief Function updating the relperm (and derivative) for the non-wetting phase in imbibition using Killough's method
     * @param[in] drainageRelPermKernelWrapper kernel wrapper storing the drainage relperm table for the non-wetting phase
     * @param[in] imbibitionRelPermKernelWrapper kernel wrapper storing the imbibition relperm table for the non-wetting phase
     * @param[in] jerauldParam_a first (modification) parameter proposed by Jerauld
     * @param[in] jerauldParam_b second (exponent) parameter proposed by Jerauld
     * @param[in] landParam Land trapping coefficient
     * @param[in] phaseVolFraction volume fraction for this phase
     * @param[in] phaseMaxHistoricalVolFraction max historical volume fraction for this phase
     * @param[in] drainageMinPhaseVolFraction min drainage volume fraction for this phase
     * @param[in] imbibitionMinPhaseVolFraction min imbibition volume fraction for this phase
     * @param[in] drainageMaxPhaseVolFraction max drainage volume fraction for this phase
     * @param[in] drainageRelPermEndPoint drainage end-point relperm for this phase
     * @param[out] phaseRelPerm relative permeability of the non-wetting phase
     * @param[out] dPhaseRelPerm_dPhaseVolFrac derivative of the relative permeability wrt phase volume fraction for the non-wetting phase
     */
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void computeImbibitionNonWettingRelPerm( TableFunction::KernelWrapper const & drainageRelPermKernelWrapper,
                                             TableFunction::KernelWrapper const & imbibitionRelPermKernelWrapper,
                                             real64 const & jerauldParam_a,
                                             real64 const & jerauldParam_b,
                                             real64 const & landParam,
                                             real64 const & phaseVolFraction,
                                             real64 const & phaseMaxHistoricalVolFraction,
                                             real64 const & drainageMinPhaseVolFraction,
                                             real64 const & imbibitionMinPhaseVolFraction,
                                             real64 const & drainageMaxPhaseVolFraction,
                                             real64 const & drainageRelPermEndPoint,
                                             real64 & phaseRelPerm,
                                             real64 & dPhaseRelPerm_dPhaseVolFrac ) const;

    /**
     * @brief Function updating all the phase relperms (and derivatives) for two-phase flow
     * @param[in] ipWetting
     * @param[in] ipNonWetting
     * @param[in] phaseVolFraction
     * @param[in] phaseMaxHistoricalVolFraction
     * @param[in] phaseMinHistoricalVolFraction
     * @param[out] phaseRelPerm
     * @param[out] dPhaseRelPerm_dPhaseVolFrac
     * @detail depending of the flow direction for a given phase, this function updates the phase relative permeability
     *         using computeDrainageRelPerm (in drainage) or using one of the imbibition update functions implementing Killough's method
     */
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void computeTwoPhase( integer const ipWetting,
                          integer const ipNonWetting,
                          arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
                          arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseMaxHistoricalVolFraction,
                          arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseMinHistoricalVolFraction,
                          arraySlice1d< real64, relperm::USD_RELPERM - 2 > const & phaseRelPerm,
                          arraySlice2d< real64, relperm::USD_RELPERM_DS - 2 > const & dPhaseRelPerm_dPhaseVolFrac ) const;

    /**
     * @brief Function updating all the phase relperms (and derivatives) for three-phase flow
     * @param[in] ipWetting index of the wetting phase
     * @param[in] ipInter index of the intermediate phase (oil)
     * @param[in] ipNonWetting index of the non-wetting phase
     * @param[in] phaseVolFraction volume fractions for the three phases
     * @param[in] phaseMaxHistoricalVolFraction max historical volume fractions for the three phases
     * @param[in] phaseMinHistoricalVolFraction min historical volume fractions for the three phases
     * @param[out] phaseRelPerm relative permeabilities for the three phases
     * @param[out] dPhaseRelPerm_dPhaseVolFrac derivatives of relative permeabilities wrt phase volume fraction for the three phases
     * @detail depending of the flow direction for a given phase, this function updates the phase relative permeability
     *         using computeDrainageRelPerm (in drainage) or using one of the imbibition update functions implementing Killough's method
     */
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void computeThreePhase( integer const ipWetting,
                            integer const ipInter,
                            integer const ipNonWetting,
                            arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
                            arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseMaxHistoricalVolFraction,
                            arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseMinHistoricalVolFraction,
                            arraySlice1d< real64, relperm::USD_RELPERM - 2 > const & phaseRelPerm,
                            arraySlice2d< real64, relperm::USD_RELPERM_DS - 2 > const & dPhaseRelPerm_dPhaseVolFrac ) const;

    /**
     * @brief Main function updating all the phase relperms (and derivatives)
     * @param[in] phaseVolFraction volume fractions for all the phases
     * @param[in] phaseMaxHistoricalVolFraction max historical volume fractions for all the phases
     * @param[in] phaseMinHistoricalVolFraction min historical volume fractions for all the phases
     * @param[out] phaseRelPerm relative permeabilities for all the phases
     * @param[out] dPhaseRelPerm_dPhaseVolFrac derivatives of relative permeabilities wrt phase volume fraction for all the phases
     */
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    void compute( arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
                  arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseMaxHistoricalVolFraction,
                  arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseMinHistoricalVolFraction,
                  arraySlice1d< real64, relperm::USD_RELPERM - 2 > const & phaseRelPerm,
                  arraySlice2d< real64, relperm::USD_RELPERM_DS - 2 > const & dPhaseRelPerm_dPhaseVolFrac ) const;

    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    virtual void update( localIndex const k,
                         localIndex const q,
                         arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction ) const override;

private:

    /// Drainage kernel wrappers for relative permeabilities in the following order:
    /// Two-phase flow:
    ///  0- wetting-phase
    ///  1- non-wetting-phase
    /// Three-phase flow:
    ///  0- wetting-phase
    ///  1- intermediate phase (wetting-intermediate data)
    ///  2- non-wetting-phase
    ///  3- intermediate phase (non-wetting-intermediate data)
    arrayView1d< TableFunction::KernelWrapper const > m_drainageRelPermKernelWrappers;

    /// Imbibition kernel wrappers for relative permeabilities in the following order:
    ///  0- wetting-phase
    ///  1- non-wetting-phase
    arrayView1d< TableFunction::KernelWrapper const > m_imbibitionRelPermKernelWrappers;

    /// Parameter a introduced by Jerauld in the Land model
    real64 const m_jerauldParam_a;

    /// Parameter b introduced by Jerauld in the Land model
    real64 const m_jerauldParam_b;

    /// Curvature parameter introduced for wetting phase hysteresis in Killough
    real64 const m_killoughCurvatureParam;

    /// Flag to specify whether the phase has hysteresis or not (deduced from table input)
    arrayView1d< integer const > m_phaseHasHysteresis;

    /// Trapping parameter from the Land model (typically called C)
    arrayView1d< real64 const > m_landParam;

    /// Minimum volume fraction for each phase in drainage (deduced from the drainage table)
    arrayView1d< real64 const > m_drainagePhaseMinVolFraction;

    /// Minimum volume fraction for each phase in imbibition (deduced from the imbibition table)
    arrayView1d< real64 const > m_imbibitionPhaseMinVolFraction;

    /// Maximum volume fraction for each phase
    arrayView1d< real64 const > m_drainagePhaseMaxVolFraction;

    /// Maximum volume fraction for each phase
    arrayView1d< real64 const > m_imbibitionPhaseMaxVolFraction;

    /// Relperm endpoint for each phase in drainage (deduced from the drainage table)
    arrayView1d< real64 const > m_drainagePhaseRelPermEndPoint;

    /// Relperm endpoint for each phase in imbibition (deduced from the imbibition table)
    arrayView1d< real64 const > m_imbibitionPhaseRelPermEndPoint;

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

  virtual void saveConvergedPhaseVolFractionState( arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFraction ) const override;

  struct viewKeyStruct : RelativePermeabilityBase::viewKeyStruct
  {
    static constexpr char const * drainageRelPermKernelWrappersString() { return "drainageRelPermWrappers"; }
    static constexpr char const * imbibitionRelPermKernelWrappersString() { return "imbibitionRelPermWrappers"; }

    static constexpr char const * phaseHasHysteresisString() { return "phaseHasHysteresis"; }
    static constexpr char const * jerauldParameterAString() { return "jerauldParameterA"; }
    static constexpr char const * jerauldParameterBString() { return "jerauldParameterB"; }
    static constexpr char const * killoughCurvatureParameterString() { return "killoughCurvatureParameter"; }
    static constexpr char const * landParameterString() { return "landParameter"; }

    static constexpr char const * drainagePhaseRelPermEndPointString() { return "drainagePhaseRelPermEndPoint"; }
    static constexpr char const * imbibitionPhaseRelPermEndPointString() { return "imbibitionPhaseRelPermEndPoint"; }
    static constexpr char const * drainagePhaseMinVolumeFractionString() { return "drainagePhaseMinVolumeFraction"; }
    static constexpr char const * imbibitionPhaseMinVolumeFractionString() { return "imbibitionPhaseMinVolumeFraction"; }
    static constexpr char const * drainagePhaseMaxVolumeFractionString() { return "drainagePhaseMaxVolumeFraction"; }
    static constexpr char const * imbibitionPhaseMaxVolumeFractionString() { return "imbibitionPhaseMaxVolumeFraction"; }

    static constexpr char const * drainageWettingNonWettingRelPermTableNamesString() { return "drainageWettingNonWettingRelPermTableNames"; }
    static constexpr char const * drainageWettingIntermediateRelPermTableNamesString() { return "drainageWettingIntermediateRelPermTableNames"; }
    static constexpr char const * drainageNonWettingIntermediateRelPermTableNamesString() { return "drainageNonWettingIntermediateRelPermTableNames"; }

    static constexpr char const * imbibitionWettingRelPermTableNameString() { return "imbibitionWettingRelPermTableName"; }
    static constexpr char const * imbibitionNonWettingRelPermTableNameString() { return "imbibitionNonWettingRelPermTableName"; }

  };

private:

  virtual void postProcessInput() override;

  virtual void initializePreSubGroups() override;

  virtual void resizeFields( localIndex const size, localIndex const numPts ) override;

  /**
   * @brief Create all the table kernel wrappers needed for the simulation (for all the phases present)
   */
  void createAllTableKernelWrappers();

  /**
   * @brief Check whether the drainage tables exist and validate all of them
   */
  void checkExistenceAndValidateDrainageRelPermTables();

  /**
   * @brief Check whether the imbibition tables exist and validate all of them
   */
  void checkExistenceAndValidateImbibitionRelPermTables();

  /**
   * @brief Check whether the table exists and validate it (increasing phase vol frac and rel perm, etc)
   * @param[in] relPermTableName the name of tje relative permeability table (kr vs s) for a given phase)
   * @param[out] phaseMinVolFrac the phase minimum volume fraction read from the table
   * @param[out] phaseMaxVolFrac the phase maximum volume fraction read from the table
   * @param[out] phaseRelPermEndPoint the end-point relative permeability
   */
  void checkExistenceAndValidateRelPermTable( string const & relPermTableName,
                                              real64 & phaseMinVolFrac,
                                              real64 & phaseMaxVolFrac,
                                              real64 & phaseRelPermEndPoint ) const;

  /**
   * @brief Compute the Land coefficient for the wetting and non-wetting phases
   */
  void computeLandCoefficient();

  // Table names

  /// Drainage relative permeability table names (one for each phase in the wetting-non-wetting pair)
  array1d< string > m_drainageWettingNonWettingRelPermTableNames;

  /// Drainage relative permeability table names (one for each phase in the wetting-intermediate pair)
  array1d< string > m_drainageWettingIntermediateRelPermTableNames;

  /// Drainage relative permeability table names (one for each phase in the non-wetting-intermediate pair)
  array1d< string > m_drainageNonWettingIntermediateRelPermTableNames;

  /// Imbibition relative permeability table name for the wetting phase
  string m_imbibitionWettingRelPermTableName;

  /// Imbibition relative permeability table name for the non-wetting phase
  string m_imbibitionNonWettingRelPermTableName;

  // Kernel wrappers

  /// Drainage kernel wrappers for relative permeabilities in the following order:
  /// Two-phase flow:
  ///  0- wetting-phase
  ///  1- non-wetting-phase
  /// Three-phase flow:
  ///  0- wetting-phase
  ///  1- intermediate phase (wetting-intermediate data)
  ///  2- non-wetting-phase
  ///  3- intermediate phase (non-wetting-intermediate data)
  array1d< TableFunction::KernelWrapper > m_drainageRelPermKernelWrappers;

  /// Imbibition kernel wrappers for relative permeabilities in the following order:
  ///  0- wetting-phase
  ///  1- non-wetting-phase
  array1d< TableFunction::KernelWrapper > m_imbibitionRelPermKernelWrappers;

  // Hysteresis parameters

  /// Parameter a introduced by Jerauld in the Land model
  real64 m_jerauldParam_a;

  /// Parameter b introduced by Jerauld in the Land model
  real64 m_jerauldParam_b;

  /// Curvature parameter in Killough wetting phase hysteresis (enpoints durvatures)
  real64 m_killoughCurvatureParam;

  /// Flag to specify whether the phase has hysteresis or not (deduced from table input)
  array1d< integer > m_phaseHasHysteresis;

  /// Trapping parameter from the Land model (typically called C)
  array1d< real64 > m_landParam;

  /// Minimum volume fraction for each phase in drainage (deduced from the drainage table)
  array1d< real64 > m_drainagePhaseMinVolFraction;

  /// Minimum volume fraction for each phase in imbibition (deduced from the imbibition table)
  array1d< real64 > m_imbibitionPhaseMinVolFraction;

  /// Relperm endpoint for each phase in drainage (deduced from the drainage table)
  array1d< real64 > m_drainagePhaseRelPermEndPoint;

  /// Relperm endpoint for each phase in imbibition (deduced from the imbibition table)
  array1d< real64 > m_imbibitionPhaseRelPermEndPoint;

  /// Maximum volume fraction for each phase
  array1d< real64 > m_drainagePhaseMaxVolFraction;

  /// Maximum volume fraction for each phase
  array1d< real64 > m_imbibitionPhaseMaxVolFraction;

  // Max historical saturations

  /// Minimum historical phase volume fraction for each phase
  array2d< real64, compflow::LAYOUT_PHASE > m_phaseMinHistoricalVolFraction;

  /// Maximum historical phase volume fraction for each phase
  array2d< real64, compflow::LAYOUT_PHASE > m_phaseMaxHistoricalVolFraction;

};

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
TableRelativePermeabilityHysteresis::KernelWrapper::
  computeDrainageRelPerm( TableFunction::KernelWrapper const & drainageRelPermKernelWrapper,
                          real64 const & phaseVolFraction,
                          real64 & phaseRelPerm,
                          real64 & dPhaseRelPerm_dPhaseVolFrac ) const
{
  phaseRelPerm =
    drainageRelPermKernelWrapper.compute( &phaseVolFraction,
                                          &dPhaseRelPerm_dPhaseVolFrac );
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
TableRelativePermeabilityHysteresis::KernelWrapper::
  computeImbibitionWettingRelPerm( TableFunction::KernelWrapper const & drainageRelPermKernelWrapper,
                                   TableFunction::KernelWrapper const & imbibitionRelPermKernelWrapper,
                                   real64 const & jerauldParam_a,
                                   real64 const & jerauldParam_b,
                                   real64 const & landParam,
                                   real64 const & phaseVolFraction,
                                   real64 const & phaseMinHistoricalVolFraction,
                                   real64 const & imbibitionPhaseMinWettingVolFraction,
                                   real64 const & drainagePhaseMaxVolFraction,
                                   real64 const & imbibitionPhaseMaxVolFraction,
                                   real64 const & drainageRelPermEndPoint,
                                   real64 const & imbibitionRelPermEndPoint,
                                   real64 & phaseRelPerm,
                                   real64 & dPhaseRelPerm_dPhaseVolFrac ) const
{

  // Step 0: preparing keypoints in the (S,kr) plan
  // if consistent, S should be equal to 1 - imbibitionPhaseMinVolNonWettingFraction for two-phase flow
  // (but wetting and nonwetting phase hysteresis are implemented in a decoupled fashion)
  real64 const S = phaseVolFraction;
  real64 const Smxi = imbibitionPhaseMaxVolFraction;
  real64 const Smxd = drainagePhaseMaxVolFraction;

  // Swc is the common end min endpoint saturation for wetting curves
  real64 const Swc = imbibitionPhaseMinWettingVolFraction;

  if( S <= Swc )
  {
    phaseRelPerm = 0.0;
    dPhaseRelPerm_dPhaseVolFrac = 0.0;
  }
  else if( S >= Smxd )
  {
    phaseRelPerm = drainageRelPermEndPoint;
    dPhaseRelPerm_dPhaseVolFrac = 0.0;
  }
  else
  {
    real64 const krwei = imbibitionRelPermEndPoint;
    real64 const krwedAtSmxi = drainageRelPermKernelWrapper.compute( &Smxi );

    // Step 1: Compute the new end point

    // Step 1.a: get the value at the max non-wetting residual value
    real64 const deltak = krwei - krwedAtSmxi;

    // Step 1.b: get the trapped from wetting data
    real64 const Shy = ( phaseMinHistoricalVolFraction > Swc ) ? phaseMinHistoricalVolFraction : Swc;
    real64 const A = 1 + jerauldParam_a * ( Shy - Swc );
    real64 const numerator = Shy - Smxd;
    real64 const denom = A + landParam * pow( ( Smxd - Shy ) / ( Smxd - Swc ), 1 + jerauldParam_b/landParam );
    real64 const Scrt = Smxd + numerator / denom;

    // Step 1.c: find the new endpoint
    // this is the saturation for the scanning curve endpoint
    real64 const krwedAtScrt = drainageRelPermKernelWrapper.compute( &Scrt );
    real64 const krwieStar = krwedAtScrt
                             + deltak * pow( ( Smxd - Scrt ) / LvArray::math::max( minScriMinusScrd, ( Smxd - Smxi ) ), m_killoughCurvatureParam );

    // Step 2: get the normalized value of saturation
    real64 const ratio = ( Smxi - Swc ) / ( Scrt - Shy );
    real64 const Snorm = Smxi - ( Scrt - S ) * ratio; // normalized saturation from equation 2.166
    real64 const dSnorm_dS =  ratio;
    real64 dkri_dSnorm = 0.0;
    real64 const krwiAtSnorm = imbibitionRelPermKernelWrapper.compute( &Snorm, &dkri_dSnorm );
    real64 const dkriAtSnorm_dS = dkri_dSnorm * dSnorm_dS;

    // Step 3: Get the final value at evaluated saturation
    real64 const krdAtShy = drainageRelPermKernelWrapper.compute( &Shy );
    real64 const imbibitionRelPermRatio = (krwieStar - krdAtShy) / krwei;

    phaseRelPerm = krdAtShy + krwiAtSnorm * imbibitionRelPermRatio;
    dPhaseRelPerm_dPhaseVolFrac = dkriAtSnorm_dS * imbibitionRelPermRatio;
  }
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
TableRelativePermeabilityHysteresis::KernelWrapper::
  computeImbibitionNonWettingRelPerm( TableFunction::KernelWrapper const & drainageRelPermKernelWrapper,
                                      TableFunction::KernelWrapper const & imbibitionRelPermKernelWrapper,
                                      real64 const & jerauldParam_a,
                                      real64 const & jerauldParam_b,
                                      real64 const & landParam,
                                      real64 const & phaseVolFraction,
                                      real64 const & phaseMaxHistoricalVolFraction,
                                      real64 const & drainagePhaseMinVolFraction,
                                      real64 const & imbibitionPhaseMinVolFraction,
                                      real64 const & drainagePhaseMaxVolFraction,
                                      real64 const & drainageRelPermEndPoint,
                                      real64 & phaseRelPerm,
                                      real64 & dPhaseRelPerm_dPhaseVolFrac ) const
{
  // note: for simplicity, the notations are taken from IX documentation (although this breaks our phaseVolFrac naming convention)

  // Step 1: for a given value of the max historical saturation, Shy, compute the trapped critical saturation, Scrt,
  //         using Land's method. The calculation includes the modifications from Jerauld. This is equation 2.162 from
  //         the IX technical description.
  real64 const S = phaseVolFraction;
  real64 const Scri = imbibitionPhaseMinVolFraction;
  real64 const Scrd = drainagePhaseMinVolFraction;
  real64 const Smx = drainagePhaseMaxVolFraction;
  real64 const Shy = phaseMaxHistoricalVolFraction < Smx ? phaseMaxHistoricalVolFraction : Smx; // to make sure that Shy < Smax
  real64 const A = 1 + jerauldParam_a * ( Smx - Shy );
  real64 const numerator = Shy - Scrd;
  real64 const denom = A + landParam * pow( ( Shy - Scrd ) / ( Smx - Scrd ), 1 + jerauldParam_b / landParam );
  real64 const Scrt = Scrd + numerator / denom; // trapped critical saturation from equation 2.162

  if( S <= Scrt )  // S is below the trapped critical saturation, so the relperm is zero
  {
    phaseRelPerm = 0.0;
    dPhaseRelPerm_dPhaseVolFrac = 0.0;
  }
  else if( S >= Smx ) // S is above the max saturation, so we just skip the rest and set the relperm to the endpoint
  {
    phaseRelPerm = drainageRelPermEndPoint;
    dPhaseRelPerm_dPhaseVolFrac = 0.0;
  }
  else
  {
    // Step 2: compute the normalized saturation, S_norm, at which the imbibition relperm curve will be evaluated.
    //         This is equation 2.166 from the IX technical description.
    real64 const ratio = ( Smx - Scri ) / ( Shy - Scrt );
    real64 const Snorm = Scri + ( S - Scrt ) * ratio; // normalized saturation from equation 2.166
    real64 const dSnorm_dS = ratio;

    // Step 3: evaluate the imbibition relperm, kri(Snorm), at the normalized saturation, Snorm.
    real64 dkri_dSnorm = 0;
    real64 const kriAtSnorm = imbibitionRelPermKernelWrapper.compute( &Snorm, &dkri_dSnorm );
    real64 const dkriAtSnorm_dS = dkri_dSnorm * dSnorm_dS;

    // Step 4: evaluate the drainage relperm, krd(Shy), at the max hystorical saturation, Shy.
    real64 const krdAtShy = drainageRelPermKernelWrapper.compute( &Shy );

    // Step 5: evaluate the drainage relperm, krd(Smx), at the max drainage saturation, Smx.
    real64 const krdAtSmx = drainageRelPermEndPoint;

    // Step 6: apply the formula blending drainage and imbibition relperms from the Killough model.
    //         This equation 2.165 from the IX technical description.
    real64 const drainageRelPermRatio = krdAtShy / krdAtSmx;
    phaseRelPerm = kriAtSnorm * drainageRelPermRatio;
    dPhaseRelPerm_dPhaseVolFrac = dkriAtSnorm_dS * drainageRelPermRatio;
  }
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
TableRelativePermeabilityHysteresis::KernelWrapper::
  computeTwoPhase( integer const ipWetting,
                   integer const ipNonWetting,
                   arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
                   arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseMaxHistoricalVolFraction,
                   arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseMinHistoricalVolFraction,
                   arraySlice1d< real64, relperm::USD_RELPERM - 2 > const & phaseRelPerm,
                   arraySlice2d< real64, relperm::USD_RELPERM_DS - 2 > const & dPhaseRelPerm_dPhaseVolFrac ) const
{
  using TPT = TableRelativePermeabilityHysteresis::TwoPhasePairPhaseType;
  using IPT = TableRelativePermeabilityHysteresis::ImbibitionPhasePairPhaseType;

  // ---------- wetting rel perm
  if( !m_phaseHasHysteresis[IPT::WETTING] ||
      phaseVolFraction[ipWetting] <= phaseMinHistoricalVolFraction[ipWetting] + flowReversalBuffer )
  {
    computeDrainageRelPerm( m_drainageRelPermKernelWrappers[TPT::WETTING],
                            phaseVolFraction[ipWetting],
                            phaseRelPerm[ipWetting],
                            dPhaseRelPerm_dPhaseVolFrac[ipWetting][ipWetting] );
  }
  else
  {
    computeImbibitionWettingRelPerm( m_drainageRelPermKernelWrappers[TPT::WETTING],
                                     m_imbibitionRelPermKernelWrappers[IPT::WETTING],
                                     m_jerauldParam_a,
                                     m_jerauldParam_b,
                                     m_landParam[IPT::WETTING],
                                     phaseVolFraction[ipWetting],
                                     phaseMinHistoricalVolFraction[ipWetting],
                                     m_imbibitionPhaseMinVolFraction[IPT::WETTING],
                                     m_drainagePhaseMaxVolFraction[ipWetting],
                                     m_imbibitionPhaseMaxVolFraction[IPT::WETTING],
                                     m_drainagePhaseRelPermEndPoint[ipWetting],
                                     m_imbibitionPhaseRelPermEndPoint[IPT::WETTING],
                                     phaseRelPerm[ipWetting],
                                     dPhaseRelPerm_dPhaseVolFrac[ipWetting][ipWetting] );
  }

  // --------- non-wetting rel perm
  if( !m_phaseHasHysteresis[IPT::NONWETTING] ||
      phaseVolFraction[ipNonWetting] >= phaseMaxHistoricalVolFraction[ipNonWetting] - flowReversalBuffer )
  {
    computeDrainageRelPerm( m_drainageRelPermKernelWrappers[TPT::NONWETTING],
                            phaseVolFraction[ipNonWetting],
                            phaseRelPerm[ipNonWetting],
                            dPhaseRelPerm_dPhaseVolFrac[ipNonWetting][ipNonWetting] );
  }
  else
  {
    computeImbibitionNonWettingRelPerm( m_drainageRelPermKernelWrappers[TPT::NONWETTING],
                                        m_imbibitionRelPermKernelWrappers[IPT::NONWETTING],
                                        m_jerauldParam_a,
                                        m_jerauldParam_b,
                                        m_landParam[IPT::NONWETTING],
                                        phaseVolFraction[ipNonWetting],
                                        phaseMaxHistoricalVolFraction[ipNonWetting],
                                        m_drainagePhaseMinVolFraction[ipNonWetting],
                                        m_imbibitionPhaseMinVolFraction[IPT::NONWETTING],
                                        m_drainagePhaseMaxVolFraction[ipNonWetting],
                                        m_drainagePhaseRelPermEndPoint[ipNonWetting],
                                        phaseRelPerm[ipNonWetting],
                                        dPhaseRelPerm_dPhaseVolFrac[ipNonWetting][ipNonWetting] );
  }
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
TableRelativePermeabilityHysteresis::KernelWrapper::
  computeThreePhase( integer const ipWetting,
                     integer const ipInter,
                     integer const ipNonWetting,
                     arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
                     arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseMaxHistoricalVolFraction,
                     arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseMinHistoricalVolFraction,
                     arraySlice1d< real64, relperm::USD_RELPERM - 2 > const & phaseRelPerm,
                     arraySlice2d< real64, relperm::USD_RELPERM_DS - 2 > const & dPhaseRelPerm_dPhaseVolFrac ) const
{
  real64 interRelPerm_wi = 0; // oil rel perm using two-phase gas-oil data
  real64 dInterRelPerm_wi_dInterVolFrac = 0; // derivative w.r.t to So
  real64 interRelPerm_nwi = 0; // oil rel perm using two-phase gas-oil data
  real64 dInterRelPerm_nwi_dInterVolFrac = 0; // derivative w.r.t to So

  using TPT = TableRelativePermeabilityHysteresis::ThreePhasePairPhaseType;
  using IPT = TableRelativePermeabilityHysteresis::ImbibitionPhasePairPhaseType;

  // 1) Wetting and intermediate phase relative permeabilities using two-phase wetting-intermediate data

  // ---------- wetting rel perm
  if( !m_phaseHasHysteresis[IPT::WETTING] ||
      phaseVolFraction[ipWetting] <= phaseMinHistoricalVolFraction[ipWetting] + flowReversalBuffer )
  {
    computeDrainageRelPerm( m_drainageRelPermKernelWrappers[TPT::WETTING],
                            phaseVolFraction[ipWetting],
                            phaseRelPerm[ipWetting],
                            dPhaseRelPerm_dPhaseVolFrac[ipWetting][ipWetting] );
  }
  else
  {
    computeImbibitionWettingRelPerm( m_drainageRelPermKernelWrappers[TPT::WETTING],
                                     m_imbibitionRelPermKernelWrappers[IPT::WETTING],
                                     m_jerauldParam_a,
                                     m_jerauldParam_b,
                                     m_landParam[IPT::WETTING],
                                     phaseVolFraction[ipWetting],
                                     phaseMinHistoricalVolFraction[ipWetting],
                                     m_imbibitionPhaseMinVolFraction[IPT::WETTING],
                                     m_drainagePhaseMaxVolFraction[ipWetting],
                                     m_imbibitionPhaseMaxVolFraction[IPT::WETTING],
                                     m_drainagePhaseRelPermEndPoint[ipWetting],
                                     m_imbibitionPhaseRelPermEndPoint[IPT::WETTING],
                                     phaseRelPerm[ipWetting],
                                     dPhaseRelPerm_dPhaseVolFrac[ipWetting][ipWetting] );
  }

  // ---------- intermediate rel perm (ALWAYS DRAINAGE!)
  interRelPerm_wi =
    m_drainageRelPermKernelWrappers[TPT::INTERMEDIATE_WETTING].compute( &(phaseVolFraction)[ipInter],
                                                                        &dInterRelPerm_wi_dInterVolFrac );


  // 2) Non-wetting and intermediate phase relative permeabilities using two-phase non-wetting-intermediate data

  // ---------- non-wetting rel perm
  if( !m_phaseHasHysteresis[IPT::NONWETTING] ||
      phaseVolFraction[ipNonWetting] >= phaseMaxHistoricalVolFraction[ipNonWetting] - flowReversalBuffer )
  {
    computeDrainageRelPerm( m_drainageRelPermKernelWrappers[TPT::NONWETTING],
                            phaseVolFraction[ipNonWetting],
                            phaseRelPerm[ipNonWetting],
                            dPhaseRelPerm_dPhaseVolFrac[ipNonWetting][ipNonWetting] );
  }
  else
  {
    computeImbibitionNonWettingRelPerm( m_drainageRelPermKernelWrappers[TPT::NONWETTING],
                                        m_imbibitionRelPermKernelWrappers[IPT::NONWETTING],
                                        m_jerauldParam_a,
                                        m_jerauldParam_b,
                                        m_landParam[IPT::NONWETTING],
                                        phaseVolFraction[ipNonWetting],
                                        phaseMaxHistoricalVolFraction[ipNonWetting],
                                        m_drainagePhaseMinVolFraction[ipNonWetting],
                                        m_imbibitionPhaseMinVolFraction[IPT::NONWETTING],
                                        m_drainagePhaseMaxVolFraction[ipNonWetting],
                                        m_drainagePhaseRelPermEndPoint[ipNonWetting],
                                        phaseRelPerm[ipNonWetting],
                                        dPhaseRelPerm_dPhaseVolFrac[ipNonWetting][ipNonWetting] );
  }

  // ---------- intermediate rel perm (ALWAYS DRAINAGE!)
  interRelPerm_nwi =
    m_drainageRelPermKernelWrappers[TPT::INTERMEDIATE_NONWETTING].compute( &(phaseVolFraction)[ipInter],
                                                                           &dInterRelPerm_nwi_dInterVolFrac );

  // 3) Compute the "three-phase" oil relperm

  // use saturation-weighted interpolation
  real64 const shiftedWettingVolFrac = (phaseVolFraction[ipWetting] - m_drainagePhaseMinVolFraction[ipWetting]);

  relpermInterpolators::Baker::compute( shiftedWettingVolFrac,
                                        phaseVolFraction[ipNonWetting],
                                        m_phaseOrder,
                                        interRelPerm_wi,
                                        dInterRelPerm_wi_dInterVolFrac,
                                        interRelPerm_nwi,
                                        dInterRelPerm_nwi_dInterVolFrac,
                                        phaseRelPerm[ipInter],
                                        dPhaseRelPerm_dPhaseVolFrac[ipInter] );
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
TableRelativePermeabilityHysteresis::KernelWrapper::
  compute( arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
           arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseMaxHistoricalVolFraction,
           arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseMinHistoricalVolFraction,
           arraySlice1d< real64, relperm::USD_RELPERM - 2 > const & phaseRelPerm,
           arraySlice2d< real64, relperm::USD_RELPERM_DS - 2 > const & dPhaseRelPerm_dPhaseVolFrac ) const
{
  LvArray::forValuesInSlice( dPhaseRelPerm_dPhaseVolFrac, []( real64 & val ){ val = 0.0; } );

  using PT = RelativePermeabilityBase::PhaseType;
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
                       phaseRelPerm,
                       dPhaseRelPerm_dPhaseVolFrac );

  }
  else if( ipWater < 0 )
  {
    computeTwoPhase( ipOil, // wetting
                     ipGas, // non-wetting
                     phaseVolFraction,
                     phaseMaxHistoricalVolFraction,
                     phaseMinHistoricalVolFraction,
                     phaseRelPerm,
                     dPhaseRelPerm_dPhaseVolFrac );
  }
  else if( ipOil < 0 )
  {
    computeTwoPhase( ipWater, // wetting
                     ipGas,   // non-wetting
                     phaseVolFraction,
                     phaseMaxHistoricalVolFraction,
                     phaseMinHistoricalVolFraction,
                     phaseRelPerm,
                     dPhaseRelPerm_dPhaseVolFrac );
  }
  else if( ipGas < 0 )
  {
    computeTwoPhase( ipWater, // wetting
                     ipOil,   // non-wetting
                     phaseVolFraction,
                     phaseMaxHistoricalVolFraction,
                     phaseMinHistoricalVolFraction,
                     phaseRelPerm,
                     dPhaseRelPerm_dPhaseVolFrac );
  }
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
TableRelativePermeabilityHysteresis::KernelWrapper::
  update( localIndex const k,
          localIndex const q,
          arraySlice1d< geosx::real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction ) const
{
#if defined(GEOSX_USE_HIP) && defined(GEOSX_DEVICE_COMPILE) && defined(NDEBUG)
  GEOSX_ERROR("Can't compile this kernel with HIP yet.");
#else
  compute( phaseVolFraction,
           m_phaseMaxHistoricalVolFraction[k],
           m_phaseMinHistoricalVolFraction[k],
           m_phaseRelPerm[k][q],
           m_dPhaseRelPerm_dPhaseVolFrac[k][q] );
#endif
}

} // namespace constitutive

} // namespace geosx

#endif //GEOSX_CONSTITUTIVE_TABLERELATIVEPERMEABILITYHYSTERESIS_HPP
