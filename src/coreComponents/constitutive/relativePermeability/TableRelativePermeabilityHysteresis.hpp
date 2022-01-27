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

    KernelWrapper( arrayView1d< TableFunction::KernelWrapper const > const & drainageRelPermKernelWrappers,
                   arrayView1d< TableFunction::KernelWrapper const > const & imbibitionRelPermKernelWrappers,
                   real64 const & jerauldParam_a,
                   real64 const & jerauldParam_b,
                   real64 const & alphaParam_2,
                   arrayView1d< integer const > const & phaseHasHysteresis,
                   arrayView1d< real64 const > const & landParam,
                   arrayView1d< real64 const > const & drainageMinPhaseVolFraction,
                   arrayView1d< real64 const > const & imbibitionMinPhaseVolFraction,
                   arrayView1d< real64 const > const & drainageMaxPhaseVolFraction,
                   arrayView1d< real64 const > const & imbibitionMaxPhaseVolFraction,
                   arrayView1d< integer const > const & phaseTypes,
                   arrayView1d< integer const > const & phaseOrder,
                   arrayView2d< real64 const, compflow::USD_PHASE > const & phaseMaxHistoricalVolFraction,
                   arrayView2d< real64 const, compflow::USD_PHASE > const & phaseMinHistoricalVolFraction,
                   arrayView1d< real64 const > const & drainageKrEndPoint,
                   arrayView1d< real64 const > const & imbibitionKrEndPoint,
                   arrayView3d< real64, relperm::USD_RELPERM > const & phaseRelPerm,
                   arrayView4d< real64, relperm::USD_RELPERM_DS > const & dPhaseRelPerm_dPhaseVolFrac );

    GEOSX_HOST_DEVICE
    void computeDrainageRelPerm( TableFunction::KernelWrapper const & drainageRelPermKernelWrapper,
                                 real64 const & phaseVolFraction,
                                 real64 & phaseRelPerm,
                                 real64 & dPhaseRelPerm_dPhaseVolFrac ) const;

    GEOSX_HOST_DEVICE
    void computeImbibitionWettingRelPerm( TableFunction::KernelWrapper const & drainageRelPermKernelWrapper,
                                          TableFunction::KernelWrapper const & imbibitionRelPermKernelWrapper,
                                          real64 const & jerauldParam_a,
                                          real64 const & jerauldParam_b,
                                          real64 const & landParam,
                                          real64 const & imbibitionKrEndPoint,
                                          real64 const & phaseVolFraction,
                                          real64 const & phaseMinHistoricalVolFraction,
                                          real64 const & imbibitionMinPhaseWettingVolFraction,
                                          real64 const & drainagePhaseMaxVolFraction,
                                          real64 const & imbibitionPhaseMaxVolFraction,
                                          real64 & phaseRelPerm,
                                          real64 & dPhaseRelPerm_dPhaseVolFrac ) const;

    GEOSX_HOST_DEVICE
    void computeImbibitionNonWettingRelPerm( TableFunction::KernelWrapper const & drainageRelPermKernelWrapper,
                                             TableFunction::KernelWrapper const & imbibitionRelPermKernelWrapper,
                                             real64 const & phaseVolFraction,
                                             real64 const & phaseMaxHistoricalVolFraction,
                                             real64 const & drainageMinPhaseVolFraction,
                                             real64 const & imbibitionMinPhaseVolFraction,
                                             real64 const & maxPhaseVolFraction,
                                             real64 const & jerauldParam_a,
                                             real64 const & jerauldParam_b,
                                             real64 const & landParam,
                                             real64 & phaseRelPerm,
                                             real64 & dPhaseRelPerm_dPhaseVolFrac ) const;

    GEOSX_HOST_DEVICE
    void computeTwoPhase( localIndex const ipWetting,
                          localIndex const ipNonWetting,
                          arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
                          arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseMaxHistoricalVolFraction,
                          arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseMinHistoricalVolFraction,
                          arraySlice1d< real64, relperm::USD_RELPERM - 2 > const & phaseRelPerm,
                          arraySlice2d< real64, relperm::USD_RELPERM_DS - 2 > const & dPhaseRelPerm_dPhaseVolFrac ) const;

    GEOSX_HOST_DEVICE
    void computeThreePhase( localIndex const ipWetting,
                            localIndex const ipInter,
                            localIndex const ipNonWetting,
                            arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
                            arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseMaxHistoricalVolFraction,
                            arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseMinHistoricalVolFraction,
                            arraySlice1d< real64, relperm::USD_RELPERM - 2 > const & phaseRelPerm,
                            arraySlice2d< real64, relperm::USD_RELPERM_DS - 2 > const & dPhaseRelPerm_dPhaseVolFrac ) const;

    GEOSX_HOST_DEVICE
    void compute( arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
                  arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseMaxHistoricalVolFraction,
                  arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseMinHistoricalVolFraction,
                  arraySlice1d< real64, relperm::USD_RELPERM - 2 > const & phaseRelPerm,
                  arraySlice2d< real64, relperm::USD_RELPERM_DS - 2 > const & dPhaseRelPerm_dPhaseVolFrac ) const;

    GEOSX_HOST_DEVICE
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

    /// Parameter alpha_2 introduced for wetting phase hysteresis in Killough
    real64 const m_alphaParam_2;

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

    /// Maximum historical phase volume fraction for each phase
    arrayView2d< real64 const, compflow::USD_PHASE > m_phaseMaxHistoricalVolFraction;

    /// Minimum historical phase volume fraction for each phase
    arrayView2d< real64 const, compflow::USD_PHASE > m_phaseMinHistoricalVolFraction;

    /// Relperm endpoint for each phase in drainage (deduced from the drainage table)
    arrayView1d< real64 const > m_drainagePhaseRelPermEndPoint;

    /// Relperm endpoint for each phase in imbibition (deduced from the imbibition table)
    arrayView1d< real64 const > m_imbibitionPhaseRelPermEndPoint;

  };

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

  virtual void initializeState( arrayView2d< real64 const, compflow::USD_PHASE > const & initialPhaseVolFraction ) const override;

  virtual void saveConvergedPhaseVolFraction( arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFraction ) const override;

  struct viewKeyStruct : RelativePermeabilityBase::viewKeyStruct
  {
    static constexpr char const * drainageRelPermKernelWrappersString() { return "drainageRelPermWrappers"; }
    static constexpr char const * imbibitionRelPermKernelWrappersString() { return "imbibitionRelPermWrappers"; }

    static constexpr char const * phaseMaxHistoricalVolumeFractionString() { return "phaseMaxHistoricalVolumeFraction"; }
    static constexpr char const * phaseMinHistoricalVolumeFractionString() { return "phaseMinHistoricalVolumeFraction"; }

    static constexpr char const * phaseHasHysteresisString() { return "phaseHasHysteresis"; }
    static constexpr char const * jerauldParameterAString() { return "jerauldParameterA"; }
    static constexpr char const * jerauldParameterBString() { return "jerauldParameterB"; }
    static constexpr char const * alphaParameter2String() { return "alphaParameter2"; }
    static constexpr char const * landParameterString() { return "landParameter"; }

    static constexpr char const * imbibitionPhaseMinVolumeFractionString() { return "imbibitionPhaseMinVolumeFraction"; }
    static constexpr char const * drainagePhaseRelPermEndPointString() { return "drainagePhaseRelPermEndPoint"; }
    static constexpr char const * imbibitionPhaseRelPermEndPointString() { return "imbibitionPhaseRelPermEndPoint"; }
    static constexpr char const * drainagePhaseMinVolumeFractionString() { return "drainagePhaseMinVolumeFraction"; }
    static constexpr char const * drainagePhaseMaxVolumeFractionString() { return "drainagePhaseMaxVolumeFraction"; }
    static constexpr char const * imbibitionPhaseMaxVolumeFractionString() { return "imbibitionPhaseMaxVolumeFraction"; }

    static constexpr char const * drainageWettingNonWettingRelPermTableNamesString() { return "drainageWettingNonWettingRelPermTableNames"; }
    static constexpr char const * drainageWettingIntermediateRelPermTableNamesString() { return "drainageWettingIntermediateRelPermTableNames"; }
    static constexpr char const * drainageNonWettingIntermediateRelPermTableNamesString() { return "drainageNonWettingIntermediateRelPermTableNames"; }

    static constexpr char const * imbibitionWettingRelPermTableNameString() { return "imbibibitionWettingRelPermTableName"; }
    static constexpr char const * imbibitionNonWettingRelPermTableNameString() { return "imbibitionNonWettingImbibitionRelPermTableName"; }

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
   * @brief Validate the relative permeability table provided in input (increasing phase vol frac and rel perm, etc)
   * @param[in] relPermTable the relative permeability table (kr vs s) for a given phase)
   * @param[out] phaseMinVolFrac the phase minimum volume fraction read from the table
   * @param[out] phaseMaxVolFrac the phase maximum volume fraction read from the table
   * @param[out] phaseRelPermEndPoint the end-point relative permeability
   */
  void validateRelPermTable( TableFunction const & relPermTable,
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

  /// Parameter alpha2 in Killough wetting phase hysteresis (enpoints durvatures)
  real64 m_alphaParam_2;

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

  /// Maximum historical phase volume fraction for each phase
  array2d< real64, compflow::LAYOUT_PHASE > m_phaseMaxHistoricalVolFraction;

  /// Minimum historical phase volume fraction for each phase
  array2d< real64, compflow::LAYOUT_PHASE > m_phaseMinHistoricalVolFraction;

};

GEOSX_HOST_DEVICE
inline void
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
inline void
TableRelativePermeabilityHysteresis::KernelWrapper::
  computeImbibitionWettingRelPerm( TableFunction::KernelWrapper const & drainageRelPermKernelWrapper,
                                   TableFunction::KernelWrapper const & imbibitionRelPermKernelWrapper,
                                   real64 const & jerauldParam_a,
                                   real64 const & jerauldParam_b,
                                   real64 const & landParam,
                                   real64 const & imbibitionKrEndPoint,
                                   real64 const & phaseVolFraction,
                                   real64 const & phaseMinHistoricalVolFraction,
                                   real64 const & imbibitionMinPhaseWettingVolFraction,
                                   real64 const & drainagePhaseMaxVolFraction,
                                   real64 const & imbibitionPhaseMaxVolFraction,
                                   real64 & phaseRelPerm,
                                   real64 & dPhaseRelPerm_dPhaseVolFrac ) const
{

//0. preparing keypoints in the (S,kr) plan
//real64 const Scri = imbibitionPhaseMinVolNonWettingFraction;
//if consistify should be equal to 1. - imbibitionPhaseMinVolNonWettingFraction
//but wetting and nonwetting are implemented decoupled
  real64 const Smxi = imbibitionPhaseMaxVolFraction;
  real64 const Smxd = drainagePhaseMaxVolFraction;
//  real64 const Scri = 1. - Smxi;
//  real64 const Scrd = 1. - Smxd;
  // Swc is the common end min endpoint sat for wetting curves
  real64 const Swc = imbibitionMinPhaseWettingVolFraction;

  real64 const krwei = imbibitionKrEndPoint;
  real64 const krwedAtSmxi = drainageRelPermKernelWrapper.compute( &Smxi );

  //1. Compute the new end point
  //a. get the value at the max NW residual value
  real64 const deltak = krwei - krwedAtSmxi;

  //b. get the trapped from wetting data
  real64 const Shy = phaseMinHistoricalVolFraction;
  real64 const A = 1 + jerauldParam_a * ( Shy - Swc );
  real64 const numerator = (Shy - Smxd);
  real64 const denom = A + landParam * pow( (Smxd-Shy)/(Smxd-Swc), 1 + jerauldParam_b/landParam );
  real64 const Scrt = Smxd + numerator / denom;


  //c. Finding the new endpoint
  //this is the saturation for the scanning curve endpoint
  real64 const krwedAtScrt = drainageRelPermKernelWrapper.compute( &Scrt );
  real64 const krwieStar = krwedAtScrt + deltak * pow( (Smxd - Scrt) / (Smxd - Smxi), m_alphaParam_2 );

  //2. Get the normalized value of sat
  real64 const S = phaseVolFraction;
  real64 const ratio = ( Smxi - Swc ) / ( Scrt - Shy );
  real64 const Snorm = Smxi - ( Scrt - S ) * ratio; // normalized saturation from equation 2.166
  real64 const dSnorm_dS =  -ratio;
  real64 dkri_dSnorm = 0.0;
  real64 const krwiAtSnorm = imbibitionRelPermKernelWrapper.compute( &Snorm, &dkri_dSnorm );
  real64 const dkriAtSnorm_dS = dkri_dSnorm * dSnorm_dS;

  //3. Get the final value at evaluated sat
  real64 const krdAtShy = drainageRelPermKernelWrapper.compute( &Shy );
  real64 const imbibitionRelPermRatio = (krwieStar - krdAtShy) / krwei;

  phaseRelPerm = krdAtShy + krwiAtSnorm * imbibitionRelPermRatio;
  dPhaseRelPerm_dPhaseVolFrac = dkriAtSnorm_dS * imbibitionRelPermRatio;
}

GEOSX_HOST_DEVICE
inline void
TableRelativePermeabilityHysteresis::KernelWrapper::
  computeImbibitionNonWettingRelPerm( TableFunction::KernelWrapper const & drainageRelPermKernelWrapper,
                                      TableFunction::KernelWrapper const & imbibitionRelPermKernelWrapper,
                                      real64 const & jerauldParam_a,
                                      real64 const & jerauldParam_b,
                                      real64 const & landParam,
                                      real64 const & imbibitionPhaseMinVolFraction,
                                      real64 const & drainagePhaseMinVolFraction,
                                      real64 const & phaseMaxVolFraction,
                                      real64 const & phaseVolFraction,
                                      real64 const & phaseMaxHistoricalVolFraction,
                                      real64 & phaseRelPerm,
                                      real64 & dPhaseRelPerm_dPhaseVolFrac ) const
{
  // note: for simplicity, the notations are taken from IX documentation (although this breaks our phaseVolFrac naming convention)

  // Step 1: for a given value of the max historical saturation, Shy, compute the trapped critical saturation, Scrt,
  //         using Land's method. The calculation includes the modifications from Jerauld. This is equation 2.162 from
  //         the IX technical description.
  real64 const Shy = phaseMaxHistoricalVolFraction;
  real64 const Scrd = drainagePhaseMinVolFraction;
  real64 const Smx = phaseMaxVolFraction;
  real64 const A = 1 + jerauldParam_a * ( Smx - Shy );
  real64 const numerator = Shy - Scrd;
  real64 const denom = A + landParam * pow( ( Shy - Scrd ) / ( Smx - Scrd ), 1 + jerauldParam_b / landParam );
  real64 const Scrt = Scrd + numerator / denom; // trapped critical saturation from equation 2.162

  // Step 2: compute the normalized saturation, S_norm, at which the imbibition relperm curve will be evaluated.
  //         This is equation 2.166 from the IX technical description.
  real64 const S = phaseVolFraction;
  real64 const Scri = imbibitionPhaseMinVolFraction;
  real64 const ratio = ( Smx - Scri ) / ( Shy - Scrt );
  real64 const Snorm = Scri + ( S - Scrt ) * ratio; // normalized saturation from equation 2.166
  real64 const dSnorm_dS = ratio;


  GEOSX_ERROR_IF( Scrt > Scri,
                  GEOSX_FMT( string( "For non-wetting phase hysteresis, trapping saturation Snrt: {} " ) +
                             string( " should be less or equal to the trapping sat for imbibition Sncri: {}" ),
                             Scrt, Scri )
                  );


  // Step 3: evaluate the imbibition relperm, kri(Snorm), at the normalized saturation, Snorm.
  real64 dkri_dSnorm = 0;
  real64 const kriAtSnorm = imbibitionRelPermKernelWrapper.compute( &Snorm, &dkri_dSnorm );
  real64 const dkriAtSnorm_dS = dkri_dSnorm * dSnorm_dS;

  // Step 4: evaluate the drainage relperm, krd(Shy), at the max hystorical saturation, Shy.
  real64 const krdAtShy = drainageRelPermKernelWrapper.compute( &Shy );

  // Step 5: evaluate the drainage relperm, krd(Smx), at the max drainage saturation, Smx.
  real64 const krdAtSmx = drainageRelPermKernelWrapper.compute( &Smx );

  // Step 6: apply the formula blending drainage and imbibition relperms from the Killough model.
  //         This equation 2.165 from the IX technical description.
  real64 const drainageRelPermRatio = krdAtShy / krdAtSmx;
  phaseRelPerm = kriAtSnorm * drainageRelPermRatio;
  dPhaseRelPerm_dPhaseVolFrac = dkriAtSnorm_dS * drainageRelPermRatio;
}

GEOSX_HOST_DEVICE
inline void
TableRelativePermeabilityHysteresis::KernelWrapper::
  computeTwoPhase( localIndex const ipWetting,
                   localIndex const ipNonWetting,
                   arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction,
                   arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseMaxHistoricalVolFraction,
                   arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseMinHistoricalVolFraction,
                   arraySlice1d< real64, relperm::USD_RELPERM - 2 > const & phaseRelPerm,
                   arraySlice2d< real64, relperm::USD_RELPERM_DS - 2 > const & dPhaseRelPerm_dPhaseVolFrac ) const
{
  using TPT = TableRelativePermeabilityHysteresis::TwoPhasePairPhaseType;
  using IPT = TableRelativePermeabilityHysteresis::ImbibitionPhasePairPhaseType;

  // ---------- wetting rel perm
  if( !m_phaseHasHysteresis[IPT::WETTING] || phaseVolFraction[ipWetting] <= phaseMinHistoricalVolFraction[ipWetting] )
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
                                     m_imbibitionPhaseRelPermEndPoint[IPT::WETTING],
                                     phaseVolFraction[ipWetting],
                                     phaseMinHistoricalVolFraction[ipWetting],
                                     m_imbibitionPhaseMinVolFraction[IPT::WETTING],
                                     m_drainagePhaseMaxVolFraction[TPT::WETTING],
                                     m_imbibitionPhaseMaxVolFraction[IPT::WETTING],
                                     phaseRelPerm[ipWetting],
                                     dPhaseRelPerm_dPhaseVolFrac[ipWetting][ipWetting] );
  }

  // --------- non-wetting rel perm
  if( !m_phaseHasHysteresis[IPT::NONWETTING] || phaseVolFraction[ipNonWetting] >= phaseMaxHistoricalVolFraction[ipNonWetting] )
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
                                        m_imbibitionPhaseMinVolFraction[IPT::NONWETTING],
                                        m_drainagePhaseMinVolFraction[ipNonWetting],
                                        m_drainagePhaseMaxVolFraction[ipNonWetting],
                                        phaseVolFraction[ipNonWetting],
                                        phaseMaxHistoricalVolFraction[ipNonWetting],
                                        phaseRelPerm[ipNonWetting],
                                        dPhaseRelPerm_dPhaseVolFrac[ipNonWetting][ipNonWetting] );
  }
}

GEOSX_HOST_DEVICE
inline void
TableRelativePermeabilityHysteresis::KernelWrapper::
  computeThreePhase( localIndex const ipWetting,
                     localIndex const ipInter,
                     localIndex const ipNonWetting,
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
  if( phaseVolFraction[ipWetting] <= phaseMinHistoricalVolFraction[ipWetting] )
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
                                     m_imbibitionPhaseRelPermEndPoint[IPT::WETTING],
                                     phaseVolFraction[ipWetting],
                                     phaseMinHistoricalVolFraction[ipWetting],
                                     m_imbibitionPhaseMinVolFraction[IPT::WETTING],
                                     m_drainagePhaseMaxVolFraction[TPT::WETTING],
                                     m_imbibitionPhaseMaxVolFraction[IPT::WETTING],
                                     phaseRelPerm[ipWetting],
                                     dPhaseRelPerm_dPhaseVolFrac[ipWetting][ipWetting] );
  }

  // ---------- intermediate rel perm (ALWAYS DRAINAGE!)
  interRelPerm_wi =
    m_drainageRelPermKernelWrappers[TPT::INTERMEDIATE_WETTING].compute( &(phaseVolFraction)[ipInter],
                                                                        &dInterRelPerm_wi_dInterVolFrac );


  // 2) Non-wetting and intermediate phase relative permeabilities using two-phase non-wetting-intermediate data

  // ---------- non-wetting rel perm
  if( phaseVolFraction[ipNonWetting] >= phaseMaxHistoricalVolFraction[ipNonWetting] )
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
                                        m_drainagePhaseMinVolFraction[ipNonWetting],
                                        m_imbibitionPhaseMinVolFraction[IPT::NONWETTING],
                                        m_drainagePhaseMaxVolFraction[ipNonWetting],
                                        phaseVolFraction[ipNonWetting],
                                        phaseMaxHistoricalVolFraction[ipNonWetting],
                                        phaseRelPerm[ipNonWetting],
                                        dPhaseRelPerm_dPhaseVolFrac[ipNonWetting][ipNonWetting] );
  }

  // ---------- intermediate rel perm (ALWAYS DRAINAGE!)
  interRelPerm_nwi =
    m_drainageRelPermKernelWrappers[TPT::INTERMEDIATE_NONWETTING].compute( &(phaseVolFraction)[ipInter],
                                                                           &dInterRelPerm_nwi_dInterVolFrac );

  // 3) Compute the "three-phase" oil relperm

  // TODO FRANCOIS: double check that
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
inline void
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
inline void
TableRelativePermeabilityHysteresis::KernelWrapper::
  update( localIndex const k,
          localIndex const q,
          arraySlice1d< geosx::real64 const, compflow::USD_PHASE - 1 > const & phaseVolFraction ) const
{
  compute( phaseVolFraction,
           m_phaseMaxHistoricalVolFraction[k],
           m_phaseMinHistoricalVolFraction[k],
           m_phaseRelPerm[k][q],
           m_dPhaseRelPerm_dPhaseVolFrac[k][q] );
}

} // namespace constitutive

} // namespace geosx

#endif //GEOSX_CONSTITUTIVE_TABLERELATIVEPERMEABILITYHYSTERESIS_HPP
