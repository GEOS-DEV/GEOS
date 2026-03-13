/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file InjPipeFlowTableFunction.hpp
 */

#ifndef GEOS_FUNCTIONS_INJPIPEFLOWTABLEFUNCTION_HPP_
#define GEOS_FUNCTIONS_INJPIPEFLOWTABLEFUNCTION_HPP_
#include "functions/TableFunction.hpp"
#include "functions/MultivariableNonuniformTableFunction.hpp"

#include "common/format/EnumStrings.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{

/**
 * @class InjPipeFlowTableFunction
 *
 * An interface class for pipeflow table function (function with multiple inputs and outputs) with uniform discretization
 */

class InjPipeFlowTableFunction : public MultivariableNonuniformTableFunction
{
public:

  /**
   * @brief The constructor
   * @param[in] name the name of this object manager
   * @param[in] parent the parent Group
   */
  InjPipeFlowTableFunction( const string & name,
                            Group * const parent );

  /**
   * @brief The catalog name interface
   * @return name of the InjPipeFlowTableFunction in the FunctionBase catalog
   */
  static string catalogName() { return "InjPipeFlowTableFunction"; }

  /**
   * @name Getters / Setters
   */
  ///@{

  /**
   * @brief Get type of rate array
   * @return name of type
   */
  const string & getTableName() const { return m_tableName; }

  /**
   * @brief Get type of rate array
   * @return name of type
   */
  const string & getRateType() const { return m_rateType; }

  /**
   * @brief Get phases associated with rate type
   * @return array of phases
   */
  const string_array & getRatePhases() const { return m_ratePhases; }

  /**
   * @brief Get rates table coordinates
   * @return array of phases
   */
  const array1d< real64 > & getRates() const { return m_rate; }

  /**
   * @brief Get type of water fraction array
   * @return name of type
   */
  const string & getWaterFractionType() const { return m_waterFractionType; }

  /**
   * @brief Get type of gas fraction array
   * @return name of type
   */
  const string & getGasFractionType() const { return m_gasFractionType; }

  /**
   * @brief Get type of gas fraction array
   * @param[in] rate The flow rate to bracket
   * @param[out] b0 The lower bracket index
   * @param[out] b1 The upper bracket index
   * @return 0 if rate is below table min, 2 if rate is above table max, else 1
   */
  integer getRateBracket( real64 const & rate, integer & b0, integer & b1 ) const;

  /**
   * @brief Get dQdP for bracket containing rate
   * @param[in] phaseRates array of phase rates, indexing according to phase ordering in fluid model
   * @param[in] whp well head pressure
   * @param[out] ql0 liquid rate at lower bracket
   * @param[out] ql1 liquid rate at upper bracket
   * @param[out] bhp0 bottom hole pressure at lower bracket
   * @param[out] bhp1 bottom hole pressure at upper bracket
   * @return dQdP value
   */

  real64 calculatedPdQ( array1d< real64 > const & phaseRates, real64 const & whp, real64 & ql0, real64 & ql1, real64 & bhp0, real64 & bhp1 ) const;

  ///@}

  /**
   * @brief Struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct
  {
    /// @return String key type of table name
    static constexpr char const *tableName() { return "name"; }

    /// @return String key type of flow rate associated with the "rate" array
    static constexpr char const *rateType() { return "rateType"; }
    /// @return String key for "rate" array
    static constexpr char const *rateArray() { return "rate"; }


    /// @return String key for "whp" array
    static constexpr char const *wellHeadPressureArray() { return "wellHeadPressure"; }

    /// @return String key type of water fraction associated with the "wfr" array
    static constexpr char const *waterFractionType() { return "waterFractionType"; }
    /// @return String key for "wfr" array
    static constexpr char const *waterFractionArray() { return "wfr"; }

    /// @return String key type of gass fraction associated with the "gfr" array
    static constexpr char const *gasFractionType() { return "gasFractionType"; }
    /// @return String key for "wfr" array
    static constexpr char const *gasFractionArray() { return "gfr"; }

    /// @return String key type of gas lift rate associated with the "gasLift" array
    static constexpr char const *gasLiftType() { return "gasLiftType"; }
    /// @return String key for "gasLift" array
    static constexpr char const *gasLiftArray() { return "gasLift"; }

    /// @return String key for "bhp" array
    static constexpr char const *bottomHolePressureArray() { return "bottomHolePressure"; }
  }
  /// ViewKey struct for the Perforation class
  viewKeysInjPipeFlowTableFunction;

  void calculateBHP( array1d< real64 > const & phaseRates, real64 const & whp, real64 & bhp, integer & solveStat ) const;

  void calculateWHP( const std::string & wellName, real64 const & bhp, array1d< real64 > const & phaseRates, real64 & whp, integer & solveStat ) const;

  void writeTable() const;

  virtual void postInputInitialization() override;

protected:

  /**
   * @brief Initialize the table function after setting table coordinates and values
   */
  virtual void initializeFunction() override;
  //void initializeFunction1()  ;
private:

  /// Name of the flow table
  string m_tableName;

  /// Rate
  string m_rateType;
  string_array m_ratePhases;
  array1d< real64 > m_rate;


  /// Well head pressure
  string m_pressureType;
  array1d< real64 > m_whp;

  /// Water fraction
  string m_waterFractionType;
  array1d< real64 > m_wfr;

  /// gas fraction
  string m_gasFractionType;
  array1d< real64 > m_gfr;

  /// gas lift rate
  string m_gasLiftType;
  array1d< real64 > m_gasLift;

  /// bottom hole pressure
  array1d< real64 > m_bhp;

  // Table function for evaluating values and derivatives
  MultivariableNonuniformTableFunction * m_tableFunction;
  TableFunction * m_tableFunction0;

};


} /* namespace geos */

#endif /* GEOS_FUNCTIONS_INJPIPEFLOWTABLEFUNCTION_HPP_ */
