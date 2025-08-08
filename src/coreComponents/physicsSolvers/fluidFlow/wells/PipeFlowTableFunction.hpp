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
 * @file PipeFlowTableFunction.hpp
 */

#ifndef GEOS_FUNCTIONS_PIPEFLOWTABLEFUNCTION_HPP_
#define GEOS_FUNCTIONS_PIPEFLOWTABLEFUNCTION_HPP_

#include "functions/MultivariableNonuniformTableFunction.hpp"

#include "common/format/EnumStrings.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geos
{

/**
 * @class PipeFlowTableFunction
 *
 * An interface class for pipeflow table function (function with multiple inputs and outputs) with uniform discretization
 */

class PipeFlowTableFunction : public MultivariableNonuniformTableFunction
{
public:

  /**
   * @brief The constructor
   * @param[in] name the name of this object manager
   * @param[in] parent the parent Group
   */
  PipeFlowTableFunction( const string & name,
                         Group * const parent );

  /**
   * @brief The catalog name interface
   * @return name of the PipeFlowTableFunction in the FunctionBase catalog
   */
  static string catalogName() { return "PipeFlowTableFunction"; }

    /**
   * @name Getters / Setters
   */
  ///@{

    /**
   * @brief Get type of rate array
   * @return name of type
   */
  const string & getRateType() const { return m_rateType; }

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

    ///@}

  /**
   * @brief Struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct
  {
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

    /// @return String key for "bhp" array
    static constexpr char const *bottomHolePressureArray() { return "bottomHolePressure"; }
  }
  /// ViewKey struct for the Perforation class
  viewKeysPipeFlowTableFunction;

protected:
  virtual void postInputInitialization() override;
  /**
   * @brief Initialize the table function after setting table coordinates and values
   */
  virtual void initializeFunction() override;
private:

  /// Rate
  string m_rateType;
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

  /// bottom hole pressure
  array1d< real64 > m_bhp;

  // Table function for evaluating values and derivatives
  MultivariableNonuniformTableFunction * m_tableFunction;
};


} /* namespace geos */

#endif /* GEOS_FUNCTIONS_PipeFlowTableFunction_HPP_ */
