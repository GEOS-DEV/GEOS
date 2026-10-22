/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 Total, S.A
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file KValueFlashParameters.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_PARAMETERS_KVALUEFLASHPARAMETERS_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_PARAMETERS_KVALUEFLASHPARAMETERS_HPP_

#include "PressureTemperatureCoordinates.hpp"

namespace geos
{

class TableFunction;

namespace constitutive
{

namespace compositional
{

/**
 * @brief Manages the parameters and data structures required for K-Value based
 * flash calculations, specifically using pre-tabulated K-values.
 *
 * This class inherits from PressureTemperatureCoordinates and holds the static
 * lookup tables for K-values ($K = y_i / x_i$) which are used to model the
 * equilibrium between phases. The tables allow for non-uniform spacing of
 * pressure and temperature coordinate points.
 *
 * @tparam NUM_PHASE The total number of phases present in the system, including
 * the primary phase.
 */
template< integer NUM_PHASE >
class KValueFlashParameters : public PressureTemperatureCoordinates
{
  static constexpr integer numPhases = NUM_PHASE;
public:
  KValueFlashParameters( std::unique_ptr< ModelParameters > parameters );
  ~KValueFlashParameters() override = default;

  static std::unique_ptr< ModelParameters > create( std::unique_ptr< ModelParameters > parameters );

  /**
   * @brief Creates and registers the necessary pressure and temperature index tables
   * based on the coordinate arrays stored in this object.
   *
   * These tables allow conversion from actual P/T values to the indices needed
   * for accessing the m_kValueHyperCube.
   *
   * @param name The base name for the tables.
   * @param pressureTableName Output reference to store the name of the created pressure index table.
   * @param temperatureTableName Output reference to store the name of the created temperature index table.
   */
  void createTables( string const & name,
                     string & pressureTableName,
                     string & temperatureTableName ) const;

  /**
   * @brief Array of strings holding the names of the Functions used to generate
   * the K-values. There is one TableFunction per component and per secondary phase.
   * The functions are evaluated on the host to populate the hypercube.
   */
  string_array m_kValueTables;

  array1d< array1d< real64 > > m_pressureValues;
  array1d< array1d< real64 > > m_temperatureValues;

  /**
   * @brief The core 4D data structure containing the pre-calculated K-values.
   * The indices correspond to: (Phase Index, Component Index, Pressure Index, Temperature Index).
   */
  array4d< real64 > m_kValueHyperCube;

  struct viewKeyStruct : public PressureTemperatureCoordinates::viewKeyStruct
  {
    static constexpr char const * kValueTablesString() { return "kValueTables"; }
    static constexpr char const * kValuePressureCoordinatesString() { return "kValuePressureCoordinates"; }
    static constexpr char const * kValueTemperatureCoordinatesString() { return "kValueTemperatureCoordinates"; }
    static constexpr char const * kValueHyperCubeString() { return "kValueHyperCube"; }
  };

protected:
  void registerParametersImpl( MultiFluidBase * fluid ) override;
  void postInputInitializationImpl( MultiFluidBase const * fluid, ComponentProperties const & componentProperties ) override;

private:
  /**
   * @brief Populates the 4D m_kValueHyperCube by evaluating the TableFunctions
   * listed in m_kValueTables at all pressure/temperature grid points.
   * @param numComps The total number of components in the system.
   */
  void generateHyperCube( integer const numComps );

  /**
   * @brief Performs validation checks on the loaded or generated K-values
   * to ensure they are physically sound (e.g., non-negative).
   * @param fluid Constant pointer to the MultiFluidBase object for context.
   * @return true if K-values are valid, false otherwise.
   */
  bool validateKValues( MultiFluidBase const * fluid ) const;
};

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_PARAMETERS_KVALUEFLASHPARAMETERS_HPP_
