/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file BlackOilFluidBase.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_BLACKOILFLUIDBASE_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_BLACKOILFLUIDBASE_HPP_

#include "constitutive/fluid/MultiFluidBase.hpp"
#include "functions/TableFunction.hpp"


namespace geosx
{

namespace constitutive
{

class BlackOilFluidBase : public MultiFluidBase
{
public:

  using exec_policy = parallelDevicePolicy<>;

  struct PhaseType
  {
    static constexpr integer OIL            = 0;
    static constexpr integer GAS            = 1;
    static constexpr integer WATER          = 2;
    static constexpr integer MAX_NUM_PHASES = 3;
  };

  BlackOilFluidBase( string const & name, Group * const parent );

  virtual ~BlackOilFluidBase() override = default;

  struct viewKeyStruct : MultiFluidBase::viewKeyStruct
  {
    static constexpr char const * surfacePhaseMassDensitiesString() { return "surfaceDensities"; }
    static constexpr char const * tableFilesString() { return "tableFiles"; }
    static constexpr char const * formationVolumeFactorTableNamesString() { return "hydrocarbonFormationVolFactorTableNames"; }
    static constexpr char const * viscosityTableNamesString() { return "hydrocarbonViscosityTableNames"; }
    static constexpr char const * waterRefPressureString() { return "waterReferencePressure"; }
    static constexpr char const * waterFormationVolumeFactorString() { return "waterFormationVolumeFactor"; }
    static constexpr char const * waterCompressibilityString() { return "waterCompressibility"; }
    static constexpr char const * waterViscosityString() { return "waterViscosity"; }

  };

protected:

  virtual void postProcessInput() override;

  virtual void initializePostSubGroups() override;

  /**
   * @brief Create all the table kernel wrappers
   */
  virtual void createAllKernelWrappers();

  /**
   * @brief Use the TableFunctions provided by the user to get the PVT data
   */
  virtual void useProvidedTableFunctions() = 0;

  /**
   * @brief Read all the PVT table provided by the user in Eclipse format
   */
  virtual void readInputDataFromPVTFiles() = 0;

  /**
   * @brief Fill the water data (formation vol factor, compressibility, etc)
   * @param[in] tableValues the values in the water table
   */
  void fillWaterData( array1d< array1d< real64 > > const & tableValues );

  /**
   * @brief Fill the hydrocarbon data (pressure, formation vol factor, viscosity)
   * @param[in] ip the index of the phase
   * @param[in] tableValues the values in the oil or gas table
   */
  void fillHydrocarbonData( localIndex const ip,
                            array1d< array1d< real64 > > const & tableValues );

  /**
   * @brief Check that the table values make sense
   * @param[in] table the values in the oil or gas table
   */
  void validateTable( TableFunction const & table ) const;


  // Input data

  // Black-oil table filenames
  path_array m_tableFiles;

  // Fluid data

  /// Names of the formation volume factor tables
  array1d< string > m_formationVolFactorTableNames;

  /// Names of the viscosity tables
  array1d< string > m_viscosityTableNames;

  /// Surface densities
  array1d< real64 > m_surfacePhaseMassDensity;

  /// Water reference pressure
  real64 m_waterRefPressure;

  /// Water formation volume factor
  real64 m_waterFormationVolFactor;

  /// Water compressibility
  real64 m_waterCompressibility;

  /// Water viscosity
  real64 m_waterViscosity;

  /// Data after processing of input

  /// Phase ordering info
  array1d< integer > m_phaseTypes;
  array1d< integer > m_phaseOrder;
  array1d< integer > m_hydrocarbonPhaseOrder;

  /// Table kernel wrappers to interpolate in the oil and gas (B vs p) tables
  array1d< TableFunction::KernelWrapper > m_formationVolFactorTables;

  /// Table kernel wrappers to interpolate in the oil and gas (\mu vs p) tables
  array1d< TableFunction::KernelWrapper > m_viscosityTables;

};

} //namespace constitutive

} //namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_BLACKOILFLUIDBASE_HPP_
