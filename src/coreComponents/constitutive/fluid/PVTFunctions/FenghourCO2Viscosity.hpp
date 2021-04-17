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
 * @file FenghourCO2Viscosity.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_FENGHOURCO2VISCOSITY_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_FENGHOURCO2VISCOSITY_HPP_

#include "PVTFunctionBase.hpp"

#include "constitutive/fluid/PVTFunctions/PVTFunctionHelpers.hpp"
#include "functions/TableFunction.hpp"

namespace geosx
{

namespace constitutive
{

namespace PVTProps
{

class FenghourCO2ViscosityUpdate final : public PVTFunctionBaseUpdate
{
public:

  FenghourCO2ViscosityUpdate( arrayView1d< string const > const & componentNames,
                              arrayView1d< real64 const > const & componentMolarWeight,
                              TableFunction * CO2ViscosityTable )
    : PVTFunctionBaseUpdate( componentNames,
                             componentMolarWeight ),
    m_CO2ViscosityTable( CO2ViscosityTable->createKernelWrapper() )
  {}

  /// Default copy constructor
  FenghourCO2ViscosityUpdate( FenghourCO2ViscosityUpdate const & ) = default;

  /// Default move constructor
  FenghourCO2ViscosityUpdate( FenghourCO2ViscosityUpdate && ) = default;

  /// Deleted copy assignment operator
  FenghourCO2ViscosityUpdate & operator=( FenghourCO2ViscosityUpdate const & ) = delete;

  /// Deleted move assignment operator
  FenghourCO2ViscosityUpdate & operator=( FenghourCO2ViscosityUpdate && ) = delete;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void compute( real64 const & pressure,
                        real64 const & temperature,
                        arraySlice1d< real64 const > const & phaseComposition,
                        arraySlice1d< real64 const > const & dPhaseComposition_dPressure,
                        arraySlice1d< real64 const > const & dPhaseComposition_dTemperature,
                        arraySlice2d< real64 const > const & dPhaseComposition_dGlobalCompFraction,
                        real64 & value,
                        real64 & dValue_dPressure,
                        real64 & dValue_dTemperature,
                        arraySlice1d< real64 > const & dValue_dGlobalCompFraction,
                        bool useMass = 0 ) const override;

  /**
   * @brief Move the KernelWrapper to the given execution space, optionally touching it.
   * @param space the space to move the KernelWrapper to
   * @param touch whether the KernelWrapper should be touched in the new space or not
   * @note This function exists to enable holding KernelWrapper objects in an ArrayView
   *       and have their contents properly moved between memory spaces.
   */
  void move( LvArray::MemorySpace const space, bool const touch = false )
  {
    m_CO2ViscosityTable.move( space, touch );
  }

protected:

  /// Table with viscosity tabulated as a function (P,T)
  TableFunction::KernelWrapper m_CO2ViscosityTable;

};

class FenghourCO2Viscosity : public PVTFunctionBase
{
public:

  FenghourCO2Viscosity( string_array const & inputPara,
                        string_array const & componentNames,
                        array1d< real64 > const & componentMolarWeight );
  ~FenghourCO2Viscosity() override {}

  static string catalogName() { return "FenghourCO2Viscosity"; }

  virtual string getCatalogName() const override final { return catalogName(); }

  virtual PVTFunctionType functionType() const override
  {
    return PVTFunctionType::VISCOSITY;
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = FenghourCO2ViscosityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper();

private:

  void makeTable( string_array const & inputPara );

  void calculateCO2Viscosity( PVTProps::PTTableCoordinates const & tableCoords,
                              array1d< real64 > const & density,
                              array1d< real64 > const & viscosity );

  void fenghourCO2ViscosityFunction( real64 const & temperatureCent,
                                     real64 const & density,
                                     real64 & viscosity );

  /// Table with CO2 viscosity tabulated as a function of (P,T)
  TableFunction * m_CO2ViscosityTable;

};

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void FenghourCO2ViscosityUpdate::compute( real64 const & pressure,
                                          real64 const & temperature,
                                          arraySlice1d< real64 const > const & phaseComposition,
                                          arraySlice1d< real64 const > const & dPhaseComposition_dPressure,
                                          arraySlice1d< real64 const > const & dPhaseComposition_dTemperature,
                                          arraySlice2d< real64 const > const & dPhaseComposition_dGlobalCompFraction,
                                          real64 & value,
                                          real64 & dValue_dPressure,
                                          real64 & dValue_dTemperature,
                                          arraySlice1d< real64 > const & dValue_dGlobalCompFraction,
                                          bool useMass ) const
{
  GEOSX_UNUSED_VAR( phaseComposition );
  GEOSX_UNUSED_VAR( dPhaseComposition_dPressure );
  GEOSX_UNUSED_VAR( dPhaseComposition_dTemperature );
  GEOSX_UNUSED_VAR( dPhaseComposition_dGlobalCompFraction );
  GEOSX_UNUSED_VAR( useMass );

  real64 const input[2] = { pressure, temperature };
  real64 densityDeriv[2]{};
  m_CO2ViscosityTable.compute( input, value, densityDeriv );
  dValue_dPressure = densityDeriv[0];
  dValue_dTemperature = densityDeriv[1];

  for( real64 & val : dValue_dGlobalCompFraction )
  {
    val = 0.0;
  }
}

} // end namespace PVTProps

} // end namespace constitutive

} // end namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_FENGHOURCO2VISCOSITY_HPP_
