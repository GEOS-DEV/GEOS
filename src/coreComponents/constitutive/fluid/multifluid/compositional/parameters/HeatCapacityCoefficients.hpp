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
 * @file HeatCapacityCoefficients.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_PARAMETERS_HEATCAPACITYCOEFFICIENTS_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_PARAMETERS_HEATCAPACITYCOEFFICIENTS_HPP_

#include "ModelParameters.hpp"
#include "PhaseType.hpp"
#include "common/DataTypes.hpp"

namespace geos
{

namespace constitutive
{

namespace compositional
{

class HeatCapacityCoefficients : public ModelParameters
{
public:
  HeatCapacityCoefficients( std::unique_ptr< ModelParameters > parameters );
  ~HeatCapacityCoefficients() override = default;

  static std::unique_ptr< ModelParameters > create( std::unique_ptr< ModelParameters > parameters );

  struct viewKeyStruct
  {
    static constexpr char const * enthalpyReferenceTemperatureString() { return "enthalpyReferenceTemperature"; }
    static constexpr char const * referenceEnthalpyString() { return "referenceEnthalpy"; }
    static constexpr char const * componentHeatCapacityCoefficientsString() { return "componentHeatCapacityCoefficients"; }
  };

  real64 m_referenceTemperature{298.15};
  array2d< real64 > m_referenceEnthalpy;
  array2d< real64 > m_coefficients;
  array1d< integer > m_phaseTypes;

protected:
  void registerParametersImpl( MultiFluidBase * fluid ) override;

  void postInputInitializationImpl( MultiFluidBase const * fluid, ComponentProperties const & componentProperties ) override;

private:
  static bool isPolynomialPositive( arraySlice1d< real64 const > const a, real64 const T0, real64 const T1, real64 & T, real64 & hT );
};

} // end namespace compositional

} // end namespace constitutive

} // end namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_PARAMETERS_HEATCAPACITYCOEFFICIENTS_HPP_
