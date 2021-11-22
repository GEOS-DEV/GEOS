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
 * @file SpanWagnerCO2DensityFunction.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_SPANWAGNERCO2DENSITYFUNCTION_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_SPANWAGNERCO2DENSITYFUNCTION_HPP_

#include "constitutive/fluid/PVTFunctions/old/PVTFunctionBase.hpp"

namespace geosx
{

namespace PVTProps
{

class SpanWagnerCO2DensityFunction : public PVTFunction
{
public:


  SpanWagnerCO2DensityFunction( string_array const & inputPara,
                                string_array const & componentNames,
                                real64_array const & componentMolarWeight );

  ~SpanWagnerCO2DensityFunction() override {}


  static constexpr auto m_catalogName = "SpanWagnerCO2Density";
  static string catalogName()                    { return m_catalogName; }
  virtual string getCatalogName() const override final { return catalogName(); }


  virtual PVTFuncType functionType() const override
  {
    return PVTFuncType::DENSITY;
  }

  virtual void evaluation( EvalVarArgs const & pressure, EvalVarArgs const & temperature, arraySlice1d< EvalVarArgs const > const & phaseComposition,
                           EvalVarArgs & value, bool useMass = 0 ) const override;

  static void calculateCO2Density( real64_array const & pressure, real64_array const & temperature, real64_array2d const & density );

private:

  void makeTable( string_array const & inputPara );

  static void spanWagnerCO2Density( real64 const & T, real64 const & P, real64 & rho, real64 (*f)( real64 const & x1, real64 const & x2, real64 const & x3 ));


  TableFunctionPtr m_CO2DensityTable;
  localIndex m_CO2Index;

};

}

}

#endif //GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_SPANWAGNERCO2DENSITYFUNCTION_HPP_
