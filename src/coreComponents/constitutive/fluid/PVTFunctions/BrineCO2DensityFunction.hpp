/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file BrineCO2DensityFunction.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_BRINECO2DENSITYFUNCTION_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_BRINECO2DENSITYFUNCTION_HPP_

#include "PVTFunctionBase.hpp"

namespace geosx
{

namespace PVTProps
{

class BrineCO2DensityFunction : public PVTFunction
{
public:

  BrineCO2DensityFunction( string_array const & inputPara,
                           string_array const & componentNames,
                           real64_array const & componentMolarWeight );
  ~BrineCO2DensityFunction() override {}


  static constexpr auto m_catalogName = "BrineCO2Density";
  static string CatalogName()                    { return m_catalogName; }
  virtual string GetCatalogName() override final { return CatalogName(); }


  virtual PVTFuncType FunctionType() const override
  {
    return PVTFuncType::DENSITY;

  }

  virtual void Evaluation( EvalVarArgs const & pressure,
                           EvalVarArgs const & temperature,
                           arraySlice1d< EvalVarArgs const > const & phaseComposition,
                           EvalVarArgs & value,
                           bool useMass = 0 ) const override;


private:

  void MakeTable( string_array const & inputPara );

  void CalculateBrineDensity( real64_array const & pressure, real64_array const & temperature, real64 const & salinity, real64_array2d const & density );

  TableFunctionPtr m_BrineDensityTable;
  localIndex m_CO2Index;
  localIndex m_waterIndex;

};

}

}

#endif //GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_BRINECO2DENSITYFUNCTION_HPP_
