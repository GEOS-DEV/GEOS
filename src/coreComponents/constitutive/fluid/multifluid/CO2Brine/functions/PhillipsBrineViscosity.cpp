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
 * @file PhillipsBrineViscosity.cpp
 */

#include "constitutive/fluid/multifluid/CO2Brine/functions/PhillipsBrineViscosity.hpp"

#include "constitutive/fluid/multifluid/CO2Brine/functions/PureWaterProperties.hpp"
#include "functions/FunctionManager.hpp"

namespace geos
{

using namespace stringutilities;

namespace constitutive
{

namespace PVTProps
{

PhillipsBrineViscosity::PhillipsBrineViscosity( string const & name,
                                                string_array const & inputPara,
                                                string_array const & componentNames,
                                                array1d< real64 > const & componentMolarWeight,
                                                bool const printInCsv,
                                                bool const printInLog ):
  PVTFunctionBase( name,
                   componentNames,
                   componentMolarWeight )
{
  m_waterViscosityTable = PureWaterProperties::makeSaturationViscosityTable( m_functionName, FunctionManager::getInstance() );
  makeCoefficients( inputPara );

  checkPrint( m_waterViscosityTable, printInCsv, printInLog );
}

void PhillipsBrineViscosity::makeCoefficients( string_array const & inputPara )
{
  GEOS_THROW_IF_LT_MSG( inputPara.size(), 3,
                        GEOS_FMT( "{}: insufficient number of model parameters", m_functionName ),
                        InputError );

  real64 m;
  try
  {
    m = stod( inputPara[2] );
  }
  catch( std::invalid_argument const & e )
  {
    GEOS_THROW( GEOS_FMT( "{}: invalid model parameter value '{}'", m_functionName, e.what() ), InputError );
  }

  // these coefficients come from Phillips et al. (1981), equation (1), pages 5-6
  constexpr real64 a = 0.0816;
  constexpr real64 b = 0.0122;
  constexpr real64 c = 0.000128;
  constexpr real64 d = 0.000629;
  constexpr real64 k = -0.7;

  // precompute the model coefficients
  // (excluding water viscosity, which will multiply them in the compute function)
  m_coef0 = (1.0 + a * m + b * m * m + c * m * m * m);
  m_coef1 =  d * (1.0 - exp( k * m ));
}

void PhillipsBrineViscosity::checkTablesParameters( real64 const GEOS_UNUSED_PARAM( pressure ),
                                                    real64 const temperature ) const
{
  m_waterViscosityTable->checkCoord( temperature, 0 );
}

PhillipsBrineViscosity::KernelWrapper
PhillipsBrineViscosity::createKernelWrapper() const
{
  return KernelWrapper( m_componentMolarWeight,
                        *m_waterViscosityTable,
                        m_coef0,
                        m_coef1 );
}

REGISTER_CATALOG_ENTRY( PVTFunctionBase, PhillipsBrineViscosity, string const &, string_array const &, string_array const &, array1d< real64 > const &, bool const, bool const )

} // end namespace PVTProps

} // namespace constitutive

} // end namespace geos
