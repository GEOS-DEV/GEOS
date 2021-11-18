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

#include "constitutive/fluid/PVTFunctions/PhillipsBrineViscosity.hpp"

namespace geosx
{

using namespace stringutilities;

namespace constitutive
{

namespace PVTProps
{

PhillipsBrineViscosity::PhillipsBrineViscosity( string const & name,
                                                string_array const & inputPara,
                                                string_array const & componentNames,
                                                array1d< real64 > const & componentMolarWeight ):
  PVTFunctionBase( name,
                   componentNames,
                   componentMolarWeight )
{
  makeCoefficients( inputPara );
}

void PhillipsBrineViscosity::makeCoefficients( string_array const & inputPara )
{
  // these coefficients come from Phillips et al. (1981), equation (1), pages 5-6
  constexpr real64 a = 0.0816;
  constexpr real64 b = 0.0122;
  constexpr real64 c = 0.000128;
  constexpr real64 d = 0.000629;
  constexpr real64 k = -0.7;
  constexpr real64 waterVisc = 8.9e-4; //at 25C

  GEOSX_THROW_IF_LT_MSG( inputPara.size(), 3,
                         GEOSX_FMT( "{}: insufficient number of model parameters", m_functionName ),
                         InputError );

  real64 m;
  try
  {
    m = stod( inputPara[2] );
  }
  catch( std::invalid_argument const & e )
  {
    GEOSX_THROW( GEOSX_FMT( "{}: invalid model parameter value '{}'", m_functionName, e.what() ), InputError );
  }

  m_coef0 = (1.0 + a * m + b * m * m + c * m * m * m) * waterVisc;
  m_coef1 =  d * (1.0 - exp( k * m )) * waterVisc;
}

PhillipsBrineViscosity::KernelWrapper
PhillipsBrineViscosity::createKernelWrapper() const
{
  return KernelWrapper( m_componentMolarWeight,
                        m_coef0,
                        m_coef1 );
}

REGISTER_CATALOG_ENTRY( PVTFunctionBase, PhillipsBrineViscosity, string const &, string_array const &, string_array const &, array1d< real64 > const & )

} // end namespace PVTProps

} // namespace constitutive

} // end namespace geosx
