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
 * @file BrineViscosity.cpp
 */

#include "constitutive/fluid/PVTFunctions/BrineViscosityFunction.hpp"

namespace geosx
{

using namespace stringutilities;

namespace PVTProps
{

BrineViscosity::BrineViscosity( array1d< string > const & inputPara,
                                array1d< string > const & componentNames,
                                array1d< real64 > const & componentMolarWeight ):
  PVTFunctionBase( inputPara[1],
                   componentNames,
                   componentMolarWeight )
{
  makeCoefficients( inputPara );
}

void BrineViscosity::makeCoefficients( array1d< string > const & inputPara )
{
  // these coefficients come from Phillips et al. (1981), equation (1), pages 5-6
  constexpr real64 a = 0.0816;
  constexpr real64 b = 0.0122;
  constexpr real64 c = 0.000128;
  constexpr real64 d = 0.000629;
  constexpr real64 k = -0.7;
  constexpr real64 waterVisc = 8.9e-4; //at 25C

  real64 m = -1.0;

  GEOSX_ERROR_IF( inputPara.size() < 3, "Invalid BrineViscosity input!" );

  try
  {
    m = stod( inputPara[2] );
  }
  catch( const std::invalid_argument & e )
  {
    GEOSX_ERROR( "Invalid BrineViscosity argument:" + string( e.what() ) );
  }

  m_coef0 = (1.0 + a * m + b * m * m + c * m * m * m) * waterVisc;
  m_coef1 =  d * (1.0 - exp( k * m )) * waterVisc;
}

BrineViscosity::KernelWrapper BrineViscosity::createKernelWrapper()
{
  return KernelWrapper( m_componentNames,
                        m_componentMolarWeight,
                        m_coef0,
                        m_coef1 );
}

REGISTER_CATALOG_ENTRY( PVTFunctionBase, BrineViscosity, array1d< string > const &, array1d< string > const &, array1d< real64 > const & )

} // end namespace PVTProps

} // end namespace geosx
