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
 * @file EzrokhiBrineViscosity.cpp
 */

#include "constitutive/fluid/PVTFunctions/EzrokhiBrineViscosity.hpp"

#include "constitutive/fluid/PVTFunctions/PVTFunctionHelpers.hpp"
#include "functions/FunctionManager.hpp"

namespace geosx
{

using namespace stringutilities;

namespace constitutive
{

namespace PVTProps
{

namespace
{

TableFunction const * makeViscosityTable( string const & functionName,
                                          FunctionManager & functionManager )
{
  array1d< array1d< real64 > > temperatures;
  array1d< real64 > viscosities;

  temperatures.resize( 1 );
  temperatures[0].resize( 26 );
  viscosities.resize( 26 );

  temperatures[0][0] = 0.01;
  temperatures[0][1] = 10;
  temperatures[0][2] = 20;
  temperatures[0][3] = 25;
  temperatures[0][4] = 30;
  temperatures[0][5] = 40;
  temperatures[0][6] = 50;
  temperatures[0][7] = 60;
  temperatures[0][8] = 70;
  temperatures[0][9] = 80;
  temperatures[0][10] = 90;
  temperatures[0][11] = 100;
  temperatures[0][12] = 110;
  temperatures[0][13] = 120;
  temperatures[0][14] = 140;
  temperatures[0][15] = 160;
  temperatures[0][16] = 180;
  temperatures[0][17] = 200;
  temperatures[0][18] = 220;
  temperatures[0][19] = 240;
  temperatures[0][20] = 260;
  temperatures[0][21] = 280;
  temperatures[0][22] = 300;
  temperatures[0][23] = 320;
  temperatures[0][24] = 340;
  temperatures[0][25] = 360;


  viscosities[0] = 0.0017914;
  viscosities[1] = 0.0013060;
  viscosities[2] = 0.0010016;
  viscosities[3] = 0.0008900;
  viscosities[4] = 0.0007972;
  viscosities[5] = 0.0006527;
  viscosities[6] = 0.0005465;
  viscosities[7] = 0.0004660;
  viscosities[8] = 0.0004035;
  viscosities[9] = 0.0003540;
  viscosities[10] = 0.0003142;
  viscosities[11] = 0.0002816;
  viscosities[12] = 0.0002546;
  viscosities[13] = 0.0002320;
  viscosities[14] = 0.0001966;
  viscosities[15] = 0.0001704;
  viscosities[16] = 0.0001504;
  viscosities[17] = 0.0001346;
  viscosities[18] = 0.0001218;
  viscosities[19] = 0.0001111;
  viscosities[20] = 0.0001018;
  viscosities[21] = 0.0000936;
  viscosities[22] = 0.0000859;
  viscosities[23] = 0.0000783;
  viscosities[24] = 0.0000703;
  viscosities[25] = 0.0000603;

  string const tableName = functionName +  "_table";
  if( functionManager.hasGroup< TableFunction >( tableName ) )
  {
    return functionManager.getGroupPointer< TableFunction >( tableName );
  }
  else
  {
    TableFunction * const viscosityTable = dynamicCast< TableFunction * >( functionManager.createChild( TableFunction::catalogName(), tableName ) );
    viscosityTable->setTableCoordinates( temperatures );
    viscosityTable->setTableValues( viscosities );
    viscosityTable->setInterpolationMethod( TableFunction::InterpolationType::Linear );
    return viscosityTable;
  }
}

} // namespace

EzrokhiBrineViscosity::EzrokhiBrineViscosity( string const & name,
                                              string_array const & inputPara,
                                              string_array const & componentNames,
                                              array1d< real64 > const & componentMolarWeight ):
  PVTFunctionBase( name,
                   componentNames,
                   componentMolarWeight )
{
  string const expectedCO2ComponentNames[] = { "CO2", "co2" };
  m_CO2Index = PVTFunctionHelpers::findName( componentNames, expectedCO2ComponentNames, "componentNames" );

  string const expectedWaterComponentNames[] = { "Water", "water" };
  m_waterIndex = PVTFunctionHelpers::findName( componentNames, expectedWaterComponentNames, "componentNames" );

  makeCoefficients( inputPara );
  m_waterViscosityTable = makeViscosityTable( m_functionName, FunctionManager::getInstance() );
}

void EzrokhiBrineViscosity::makeCoefficients( string_array const & inputPara )
{
  // compute brine viscosity following Ezrokhi`s method (referenced in Eclipse TD, Aqueous phase properties)
  // Reference : Zaytsev, I.D. and Aseyev, G.G. Properties of Aqueous Solutions of Electrolytes, Boca Raton, Florida, USA CRC Press (1993).
  GEOSX_THROW_IF_LT_MSG( inputPara.size(), 5,
                         GEOSX_FMT( "{}: insufficient number of model parameters", m_functionName ),
                         InputError );

  try
  {
    // assume CO2 is the only non-water component in the brine
    m_coef0 = stod( inputPara[2] );
    m_coef1 = stod( inputPara[3] );
    m_coef2 = stod( inputPara[4] );
  }
  catch( std::invalid_argument const & e )
  {
    GEOSX_THROW( GEOSX_FMT( "{}: invalid model parameter value '{}'", m_functionName, e.what() ), InputError );
  }
}

EzrokhiBrineViscosity::KernelWrapper
EzrokhiBrineViscosity::createKernelWrapper() const
{
  return KernelWrapper( m_componentMolarWeight,
                        *m_waterViscosityTable,
                        m_CO2Index,
                        m_waterIndex,
                        m_coef0,
                        m_coef1,
                        m_coef2 );
}

REGISTER_CATALOG_ENTRY( PVTFunctionBase, EzrokhiBrineViscosity, string const &, string_array const &, string_array const &, array1d< real64 > const & )

} // end namespace PVTProps

} // namespace constitutive

} // end namespace geosx
