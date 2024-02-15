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
#include "codingUtilities/Table.hpp"

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
                                                bool const printTable ):
  PVTFunctionBase( name,
                   componentNames,
                   componentMolarWeight )
{
  m_waterViscosityTable = PureWaterProperties::makeSaturationViscosityTable( m_functionName, FunctionManager::getInstance() );
  if( printTable )
    m_waterViscosityTable->print( m_waterViscosityTable->getName() );
  makeCoefficients( inputPara );

}

void PhillipsBrineViscosity::debugViscosityTable() const
{
  if( MpiWrapper::commRank( MPI_COMM_GEOSX ) != 0 )
  {
    return;
  }

  Table tablePerforation = Table( {"Perforation no." } );

  std::cout << " numCoord " <<  m_waterViscosityTable->numDimensions() << std::endl;
  std::cout << " dimUnit " <<  m_waterViscosityTable->getDimUnit( 0 ) << std::endl;
  std::cout << " Unit " << units::getDescription( units::Unit::Viscosity ) << std::endl;

  arrayView1d< real64 const > viscosity = m_waterViscosityTable->getValues();
  ArrayOfArraysView< real64 const > coords = m_waterViscosityTable->getCoordinates();
  arraySlice1d< real64 const > tempVar = coords[0];
  arraySlice1d< real64 const > pressure = coords[1];
  // for( auto value : m_waterViscosityTable->getValues() )
  // {
  //   std::cout << " value m_water : " <<  value << std::endl;

  // }
  for( localIndex i = 1; i < coords.sizeOfArray( 0 ); ++i )
  {
    tablePerforation.addRow< 1 >( tempVar[i] );
  }

 tablePerforation.draw(std::cout);
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

REGISTER_CATALOG_ENTRY( PVTFunctionBase, PhillipsBrineViscosity, string const &, string_array const &, string_array const &, array1d< real64 > const &, bool const )

} // end namespace PVTProps

} // namespace constitutive

} // end namespace geos
