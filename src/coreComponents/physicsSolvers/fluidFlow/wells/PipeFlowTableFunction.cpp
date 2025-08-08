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
 * @file PipeFlowTableFunction.cpp
 */

#include "PipeFlowTableFunction.hpp"
#include "functions/FunctionManager.hpp"
#include "common/DataTypes.hpp"
#include <algorithm>

namespace geos
{

using namespace dataRepository;

PipeFlowTableFunction::PipeFlowTableFunction( const string & name,
                                              Group * const parent ):
  MultivariableNonuniformTableFunction( name, parent ),
  m_tableFunction( nullptr )
{

  registerWrapper( viewKeyStruct::rateType(), &m_rateType ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Type of rate entered in the rates array. Valid entires are ..." );

  registerWrapper( viewKeyStruct::rateArray(), &m_rate ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Array of rates" );

  registerWrapper( viewKeyStruct::waterFractionType(), &m_waterFractionType ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Type of water fraction  entered in the wfr array. Valid entires are ..." );

  registerWrapper( viewKeyStruct::waterFractionArray(), &m_wfr ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Array of water fractions " );

  registerWrapper( viewKeyStruct::gasFractionType(), &m_gasFractionType ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Type of gas fraction  entered in the wfr array. Valid entires are ..." );

  registerWrapper( viewKeyStruct::gasFractionArray(), &m_gfr ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Array of gas fractions " );

  registerWrapper( viewKeyStruct::wellHeadPressureArray(), &m_whp ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Array of well head pressures " );


  registerWrapper( viewKeyStruct::bottomHolePressureArray(), &m_bhp ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Array of bottom hole pressures representing the dependent variable for the table function.\n"
                    "The quantities must be entered such that the rate index increases first, followed by whp, waterFraction, and finally gasFraction.\n"
                    " For example  bhp[1]= f(rate[1],whp[1],,waterFraction[1],gasFraction[1]) , bhp[2]= f(rate[2],whp[1],waterFraction[1],gasFraction[1])" );

}

void PipeFlowTableFunction::postInputInitialization()
{
  // Validate independent table values are not decreasing
  auto checkNotDecreasing = []( const auto & arr, const std::string & name )
  {
    for( auto i = 1; i < arr.size(); ++i )
    {
      if( arr[i] < arr[i-1] )
      {
        GEOS_ERROR( name << " array values must be non-decreasing, but arr[" << i << "] < arr[" << i-1 << "]" );
      }
    }
  };

  checkNotDecreasing( m_rate, getRateType());
  checkNotDecreasing( m_whp, "wellHeadPressure" );
  checkNotDecreasing( m_wfr, getWaterFractionType());
  checkNotDecreasing( m_gfr, getGasFractionType());

  initializeFunction();
}

void PipeFlowTableFunction::initializeFunction()
{
  return;
  localIndex constexpr nDims = 4;
  localIndex constexpr nOps = 1;

  FunctionManager * functionManager = &FunctionManager::getInstance();
  // Create nonuniform version of uniformly spaced table
  m_tableFunction = &(dynamicCast< MultivariableNonuniformTableFunction & >( *functionManager->createChild( "MultivariableNonuniformTableFunction", "PipeFlowModel" ) ));

  // find max number if independent vars
  int maxVar = std::max( m_rate.size(), m_whp.size());
  maxVar = std::max( maxVar, m_wfr.size());
  maxVar = std::max( maxVar, m_gfr.size());

  // Copy inputs into formats needed by table function

  integer_array axisPoints( nDims );
  axisPoints[0] = m_rate.size();
  axisPoints[1] = m_whp.size();
  axisPoints[2] = m_wfr.size();
  axisPoints[3] = m_gfr.size();


  array2d< real64 > axisCoordinates( nDims, maxVar );
  for( int i=0; i< m_rate.size(); i++ )
  {
    axisCoordinates[0][i] = m_rate[i];
  }
  for( int i=0; i<m_whp.size(); i++ )
  {
    axisCoordinates[1][i] = m_whp[i];
  }
  for( int i=0; i<m_wfr.size(); i++ )
  {
    axisCoordinates[2][i] = m_wfr[i];
  }
  for( int i=0; i<m_gfr.size(); i++ )
  {
    axisCoordinates[3][i] = m_gfr[i];
  }
  m_tableFunction->setTableCoordinates(nDims,nOps,axisCoordinates,axisPoints);
  m_tableFunction->setTableValues( m_bhp );
  m_tableFunction->initializeFunction();

}
REGISTER_CATALOG_ENTRY( FunctionBase, PipeFlowTableFunction, string const &, Group * const )

} // end of namespace geos
